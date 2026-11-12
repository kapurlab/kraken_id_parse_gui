#!/usr/bin/env python
"""Which sequencer made the reads, and the tools and options that suit them.

The pipeline was built for Illumina reads. Long reads (Oxford Nanopore,
PacBio) need a long-read assembler (Flye, not SPAdes), a long-read mapper
(minimap2, not bwa mem) and variant calling that leaves out indels, because
their errors are mostly small indels: freebayes' own -e 1 (drop reads with
more than one gap) would drop nearly every nanopore read. Every step that
assembles, maps or calls asks this module, so one detection drives them all
and is recorded in the run's provenance.
"""

import gzip
import re
from collections import namedtuple

__all__ = ['Platform', 'detect', 'map_command', 'needs_bwa_index',
           'freebayes_read_options', 'LONG_READ_MEAN']

#: Mean read length above which single-end reads are long reads. Illumina
#: tops out at 301 bp and Ion Torrent near 600; nanopore reads that short are
#: rare, and nanopore amplicon reads (~400-500 bp) carry nanopore read names,
#: which are checked first.
LONG_READ_MEAN = 700

Platform = namedtuple('Platform', 'name long_reads evidence flye_mode minimap2_preset read_group_pl')

_UUID = re.compile(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')
_ION = re.compile(r'[A-Z0-9]{5}:\d+:\d+')
_PACBIO = re.compile(r'm\d+[A-Za-z]?_\d+_\d+/\d+/(ccs|\d+_\d+)\S*')


def _sample(fastq, n=500):
    '''First n records of a FASTQ: (header of the first, read lengths, per-read mean Phred).'''
    header, lengths, quals = '', [], []
    opener = gzip.open if fastq.endswith('.gz') else open
    with opener(fastq, 'rt', errors='replace') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                header = line.rstrip('\n')
            if i % 4 == 1:
                lengths.append(len(line.rstrip('\n')))
            elif i % 4 == 3:
                q = line.rstrip('\n')
                if q:
                    quals.append(sum(ord(c) - 33 for c in q) / len(q))
                if len(quals) >= n:
                    break
    return header, lengths, quals


def detect(fastq_r1, fastq_r2=None):
    '''Platform of the reads, with the evidence for the call.

    Paired reads are Illumina. Single-end reads are named by their headers
    where those say (nanopore UUIDs or runid=, PacBio movie names, Ion Torrent
    RUNID:ROW:COL, also inside an SRA defline), and otherwise by length.
    '''
    if fastq_r2:
        return Platform('illumina', False, 'paired-end reads', None, None, 'ILLUMINA')
    try:
        header, lengths, quals = _sample(fastq_r1)
    except (OSError, EOFError, UnicodeError):
        header, lengths, quals = '', [], []
    tokens = header[1:].split()
    mean_len = sum(lengths) / len(lengths) if lengths else 0
    mean_q = sum(quals) / len(quals) if quals else 0
    if any(_PACBIO.fullmatch(t) for t in tokens):
        hifi = any(_PACBIO.fullmatch(t) and t.split('/')[2].startswith('ccs') for t in tokens)
        return Platform('pacbio', True, f'PacBio read names ({"HiFi/CCS" if hifi else "CLR"})',
                        '--pacbio-hifi' if hifi else '--pacbio-raw',
                        'map-hifi' if hifi else 'map-pb', 'PACBIO')
    nanopore_names = 'runid=' in header or any(_UUID.fullmatch(t) for t in tokens)
    if nanopore_names or mean_len > LONG_READ_MEAN:
        evidence = ('nanopore read names' if nanopore_names
                    else f'mean read length {mean_len:,.0f} bp (long reads; treated as nanopore)')
        # Flye: --nano-hq for reads under ~5% error (Guppy 5+ / Dorado SUP or
        # HAC, Q20 chemistry); --nano-raw for older, noisier reads.
        hq = mean_q >= 13
        return Platform('nanopore', True, f'{evidence}, mean read quality Q{mean_q:.0f}',
                        '--nano-hq' if hq else '--nano-raw', 'map-ont', 'ONT')
    if any(_ION.fullmatch(t) for t in tokens):
        return Platform('iontorrent', False, 'Ion Torrent read names', None, None, 'IONTORRENT')
    return Platform('illumina', False, f'single-end short reads (mean {mean_len:,.0f} bp)',
                    None, None, 'ILLUMINA')


def needs_bwa_index(platform):
    '''bwa needs an index built first; minimap2 indexes on the fly.'''
    return not platform.long_reads


def map_command(platform, reference, fastq_r1, fastq_r2=None, threads=4, sample=None, mark_secondary=False):
    '''argv mapping the reads to reference, SAM on stdout.

    sample adds a read group (ID/SM = sample, PL = the platform);
    mark_secondary is bwa's -M (split hits marked secondary, for Picard).
    '''
    rg = f'@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform.read_group_pl}' if sample else None
    if platform.long_reads:
        cmd = ['minimap2', '-a', '-x', platform.minimap2_preset, '-t', str(threads)]
        if rg:
            cmd += ['-R', rg]
        return cmd + [reference, fastq_r1]
    cmd = ['bwa', 'mem', '-t', str(threads)]
    if mark_secondary:
        cmd.append('-M')
    if rg:
        cmd += ['-R', rg + ('\\tPI:250' if platform.name == 'illumina' else '')]
    return cmd + [reference, fastq_r1] + ([fastq_r2] if fastq_r2 else [])


def freebayes_read_options(platform, short_read_options):
    '''freebayes options for these reads. Long reads: SNPs only (-i no indels,
    -X no MNPs, -u no complex events), since their errors are mostly small
    indels; no per-read gap limit (-e), which would drop nearly every long
    read; and a base-quality floor of at most Q10 (-q), as older nanopore
    reads carry many bases under Q20. Short reads keep the caller's options.'''
    if not platform.long_reads:
        return list(short_read_options)
    keep, opts = [], list(short_read_options)
    while opts:
        opt = opts.pop(0)
        if opt in ('-e', '--read-indel-limit'):
            opts.pop(0)
        elif opt in ('-q', '--min-base-quality'):
            keep += [opt, str(min(10, int(opts.pop(0))))]
        else:
            keep.append(opt)
    return keep + [o for o in ('-i', '-X', '-u') if o not in keep]
