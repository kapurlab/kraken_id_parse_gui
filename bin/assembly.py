#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import json
import hashlib
import heapq
import tempfile
from datetime import datetime
import shutil
import re
import glob
import subprocess
import argparse
import textwrap
import numpy as np
from Bio import SeqIO

from file_setup import Setup, bcolors, Excel_Stats, safe_move

from fastq_stats_seqkit import FASTQ_Stats
import read_platform

__all__ = ['SPAdesDidNotAssembleFASTA', 'Assemble']

class SPAdesDidNotAssembleFASTA(Exception):
    """SPAdes assembly failed. str(e) says why; .reason is 'out_of_memory',
    'timeout', 'no_scaffolds', 'no_contigs' (too few reads to assemble) or
    'error'."""
    def __init__(self, message='', reason='error'):
        super().__init__(message)
        self.reason = reason

#: Hard ceiling on a single SPAdes invocation. SPAdes' internal k-mer coverage
#: model fit (used for graph simplification, independent of --only-assembler)
#: can converge very slowly on read sets with an unusually uniform coverage
#: distribution — observed with shallow, low-diversity viral/Orbivirus
#: segments — doubling its iteration count each retry (32, 64, 128, 256...).
#: It normally still converges, just slowly, so this ceiling is generous
#: rather than tight; it exists only to bound the rare case where it never
#: converges at all, so the pipeline fails with a clear message instead of
#: hanging indefinitely.
SPADES_TIMEOUT_SECONDS = 7200

#: Most read bases handed to SPAdes; a deeper read set is randomly downsampled
#: first. The genome size is not known here (Viruses can be a 10 kb or a 190 kb
#: genome; Mycobacterium tuberculosis complex is 4.4 Mb), so the cap is a base
#: count, not a depth. 500 Mbp is still ~110x of a 4.4 Mb bacterial genome,
#: past the depth where SPAdes assemblies stop improving, so ordinary bacterial
#: runs are never downsampled. Above it, read errors and adapter read-through
#: add false k-mers in proportion to the bases, not the genome, and they fill
#: memory: 1.34 Gbp of ASFV reads (ERR4507142, ~7,000x of 190 kb) ran SPAdes
#: out of 8 GB at K55. 500 Mbp of simulated ASFV reads like them (251 bp, 130-
#: 300 bp inserts, ~2,600x) assembled with a 3.7 GB peak (2.2 GB in fast
#: mode); all 1.34 Gbp, given unlimited memory, gave only sub-210 bp contigs.
SPADES_MAX_BASES = 500_000_000

#: What SPAdes, or the allocator under it, writes when it runs out of memory.
SPADES_OOM_MARKERS = ('Cannot allocate memory', 'unable to allocate OS memory',
                      'mimalloc: error', 'std::bad_alloc')

def spades_failure(returncode, output, log_path, assembler='SPAdes'):
    '''Say why the assembler (SPAdes or Flye) stopped without contigs, from
    its exit code and log.

    Returns (reason, message) for SPAdesDidNotAssembleFASTA. Too few reads and
    too many (out of memory) both end here, and they need opposite fixes.
    '''
    try:
        with open(log_path, errors='replace') as fh:
            output = fh.read() + output
    except OSError:
        pass
    if any(marker in output for marker in SPADES_OOM_MARKERS):
        return ('out_of_memory', f'{assembler} ran out of memory: the read set is too large '
                'for the memory limit, not too small')
    if returncode in (-9, 137) or 'OS return value: -9' in output:
        return ('out_of_memory', f'{assembler} was killed by the system (signal 9), which '
                'usually means it ran out of memory')
    # Exit 0 without scaffolds: SPAdes could not pair the reads. Seen on 1.34
    # Gbp of simulated ASFV reads with no memory limit: 2.2 million contigs
    # under 210 bp, so this is no sign of too few reads.
    if returncode == 0 and 'Insert size was not estimated' in output:
        return ('no_scaffolds', 'SPAdes wrote no scaffolds: it could not estimate '
                'the insert size from the read pairs')
    # SPAdes 4 exits 64 with this error when a few reads give it no k-mer
    # coverage to model; exit 0 without scaffolds is the same shortfall.
    # Flye stops with these when the reads do not overlap enough to assemble
    # (too few, or shorter than its 1 kb minimum overlap).
    if (returncode == 0 or 'Invalid kmer coverage histogram' in output
            or 'No disjointigs were assembled' in output or 'No contigs were assembled' in output):
        return ('no_contigs', f'{assembler} assembled no contigs')
    if 'No space left on device' in output:
        return ('error', f'{assembler} ran out of disk space in the sample folder')
    return ('error', f'{assembler} stopped with exit code {returncode} '
            f'(see {os.path.basename(log_path)})')

#: SPAdes RAM ceiling in GB per mode. -m is the lesser of this and three
#: quarters of the memory the run may use (usable_memory_gb), so a laptop, a
#: container or a SLURM job never promises SPAdes more than it has. Linux
#: enforces -m; macOS refuses SPAdes' setrlimit call, so there it is advice.
SPADES_MEMORY_GB = {'standard': 16, 'fast': 4}

def usable_memory_gb(root='/'):
    '''Memory this run may use, in GB: the least of physical RAM, Linux's
    MemAvailable, the cgroup limit (containers, SLURM jobs; on the job's
    cgroup or any ancestor, v1 or v2) and the SLURM allocation. None when
    none of these can be read. root is for tests.'''
    limits = []
    if root == '/':
        try:
            limits.append(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES'))
        except (AttributeError, ValueError, OSError):
            pass

    def read(path):
        try:
            with open(os.path.join(root, path.lstrip('/'))) as fh:
                return fh.read()
        except OSError:
            return ''

    for line in read('/proc/meminfo').splitlines():
        if line.startswith('MemAvailable:'):
            limits.append(int(line.split()[1]) * 1024)
    for line in read('/proc/self/cgroup').splitlines():
        _, controllers, path = line.split(':', 2)
        if controllers == '':
            base, name = '/sys/fs/cgroup', 'memory.max'
        elif 'memory' in controllers.split(','):
            base, name = '/sys/fs/cgroup/memory', 'memory.limit_in_bytes'
        else:
            continue
        while True:
            value = read(f"{base}{path.rstrip('/')}/{name}").strip()
            if value.isdigit() and int(value) < 2**60:   # v1 "unlimited" is ~2**63
                limits.append(int(value))
            if path in ('', '/'):
                break
            path = os.path.dirname(path.rstrip('/'))
    try:
        limits.append(int(os.environ['SLURM_MEM_PER_NODE']) * 2**20)
    except (KeyError, ValueError):
        pass
    try:
        cpus = int(os.environ.get('SLURM_CPUS_ON_NODE') or os.environ.get('SLURM_CPUS_PER_TASK') or 1)
        limits.append(int(os.environ['SLURM_MEM_PER_CPU']) * 2**20 * cpus)
    except (KeyError, ValueError):
        pass
    return min(limits) / 2**30 if limits else None

def _text_lines(path):
    '''Lines of a plain or gzipped text file, without their newlines.'''
    import gzip
    with (gzip.open if path.endswith('.gz') else open)(path, 'rt', errors='replace') as fh:
        for line in fh:
            yield line.rstrip('\n')

def tool_version(cmd):
    '''First line a tool prints for its version, for the provenance record.'''
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return 'not found'
    lines = (out.stdout.strip() or out.stderr.strip()).splitlines()
    return lines[0] if lines else 'unknown'

def sha256(path):
    '''SHA-256 of a file, for the provenance record.'''
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()

class Assemble(Setup):
    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, debug=False, fast_mode=False):
        # First call parent class initialization
        super().__init__(FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)

        self.fast_mode = fast_mode
        # Set by cap_depth() when SPAdes gets a subsample of the reads, and by
        # trim_reads() to say how (or whether) fastp trimmed them.
        self.downsample_note = None
        self.trim_note = None
        self.temp_reads = []
        # What run() did, for audit: platform, tool versions, every command,
        # each attempt's outcome. Written to <sample>_assembly_provenance.json.
        self.platform = None
        self.method_note = None
        self.provenance = {}
        # Long reads: contigs plus the reads no contig explains, for BLAST
        # (unassembled_reads()); None means BLAST the contigs alone.
        self.blast_query = None
        self.read_rescue_note = None

        # Ensure sample_name is set
        if not hasattr(self, 'sample_name'):
            if FASTQ_R1:
                self.sample_name = re.sub('[_.].*', '', os.path.basename(FASTQ_R1))
            elif FASTA:
                self.sample_name = re.sub('[_.].*', '', os.path.basename(FASTA))
            else:
                self.sample_name = 'sample'

        # Get FASTQ stats if FASTQ files provided
        if FASTQ_R1:
            fastq_stats = FASTQ_Stats(FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
            fastq_stats.run()
            self.R1 = fastq_stats.R1
            self.R2 = fastq_stats.R2

    def read_bases(self):
        '''Bases in the reads handed to Assemble (R1 + R2).'''
        reads = [r for r in (getattr(self, 'R1', None), getattr(self, 'R2', None)) if r]
        return sum(int(r.sum_len.replace(',', '')) for r in reads)

    def cap_depth(self, threads, max_bases=SPADES_MAX_BASES):
        '''
        Return the reads to assemble: all of them, or a random subsample when
        they total more than max_bases. seqkit sample -p keeps or drops
        each record from one random stream, so the same seed on R1 and R2
        keeps the same records in each and mates stay paired.
        '''
        total_bases = self.read_bases()
        if total_bases <= max_bases:
            return self.FASTQ_list
        fraction = max_bases / total_bases
        subsample = [f'{self.cwd}/{self.sample_name}_spades_input_R{i}.fastq.gz'
                     for i in range(1, len(self.FASTQ_list) + 1)]
        self.temp_reads += subsample
        for fastq, out in zip(self.FASTQ_list, subsample):
            cmd = ['seqkit', 'sample', '-p', f'{fraction:.6f}', '-s', '11', '-j', str(threads), '-o', out, fastq]
            self.commands.append(cmd)
            subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
        rows = subprocess.run(['seqkit', 'stats', '-T', *subsample],
                              check=True, capture_output=True, text=True).stdout.splitlines()
        kept = [dict(zip(rows[0].split('\t'), row.split('\t'))) for row in rows[1:]]
        unit = 'pairs' if self.paired else 'reads'
        self.downsample_note = (
            f'assembly used a random {fraction:.1%} of the reads: {int(kept[0]["num_seqs"]):,} of '
            f'{self.R1.num_seqs} {unit}, {sum(int(k["sum_len"]) for k in kept) / 1e6:,.0f} of '
            f'{total_bases / 1e6:,.0f} Mbp (cap {max_bases / 1e6:,.0f} Mbp)')
        print(f'  Deep read set: {self.downsample_note}.', flush=True)
        return subsample

    def trim_reads(self, FASTQ_list, threads):
        '''
        Adapter- and quality-trim the reads for SPAdes with fastp, as SPAdes
        recommends. Paired reads have their adapters found by mate overlap,
        which catches the read-through of inserts shorter than the read
        (ERR4507142: 134 bp inserts on 251 bp reads); 3' tails are cut where
        the mean quality of 4 bases falls under Q20, and reads left under
        30 bp are dropped. Without fastp, or if it fails, SPAdes gets the
        reads untrimmed and the note says so.
        '''
        if self.platform and self.platform.long_reads:
            self.trim_note = (f'{self.platform.name} reads not trimmed here: fastp is for short '
                              'reads, and long-read adapters are removed at basecalling')
            print(f'  {self.trim_note}.', flush=True)
            return FASTQ_list
        if not shutil.which('fastp'):
            self.trim_note = ('reads not trimmed: fastp is not in this environment '
                              '(conda_setup/environment.yml lists it)')
            print(f'  {self.trim_note}.', flush=True)
            return FASTQ_list
        trimmed = [f'{self.cwd}/{self.sample_name}_spades_trimmed_R{i}.fastq.gz'
                   for i in range(1, len(FASTQ_list) + 1)]
        # fastp's JSON QC report is kept in the sample folder as a record.
        report = f'{self.cwd}/{self.sample_name}_fastp.json'
        self.temp_reads += trimmed
        cmd = ['fastp', '-i', FASTQ_list[0], '-o', trimmed[0]]
        if len(FASTQ_list) == 2:
            cmd += ['-I', FASTQ_list[1], '-O', trimmed[1], '--detect_adapter_for_pe']
        cmd += ['--cut_tail', '--length_required', '30', '--thread', str(min(threads, 16)),
                '--json', report, '--html', os.devnull]
        self.commands.append(cmd)
        result = subprocess.run(cmd, capture_output=True, text=True)
        try:
            with open(report) as fh:
                summary = json.load(fh)
            before = summary['summary']['before_filtering']
            after = summary['summary']['after_filtering']
        except (OSError, ValueError, KeyError):
            before = None
        if result.returncode != 0 or not before or not after['total_reads']:
            self.trim_note = f'reads not trimmed: fastp failed (exit {result.returncode})'
            print(f'  {self.trim_note}:', *result.stderr.strip().splitlines()[-5:],
                  sep='\n    ', flush=True)
            return FASTQ_list
        adapters = summary.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0)
        self.trim_note = (
            f'fastp trimmed adapters from {adapters / before["total_reads"]:.1%} of reads and kept '
            f'{after["total_reads"] / before["total_reads"]:.1%} of them '
            f'({after["total_bases"] / 1e6:,.1f} Mbp)')
        print(f'  {self.trim_note}.', flush=True)
        return trimmed

    def unassembled_reads(self, threads, max_reads=5000, min_support=3):
        '''
        Long reads: the reads no contig explains, one per group, for BLAST.

        Flye drops replicons shorter than about 2 kb (its minimum read overlap
        is 1 kb), which for a segmented virus is most segments (BTV/EHDV
        segments 4-10 are 0.8-2 kb). Their nanopore reads are often full
        length, so the reads no contig explains (aligned over less than half
        their length) are grouped by mutual overlap, and the longest read of
        each group of at least min_support stands in for the missing contig
        in the BLAST search that picks references. Writes
        <sample>_unassembled_reads.fasta and <sample>_blast_input.fasta
        (contigs + those reads) and returns the latter, or None.
        '''
        reads_fastq, cwd = self.FASTQ_list[0], self.cwd
        preset = self.platform.minimap2_preset
        explained = set()
        cmd = ['minimap2', '-x', preset, '-t', str(threads), '--secondary=no', self.FASTA, reads_fastq]
        self.provenance.setdefault('read_rescue_commands', []).append(cmd)
        paf = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        for line in paf.splitlines():
            f = line.split('\t')
            if int(f[3]) - int(f[2]) >= 0.5 * int(f[1]):
                explained.add(f[0])
        # Keep only the max_reads longest unexplained reads in memory.
        heap, unexplained, lines = [], 0, _text_lines(reads_fastq)
        for header in lines:              # FASTQ: header, sequence, +, qualities
            seq = next(lines)
            next(lines), next(lines)
            name = header[1:].split()[0]
            if name in explained:
                continue
            unexplained += 1
            if len(heap) < max_reads:
                heapq.heappush(heap, (len(seq), name, seq))
            elif len(seq) > heap[0][0]:
                heapq.heapreplace(heap, (len(seq), name, seq))
        if unexplained < min_support:
            return None
        seqs = {name: seq for _, name, seq in heap}
        longest = sorted(seqs, key=lambda n: len(seqs[n]), reverse=True)
        pool = f'{cwd}/{self.sample_name}_unexplained_reads.fasta'
        with open(pool, 'w') as fh:
            fh.writelines(f'>{n}\n{seqs[n]}\n' for n in longest)
        cmd = ['minimap2', '-x', 'ava-ont' if self.platform.name == 'nanopore' else 'ava-pb',
               '-t', str(threads), pool, pool]
        self.provenance['read_rescue_commands'].append(cmd)
        overlaps = {n: set() for n in longest}
        for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
            f = line.split('\t')
            if f[0] != f[5] and int(f[3]) - int(f[2]) >= 0.5 * min(int(f[1]), int(f[6])):
                overlaps[f[0]].add(f[5])
                overlaps[f[5]].add(f[0])
        if not self.debug:
            os.remove(pool)
        covered, groups = set(), []
        for name in longest:
            if name in covered:
                continue
            group = ({name} | overlaps[name]) - covered
            covered |= group
            if len(group) >= min_support:
                groups.append((name, len(group)))
        if not groups:
            return None
        rescued = f'{cwd}/{self.sample_name}_unassembled_reads.fasta'
        with open(rescued, 'w') as fh:
            for i, (name, support) in enumerate(groups, 1):
                fh.write(f'>unassembled_read_{i} {len(seqs[name])} bp, longest of {support} '
                         f'overlapping {self.platform.name} reads no contig explains ({name})\n{seqs[name]}\n')
        query = f'{cwd}/{self.sample_name}_blast_input.fasta'
        with open(query, 'w') as out:
            for path in (self.FASTA, rescued):
                with open(path) as fh:
                    shutil.copyfileobj(fh, out)
        self.read_rescue_note = (
            f'{len(groups)} unassembled {self.platform.name} read(s) were searched with the contigs '
            f'(labelled unassembled_read_N, each the longest of at least {min_support} overlapping reads '
            f'no contig explains): Flye does not assemble replicons or segments under about 2 kb')
        self.provenance.update(read_rescue=self.read_rescue_note, unassembled_reads={
            'path': rescued, 'sha256': sha256(rescued), 'unexplained_reads': unexplained,
            'groups': [{'read': n, 'support': k} for n, k in groups]})
        print(f'  {self.read_rescue_note}.', flush=True)
        return query

    def run(self):
        '''
        Assemble the reads: SPAdes for short reads (Illumina, Ion Torrent),
        Flye for long reads (nanopore, PacBio). Depth is capped and short
        reads trimmed first; every step is recorded in self.provenance.
        '''
        self.print_run_time('Assembly')
        if len(self.FASTQ_list) not in (1, 2):
            print(f'\n### Must have either single or paired read set.\n')
            sys.exit(0)

        self.platform = platform = read_platform.detect(
            self.FASTQ_list[0], self.FASTQ_list[1] if self.paired else None)
        usable = usable_memory_gb()
        if platform.long_reads:
            # Flye: --meta tolerates the uneven coverage of segmented and
            # amplicon genomes; it sets no memory limit of its own.
            assembler = 'Flye'
            max_threads = min(8, max(1, self.cpus // 2))
            memory_limit = None
            options = [platform.flye_mode, '--meta']
        elif self.fast_mode:
            # Slower systems: fewer threads, less RAM, two k-mer sizes and no
            # read error correction, trading some contiguity for time.
            assembler = 'SPAdes'
            max_threads = min(2, max(1, self.cpus // 4))
            ceiling = SPADES_MEMORY_GB['fast']
            options = ['-k', '21,33', '--only-assembler']
        else:
            # SPAdes' recommendations: --isolate ("highly recommended for
            # high-coverage isolate data", which target reads extracted by
            # Kraken are), BayesHammer read error correction, and k-mer sizes
            # picked from the read length (21-127 for 250 bp reads; 55 was the
            # largest before). On a BTV sample whose segments ran 7X to 4,000X,
            # this kept all 10 segments as contigs >= 300 bp; plain error
            # correction lost the 7X segment 6 (296 bp), the old settings 2.
            assembler = 'SPAdes'
            max_threads = min(4, max(1, self.cpus // 2))
            ceiling = SPADES_MEMORY_GB['standard']
            options = ['--isolate']
        if assembler == 'SPAdes':
            memory_limit = max(1, min(ceiling, int(usable * 0.75))) if usable else ceiling
            if platform.name == 'iontorrent':
                options.append('--iontorrent')   # SPAdes: required for Ion Torrent data
        version = (f"Flye {tool_version(['flye', '--version'])}" if assembler == 'Flye'
                   else tool_version(['spades.py', '--version']))
        print(f'  {assembler} assembling filtered reads (this usually takes a few minutes)...', flush=True)
        print(f'  Reads: {platform.name} ({platform.evidence}).', flush=True)
        print(f"  {'Fast' if self.fast_mode else 'Standard'} mode: {version} with {max_threads} threads"
              + (f', max {memory_limit}GB RAM' if memory_limit else '')
              + (f' of {usable:.0f} GB usable' if usable else '')
              + f", options: {' '.join(options) or 'SPAdes defaults'}", flush=True)

        self.provenance = {
            'started_at': datetime.now().isoformat(timespec='seconds'),
            'read_platform': {'name': platform.name, 'evidence': platform.evidence},
            'input_reads': [{'path': f, 'sha256': sha256(f)} for f in self.FASTQ_list],
            'input_bases': self.read_bases(),
            'assembler': assembler,
            'assembler_version': version,
            'fast_mode': self.fast_mode,
            'threads': max_threads,
            'memory_limit_gb': memory_limit,
            'usable_memory_gb': round(usable, 1) if usable else None,
            'max_bases': SPADES_MAX_BASES,
            'subsample_seed': 11,
            'tool_versions': {name: tool_version(cmd) for name, cmd in (
                ('seqkit', ['seqkit', 'version']), ('fastp', ['fastp', '--version']))},
            'attempts': [],
        }
        max_bases = SPADES_MAX_BASES
        try:
            while True:
                try:
                    self.assemble_once(assembler, max_bases, max_threads, memory_limit, options)
                    self.provenance['attempts'][-1]['outcome'] = 'assembled'
                    if platform.long_reads:
                        self.blast_query = self.unassembled_reads(max_threads)
                    break
                except SPAdesDidNotAssembleFASTA as e:
                    if self.provenance['attempts']:
                        self.provenance['attempts'][-1].update(outcome=e.reason, message=str(e))
                    # Out of memory is too much data for this machine: one more
                    # try with half the reads, never more.
                    if e.reason != 'out_of_memory' or max_bases < SPADES_MAX_BASES:
                        raise
                    max_bases = min(max_bases, self.read_bases()) // 2
                    print(f'  {e}; retrying with {max_bases / 1e6:,.0f} Mbp of reads.', flush=True)
        finally:
            self.method_note = (
                f'{version} ({platform.name} reads: {platform.evidence}; '
                f"options: {' '.join(options) or 'SPAdes defaults'}"
                + (f', -m {memory_limit}' if memory_limit else '') + ')')
            self.provenance.update(
                finished_at=datetime.now().isoformat(timespec='seconds'),
                downsampling=self.downsample_note, trimming=self.trim_note, method=self.method_note)
            if self.FASTA and os.path.exists(self.FASTA):
                self.provenance['contigs'] = {'path': self.FASTA, 'sha256': sha256(self.FASTA)}
            with open(f'{self.cwd}/{self.sample_name}_assembly_provenance.json', 'w') as fh:
                json.dump(self.provenance, fh, indent=2)

        if not self.debug:
            shutil.rmtree(f'{self.cwd}/{self.assembly_dir}')

    def assemble_once(self, assembler, max_bases, threads, memory_limit, options):
        '''One assembler run on the reads, capped at max_bases and (short
        reads) trimmed. Its log is copied to <sample>_<assembler>.log.'''
        cwd = self.cwd
        self.assembly_dir = 'flye_assembly' if assembler == 'Flye' else 'spades_assembly'
        log = f'{cwd}/{self.assembly_dir}/' + ('flye.log' if assembler == 'Flye' else 'spades.log')
        # Contigs left by an earlier run or attempt must not pass for these.
        shutil.rmtree(f'{cwd}/{self.assembly_dir}', ignore_errors=True)
        self.temp_reads = []
        self.commands = []
        self.provenance['attempts'].append({'max_bases': max_bases, 'commands': self.commands})
        workdir = cwd
        try:
            FASTQ_list = self.trim_reads(self.cap_depth(threads, max_bases), threads)
            # SPAdes stops (exit 255) when its folder or reads have whitespace in
            # their paths, as macOS paths often do: assemble in a space-free
            # temporary folder with the reads linked in, and move the result back.
            if any(c.isspace() for c in cwd + ''.join(FASTQ_list)):
                workdir = tempfile.mkdtemp(prefix='assembly_')
                links = []
                for i, fastq in enumerate(FASTQ_list, 1):
                    links.append(os.path.join(workdir, f'R{i}.fastq' + ('.gz' if fastq.endswith('.gz') else '')))
                    os.symlink(os.path.abspath(fastq), links[-1])
                FASTQ_list = links
                print(f'  Assembling in {workdir} (the sample folder path has spaces).', flush=True)
            out_dir = os.path.join(workdir, self.assembly_dir)
            if assembler == 'Flye':
                cmd = ['flye', options[0], FASTQ_list[0], *options[1:], '--out-dir', out_dir,
                       '--threads', str(threads)]
                contigs = f'{cwd}/{self.assembly_dir}/assembly.fasta'
            else:
                if len(FASTQ_list) == 2:
                    reads = ['-1', FASTQ_list[0], '-2', FASTQ_list[1]]
                else:
                    reads = ['-s', FASTQ_list[0]]
                cmd = ['spades.py', *reads, '-o', out_dir, '-t', str(threads),
                       '-m', str(memory_limit), *options]
                contigs = f'{cwd}/{self.assembly_dir}/scaffolds.fasta'
            self.commands.append(cmd)
            # text=True: the tail of the assembler's own output is printed if it fails.
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=SPADES_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            print(f'\n### {assembler} did not finish within {SPADES_TIMEOUT_SECONDS // 60} minutes and was '
                  f'terminated. This usually means its internal coverage-model fit failed to converge '
                  f'on this read set rather than a resource problem — see {os.path.basename(log)}.\n')
            raise SPAdesDidNotAssembleFASTA(
                f'{assembler} did not finish within {SPADES_TIMEOUT_SECONDS // 60} minutes and was stopped '
                f'(its coverage model likely failed to converge)', reason='timeout')
        except FileNotFoundError as e:
            raise SPAdesDidNotAssembleFASTA(f'{e.filename} is not installed in this environment')
        finally:
            if workdir != cwd:
                if os.path.isdir(os.path.join(workdir, self.assembly_dir)):
                    shutil.move(os.path.join(workdir, self.assembly_dir), f'{cwd}/{self.assembly_dir}')
                shutil.rmtree(workdir, ignore_errors=True)
            # Subsampled and trimmed reads are only assembler input.
            if not self.debug:
                for path in self.temp_reads:
                    if os.path.exists(path):
                        os.remove(path)
            # The assembler's own log stays with the results, for the record.
            if os.path.exists(log):
                shutil.copy2(log, f'{cwd}/{self.sample_name}_{assembler.lower()}.log')

        if os.path.exists(contigs) and os.path.getsize(contigs) > 0:
            shutil.copy2(contigs, f'{cwd}/{self.sample_name}.fasta')
            self.FASTA = f'{cwd}/{self.sample_name}.fasta'
            print(f'  {assembler} assembly complete.', flush=True)
        else:
            # Surface why the assembler failed instead of swallowing it.
            print(f'\n### {assembler} did not complete (exit {result.returncode}). Last output:', flush=True)
            tail = (result.stderr or result.stdout or '').strip().splitlines()[-15:]
            for line in tail:
                print(f'    {line}', flush=True)
            reason, message = spades_failure(result.returncode, result.stdout + result.stderr, log, assembler)
            raise SPAdesDidNotAssembleFASTA(message, reason=reason)

    def stats(self, FASTA=None):
        '''
        description
        '''
        records = list(SeqIO.parse(FASTA, "fasta"))
        cov_length_list=[]
        contig_count = 0
        coverage_list=[]
        length_list=[]
        small_contigs=[]
        greater_one_kb=[]
        mid_size = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
            except IndexError:
                coverage_value = 1
            coverage_value = int(float(coverage_value))
            cov_length_list.append({'name': rec.description, 'cov': coverage_value, 'length': len(rec)})
            coverage_list.append(coverage_value)
            length_list.append(len(rec))
            contig_count += 1
            if len(rec) <= 300:
                small_contigs.append(len(rec))
            elif len(rec) >= 1000:
                greater_one_kb.append(len(rec))
            else:
                mid_size.append(len(rec))
        total_contig_lengths = int(sum(length_list))

        # No usable contigs (empty or headers-only scaffolds). Treat this the
        # same as "SPAdes did not assemble" — the caller already handles that
        # exception gracefully (writes the report, shows detected taxa) rather
        # than crashing on the divisions / empty-array reductions below.
        if not length_list or total_contig_lengths <= 0:
            raise SPAdesDidNotAssembleFASTA(
                'Assembly produced no usable contigs', reason='no_contigs')

        if self.FASTQ_R1:
            # All extracted reads, even when SPAdes got a subsample (cap_depth):
            # this is the sample's depth of target reads, and a capped value
            # would read the same for every deep sample. downsample_note,
            # printed and in its own Excel column, says what SPAdes used.
            # Read bases (R1 + R2, exact) / total assembly length. Reads x mean
            # length crashed on long reads: seqkit prints "4,968" from 1 kb up.
            self.coverage_title = 'FASTQ calculated mean coverage'
            mean_coverage = self.read_bases() / total_contig_lengths
        else:
            # No FASTQs so calculate mean coverage via SPAdes reportings.
            normalized_list = []
            for rec in records:
                header = rec.description
                try:
                    coverage_value = header.split('_')[5]
                except IndexError:
                    coverage_value = 1
                coverage_value = int(float(coverage_value))
                normalized_list.append((len(rec) / total_contig_lengths) * coverage_value)
            self.coverage_title = 'SPAdes calculated mean coverage' 
            mean_coverage = sum(normalized_list)

        #N50 calculation
        all_len = sorted(length_list, reverse=True)
        csum = np.cumsum(all_len)
        n2 = int(sum(length_list)/2)
        csumn2 = min(csum[csum >= n2])
        ind = np.where(csum == csumn2)
        self.n50 = all_len[int(ind[0])] # n50 smallest size contig which, along with the larger contigs, contain half of sequence of a particular genome
        self.l50 = int(ind[0][0]) + 1 # l50 smallest number of contigs whose length sum makes up half of genome
        self.greater_one_kb_count = len(greater_one_kb)
        self.longest_contig = int(max(length_list))
        self.small_contigs_count = len(small_contigs)
        self.mid_size = len(mid_size)
        self.contig_count = contig_count
        self.total_contig_lengths = total_contig_lengths
        self.mean_coverage = mean_coverage


        print(f'\t     Contig count: {bcolors.YELLOW}{self.contig_count:,}{bcolors.ENDC}, \n \
            Contig length counts <|301-999bp|>: {bcolors.RED}{self.small_contigs_count:,}{bcolors.ENDC}|{bcolors.BLUE}{self.mid_size:,}{bcolors.ENDC}|{bcolors.GREEN}{self.greater_one_kb_count:,}{bcolors.ENDC}, \n \
            Longest contig: {bcolors.GREEN}{self.longest_contig:,}{bcolors.ENDC}, \n \
            Total length: {bcolors.WHITE}{self.total_contig_lengths:,}{bcolors.ENDC}, \n \
            N50: {bcolors.PURPLE}{self.n50:,}{bcolors.ENDC}, \n \
            {self.coverage_title}: {bcolors.YELLOW}{self.mean_coverage:,.1f}X{bcolors.ENDC}\n')
        if self.downsample_note:
            print(f'\t     Coverage counts all reads; {self.downsample_note}\n')

    def excel(self, excel_dict):
        excel_dict['Contig count'] = f'{self.contig_count:,}'
        excel_dict['Contig length counts <|301-999bp|>'] = f'{self.small_contigs_count:,}|{self.mid_size:,}|{self.greater_one_kb_count:,}'
        excel_dict['Longest contig'] = f'{self.longest_contig:,}'
        excel_dict['Total length'] = f'{self.total_contig_lengths:,}'
        excel_dict['N50'] = f'{self.n50:,}'
        excel_dict[f'{self.coverage_title}'] = f'{self.mean_coverage:,.1f}X'
        spades_input = '; '.join(n for n in (self.downsample_note, self.trim_note) if n)
        if spades_input:
            excel_dict['Assembly input reads'] = spades_input

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')
    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=False, default=None, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-f', '--fasta', action='store', dest='FASTA', default=None, help='provide assembly if just stats are needed.  Assumed SPAdes assembly')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    assemble = Assemble(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2)
    if args.FASTQ_R1:
        assemble.run()
        assemble.stats(assemble.FASTA)
    elif args.FASTA:
        assemble.stats(args.FASTA)
    else:
        print('### Error: Provide FASTQ or FASTA file.  See usda_assembly.py -h for option')

    #Excel Stats
    excel_stats = Excel_Stats(assemble.sample_name)
    assemble.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        safe_move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created 2021 by Tod Stuber
