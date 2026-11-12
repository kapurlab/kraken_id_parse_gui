#!/usr/bin/env python3
"""Regression tests for assembly input handling and read-platform routing.

Pins what assembly.py and read_platform.py decide, without running SPAdes or
Flye:

  1. Why an assembly failed (spades_failure): out of memory, killed, too few
     reads, no scaffolds, disk full. Too few reads and too many used to read
     the same ("not enough reads").
  2. How much memory a run may give the assembler (usable_memory_gb) on
     Linux, in a cgroup v1 or v2 container or SLURM job, and from SLURM's
     environment, using fake /proc and /sys trees.
  3. Which sequencer made the reads (read_platform.detect), and the mapper
     and freebayes options each platform gets.

Run:  python bin/test_assembly_platform.py
Exits non-zero on any failure.
"""
import gzip
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import assembly            # noqa: E402
import read_platform       # noqa: E402

FAILED = []


def check(label, got, want):
    if got == want:
        print(f'[PASS] {label}')
    else:
        FAILED.append(label)
        print(f'[FAIL] {label}: got {got!r}, want {want!r}')


def failure_cases():
    oom = ('mimalloc: error: unable to allocate OS memory (67108864 bytes, error code: 12 '
           '[Cannot allocate memory])\n== Error ==  system call for: ... OS return value: 68\n')
    cases = [
        ('mimalloc OOM, exit 68', (68, oom), 'out_of_memory'),
        ('std::bad_alloc', (1, "terminate called after throwing an instance of 'std::bad_alloc'"), 'out_of_memory'),
        ('killed by signal 9', (-9, ''), 'out_of_memory'),
        ('spades.py reports child killed', (255, 'finished abnormally, OS return value: -9'), 'out_of_memory'),
        ('pairs unusable, exit 0', (0, 'WARN Insert size was not estimated for any of the paired libraries'), 'no_scaffolds'),
        ('too few reads (SPAdes 4 exit 64)', (64, 'ERROR Invalid kmer coverage histogram'), 'no_contigs'),
        ('Flye: no disjointigs', (1, 'ERROR: No disjointigs were assembled - please check if the read type and genome size parameters are correct'), 'no_contigs'),
        ('exit 0, nothing written', (0, ''), 'no_contigs'),
        ('disk full', (1, 'No space left on device'), 'error'),
        ('anything else', (12, 'boom'), 'error'),
    ]
    for label, (rc, out), want in cases:
        check(f'spades_failure: {label}', assembly.spades_failure(rc, out, '/nonexistent')[0], want)
    check('spades_failure names the assembler',
          assembly.spades_failure(-9, '', '/nonexistent', 'Flye')[1].startswith('Flye '), True)


def write(root, path, text):
    full = os.path.join(root, path.lstrip('/'))
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, 'w') as fh:
        fh.write(text)


def memory_cases():
    gib = 2**30
    saved = {k: os.environ.pop(k, None) for k in ('SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU',
                                                  'SLURM_CPUS_ON_NODE', 'SLURM_CPUS_PER_TASK')}
    try:
        with tempfile.TemporaryDirectory() as root:
            check('memory: nothing readable', assembly.usable_memory_gb(root), None)
            write(root, '/proc/meminfo', f'MemTotal: {64 * gib // 1024} kB\nMemAvailable: {40 * gib // 1024} kB\n')
            check('memory: Linux MemAvailable', assembly.usable_memory_gb(root), 40)
            # cgroup v2: the job's own limit is unlimited, its parent's is 12 GiB.
            write(root, '/proc/self/cgroup', '0::/slurm/job_7/step_0\n')
            write(root, '/sys/fs/cgroup/slurm/job_7/step_0/memory.max', 'max\n')
            write(root, '/sys/fs/cgroup/slurm/job_7/memory.max', f'{12 * gib}\n')
            check('memory: cgroup v2 limit on an ancestor', assembly.usable_memory_gb(root), 12)
        with tempfile.TemporaryDirectory() as root:
            write(root, '/proc/meminfo', f'MemAvailable: {40 * gib // 1024} kB\n')
            write(root, '/proc/self/cgroup', '5:cpuacct,cpu:/docker/abc\n4:memory:/docker/abc\n')
            write(root, '/sys/fs/cgroup/memory/docker/abc/memory.limit_in_bytes', f'{6 * gib}\n')
            write(root, '/sys/fs/cgroup/memory/memory.limit_in_bytes', '9223372036854771712\n')
            check('memory: cgroup v1 container, root unlimited', assembly.usable_memory_gb(root), 6)
            os.environ['SLURM_MEM_PER_CPU'] = '1024'
            os.environ['SLURM_CPUS_ON_NODE'] = '4'
            check('memory: SLURM --mem-per-cpu x CPUs', assembly.usable_memory_gb(root), 4)
            os.environ['SLURM_MEM_PER_NODE'] = '2048'
            check('memory: SLURM --mem', assembly.usable_memory_gb(root), 2)
        check('memory: this machine is readable', assembly.usable_memory_gb() is not None, True)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def fastq(directory, name, header, seq_len, n=20):
    path = os.path.join(directory, name)
    with gzip.open(path, 'wt') as fh:
        for i in range(n):
            fh.write(f'@{header.format(i=i)}\n{"A" * seq_len}\n+\n{"5" * seq_len}\n')
    return path


def platform_cases():
    with tempfile.TemporaryDirectory() as d:
        illumina = 'M01234:12:000000000-ABCDE:1:1101:15589:{i} 1:N:0:1'
        cases = [
            ('paired Illumina', fastq(d, 'p1.fq.gz', illumina, 251), fastq(d, 'p2.fq.gz', illumina, 251), 'illumina'),
            ('single-end Illumina', fastq(d, 's.fq.gz', illumina, 150), None, 'illumina'),
            ('Ion Torrent names', fastq(d, 'ion.fq.gz', 'ZJ7HH:00009:{i:05d}', 220), None, 'iontorrent'),
            ('Ion Torrent in an SRA defline', fastq(d, 'sra.fq.gz', 'SRR1.{i} ZJ7HH:00009:00029 length=220', 220), None, 'iontorrent'),
            ('nanopore UUID + runid=', fastq(d, 'ont.fq.gz', '0d6f0f5e-1c2b-4f6e-9a8b-0123456789ab runid=abc read={i} ch=5', 4000), None, 'nanopore'),
            ('nanopore amplicons (short, named)', fastq(d, 'amp.fq.gz', '0d6f0f5e-1c2b-4f6e-9a8b-0123456789ab runid=abc', 450), None, 'nanopore'),
            ('long reads without names', fastq(d, 'long.fq.gz', 'read{i}', 3000), None, 'nanopore'),
            ('PacBio HiFi', fastq(d, 'hifi.fq.gz', 'm64011_190830_220126/{i}/ccs', 15000), None, 'pacbio'),
        ]
        for label, r1, r2, want in cases:
            check(f'detect: {label}', read_platform.detect(r1, r2).name, want)
        hifi = read_platform.detect(os.path.join(d, 'hifi.fq.gz'))
        check('PacBio HiFi -> Flye --pacbio-hifi, minimap2 map-hifi', (hifi.flye_mode, hifi.minimap2_preset),
              ('--pacbio-hifi', 'map-hifi'))
        ont = read_platform.detect(os.path.join(d, 'ont.fq.gz'))
        check('nanopore Q20 -> Flye --nano-hq', ont.flye_mode, '--nano-hq')
        check('nanopore maps with minimap2 map-ont, read group PL:ONT',
              read_platform.map_command(ont, 'ref.fa', 'r.fq', threads=8, sample='S1'),
              ['minimap2', '-a', '-x', 'map-ont', '-t', '8', '-R', '@RG\\tID:S1\\tSM:S1\\tPL:ONT', 'ref.fa', 'r.fq'])
        ilm = read_platform.detect(os.path.join(d, 'p1.fq.gz'), os.path.join(d, 'p2.fq.gz'))
        check('Illumina maps with bwa mem as before',
              read_platform.map_command(ilm, 'ref.fa', 'r1.fq', 'r2.fq', threads=8, sample='S1', mark_secondary=True),
              ['bwa', 'mem', '-t', '8', '-M', '-R', '@RG\\tID:S1\\tSM:S1\\tPL:ILLUMINA\\tPI:250', 'ref.fa', 'r1.fq', 'r2.fq'])
        check('Illumina freebayes options unchanged',
              read_platform.freebayes_read_options(ilm, ['-E', '-1', '-e', '1', '-u']), ['-E', '-1', '-e', '1', '-u'])
        check('nanopore freebayes: SNPs only, no -e gap limit',
              read_platform.freebayes_read_options(ont, ['-E', '-1', '-e', '1', '-u']), ['-E', '-1', '-u', '-i', '-X'])
        check('nanopore freebayes: base quality floor Q10',
              read_platform.freebayes_read_options(ont, ['-m', '30', '-q', '20']), ['-m', '30', '-q', '10', '-i', '-X', '-u'])


if __name__ == '__main__':
    failure_cases()
    memory_cases()
    platform_cases()
    print(f'\n{len(FAILED)} failed' if FAILED else '\nall passed')
    sys.exit(1 if FAILED else 0)
