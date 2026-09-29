#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import shutil
import re
import glob
import subprocess
import argparse
import textwrap
import numpy as np
from Bio import SeqIO

from file_setup import Setup, bcolors, Excel_Stats

from fastq_stats_seqkit import FASTQ_Stats

__all__ = ['SPAdesDidNotAssembleFASTA', 'Assemble']

class SPAdesDidNotAssembleFASTA(Exception):
    """Custom exception for when SPAdes assembly fails - Likely due to input reads."""
    pass

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

class Assemble(Setup):
    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, debug=False, fast_mode=False):
        # First call parent class initialization
        super().__init__(FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)

        self.fast_mode = fast_mode
        
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

    def run(self):
        '''
        Run SPAdes assembly, single read assumes ion torrent
        '''
        FASTQ_list = self.FASTQ_list
        cwd = self.cwd
        debug = self.debug

        self.print_run_time('SPAdes')

        # Optimized SPAdes parameters to reduce system load
        if self.fast_mode:
            # Aggressive optimization for slower systems
            max_threads = min(2, max(1, self.cpus // 4))  # Use fewer threads
            memory_limit = 4                              # Use less RAM
            k_mers = "21,33"                             # Use only 2 k-mer sizes (faster)
            print(f"Fast mode: SPAdes will use {max_threads} threads, max {memory_limit}GB RAM, k-mers: {k_mers}")
        else:
            # Standard optimization - balance speed vs system responsiveness
            max_threads = min(4, max(1, self.cpus // 2))  # Use at most 4 threads or half available CPUs
            memory_limit = 8                              # Limit to 8GB RAM
            k_mers = "21,33,55"                          # Standard k-mer sizes
            print(f"Standard mode: SPAdes will use {max_threads} threads, max {memory_limit}GB RAM, k-mers: {k_mers}")

        if len(FASTQ_list) == 2:
            cmd = [
                "spades.py",
                "-1", FASTQ_list[0],
                "-2", FASTQ_list[1],
                "-o", "spades_assembly",
                "-t", str(max_threads),           # Limit CPU threads
                "-m", str(memory_limit),          # Limit memory (GB)
                "-k", k_mers,                     # Use optimized k-mer sizes
                "--only-assembler"                # Skip error correction (faster, less RAM)
            ]
        elif len(FASTQ_list) == 1:
            cmd = [
                "spades.py",
                "-s", FASTQ_list[0],
                "-o", "spades_assembly",
                "-t", str(max_threads),           # Limit CPU threads
                "-m", str(memory_limit),          # Limit memory (GB)
                "-k", k_mers,                     # Use optimized k-mer sizes
                "--only-assembler"                # Skip error correction (faster, less RAM)
            ]
        else:
            print(f'\n### Must have either single or paired read set.\n')
            sys.exit(0)

        try:
            subprocess.run(cmd, capture_output=True, timeout=SPADES_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            print(f'\n### SPAdes did not finish within {SPADES_TIMEOUT_SECONDS // 60} minutes and was '
                  f'terminated. This usually means its internal coverage-model fit failed to converge '
                  f'on this read set rather than a resource problem — see spades_assembly/spades.log.\n')
            raise SPAdesDidNotAssembleFASTA(
                f'SPAdes timed out after {SPADES_TIMEOUT_SECONDS} seconds (coverage model likely failed to converge)')

        if os.path.exists(f'{cwd}/spades_assembly/scaffolds.fasta'):
            shutil.copy2(f'{cwd}/spades_assembly/scaffolds.fasta', f'{cwd}/{self.sample_name}.fasta')
            self.FASTA = f'{cwd}/{self.sample_name}.fasta'
        else:
            print(f'\n### SPAdes did not complete, see log\n')
            raise SPAdesDidNotAssembleFASTA(f'SPAdes assembly fails - Likely due to input reads')
        
        if not debug:
            shutil.rmtree(f'{cwd}/spades_assembly')

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
                'Assembly produced no usable contigs')

        if self.FASTQ_R1:
            self.coverage_title = 'FASTQ calculated mean coverage' #read count * read size / total assembly length
            read_multiplier = 2 if self.paired else 1
            mean_coverage = ((int(self.R1.num_seqs.replace(',', '')) * float(self.R1.avg_len)) * read_multiplier) / total_contig_lengths
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
    
    def excel(self, excel_dict):
        excel_dict['Contig count'] = f'{self.contig_count:,}'
        excel_dict['Contig length counts <|301-999bp|>'] = f'{self.small_contigs_count:,}|{self.mid_size:,}|{self.greater_one_kb_count:,}'
        excel_dict['Longest contig'] = f'{self.longest_contig:,}'
        excel_dict['Total length'] = f'{self.total_contig_lengths:,}'
        excel_dict['N50'] = f'{self.n50:,}'
        excel_dict[f'{self.coverage_title}'] = f'{self.mean_coverage:,.1f}X'

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
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created 2021 by Tod Stuber
