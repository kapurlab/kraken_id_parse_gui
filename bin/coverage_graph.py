#!/usr/bin/env python3

__version__ = "0.0.1"

import os
import sys
import random
import string
import re
import glob
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')
from matplotlib import cm, colormaps
import pysam
from Bio import SeqIO
import subprocess
from pathlib import Path

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats


class Coverage_Graph(Setup):
    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        super().__init__(FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.reference = FASTA
        try:
            self.reference_name = re.sub('[._].*', '', os.path.basename(FASTA))
        except TypeError:
            print("""
        ########### ERROR ###########
        Does the selected BLAST database contain your target sample?

        If not, you'll need to choose a different BLAST database that includes your target sample.
        For example, if 'nt_viruses' is selected but your target is a bacterium, you'll need to select a different database.
        #############################
        """)
            sys.exit(1)
        self.output_pdf = None
        self.alignment_stats = {}
        self.cmap = colormaps['jet']

    def _process_references(self):
        """Process reference FASTA files for alignment"""
        if self.debug:
            print(f"Processing reference: {self.reference}")
        
        record_list = []
        with open(self.reference) as handle:
            records = SeqIO.parse(handle, "fasta")
            for rec in records:
                rec.id = rec.id.replace(":", "_")
                record_list.append(rec)
        
        temp_reference = f'{self.reference}.temp.fasta'
        SeqIO.write(record_list, temp_reference, "fasta")
        
        if self.debug:
            print(f"Created temporary reference: {temp_reference}")
        return temp_reference

    def _align_reads(self, temp_reference):
        """Align reads to reference using BWA"""
        if self.debug:
            print("Starting read alignment")
            
        sample_name = re.sub('[._].*', '', os.path.basename(self.FASTQ_R1))
        
        # Index reference
        subprocess.run(["bwa", "index", temp_reference], check=True)
        
        # Align reads
        sam_file = f"{sample_name}.sam"
        bam_file = f"{sample_name}.bam"
        sorted_bam = f"{sample_name}.sorted.bam"
        
        print(f"Aligning reads to reference...")
            
        with open(sam_file, 'w') as sam:
            subprocess.run(["bwa", "mem", "-t", str(self.cpus), temp_reference, 
                          self.FASTQ_R1, self.FASTQ_R2],
                         stdout=sam, check=True)
        
        # Convert to BAM and sort
        subprocess.run(["samtools", "view", "-bS", sam_file], 
                      stdout=open(bam_file, 'wb'), check=True)
        subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
        subprocess.run(["samtools", "index", sorted_bam], check=True)
        mapped_reads = subprocess.check_output(["samtools", "view", "-c", "-F", "4", sorted_bam])
        print(f"Number of mapped reads: {mapped_reads.decode().strip()}")
        bam_stats = subprocess.check_output(["samtools", "flagstat", sorted_bam])
        print(f"BAM statistics:\n{bam_stats.decode()}")
        
        if self.debug:
            print(f"Created sorted BAM: {sorted_bam}")
            
        # Clean up intermediate files
        if not self.debug:
            os.remove(sam_file)
            os.remove(bam_file)
            
        return sorted_bam

    def _get_reference_headers(self):
        """Get original headers from reference FASTA"""
        headers = {}
        for record in SeqIO.parse(self.reference, "fasta"):
            headers[record.id.replace(":", "_")] = record.description
        return headers

    def _calculate_alignment_stats(self, sorted_bam, ref_lengths):
        """Calculate alignment statistics for each reference"""
        if self.debug:
            print("Calculating alignment statistics...")
            
        stats = {}
        ref_headers = self._get_reference_headers()
        
        bam = pysam.AlignmentFile(sorted_bam, "rb")
        for ref_id, ref_len in ref_lengths.items():
            # Get coverage array for this reference
            coverage_array = np.zeros(ref_len)
            try:
                for pileup_column in bam.pileup(ref_id):
                    coverage_array[pileup_column.pos] = pileup_column.n
            except Exception as e:
                if self.debug:
                    print(f"Warning: Error processing reference {ref_id}: {str(e)}")
                continue
                
            # Calculate statistics
            covered_bases = np.count_nonzero(coverage_array)
            mean_coverage = np.mean(coverage_array)
            percent_covered = (covered_bases / ref_len) * 100
            
            stats[ref_id] = {
                'header': ref_headers.get(ref_id, ref_id),
                'length': ref_len,
                'mean_coverage': mean_coverage,
                'percent_covered': percent_covered
            }
            
            if self.debug:
                print(f"Stats for {ref_id}:")
                print(f"  Mean coverage: {mean_coverage:.1f}X")
                print(f"  Percent covered: {percent_covered:.1f}%")
                
        bam.close()
        return stats

    def get_coverage_graph(self):
        """Generate coverage graph with alignment statistics"""
        self.print_run_time('Coverage Graph')
        
        # Process references and get headers
        temp_reference = self._process_references()
        sample_name = re.sub('[._].*', '', os.path.basename(self.FASTQ_R1))
        
        # Get reference lengths
        ref_lengths = {}
        for record in SeqIO.parse(temp_reference, "fasta"):
            ref_lengths[record.id] = len(record.seq)
        
        # Align reads
        sorted_bam = self._align_reads(temp_reference)
        
        # Calculate alignment statistics
        self.alignment_stats = self._calculate_alignment_stats(sorted_bam, ref_lengths)
        
        # Get coverage data
        coverage_list = pysam.depth(sorted_bam, split_lines=True)
        if self.debug:
            print(f"Number of coverage lines: {len(list(coverage_list))}")
        list_of_lists = []
        for line in coverage_list:
            chrom, position, depth = line.split('\t')
            list_of_lists.append([chrom, position, depth])
        
        df = pd.DataFrame(list_of_lists, columns=['chrom', 'position', 'depth'])
        
        # Format each segment
        df_list = []
        for chrom, df_chrom in df.groupby('chrom'):
            df = df_chrom[['position', 'depth']]
            df = df.set_index('position')
            df = df.rename(columns={"depth": f'{chrom}'})
            df_list.append(df)
        
        table = pd.concat(df_list, axis=1, sort=False)
        table = table.astype(float)
        
        # Calculate mean coverage
        mean_coverage_dict = {}
        for column in table:
            mean_coverage_dict[column] = table[column].mean()
        
        # Create figure with 3 subplots (2 for coverage, 1 for stats)
        fig = plt.figure(figsize=(12, 10))
        gs = fig.add_gridspec(3, 1, height_ratios=[2, 2, 1.5])
        
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])
        
        # Plot high coverage
        gt_hundredX = [k for k, v in mean_coverage_dict.items() if v > 99]
        gt_hundredX = sorted(gt_hundredX)
        
        rolling_window = 1 if table.shape[0] <= 50000 else int(table.shape[0] / 1500)
        
        try:
            table[table.columns.intersection(gt_hundredX)].rolling(rolling_window).mean().plot(
                logy=True, ax=ax1, cmap=self.cmap,  # Using self.cmap here
                title='Average Depth of Coverage (>100X)'
            )
            ax1.set_ylabel("log depth coverage")
            ax1.legend(loc='lower center', ncol=3, prop={'size': 6})
            ax1.axhline(100, color='r', linestyle='--')
        except TypeError:
            ax1.text(.05, .1, 'No Coverage Data Above 100X', dict(size=12))
        
        # Plot low coverage
        try:
            table[table.columns.difference(gt_hundredX)].rolling(rolling_window).mean().plot(
                logy=False, ax=ax2, cmap=self.cmap,  # Using self.cmap here
                title='Average Depth of Coverage (<100X)'
            )
            ax2.set_xlabel("position")
            ax2.set_ylabel("linear depth coverage")
            ax2.legend(loc='upper right', ncol=3, prop={'size': 6})
            ax2.axhline(100, color='r', linestyle='--')
        except TypeError:
            ax2.text(0.05, .1, 'No Coverage Data Below 100X', dict(size=12))
              
        # Overall title
        fig.suptitle(f'{sample_name} against reference sequences', fontsize=12)
        plt.tight_layout()
        
        # Save plot
        random_part = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
        output_pdf = f'{sample_name}-{random_part}-coverage_graph.pdf'
        plt.savefig(output_pdf, bbox_inches='tight')
        self.output_pdf = output_pdf

        if self.debug:
            print(f"\nCreated coverage graph: {output_pdf}")
        
        # Cleanup temp files if not in debug mode
        if not self.debug:
            for pattern in ['*.fai', '*.sam', '*.bam', '*.bai', '*.amb', '*.ann', 
                          '*.bwt', '*.pac', '*.sa', temp_reference]:
                for f in glob.glob(pattern):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

    def latex(self, tex):
        """Add coverage graph to LaTeX report"""
        if not hasattr(self, 'output_pdf') or self.output_pdf is None:
            print("Warning: Coverage graph not generated. Running analysis...")
            self.get_coverage_graph()
        
        if not self.output_pdf:
            print("Error: Could not generate coverage graph")
            return

        # Add banner and start figure environment
        coverage_banner = Banner("Coverage Analysis")
        print(r'\newpage', file=tex)

        # Coverage banner
        print(r'\begin{figure}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(f'\\includegraphics[scale=1]{{{coverage_banner.banner}}}', file=tex)
        print(r'\end{adjustbox}', file=tex)

        # Add coverage graph
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(f'\\includegraphics[width=\\textwidth]{{{self.output_pdf}}}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\caption{Coverage depth analysis showing read alignment depth across reference sequences.}', file=tex)
        print(r'\end{figure}', file=tex)

        print(r'\vspace{2cm}', file=tex)  # Add more vertical space between figure and table

        # Add alignment statistics table
        print(r'\begin{table}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{l|r|r|r}', file=tex)
        print(r'\hline', file=tex)
        print(r'Reference & Length & Mean Coverage & \% Genome Covered \\', file=tex)
        print(r'\hline', file=tex)

        # Add reference statistics
        for ref_id, stats in self.alignment_stats.items():
            header = stats["header"].replace("_", "\\_").replace("&", "\\&")
            print(f'{header} & {stats["length"]:,} & {stats["mean_coverage"]:.1f}X & {stats["percent_covered"]:.1f}\\% \\\\', file=tex)

        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\caption{Alignment statistics for each reference sequence.}', file=tex)
        print(r'\end{table}', file=tex)

        print(r'\clearpage', file=tex)  # Force a page break after each graph/table pair

def main():
    import argparse
    import textwrap
    
    parser = argparse.ArgumentParser(
        prog='Coverage Graph Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        ---------------------------------------------------------
        Generate coverage graphs and alignment statistics from 
        reference FASTA and FASTQ reads.
        
        Outputs:
        - Coverage depth plots
        - Alignment statistics
        - Per-reference metrics
        ---------------------------------------------------------
        '''))
    
    # Required arguments
    parser.add_argument('-f', '--fasta',
                       required=True,
                       help='Reference FASTA file')
    parser.add_argument('-r1', '--fastq1',
                       required=True,
                       help='Input R1 FASTQ file (can be gzipped)')
    parser.add_argument('-r2', '--fastq2',
                       required=True,
                       help='Input R2 FASTQ file (can be gzipped)')
    
    # Optional arguments
    parser.add_argument('-d', '--debug',
                       action='store_true',
                       help='Keep intermediate files and print debug info')
    parser.add_argument('-o', '--output',
                       help='Output directory (default: current directory)')
    parser.add_argument('-v', '--version',
                       action='version',
                       version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    # Print settings
    print(f"\n{os.path.basename(__file__)} Settings:")
    print(f"  Reference FASTA: {args.fasta}")
    print(f"  FASTQ R1: {args.fastq1}")
    print(f"  FASTQ R2: {args.fastq2}")
    print(f"  Debug mode: {args.debug}")
    print(f"  Output directory: {args.output if args.output else 'current directory'}\n")
    
    # Create output directory if specified
    if args.output:
        os.makedirs(args.output, exist_ok=True)
        os.chdir(args.output)
    
    # Run coverage analysis
    try:
        coverage_graph = Coverage_Graph(
            FASTA=args.fasta,
            FASTQ_R1=args.fastq1,
            FASTQ_R2=args.fastq2,
            debug=args.debug
        )
        coverage_graph.get_coverage_graph()
        
        print(f"\nAnalysis complete!")
        print(f"Output PDF: {coverage_graph.output_pdf}")
        
    except Exception as e:
        print(f"\nError during analysis: {str(e)}", file=sys.stderr)
        if args.debug:
            raise
        sys.exit(1)

if __name__ == "__main__":
    main()