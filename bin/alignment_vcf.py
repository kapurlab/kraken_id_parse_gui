#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import glob
import shutil
import argparse
import textwrap
import zipfile
import pandas as pd
from Bio import SeqIO

from pathlib import Path

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats
from fastq_stats_seqkit import FASTQ_Stats
# from vsnp3_vcf_annotation import VCF_Annotation
from assembly import Assemble
from zero_coverage import Zero_Coverage


def _safe_move(src, dst):
    """shutil.move that silently overwrites an existing destination."""
    src_path = Path(src)
    dst_path = Path(dst)
    actual_dst = dst_path / src_path.name if dst_path.is_dir() else dst_path
    if actual_dst.exists():
        shutil.rmtree(actual_dst) if actual_dst.is_dir() else actual_dst.unlink()
    shutil.move(str(src), str(dst))



class Alignment(Setup):
    ''' 
    '''

    def __init__(self, FASTQ_R1=None, FASTQ_R2=None, reference=None, gbk=None, skip_assembly=None, debug=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, FASTA=reference, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.print_run_time('Align and make VCF file')
        self.reference = self.FASTA
        self.gbk = gbk
        self.skip_assembly = skip_assembly
        self.reference_name = re.sub('[_.].*', '', os.path.basename(reference))

    def run(self,):
        '''
        description
        '''
        fq = FASTQ_Stats(self.FASTQ_R1, self.FASTQ_R2)
        fq.run()
        sample_name = self.sample_name
        reference = self.reference
        samfile = f'{sample_name}.sam'
        all_bamfile = f'{sample_name}_all.bam'
        sorted_bamfile = f'{sample_name}_sorted.bam'
        nodup_bamfile = f'{sample_name}_nodup.bam'
        unfiltered_hapall = f'{sample_name}_unfiltered_hapall.vcf'
        mapfix_hapall = f'{sample_name}_mapfix_hapall.vcf'
        filtered_hapall = f'{sample_name}_filtered_hapall.vcf'
        unmapped_read1 = f'{sample_name}_unmapped_R1.fastq'
        unmapped_read2 = f'{sample_name}_unmapped_R2.fastq'
        unmapped_read = f'{sample_name}_unmapped.fastq'
        zero_coverage_vcf = f'{sample_name}_zc.vcf'
        os.system(f'samtools faidx {reference}')
        os.system(f'picard CreateSequenceDictionary REFERENCE={reference} OUTPUT={reference.rsplit(".", 1)[0]}.dict 2> /dev/null')
        os.system(f'bwa index {reference} 2> /dev/null')
        if self.paired:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {self.FASTQ_R1} {self.FASTQ_R2} > {samfile}')
        else:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {self.FASTQ_R1} > {samfile}')
        os.system(f'samtools view -Sb {samfile} -o {all_bamfile}')
        os.system(f'samtools sort {all_bamfile} -o {sorted_bamfile}')
        os.system(f'samtools index {sorted_bamfile}')
        os.system(f'picard MarkDuplicates INPUT={sorted_bamfile} OUTPUT={nodup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv') # 2> /dev/null')
        dup_metrics_df = pd.read_csv('dup_metrics.csv', delimiter='\t', skiprows=6, nrows=1)
        self.UNPAIRED_READS_EXAMINED = int(dup_metrics_df['UNPAIRED_READS_EXAMINED'])
        self.READ_PAIRS_EXAMINED = int(dup_metrics_df['READ_PAIRS_EXAMINED'])
        self.SECONDARY_OR_SUPPLEMENTARY_RDS = int(dup_metrics_df['SECONDARY_OR_SUPPLEMENTARY_RDS'])
        self.UNMAPPED_READS = int(dup_metrics_df['UNMAPPED_READS'])
        self.UNPAIRED_READ_DUPLICATES = int(dup_metrics_df['UNPAIRED_READ_DUPLICATES'])
        self.READ_PAIR_DUPLICATES = int(dup_metrics_df['READ_PAIR_DUPLICATES'])
        self.READ_PAIR_OPTICAL_DUPLICATES = int(dup_metrics_df['READ_PAIR_OPTICAL_DUPLICATES'])
        self.PERCENT_DUPLICATION	 = float(dup_metrics_df['PERCENT_DUPLICATION'])
        os.system(f'samtools index {nodup_bamfile}')
        chrom_ranges = open("chrom_ranges.txt", 'w')
        for record in SeqIO.parse(reference, "fasta"):
            chrom = record.id
            total_len = len(record.seq)
            min_number = 0
            step = 100000
            if step < total_len:
                for chunk in range(min_number, total_len, step)[1:]:
                    print("{}:{}-{}".format(chrom, min_number, chunk), file=chrom_ranges)
                    min_number = chunk
            print("{}:{}-{}".format(chrom, min_number, total_len), file=chrom_ranges)
        chrom_ranges.close()
        os.system(f'freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {reference} {nodup_bamfile} > {unfiltered_hapall}')
        write_fix = open(mapfix_hapall, 'w+')
        with open(unfiltered_hapall, 'r') as unfiltered:
            for line in unfiltered:
                line = line.strip()
                new_line = re.sub(r';MQM=', r';MQ=', line)
                new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
                print(new_line, file=write_fix)
            write_fix.close()
        # remove clearly poor positions
        os.system(f'vcffilter -f "QUAL > 20" {mapfix_hapall} > {filtered_hapall}')

        unmapped_dir = 'unmapped_reads'
        if not os.path.exists(unmapped_dir):
            os.makedirs(unmapped_dir)
        if self.paired:
            os.system(f'samtools fastq -f4 -1 {unmapped_read1} -2 {unmapped_read2} --reference {reference} --threads 8 {nodup_bamfile} 2> /dev/null' )
            os.system(f'pigz {unmapped_read1}')
            os.system(f'pigz {unmapped_read2}')
            _safe_move(f'{unmapped_read1}.gz', unmapped_dir)
            _safe_move(f'{unmapped_read2}.gz', unmapped_dir)
            self.unmapped_read1 = f'{self.cwd}/{unmapped_dir}/{unmapped_read1}.gz'
            self.unmapped_read2 = f'{self.cwd}/{unmapped_dir}/{unmapped_read2}.gz'
            self.unmapped_read_list = [self.unmapped_read1, self.unmapped_read2]
            if not self.skip_assembly:
                assemble = Assemble(FASTQ_R1=self.unmapped_read1, FASTQ_R2=self.unmapped_read2, debug=True)
        else:
            os.system(f'samtools fastq -f4 -0 {unmapped_read} --reference {reference} --threads 8 {nodup_bamfile} 2> /dev/null' )
            os.system(f'pigz {unmapped_read}')
            _safe_move(f'{unmapped_read}.gz', unmapped_dir)
            self.unmapped_read = f'{self.cwd}/{unmapped_dir}/{unmapped_read}.gz'
            self.unmapped_read_list = [self.unmapped_read]
            if not self.skip_assembly:
                assemble = Assemble(FASTQ_R1=self.unmapped_read, debug=True)

        if not self.skip_assembly:
            assemble.run()
            try:
                _safe_move(assemble.FASTA, f'{assemble.fastq_name}_unmapped.fasta')
                _safe_move(f'{assemble.fastq_name}_unmapped.fasta', unmapped_dir)
                assemble.FASTA = f'{self.cwd}/{unmapped_dir}/{assemble.fastq_name}_unmapped.fasta'
                assemble.stats(assemble.FASTA)
                self.assemble = assemble
                self.unmapped_assemble = assemble.FASTA
                if not self.debug:
                    shutil.rmtree('spades_assembly')
            except (TypeError):
                self.unmapped_assemble = None
                with open(f'{self.cwd}/{unmapped_dir}/failed_assembly', "w") as opened_file:
                    print("see unmapped FASTQ files for troubleshooting", file=opened_file)
                _safe_move('spades_assembly/spades.log', unmapped_dir)
                if not self.debug:
                    shutil.rmtree('spades_assembly')

        if self.gbk:
            self.vcf_annotation = VCF_Annotation(gbk_dir=self.gbk, vcf_file=filtered_hapall)

        self.zero_coverage = Zero_Coverage(FASTA=reference, bam=nodup_bamfile, vcf=filtered_hapall,)
        

        alignment = f'alignment_{self.reference_name}'
        if not os.path.exists(alignment):
            os.makedirs(alignment)
        if os.path.exists(unmapped_dir):
            _safe_move(unmapped_dir, alignment)
        files_grab = []
        for files in ('*_nodup.bam', '*_zc.vcf', '*_nodup.bam.bai', '*_annotated.vcf', '*_filtered_hapall.vcf'):
            files_grab.extend(glob.glob(files))
        for each in files_grab:
            _safe_move(each, alignment)
        _safe_move(reference, alignment)
        _safe_move(f'{reference}.fai', alignment)
        self.reference = f'{self.cwd}/{alignment}/{os.path.basename(reference)}'
        self.nodup_bamfile = f'{self.cwd}/{alignment}/{nodup_bamfile}'
        self.filtered_hapall = f'{self.cwd}/{alignment}/{filtered_hapall}'
        self.zc_vcf = f'{self.cwd}/{alignment}/{self.zero_coverage.zero_coverage_vcf}'
        if self.gbk:
            self.annotated_vcf = f'{self.cwd}/{alignment}/{os.path.basename(self.vcf_annotation.vcf)}'

        temp_dir = './temp'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        files_grab = []
        for files in ('*_unmapped*.fastq.gz', '*_all.bam', '*.bai', '*_mapfix_hapall.vcf', '*_unfiltered_hapall.vcf', '*.sam', '*.amb', '*.ann', '*.bwt', '*.pac', '*.fasta.sa', '*_sorted.bam', '*.dict', 'chrom_ranges.txt', '*.fai', 'dup_metrics.csv'):
            files_grab.extend(glob.glob(files))
        for each in files_grab:
            _safe_move(each, temp_dir)

        if self.debug is False:
            shutil.rmtree(temp_dir)


    def latex(self, tex):
        blast_banner = Banner(f'Read Mapping against {self.reference_name}')
        print(r'\begin{table}[H]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        
        print(r'\begin{tabular}{ l | l | l | l | l | l }', file=tex)
        print(r'Mapped Paired Reads & Mapped Single Reads & Unmapped Reads & Unmapped Percent & \multicolumn{2}{l}{Unmapped Assembled Contigs} \\', file=tex)
        print(r'\hline', file=tex) 
        mapped_reads = self.READ_PAIRS_EXAMINED + self.UNPAIRED_READS_EXAMINED
        total_reads = mapped_reads + self.UNMAPPED_READS
        self.freq_unmapped_reads = self.UNMAPPED_READS / total_reads
        # if self.unmapped_assemble:
        #     print(f'{self.READ_PAIRS_EXAMINED:,} & {self.UNPAIRED_READS_EXAMINED:,} & {self.UNMAPPED_READS:,} & {(self.freq_unmapped_reads*100):,.1f}' + r'\%' + f' & ' + r'\multicolumn{2}{l}{' + f'{self.assemble.contig_count:,}' + r' } \\', file=tex)
        # else:
        print(f'{self.READ_PAIRS_EXAMINED:,} & {self.UNPAIRED_READS_EXAMINED:,} & {self.UNMAPPED_READS:,} & {(self.freq_unmapped_reads*100):,.1f}' + r'\%' + f' & ' + r'\multicolumn{2}{l}{' + f'n/a' + r' } \\', file=tex)
        print(r'\hline', file=tex)
        print(r'\hline', file=tex)
        
        print(r'Duplicate Paired Reads & Duplicate Single Reads & \multicolumn{4}{l}{Duplicate Percent of Mapped Reads} \\', file=tex)
        print(r'\hline', file=tex)
        print(f'{self.READ_PAIR_DUPLICATES:,} & {self.UNPAIRED_READ_DUPLICATES:,} & ' + r'\multicolumn{4}{l}{' + f'{(self.PERCENT_DUPLICATION*100):,.1f}' + r'\%} \\', file=tex)
        print(r'\hline', file=tex)
        print(r'\hline', file=tex)

        print(f'BAM File & Reference Length & Genome with Coverage & Average Depth & No Coverage Bases & Quality SNPs \\\\', file=tex)
        print(r'\hline', file=tex)
        bam = self.zero_coverage.bam.replace('_', '\_')
        print(f'{bam} & {self.zero_coverage.reference_length:,} & {(self.zero_coverage.genome_coverage*100):,.2f}\% & {self.zero_coverage.ave_coverage:,.1f}X & {self.zero_coverage.total_zero_coverage:,} & {self.zero_coverage.good_snp_count:,} \\\\', file=tex)
        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)

        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        excel_dict['Mapped Paired Reads'] = f'{self.READ_PAIRS_EXAMINED:,}'
        excel_dict['Mapped Single Reads'] = f'{self.UNPAIRED_READS_EXAMINED:,}'
        excel_dict['Unmapped Reads'] = f'{self.UNMAPPED_READS:,}'
        excel_dict['Unmapped Percent'] = f'{(self.freq_unmapped_reads*100):,.1f}%'
        # if self.unmapped_assemble:
        #     excel_dict['Unmapped Assembled Contigs'] = f'{self.assemble.contig_count:,}'
        # else:
        excel_dict['Unmapped Assembled Contigs'] = 'n/a'
        excel_dict['Duplicate Paired Reads'] = f'{self.READ_PAIR_DUPLICATES:,}'
        excel_dict['Duplicate Single Reads'] = f'{self.UNPAIRED_READ_DUPLICATES:,}'
        excel_dict['Duplicate Percent of Mapped Reads'] = f'{(self.PERCENT_DUPLICATION*100):,.1f}%'
        excel_dict['BAM/Reference File'] = f'{self.zero_coverage.bam} made with {self.reference_name}'
        excel_dict['Reference Length'] = f'{self.zero_coverage.reference_length:,}'
        excel_dict['Genome with Coverage'] = f'{(self.zero_coverage.genome_coverage*100):,.2f}%'
        excel_dict['Average Depth'] = f'{self.zero_coverage.ave_coverage:,.1f}X'
        excel_dict['No Coverage Bases'] = f'{self.zero_coverage.total_zero_coverage:,}'
        excel_dict['Quality SNPs'] = f'{self.zero_coverage.good_snp_count:,}'


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Usage:
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta
    alignment_vcf.py -r1 *fastq.gz -r *fasta
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta -g *gbk

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=True, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-r', '--reference', action='store', dest='FASTA', required=True, default=None, help="Optional: Provide reference option or FASTA file.  If neither are given, no -r option, then a TB/Brucella/paraTB best reference are searched")
    parser.add_argument('-g', '--gbk', action='store', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file')
    parser.add_argument('-skip_assembly', '--skip_assembly', action='store_true', dest='skip_assembly', help='skip assembly of unmapped reads')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference=args.FASTA, gbk=args.gbk, skip_assembly=args.skip_assembly, debug=args.debug)
    alignment.run()

    #Latex report
    latex_report = Latex_Report(alignment.sample_name)
    alignment.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(alignment.sample_name)
    alignment.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out', '*_all.bam', '*.bai', '*_mapfix_hapall.vcf', '*_unfiltered_hapall.vcf', '*.sam', '*.amb', '*.ann', '*.bwt', '*.pac', '*.fasta.sa', '*_sorted.bam', '*.dict', 'chrom_ranges.txt', 'dup_metrics.csv', '*.fai'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        _safe_move(each, temp_dir)

    # if args.debug is False:
    #     shutil.rmtree(temp_dir)

# Created March 2021 by Tod Stuber
