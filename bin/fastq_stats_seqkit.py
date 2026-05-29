#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import subprocess
import shutil
import glob
import argparse
import textwrap
import humanize

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats, safe_move

class FASTQ_Container:
    """Provide nested dot notation to object for each read with stats"""
    def __init__(self, file_name, file_size, read_format, seq_type, num_seqs, sum_len, min_len, avg_len, max_len, Q1, Q2, Q3, sum_gap, N50, passQ20, passQ30, read_quality_average):
        self.file_name = file_name
        self.file_size = file_size
        self.read_format = read_format
        self.seq_type = seq_type
        self.num_seqs = num_seqs
        self.sum_len = sum_len
        self.min_len = min_len
        self.avg_len = avg_len
        self.max_len = max_len
        self.Q1 = Q1
        self.Q2 = Q2
        self.Q3 = Q3
        self.sum_gap = sum_gap
        self.N50 = N50
        self.passQ20 = passQ20
        self.passQ30 = passQ30
        self.read_quality_average = read_quality_average

class FASTQ_Stats(Setup):
    def __init__(self, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        # Initialize the parent class
        super().__init__(FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        
        # Initialize R1 and R2 to None
        self.R1 = None
        self.R2 = None

        # Store file paths directly
        self.fastq_r1_path = FASTQ_R1
        self.fastq_r2_path = FASTQ_R2

        if self.debug:
            print(f"Debug: FASTQ_Stats initialized with:")
            print(f"R1 path: {self.fastq_r1_path}")
            print(f"R2 path: {self.fastq_r2_path}")

    def run(self):
        """Process FASTQ files and collect statistics"""
        files_to_process = {}
        
        if self.fastq_r1_path:
            files_to_process['R1'] = self.fastq_r1_path
        if self.fastq_r2_path:
            files_to_process['R2'] = self.fastq_r2_path

        for read_type, filepath in files_to_process.items():
            if not os.path.exists(filepath):
                print(f"Warning: File not found {filepath}")
                continue

            file_size = humanize.naturalsize(os.path.getsize(filepath))
            
            # Run seqkit stats
            subprocess.run(["seqkit", "stats", "-a", "-o", "temp_fastq_seqkit_stats.txt", filepath], 
                         stderr=subprocess.DEVNULL)
            
            try:
                with open('temp_fastq_seqkit_stats.txt', 'r') as fopen:
                    last_line = fopen.readlines()[-1].split()
                    file_name = last_line[0]
                    read_format = last_line[1]
                    seq_type = last_line[2]
                    num_seqs = last_line[3]
                    sum_len = last_line[4]
                    min_len = last_line[5]
                    avg_len = last_line[6]
                    max_len = last_line[7]
                    Q1 = last_line[8]
                    Q2 = last_line[9]
                    Q3 = last_line[10]
                    sum_gap = last_line[11]
                    N50 = last_line[12]
                    passQ20 = last_line[13]
                    passQ30 = last_line[14]

                os.remove('temp_fastq_seqkit_stats.txt')

                # Average read quality
                cmd = f'seqkit fx2tab {filepath} -l -q -n -i -H | ' + r"awk '{sum+=$3} END{print sum/(NR-1)}'"
                ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                output = ps.communicate()[0]
                read_quality_average = output.decode("utf-8").strip()

                # Create container and set attribute
                container = FASTQ_Container(
                    file_name, file_size, read_format, seq_type, num_seqs,
                    sum_len, min_len, avg_len, max_len, Q1, Q2, Q3,
                    sum_gap, N50, passQ20, passQ30, read_quality_average
                )
                
                setattr(self, read_type, container)

            except Exception as e:
                print(f"Error processing {filepath}: {str(e)}")
                if self.debug:
                    raise
        
        # Ensure R2 is None for single-end data
        if len(files_to_process) == 1:
            self.R2 = None

    def latex(self, tex):
        blast_banner = Banner("FASTQ Quality")
        print(r'\begin{table}[H]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\small', file=tex)
        print(r'\begin{tabular}{ l | l | l }', file=tex)
        
        if self.R2:
            print('Filename & ' + os.path.basename(self.R1.file_name).replace("_", "\_") + ' & ' + os.path.basename(self.R2.file_name).replace("_", "\_") + ' \\\\', file=tex)
            print(r'\hline', file=tex)
            print(f'File Size & {self.R1.file_size} & {self.R2.file_size} \\\\', file=tex)
            print(f'Q30 Passing & {self.R1.passQ30}\% & {self.R2.passQ30}\% \\\\', file=tex)
            print(f'Mean Read Score & {float(self.R1.read_quality_average):0.1f} & {float(self.R2.read_quality_average):0.1f} \\\\', file=tex)
            print(f'Average Read Length & {self.R1.avg_len} & {self.R1.avg_len} \\\\', file=tex)
        else:
            print('Filename & ' + os.path.basename(self.R1.file_name).replace("_", "\_") + ' & Read 2 \\\\', file=tex)
            print(r'\hline', file=tex)
            print(f'File Size & {self.R1.file_size} & N/A \\\\', file=tex)
            print(f'Q30 Passing & {self.R1.passQ30}\% & N/A \\\\', file=tex)
            print(f'Mean Read Score & {float(self.R1.read_quality_average):0.1f} & N/A \\\\', file=tex)
            print(f'Average Read Length & {self.R1.avg_len} & N/A \\\\', file=tex)
        
        print(r'\hline', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\vspace{0.1 mm}', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        excel_dict['FASTQ_R1'] = os.path.basename(self.R1.file_name)
        excel_dict['R1 File Size'] = self.R1.file_size
        excel_dict['R1 Read Count'] = self.R1.num_seqs
        excel_dict['R1 Length Sum'] = self.R1.sum_len
        excel_dict['R1 Min Length'] = self.R1.min_len
        excel_dict['R1 Ave Length'] = self.R1.avg_len
        excel_dict['R1 Max Length'] = self.R1.max_len
        excel_dict['R1 Passing Q20'] = f'{self.R1.passQ20}%'
        excel_dict['R1 Passing Q30'] = f'{self.R1.passQ30}%'
        excel_dict['R1 Read Quality Ave'] = f'{self.R1.read_quality_average}'
        
        if self.R2:
            excel_dict['FASTQ_R2'] = os.path.basename(self.R2.file_name)
            excel_dict['R2 File Size'] = self.R2.file_size
            excel_dict['R2 Read Count'] = self.R2.num_seqs
            excel_dict['R2 Length Sum'] = self.R2.sum_len
            excel_dict['R2 Min Length'] = self.R2.min_len
            excel_dict['R2 Ave Length'] = self.R2.avg_len
            excel_dict['R2 Max Length'] = self.R2.max_len
            excel_dict['R2 Passing Q20'] = f'{self.R2.passQ20}%'
            excel_dict['R2 Passing Q30'] = f'{self.R2.passQ30}%'
            excel_dict['R2 Read Quality Ave'] = f'{self.R2.read_quality_average}'

def main():
    parser = argparse.ArgumentParser(
        prog='PROG', 
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=textwrap.dedent('''
        ---------------------------------------------------------
        seqkit used to calculate FASTQ stats
        https://bioinf.shenwei.me/seqkit/

        Usage:
        usda_fastq_stats_seqkit.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz #paired
        usda_fastq_stats_seqkit.py -r1 *fastq.gz #single
        '''),
        epilog='---------------------------------------------------------'
    )

    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=False, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    fastq_stats = FASTQ_Stats(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=args.debug)
    fastq_stats.run()

    print(f'\t R1 File Size: {bcolors.WHITE}{fastq_stats.R1.file_size}{bcolors.ENDC} \n \
        R1 Passing Q30: {bcolors.WHITE}{fastq_stats.R1.passQ30}{bcolors.ENDC} \n \
        R1 Mean Read Score: {bcolors.WHITE}{fastq_stats.R1.read_quality_average}{bcolors.ENDC} \n \
        R1 Average Read Length: {bcolors.WHITE}{fastq_stats.R1.avg_len}{bcolors.ENDC} \n')

    if args.FASTQ_R2:
        print(f'\t R2 File Size: {bcolors.WHITE}{fastq_stats.R2.file_size}{bcolors.ENDC} \n \
        R2 Passing Q30: {bcolors.WHITE}{fastq_stats.R2.passQ30}{bcolors.ENDC} \n \
        R2 Mean Read Score: {bcolors.WHITE}{fastq_stats.R2.read_quality_average}{bcolors.ENDC} \n \
        R2 Average Read Length: {bcolors.WHITE}{fastq_stats.R2.avg_len}{bcolors.ENDC} \n')

    #Latex report
    latex_report = Latex_Report(fastq_stats.sample_name)
    fastq_stats.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(fastq_stats.sample_name)
    fastq_stats.excel(excel_stats.excel_dict)
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

if __name__ == "__main__":
    main()