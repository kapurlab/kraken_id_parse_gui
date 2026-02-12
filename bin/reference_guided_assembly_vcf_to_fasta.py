#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import shutil
import glob
import argparse
import textwrap
import pandas as pd
import allel
from Bio import SeqIO
from Bio.Seq import Seq

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats

class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'


class Reference_Guided_Assembly():
    ''' 
    '''

    def __init__(self, FASTA=None, vcf=None, output_name=None, qual=150, map_quality=58, iupac=True, depth=2):
        '''
        Start at class call
        '''

        if output_name:
            pass
        else:
            fasta_name = re.sub('[.].*', '', FASTA)
            vcf_name = re.sub('[.].*', '', vcf)
            output_name = f'{fasta_name}_{vcf_name}_consensus.fasta'

        ambiguous_lookup = {}
        ambiguous_lookup['AG'] = 'R'
        ambiguous_lookup['CT'] = 'Y'
        ambiguous_lookup['GC'] = 'S'
        ambiguous_lookup['AT'] = 'W'
        ambiguous_lookup['GT'] = 'K'
        ambiguous_lookup['AC'] = 'M'
        ambiguous_lookup['GA'] = 'R'
        ambiguous_lookup['TC'] = 'Y'
        ambiguous_lookup['CG'] = 'S'
        ambiguous_lookup['TA'] = 'W'
        ambiguous_lookup['TG'] = 'K'
        ambiguous_lookup['CA'] = 'M'

        df = allel.vcf_to_dataframe(vcf, fields=['variants/CHROM', 'variants/POS', 'variants/QUAL',
                                    'variants/REF', 'variants/ALT', 'variants/AC', 'variants/DP', 'variants/MQ',], alt_number=1)
        try:  # change AC=1 to ambiguous
            for index, row in df.loc[df['AC'] == 1].iterrows():
                df.at[index, 'ALT'] = ambiguous_lookup[row['REF'] + row['ALT']]
        except (KeyError, TypeError) as e:
            print(f'\n\t#####\n\t##### {e}\n')

        good_snps_df = df[(df.QUAL > qual) & (df.MQ > map_quality) & (df.DP > depth) & (
            df.AC == 2) & (df.ALT.astype(str).str.len() == 1)]  # select all to update SNP
        ambiguous_snps_df = df[(df.QUAL > qual) & (df.MQ > map_quality) & (
            df.DP > depth) & (df.AC == 1) & (df.ALT.astype(str).str.len() == 1)]
        caution_df = df[(df.REF.str.len() > 1) | (
            df.ALT.astype(str).str.len() > 1)]
        if not caution_df.empty:
            caution_df.to_excel('CAUTION_SITES.xlsx')
        N_calls_df = df[(df.REF == 'N')]
        N_calls_df['ALT'] = 'N'
        print("")  # debug breakpoint line to check dataframes calls in debug console
        if iupac:
            dfs = pd.concat([good_snps_df, ambiguous_snps_df, N_calls_df])
        else:
            dfs = pd.concat([good_snps_df, N_calls_df])
        updates_df = dfs[['CHROM', 'POS', 'ALT']]

        record_list = []
        for record in SeqIO.parse(FASTA, "fasta"):
            record_df = updates_df[updates_df['CHROM'] == record.id]
            record_dict = dict(zip(record_df['POS'], record_df['ALT']))
            seq_list = list(record.seq)
            for pos, update in record_dict.items():
                seq_list[pos - 1] = update
            record.seq = Seq("".join(seq_list))
            record_list.append(record)
        SeqIO.write(record_list, output_name, "fasta")
        self.output_name = output_name

        print(f'{bcolors.WHITE}Total VCF calls: {df.shape[0]}{bcolors.ENDC}')
        print(
            f'{bcolors.GREEN}Good SNPs: {good_snps_df.shape[0]}{bcolors.ENDC}')
        if iupac:
            print(
                f'{bcolors.BLUE}Ambiguous: {ambiguous_snps_df.shape[0]}{bcolors.ENDC}')
        else:
            print(
                f'Ambiguous (AC=1) calls not applied.  Updated reference not changed at these locations.')
        print(
            f'{bcolors.YELLOW}N positions: {N_calls_df.shape[0]}{bcolors.ENDC}')
        if not caution_df.empty:
            print(
                f'{bcolors.RED}Unreported indels: {caution_df.shape[0]}{bcolors.ENDC}\n')
        else:
            print(f'{bcolors.RED}0 found, which is good{bcolors.ENDC}\n')

        caution_df['QUAL'] = caution_df['QUAL'].map('{:,.1f}'.format)
        caution_df['MQ'] = caution_df['MQ'].map('{:,.1f}'.format)
        caution_df['LEN_DIFF'] = caution_df['ALT'].astype(str).str.len() - \
            caution_df['REF'].str.len()
        caution_df = caution_df[['CHROM', 'POS', 'QUAL',
                                 'REF', 'ALT', 'LEN_DIFF', 'DP', 'AC', 'MQ']]
        self.df = df
        self.good_snps_df = good_snps_df
        self.ambiguous_snps_df = ambiguous_snps_df
        self.caution_df = caution_df

    def latex(self, tex):
        blast_banner = Banner(
            "Sites Not Applied to Consensus -- Additional Verification Required")
        print(r'\begin{table}[H]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' +
              blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(self.caution_df.to_latex(), file=tex)
        print(r'\\', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)

    def excel(self, excel_dict):
        if not self.caution_df.empty:
            excel_dict[
                'Sites Not Applied to Consensus -- Additional Verification Required'] = f'{self.caution_df.shape[0]}'


if __name__ == "__main__":  # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='FASTA',
                        required=True, help='FASTA used to create VCF file')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        required=True, help='Required: VCF file corresponding to FASTA')
    parser.add_argument('-o', '--output_name', action='store', dest='output_name', required=False,
                        help='Output name of updated FASTA.  If not provided combined name of fasta and vcf input will be output.')
    parser.add_argument('-q', '--qual', action='store', dest='qual',
                        required=False, default=150, help='Quality threshold to call a SNP')
    parser.add_argument('-m', '--map_quality', action='store', dest='map_quality',
                        required=False, default=58, help='Map quality threshold to call a SNP')
    parser.add_argument('-d', '--depth', action='store', dest='depth',
                        required=False, default=2, help='Depth threshold to call a SNP')
    parser.add_argument('-i', '--iupac', action='store_true', dest='iupac',
                        default=True, help='report AC=1 as iupac ambiguity')
    args = parser.parse_args()
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    reference_guided_assembly = Reference_Guided_Assembly(
        FASTA=args.FASTA, vcf=args.vcf, output_name=args.output_name, qual=args.qual, map_quality=args.map_quality, depth=args.depth, iupac=args.iupac)

    # Latex report
    latex_report = Latex_Report(reference_guided_assembly.output_name)
    reference_guided_assembly.latex(latex_report.tex)
    latex_report.latex_ending()

    # Excel Stats
    excel_stats = Excel_Stats(reference_guided_assembly.output_name)
    reference_guided_assembly.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    shutil.rmtree(temp_dir)

# Created September 2021 by Tod Stuber
