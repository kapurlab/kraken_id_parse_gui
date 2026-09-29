#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import shutil
import glob
import argparse
import textwrap
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('spawn', True)
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')
from matplotlib import pyplot as plt

from file_setup import Setup, bcolors, Excel_Stats, move_overwrite


class BrackenNoReadsError(Exception):
    """Raised when Bracken produces no output because essentially no reads were
    classified at the species level — usually the sample does not match the
    Kraken database (e.g. the wrong --taxon/--kraken_db for this sample).
    The pipeline treats this as a non-fatal 'nothing found' condition."""


class Kraken_Identification(Setup):
    ''' 
    Assemble reads using Spades assembler.
    Paired or single reads
    '''

    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, directory='kraken', kraken_db=None, influenza=None, debug=None):

        Setup.__init__(self, FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.directory = directory
        self.influenza = influenza
        if influenza:
            self.kraken_db = "/project/bioinformatic_databases/databases/kraken/flu_jhu"
        if kraken_db:
            self.kraken_db = kraken_db
        
    def run(self,):
        self.print_run_time('Kraken')
        kraken_db = self.kraken_db
        cpus = self.cpus
        sample_name = self.sample_name
        FASTQ_list = self.FASTQ_list
        FASTA =  self.FASTA
        cwd = self.cwd
        if self.influenza: #Need to run JHU database using Kraken1
            if len(FASTQ_list) == 2:
                os.system(f'kraken --db {self.kraken_db} --paired {FASTQ_list[0]} {FASTQ_list[1]} > {sample_name}_outputkraken.txt')
            elif len(FASTQ_list) == 1:
                os.system(f'kraken --db {self.kraken_db} {FASTQ_list[0]} > {sample_name}_outputkraken.txt')
            else:
                os.system(f'kraken --db {self.kraken_db} {FASTA} > {sample_name}_outputkraken.txt')
            os.system(f'kraken-report --db {self.kraken_db} {sample_name}_outputkraken.txt > {sample_name}_reportkraken.txt')
            os.system(f'dvl_krakenreport2krona.sh -i {sample_name}_reportkraken.txt -k {self.kraken_db} -t {sample_name}-jhu-output.txt -o {sample_name}-jhu-Krona_id_graphic.html')
            if not os.path.exists(self.directory):
                os.mkdir(self.directory)
            move_overwrite(f'{sample_name}-jhu-Krona_id_graphic.html', self.directory)
            os.remove(f'{sample_name}-jhu-output.txt')
        else:
            if len(FASTQ_list) == 2:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} --paired {FASTQ_list[0]} {FASTQ_list[1]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')
            elif len(FASTQ_list) == 1:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} {FASTQ_list[0]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')
            else:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} {FASTA} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')

        # Nothing downstream can run without Kraken output, so stop with a
        # non-zero status (exit 0 made a failed run look successful to SLURM).
        if os.path.exists(f'{cwd}/{sample_name}_outputkraken.txt'):
            output = f'{cwd}/{sample_name}_outputkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete - check that the Kraken '
                  f'database exists and is valid: {kraken_db}')
            sys.exit(1)
        if os.path.exists(f'{cwd}/{sample_name}_reportkraken.txt'):
            report = f'{cwd}/{sample_name}_reportkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete - check that the Kraken '
                  f'database exists and is valid: {kraken_db}')
            sys.exit(1)

        if self.directory:
            if not os.path.exists(self.directory):
                os.mkdir(self.directory)
            move_overwrite(report, self.directory)
            move_overwrite(output, self.directory)
            self.report = f'{cwd}/{self.directory}/{sample_name}_reportkraken.txt'
            self.output = f'{cwd}/{self.directory}/{sample_name}_outputkraken.txt'
            log_file = open("kraken_log.txt", "a")
            try:
                log_file.write(f'DB used: {os.readlink(self.kraken_db)}')
            except OSError:
                log_file.write(f'DB used: {self.kraken_db}')
            log_file.close()
            move_overwrite("kraken_log.txt", self.directory)

    def krona_make_graph(self, report):
        '''
        kreport2krona.py -r 20-037580-001s_reportkraken.txt -o sample.krona 
        kreport2krona.py --intermediate-ranks -r 20-037580-001s_reportkraken.txt -o sample.krona
        '''
        os.system(f'kreport2krona.py --intermediate-ranks -r {report} -o {self.sample_name}.krona')
        os.system(f'ktImportText {self.sample_name}.krona -o {self.sample_name}_{self.date_stamp}_krona.html')
        os.remove(f'{self.sample_name}.krona')

        if os.path.exists(f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'):
            self.krona_html = f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'
        else:
            print(f'\n### Error: Krona HTML did not complete')
            sys.exit(0)
        if self.directory:
            move_overwrite(f'{self.sample_name}_{self.date_stamp}_krona.html', self.directory)
            self.krona_html = f'{self.cwd}/{self.directory}/{self.sample_name}_{self.date_stamp}_krona.html'
            
    @staticmethod
    def _bracken_prefix():
        """Bracken 2.x can't live in the modern main env (old-zlib conflict), so it
        ships in a sibling 'bracken' conda env. Prepend that env's bin to PATH for
        the call so its own python runs est_abundance.py. Falls back to PATH, and
        honors an explicit $BRACKEN_BIN override."""
        override = os.environ.get('BRACKEN_BIN')
        if override and os.path.exists(override):
            return f'PATH="{os.path.dirname(os.path.abspath(override))}:$PATH" '
        cand = os.path.join(os.path.dirname(os.path.dirname(sys.prefix)), 'envs', 'bracken', 'bin')
        if os.path.exists(os.path.join(cand, 'bracken')):
            return f'PATH="{cand}:$PATH" '
        if not shutil.which('bracken'):
            print('\n### Error: bracken not found. Create it with:\n'
                  '    CONDA_SUBDIR=osx-64 mamba env create -f conda_setup/environment.bracken.yml',
                  file=sys.stderr)
        return ''

    def bracken(self, report, output):
        bracken_txt = f'{self.sample_name}-bracken.txt'
        if os.path.exists(bracken_txt):
            os.remove(bracken_txt)
        os.system(f'{self._bracken_prefix()}bracken -d {self.kraken_db} -i {report} -o {bracken_txt} -r 250')
        # Bracken writes no output file when no reads are classified at the
        # species level. Surface that as a clear, catchable condition instead of
        # a downstream pandas FileNotFoundError.
        if not os.path.exists(bracken_txt):
            raise BrackenNoReadsError(
                'Bracken produced no output — essentially no reads were classified '
                'at the species level. The sample likely does not match the Kraken '
                'database (check --taxon and --kraken_db for this sample).')
        df = pd.read_csv(bracken_txt, sep='\t')
        df.to_excel(f'{self.sample_name}-bracken.xlsx', index=False)
        os.remove(f'{self.sample_name}-bracken.txt')
        self.bracken_excel = f'{os.getcwd()}/{self.sample_name}-bracken.xlsx'
        if self.directory:
            move_overwrite(f'{self.sample_name}-bracken.xlsx', self.directory)
            self.bracken_excel = f'{os.getcwd()}/{self.directory}/{self.sample_name}-bracken.xlsx'

class Bracken_Pie_Charts:
    # Categorical USDA / USWDS palette (official USDA Dark Blue and Dark Green
    # first, then USWDS color tokens and the usda.gov gold) so neighbouring slices
    # stay distinct. Each slice's percentage uses whichever of white / ink text
    # has >= 4.5:1 contrast on it, and the legend lists every taxon with its
    # percentage, so the chart never relies on color alone (USDA Design and Brand
    # guidance).
    PIE_COLORS = ['#002D72', '#005440', '#CEA467', '#005EA2', '#73B3E7',
                  '#C05600', '#3D4551', '#97D4EA', '#936F38', '#A9AEB1']
    WHITE_TEXT = {'#002D72', '#005440', '#005EA2', '#C05600', '#3D4551', '#936F38'}
    UNCLASSIFIED_COLOR = '#DFE1E2'   # USWDS base-lighter
    INK = '#0B2437'                  # usda.gov text color
    # Bracken fractions are shares of the species-level reads, so the remainder
    # is every species under 1% - not unclassified reads (Kraken's unclassified
    # reads are not part of Bracken's species table at all).
    OTHER_LABEL = 'other species (<1% each)'

    def __init__(self, FASTA=False):
        self.FASTA = FASTA

    def run(self, bracken_excel,):
        from coverage_plot import mpl_report_font
        df = pd.read_excel(bracken_excel)
        df = df[df['fraction_total_reads'] > 0.01 ]
        slices = list(zip(df['name'], df['fraction_total_reads']))
        remainder = 1 - df['fraction_total_reads'].sum()
        if remainder > 0:   # never a negative wedge from float rounding
            slices.append((self.OTHER_LABEL, remainder))
        values = [v for _, v in slices]
        total = sum(values)
        colors = [self.UNCLASSIFIED_COLOR if name == self.OTHER_LABEL
                  else self.PIE_COLORS[i % len(self.PIE_COLORS)]
                  for i, (name, _) in enumerate(slices)]
        title = 'Identification of Assembled Scaffolds' if self.FASTA else 'FASTQ Read Identification'

        with plt.rc_context({'font.family': mpl_report_font(), 'text.color': self.INK}):
            fig, ax = plt.subplots(figsize=(9, 5))
            wedges, _labels, pct_texts = ax.pie(
                values, colors=colors, startangle=90, counterclock=False,
                autopct=lambda p: f'{p:.1f}%' if p >= 3 else '', pctdistance=0.72,
                wedgeprops=dict(edgecolor='white', linewidth=1.2), textprops=dict(fontsize=9))
            for color, text in zip(colors, pct_texts):
                text.set_color('white' if color in self.WHITE_TEXT else self.INK)
            ax.legend(wedges, [f'{name} ({v / total:.1%})' for name, v in slices],
                      loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
            ax.set_title(title, fontsize=12, fontweight='bold', color=self.INK)
            ax.axis('equal')
            fig.savefig(f'{os.getcwd()}/bracken_pie.png', format='png', bbox_inches='tight', dpi=150)
            plt.close(fig)
        self.pie_chart = f'{os.getcwd()}/bracken_pie.png'


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        Provide either a single FASTA file, single FASTQ or Paired files.
        Usage:
            kraken_run.py -r1 *_R1*fastg.gz
            kraken_run.py -r1 *_R1*fastg.gz -r2 *_R2*fastq.gz -d
            kraken_run.py -f *fasta

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--FASTA', action='store', dest='FASTA', required=False, help='Provide FASTA file')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, help='Provide R1 FASTQ gz file, or single read')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Provide R2 FASTQ gz file')
    parser.add_argument('-y', '--directory', action='store', dest='directory', required=False, default="kraken", help='Put output to directory')
    parser.add_argument('-i', '--influenza', action='store_true', dest='influenza', default=False, help='Use JHU influenza specific database')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    kraken = Kraken_Identification(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, directory=args.directory, influenza=args.influenza, debug=args.debug)
    kraken.run()
    if args.influenza is False:
        krona_html = kraken.krona_make_graph(kraken.report)
        kraken.bracken(kraken.report, kraken.output)

    print('done')
# Created February 2021 by Tod Stuber