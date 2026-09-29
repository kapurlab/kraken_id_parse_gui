#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import shutil
import glob
import argparse
import subprocess
import textwrap
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('spawn', True)

from file_setup import Setup, bcolors, Excel_Stats, apply_mpl_style, move_overwrite


class BrackenNoReadsError(Exception):
    """Raised when Bracken produces no output because essentially no reads were
    classified at the species level — usually the sample does not match the
    Kraken database (e.g. the wrong --taxon/--kraken_db for this sample).
    The pipeline treats this as a non-fatal 'nothing found' condition."""


def _available_ram_bytes():
    """How much memory a kraken2 in THIS process's world can actually use.

    Inside a batch/OOD session that is the cgroup limit, not the host's RAM —
    loading an 8 GB database into a 4 GB session doesn't fail politely, the
    kernel SIGKILLs kraken2 ('Killed' with no explanation). Returns the
    smallest applicable bound, or None when nothing is readable (macOS)."""
    bounds = []
    for path in ("/sys/fs/cgroup/memory.max",                 # cgroup v2
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes"):  # cgroup v1
        try:
            raw = open(path).read().strip()
            if raw and raw != "max":
                val = int(raw)
                if 0 < val < 1 << 60:      # v1 reports "unlimited" as a huge number
                    bounds.append(val)
        except (OSError, ValueError):
            pass
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    bounds.append(int(line.split()[1]) * 1024)
                    break
    except (OSError, ValueError):
        pass
    return min(bounds) if bounds else None


def _report_percent_classified(report_path):
    """Percent of reads classified, straight from the kraken2 report.

    Report line: pct, clade reads, taxon reads, rank, taxid, name — the
    'unclassified' line has rank U. Returns None if unparsable."""
    try:
        with open(report_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6 and parts[3].strip() == "U":
                    return 100.0 - float(parts[0].strip())
        return 100.0   # no U line at all -> nothing unclassified
    except (OSError, ValueError):
        return None

plt = apply_mpl_style()

class Kraken_Identification(Setup):
    ''' 
    Assemble reads using Spades assembler.
    Paired or single reads
    '''

    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, directory='kraken', kraken_db=None, influenza=None, debug=None):

        Setup.__init__(self, FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.directory = directory
        if influenza:
            raise ValueError(
                "the legacy Kraken1/JHU influenza mode is no longer supported; "
                "use the IRMA or GenoFLU GUI for influenza analysis"
            )
        self.influenza = False
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

        # A Kraken2 database is exactly three .k2d files. Name what's missing
        # up front — a half-downloaded DB otherwise surfaces as a cryptic
        # kraken2 failure (or worse, a silently empty classification).
        missing = [f for f in ("hash.k2d", "opts.k2d", "taxo.k2d")
                   if not os.path.isfile(os.path.join(kraken_db, f))]
        if missing:
            print(f'\n{bcolors.RED}### Error: {kraken_db} is not a complete Kraken2 database — '
                  f'missing {", ".join(missing)}.{bcolors.ENDC}')
            print('    A DB directory must contain hash.k2d, opts.k2d and taxo.k2d. '
                  'If this DB was downloaded, the download may have been interrupted — re-fetch it.')
            sys.exit(1)
        hash_bytes = os.path.getsize(os.path.join(kraken_db, "hash.k2d"))
        print(f'Kraken2 DB: {kraken_db}  (hash.k2d {hash_bytes/1e9:.1f} GB)')

        # kraken2 loads hash.k2d into RAM. Inside a memory-capped session (an
        # OOD allocation) a too-big DB gets the process SIGKILLed mid-load —
        # seen as a bare 'Killed', or as an inexplicably unclassified run.
        # --memory-mapping reads the DB from disk instead: slower, but correct.
        mem_flag = ''
        avail = _available_ram_bytes()
        if avail is not None and hash_bytes > avail * 0.8:
            mem_flag = ' --memory-mapping'
            print(f'{bcolors.YELLOW}NOTE: DB needs ~{hash_bytes/1e9:.1f} GB but only '
                  f'{avail/1e9:.1f} GB memory is available to this session — running '
                  f'kraken2 with --memory-mapping (slower, avoids the out-of-memory kill).{bcolors.ENDC}')

        if len(FASTQ_list) == 2:
            cmd = f'kraken2 --db {kraken_db}{mem_flag} --threads {cpus} --paired {FASTQ_list[0]} {FASTQ_list[1]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt'
        elif len(FASTQ_list) == 1:
            cmd = f'kraken2 --db {kraken_db}{mem_flag} --threads {cpus} {FASTQ_list[0]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt'
        else:
            cmd = f'kraken2 --db {kraken_db}{mem_flag} --threads {cpus} {FASTA} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt'
        print(f'$ {cmd}')
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            # -9/137 is the kernel OOM/stall killer, not a kraken2 bug.
            hint = (' (SIGKILL — almost always the session memory limit; '
                    'use a bigger allocation)' if rc in (-9, 137) else '')
            print(f'\n{bcolors.RED}### Error: kraken2 exited {rc}{hint}. '
                  f'See its messages above.{bcolors.ENDC}')
            sys.exit(1)

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

        # Say out loud how much was classified — a near-zero number with a
        # plausible-looking Krona is exactly the failure users should not have
        # to discover by squinting at a pie chart.
        pct = _report_percent_classified(report)
        self.percent_classified = pct
        if pct is not None:
            print(f'\nkraken2 classification: {pct:.1f}% of reads classified against this DB')
            if pct < 5.0:
                print(f'{bcolors.YELLOW}WARNING: almost nothing classified. For real samples this '
                      f'usually means the database did not load fully (memory kill / truncated '
                      f'download) or the wrong database was selected — verify hash.k2d\'s size '
                      f'against its source and re-run.{bcolors.ENDC}')

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
        Text-mode Krona from the kraken2 report (kreport2krona.py + ktImportText).
        Deliberately NOT ktImportTaxonomy: text mode needs no Krona taxonomy
        database, so there is nothing to download and nothing to go stale —
        the chart always shows exactly what the kraken2 report says.
        '''
        for cmd in (f'kreport2krona.py --intermediate-ranks -r {report} -o {self.sample_name}.krona',
                    f'ktImportText {self.sample_name}.krona -o {self.sample_name}_{self.date_stamp}_krona.html'):
            print(f'$ {cmd}')
            rc = subprocess.call(cmd, shell=True)
            if rc != 0:
                print(f'\n{bcolors.RED}### Error: Krona step exited {rc}: {cmd.split()[0]}{bcolors.ENDC}')
                sys.exit(1)
        os.remove(f'{self.sample_name}.krona')

        if os.path.exists(f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'):
            self.krona_html = f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'
        else:
            print(f'\n### Error: Krona HTML did not complete')
            sys.exit(1)
        if self.directory:
            move_overwrite(f'{self.sample_name}_{self.date_stamp}_krona.html', self.directory)
            self.krona_html = f'{self.cwd}/{self.directory}/{self.sample_name}_{self.date_stamp}_krona.html'
        return self.krona_html

    @staticmethod
    def _bracken_bin_dir():
        """Where the bracken to run lives: an explicit $BRACKEN_BIN, else a
        sibling 'bracken' conda env (Bracken 2.x can't share the main env on
        osx-64 — old-zlib conflict), else None for plain PATH."""
        override = os.environ.get('BRACKEN_BIN')
        if override and os.path.exists(override):
            return os.path.dirname(os.path.abspath(override))
        cand = os.path.join(os.path.dirname(os.path.dirname(sys.prefix)), 'envs', 'bracken', 'bin')
        if os.path.exists(os.path.join(cand, 'bracken')):
            return cand
        return None

    @classmethod
    def bracken_available(cls):
        """True when bracken can run: $BRACKEN_BIN, the sibling env, or PATH.
        It is optional — it has no macOS/arm64 conda build — so the pipeline
        checks this and skips Bracken instead of failing without it."""
        return cls._bracken_bin_dir() is not None or shutil.which('bracken') is not None

    @classmethod
    def _bracken_prefix(cls):
        """Prepend the sibling env's bin to PATH for the call so its own python
        runs est_abundance.py. Falls back to PATH."""
        bin_dir = cls._bracken_bin_dir()
        if bin_dir:
            return f'PATH="{bin_dir}:$PATH" '
        if not shutil.which('bracken'):
            print('\n### Error: bracken not found. Create it with:\n'
                  '    CONDA_SUBDIR=osx-64 mamba env create -f conda_setup/environment.bracken.yml',
                  file=sys.stderr)
        return ''

    def bracken(self, report, output):
        bracken_txt = f'{self.sample_name}-bracken.txt'
        self.bracken_excel = None
        if os.path.exists(bracken_txt):
            os.remove(bracken_txt)
        rc = subprocess.call(f'{self._bracken_prefix()}bracken -d {self.kraken_db} -i {report} -o {bracken_txt} -r 250', shell=True)
        # Bracken writes no output file when no reads are classified at the
        # species level (or it could not run at all, e.g. a DB without its
        # kmer distribution files). Surface that as a clear, catchable
        # condition instead of a downstream pandas FileNotFoundError — it is
        # abundance refinement, not the identification itself.
        if not os.path.exists(bracken_txt):
            raise BrackenNoReadsError(
                f'Bracken produced no output (exit {rc}) — essentially no reads were classified '
                'at the species level, or the database has no Bracken kmer distribution for '
                '250 bp reads. The Kraken report and Krona graph stand on their own.')
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
    parser.add_argument(
        '-i', '--influenza', action='store_true', dest='influenza', default=False,
        help='deprecated: use the IRMA or GenoFLU GUI for influenza analysis',
    )
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    if args.influenza:
        parser.error(
            "the legacy Kraken1/JHU influenza mode is no longer supported; "
            "use the IRMA or GenoFLU GUI"
        )

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    kraken = Kraken_Identification(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, directory=args.directory, influenza=args.influenza, debug=args.debug)
    kraken.run()
    if args.influenza is False:
        krona_html = kraken.krona_make_graph(kraken.report)
        kraken.bracken(kraken.report, kraken.output)

    print('done')
# Created February 2021 by Tod Stuber
