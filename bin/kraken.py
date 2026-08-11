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

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats, apply_mpl_style


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
cmap = plt.get_cmap('coolwarm')

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

        if os.path.exists(f'{cwd}/{sample_name}_outputkraken.txt'):
            output = f'{cwd}/{sample_name}_outputkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(1)
        if os.path.exists(f'{cwd}/{sample_name}_reportkraken.txt'):
            report = f'{cwd}/{sample_name}_reportkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(1)

        # Say out loud how much was classified — a near-zero number with a
        # plausible-looking Krona is exactly the failure users should not have
        # to discover by squinting at a pie chart.
        pct = _report_percent_classified(report)
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
            for src in (report, output):
                dst = os.path.join(self.directory, os.path.basename(src))
                if os.path.exists(dst):
                    os.remove(dst)
                shutil.move(src, self.directory)
            self.report = f'{cwd}/{self.directory}/{sample_name}_reportkraken.txt'
            self.output = f'{cwd}/{self.directory}/{sample_name}_outputkraken.txt'
            log_file = open("kraken_log.txt", "a")
            try:
                log_file.write(f'DB used: {os.readlink(self.kraken_db)}')
            except OSError:
                log_file.write(f'DB used: {self.kraken_db}')
            log_file.close()
            dst_log = os.path.join(self.directory, "kraken_log.txt")
            if os.path.exists(dst_log):
                os.remove(dst_log)
            shutil.move("kraken_log.txt", self.directory)

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
            dst = os.path.join(self.directory, f'{self.sample_name}_{self.date_stamp}_krona.html')
            if os.path.exists(dst):
                os.remove(dst)
            shutil.move(f'{self.sample_name}_{self.date_stamp}_krona.html', self.directory)
            self.krona_html = f'{self.cwd}/{self.directory}/{self.sample_name}_{self.date_stamp}_krona.html'
            
    def bracken(self, report, output):
        rc = subprocess.call(f'bracken -d {self.kraken_db} -i {report} -o {self.sample_name}-bracken.txt -r 250', shell=True)
        if rc != 0 or not os.path.exists(f'{self.sample_name}-bracken.txt'):
            # Bracken refuses when the DB lacks its kmer distribution files or
            # when (near-)nothing was classified. It's an abundance refinement,
            # not the identification itself — skip it loudly rather than die.
            print(f'{bcolors.YELLOW}WARNING: bracken exited {rc} — skipping abundance '
                  f're-estimation (the Kraken report and Krona graph stand on their own).{bcolors.ENDC}')
            self.bracken_excel = None
            return
        df = pd.read_csv(f'{self.sample_name}-bracken.txt', sep='\t')
        df.to_excel(f'{self.sample_name}-bracken.xlsx', index=False)
        os.remove(f'{self.sample_name}-bracken.txt')
        self.bracken_excel = f'{os.getcwd()}/{self.sample_name}-bracken.xlsx'
        if self.directory:
            dst = os.path.join(self.directory, f'{self.sample_name}-bracken.xlsx')
            if os.path.exists(dst):
                os.remove(dst)
            shutil.move(f'{self.sample_name}-bracken.xlsx', self.directory)
            self.bracken_excel = f'{os.getcwd()}/{self.directory}/{self.sample_name}-bracken.xlsx'

class Bracken_Pie_Charts:

    def __init__(self, FASTA=False):
        self.FASTA = FASTA

    def run(self, bracken_excel,):
        df = pd.read_excel(bracken_excel)
        df = df[df['fraction_total_reads'] > 0.01 ]
        df2 = pd.DataFrame([['unclassified', 0, 'na', 0, 0, 0, 1 - df['fraction_total_reads'].sum()]], columns=['name','taxonomy_id','taxonomy_lvl', 'kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads'])
        df3 = pd.concat([df, df2])
        df3 = df3.set_index('name')
        if self.FASTA:
            plot = df3.plot.pie(y='fraction_total_reads', title='Identification of Assembled Scaffolds', figsize=(9, 5), cmap=cmap, labeldistance=None, legend=True, autopct='%1.1f%%')
        else:  #default FASTQ
            plot = df3.plot.pie(y='fraction_total_reads', title='FASTQ Read Identification', figsize=(9, 5), cmap=cmap, labeldistance=None, legend=True, autopct='%1.1f%%')
        plot.axis('off')
        plot.legend(bbox_to_anchor=(0.9, 0.9))
        plot.yaxis.label.set_visible(False)
        plot.get_figure().savefig(f'{os.getcwd()}/bracken_pie.png', format='png', bbox_inches='tight')
        self.pie_chart = f'{os.getcwd()}/bracken_pie.png'

    def latex(self, build_latex=False):
            tex = build_latex
            print(r'\begin{table}[H]', file=tex)
            print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
            if self.FASTA:
                bracken_pie_banner = Banner("FASTA Identifications")
            else:
                bracken_pie_banner = Banner("FASTQ Identifications")
            print(r'\begin{center}', file=tex)
            print(r'\includegraphics[scale=1]{' + bracken_pie_banner.banner + r'}', file=tex)
            print(r'\end{center}', file=tex)
            print(r'\end{adjustbox}', file=tex)
            print(r'\begin{center}', file=tex)
            print(r'\includegraphics[scale=0.8]{' + self.pie_chart + r'}', file=tex)
            print(r'\end{center}', file=tex)
            print(r'Identified using: \href{https://ccb.jhu.edu/software/kraken2/}{Kraken} and \href{https://ccb.jhu.edu/software/bracken/}{Bracken}', file=tex)
            print(r'\end{table}', file=tex)


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
