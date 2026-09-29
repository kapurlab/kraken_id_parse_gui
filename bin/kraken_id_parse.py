#!/usr/bin/env python


# --- provenance: log every external command this pipeline runs (best-effort) ---
# Captures commands launched directly by this orchestrator. Child executables
# remain responsible for their own nested-command provenance. Kraken's output
# contract is its working directory, so capture remains rooted at cwd here.
def _install_provenance_capture():
    import os as _o, subprocess as _s, shlex as _sh
    from pathlib import Path as _P
    from datetime import datetime as _dt
    _tool = _P(__file__).resolve().parents[1].name
    _out = _P.cwd() / ".provenance"
    _f = _out / (_tool + "_commands.txt")
    def _log(_cmd, _cwd=None):
        try:
            _out.mkdir(parents=True, exist_ok=True)
            _ln = _cmd if isinstance(_cmd, str) else _sh.join(str(c) for c in _cmd)
            _ts = _dt.now().astimezone().strftime("%H:%M:%S")
            with open(_f, "a", encoding="utf-8") as _h:
                _where = _P(_cwd).resolve() if _cwd else _P.cwd().resolve()
                _h.write(_ts + "  cwd=" + _sh.quote(str(_where)) + "  " + _ln + "\n")
        except Exception:
            pass
    try:
        _out.mkdir(parents=True, exist_ok=True)
        with open(_f, "a", encoding="utf-8") as _h:
            _h.write("\n# === %s run %s — external commands that produced results in this folder ===\n"
                     % (_tool, _dt.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")))
    except Exception:
        pass
    _orig_popen = _s.Popen
    class _Popen(_orig_popen):
        def __init__(self, args, *a, **k):
            _log(args, k.get("cwd"))
            super().__init__(args, *a, **k)
    _s.Popen = _Popen
    _osys = _o.system
    def _sysw(_cmd):
        _log(_cmd)
        return _osys(_cmd)
    _o.system = _sysw
try:
    _install_provenance_capture()
except Exception:
    pass
# --- end provenance ------------------------------------------------------------

__version__ = "0.0.1"

import os
import shutil
import sys
import glob
import re
import datetime
import random
import time
import logging
import subprocess
from pathlib import Path
import argparse
import textwrap
import pandas as pd
from collections import OrderedDict, defaultdict
from random import randint
from time import sleep
from Bio import SeqIO

from file_setup import bcolors, Excel_Stats, SummaryLog, UI
from report_html import HtmlReport
import coverage_plot

from fastq_stats_seqkit import FASTQ_Stats
from alignment_vcf import Alignment
from blast_fasta_and_search import Blast_Fasta, BLAST_TIMEOUT_SECONDS
# Same downloader the segmented engine uses: honors NCBI_EMAIL / NCBI_API_KEY,
# which the reference-download failure message tells users to set.
from download_fasta_by_acc import Downloader
from kraken import Kraken_Identification
from kraken import Bracken_Pie_Charts
from kraken import BrackenNoReadsError
from reference_guided_assembly_vcf_to_fasta import Reference_Guided_Assembly
from parse_reads import ParseReads, TaxonNotFoundError
from assembly import Assemble, SPAdesDidNotAssembleFASTA
# from blast_hpc import Blast
from coverage_graph import Coverage_Graph
from blast_to_coverage import BlastCoverageBridge
from orbivirus_specific import Orbivirus_Specific, resolve_segment_number
from orbivirus_specific import orbivirus_species as _orbivirus_species
from isav_specific import ISAV_Specific
from apicomplexa_specific import Apicomplexa
from organism_profiles import SegmentedEngine
from organism_registry import PROFILES, BESPOKE_TAXA
from coverage_graph_generator import CoverageGraphGenerator
from btv_serotyping import BTVSerotyping
from reporting import build_run_manifest, render_html_report, render_pdf_report, write_manifest


def _is_btv(description):
    """True if a BLAST/FASTA description names Bluetongue virus."""
    return _orbivirus_species(description) == 'BTV'


def _is_ehd(description):
    """True if a BLAST/FASTA description names Epizootic Hemorrhagic Disease virus."""
    return _orbivirus_species(description) == 'EHD'


def _segment_label(ref_id, fallback):
    """Excel label for a reference-guided sequence ID: 'seg2', or 'BTV seg2' /
    'EHD seg2' for the species-tagged IDs of a mixed BTV + EHD sample."""
    m = re.search(r'(?:(BTV|EHD)_)?segment(\d+)', str(ref_id), re.IGNORECASE)
    if not m:
        return fallback
    return f"{m.group(1) + ' ' if m.group(1) else ''}seg{m.group(2)}"


def tentative_segments_by_species(segment_mapping):
    """Group the segments assigned by the cautious low-confidence fallback by
    organism: {'BTV': {6}, 'EHD': {...}, 'other': {...}}. Keeping the species
    apart means a tentative BTV segment never badges the same-numbered EHD
    segment in a mixed sample."""
    groups = defaultdict(set)
    for info in (segment_mapping or {}).values():
        if info.get('confidence') == 'tentative':
            species = _orbivirus_species(info.get('description', '')) or 'other'
            groups[species].add(info['segment_number'])
    return groups


# ----------------------------------------------------------------------------
# Coverage-graph report helpers (module-level so they are unit-testable)
# ----------------------------------------------------------------------------

def get_segment_status(png_dict, descriptions_dict, total_expected=10, tentative_segments=None):
    """Determine which segments were found and return status string.

    Uses the canonical resolver (orbivirus_specific.resolve_segment_number)
    so VP2/VP5-labeled references are recognized as segments 2/6 even when
    the header has no explicit "segment N" text — keeping this banner
    consistent with serotyping and the coverage graphs.

    tentative_segments (optional set of ints) are segments that were assigned
    by the cautious low-confidence fallback; they are noted in the banner.
    """
    tentative_segments = tentative_segments or set()
    found_segments = set()
    for var_name in png_dict.keys():
        desc = descriptions_dict.get(var_name, var_name)
        seg_num = resolve_segment_number(desc)
        if seg_num is None:
            seg_num = resolve_segment_number(var_name)
        if seg_num is not None:
            found_segments.add(seg_num)
    expected = set(range(1, total_expected + 1))
    missing = expected - found_segments
    if not missing:
        status = f"All {total_expected} segments found"
    else:
        missing_sorted = sorted(missing)
        status = (f"Segment{'s' if len(missing_sorted) > 1 else ''} "
                  f"{', '.join(str(s) for s in missing_sorted)} not found")
    tent = sorted(s for s in tentative_segments if s in found_segments)
    if tent:
        plural = 's' if len(tent) > 1 else ''
        status += (f" — segment{plural} {', '.join(str(s) for s in tent)} "
                   f"tentative (low-confidence assignment)")
    return status


def _is_no_reference_placeholder(text):
    """True if a reference id/header is one of the synthetic placeholders the
    pipeline writes when no real reference could be obtained (all downloads
    failed, or none available): 'No_Reference_Available', 'No_Segments_Available',
    'No_Sequences_Available'. Matches the truncated id form too (the consensus id
    is cut to ~20 chars, e.g. '..._No_Reference_Availab'), so we match on the
    surviving prefix rather than the full token."""
    t = str(text).lower()
    return any(tok in t for tok in ('no_reference', 'no_segments', 'no_sequences'))


def add_virus_section(report, banner_title, png_generator, group_keys, descriptions_dict,
                      show_segment_status=True, tentative_segments=None):
    """Add a coverage section to the HTML/PDF report for the given group of
    reference sequences. Builds an interactive Plotly graph (HTML) and a static
    PNG (PDF; SNP-density track for large genomes) from the generator's data.

    group_keys are CoverageGraphGenerator png_file_paths keys (one per reference).
    show_segment_status adds the "All N segments found" banner (Orbivirus only).
    """
    if not group_keys:
        return
    tentative_segments = tentative_segments or set()
    align = getattr(png_generator, '_last_alignment_stats', {}) or {}
    snps_all = getattr(png_generator, '_last_snps', {}) or {}
    nocov_all = getattr(png_generator, '_last_no_coverage', {}) or {}

    # Map png variable-name keys back to reference ids via matching header.
    items = []
    for var_name in group_keys:
        desc = descriptions_dict.get(var_name, var_name)
        # find the alignment-stats entry whose header matches this description
        ref_id, stats = None, None
        for rid, st in align.items():
            if st.get('header') == desc or rid == desc:
                ref_id, stats = rid, st
                break
        if stats is None:
            continue
        snps = snps_all.get(ref_id, [])
        nocov = nocov_all.get(ref_id, [])
        div = coverage_plot.interactive_div(ref_id, stats, snps, nocov, include_plotlyjs=False)
        png = coverage_plot.static_png(ref_id, stats, snps, nocov,
                                       out_dir='coverage_graph_alignment')
        caption = f'Coverage analysis with SNPs for: {desc}'
        seg_num = resolve_segment_number(desc)
        if seg_num is None:
            seg_num = resolve_segment_number(var_name)
        if seg_num in tentative_segments:
            caption += '  (tentative — low-confidence segment assignment; verify manually)'
        items.append({'caption': caption, 'div': div, 'img': png})
    if not items:
        return
    status = (get_segment_status({k: 1 for k in group_keys}, descriptions_dict,
                                 tentative_segments=tentative_segments)
              if show_segment_status else None)
    report.add_coverage(banner_title, items, status=status)


def add_alignment_stats_table(report, stats_dict, tentative_by_species=None):
    """Add the 'Alignment Statistics' (Table 1) section from generator stats.
    Synthetic no-reference placeholders are excluded so the table never shows a
    bogus all-N reference row.

    References whose segment was assigned by the cautious low-confidence
    fallback get a '(tentative)' badge and trigger an explanatory footnote.
    tentative_by_species comes from tentative_segments_by_species(), so each
    row is checked against its own organism's tentative segments."""
    tentative_by_species = tentative_by_species or {}
    rows = []
    any_tentative = False
    for rid, st in (stats_dict or {}).items():
        if _is_no_reference_placeholder(rid) or _is_no_reference_placeholder(st.get('header', '')):
            continue
        label = st.get('header', rid)
        seg_num = resolve_segment_number(rid)
        if seg_num is None:
            seg_num = resolve_segment_number(st.get('header', ''))
        species = _orbivirus_species(st.get('header', rid)) or 'other'
        if seg_num is not None and seg_num in tentative_by_species.get(species, set()):
            label = f"{label} (tentative)"
            any_tentative = True
        rows.append((label, st['length'], st['mean_coverage'], st['percent_covered'],
                     st.get('snp_count')))
    if not rows:
        return
    report.add_alignment_stats(rows, tentative_note=any_tentative)


def render_kraken_overview_table(report, report_file, top_n=25, min_percent=0.0):
    """Add a 'Kraken Detailed Classification' table of the taxa actually detected
    (from the Kraken report) so the report shows what IS present when the target
    taxon was not found. Kraken columns: percent, clade_reads, taxon_reads, rank,
    taxid, name."""
    if not report_file or not os.path.exists(report_file):
        return
    rows = []
    try:
        with open(report_file) as fh:
            for order, line in enumerate(fh):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 6:
                    continue
                try:
                    percent = float(parts[0])
                    clade_reads = int(parts[1])
                except ValueError:
                    continue
                rows.append({'name': parts[5].strip(), 'level': parts[3],
                             'reads': clade_reads, 'percent': percent, 'order': order})
    except Exception:
        return
    if not rows:
        return
    rows = [r for r in rows if r['percent'] >= min_percent]
    rows.sort(key=lambda r: r['reads'], reverse=True)
    rows = sorted(rows[:top_n], key=lambda r: r['order'])
    disp = [[r['name'], r['level'], f"{r['reads']:,}", f"{r['percent']:.2f}"] for r in rows]
    report.add_table('Kraken Detailed Classification',
                     ['Taxonomic Name', 'Level', 'Reads', 'Percent'], disp, num_cols={2, 3},
                     compact=True)


# Below this many extracted target reads, an assembly/consensus cannot be built;
# it almost always means the sample doesn't match the chosen taxon/database.
MIN_TARGET_READS = 100


def _taxon_db_mismatch_hint(taxon, kraken_db, extracted=None):
    """A clear, actionable message for the common 'wrong --taxon/--kraken_db for
    this sample' mistake, shown when almost nothing matched the target."""
    db = os.path.basename(str(kraken_db).rstrip('/')) if kraken_db else 'the Kraken database'
    got = f'Only {extracted:,} reads matched' if extracted is not None else 'Almost no reads matched'
    return (f"{got} target taxon '{taxon}' in database '{db}'. "
            f"This usually means --taxon and/or --kraken_db do not match this "
            f"sample. Check that the target organism, its Kraken database, and "
            f"--blast_db are the right ones for these reads (see the taxa actually "
            f"detected below).")


def _write_reports(report, logger=None, summary_log=None, write_pdf=True):
    """Write the interactive HTML report and render the matching PDF via
    WeasyPrint (write_pdf=False, i.e. --no-pdf, skips the PDF). Returns
    (html_path, pdf_path)."""
    html_path = pdf_path = None
    if summary_log:
        summary_log.start_section('Report Generation (HTML + PDF)')
    try:
        html_path = report.write_html()
        if logger:
            logger.info(f"HTML report written: {html_path}")
    except Exception as e:
        if logger:
            logger.error(f"HTML report generation failed: {e}")
        if summary_log:
            summary_log.log_error(f"HTML report: {e}")
    if write_pdf:
        try:
            pdf_path = report.write_pdf()
            if logger:
                logger.info(f"PDF report written: {pdf_path}")
        except Exception as e:
            if logger:
                logger.error(f"PDF report generation failed: {e}")
            if summary_log:
                summary_log.log_error(f"PDF report: {e}")
    elif logger:
        logger.info("--no-pdf: skipping the PDF rendering of the HTML report")
    if summary_log:
        summary_log.end_section() if (html_path or pdf_path) else summary_log.end_section('FAILED')
    return html_path, pdf_path


def setup_logging(debug_mode=False, sample_name="kraken_analysis"):
    """
    Set up comprehensive logging system that writes to both console and log file.
    """
    # Create formatters
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s')

    # Create logger
    logger = logging.getLogger('kraken_pipeline')
    logger.setLevel(logging.DEBUG if debug_mode else logging.INFO)

    # Clear any existing handlers
    logger.handlers.clear()

    # Console handler - always show INFO and above
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    if debug_mode:
        # File handler - show DEBUG and above when in debug mode
        log_filename = f"{sample_name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}_debug.log"
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

        # Also show DEBUG on console when in debug mode
        console_handler.setLevel(logging.DEBUG)

        logger.info(f"Debug mode enabled - extensive logging to {log_filename}")
        logger.debug(f"Logger initialized with {len(logger.handlers)} handlers")

    return logger


def debug_command(cmd, logger, description="Running command"):
    """
    Execute a command with extensive logging when in debug mode.
    """
    logger.debug(f"{description}: {cmd}")

    try:
        # Use subprocess for better control and logging
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=3600)

        if result.stdout:
            logger.debug(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"STDERR:\n{result.stderr}")

        logger.debug(f"Command exit code: {result.returncode}")

        if result.returncode != 0:
            logger.error(f"Command failed with exit code {result.returncode}")
            logger.error(f"Failed command: {cmd}")
            if result.stderr:
                logger.error(f"Error output: {result.stderr}")
        else:
            logger.debug(f"Command completed successfully")

        return result.returncode == 0, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after 1 hour: {cmd}")
        return False, "", "Command timed out"
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return False, "", str(e)


if __name__ == "__main__": # execute if directly access by the interpreter
    start_time = datetime.datetime.now()
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    cp /project/bioinformatic_databases/bluetongue_virus/ehd_refseq.fasta .

    ehd_report.py -r1 *fastq.gz -f *fasta
    ehd_report.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -f *fasta
    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=False, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-l', '--logo', action='store', dest='logo', required=False, help='Logo image (.png or .svg) for the report header. For USDA reports use the official USDA signature lockup file from the USDA Style Guide (Office of Communications); it is shown unmodified at the top left.')
    parser.add_argument("-t", "--taxon", action='store', dest='taxon', help='Target Taxon')
    parser.add_argument("-k", "--kraken_db", action='store', dest='kraken_db', help='Specify Kraken db to use')
    parser.add_argument("-b", "--blast_db", action='store', dest='blast_db', default="nt", help='Specify BLAST db to use')
    parser.add_argument("-s", "--specific", action='store', dest='specific', default=None, help='Specify custom script/function for the target being used.  Often just default to the taxon name')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Enable extensive debugging with log file creation and preserve temporary files')
    parser.add_argument('--fast', action='store_true', dest='fast_mode', default=False, help='Enable fast assembly mode - uses fewer CPU threads and less RAM to reduce system lag')
    parser.add_argument('--threads', action='store', dest='threads', type=int, default=None, help='Worker threads for kraken2/BLAST/SPAdes/bwa. Default: auto-detect (leaves 2 cores free). Set this (or KIP_THREADS / sbatch --cpus-per-task) to pack several samples onto one node.')
    parser.add_argument('--kraken-only', action='store_true', dest='kraken_only', default=False, help='Run only Kraken2 and produce the Krona graph; skip read parsing, assembly, and BLAST. Requires -k.')
    parser.add_argument('--no-blast', action='store_true', dest='no_blast', default=False, help='Run Kraken2 and taxonomic read parsing, then stop: skip assembly, BLAST, and coverage. Leaves the parsed FASTQ.gz reads for the target taxon. Requires -k and -t.')
    parser.add_argument('--no-pdf', action='store_true', dest='no_pdf', default=False, help='Skip the PDF renderings of the reports. The HTML reports (and run_manifest.json) are always written; each PDF is just WeasyPrint\'s print rendering of its HTML and is not always needed.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    # Publish an explicit thread count to every downstream module (they read it
    # via file_setup.resolve_cpus()). Setting the env var here means a single
    # --threads N (or --override threads=N via run_with_config.py) is honored by
    # the whole pipeline.
    if args.threads and args.threads > 0:
        os.environ['KIP_THREADS'] = str(args.threads)

    # Pre-flight: verify all required external tools are on PATH
    REQUIRED_TOOLS = [
        'kraken2', 'kreport2krona.py', 'ktImportText',
        'spades.py', 'blastn', 'bwa', 'samtools', 'picard',
        'freebayes', 'freebayes-parallel', 'vcffilter', 'pigz', 'seqkit',
    ]
    print(f'\n{"="*55}')
    print(f'  PRE-FLIGHT TOOL CHECK')
    print(f'{"="*55}')
    missing_tools = []
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool)
        if path:
            print(f'  OK      {tool}: {path}')
        else:
            print(f'  MISSING {tool}  <-- not found on PATH')
            missing_tools.append(tool)
    # Bracken is optional (no macOS/arm64 build): found on PATH, in a sibling
    # 'bracken' env, or via $BRACKEN_BIN — otherwise its pie chart is skipped.
    if Kraken_Identification.bracken_available():
        print(f'  OK      bracken (optional)')
    else:
        print(f'  SKIP    bracken (optional) — not found; the Bracken abundance chart will be skipped')
    if missing_tools:
        print(f'\n  WARNING: {len(missing_tools)} tool(s) missing: {", ".join(missing_tools)}')
        print(f'  Pipeline will fail when these are reached.\n')
    else:
        print(f'  All tools found.\n')
    print(f'{"="*55}\n')

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    # Validate inputs up front so a missing/mistyped FASTQ path fails with a
    # clear message instead of a confusing traceback deeper in the pipeline.
    input_errors = []
    if not args.FASTQ_R1:
        input_errors.append('No R1 FASTQ provided (use -r1/--read1).')
    elif not os.path.isfile(args.FASTQ_R1):
        input_errors.append(f'R1 FASTQ not found: {args.FASTQ_R1}')
    if args.FASTQ_R2 and not os.path.isfile(args.FASTQ_R2):
        input_errors.append(f'R2 FASTQ not found: {args.FASTQ_R2}')
    if input_errors:
        for msg in input_errors:
            print(f'ERROR: {msg}', file=sys.stderr)
        print('Run from the sample\'s working directory and pass the actual FASTQ '
              'filenames (e.g. -r1 mysample_R1.fastq.gz -r2 mysample_R2.fastq.gz).',
              file=sys.stderr)
        sys.exit(2)

    # The two short run modes need their inputs up front, not after minutes of
    # FASTQ statistics.
    if args.kraken_only and not args.kraken_db:
        print(f"{bcolors.RED}ERROR: --kraken-only requires a Kraken DB (-k/--kraken_db).{bcolors.ENDC}")
        sys.exit(2)
    if args.no_blast:
        if not args.kraken_db:
            print(f"{bcolors.RED}ERROR: --no-blast requires a Kraken DB (-k/--kraken_db).{bcolors.ENDC}")
            sys.exit(2)
        if not args.taxon:
            print(f"{bcolors.RED}ERROR: --no-blast requires a target taxon (-t/--taxon).{bcolors.ENDC}")
            sys.exit(2)

    # Extract sample name early for logger setup
    sample_name = re.sub('[._].*', '', os.path.basename(args.FASTQ_R1))

    # Set up comprehensive logging system
    logger = setup_logging(args.debug, sample_name)
    logger.info(f"Starting {os.path.basename(__file__)} version {__version__}")
    logger.info(f"Sample name: {sample_name}")
    if args.debug:
        logger.debug(f"Debug mode enabled - preserving temporary files and extensive logging")

    if args.fast_mode:
        logger.info("Fast assembly mode enabled - using reduced CPU/RAM to minimize system lag")

    logger.debug(f"Arguments: {vars(args)}")

    if args.debug:
        logger.debug("Environment variables:")
        for key, value in os.environ.items():
            if any(term in key.upper() for term in ['PATH', 'PYTHON', 'CONDA', 'BLAST', 'KRAKEN']):
                logger.debug(f"  {key}={value}")

    logger.debug(f"Working directory: {os.getcwd()}")
    logger.debug(f"Current user: {os.getenv('USER', 'unknown')}")
    logger.debug(f"Available disk space: {shutil.disk_usage('.').free / (1024**3):.2f} GB")

    # Initialize reporting systems
    logger.info("Initializing report generation systems...")
    logger.debug(f"Creating HTML report with logo: {args.logo}")
    # Report colors are the official USDA palette (USDA visual standards), not
    # derived from the logo image.
    report = HtmlReport(sample_name=sample_name, logo=args.logo, out_dir=os.getcwd())

    logger.debug("Creating Excel stats report")
    excel_stats = Excel_Stats(sample_name)
    logger.debug(f"Excel report file: {excel_stats.excel_filename}")

    summary_log = SummaryLog(sample_name=sample_name, script_version=__version__, args=args)
    summary_log.collect_versions()
    logger.info(f"Git branch: {summary_log.git_branch} (commit {summary_log.git_commit})")
    logger.debug(f"Summary log file: {summary_log.log_file}")

    logger.info("Starting FASTQ quality statistics analysis...")
    logger.debug(f"FASTQ R1: {args.FASTQ_R1}")
    logger.debug(f"FASTQ R2: {args.FASTQ_R2}")

    if args.FASTQ_R1 and os.path.exists(args.FASTQ_R1):
        file_size_r1 = os.path.getsize(args.FASTQ_R1) / (1024**2)
        logger.debug(f"R1 file size: {file_size_r1:.2f} MB")
    else:
        logger.error(f"R1 file not found or not provided: {args.FASTQ_R1}")

    if args.FASTQ_R2 and os.path.exists(args.FASTQ_R2):
        file_size_r2 = os.path.getsize(args.FASTQ_R2) / (1024**2)
        logger.debug(f"R2 file size: {file_size_r2:.2f} MB")

    summary_log.start_section('FASTQ Quality Statistics')
    try:
        fastq_stats = FASTQ_Stats(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=args.debug)
        logger.debug("FASTQ_Stats object created successfully")

        logger.debug("Running FASTQ statistics analysis...")
        fastq_stats.run()
        logger.info("FASTQ statistics completed successfully")

        logger.debug("Adding FASTQ stats to report...")
        report.add_fastq_quality(fastq_stats.R1, fastq_stats.R2)

        logger.debug("Adding FASTQ stats to Excel report...")
        fastq_stats.excel(excel_stats.excel_dict)
        summary_log.end_section()

    except Exception as e:
        summary_log.log_error(f"FASTQ statistics: {e}")
        summary_log.end_section('FAILED')
        logger.error(f"Error in FASTQ statistics: {e}")
        raise

    # Glanceable run overview strip at the top of the report (always available,
    # even if a later step fails, since it only depends on args + FASTQ stats).
    try:
        total_reads = int(str(fastq_stats.R1.num_seqs).replace(',', ''))
        if fastq_stats.R2:
            total_reads += int(str(fastq_stats.R2.num_seqs).replace(',', ''))
        reads_display = f'{total_reads:,}'
    except Exception:
        reads_display = 'N/A'
    report.add_overview([
        ('Target Taxon', args.taxon or 'N/A'),
        ('Total Reads', reads_display),
        ('Kraken Database', os.path.basename(str(args.kraken_db).rstrip('/')) if args.kraken_db else 'N/A'),
        ('Run Date', datetime.datetime.now().strftime('%B %d, %Y')),
    ])

    if args.kraken_db:
        summary_log.start_section('Kraken Classification & Bracken')
        logger.info("Running Kraken taxonomic classification...")
        logger.debug(f"Kraken database: {args.kraken_db}")

        # Verify Kraken database exists
        if os.path.exists(args.kraken_db):
            db_size = sum(os.path.getsize(os.path.join(args.kraken_db, f))
                         for f in os.listdir(args.kraken_db) if os.path.isfile(os.path.join(args.kraken_db, f)))
            logger.debug(f"Kraken database size: {db_size / (1024**3):.2f} GB")
        else:
            logger.error(f"Kraken database not found: {args.kraken_db}")

        try:
            kraken = Kraken_Identification(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, kraken_db=args.kraken_db, directory='kraken')
            logger.debug("Kraken_Identification object created successfully")

            logger.debug("Running Kraken classification...")
            kraken.run()
            logger.info("Kraken classification completed")

            logger.debug("Creating Krona visualization...")
            krona_html = kraken.krona_make_graph(kraken.report)
            logger.debug(f"Krona HTML created: {krona_html}")

            if args.kraken_only:
                # Kraken-only mode: the Krona graph is the deliverable. Skip Bracken,
                # read parsing, assembly, and BLAST entirely and finish here.
                print(f"\n{bcolors.GREEN}Kraken-only mode: Krona graph generated at {krona_html}.{bcolors.ENDC}")
                print(f"{bcolors.GREEN}Skipping Bracken, read parsing, assembly, and BLAST.{bcolors.ENDC}")
                summary_log.end_section()
                summary_log.write()
                sys.exit(0)

            logger.debug("Running Bracken abundance estimation...")
            # Bracken re-estimates abundances but has no macOS/arm64 conda build,
            # so it may be absent. It's optional — when it can't be found, warn
            # and skip Bracken + its pie chart so read parsing/identification
            # still run.
            if not Kraken_Identification.bracken_available():
                logger.warning("bracken not found ($BRACKEN_BIN, a sibling 'bracken' env, or PATH) "
                               "— skipping Bracken abundance re-estimation and its pie chart. "
                               "Read parsing and identification continue.")
                summary_log.log_error("Bracken: not installed — abundance chart skipped (continuing)")
            else:
                # A near-empty classification (sample doesn't match the DB) is a
                # legitimate 'nothing found' outcome, not a crash: skip the pie chart
                # and let the pipeline flow to the graceful 'taxon not found' report.
                try:
                    kraken.bracken(kraken.report, kraken.output)
                    logger.debug("Creating Bracken pie charts...")
                    bracken_pie_charts = Bracken_Pie_Charts()
                    bracken_pie_charts.run(kraken.bracken_excel)
                    report.add_pie(getattr(bracken_pie_charts, 'pie_chart', None))
                    logger.info("Bracken analysis and visualization completed")
                except BrackenNoReadsError as e:
                    logger.warning(f"Bracken abundance estimation skipped: {e}")
                    classified = getattr(kraken, 'percent_classified', None)
                    if classified is not None and classified >= 5.0:
                        # Plenty was classified, so this is not a sample/database
                        # mismatch: Bracken itself could not produce an estimate.
                        summary_log.log_error("Bracken: no output (continuing)")
                        report.add_message(
                            'Bracken abundance not estimated',
                            f'Bracken exited without an abundance estimate, so no abundance '
                            f'chart was produced. Kraken classified {classified:.1f}% of the '
                            f'reads, so this is not a database mismatch; the usual cause is a '
                            f'Kraken database without Bracken\'s kmer distribution file for '
                            f'250 bp reads. The Kraken classification below is unaffected.',
                            level='warn')
                    else:
                        summary_log.log_error("Bracken: no species-level reads (continuing)")
                        report.add_message(
                            'Very low classification',
                            'Almost no reads were classified against the Kraken database, so '
                            'no abundance chart was produced. If the target taxon is not '
                            'found below, the sample likely does not match this database '
                            '(check the --taxon and --kraken_db choices).', level='warn')
            summary_log.end_section()

        except Exception as e:
            summary_log.log_error(f"Kraken/Bracken: {e}")
            summary_log.end_section('FAILED')
            logger.error(f"Error in Kraken/Bracken analysis: {e}")
            raise
    else:
        logger.info("No Kraken database specified - looking for existing results...")

        class Kraken:
            def __init__(self):
                self.output = None
                self.report = None

        # Create kraken instance
        kraken = Kraken()

        # Check if either kraken or kraken2 folders exist
        kraken_folder = None
        logger.debug("Checking for existing Kraken output directories...")

        if os.path.isdir("kraken"):
            kraken_folder = "kraken"
            logger.debug("Found 'kraken' directory")
        elif os.path.isdir("kraken2"):
            kraken_folder = "kraken2"
            logger.debug("Found 'kraken2' directory")
        else:
            logger.debug("No kraken or kraken2 directory found")

        if kraken_folder:
            logger.debug(f"Searching for Kraken files in {kraken_folder}/")

            # Find the files using glob
            output_files = glob.glob(os.path.join(kraken_folder, "*_outputkraken.txt"))
            report_files = glob.glob(os.path.join(kraken_folder, "*_reportkraken.txt"))

            logger.debug(f"Found {len(output_files)} output files: {output_files}")
            logger.debug(f"Found {len(report_files)} report files: {report_files}")

            # Check if files were found and assign to kraken object
            if output_files and report_files:
                kraken.output = output_files[0]  # Get the first matching output file
                kraken.report = report_files[0]  # Get the first matching report file

                # Check file sizes
                output_size = os.path.getsize(kraken.output) / (1024**2)
                report_size = os.path.getsize(kraken.report) / (1024**2)
                logger.info(f"Using existing Kraken results:")
                logger.info(f"  Output: {kraken.output} ({output_size:.2f} MB)")
                logger.info(f"  Report: {kraken.report} ({report_size:.2f} MB)")
            else:
                logger.error("Required Kraken files not found in the kraken folder")
                logger.debug("Expected files: *_outputkraken.txt and *_reportkraken.txt")
        else:
            logger.warning("Neither kraken nor kraken2 folder found in current directory.")
            logger.info("Run script with -k option to first run Kraken.")

############################
    def cleanup_artifacts(debug: bool = False, logger=None):
        """Clean up temporary artifacts - preserve everything in debug mode"""
        if logger:
            logger.debug(f"Cleanup artifacts called with debug={debug}")

        if debug:
            if logger:
                logger.info("Debug mode: Preserving all temporary files for troubleshooting")
                logger.debug("Files preserved include:")
                logger.debug("  BLAST output files")
                logger.debug("  Coverage files")
                logger.debug("  Alignment directories")
                logger.debug("  PNG images and other temporary files")
            return  # Exit early - don't clean anything in debug mode

        if logger:
            logger.info("Cleaning up temporary artifacts...")

        temp_dir = './temp'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        files_grab = []
        # NOTE: do not add '*.log' here — it would delete the '*_summary.log'
        # deliverable (and any user run.log). The old '*.aux'/'*.log'/'*tex'/
        # '*out' patterns were pdflatex intermediates and are gone with LaTeX.
        # The BLAST summary, coverage-graph PDF/PNGs, the Bracken pie,
        # combined_references.fasta and CAUTION_SITES.xlsx are kept: the
        # structured report (run_manifest.json + report.html) is built from them
        # after this cleanup, and the GUI's Results pane serves them.
        file_patterns = ('batch.sh', '*_blast_all.txt', '*_blast_out.txt',
                         'coverage_list.txt', 'coverage.txt')

        for pattern in file_patterns:
            found_files = glob.glob(pattern)
            files_grab.extend(found_files)
            if logger and found_files:
                logger.debug(f"Found {len(found_files)} files matching {pattern}: {found_files}")

        if logger:
            logger.debug(f"Total files to clean up: {len(files_grab)}")

        if files_grab:
            print("Removing:")
            for each in files_grab:
                print(f'\t{each}')
                if logger:
                    logger.debug(f"Moving to temp: {each}")

                dest_path = os.path.join(temp_dir, os.path.basename(each))
                # Remove existing file if it exists to prevent shutil.Error
                if os.path.exists(dest_path):
                    os.remove(dest_path)
                    if logger:
                        logger.debug(f"Removed existing file: {dest_path}")

                shutil.move(each, temp_dir)

        shutil.rmtree(temp_dir)
        if logger:
            logger.debug("Temporary directory removed")

        # Clean up alignment directories
        for dir_name in ['alignment_all', 'alignment_top']:
            try:
                if os.path.exists(dir_name):
                    shutil.rmtree(dir_name)
                    if logger:
                        logger.debug(f"Removed directory: {dir_name}")
            except Exception as e:
                if logger:
                    logger.debug(f"Could not remove {dir_name}: {e}")

    def write_structured_reports(status: str, warnings=None):
        warnings = warnings or []
        try:
            output_dir = Path.cwd()
            manifest = build_run_manifest(
                sample_id=sample_name,
                status=status,
                warnings=warnings,
                inputs={"r1": args.FASTQ_R1, "r2": args.FASTQ_R2},
                parameters={
                    "taxon": args.taxon,
                    "kraken_db": args.kraken_db,
                    "blast_db": args.blast_db,
                },
                output_dir=output_dir,
                started_at=start_time,
            )
            manifest_path = write_manifest(manifest, output_dir)
            html_path = render_html_report(manifest, output_dir)
            # The PDF is nothing more than WeasyPrint's print rendering of
            # report.html — skippable because it is not always wanted.
            if args.no_pdf:
                print("--no-pdf: skipping the PDF rendering of report.html.")
            else:
                pdf_path, pdf_warning = render_pdf_report(html_path, output_dir)
                if pdf_warning:
                    manifest.setdefault("warnings", []).append(pdf_warning)
                    print(f"WARNING: {pdf_warning}")
                else:
                    print(f"PDF report generated: {pdf_path}")
                # Re-scan so the manifest's artifact list includes (or the
                # warning mentions) the PDF outcome, then re-render the HTML
                # from that final manifest.
                manifest = build_run_manifest(
                    sample_id=sample_name,
                    status=status,
                    warnings=manifest.get("warnings", warnings),
                    inputs={"r1": args.FASTQ_R1, "r2": args.FASTQ_R2},
                    parameters={
                        "taxon": args.taxon,
                        "kraken_db": args.kraken_db,
                        "blast_db": args.blast_db,
                    },
                    output_dir=output_dir,
                    started_at=start_time,
                )
                manifest_path = write_manifest(manifest, output_dir)
                html_path = render_html_report(manifest, output_dir)
            print(f"Structured report manifest generated: {manifest_path}")
            print(f"HTML report generated: {html_path}")
        except Exception as exc:
            print(f"WARNING: structured HTML report generation failed: {exc}")

    """Execute the complete analysis pipeline"""
    logger.info(UI.banner("Starting bioinformatics analysis pipeline"))
    logger.debug("Pipeline overview:")
    logger.debug("  1. Taxonomic read parsing")
    logger.debug("  2. Assembly of filtered reads")
    logger.debug("  3. BLAST analysis")
    logger.debug("  4. Coverage analysis")
    logger.debug("  5. Organism-specific analysis (if applicable)")
    logger.debug("  6. Reference-guided assembly")
    logger.debug("  7. Report generation")

    # Initialize global segment mapping dictionary to preserve segment assignments throughout pipeline
    global_segment_mapping = {}
    # Number of segments/chromosomes the mapping is checked against for the
    # Excel 'Segments Missing' summary (Orbivirus = 10; the generic segmented
    # handler uses its organism profile's count, e.g. ISAV = 8).
    segment_total_expected = 10
    logger.debug("Initialized global segment mapping dictionary")

    # Initialize Orbivirus dictionaries to preserve throughout pipeline
    btv_dict = {}
    ehv_dict = {}
    logger.debug("Initialized Orbivirus dictionaries (btv_dict, ehv_dict)")

    # Step 1: Parse reads using Kraken results
    summary_log.start_section('Taxonomic Read Parsing')
    logger.info(UI.step(1, 4, "Running taxonomic read parsing..."))
    logger.debug(f"Target taxon: {args.taxon}")
    logger.debug(f"Kraken output file: {kraken.output}")
    logger.debug(f"Kraken report file: {kraken.report}")

    if kraken.output and os.path.exists(kraken.output):
        output_lines = sum(1 for _ in open(kraken.output, 'r'))
        logger.debug(f"Kraken output contains {output_lines} classified reads")

    try:
        parser = ParseReads(
            R1=args.FASTQ_R1,
            R2=args.FASTQ_R2,
            kraken_output=kraken.output,
            kraken_report=kraken.report,
            taxon=args.taxon,
            output_prefix=args.taxon,
            debug=args.debug
        )
        logger.debug("ParseReads object created successfully")

        logger.debug("Running taxonomic read parsing...")
        parser.run()
        logger.info("Taxonomic read parsing completed successfully")

        logger.debug("Adding parsing results to report...")
        report.add_kv('Kraken Classification', [
            ('Target Taxon', parser.metrics['target_taxon']),
            ('Taxon ID', parser.metrics['taxid']),
            ('Total Input Reads', f"{parser.metrics['total_input_reads']:,}"),
            ('Extracted Reads', f"{parser.metrics['extracted_reads']:,}"),
            ('Extraction Rate', f"{parser.metrics['extraction_rate']}%"),
        ])
        _kd = [[c['name'], c['level'], f"{c['reads']:,}", f"{c['percent']:.2f}"]
               for c in parser.classifications]
        report.add_table('Kraken Detailed Classification',
                         ['Taxonomic Name', 'Level', 'Reads', 'Percent'], _kd, num_cols={2, 3},
                         compact=True)

        # Catch a taxon/database mismatch: if almost nothing matched the target,
        # an assembly/consensus can't be built. Stop cleanly with clear, colored
        # guidance and the taxa actually detected, instead of failing cryptically
        # deep in assembly/coverage.
        extracted = int(parser.metrics.get('extracted_reads') or 0)
        if extracted < MIN_TARGET_READS:
            hint = _taxon_db_mismatch_hint(args.taxon, args.kraken_db, extracted)
            print(UI.error('Possible taxon / database mismatch'))
            print(UI.warn(hint))
            logger.error(f"Only {extracted} reads matched '{args.taxon}' "
                         f"(< {MIN_TARGET_READS}) — likely taxon/database mismatch")
            summary_log.log_error(
                f"Only {extracted} reads matched '{args.taxon}' — likely taxon/database mismatch")
            summary_log.end_section('FAILED — likely taxon/database mismatch')
            report.add_message('Possible taxon / database mismatch', hint, level='error')
            render_kraken_overview_table(report, kraken.report)
            parser.excel(excel_stats.excel_dict)
            _write_reports(report, logger, summary_log, write_pdf=not args.no_pdf)
            excel_stats.post_excel()
            summary_log.write()
            cleanup_artifacts(debug=args.debug, logger=logger)
            write_structured_reports("completed_with_warnings", [hint])
            sys.exit(0)

        logger.debug("Adding parsing results to Excel report...")
        parser.excel(excel_stats.excel_dict)
        summary_log.end_section()

    except TaxonNotFoundError as e:
        summary_log.log_error(f"Target taxon '{args.taxon}' not found by Kraken")
        summary_log.end_section('FAILED — taxon not found')
        logger.error(f"Target taxon '{args.taxon}' not found by Kraken")
        logger.info("Terminating target-specific analysis - showing taxa detected in the sample")

        report.add_message(
            f'Target {args.taxon} taxon not found by Kraken',
            'The target taxon was not detected. The most abundant taxa identified in '
            'this sample are shown below; no target-specific assembly or coverage '
            'analysis was performed.', level='error')
        # Still show what IS present in the sample
        render_kraken_overview_table(report, kraken.report)
        _write_reports(report, logger, summary_log, write_pdf=not args.no_pdf)
        excel_stats.excel_dict['Target Taxon'] = args.taxon
        excel_stats.excel_dict['Extracted Reads'] = "No Reads Found"
        excel_stats.post_excel()
        summary_log.write()
        cleanup_artifacts(debug=args.debug, logger=logger)
        write_structured_reports(
            "completed_with_warnings",
            [f"Target taxon not found by Kraken: {args.taxon}"],
        )
        sys.exit(0)
    except Exception as e:
        summary_log.log_error(f"Taxonomic read parsing: {e}")
        summary_log.end_section('FAILED')
        logger.error(f"Error in taxonomic read parsing: {e}")
        raise

    # --no-blast: the taxonomically-parsed reads are the deliverable. Stop here
    # before assembly/BLAST/coverage. The gzipped reads for the target taxon
    # (parser.r1_out.gz / parser.r2_out.gz) are left in place for reuse (e.g.
    # re-running them through vSNP).
    if args.no_blast:
        _parsed = ', '.join(f'{r}.gz' for r in (parser.r1_out, parser.r2_out) if r)
        print(f"\n{bcolors.GREEN}--no-blast mode: extracted target reads ({_parsed}).{bcolors.ENDC}")
        print(f"{bcolors.GREEN}Skipping assembly, BLAST, and coverage analysis.{bcolors.ENDC}")
        _write_reports(report, logger, summary_log, write_pdf=not args.no_pdf)
        excel_stats.post_excel()
        summary_log.write()
        cleanup_artifacts(debug=args.debug, logger=logger)
        write_structured_reports("completed")
        print(f"{bcolors.WHITE}Total runtime: {datetime.datetime.now() - start_time}{bcolors.ENDC}\n")
        sys.exit(0)

    # Step 2: Assembly of filtered reads
    summary_log.start_section('De Novo Assembly (SPAdes)')
    logger.info(UI.step(2, 4, "Running assembly of filtered reads..."))

    parsed_r1 = f"{parser.r1_out}.gz"
    # Single-end support (from master): no R2 when the parser produced none.
    parsed_r2 = f"{parser.r2_out}.gz" if parser.r2_out else None

    logger.debug(f"Assembly input files:")
    logger.debug(f"  R1: {parsed_r1}")
    logger.debug(f"  R2: {parsed_r2}")

    # Check parsed read file sizes
    if os.path.exists(parsed_r1):
        r1_size = os.path.getsize(parsed_r1) / (1024**2)
        logger.debug(f"  R1 size: {r1_size:.2f} MB")
    else:
        logger.error(f"Parsed R1 file not found: {parsed_r1}")

    if parsed_r2 and os.path.exists(parsed_r2):
        r2_size = os.path.getsize(parsed_r2) / (1024**2)
        logger.debug(f"  R2 size: {r2_size:.2f} MB")
    else:
        logger.debug(f"Parsed R2 file not found (single-end run): {parsed_r2}")

    try:
        assembler = Assemble(
            FASTQ_R1=parsed_r1,
            FASTQ_R2=parsed_r2,
            debug=args.debug,
            fast_mode=args.fast_mode
        )
        logger.debug("Assemble object created successfully")

        logger.debug("Starting SPAdes assembly...")
        assembler.run()
        logger.info("Assembly completed successfully")

        logger.debug("Generating assembly statistics...")
        assembler.stats(assembler.FASTA)

        if os.path.exists(assembler.FASTA):
            fasta_size = os.path.getsize(assembler.FASTA)
            logger.debug(f"Assembly output file: {assembler.FASTA} ({fasta_size} bytes)")
        summary_log.end_section()

    except SPAdesDidNotAssembleFASTA as e:
        summary_log.log_error(f"SPAdes assembly failed — not enough {args.taxon} reads")
        summary_log.end_section('FAILED')
        logger.error("SPAdes assembly failed")
        logger.error(f"Likely not enough {args.taxon} reads to complete assembly")
        logger.debug(f"Assembly error details: {e}")

        print('Assembly failed, Check reads')
        report.add_message(
            'Assembly failed',
            f'Assembly failed — check reads. Likely not enough {args.taxon} reads to '
            'complete assembly. Script terminated.', level='error')
        # Show what was detected even though assembly could not proceed.
        if 'kraken' in dir() and getattr(kraken, 'report', None):
            render_kraken_overview_table(report, kraken.report)
        _write_reports(report, logger, summary_log, write_pdf=not args.no_pdf)
        excel_stats.excel_dict['Contig count'] = "Assembly did not complete"
        excel_stats.post_excel()
        summary_log.write()
        cleanup_artifacts(debug=args.debug, logger=logger)
        write_structured_reports(
            "completed_with_warnings",
            [f"Assembly did not complete for target taxon: {args.taxon}"],
        )
        sys.exit(0)
    except Exception as e:
        summary_log.log_error(f"Assembly: {e}")
        summary_log.end_section('FAILED')
        logger.error(f"Error in assembly step: {e}")
        raise
    
    # Step 3: BLAST analysis
    summary_log.start_section('BLAST Analysis')
    logger.info(UI.step(3, 4, f"Running BLAST analysis using {args.blast_db}..."))
    logger.debug(f"BLAST database: {args.blast_db}")
    logger.debug(f"Assembly FASTA for BLAST: {assembler.FASTA}")

    try:
        # Segmented organisms (Orbivirus, ISAV, ...): SPAdes can fuse segments that
        # share conserved termini into one contig, hiding a segment behind the
        # contig's top hit, so search the unexplained stretch of such contigs too.
        blast = Blast_Fasta(FASTA=assembler.FASTA, blast_db=args.blast_db,
                            split_chimeras=args.taxon in PROFILES)
        logger.debug("Blast_Fasta object created successfully")

        logger.info("BLAST analysis completed")
        logger.debug(f"BLAST summary file: {blast.blast_summary_file}")

        if os.path.exists(blast.blast_summary_file):
            summary_size = os.path.getsize(blast.blast_summary_file)
            logger.debug(f"BLAST summary file size: {summary_size} bytes")

        logger.debug("Processing BLAST results for coverage analysis...")
        blast_to_coverage = BlastCoverageBridge(blast.blast_summary_file)
        blast_to_coverage.process()
        logger.debug(f"Combined FASTA for coverage: {blast_to_coverage.combined_fasta}")
        summary_log.end_section()

    except Exception as e:
        summary_log.log_error(f"BLAST analysis: {e}")
        summary_log.end_section('FAILED')
        logger.error(f"Error in BLAST analysis: {e}")
        raise

    # Step 4: Coverage analysis
    summary_log.start_section('Initial Coverage Analysis')
    logger.info(UI.step(4, 4, "Running coverage analysis..."))
    logger.debug("Getting top 5 unique accessions from BLAST results")

    try:
        # Get top 5 unique accessions from BLAST results
        coverage_graph = Coverage_Graph(FASTA=blast_to_coverage.combined_fasta, FASTQ_R1=assembler.FASTQ_R1, FASTQ_R2=assembler.FASTQ_R2, debug=args.debug)
        logger.debug("Coverage_Graph object created successfully")

        logger.debug("Generating coverage graphs...")
        coverage_graph.get_coverage_graph()
        logger.info("Coverage analysis completed successfully")

        if hasattr(coverage_graph, 'alignment_stats') and coverage_graph.alignment_stats:
            logger.debug(f"Coverage stats for {len(coverage_graph.alignment_stats)} sequences")
            for seq_id, stats in coverage_graph.alignment_stats.items():
                logger.debug(f"  {seq_id}: {stats.get('header', 'No header')}")
        summary_log.end_section()

    except Exception as e:
        summary_log.log_error(f"Coverage analysis: {e}")
        summary_log.end_section('FAILED')
        logger.error(f"Error in coverage analysis: {e}")
        raise

    # Update the assembled FASTA file with file name containing "denovo"
    path, filename = os.path.split(assembler.FASTA)
    name, ext = os.path.splitext(filename)
    new_filename = f"{name}_denovo{ext}"
    new_FASTA = os.path.join(path, new_filename)
    os.rename(assembler.FASTA, new_FASTA)
    assembler.FASTA = new_FASTA
    print(f"Updated FASTA variable: {assembler.FASTA}")

    report.add_assembly(assembler)
    _blast_base = os.path.basename(blast.blast_db)
    _blast_rows = [[r[0], r[1], (r[2] if not isinstance(r[2], (list, tuple)) else ' '.join(map(str, r[2])))]
                   for r in (blast.summary_list[::-1] if getattr(blast, 'summary_list', None) else [])]
    report.add_table(f'BLAST {_blast_base} - Assembly Identification',
                     ['nt base count', 'contigs', 'Description'], _blast_rows, num_cols={0, 1},
                     footer=f'Results provided by: BLAST {_blast_base} database', compact=True)
    # (BTV serotyping section, if any, is added right here in display order below.)
    summary_log.start_section('Organism-Specific Analysis')
    # Generic, registry-driven handler for segmented / multi-chromosome
    # organisms that do NOT need sub-species splitting or serotyping (Isavirus,
    # Apicomplexa, and any future organism). Adding one = declare an
    # OrganismProfile and register it in organism_registry.PROFILES; no new
    # dispatch code. Orbivirus is intentionally excluded (BESPOKE_TAXA) and kept
    # on its dedicated branch below because it needs BTV/EHD splitting + serotyping.
    if args.taxon in PROFILES and args.taxon not in BESPOKE_TAXA:
        profile = PROFILES[args.taxon]
        logger.info(f"Running {profile.name}-specific analysis (generic segmented handler)...")
        args.specific = args.taxon
        segment_total_expected = profile.n_segments

        try:
            engine = SegmentedEngine(profile)
            hits = {}
            dest = open('all_downloaded.fasta', 'a')

            for seq_id, seq_info in coverage_graph.alignment_stats.items():
                header = seq_info['header'].lower()
                if any(kw in header for kw in profile.match_keywords):
                    logger.debug(f"Found {profile.name} match: {seq_id}")
                    hits[seq_id] = seq_info

            logger.debug(f"Found {len(hits)} {profile.name} sequences")

            if hits:
                ordered = engine.run(alignment_stats=hits)

                # Capture segment mapping (+ confidence) so the reference-guided
                # naming, alignment table, coverage banner, and Excel can badge
                # tentative slots — mirroring the Orbivirus path.
                for seq_id, info in (ordered or {}).items():
                    accession = str(seq_id).split('.')[0]
                    global_segment_mapping[accession] = {
                        'segment_number': info['segment_number'],
                        'full_accession': str(seq_id),
                        'description': info.get('header', str(seq_id)),
                        'confidence': info.get('confidence', 'high'),
                    }

                logger.debug("Generating alignment statistics for concatenated sequences...")
                coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
                coverage_graph.get_coverage_graph()
                pass  # alignment stats table is rendered later from CoverageGraphGenerator

                logger.debug("Adding sequences to combined FASTA...")
                with open('concatenated_specific.fasta', 'r') as source:
                    dest.write(source.read())

                if not args.debug:
                    os.remove('concatenated_specific.fasta')
                else:
                    logger.debug("Debug mode: Preserving concatenated_specific.fasta")
            else:
                logger.info(f"No {profile.name} sequences found for specific analysis")

            dest.close()

            if not hits:
                # Same guard as the Orbivirus branch: with nothing to download,
                # all_downloaded.fasta would stay empty and the reference-guided
                # alignment below would run against an empty FASTA. Write a
                # placeholder instead (recognized and hidden from tables/Excel).
                logger.warning(f"No {profile.name} sequences found among the BLAST hits — writing a "
                               "placeholder reference so the pipeline can finish")
                report.add_message(
                    f'{profile.name}-Specific Analysis',
                    f'No {profile.name} sequences were identified among the BLAST hits for this '
                    'sample. Reference-guided assembly and coverage graphs were skipped for the '
                    'organism-specific step.', level='warn')
                with open('all_downloaded.fasta', 'w') as f:
                    f.write('>No_Reference_Available\nNNNNNNNNNN\n')

            logger.info(f"{profile.name}-specific analysis completed")

        except Exception as e:
            logger.error(f"Error in {profile.name}-specific analysis: {e}")
            raise
    # if a specific test, Orbivirus, is an option
    elif args.taxon == 'Orbivirus':
        logger.info("Running Orbivirus-specific analysis...")
        logger.debug("Searching for Bluetongue virus and Epizootic hemorrhagic disease sequences")
        args.specific = args.taxon

        try:
            orbivirus_specific = Orbivirus_Specific()
            # Clear and populate the globally initialized dictionaries
            btv_dict.clear()
            ehv_dict.clear()
            dest = open('all_downloaded.fasta', 'a')

            for seq_id, seq_info in coverage_graph.alignment_stats.items():
                header = seq_info['header'].lower()
                logger.debug(f"Checking Orbivirus header: {header}")
                # Same species rule as the report sections and serotyping.
                species = _orbivirus_species(header)
                if species == 'BTV':
                    logger.debug(f"Found Bluetongue virus match: {seq_id}")
                    btv_dict[seq_id] = seq_info
                elif species == 'EHD':
                    logger.debug(f"Found Epizootic hemorrhagic disease match: {seq_id}")
                    ehv_dict[seq_id] = seq_info

            logger.debug(f"Found {len(btv_dict)} Bluetongue virus and {len(ehv_dict)} Epizootic hemorrhagic disease sequences")
            if btv_dict:
                logger.info("Processing Bluetongue virus sequences...")
                btv_ordered = orbivirus_specific.run(alignment_stats=btv_dict)
                # Per-segment confidence ('high' strict / 'tentative' cautious
                # fallback), keyed by base accession, to carry into the mapping.
                btv_confidence = {str(k).split('.')[0]: v.get('confidence', 'high')
                                  for k, v in (btv_ordered or {}).items()}

                # Capture segment mapping after processing to preserve throughout pipeline
                print(f"\n{bcolors.YELLOW}=== CAPTURING BLUETONGUE VIRUS SEGMENT ASSIGNMENTS ==={bcolors.ENDC}")
                logger.debug("Capturing segment mapping for BTV sequences...")

                with open('concatenated_specific.fasta', 'r') as fasta_file:
                    from Bio import SeqIO
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        # Extract accession from FASTA ID
                        accession = record.id.split('.')[0]  # Get base accession without version

                        # Extract segment number from description (canonical resolver:
                        # explicit "segment N" notation, falling back to protein codes).
                        segment_num = resolve_segment_number(record.description)
                        if segment_num is not None:
                            confidence = btv_confidence.get(accession, 'high')
                            global_segment_mapping[accession] = {
                                'segment_number': segment_num,
                                'full_accession': record.id,
                                'description': record.description,
                                'confidence': confidence
                            }
                            _conf_tag = '' if confidence == 'high' else f' {bcolors.YELLOW}(tentative){bcolors.ENDC}'
                            print(f"  📍 {bcolors.GREEN}{accession}{bcolors.ENDC} → {bcolors.BLUE}Segment {segment_num}{bcolors.ENDC}{_conf_tag}")
                            logger.debug(f"  Mapped {accession} → segment {segment_num} ({confidence})")
                        else:
                            print(f"  ⚠️  Could not extract segment for {accession}: {record.description}")
                            logger.debug(f"  Could not extract segment for {accession}: {record.description}")

                print(f"\n{bcolors.YELLOW}SEGMENT MAPPING DATA STRUCTURE:{bcolors.ENDC}")
                for acc, info in sorted(global_segment_mapping.items(), key=lambda x: x[1]['segment_number']):
                    print(f"  {bcolors.WHITE}{acc}{bcolors.ENDC}: {{'segment_number': {bcolors.BLUE}{info['segment_number']}{bcolors.ENDC}, 'full_accession': '{info['full_accession']}'}}")
                print(f"{bcolors.YELLOW}==================================================={bcolors.ENDC}\n")

                logger.debug("Generating alignment statistics for Bluetongue virus...")
                coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
                coverage_graph.get_coverage_graph()
                pass  # alignment stats table is rendered later from CoverageGraphGenerator

                logger.debug("Adding Bluetongue virus sequences to combined FASTA...")
                with open('concatenated_specific.fasta', 'r') as source:
                    dest.write(source.read())

                # Only remove if not in debug mode
                if not args.debug:
                    os.remove('concatenated_specific.fasta')
                else:
                    logger.debug("Debug mode: Preserving BTV concatenated_specific.fasta")

            if ehv_dict:
                logger.info("Processing Epizootic hemorrhagic disease sequences...")
                ehv_ordered = orbivirus_specific.run(alignment_stats=ehv_dict)
                ehv_confidence = {str(k).split('.')[0]: v.get('confidence', 'high')
                                  for k, v in (ehv_ordered or {}).items()}

                # Capture segment mapping for EHV sequences
                print(f"\n{bcolors.YELLOW}=== CAPTURING EHV SEGMENT ASSIGNMENTS ==={bcolors.ENDC}")
                logger.debug("Capturing segment mapping for EHV sequences...")

                with open('concatenated_specific.fasta', 'r') as fasta_file:
                    from Bio import SeqIO
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        # Extract accession from FASTA ID
                        accession = record.id.split('.')[0]  # Get base accession without version

                        # Extract segment number from description (canonical resolver:
                        # explicit "segment N" notation, falling back to protein codes).
                        segment_num = resolve_segment_number(record.description)
                        if segment_num is not None:
                            confidence = ehv_confidence.get(accession, 'high')
                            global_segment_mapping[accession] = {
                                'segment_number': segment_num,
                                'full_accession': record.id,
                                'description': record.description,
                                'confidence': confidence
                            }
                            _conf_tag = '' if confidence == 'high' else f' {bcolors.YELLOW}(tentative){bcolors.ENDC}'
                            print(f"  📍 {bcolors.GREEN}{accession}{bcolors.ENDC} → {bcolors.BLUE}Segment {segment_num}{bcolors.ENDC}{_conf_tag}")
                            logger.debug(f"  Mapped {accession} → segment {segment_num} ({confidence})")
                        else:
                            print(f"  ⚠️  Could not extract segment for {accession}: {record.description}")
                            logger.debug(f"  Could not extract segment for {accession}: {record.description}")

                print(f"\n{bcolors.YELLOW}UPDATED SEGMENT MAPPING DATA STRUCTURE:{bcolors.ENDC}")
                for acc, info in sorted(global_segment_mapping.items(), key=lambda x: x[1]['segment_number']):
                    print(f"  {bcolors.WHITE}{acc}{bcolors.ENDC}: {{'segment_number': {bcolors.BLUE}{info['segment_number']}{bcolors.ENDC}, 'full_accession': '{info['full_accession']}'}}")
                print(f"{bcolors.YELLOW}============================================={bcolors.ENDC}\n")

                logger.debug("Generating alignment statistics for Epizootic hemorrhagic disease...")
                coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
                coverage_graph.get_coverage_graph()
                pass  # alignment stats table is rendered later from CoverageGraphGenerator

                logger.debug("Adding EHD sequences to combined FASTA...")
                with open('concatenated_specific.fasta', 'r') as source:
                    dest.write(source.read())

                # Only remove if not in debug mode
                if not args.debug:
                    os.remove('concatenated_specific.fasta')
                else:
                    logger.debug("Debug mode: Preserving EHD concatenated_specific.fasta")

            dest.close()

            if not btv_dict and not ehv_dict:
                logger.warning("No Bluetongue virus or Epizootic hemorrhagic disease sequences found "
                                "among the BLAST hits — the recovered Orbivirus reads did not match "
                                "either known species header. Writing a placeholder reference so the "
                                "pipeline can finish instead of aligning against an empty FASTA.")
                report.add_message(
                    'Orbivirus-Specific Analysis',
                    'No Bluetongue virus or Epizootic hemorrhagic disease sequences were identified '
                    'among the BLAST hits for this Orbivirus-classified sample. Reference-guided '
                    'assembly and coverage graphs were skipped for the organism-specific step.',
                    level='warn')
                with open('all_downloaded.fasta', 'w') as f:
                    f.write('>No_Reference_Available\nNNNNNNNNNN\n')

            logger.info("Orbivirus-specific analysis completed")

        except Exception as e:
            logger.error(f"Error in Orbivirus-specific analysis: {e}")
            raise
    else:
        def get_last_accession(filename):
            """
            Extract the NCBI accession number from the last line of a BLAST summary file.
            
            Args:
                filename (str): Path to the BLAST summary file
                
            Returns:
                str: The NCBI accession number from the last line
            """
            try:
                with open(filename, 'r') as file:
                    # Read all lines and get the last one
                    lines = file.readlines()
                    if not lines:
                        raise ValueError(f"BLAST summary file {filename} is empty")
                    
                    last_line = lines[-1].strip()
                    if not last_line:
                        raise ValueError(f"Last line of {filename} is empty")
                    
                    # Split the line by spaces and get the third element
                    # which contains the accession number
                    parts = last_line.split()
                    if len(parts) < 3:
                        raise ValueError(f"Invalid format in last line of {filename}: {last_line}")
                    
                    accession = parts[2]
                    print(f"Extracted accession from BLAST results: {accession}")
                    return accession
                    
            except (FileNotFoundError, IndexError, ValueError) as e:
                print(f"Error reading BLAST summary file {filename}: {e}")
                raise
        
        # Get the top accession from BLAST results
        logger.info("Downloading reference sequence for standard analysis...")
        try:
            accession = get_last_accession(blast.blast_summary_file)
            logger.info(f'Downloading reference: {accession}')
            logger.debug(f"BLAST summary file: {blast.blast_summary_file}")

            # Download with proper error handling and rate limiting
            download = Downloader(accession)
            max_retries = 3
            success = False

            for attempt in range(max_retries):
                logger.debug(f"Download attempt {attempt + 1}/{max_retries} for {accession}")
                try:
                    download.fasta()
                    success = True
                    logger.info(f"✓ Successfully downloaded {accession}")
                    break
                except Exception as e:
                    logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
                    if attempt < max_retries - 1:
                        wait_time = (attempt + 1) * 5  # 5, 10, 15 seconds
                        logger.debug(f"Waiting {wait_time} seconds before retry...")
                        sleep(wait_time)
            
            if success:
                # Rename the downloaded file
                downloaded_file = f'{download.accession}.fasta'
                if os.path.exists(downloaded_file):
                    file_size = os.path.getsize(downloaded_file)
                    logger.debug(f"Downloaded file size: {file_size} bytes")
                    os.rename(downloaded_file, "all_downloaded.fasta")
                    logger.info(f"Renamed {downloaded_file} to all_downloaded.fasta")
                else:
                    logger.error(f"Downloaded file {downloaded_file} not found after successful download")
                    raise FileNotFoundError(f"Downloaded file {downloaded_file} not found")
            else:
                logger.error(f"Failed to download {accession} after {max_retries} attempts")
                raise Exception(f"Failed to download {accession} after {max_retries} attempts")

        except Exception as e:
            logger.error(f"Error in reference download process: {e}")
            logger.warning("Creating empty reference file to prevent pipeline failure")
            with open("all_downloaded.fasta", 'w') as f:
                f.write(f">No_Reference_Available\nNNNNNNNNNN\n")
            logger.debug("Empty reference file created")
        
        args.specific = args.taxon
        
    summary_log.end_section()
    logger.info(f"\n{bcolors.BLUE}Generating analysis reports...{bcolors.ENDC}")

    try:
        logger.debug("Adding assembly stats to Excel report...")
        assembler.excel(excel_stats.excel_dict)

        logger.debug("Adding BLAST results to Excel report...")
        blast.excel(excel_stats.excel_dict)

        logger.debug("Excel report data compilation completed")

    except Exception as e:
        logger.error(f"Error adding data to Excel report: {e}")
        raise

    # Clean up temporary directory if it exists
    if not args.debug:
        temp_dir = Path('./temp')
        if temp_dir.exists():
            logger.debug(f"Removing temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir)
    else:
        logger.debug("Debug mode: Preserving temporary directory for troubleshooting")

    pipeline_duration = datetime.datetime.now() - start_time
    logger.info(f"\n{bcolors.GREEN}Initial pipeline phase completed successfully!{bcolors.ENDC}")
    logger.info(f"Duration so far: {pipeline_duration}")
    logger.debug("Proceeding to organism-specific analysis phase...")

####################################################################################################################
    btv_serotyping = None
    if args.specific:
        summary_log.start_section('Reference-Guided Assembly')
        logger.info("Starting reference-guided assembly phase...")
        logger.debug("Running alignment against downloaded reference...")

        try:
            alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='all_downloaded.fasta', skip_assembly=True, debug=args.debug)
            alignment.run()
            logger.debug(f"Alignment completed - reference: {alignment.reference}")
            logger.debug(f"VCF file: {alignment.zc_vcf}")

            logger.debug("Running reference-guided assembly...")
            reference_guided_assembly = Reference_Guided_Assembly(FASTA=f'{alignment.reference}', vcf=f'{alignment.zc_vcf}', iupac=True, output_name='consensus.fasta')
            logger.debug("Reference-guided assembly completed")

            logger.debug("Running BLAST on consensus sequence...")
            blast = Blast_Fasta(FASTA='consensus.fasta', num_alignment=1, blast_db=args.blast_db)
            logger.debug("Consensus BLAST completed")

        except Exception as e:
            logger.error(f"Error in reference-guided assembly phase: {e}")
            raise

        # IMPROVED DOWNLOAD SECTION - REPLACE THE OLD LOOP WITH THIS:
        logger.info("Downloading top hit reference genomes...")
        logger.info(f"Found {len(blast.top_hit_acc)} accessions to download")
        logger.debug(f"Top hit accessions: {blast.top_hit_acc}")
        
        downloaded_count = 0
        failed_downloads = []          # list of (accession, error_message)
        reference_download_failed = False

        with open('top_downloaded_genomes.fasta', 'w') as outfile:
            for i, acc in enumerate(blast.top_hit_acc, 1):
                logger.info(f"Downloading {i}/{len(blast.top_hit_acc)}: {acc}")

                # PROACTIVE rate limiting - wait BEFORE each request (except first)
                if i > 1:
                    wait_time = 3  # Fixed 3-second delay between downloads
                    logger.debug(f"Waiting {wait_time} seconds (NCBI rate limiting)...")
                    time.sleep(wait_time)
                
                download = Downloader(acc)
                max_retries = 3
                success = False
                last_error = None

                # Retry with exponential backoff
                for attempt in range(max_retries):
                    try:
                        download.fasta()
                        success = True
                        downloaded_count += 1
                        logger.info(f"  ✓ Successfully downloaded {acc}")
                        break
                    except Exception as e:
                        last_error = e
                        logger.warning(f"  ⚠ Attempt {attempt + 1}/{max_retries} failed for {acc}: {e}")
                        if attempt < max_retries - 1:
                            # Exponential backoff: 10, 20, 30 seconds
                            backoff_time = (attempt + 1) * 10
                            logger.debug(f"    Waiting {backoff_time} seconds before retry...")
                            time.sleep(backoff_time)
                
                # Handle the downloaded file
                if success:
                    fasta_file = f'{acc}.fasta'
                    try:
                        if os.path.exists(fasta_file):
                            with open(fasta_file, 'r') as infile:
                                content = infile.read().rstrip()
                                if content:  # Only write if file has content
                                    outfile.write(content)
                                    outfile.write('\n')
                            os.remove(fasta_file)  # Clean up individual file
                        else:
                            logger.warning(f"  ⚠ Warning: {fasta_file} not found after download")
                            success = False
                    except Exception as e:
                        last_error = e
                        logger.error(f"  ⚠ Error processing {fasta_file}: {e}")
                        success = False

                if not success:
                    failed_downloads.append((acc, str(last_error) if last_error else 'unknown error'))
                    logger.error(f"  ✗ Failed to download {acc} after {max_retries} attempts")
        
        # Summary and error handling
        logger.info(f"\nDownload Summary:")
        logger.info(f"  Successfully downloaded: {downloaded_count}/{len(blast.top_hit_acc)}")

        if failed_downloads:
            logger.warning(f"  Failed downloads: {len(failed_downloads)}")
            logger.debug(f"  Failed accessions: {', '.join(acc for acc, _ in failed_downloads)}")

        # Check if we have any successful downloads
        if downloaded_count == 0 and not blast.top_hit_acc:
            # Nothing was requested at all: the BLAST of the reference-guided
            # consensus returned no hits (e.g. it reached BLAST_TIMEOUT_SECONDS), so
            # this is not an NCBI download problem and must not be reported as one.
            reference_download_failed = True
            logger.error("ERROR: The consensus BLAST returned no hits - no reference genome to download")
            title = 'Reference-Guided Coverage Not Performed'
            log_msg = ('Consensus BLAST returned no hits — reference-guided coverage skipped '
                       '(check the run log for "BLAST timed out").')
            detail = (
                "The BLAST search of the reference-guided consensus returned no hits, so no reference genome "
                "could be selected for the final alignment. Reference-guided alignment and coverage analysis were "
                "therefore SKIPPED. Taxonomic identification from Kraken, the de novo assembly, and BLAST above "
                "is unaffected.\n\n"
                "Check the run log for \"BLAST timed out\": a large consensus can exceed BLAST_TIMEOUT_SECONDS "
                f"({BLAST_TIMEOUT_SECONDS} s). Raise that limit, or confirm the BLAST database ({args.blast_db}) "
                "contains this organism, and re-run."
            )
        elif downloaded_count == 0:
            reference_download_failed = True
            logger.error("ERROR: No reference genomes could be downloaded from NCBI!")

            # Capture the ACTUAL error(s) so the failure is diagnosable instead of
            # being silently replaced by a fake all-N reference. BLAST already
            # identified the organism (see the BLAST section) — only the reference
            # fetch failed, so identification remains valid; we just cannot do
            # reference-guided coverage.
            attempted = [acc for acc, _ in failed_downloads] or list(blast.top_hit_acc)
            sample_err = next((err for _, err in failed_downloads if err), 'unknown error')
            acc_lines = "\n".join(f"    • {acc}: {err}" for acc, err in failed_downloads[:5])
            if len(failed_downloads) > 5:
                acc_lines += f"\n    • …and {len(failed_downloads) - 5} more"

            title = 'Reference Genome Download Failed'
            log_msg = ("Reference genome download failed — reference-guided coverage skipped. "
                       f"First error: {sample_err}. See report for remediation steps.")
            detail = (
                "A reference genome was identified by BLAST (see the BLAST section above), but none of the "
                f"candidate reference sequences could be downloaded from NCBI ({len(attempted)} attempted, "
                "3 retries each). Reference-guided alignment and coverage analysis were therefore SKIPPED. "
                "Taxonomic identification from Kraken, the de novo assembly, and BLAST above is unaffected.\n\n"
                "Accessions attempted (with the error returned):\n"
                f"{acc_lines}\n\n"
                "This is almost always a transient NCBI E-utilities outage, network interruption, or rate\n"
                "limiting during the run — the download step is normally reliable. To resolve:\n"
                "  1. Re-run the sample; transient NCBI/network issues usually clear on their own.\n"
                "  2. Confirm the compute node can reach NCBI:\n"
                "       python -c \"from Bio import Entrez; Entrez.email='you@org'; "
                "print(Entrez.efetch(db='nucleotide', id='NC_000962.3', rettype='fasta', retmode='text').read()[:60])\"\n"
                "  3. Set a real NCBI identity to avoid throttling (raises the rate limit 3→10/s):\n"
                "       export NCBI_EMAIL='you@org'\n"
                "       export NCBI_API_KEY='<key from https://www.ncbi.nlm.nih.gov/account/settings/>'\n"
                "  4. If NCBI is unreachable from this network, download the accession(s) elsewhere, save them\n"
                "     as 'top_downloaded_genomes.fasta' in the sample folder, and re-run."
            )

        if downloaded_count == 0:
            logger.error(detail)
            try:
                summary_log.log_error(log_msg)
            except Exception:
                pass
            report.add_message(title, detail, level='error')

            # Write a minimal placeholder ONLY to keep the downstream mechanics
            # from crashing; it is recognized and suppressed everywhere it would
            # otherwise be shown (no fake coverage numbers reach the report/Excel).
            with open('top_downloaded_genomes.fasta', 'w') as f:
                f.write(">No_Reference_Available\n")
                f.write("N" * 1000 + "\n")
            logger.debug("Placeholder reference written (suppressed downstream)")
        else:
            logger.info(f"Proceeding with {downloaded_count} successfully downloaded references...")
        # END OF IMPROVED DOWNLOAD SECTION
        
        # Clean up temp directory unless in debug mode
        if not args.debug:
            if os.path.exists('./temp'):
                shutil.rmtree('./temp')
                logger.debug("Removed temp directory after downloads")
        else:
            logger.debug("Debug mode: Preserving temp directory after downloads")

        logger.debug("Running final alignment against top downloaded genomes...")
        alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='top_downloaded_genomes.fasta', skip_assembly=True, debug=args.debug)
        alignment.run()
        logger.debug("Final alignment completed")

        logger.debug("Running final reference-guided assembly...")
        reference_guided_assembly = Reference_Guided_Assembly(FASTA=f'{alignment.reference}', vcf=f'{alignment.zc_vcf}', iupac=True, output_name='consensus.fasta')
        logger.debug("Final reference-guided assembly completed")

        final_consensus = f'{assembler.sample_name}_reference_guided.fasta'
        logger.debug(f"Final consensus output: {final_consensus}")
        print(f"\n{bcolors.PURPLE}=== APPLYING PRESERVED SEGMENT ASSIGNMENTS TO REFERENCE-GUIDED FASTA ==={bcolors.ENDC}")
        logger.debug("Processing consensus sequences and renaming headers...")

        print(f"Using preserved segment mapping with {bcolors.BLUE}{len(global_segment_mapping)}{bcolors.ENDC} entries:")
        for acc, info in sorted(global_segment_mapping.items(), key=lambda x: x[1]['segment_number']):
            print(f"  {bcolors.WHITE}{acc}{bcolors.ENDC} → {bcolors.BLUE}segment {info['segment_number']}{bcolors.ENDC}")
            logger.debug(f"  {acc} → segment {info['segment_number']}")

        renamed_fastas=[]
        # A mixed BTV + EHD sample has a segment N from each virus; tag the species
        # in those IDs ({sample}_BTV_segment2 / {sample}_EHD_segment2) so the
        # reference-guided FASTA never carries duplicate sequence IDs.
        mixed_orbivirus = args.taxon == 'Orbivirus' and bool(btv_dict) and bool(ehv_dict)
        record_dict = SeqIO.to_dict(SeqIO.parse(f'consensus.fasta', "fasta"))
        print(f"\nProcessing {bcolors.YELLOW}{len(record_dict)}{bcolors.ENDC} consensus sequences for renaming:")
        logger.debug(f"Processing {len(record_dict)} consensus sequences")

        for acc, seq in record_dict.items():
            logger.debug(f"Renaming sequence: {acc}")

            # Extract segment information using preserved mapping first
            segment_info = ""
            original_description = seq.description

            # First try to use the preserved segment mapping
            base_accession = acc.split('.')[0]  # Remove version number if present
            if base_accession in global_segment_mapping:
                segment_num = global_segment_mapping[base_accession]['segment_number']
                segment_info = f"_segment{segment_num}"
                print(f"  ✅ {bcolors.GREEN}{base_accession}{bcolors.ENDC} → Using preserved mapping → {bcolors.BLUE}segment {segment_num}{bcolors.ENDC}")
                logger.debug(f"  Using preserved segment mapping: {base_accession} → segment {segment_num}")
            else:
                # Fallback to extracting from description if not in mapping
                print(f"  ⚠️  {bcolors.YELLOW}{base_accession}{bcolors.ENDC} → Not in preserved mapping, extracting from description")
                logger.debug(f"  Segment mapping not found for {base_accession}, extracting from description")

                # Try the canonical resolver first (explicit "segment N"/"seg N"
                # notation, then bare protein codes) so this matches the same
                # segment number the coverage-graph banner and serotyping use.
                canonical_segment = resolve_segment_number(original_description)
                if canonical_segment is not None:
                    segment_info = f"_segment{canonical_segment}"
                    print(f"    → Resolved via canonical resolver: segment {canonical_segment}")
                else:
                    # Try alternative patterns like "seg 1", "S1", etc.
                    seg_match = re.search(r'seg(?:ment)?\s*[:\-]?\s*([0-9]+)', original_description, re.IGNORECASE)
                    if seg_match:
                        segment_info = f"_segment{seg_match.group(1)}"
                        print(f"    → Extracted from 'seg X' pattern: segment {seg_match.group(1)}")
                    else:
                        # Map protein names to segment numbers (critical for NS proteins!)
                        protein_to_segment = {
                            'VP1': 1, 'RNA-dependent RNA polymerase': 1,
                            'VP2': 2,
                            'VP3': 3, 'T2': 3,
                            'VP4': 4,
                            'NS1': 5, 'nonstructural protein NS1': 5,  # ← This catches X17041!
                            'VP5': 6,
                            'VP7': 7,
                            'NS2': 8, 'nonstructural protein NS2': 8,
                            'VP6': 9,
                            'NS3': 10, 'nonstructural protein NS3': 10
                        }

                        found_segment = None
                        for protein, seg_num in protein_to_segment.items():
                            if protein.lower() in original_description.lower():
                                segment_info = f"_segment{seg_num}"
                                found_segment = seg_num
                                print(f"    → Mapped protein '{protein}' to segment {seg_num}")
                                break

                        if not found_segment:
                            # Try pattern like "S1", "L1", etc.
                            s_match = re.search(r'\b[SLM](\d+)\b', original_description, re.IGNORECASE)
                            if s_match:
                                segment_info = f"_segment{s_match.group(1)}"
                                print(f"    → Extracted from 'SX' pattern: segment {s_match.group(1)}")
                            else:
                                # FALLBACK: Check if this accession exists anywhere in the original captured mapping
                                fallback_segment = None
                                for mapped_acc, mapping_info in global_segment_mapping.items():
                                    # Check if current accession matches any in the original mapping
                                    if base_accession == mapped_acc or acc == mapped_acc:
                                        fallback_segment = mapping_info['segment_number']
                                        segment_info = f"_segment{fallback_segment}"
                                        print(f"    → {bcolors.GREEN}Found in original data structure: {mapped_acc} → segment {fallback_segment}{bcolors.ENDC}")
                                        break

                                if not fallback_segment:
                                    # ABSOLUTE LAST RESORT: use accession but warn user
                                    clean_acc = re.sub(r'[^\w\-]', '_', acc)[:20]
                                    segment_info = f"_{clean_acc}"
                                    print(f"    → {bcolors.RED}FINAL WARNING: Could not determine segment number! Using accession: {clean_acc}{bcolors.ENDC}")
                                    print(f"    → Description was: {original_description}")
                                    print(f"    → This accession was not found in original captured assignments:")
                                    for mapped_acc, mapping_info in sorted(global_segment_mapping.items(), key=lambda x: x[1]['segment_number']):
                                        print(f"      {mapped_acc} → segment {mapping_info['segment_number']}")
                                    print(f"    → Please check if this should have a segment number!")

            # Create unique, descriptive sequence ID
            species = _orbivirus_species(original_description) if mixed_orbivirus else None
            if species and segment_info.startswith('_segment'):
                segment_info = f'_{species}{segment_info}'
            unique_id = f'{assembler.sample_name}{segment_info}'

            seq.name = unique_id
            seq.id = unique_id
            seq.description = f'guided by {original_description}'
            renamed_fastas.append(seq)

            print(f"    → Final ID: {bcolors.PURPLE}{unique_id}{bcolors.ENDC}")
            logger.debug(f"  Renamed to: {unique_id}")
            logger.debug(f"  Full description: {seq.description}")

        print(f"\n{bcolors.PURPLE}FINAL REFERENCE-GUIDED FASTA SEQUENCES:{bcolors.ENDC}")
        for seq in renamed_fastas:
            print(f"  {bcolors.BLUE}>{seq.id}{bcolors.ENDC} {seq.description}")
        print(f"{bcolors.PURPLE}================================================================{bcolors.ENDC}\n")

        SeqIO.write(renamed_fastas, final_consensus, "fasta")
        logger.debug(f"Final consensus written to {final_consensus}")

        # Only remove consensus.fasta if not in debug mode
        if not args.debug:
            os.remove('consensus.fasta')
            logger.debug("Removed intermediate consensus.fasta")
        else:
            logger.debug("Debug mode: Preserving intermediate consensus.fasta")

        # Clean up temp directory unless in debug mode
        if not args.debug:
            if os.path.exists('./temp'):
                shutil.rmtree('./temp')
                logger.debug("Removed temp directory after consensus processing")
        else:
            logger.debug("Debug mode: Preserving temp directory after consensus processing")
        
        # BTV Serotyping - run when sample is Orbivirus BTV
        btv_serotyping = None
        if args.taxon == 'Orbivirus' and btv_dict:
            try:
                if os.path.exists(final_consensus):
                    logger.info("Running BTV serotyping analysis...")
                    btv_serotyping = BTVSerotyping(
                        consensus_fasta=final_consensus,
                        output_dir='btv_serotype',
                        debug=args.debug
                    )
                    btv_serotyping.run()

                    # Add serotyping section (display order is naturally preserved:
                    # it lands after the BLAST section and before coverage graphs).
                    # Flag the serotyping protein(s) if their BTV segment was assigned
                    # by the cautious low-confidence fallback (VP2=seg2, VP5=seg6).
                    tentative_segments = tentative_segments_by_species(global_segment_mapping).get('BTV', set())
                    report.add_serotype(btv_serotyping.get_consensus_serotype(),
                                        btv_serotyping.interpretation,
                                        btv_serotyping.predictions,
                                        tentative_segments=tentative_segments)
                    logger.info(f"BTV serotype prediction: {btv_serotyping.get_consensus_serotype()}")
                else:
                    logger.warning(f"Skipping BTV serotyping - consensus FASTA not found: {final_consensus}")
            except Exception as e:
                logger.error(f"BTV serotyping failed: {e}")

        def clean_fasta_headers(input_file, output_file):
            """
            Replace only the first 4 spaces with underscores in FASTA headers while preserving sequence data.
            
            Args:
                input_file (str): Path to input FASTA file
                output_file (str): Path to output FASTA file
            """
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    if line.startswith('>'):  # Header line
                        # Replace only first 4 spaces with underscores in header
                        cleaned_header = line.strip()
                        for _ in range(4):
                            cleaned_header = cleaned_header.replace(' ', '_', 1)
                        fout.write(cleaned_header + '\n')
                    else:  # Sequence line
                        fout.write(line)
                        
        input_file = final_consensus
        output_file = "output.fasta"
        clean_fasta_headers(input_file, output_file)
        
        logger.debug("Running final coverage analysis on cleaned consensus...")
        coverage_graph = Coverage_Graph(FASTA="output.fasta", FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=args.debug)
        coverage_graph.get_coverage_graph()
        logger.debug("Final coverage analysis completed")

        # Only remove output.fasta if not in debug mode
        if not args.debug:
            os.remove("output.fasta")
            logger.debug("Removed intermediate output.fasta")
        else:
            logger.debug("Debug mode: Preserving intermediate output.fasta")

        summary_log.end_section()  # end Reference-Guided Assembly

        # Generate additional PNG coverage graphs using reference-guided assembly results
        summary_log.start_section('Coverage Graph Generation (PNG)')
        logger.info("Generating additional PNG coverage graphs...")
        try:
            # Define input files based on what the pipeline created. parse_reads
            # sanitizes the taxon (spaces -> underscores) when naming the extracted
            # reads, so reuse those exact names (parser.r1_out/.r2_out) instead of
            # rebuilding from the raw args.taxon. Rebuilding with spaces made
            # multi-word taxa (e.g. "Mycobacterium tuberculosis complex") never match,
            # which silently skipped the coverage graph. parser.r2_out is None for
            # single-end runs.
            reference_guided_fasta = f"{assembler.sample_name}_reference_guided.fasta"
            taxon_r1 = f"{parser.r1_out}.gz"
            taxon_r2 = f"{parser.r2_out}.gz" if parser.r2_out else None

            logger.debug(f"Reference-guided FASTA: {reference_guided_fasta}")
            logger.debug(f"Taxon FASTQ R1: {taxon_r1}")
            logger.debug(f"Taxon FASTQ R2: {taxon_r2}")

            # Check if all required files exist (R2 optional for single-end runs)
            required_files = [reference_guided_fasta, taxon_r1] + ([taxon_r2] if taxon_r2 else [])
            if all(os.path.exists(f) for f in required_files):
                logger.debug("All required files found, generating PNG coverage graphs...")

                # Generate coverage graphs
                coverage_png_generator = CoverageGraphGenerator(
                    fasta_file=reference_guided_fasta,
                    fastq_r1=taxon_r1,
                    fastq_r2=taxon_r2,
                    debug=args.debug
                )

                # Run the coverage graph generation
                coverage_png_generator.generate_coverage_graphs()
                logger.info("PNG coverage graphs generated successfully")

                # Drop synthetic "no reference available" placeholders (written
                # upstream when reference download/selection fails) from the stats
                # at the source, so they never pollute the alignment table, coverage
                # graphs, Excel metrics, or Coverage QC. If nothing real remains,
                # surface a clear notice instead of a bogus 1000-N reference row.
                _raw_align = getattr(coverage_png_generator, '_last_alignment_stats', {}) or {}
                _placeholder_ids = [rid for rid, st in _raw_align.items()
                                    if _is_no_reference_placeholder(rid)
                                    or _is_no_reference_placeholder(st.get('header', ''))]
                for _rid in _placeholder_ids:
                    _raw_align.pop(_rid, None)
                no_reference_available = bool(_placeholder_ids) and not _raw_align
                # A detailed, actionable message was already added at the point of
                # failure when the download step failed (reference_download_failed).
                # Only add the generic fallback notice for other placeholder origins.
                if no_reference_available and not reference_download_failed:
                    logger.warning("Reference-guided coverage skipped: no reference genome was available")
                    report.add_message(
                        'Reference-Guided Coverage Analysis Not Performed',
                        'A candidate reference genome was identified by BLAST (see above), but the '
                        'reference sequence could not be downloaded, so reference-guided alignment and '
                        'coverage analysis were skipped. Identification is based on the de novo assembly '
                        'and BLAST results.', level='warn')

                # Variant calling failed (e.g. freebayes broken) — surface it loudly
                # instead of silently reporting "0 SNPs". Coverage is still valid;
                # only the SNP column/graph track is unavailable.
                _snp_err = getattr(coverage_png_generator, '_last_snp_call_error', None)
                if _snp_err and not no_reference_available:
                    logger.error(f"SNP calling unavailable: {_snp_err}")
                    try:
                        summary_log.log_error(f"Variant calling failed — SNP counts unavailable. {_snp_err}")
                    except Exception:
                        pass
                    report.add_message(
                        'SNP Calling Unavailable',
                        'Coverage was computed successfully, but variant (SNP) calling could not run, so '
                        'SNP counts are shown as “not called” (not zero) throughout this report.\n\n'
                        f'Reason: {_snp_err}\n\n'
                        'To fix: repair the variant-caller install in the pipeline environment, then '
                        're-run — coverage is unaffected, only the SNP results are missing.',
                        level='error')

                # Add coverage graphs to the report, grouped by virus type into their own sections
                if hasattr(coverage_png_generator, 'png_file_paths') and coverage_png_generator.png_file_paths:
                    logger.debug(f"Adding {len(coverage_png_generator.png_file_paths)} PNG coverage graphs to the report...")

                    # Group PNG coverage graphs by virus type (BTV vs EHD)
                    btv_pngs = {}
                    ehd_pngs = {}
                    other_pngs = {}

                    for variable_name, png_path in coverage_png_generator.png_file_paths.items():
                        description = coverage_png_generator.png_file_descriptions.get(variable_name, variable_name)
                        if _is_btv(description):
                            btv_pngs[variable_name] = png_path
                        elif _is_ehd(description):
                            ehd_pngs[variable_name] = png_path
                        else:
                            other_pngs[variable_name] = png_path

                    # Segments assigned by the cautious low-confidence fallback,
                    # captured in global_segment_mapping and grouped per organism —
                    # badged '(tentative)' in the alignment table, coverage banner,
                    # and serotyping.
                    tentative_by_species = tentative_segments_by_species(global_segment_mapping)

                    # 'Table 1: Alignment statistics' for every recovered reference.
                    add_alignment_stats_table(
                        report, getattr(coverage_png_generator, '_last_alignment_stats', {}),
                        tentative_by_species=tentative_by_species)

                    # Coverage graph sections (interactive HTML + static PDF), grouped by virus type.
                    descs = coverage_png_generator.png_file_descriptions
                    if btv_pngs:
                        logger.info(f"Adding {len(btv_pngs)} Bluetongue virus coverage graphs")
                        add_virus_section(report, 'Bluetongue Virus Coverage Graphs',
                                          coverage_png_generator, list(btv_pngs.keys()), descs,
                                          tentative_segments=tentative_by_species.get('BTV', set()))
                    if ehd_pngs:
                        logger.info(f"Adding {len(ehd_pngs)} Epizootic Hemorrhagic Disease coverage graphs")
                        add_virus_section(report, 'Epizootic Hemorrhagic Disease Coverage Graphs',
                                          coverage_png_generator, list(ehd_pngs.keys()), descs,
                                          tentative_segments=tentative_by_species.get('EHD', set()))
                    if other_pngs:
                        logger.info(f"Adding {len(other_pngs)} non-Orbivirus coverage graphs (new style) + alignment table")
                        add_virus_section(report, 'Coverage Graphs',
                                          coverage_png_generator, list(other_pngs.keys()), descs,
                                          show_segment_status=False,
                                          tentative_segments=tentative_by_species.get('other', set()))

                    total_pngs = len(btv_pngs) + len(ehd_pngs) + len(other_pngs)
                    logger.info(f"Added {total_pngs} PNG coverage graphs to report (BTV: {len(btv_pngs)}, EHD: {len(ehd_pngs)}, other: {len(other_pngs)})")
                else:
                    logger.warning("No PNG files were generated by coverage graph generator")

            else:
                missing_files = [f for f in required_files if not os.path.exists(f)]
                logger.warning(f"Skipping PNG coverage graph generation - missing files: {missing_files}")
            summary_log.end_section()

        except Exception as e:
            summary_log.log_error(f"PNG coverage graph generation: {e}")
            summary_log.end_section('FAILED')
            logger.error(f"Error generating additional PNG coverage graphs: {e}")
            logger.error("PNG coverage graph generation FAILED - stopping pipeline execution")
            logger.error("This error indicates a problem with BWA alignment or samtools processing")
            logger.error("Check the debug output above for detailed error messages")
            raise SystemExit(f"PNG coverage graph generation failed: {e}")

####################################################################################################################

    # Final report generation: interactive HTML + PDF (WeasyPrint).
    logger.info("Finalizing reports (HTML + PDF)...")
    html_path, pdf_path = _write_reports(report, logger=logger, summary_log=summary_log,
                                         write_pdf=not args.no_pdf)
    if pdf_path or args.no_pdf:
        logger.info("Report generation completed successfully")
    else:
        logger.error("PDF report generation failed - check the summary log for diagnostics")

    try:
        logger.debug("Generating Excel statistics report...")

        # --- Add additional detail to the Excel summary ---

        # Virus type identification (BTV, EHD, or mix)
        if 'renamed_fastas' in dir() and renamed_fastas:
            has_btv = False
            has_ehd = False
            for seq in renamed_fastas:
                if _is_btv(seq.description):
                    has_btv = True
                elif _is_ehd(seq.description):
                    has_ehd = True
            if has_btv and has_ehd:
                excel_stats.excel_dict['Virus Type'] = 'mix'
            elif has_btv:
                excel_stats.excel_dict['Virus Type'] = 'BTV'
            elif has_ehd:
                excel_stats.excel_dict['Virus Type'] = 'EHD'

        # Segment mapping summary (for segmented viruses like Orbivirus)
        if global_segment_mapping:
            # Per organism, so a mixed BTV + EHD sample reports each virus's
            # segments (a pooled list read "1, 1, 2, 2, ..." and hid a segment
            # missing from only one of them). One organism keeps the plain format.
            by_species = defaultdict(list)
            for info in global_segment_mapping.values():
                species = _orbivirus_species(info.get('description', '')) or 'other'
                by_species[species].append(info['segment_number'])

            def _per_species(values):
                if len(values) == 1:
                    return next(iter(values.values()))
                return '; '.join(f'{sp}: {text}' for sp, text in sorted(values.items()))

            # Check which segments are missing from the organism's full set
            # (1-10 for Orbivirus; the profile's count for other segmented organisms)
            expected = set(range(1, segment_total_expected + 1))
            identified, counts, missing = {}, {}, {}
            for species, segment_numbers in by_species.items():
                segment_numbers = sorted(segment_numbers)
                identified[species] = ', '.join(str(s) for s in segment_numbers)
                counts[species] = str(len(segment_numbers))
                absent = sorted(expected - set(segment_numbers))
                missing[species] = ', '.join(str(s) for s in absent) if absent else 'None'
            excel_stats.excel_dict['Segments Identified'] = _per_species(identified)
            excel_stats.excel_dict['Segment Count'] = _per_species(counts)
            excel_stats.excel_dict['Segments Missing'] = _per_species(missing)

            # Segments assigned by the cautious low-confidence fallback (labelled by
            # species when a mixed sample has tentative segments in more than one).
            _tentative = tentative_segments_by_species(global_segment_mapping)
            if not _tentative:
                _tentative_text = 'None'
            elif len(_tentative) == 1:
                _tentative_text = ', '.join(str(s) for s in sorted(next(iter(_tentative.values()))))
            else:
                _tentative_text = '; '.join(f"{sp}: {', '.join(str(s) for s in sorted(segs))}"
                                            for sp, segs in sorted(_tentative.items()))
            excel_stats.excel_dict['Segments Tentative (low-confidence)'] = _tentative_text

        # Reference-guided assembly details (skip synthetic no-reference placeholders)
        _real_fastas = [seq for seq in renamed_fastas
                        if not (_is_no_reference_placeholder(seq.id)
                                or _is_no_reference_placeholder(seq.description))] \
            if 'renamed_fastas' in dir() and renamed_fastas else []
        if _real_fastas:
            # Per-segment summary with accession references
            segment_details = []
            for seq in _real_fastas:
                # Extract accession from description
                acc_match = re.search(r'guided by (\S+)', seq.description)
                acc = acc_match.group(1) if acc_match else 'unknown'
                segment_details.append(f"{_segment_label(seq.id, 'seg?')}:{acc}")
            excel_stats.excel_dict['Reference Accessions'] = '; '.join(segment_details)
            excel_stats.excel_dict['Total Reference Segments'] = str(len(_real_fastas))

        # Alignment duplication metrics (from final alignment)
        if 'alignment' in dir() and alignment:
            excel_stats.excel_dict['Read Pairs Examined'] = f'{alignment.READ_PAIRS_EXAMINED:,}'
            excel_stats.excel_dict['Unmapped Reads'] = f'{alignment.UNMAPPED_READS:,}'
            excel_stats.excel_dict['Percent Duplication'] = f'{alignment.PERCENT_DUPLICATION:.4f}'

        # Coverage graph generator stats (per-segment coverage)
        if 'coverage_png_generator' in dir() and coverage_png_generator and hasattr(coverage_png_generator, '_last_alignment_stats'):
            stats = coverage_png_generator._last_alignment_stats
            if stats:
                coverages = []
                genome_coverages = []
                snp_counts = []
                for ref_id, ref_stats in stats.items():
                    seg_label = _segment_label(ref_id, ref_id[:15])
                    coverages.append(f"{seg_label}:{ref_stats['mean_coverage']:.1f}X")
                    genome_coverages.append(f"{seg_label}:{ref_stats['percent_covered']:.1f}%")
                    # None = variant calling did not run (never reported as a real 0).
                    _snps = ref_stats.get('snp_count', 0)
                    snp_counts.append(f"{seg_label}:{'not called' if _snps is None else _snps}")
                excel_stats.excel_dict['Per-Segment Coverage'] = '; '.join(coverages)
                excel_stats.excel_dict['Per-Segment Depth of Coverage'] = '; '.join(coverages)
                excel_stats.excel_dict['Per-Segment Genome with Coverage'] = '; '.join(genome_coverages)
                excel_stats.excel_dict[f'Per-Segment SNPs (Q>{coverage_plot.GOOD_SNP_MIN_QUAL})'] = '; '.join(snp_counts)

                # Summary stats across all segments
                all_coverages = [s['mean_coverage'] for s in stats.values()]
                all_pct_covered = [s['percent_covered'] for s in stats.values()]
                # snp_count is None when variant calling did not run; int(None) here
                # used to abort the whole Excel report, so report 'not called' instead.
                all_snps = [s.get('snp_count', 0) for s in stats.values()]
                excel_stats.excel_dict['Min Segment Coverage'] = f'{min(all_coverages):.1f}X'
                excel_stats.excel_dict['Max Segment Coverage'] = f'{max(all_coverages):.1f}X'
                excel_stats.excel_dict['Mean Segment Coverage'] = f'{sum(all_coverages)/len(all_coverages):.1f}X'
                excel_stats.excel_dict['Min Genome Covered'] = f'{min(all_pct_covered):.1f}%'
                excel_stats.excel_dict[f'High-Quality SNPs (Q>{coverage_plot.GOOD_SNP_MIN_QUAL})'] = (
                    'not called' if any(n is None for n in all_snps)
                    else f'{sum(int(n) for n in all_snps):,}')

        # BTV Serotyping results
        if btv_serotyping and btv_serotyping.predictions:
            btv_serotyping.excel(excel_stats.excel_dict)

        # Coverage quality assessment (Pass/Marginal/Fail)
        if 'coverage_png_generator' in dir() and coverage_png_generator and hasattr(coverage_png_generator, '_last_alignment_stats'):
            stats = coverage_png_generator._last_alignment_stats
            if stats:
                fail_count = 0
                for ref_id, ref_stats in stats.items():
                    mean_cov = ref_stats['mean_coverage']
                    pct_covered = ref_stats['percent_covered']
                    if mean_cov < 10 or pct_covered < 98:
                        fail_count += 1
                if fail_count == 0:
                    excel_stats.excel_dict['Coverage QC'] = 'Pass'
                elif fail_count == 1:
                    excel_stats.excel_dict['Coverage QC'] = 'Marginal'
                else:
                    excel_stats.excel_dict['Coverage QC'] = 'Fail'

        # Pipeline runtime
        total_runtime = datetime.datetime.now() - start_time
        excel_stats.excel_dict['Pipeline Runtime'] = str(total_runtime).split('.')[0]

        excel_stats.post_excel()
        logger.info(f"Excel report saved as: {excel_stats.excel_filename}")
    except Exception as e:
        logger.error(f"Excel report generation failed: {e}")

    total_runtime = datetime.datetime.now() - start_time
    current_directory = os.getcwd()
    logger.info(f"Analysis completed in directory: {current_directory}")
    logger.info(f"Total pipeline runtime: {total_runtime}")

    logger.debug("Performing final cleanup...")
    cleanup_artifacts(debug=args.debug, logger=logger)

    summary_log.write()
    logger.info(f"Summary log saved as: {summary_log.log_file}")
    write_structured_reports("completed")

    logger.info(f"\n{bcolors.GREEN}Bioinformatics pipeline completed successfully!{bcolors.ENDC}")
    logger.info(f"{bcolors.WHITE}Final runtime: {total_runtime}{bcolors.ENDC}\n")
    


# Created December 2024 by Tod Stuber