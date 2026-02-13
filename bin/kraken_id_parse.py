#!/usr/bin/env python

__version__ = "0.0.2"

import os
import shutil
import sys
import glob
import re
import datetime
import random
import time
from pathlib import Path
import argparse
import textwrap
import yaml
import pandas as pd
from collections import OrderedDict, defaultdict
from random import randint
from time import sleep
from Bio import SeqIO

from file_setup import bcolors, Banner, Latex_Report, Excel_Stats


# ---------------------------------------------------------------------------
# Taxon → BLAST database mapping for automatic database resolution
# ---------------------------------------------------------------------------
# Categories: viral taxa use nt_viruses, bacterial use ref_prok_rep_genomes,
# protozoan/other eukaryotic use nt.  Fallback is nt (full).
TAXON_DB_MAP = {
    # Viral families / genera
    "Orbivirus":        "nt_viruses",
    "Flaviviridae":     "nt_viruses",
    "Coronaviridae":    "nt_viruses",
    "Retroviridae":     "nt_viruses",
    "Poxviridae":       "nt_viruses",
    "Adenoviridae":     "nt_viruses",
    "Vesiculovirus":    "nt_viruses",
    "Paramyxoviridae":  "nt_viruses",
    "Metapneumovirus":  "nt_viruses",
    "Asfivirus":        "nt_viruses",
    "Arteriviridae":    "nt_viruses",
    "Novirhabdovirus":  "nt_viruses",
    "Isavirus salaris": "nt_viruses",
    "Viruses":          "nt_viruses",
    "Bluetongue virus": "nt_viruses",
    "Epizootic hemorrhagic disease": "nt_viruses",
    # Bacterial
    "Mycobacterium tuberculosis complex": "ref_prok_rep_genomes",
    "Mycobacterium bovis":                "ref_prok_rep_genomes",
    "Brucellaceae":                       "ref_prok_rep_genomes",
    "Leptospirales":                      "ref_prok_rep_genomes",
    # Protozoan / eukaryotic
    "Apicomplexa":      "nt",
}

# Broad category keywords as fallback when exact taxon isn't in the map
TAXON_CATEGORY_HINTS = [
    # (keyword_in_taxon,  db_name)
    ("virus",   "nt_viruses"),
    ("viridae", "nt_viruses"),
    ("virales", "nt_viruses"),
    ("phage",   "nt_viruses"),
]


def resolve_blast_db(taxon: str, database_root: str = None, explicit_db: str = None) -> str:
    """Resolve the BLAST database path from taxon name and optional database_root.

    Priority:
      1. If ``explicit_db`` is a real path (absolute) → use it verbatim.
      2. Look up taxon in TAXON_DB_MAP (exact match, then case-insensitive).
      3. Try keyword hints (viral family suffixes).
      4. Fall back to "nt".

    If ``database_root`` is provided, the resolved DB name is joined as
    ``database_root / db_name``.  Otherwise the bare name is returned
    (user is expected to have it on BLASTDB path or give an absolute path).
    """
    # 1. Explicit absolute path — pass through
    if explicit_db and explicit_db != "nt":
        if os.path.isabs(explicit_db):
            return explicit_db

    # 2. Exact lookup
    db_name = TAXON_DB_MAP.get(taxon)

    # 3. Case-insensitive lookup
    if db_name is None:
        for key, val in TAXON_DB_MAP.items():
            if key.lower() == taxon.lower():
                db_name = val
                break

    # 4. Keyword hints
    if db_name is None:
        taxon_lower = taxon.lower()
        for keyword, val in TAXON_CATEGORY_HINTS:
            if keyword in taxon_lower:
                db_name = val
                break

    # 5. Fallback
    if db_name is None:
        db_name = "nt"

    # Resolve against database_root if provided
    if database_root:
        candidate = os.path.join(database_root, db_name)
        # Check for BLAST DB index files (e.g. nt_viruses.nal or nt_viruses.00.nhr)
        if glob.glob(f"{candidate}*"):
            return candidate
        else:
            print(f"{bcolors.YELLOW}Warning: No BLAST index files found at "
                  f"{candidate}* — falling back to bare name '{db_name}'{bcolors.ENDC}")

    # If an explicit non-absolute value was given (like a preset db name), prefer it
    if explicit_db and explicit_db != "nt":
        return explicit_db

    return db_name


def load_database_root() -> str:
    """Load database_root from ~/.kraken_id_parse.yaml if it exists."""
    config_path = os.path.expanduser("~/.kraken_id_parse.yaml")
    if os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                cfg = yaml.safe_load(f) or {}
            return cfg.get('database_root', '')
        except Exception as e:
            print(f"Warning: Could not read {config_path}: {e}")
    return ''


def load_reference_cache() -> str:
    """Load reference_cache directory from ~/.kraken_id_parse.yaml if it exists."""
    config_path = os.path.expanduser("~/.kraken_id_parse.yaml")
    if os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                cfg = yaml.safe_load(f) or {}
            return cfg.get('reference_cache', '')
        except Exception as e:
            print(f"Warning: Could not read {config_path}: {e}")
    return ''

import gzip


def detect_platform(fastq_path: str, n_reads: int = 100) -> str:
    """Auto-detect sequencing platform from FASTQ read lengths.

    Reads the first *n_reads* records. If the max observed read length
    exceeds 701 bp the data is assumed to be Oxford Nanopore (ONT);
    otherwise Illumina short-read.  Threshold follows vSNP3 convention.
    """
    opener = gzip.open if fastq_path.endswith('.gz') else open
    max_len = 0
    count = 0
    try:
        with opener(fastq_path, 'rt') as fh:
            while count < n_reads:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()  # +
                fh.readline()  # qual
                if seq:
                    max_len = max(max_len, len(seq))
                    count += 1
    except Exception as e:
        print(f"Warning: could not auto-detect platform from {fastq_path}: {e}")
        return 'illumina'
    platform = 'ont' if max_len > 701 else 'illumina'
    print(f"Platform detected: {platform} (max read length: {max_len:,} bp from first {count} reads)")
    return platform


from fastq_stats_seqkit import FASTQ_Stats
from alignment_vcf import Alignment
from blast_fasta_and_search import Blast_Fasta
from fasta_gbk_gff_by_acc import Downloader
from kraken import Kraken_Identification
from kraken import Bracken_Pie_Charts
from reference_guided_assembly_vcf_to_fasta import Reference_Guided_Assembly
from parse_reads import ParseReads, TaxonNotFoundError
from assembly import Assemble, SPAdesDidNotAssembleFASTA
# from blast_hpc import Blast
from coverage_graph import Coverage_Graph
from blast_to_coverage import BlastCoverageBridge
from orbivirus_specific import Orbivirus_Specific
from isav_specific import ISAV_Specific
from apicomplexa_specific import Apicomplexa


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
    parser.add_argument('-l', '--logo', action='store', dest='logo', required=False, help='Logo for the Latex PDF report')
    parser.add_argument("-t", "--taxon", action='store', dest='taxon', help='Target Taxon')
    parser.add_argument("-k", "--kraken_db", action='store', dest='kraken_db', help='Specify Kraken db to use')
    parser.add_argument("-b", "--blast_db", action='store', dest='blast_db', default="nt", help='Specify BLAST db to use (overrides auto-resolution)')
    parser.add_argument("--database-root", action='store', dest='database_root', default=None,
                        help='Root directory containing BLAST databases (enables auto-resolution from taxon name)')
    parser.add_argument("-s", "--specific", action='store', dest='specific', default=None, help='Specify custom script/function for the target being used.  Often just default to the taxon name')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('--keep-extracted-reads', action='store_true', dest='keep_extracted_reads', default=False, help='Keep taxon-extracted FASTQ files (default: remove to save space)')
    parser.add_argument('--reference-cache', action='store', dest='reference_cache', default=None,
                        help='Directory to cache downloaded NCBI reference FASTAs (skips re-download on cache hit)')
    parser.add_argument('--platform', action='store', dest='platform', default='auto',
                        choices=['illumina', 'ont', 'auto'],
                        help='Sequencing platform: illumina (BWA+SPAdes), ont (minimap2+Flye), auto (detect from read length >701bp)')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    # --- Resolve reference cache directory ---
    reference_cache = args.reference_cache or load_reference_cache() or None
    if reference_cache:
        print(f"{bcolors.BLUE}Reference cache directory: {reference_cache}{bcolors.ENDC}")

    # --- Resolve sequencing platform ---
    if args.platform == 'auto':
        platform = detect_platform(args.FASTQ_R1)
    else:
        platform = args.platform
    print(f"{bcolors.BLUE}Sequencing platform: {platform}{bcolors.ENDC}")

    # --- Auto-resolve BLAST database from taxon ---
    database_root = args.database_root or load_database_root()
    user_gave_explicit_blast = (args.blast_db != "nt")  # user explicitly set -b
    if args.taxon and not user_gave_explicit_blast:
        resolved = resolve_blast_db(
            taxon=args.taxon,
            database_root=database_root if database_root else None,
            explicit_db=args.blast_db
        )
        if resolved != args.blast_db:
            print(f"{bcolors.BLUE}Auto-resolved BLAST database: {resolved} "
                  f"(from taxon '{args.taxon}'){bcolors.ENDC}")
            args.blast_db = resolved

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    # id_parse = ID_Parse(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2,)
    sample_name = re.sub('[._].*', '', os.path.basename(args.FASTQ_R1))

    #Latex report
    latex_report = Latex_Report(sample_name=sample_name, logo=args.logo)
    #Excel Stats
    excel_stats = Excel_Stats(sample_name)

    fastq_stats = FASTQ_Stats(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=True)
    fastq_stats.run()
    fastq_stats.latex(latex_report.tex)
    fastq_stats.excel(excel_stats.excel_dict)

    if args.kraken_db:
        kraken = Kraken_Identification(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, kraken_db=args.kraken_db, directory='kraken')
        #Bracken Pie
        kraken.run()
        krona_html = kraken.krona_make_graph(kraken.report)
        kraken.bracken(kraken.report, kraken.output)
        bracken_pie_charts = Bracken_Pie_Charts()
        bracken_pie_charts.run(kraken.bracken_excel)
        bracken_pie_charts.latex(build_latex=latex_report.tex)
    else:
        class Kraken:
            def __init__(self):
                self.output = None
                self.report = None

        # Create kraken instance
        kraken = Kraken()

        # Check if either kraken or kraken2 folders exist
        kraken_folder = None
        if os.path.isdir("kraken"):
            kraken_folder = "kraken"
        elif os.path.isdir("kraken2"):
            kraken_folder = "kraken2"

        if kraken_folder:
            # Find the files using glob
            output_files = glob.glob(os.path.join(kraken_folder, "*_outputkraken.txt"))
            report_files = glob.glob(os.path.join(kraken_folder, "*_reportkraken.txt"))
            
            # Check if files were found and assign to kraken object
            if output_files and report_files:
                kraken.output = output_files[0]  # Get the first matching output file
                kraken.report = report_files[0]  # Get the first matching report file
                print(f"Found files:\nOutput: {kraken.output}\nReport: {kraken.report}")
            else:
                print("Required files not found in the kraken folder")
        else:
            print("Neither kraken nor kraken2 folder found in current directory.  Run script with -k option to first run Kraken.")

############################
    def cleanup_artifacts(keep_extracted_reads: bool = False):
        """Clean up intermediate files while preserving final outputs.

        Preserved (always):
            *_report.pdf          - PDF report
            *_stats.xlsx          - Excel statistics
            *_denovo.fasta        - De novo assembly
            *_reference_guided.fasta - Consensus sequence
            *_blast_summary.txt   - BLAST summaries
            consensus_blast_summary.txt
            CAUTION_SITES.xlsx    - Ambiguous sites QC
            kraken/ folder        - Kraken report, bracken xlsx, krona HTML
                                   (but NOT the large _outputkraken.txt)
            Input FASTQ symlinks/files

        Preserved (with --keep-extracted-reads):
            *_<taxon>_R1.fastq.gz - Extracted taxon reads
            *_<taxon>_R2.fastq.gz

        Removed:
            *.aux, *.log, *.out   - LaTeX intermediates
            *_report.tex          - LaTeX source
            *-banner.png          - Section banners (embedded in PDF)
            bracken_pie.png       - Pie chart (embedded in PDF)
            *-coverage_graph.pdf  - Coverage graphs (embedded in PDF)
            combined_references.fasta - Downloaded references
            *_blast_all.txt, *_blast_out.txt - Raw BLAST output
            coverage_list.txt, coverage.txt - Coverage intermediates
            batch.sh              - Temp batch script
            kraken/_outputkraken.txt - Large per-read classification (32MB+)
        """
        # --- Files to always remove ---
        temp_patterns = [
            '*.aux', '*.log', '*.out', 'batch.sh',
            '*_report.tex',
            '*_blast_all.txt', '*_blast_out.txt',
            '*-banner.png',
            'bracken_pie.png',
            '*-coverage_graph.pdf',
            'combined_references.fasta',
            'coverage_list.txt', 'coverage.txt',
        ]

        files_to_remove = []
        for pattern in temp_patterns:
            files_to_remove.extend(glob.glob(pattern))

        # Remove large kraken output file (per-read classifications)
        for f in glob.glob('kraken/*_outputkraken.txt'):
            files_to_remove.append(f)

        # Optionally remove extracted taxon reads (default: remove)
        if not keep_extracted_reads:
            for f in glob.glob('*_R1.fastq.gz') + glob.glob('*_R2.fastq.gz'):
                # Only remove extracted reads (contain taxon name), not original inputs
                # Extracted reads have pattern: SAMPLE_Taxon_R1.fastq.gz
                basename = os.path.basename(f)
                parts = basename.replace('.fastq.gz', '').rsplit('_', 1)  # split off R1/R2
                if len(parts) == 2:
                    prefix = parts[0]  # e.g. SRR9598511_Orbivirus
                    # If prefix contains an underscore, it's likely an extracted read
                    # (original inputs are SAMPLE_R1.fastq.gz, extracted are SAMPLE_Taxon_R1.fastq.gz)
                    if prefix.count('_') >= 1:
                        # Check it's not the original input by seeing if a simpler name exists
                        sample_base = prefix.split('_')[0]
                        r_suffix = parts[1]  # R1 or R2
                        original = f'{sample_base}_{r_suffix}.fastq.gz'
                        if os.path.exists(original) and f != original and os.path.abspath(f) != os.path.abspath(original):
                            files_to_remove.append(f)

        if files_to_remove:
            print("Cleaning up intermediate files:")
            for f in files_to_remove:
                print(f'\t{f}')
                try:
                    os.remove(f)
                except OSError as e:
                    print(f'\t  Warning: could not remove {f}: {e}')

        # Remove temporary alignment directories
        for d in ['alignment_all', 'alignment_top']:
            if os.path.isdir(d):
                try:
                    shutil.rmtree(d)
                except OSError:
                    pass

        # Remove ./temp directory if it exists from earlier steps
        if os.path.isdir('./temp'):
            shutil.rmtree('./temp')

        # Print summary of preserved files
        preserved = [f for f in os.listdir('.') if os.path.isfile(f)]
        preserved_dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
        print(f"\nFinal output files ({len(preserved)} files, {len(preserved_dirs)} directories):")
        for f in sorted(preserved):
            size = os.path.getsize(f)
            if size > 1024*1024:
                print(f'\t{f}  ({size/(1024*1024):.1f} MB)')
            elif size > 1024:
                print(f'\t{f}  ({size/1024:.1f} KB)')
            else:
                print(f'\t{f}  ({size} B)')
        for d in sorted(preserved_dirs):
            print(f'\t{d}/')
        
    """Execute the complete analysis pipeline"""
    print(f"\n{bcolors.GREEN}Starting bioinformatics analysis pipeline...{bcolors.ENDC}\n")
    
    # Step 1: Parse reads using Kraken results
    print(f"{bcolors.BLUE}Step 1: Running taxonomic read parsing...{bcolors.ENDC}")
    parser = ParseReads(
        R1=args.FASTQ_R1,
        R2=args.FASTQ_R2,
        kraken_output=kraken.output,
        kraken_report=kraken.report,
        taxon=args.taxon,
        output_prefix=args.taxon,
        debug=args.debug
    )
    try:
        parser.run()
        parser.latex(latex_report.tex)
        parser.excel(excel_stats.excel_dict)
    except TaxonNotFoundError:
        print(f'\\textcolor{{red}}{{\\textbf{{\\Large Target {args.taxon} taxon not found by Kraken}}}}\\\\[1em]\n\nScript terminated', file=latex_report.tex)
        latex_report.latex_ending()
        excel_stats.excel_dict['Target Taxon'] = args.taxon
        excel_stats.excel_dict['Extracted Reads'] = "No Reads Found"
        excel_stats.post_excel()
        cleanup_artifacts()
        sys.exit(0)
        
    
    # Step 2: Assembly of filtered reads
    print(f"\n{bcolors.BLUE}Step 2: Running assembly of filtered reads...{bcolors.ENDC}")
    assembler = Assemble(
        FASTQ_R1=f"{parser.r1_out}.gz",
        FASTQ_R2=f"{parser.r2_out}.gz",
        debug=args.debug,
        platform=platform
    )
    try:
        assembler.run()
    except SPAdesDidNotAssembleFASTA:
        print('Assembly failed, Check reads')
        print(f'\\textcolor{{red}}{{\\textbf{{\\Large Assembly failed, Check reads - Likely not enough {args.taxon} reads to complete assembly}}}}\\\\[1em]\n\nScript terminated', file=latex_report.tex)
        latex_report.latex_ending()
        excel_stats.excel_dict['Contig count'] = "Assembly did not complete"
        excel_stats.post_excel()
        cleanup_artifacts()
        sys.exit(0)
        
    assembler.stats(assembler.FASTA)
    
    # Step 3: BLAST analysis
    print(f"\n{bcolors.BLUE}Step 3: Running BLAST analysis using {args.blast_db}...{bcolors.ENDC}")
    blast = Blast_Fasta(FASTA=assembler.FASTA, blast_db=args.blast_db)

    blast_to_coverage = BlastCoverageBridge(blast.blast_summary_file, cache_dir=reference_cache)
    blast_to_coverage.process()

    # Step 4: Coverage analysis
    print(f"\n{bcolors.BLUE}Step 4: Running coverage analysis...{bcolors.ENDC}")
    # Get top 5 unique accessions from BLAST results        
    coverage_graph = Coverage_Graph(FASTA=blast_to_coverage.combined_fasta, FASTQ_R1=assembler.FASTQ_R1, FASTQ_R2=assembler.FASTQ_R2, debug=args.debug, platform=platform)
    coverage_graph.get_coverage_graph()

    # Update the assembled FASTA file with file name containing "denovo"
    path, filename = os.path.split(assembler.FASTA)
    name, ext = os.path.splitext(filename)
    new_filename = f"{name}_denovo{ext}"
    new_FASTA = os.path.join(path, new_filename)
    os.rename(assembler.FASTA, new_FASTA)
    assembler.FASTA = new_FASTA
    print(f"Updated FASTA variable: {assembler.FASTA}")

    # latex_report = Latex_Report(assembler.sample_name)
    assembler.latex(latex_report.tex)
    blast.latex(latex_report.tex)
    if args.taxon == 'Apicomplexa':
        args.specific = args.taxon
        apicomplexa = Apicomplexa()
        isav_dict = {}
        dest = open('all_downloaded.fasta', 'a')
        for seq_id, seq_info in coverage_graph.alignment_stats.items():
            header = seq_info['header'].lower()
            if any(x in header for x in ["Theileria equi"]):
                isav_dict[seq_id] = seq_info
        if isav_dict:
            apicomplexa.run(alignment_stats=isav_dict)
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug, platform=platform)
            coverage_graph.get_coverage_graph()
            coverage_graph.latex(latex_report.tex)
            source = open('concatenated_specific.fasta', 'r')
            dest.write(source.read())
            source.close()
            os.remove('concatenated_specific.fasta')
        dest.close()
    # if a specific test, Orbivirus, is an option
    elif args.taxon == 'Isavirus salaris':
        args.specific = args.taxon
        isav_specific = ISAV_Specific()
        isav_dict = {}
        dest = open('all_downloaded.fasta', 'a')
        for seq_id, seq_info in coverage_graph.alignment_stats.items():
            header = seq_info['header'].lower()
            if any(x in header for x in ["infectious salmon anemia virus", "Isavirus", "ISAV"]):
                isav_dict[seq_id] = seq_info
        if isav_dict:
            isav_specific.run(alignment_stats=isav_dict)
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug, platform=platform)
            coverage_graph.get_coverage_graph()
            coverage_graph.latex(latex_report.tex)
            source = open('concatenated_specific.fasta', 'r')
            dest.write(source.read())
            source.close()
            os.remove('concatenated_specific.fasta')
        dest.close()
    # if a specific test, Orbivirus, is an option
    elif args.taxon == 'Orbivirus':
        args.specific = args.taxon
        orbivirus_specific = Orbivirus_Specific()
        btv_dict = {}
        ehv_dict = {}
        dest = open('all_downloaded.fasta', 'a')
        for seq_id, seq_info in coverage_graph.alignment_stats.items():
            header = seq_info['header'].lower()
            if "bluetongue virus" in header:
                btv_dict[seq_id] = seq_info
            elif "epizootic hemorrhagic disease" in header:
                ehv_dict[seq_id] = seq_info
        if btv_dict:
            orbivirus_specific.run(alignment_stats=btv_dict)
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug, platform=platform)
            coverage_graph.get_coverage_graph()
            coverage_graph.latex(latex_report.tex)
            source = open('concatenated_specific.fasta', 'r')
            dest.write(source.read())
            source.close()
            os.remove('concatenated_specific.fasta')

        if ehv_dict:
            orbivirus_specific.run(alignment_stats=ehv_dict)
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug, platform=platform)
            coverage_graph.get_coverage_graph()
            coverage_graph.latex(latex_report.tex)
            source = open('concatenated_specific.fasta', 'r')
            dest.write(source.read())
            source.close()
            os.remove('concatenated_specific.fasta')
        dest.close()
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
        try:
            accession = get_last_accession(blast.blast_summary_file)
            print(f'Downloading reference: {accession}')
            
            # Download with proper error handling and rate limiting
            download = Downloader(accession, cache_dir=reference_cache)
            max_retries = 3
            success = False
            
            for attempt in range(max_retries):
                try:
                    download.fasta()
                    success = True
                    print(f"✓ Successfully downloaded {accession}")
                    break
                except Exception as e:
                    print(f"⚠ Attempt {attempt + 1} failed for {accession}: {e}")
                    if attempt < max_retries - 1:
                        wait_time = (attempt + 1) * 5  # 5, 10, 15 seconds
                        print(f"  Waiting {wait_time} seconds before retry...")
                        sleep(wait_time)
            
            if success:
                # Rename the downloaded file
                if os.path.exists(f'{download.accession}.fasta'):
                    os.rename(f'{download.accession}.fasta', "all_downloaded.fasta")
                    print(f"Renamed {download.accession}.fasta to all_downloaded.fasta")
                else:
                    raise FileNotFoundError(f"Downloaded file {download.accession}.fasta not found")
            else:
                raise Exception(f"Failed to download {accession} after {max_retries} attempts")
                
        except Exception as e:
            print(f"Error in reference download process: {e}")
            print("Creating empty reference file to prevent pipeline failure")
            with open("all_downloaded.fasta", 'w') as f:
                f.write(f">No_Reference_Available\nNNNNNNNNNN\n")
        
        args.specific = args.taxon
        
    print(f"\n{bcolors.BLUE}Generating analysis reports...{bcolors.ENDC}")
    # latex_report.latex_ending()

    # excel_stats = Excel_Stats(assembler.sample_name)
    assembler.excel(excel_stats.excel_dict)
    blast.excel(excel_stats.excel_dict)
    # excel_stats.post_excel()

    # Clean up temporary directory if it exists
    if not args.debug:
        temp_dir = Path('./temp')
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    print(f"\n{bcolors.GREEN}Pipeline completed successfully!{bcolors.ENDC}")
    print(f"{bcolors.WHITE}Total runtime: {datetime.datetime.now() - start_time}{bcolors.ENDC}\n")

####################################################################################################################
    if args.specific:
        alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='all_downloaded.fasta', skip_assembly=True, debug=True, platform=platform)
        alignment.run()

        # id_parse.parse_fastas(bam=alignment.nodup_bamfile, fasta=alignment.reference)

        reference_guided_assembly = Reference_Guided_Assembly(FASTA=f'{alignment.reference}', vcf=f'{alignment.zc_vcf}', iupac=True, output_name='consensus.fasta')

        blast = Blast_Fasta(FASTA='consensus.fasta', num_alignment=1, blast_db=args.blast_db)

        # IMPROVED DOWNLOAD SECTION - REPLACE THE OLD LOOP WITH THIS:
        print("Downloading top hit reference genomes...")
        print(f"Found {len(blast.top_hit_acc)} accessions to download")
        
        downloaded_count = 0
        failed_downloads = []
        
        with open('top_downloaded_genomes.fasta', 'w') as outfile:
            for i, acc in enumerate(blast.top_hit_acc, 1):
                print(f"Downloading {i}/{len(blast.top_hit_acc)}: {acc}")
                
                # PROACTIVE rate limiting - wait BEFORE each request (except first)
                if i > 1:
                    wait_time = 3  # Fixed 3-second delay between downloads
                    print(f"  Waiting {wait_time} seconds (NCBI rate limiting)...")
                    time.sleep(wait_time)
                
                download = Downloader(acc, cache_dir=reference_cache)
                max_retries = 3
                success = False

                # Retry with exponential backoff
                for attempt in range(max_retries):
                    try:
                        download.fasta()
                        success = True
                        downloaded_count += 1
                        print(f"  ✓ Successfully downloaded {acc}")
                        break
                    except Exception as e:
                        print(f"  ⚠ Attempt {attempt + 1}/{max_retries} failed for {acc}: {e}")
                        if attempt < max_retries - 1:
                            # Exponential backoff: 10, 20, 30 seconds
                            backoff_time = (attempt + 1) * 10
                            print(f"    Waiting {backoff_time} seconds before retry...")
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
                            print(f"  ⚠ Warning: {fasta_file} not found after download")
                            success = False
                    except Exception as e:
                        print(f"  ⚠ Error processing {fasta_file}: {e}")
                        success = False
                
                if not success:
                    failed_downloads.append(acc)
                    print(f"  ✗ Failed to download {acc} after {max_retries} attempts")
        
        # Summary and error handling
        print(f"\nDownload Summary:")
        print(f"  Successfully downloaded: {downloaded_count}/{len(blast.top_hit_acc)}")
        
        if failed_downloads:
            print(f"  Failed downloads: {len(failed_downloads)}")
            print(f"  Failed accessions: {', '.join(failed_downloads)}")
        
        # Check if we have any successful downloads
        if downloaded_count == 0:
            print("ERROR: No reference genomes could be downloaded!")
            print("Creating minimal reference file to prevent pipeline crash...")
            
            # Create a minimal reference to prevent downstream errors
            with open('top_downloaded_genomes.fasta', 'w') as f:
                f.write(">No_Reference_Available\n")
                f.write("N" * 1000 + "\n")  # 1000 N's as placeholder
        else:
            print(f"Proceeding with {downloaded_count} successfully downloaded references...")
        # END OF IMPROVED DOWNLOAD SECTION

        if os.path.isdir('./temp'):
            shutil.rmtree('./temp')
        alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='top_downloaded_genomes.fasta', skip_assembly=True, debug=True, platform=platform)
        alignment.run()

        reference_guided_assembly = Reference_Guided_Assembly(FASTA=f'{alignment.reference}', vcf=f'{alignment.zc_vcf}', iupac=True, output_name='consensus.fasta')
        
        final_consensus = f'{assembler.sample_name}_reference_guided.fasta'
        renamed_fastas=[]
        record_dict = SeqIO.to_dict(SeqIO.parse(f'consensus.fasta', "fasta"))
        for acc, seq in record_dict.items():
            seq.name = f'{assembler.sample_name}'
            seq.id = f'{assembler.sample_name}'
            seq.description = f'guided by {seq.description}'
            renamed_fastas.append(seq)
        SeqIO.write(renamed_fastas, final_consensus, "fasta")
        os.remove('consensus.fasta')

        if os.path.isdir('./temp'):
            shutil.rmtree('./temp')
        
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
        
        coverage_graph = Coverage_Graph(FASTA="output.fasta", FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, platform=platform)
        coverage_graph.get_coverage_graph()
        os.remove("output.fasta")
        coverage_graph.latex(latex_report.tex)

####################################################################################################################

    # Write methodology appendix before closing the report
    def write_report_appendix(tex, taxon, kraken_db, blast_db, sample_name,
                              parser_obj=None, assembler_obj=None, blast_obj=None,
                              platform='illumina'):
        """Write a hybrid methodology appendix to the LaTeX report.

        Dynamically populates method descriptions with actual parameters
        used during this pipeline run.
        """
        # Escape underscores for LaTeX
        def esc(s):
            return str(s).replace('_', r'\_').replace('&', r'\&').replace('%', r'\%')

        # Platform-aware tool names
        assembler_name = 'Flye (ONT long-read)' if platform == 'ont' else 'SPAdes'
        aligner_name = 'minimap2 (ONT)' if platform == 'ont' else 'BWA-MEM'

        kraken_db_name = esc(os.path.basename(kraken_db)) if kraken_db else 'default'
        blast_db_name = esc(os.path.basename(blast_db)) if blast_db else 'nt'
        taxon_esc = esc(taxon) if taxon else 'target organism'

        appendix_banner = Banner("Report Description")

        print(r'\clearpage', file=tex)
        print(r'\newpage', file=tex)
        print(r'\begin{table}[H]', file=tex)
        print(r'\centering', file=tex)
        print(r'\includegraphics[width=\textwidth]{' + appendix_banner.banner + r'}', file=tex)
        print(r'\end{table}', file=tex)
        print(r'\vspace{0.3cm}', file=tex)

        # Pipeline summary
        print(r'\noindent\textbf{\large Pipeline Summary}\\[0.5em]', file=tex)
        print(r'\noindent Input paired-end FASTQ reads were processed through an automated bioinformatics pipeline. '
              f'Reads were classified, filtered for \\textit{{{taxon_esc}}}, assembled \\textit{{de novo}}, '
              f'identified by BLAST, aligned to reference sequences, and analyzed for coverage depth. '
              r'A reference-guided consensus assembly was also generated.\\[1em]', file=tex)

        # Test summary table
        print(r'\noindent\textbf{Analyses Performed}\\[0.3em]', file=tex)
        print(r'\begin{table}[H]', file=tex)
        print(r'\begin{adjustbox}{width=0.85\textwidth}', file=tex)
        print(r'\begin{tabular}{l|l}', file=tex)
        print(r'Analysis & Description \\', file=tex)
        print(r'\hline', file=tex)
        print(r'FASTQ Quality & Read quality metrics (Q30, mean Phred score, read length) \\', file=tex)
        print(f'Taxonomic Classification & Kraken2 / Bracken using \\texttt{{{kraken_db_name}}} database \\\\', file=tex)
        print(f'Read Extraction & Reads classified as \\textit{{{taxon_esc}}} extracted for downstream analysis \\\\', file=tex)
        print(f'\\textit{{De novo}} Assembly & {esc(assembler_name)} assembler \\\\', file=tex)
        print(f'BLAST Identification & BLASTn against \\texttt{{{blast_db_name}}} database \\\\', file=tex)
        print(f'Coverage Analysis & {esc(aligner_name)} alignment to top BLAST reference hits \\\\', file=tex)
        print(r'Reference-Guided Assembly & Consensus sequence from guided alignment \\', file=tex)
        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)

        # Detailed method descriptions
        print(r'\vspace{0.5cm}', file=tex)
        print(r'\noindent\textbf{\large Method Descriptions}\\[0.5em]', file=tex)

        # FASTQ Quality
        print(r'\noindent\textbf{FASTQ Quality}\\', file=tex)
        print(r'\noindent Compressed FASTQ files were evaluated for read quality. '
              r'\textbf{Q30 Passing}: percentage of reads with average Phred score $>$30. '
              r'\textbf{Mean Read Score}: average Phred score across all base calls. '
              r'\textbf{Average Read Length}: mean read length in base pairs.\\[0.8em]', file=tex)

        # Kraken/Bracken
        print(r'\noindent\textbf{Taxonomic Classification (Kraken2/Bracken)}\\', file=tex)
        print(r'\noindent Reads were classified using '
              r'\href{https://ccb.jhu.edu/software/kraken2/}{Kraken2} with the '
              f'\\texttt{{{kraken_db_name}}} database. '
              r'Species-level abundance was estimated using '
              r'\href{https://ccb.jhu.edu/software/bracken/}{Bracken}. '
              r'Species comprising $>$1\% of total reads are shown in the pie chart. '
              r'Unclassified reads had no exact k-mer matches in the database.\\[0.8em]', file=tex)

        # Assembly
        print(r'\noindent\textbf{\textit{De novo} Assembly}\\', file=tex)
        if platform == 'ont':
            print(r'\noindent Extracted long reads were assembled using '
                  r'\href{https://github.com/fenderglass/Flye}{Flye} (\texttt{-{}-nano-raw}). '
                  r'\textbf{N50}: half of the total assembled length is from contigs $\geq$ this value. '
                  r'\textbf{Mean Coverage}: estimated average read depth across total assembly length '
                  r'(read count $\times$ read length $\div$ total assembly length).\\[0.8em]', file=tex)
        else:
            print(r'\noindent Extracted reads were assembled using '
                  r'\href{https://github.com/ablab/spades}{SPAdes}. '
                  r'\textbf{N50}: half of the total assembled length is from contigs $\geq$ this value. '
                  r'\textbf{Mean Coverage}: estimated average read depth across total assembly length '
                  r'(read count $\times$ read length $\div$ total assembly length).\\[0.8em]', file=tex)

        # BLAST
        print(r'\noindent\textbf{BLAST Identification}\\', file=tex)
        print(r'\noindent Assembled contigs were searched using '
              r'\href{https://blast.ncbi.nlm.nih.gov/Blast.cgi}{BLASTn} against the '
              f'\\texttt{{{blast_db_name}}} database. '
              r'Results are sorted by total nucleotide representation (summed contig lengths for each top hit). '
              r'Top unique accessions were used as references for coverage analysis.\\[0.8em]', file=tex)

        # Coverage
        print(r'\noindent\textbf{Coverage Analysis}\\', file=tex)
        if platform == 'ont':
            print(r'\noindent Reads were aligned to reference sequences using '
                  r'\href{https://github.com/lh3/minimap2}{minimap2} (\texttt{-x map-ont}). '
                  r'Duplicates were removed with Picard MarkDuplicates. '
                  r'Depth-of-coverage was calculated using samtools. '
                  r'VCF QUAL scores were adjusted +100 to compensate for lower nanopore base quality. '
                  r'When average depth exceeds 100X, log-scale depth is plotted; otherwise linear scale is used. '
                  r'The red dashed line indicates the 100X threshold. '
                  r'\textbf{\% Genome Covered}: percentage of the reference with $\geq$1X read depth.\\[0.8em]', file=tex)
        else:
            print(r'\noindent Reads were aligned to reference sequences using '
                  r'\href{https://github.com/lh3/bwa}{BWA-MEM}. '
                  r'Duplicates were removed with Picard MarkDuplicates. '
                  r'Depth-of-coverage was calculated using samtools. '
                  r'When average depth exceeds 100X, log-scale depth is plotted; otherwise linear scale is used. '
                  r'The red dashed line indicates the 100X threshold. '
                  r'\textbf{\% Genome Covered}: percentage of the reference with $\geq$1X read depth.\\[0.8em]', file=tex)

        # Reference-guided assembly
        print(r'\noindent\textbf{Reference-Guided Consensus}\\', file=tex)
        if platform == 'ont':
            print(r'\noindent A reference-guided consensus assembly was generated by aligning reads to the '
                  r'top BLAST reference hits using minimap2, calling variants with freebayes, and producing '
                  r'a consensus FASTA. This consensus was then re-evaluated for coverage depth.\\[0.8em]', file=tex)
        else:
            print(r'\noindent A reference-guided consensus assembly was generated by aligning reads to the '
                  r'top BLAST reference hits using BWA-MEM, calling variants with freebayes, and producing '
                  r'a consensus FASTA. This consensus was then re-evaluated for coverage depth.\\[0.8em]', file=tex)

    write_report_appendix(
        tex=latex_report.tex,
        taxon=args.taxon,
        kraken_db=args.kraken_db,
        blast_db=args.blast_db,
        sample_name=sample_name,
        parser_obj=parser,
        assembler_obj=assembler,
        blast_obj=blast,
        platform=platform,
    )

    latex_report.latex_ending()
    pdf_file = latex_report.tex_file.replace('.tex', '.pdf')
    if os.path.exists(pdf_file):
        print(f"\n{bcolors.GREEN}Report generated: {pdf_file}{bcolors.ENDC}")
    else:
        print(f"\n{bcolors.RED}Warning: PDF report was not generated. Check pdflatex installation.{bcolors.ENDC}")

    excel_stats.post_excel()

    current_directory = os.getcwd()
    print(current_directory)

    cleanup_artifacts(keep_extracted_reads=args.keep_extracted_reads)



# Created December 2024 by Tod Stuber