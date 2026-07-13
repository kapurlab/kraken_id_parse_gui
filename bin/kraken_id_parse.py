#!/usr/bin/env python

__version__ = "0.0.1"

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
import pandas as pd
from collections import OrderedDict, defaultdict
from random import randint
from time import sleep
from Bio import SeqIO

from file_setup import bcolors, Latex_Report, Excel_Stats, safe_move

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
from reporting import build_run_manifest, render_html_report, render_pdf_report, write_manifest


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
    parser.add_argument("-b", "--blast_db", action='store', dest='blast_db', default="nt", help='Specify BLAST db to use')
    parser.add_argument("-s", "--specific", action='store', dest='specific', default=None, help='Specify custom script/function for the target being used.  Often just default to the taxon name')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('--kraken-only', action='store_true', dest='kraken_only', default=False, help='Run only Kraken2 and produce the Krona graph; skip read parsing, assembly, and BLAST. Requires -k.')
    parser.add_argument('--no-blast', action='store_true', dest='no_blast', default=False, help='Run Kraken2 and taxonomic read parsing, then stop: skip assembly, BLAST, and coverage. Leaves the parsed FASTQ.gz reads for the target taxon. Requires -k and -t.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    # Pre-flight: verify all required external tools are on PATH
    REQUIRED_TOOLS = [
        'kraken2', 'bracken', 'kreport2krona.py', 'ktImportText',
        'spades.py', 'blastn', 'bwa', 'samtools', 'picard',
        'freebayes', 'freebayes-parallel', 'vcffilter', 'pigz', 'seqkit',
        'tectonic',
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
    if missing_tools:
        print(f'\n  WARNING: {len(missing_tools)} tool(s) missing: {", ".join(missing_tools)}')
        print(f'  Pipeline will fail when these are reached.\n')
    else:
        print(f'  All tools found.\n')
    print(f'{"="*55}\n')

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

    if args.kraken_db:
        kraken = Kraken_Identification(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, kraken_db=args.kraken_db, directory='kraken')
        #Bracken Pie
        kraken.run()
        krona_html = kraken.krona_make_graph(kraken.report)
        if args.kraken_only:
            # Kraken-only mode: the Krona graph is the deliverable. Skip Bracken,
            # read parsing, assembly, and BLAST entirely and finish here.
            print(f"\n{bcolors.GREEN}Kraken-only mode: Krona graph generated at {krona_html}.{bcolors.ENDC}")
            print(f"{bcolors.GREEN}Skipping Bracken, read parsing, assembly, and BLAST.{bcolors.ENDC}")
            sys.exit(0)
        # Bracken re-estimates abundances but has no macOS/arm64 conda build, so
        # the binary may be absent. It's optional (kraken-only mode already skips
        # it) — when it's not on PATH, warn and skip Bracken + its pie charts
        # instead of crashing, so read parsing/identification still runs.
        if shutil.which("bracken"):
            kraken.bracken(kraken.report, kraken.output)
            bracken_pie_charts = Bracken_Pie_Charts()
            bracken_pie_charts.run(kraken.bracken_excel)
            bracken_pie_charts.latex(build_latex=latex_report.tex)
        else:
            print(f"{bcolors.YELLOW}WARNING: bracken not found on PATH — skipping "
                  f"Bracken abundance re-estimation and its pie charts "
                  f"(no macOS build). Read parsing and identification continue."
                  f"{bcolors.ENDC}")
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
    def cleanup_artifacts(debug: bool = False):
        temp_dir = './temp'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        files_grab = []
        for files in ('*.aux', '*.log', '*tex', '*out', 'batch.sh', '*_blast_all.txt', '*_blast_out.txt', 'coverage_list.txt', 'coverage.txt'):
            files_grab.extend(glob.glob(files))
        print("Removing:")
        for each in files_grab:
            print(f'\t{each}')
            safe_move(each, temp_dir)

        shutil.rmtree(temp_dir)
        try: 
            shutil.rmtree('alignment_all')
        except: 
            pass
        try:
            shutil.rmtree('alignment_top')
        except: 
            pass

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
            pdf_path, pdf_warning = render_pdf_report(html_path, output_dir)
            if pdf_warning:
                manifest.setdefault("warnings", []).append(pdf_warning)
                manifest_path = write_manifest(manifest, output_dir)
                html_path = render_html_report(manifest, output_dir)
                print(f"WARNING: {pdf_warning}")
            else:
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
                print(f"PDF report generated: {pdf_path}")
            print(f"Structured report manifest generated: {manifest_path}")
            print(f"HTML report generated: {html_path}")
        except Exception as exc:
            print(f"WARNING: structured HTML report generation failed: {exc}")
        
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
        write_structured_reports(
            "completed_with_warnings",
            [f"Target taxon not found by Kraken: {args.taxon}"],
        )
        sys.exit(0)

    # --no-blast: the taxonomically-parsed reads are the deliverable. Stop here
    # before assembly/BLAST/coverage. The gzipped R1/R2 reads for the target
    # taxon (parser.r1_out.gz / parser.r2_out.gz) are left in place for reuse
    # (e.g. re-running them through vSNP).
    if args.no_blast:
        print(f"\n{bcolors.GREEN}--no-blast mode: extracted target reads "
              f"({parser.r1_out}.gz, {parser.r2_out}.gz).{bcolors.ENDC}")
        print(f"{bcolors.GREEN}Skipping assembly, BLAST, and coverage analysis.{bcolors.ENDC}")
        latex_report.latex_ending()
        excel_stats.post_excel()
        cleanup_artifacts()
        write_structured_reports("completed")
        print(f"{bcolors.WHITE}Total runtime: {datetime.datetime.now() - start_time}{bcolors.ENDC}\n")
        sys.exit(0)

    # Step 2: Assembly of filtered reads
    print(f"\n{bcolors.BLUE}Step 2: Running assembly of filtered reads...{bcolors.ENDC}")
    assembler = Assemble(
        FASTQ_R1=f"{parser.r1_out}.gz",
        FASTQ_R2=f"{parser.r2_out}.gz",
        debug=args.debug
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
        write_structured_reports(
            "completed_with_warnings",
            [f"Assembly did not complete for target taxon: {args.taxon}"],
        )
        sys.exit(0)
        
    assembler.stats(assembler.FASTA)
    
    # Step 3: BLAST analysis
    print(f"\n{bcolors.BLUE}Step 3: Running BLAST analysis using {args.blast_db}...{bcolors.ENDC}")
    blast = Blast_Fasta(FASTA=assembler.FASTA, blast_db=args.blast_db)

    blast_to_coverage = BlastCoverageBridge(blast.blast_summary_file)
    blast_to_coverage.process()

    # Step 4: Coverage analysis
    print(f"\n{bcolors.BLUE}Step 4: Running coverage analysis...{bcolors.ENDC}")
    # Get top 5 unique accessions from BLAST results        
    coverage_graph = Coverage_Graph(FASTA=blast_to_coverage.combined_fasta, FASTQ_R1=assembler.FASTQ_R1, FASTQ_R2=assembler.FASTQ_R2, debug=args.debug)
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
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
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
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
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
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
            coverage_graph.get_coverage_graph()
            coverage_graph.latex(latex_report.tex)
            source = open('concatenated_specific.fasta', 'r')
            dest.write(source.read())
            source.close()
            os.remove('concatenated_specific.fasta')

        if ehv_dict:
            orbivirus_specific.run(alignment_stats=ehv_dict)
            coverage_graph = Coverage_Graph(FASTA="concatenated_specific.fasta", FASTQ_R1=parser.FASTQ_R1, FASTQ_R2=parser.FASTQ_R2, debug=args.debug)
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
            download = Downloader(accession)
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
        alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='all_downloaded.fasta', skip_assembly=True, debug=True)
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
                
                download = Downloader(acc)
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
        
        shutil.rmtree('./temp')
        alignment = Alignment(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference='top_downloaded_genomes.fasta', skip_assembly=True, debug=True)
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
        
        coverage_graph = Coverage_Graph(FASTA="output.fasta", FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2)
        coverage_graph.get_coverage_graph()
        os.remove("output.fasta")
        coverage_graph.latex(latex_report.tex)

####################################################################################################################

    latex_report.latex_ending()
    excel_stats.post_excel()
    
    current_directory = os.getcwd()
    print(current_directory)
    
    cleanup_artifacts()
    write_structured_reports("completed")
    


# Created December 2024 by Tod Stuber
