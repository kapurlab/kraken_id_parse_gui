#!/usr/bin/env python3

__version__ = "0.0.1"

import argparse
import textwrap
import sys
import os
from pathlib import Path
import subprocess
import re
import gzip
import pandas as pd
from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats

__all__ = ['TaxonNotFoundError', 'ParseReads']

def get_name_indent(line: str) -> int:
    """Get the number of spaces before the taxonomic name"""
    parts = line.split('\t')
    if len(parts) >= 6:
        return len(parts[5]) - len(parts[5].lstrip())
    parts = re.split(r'\s+', line.rstrip())
    if len(parts) >= 6:
        match = re.search(r'(\s+)\S+$', line)
        if match:
            return len(match.group(1))
    return 0

def get_taxonomy_info(report_file: str, target_taxon: str, debug: bool = False) -> tuple:
    """
    Extract taxonomy ID and relevant lines for the target taxon from Kraken report.
    Returns (taxid, list of relevant report lines)
    """
    target_taxon = target_taxon.lower()
    taxid = None
    relevant_lines = []
    
    if debug:
        print("\nDEBUG: Starting taxonomy search...")
    
    with open(report_file) as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if not line.strip():
            continue
            
        parts = re.split(r'\s+', line.strip())
        if len(parts) >= 6:
            name = ' '.join(parts[5:]).strip().lower()
            if name == target_taxon:
                taxid = parts[4]
                target_indent = get_name_indent(line)
                relevant_lines.append(line)
                
                if debug:
                    print(f"\nDEBUG: Found target line: {line.strip()}")
                    print(f"DEBUG: Target name indent level: {target_indent}")
                
                for next_line in lines[i+1:]:
                    if not next_line.strip():
                        continue
                        
                    next_indent = get_name_indent(next_line)
                    
                    if debug:
                        print(f"DEBUG: Checking line: {next_line.strip()}")
                        print(f"DEBUG: Name indent level: {next_indent}")
                        print(f"DEBUG: Include? {next_indent > target_indent}")
                        
                    if next_indent > target_indent:
                        relevant_lines.append(next_line)
                    else:
                        break
                break
    
    return taxid, relevant_lines

def count_fastq_reads(filename: str) -> int:
    """Count reads in a FASTQ file"""
    count = 0
    opener = gzip.open if filename.endswith('.gz') else open
    with opener(filename, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count

def gzip_file(input_file: str, remove_original: bool = True):
    """Gzip a file"""
    with open(input_file, 'rb') as f_in:
        with gzip.open(f"{input_file}.gz", 'wb') as f_out:
            f_out.write(f_in.read())
    if remove_original:
        Path(input_file).unlink()

class TaxonNotFoundError(Exception):
    """Custom exception for when a taxon ID cannot be found."""
    pass

class ParseReads(Setup):
    def __init__(self, R1=None, R2=None, kraken_output=None, kraken_report=None, 
                 taxon=None, output_prefix=None, debug=False):
        # Call parent class (Setup) initialization first
        super().__init__(FASTQ_R1=R1, FASTQ_R2=R2, debug=debug)
        
        # Store R1 and R2 paths directly
        self.R1 = R1
        self.R2 = R2
        
        # Extract sample name from R1
        self.sample_name = re.sub('[_.].*', '', os.path.basename(R1))
            
        # Initialize ParseReads-specific attributes
        self.kraken_output = kraken_output
        self.kraken_report = kraken_report
        self.taxon = taxon
        self.output_prefix = output_prefix.replace(' ', '_')
        self.metrics = {}
        
        if self.debug:
            print(f"Debug: ParseReads initialized with:")
            print(f"  Sample name: {self.sample_name}")
            print(f"  R1: {self.R1}")
            print(f"  R2: {self.R2}")

    def run(self):
        # Get taxonomy information
        self.taxid, self.relevant_report_lines = get_taxonomy_info(
            self.kraken_report, self.taxon, self.debug)
        
        if not self.taxid:
            print(f"Error: No taxid found for {self.taxon}", file=sys.stderr)
            # Create the file to indicate the Kraken did not find taxon
            notification_file =f'KRAKEN-NO_TAXID_FOUND_FOR_{self.taxon}'  
            with open(notification_file, 'w') as f:
                f.write("Script is still running.")
            print(f"Created file: {notification_file}")
            raise TaxonNotFoundError(f"No taxid found for {self.taxon}")
            
        # Save relevant report lines
        self.report_output = f"{self.sample_name}_{self.output_prefix}_kraken.report"
        
        # Setup output filenames
        self.r1_out = f"{self.sample_name}_{self.output_prefix}_R1.fastq"
        self.r2_out = f"{self.sample_name}_{self.output_prefix}_R2.fastq"
        
        # Run KrakenTools
        cmd = [
            "extract_kraken_reads.py",
            "-k", self.kraken_output,
            "-s1", self.R1,
            "-s2", self.R2,
            "-t", self.taxid,
            "-o", self.r1_out,
            "-o2", self.r2_out,
            "--include-children",
            "--fastq-output",
            "--report", self.kraken_report
        ]
        
        if self.debug:
            print("Debug: Running command:")
            print(" ".join(cmd))
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Count reads and store metrics
            self.r1_count = count_fastq_reads(self.r1_out)
            self.r2_count = count_fastq_reads(self.r2_out)
            self.total_input_reads = count_fastq_reads(self.R1)
            
            # Parse relevant metrics from report
            self.parse_report_metrics()
            
            # Gzip output files
            gzip_file(self.r1_out)
            gzip_file(self.r2_out)
            
            if result.stderr:
                self.stderr = result.stderr
                
        except subprocess.CalledProcessError as e:
            print("Error running KrakenTools:", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
            sys.exit(1)
    
    def parse_report_metrics(self):
        """Parse key metrics from the Kraken report"""
        self.metrics = {
            'target_taxon': self.taxon,
            'taxid': self.taxid,
            'total_input_reads': self.total_input_reads,
            'extracted_reads': self.r1_count,
            'extraction_rate': round(self.r1_count / self.total_input_reads * 100, 2)
        }
        
        # Parse classification details from report lines
        classifications = []
        for line in self.relevant_report_lines:
            parts = re.split(r'\s+', line.strip())
            if len(parts) >= 6:
                percent = float(parts[0])
                reads = int(parts[1])
                level = parts[3]
                name = ' '.join(parts[5:])
                classifications.append({
                    'name': name,
                    'reads': reads,
                    'percent': percent,
                    'level': level
                })
        self.classifications = classifications

    def latex(self, tex):
        """Generate LaTeX report"""
        # First table banner and table
        kraken_banner = Banner("Kraken Classification")
        
        # First table - use table environment to keep banner and table together
        print(r'\begin{table}[H]', file=tex)
        print(r'\centering', file=tex)
        print(f'\\includegraphics[width=\\textwidth]{{{kraken_banner.banner}}}', file=tex)
        
        # First table with matching width
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{p{4cm}|p{8cm}}', file=tex)
        
        # Table header
        # print(r'\hline', file=tex)
        print(r'Metric & Value \\', file=tex)
        print(r'\hline', file=tex)
        
        # Table content
        print(f'Target Taxon & {self.metrics["target_taxon"]}\\\\', file=tex)
        print(f'Taxon ID & {self.metrics["taxid"]}\\\\', file=tex)
        print(f'Total Input Reads & {self.metrics["total_input_reads"]:,}\\\\', file=tex)
        print(f'Extracted Reads & {self.metrics["extracted_reads"]:,}\\\\', file=tex)
        print(f'Extraction Rate & {self.metrics["extraction_rate"]}\%\\\\', file=tex)
        print(r'\hline', file=tex)
        
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)

        # Second table banner and table
        detail_banner = Banner("Kraken Detailed Classification")
        continuation_banner = Banner("Kraken Detailed Classification (Continued from previous page)")

        # Get all classifications to determine table length
        all_classifications = self.classifications
        max_rows_per_page = 35  # Based on article class with 1-inch margins
        total_rows = len(all_classifications)
        pages_needed = (total_rows + max_rows_per_page - 1) // max_rows_per_page  # Ceiling division

        for page in range(pages_needed):
            start_idx = page * max_rows_per_page
            end_idx = min((page + 1) * max_rows_per_page, total_rows)
            page_classifications = all_classifications[start_idx:end_idx]
            
            # Table for this page
            print(r'\begin{table}[H]', file=tex)
            
            # Choose appropriate banner based on whether this is a continuation
            current_banner = continuation_banner.banner if page > 0 else detail_banner.banner
            
            # Use a zero-spaced vbox to eliminate spacing
            print(r'\noindent\makebox[\textwidth]{%', file=tex)
            print(f'\\includegraphics[width=\\textwidth]{{{current_banner}}}%', file=tex)
            print(r'}', file=tex)
            
            # No vertical space between banner and table
            print(r'\vspace{-0.5pt}%', file=tex)  # Tiny negative space to ensure flush alignment
            
            # Ensure table has exactly the same width as the banner (textwidth)
            print(r'\begin{adjustbox}{width=\textwidth}', file=tex)
            print(r'\begin{tabular}{p{8cm}|c|r|r}', file=tex)
            
            # Table header on each page
            print(r'Taxonomic Name & Level & Reads & Percent \\', file=tex)
            print(r'\hline', file=tex)
            
            # Table content for this page
            for classif in page_classifications:
                name = classif['name'].replace('_', '\\_')
                print(f'{name} & {classif["level"]} & {classif["reads"]:,} & {classif["percent"]:.2f}\\\\', file=tex)
            print(r'\hline', file=tex)
            
            # If this is not the last page, add continuation notice
            if page < pages_needed - 1:
                print(r'\multicolumn{4}{r}{\textit{Continued on next page...}} \\', file=tex)
            
            print(r'\end{tabular}', file=tex)
            print(r'\end{adjustbox}', file=tex)
            print(r'\end{table}', file=tex)

    def excel(self, excel_dict):
        """Add metrics to provided excel dictionary"""
        # Basic metrics
        excel_dict['Target Taxon'] = self.metrics['target_taxon']
        excel_dict['Taxon ID'] = self.metrics['taxid']
        excel_dict['Total Input Reads'] = self.metrics['total_input_reads']
        excel_dict['Extracted Reads'] = self.metrics['extracted_reads']
        excel_dict['Extraction Rate (%)'] = self.metrics['extraction_rate']
        
        # Create a condensed classification summary
        classification_summary = []
        for c in self.classifications:
            if c['percent'] > 1.0:  # Only include classifications > 1%
                classification_summary.append(
                    f"{c['name']} ({c['level']}): {c['percent']:.1f}%"
                )
        
        excel_dict['Classification Summary'] = '; '.join(classification_summary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='PROG',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        ---------------------------------------------------------
        Separate FASTQ reads based on Kraken taxonomic classification using KrakenTools.
        Processes paired-end reads and generates filtered output files based on 
        taxonomic classification.
        '''),
        epilog='---------------------------------------------------------'
    )
    
    parser.add_argument('-r1', '--read1', action='store', dest='R1', required=True,
                       help='Input R1 FASTQ file (can be gzipped)')
    parser.add_argument('-r2', '--read2', action='store', dest='R2', required=True,
                       help='Input R2 FASTQ file (can be gzipped)')
    parser.add_argument('-k', '--kraken-output', action='store', dest='kraken_output',
                       required=True, help='Kraken output file')
    parser.add_argument('-r', '--kraken-report', action='store', dest='kraken_report',
                       required=True, help='Kraken report file')
    parser.add_argument('-t', '--taxon', action='store', dest='taxon', required=True,
                       help='Taxonomic name to extract (e.g., "Orbivirus")')
    parser.add_argument('-o', '--output-prefix', action='store', dest='output_prefix',
                       required=True, help='Prefix for output FASTQ files')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug',
                       default=False, help='Print debug information')
    parser.add_argument('-l', '--latex', action='store_true', dest='build_latex',
                       default=False, help='Generate LaTeX report')
    parser.add_argument('-e', '--excel', action='store_true', dest='build_excel',
                       default=False, help='Generate Excel report')
    parser.add_argument('-v', '--version', action='version',
                       version=f'{os.path.basename(__file__)}: version {__version__}')
    
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    # Initialize ParseReads class
    parser = ParseReads(
        R1=args.R1,
        R2=args.R2,
        kraken_output=args.kraken_output,
        kraken_report=args.kraken_report,
        taxon=args.taxon,
        output_prefix=args.output_prefix,
        debug=args.debug
    )
    
    # Run the main processing
    parser.run()
    
    # Generate reports if requested
    if args.build_latex:
        latex_report = Latex_Report(parser.sample_name)
        parser.latex(latex_report.tex)
        latex_report.latex_ending()
        
    if args.build_excel:
        excel_stats = Excel_Stats(parser.sample_name)
        parser.excel(excel_stats.excel_dict)
        excel_stats.post_excel()