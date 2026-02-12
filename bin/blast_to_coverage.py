#!/usr/bin/env python3

__version__ = "0.0.1"

import os
import sys
import re
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Optional

from download_fasta_by_acc import Downloader
from file_setup import Setup, bcolors

class BlastCoverageBridge(Setup):
    """Bridge between BLAST results and coverage analysis"""
    
    def __init__(self, blast_summary: str, max_refs: int = 100, debug: bool = False):
        """
        Initialize bridge processing
        
        Args:
            blast_summary: Path to BLAST summary file
            max_refs: Maximum number of references to process
            debug: Enable debug output
        """
        super().__init__(debug=debug)
        self.blast_summary = blast_summary
        self.max_refs = max_refs
        self.debug = debug
        self.downloaded_fastas = []
        self.combined_fasta = None

    def _parse_blast_summary(self) -> List[str]:
        """Parse BLAST summary file to get accession numbers"""
        accessions = []
        seen = set()
        
        try:
            # Read the summary file and parse lines
            with open(self.blast_summary) as f:
                lines = f.readlines()
            
            # Convert to dataframe for easier sorting
            parsed_data = []
            for line in lines:
                count, _, description = line.strip().split('\t')
                count = int(count.replace(',', ''))  # Handle comma in numbers
                acc = description.split()[0]  # First word is accession
                parsed_data.append({'count': count, 'accession': acc, 'description': description})
            
            df = pd.DataFrame(parsed_data)
            # Sort by count in descending order
            df = df.sort_values('count', ascending=False)
            
            # Get unique accessions maintaining order
            for _, row in df.iterrows():
                acc = row['accession']
                if acc not in seen and not acc.startswith('Query_'):
                    seen.add(acc)
                    accessions.append(acc)
                    if len(accessions) >= self.max_refs:
                        break
                    
            if self.debug:
                print("\nDebug: Top accessions by count:")
                for acc in accessions:
                    print(f"  {acc}")
                    
        except Exception as e:
            print(f"Error parsing BLAST summary: {e}", file=sys.stderr)
            return []

        return accessions

    def _download_references(self, accessions: List[str]) -> List[Tuple[str, int]]:
        """Download reference sequences and return list of (file, size) tuples"""
        downloaded = []
        
        for acc in accessions:
            try:
                downloader = Downloader(acc,)
                description = downloader.fasta()  # Downloads file and returns description
                
                fasta_file = f"{acc}.fasta"
                if os.path.exists(fasta_file):
                    # Get sequence length
                    with open(fasta_file) as handle:
                        record = next(SeqIO.parse(handle, "fasta"))
                        seq_len = len(record.seq)
                        downloaded.append((fasta_file, seq_len))
                        self.downloaded_fastas.append(fasta_file)
                
                if self.debug:
                    print(f"Debug: Downloaded {acc}, length: {seq_len:,}")
                    
            except Exception as e:
                print(f"Warning: Failed to download {acc}: {e}", file=sys.stderr)
                continue
                
        return downloaded

    def _combine_fastas(self, fasta_files: List[str]) -> str:
        """
        Combine multiple FASTA files into one and cleanup individual files
        Returns path to combined file or None if failed
        """
        combined_records = []
        output_file = "combined_references.fasta"
        
        # First collect all records
        for fasta in fasta_files:
            try:
                records = list(SeqIO.parse(fasta, "fasta"))
                combined_records.extend(records)
                if self.debug:
                    print(f"Added {len(records)} records from {fasta}")
            except Exception as e:
                print(f"Warning: Error processing {fasta}: {e}", file=sys.stderr)
        
        # Write combined file if we have records
        if combined_records:
            SeqIO.write(combined_records, output_file, "fasta")
            
            # Clean up individual FASTA files
            for fasta in fasta_files:
                try:
                    os.remove(fasta)
                    if self.debug:
                        print(f"Removed individual FASTA: {fasta}")
                except OSError as e:
                    print(f"Warning: Could not remove {fasta}: {e}", file=sys.stderr)
            
            return output_file
        return None

    def process(self) -> Optional[str]:
        """
        Process BLAST results and prepare references for coverage analysis
        
        Returns:
            Path to combined FASTA file or None if processing failed
        """
        # Get accession numbers from BLAST summary
        accessions = self._parse_blast_summary()
        if not accessions:
            print("No valid accessions found in BLAST summary", file=sys.stderr)
            return None
            
        if self.debug:
            print(f"Debug: Processing accessions in order: {accessions}")
        
        # Download and get sizes
        downloaded = self._download_references(accessions)
        if not downloaded:
            print("No references could be downloaded", file=sys.stderr)
            return None
            
        # Sort by size (largest first)
        downloaded.sort(key=lambda x: x[1], reverse=True)
        
        if self.debug:
            print("\nDebug: References ordered by size:")
            for fasta, size in downloaded:
                print(f"  {fasta}: {size:,} bp")
        
        # Combine ordered FASTAs
        ordered_fastas = [f[0] for f in downloaded]
        combined_file = self._combine_fastas(ordered_fastas)

        if combined_file and os.path.exists(combined_file):
            self.combined_fasta = combined_file
            return combined_file

        
        return None

    def cleanup(self):
        """Remove individual FASTA files after combination"""
        if not self.debug:
            for fasta in self.downloaded_fastas:
                try:
                    os.remove(fasta)
                except OSError:
                    pass

def main():
    """Run bridge processing from command line"""
    import argparse
    parser = argparse.ArgumentParser(description="Process BLAST results for coverage analysis")
    parser.add_argument("blast_summary", help="Path to BLAST summary file")
    parser.add_argument("-n", "--max-refs", type=int, default=15, help="Maximum number of references to process")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()
    
    bridge = BlastCoverageBridge(args.blast_summary, args.max_refs, args.debug)
    combined_fasta = bridge.process()
    
    if combined_fasta:
        print(f"Combined reference file created: {combined_fasta}")
        if not args.debug:
            bridge.cleanup()
    else:
        print("Failed to create combined reference file", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
