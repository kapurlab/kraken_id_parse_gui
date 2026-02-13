#!/usr/bin/env python

__version__ = "3.27"

import os
import shutil
import argparse
import textwrap
import time
from Bio import SeqIO
from Bio import Entrez


class Downloader:

    def __init__(self, accession, cleanup=False, cache_dir=None):
        self.entrezDbName = 'nucleotide'
        self.email = 'mickey_mouse@gmail.com'  # Consider using a real email
        self.accession = accession
        self.cleanup = cleanup
        self.cache_dir = cache_dir
        self.downloaded_files = []

    def gbk(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gbk")
        try:
            # Add rate limiting to prevent HTTP 400 errors
            time.sleep(1.0)  # NCBI recommends 1 second between requests
            
            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gbwithparts')
            writeFile = self.accession + ".gbk"
            
            with open(writeFile, "w") as local_file:
                local_file.write(entryData.read())
            entryData.close()
            
            self.downloaded_files.append(writeFile)
            print(f"✓ Successfully downloaded {writeFile}")
            
        except Exception as e:
            print(f"✗ Error downloading gbk for {self.accession}: {e}")
            raise

    def gff(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gff3")
        try:
            # Add rate limiting to prevent HTTP 400 errors
            time.sleep(1.0)  # NCBI recommends 1 second between requests
            
            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gff3')
            writeFile = self.accession + ".gff"
            
            with open(writeFile, "w") as local_file:
                local_file.write(entryData.read())
            entryData.close()
            
            self.downloaded_files.append(writeFile)
            print(f"✓ Successfully downloaded {writeFile}")
            
        except Exception as e:
            print(f"✗ Error downloading gff for {self.accession}: {e}")
            raise

    # ---- reference cache helpers ----

    def _cache_path(self):
        """Return the cached FASTA path, or None if caching is disabled."""
        if not self.cache_dir:
            return None
        return os.path.join(self.cache_dir, f"{self.accession}.fasta")

    def _try_cache_hit(self, write_file):
        """Check local cache; copy to *write_file* on hit.  Returns description or None."""
        cached = self._cache_path()
        if not cached or not os.path.isfile(cached) or os.path.getsize(cached) == 0:
            return None  # cache miss
        shutil.copy2(cached, write_file)
        self.downloaded_files.append(write_file)
        print(f"Cache hit: {self.accession}")
        description = None
        with open(write_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                print(f"{record.description}")
                print(f'Sequence length: {len(record.seq):,}\n')
                description = record.description
        return description

    def _write_to_cache(self, write_file):
        """Copy a successfully-downloaded FASTA into the cache directory."""
        cached = self._cache_path()
        if not cached:
            return
        os.makedirs(self.cache_dir, exist_ok=True)
        shutil.copy2(write_file, cached)

    # ---- download methods ----

    def fasta(self):
        writeFile = self.accession + ".fasta"

        # --- try local cache first ---
        desc = self._try_cache_hit(writeFile)
        if desc is not None:
            return desc

        # --- download from NCBI ---
        Entrez.email = self.email
        print(f"Downloading {self.accession} FASTA")
        try:
            # Add rate limiting to prevent HTTP 400 errors
            time.sleep(1.0)  # NCBI recommends 1 second between requests

            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='fasta')

            # Read data first
            fasta_data = entryData.read()
            entryData.close()

            # Write to file
            with open(writeFile, "w") as local_file:
                local_file.write(fasta_data)

            self.downloaded_files.append(writeFile)

            # Populate cache for next time
            self._write_to_cache(writeFile)

            # Parse and display sequence information
            description = None
            with open(writeFile, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    print(f"{record.description}")
                    print(f'Sequence length: {len(record.seq):,}\n')
                    description = record.description

            print(f"✓ Successfully downloaded {writeFile}")
            return description

        except Exception as e:
            print(f"✗ Error downloading fasta for {self.accession}: {e}")
            raise

    def fasta_with_retry(self, max_retries=3):
        """
        Download FASTA with retry logic and exponential backoff.
        This is useful when called from scripts that download many sequences.
        """
        for attempt in range(max_retries):
            try:
                return self.fasta()
            except Exception as e:
                print(f"⚠ Attempt {attempt + 1}/{max_retries} failed for {self.accession}: {e}")
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 5  # 5, 10, 15 seconds exponential backoff
                    print(f"  Waiting {wait_time} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"✗ All {max_retries} attempts failed for {self.accession}")
                    raise

    def cleanup_files(self):
        """Remove downloaded files"""
        if self.cleanup:
            for file in self.downloaded_files:
                try:
                    if os.path.exists(file):
                        os.remove(file)
                        print(f"Removed: {file}")
                except OSError as e:
                    print(f"Error removing {file}: {e}")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Usage:
        download_fasta_by_acc.py -a NC_002945 -f
        download_fasta_by_acc.py -a NC_006932 -fg
        download_fasta_by_acc.py -a NC_006933 -fbg
        download_fasta_by_acc.py -a CP023243 -fbg
        download_fasta_by_acc.py -a NZ_CP023243 -fbg
        download_fasta_by_acc.py -a KM580418 -f  # segmented virus example
        **bad request: download_fasta_by_acc.py -a NZ_AQME0000000 -fbg # must be complete chromosome

    For segmented viruses, this script will download each segment individually.
    Use with caution regarding NCBI rate limits when downloading many segments.

    Search genomes: https://www.ncbi.nlm.nih.gov/genome

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--accession', action='store', dest='accession', required=True, help='NCBI accession number')
    parser.add_argument('-f', '--fasta', action='store_true', dest='fasta', help='get FASTA file')
    parser.add_argument('-b', '--gbk', action='store_true', dest='gbk', help='get gbk file')
    parser.add_argument('-g', '--gff', action='store_true', dest='gff', help='get gff file')
    parser.add_argument('-c', '--cleanup', action='store_true', dest='cleanup', help='remove downloaded files after completion')
    parser.add_argument('-r', '--retry', action='store_true', dest='use_retry', help='use retry logic with exponential backoff')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    
    # Check that at least one format is requested
    if not (args.fasta or args.gbk or args.gff):
        print("ERROR: Must specify at least one format (-f, -b, or -g)")
        parser.print_help()
        exit(1)
    
    download = Downloader(args.accession, cleanup=args.cleanup)
    
    try:
        if args.fasta:
            if args.use_retry:
                download.fasta_with_retry()
            else:
                download.fasta()
        if args.gbk:
            download.gbk()
        if args.gff:
            download.gff()
            
        print("All requested downloads completed successfully!")
        
    except Exception as e:
        print(f"Download failed: {e}")
        exit(1)
    finally:
        download.cleanup_files()