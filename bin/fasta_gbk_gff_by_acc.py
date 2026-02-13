#!/usr/bin/env python

__version__ = "2.04"

import os
import argparse
import textwrap
import time
from Bio import SeqIO
from Bio import Entrez

class Downloader:

    def __init__(self, accession):
        self.entrezDbName = 'nucleotide'
        self.email = 'mickey_mouse@gmail.com'  # Consider using a real email
        self.accession = accession

    def gbk(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gbk")
        try:
            # Add small delay to be nice to NCBI servers
            time.sleep(0.3)
            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gb')
            writeFile = self.accession + ".gbk"
            
            with open(writeFile, "w") as local_file:
                local_file.write(entryData.read())
            entryData.close()
            
            print(f"Successfully downloaded {writeFile}")
            
        except Exception as e:
            print(f"Error downloading gbk for {self.accession}: {e}")
            raise

    def gff(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gff3")
        try:
            # Add small delay to be nice to NCBI servers
            time.sleep(0.3)
            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gff3')
            writeFile = self.accession + ".gff"
            
            with open(writeFile, "w") as local_file:
                local_file.write(entryData.read())
            entryData.close()
            
            print(f"Successfully downloaded {writeFile}")
            
        except Exception as e:
            print(f"Error downloading gff for {self.accession}: {e}")
            raise

    def fasta(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} FASTA")
        try:
            # Add small delay to be nice to NCBI servers
            time.sleep(0.3)
            entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='fasta')
            writeFile = self.accession + ".fasta"
            
            # Read the data first
            fasta_content = entryData.read()
            entryData.close()
            
            # Write to file
            with open(writeFile, "w") as local_file:
                local_file.write(fasta_content)
            
            print(f"Successfully downloaded {writeFile}")
            
            # Parse and display sequence information
            with open(writeFile, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    print(f"{record.description}")
                    print(f'Sequence length: {len(record.seq):,}\n')
                    return record.description
                    
        except Exception as e:
            print(f"Error downloading fasta for {self.accession}: {e}")
            raise

    def verify_accession(self):
        """Verify that the accession exists before attempting download"""
        Entrez.email = self.email
        try:
            print(f"Verifying accession {self.accession}...")
            search_handle = Entrez.esearch(db=self.entrezDbName, term=self.accession)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            count = int(search_results["Count"])
            if count == 0:
                print(f"ERROR: Accession {self.accession} not found in {self.entrezDbName} database")
                return False
            else:
                print(f"Accession {self.accession} found in {self.entrezDbName} database")
                return True
                
        except Exception as e:
            print(f"Error verifying accession {self.accession}: {e}")
            return False
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Usage:
        fasta_gbk_gff_by_acc.py -a NC_002945 -f
        fasta_gbk_gff_by_acc.py -a NC_006932 -fg
        fasta_gbk_gff_by_acc.py -a NC_006933 -fbg
        fasta_gbk_gff_by_acc.py -a CP023243 -fbg
        fasta_gbk_gff_by_acc.py -a NZ_CP023243 -fbg
        fasta_gbk_gff_by_acc.py -a JX121104 -fbg
        **bad request: fasta_gbk_gff_by_acc.py -a NZ_AQME0000000 -fbg # must be complete chromosome

    vSNP requires multi-chromosome genomes to be concatenated to single file

    Search genomes: https://www.ncbi.nlm.nih.gov/genome

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--accession', action='store', dest='accession', required=True, help='NCBI accession number')
    parser.add_argument('-f', '--fasta', action='store_true', dest='fasta', help='get FASTA file')
    parser.add_argument('-b', '--gbk', action='store_true', dest='gbk', help='get gbk file')
    parser.add_argument('-g', '--gff', action='store_true', dest='gff', help='get gff file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')
    parser.add_argument('--verify', action='store_true', dest='verify', help='verify accession exists before download')

    args = parser.parse_args()
    
    # Check that at least one format is requested
    if not (args.fasta or args.gbk or args.gff):
        print("ERROR: Must specify at least one format (-f, -b, or -g)")
        parser.print_help()
        exit(1)
    
    download = Downloader(args.accession)
    
    # Optional verification step
    if args.verify:
        if not download.verify_accession():
            print("Exiting due to verification failure")
            exit(1)
    
    # Download requested formats
    try:
        if args.fasta:
            download.fasta()
        if args.gbk:
            download.gbk()
        if args.gff:
            download.gff()
            
        print("All requested downloads completed successfully!")
        
    except Exception as e:
        print(f"Download failed: {e}")
        exit(1)