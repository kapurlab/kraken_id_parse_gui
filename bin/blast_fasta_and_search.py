#!/usr/bin/env python

import os
import re
import shutil
import glob
import operator
from collections import defaultdict
from collections import Counter
import argparse
import textwrap
from datetime import datetime

from Bio import SeqIO

from file_setup import Setup, bcolors, Banner, Latex_Report, Excel_Stats, safe_move


class Blast_Fasta(Setup, bcolors):
    def __init__(self, FASTA=None, search=None, blast_out=None, format="6 qseqid sacc bitscore pident stitle", num_alignment=3, blast_db="nt"):
        Setup.__init__(self, FASTA=FASTA)
        sample_name = self.sample_name
        self.blast_db = blast_db
        blastout_file = f'{sample_name}_blast_out.txt'
        self.split_file_list = []  # Initialize split_file_list here

        # Detect HPC system
        self.hpc_system = self.detect_hpc_system()
        print(f'{bcolors.BLUE}Detected HPC system: {self.hpc_system}{bcolors.ENDC}')

        if blast_out:
            blastout_file = blast_out
        else:
            fasta_size = sum([len(seq_record.seq) for seq_record in SeqIO.parse(FASTA, "fasta")])
            
            # Check if we're in a SLURM environment
            is_slurm = shutil.which('sbatch') is not None
            
            if fasta_size < 50000:
                print(f'\n{bcolors.YELLOW}Running BLAST...{bcolors.ENDC}\n')
                if is_slurm and (os.path.isdir("/project") or os.path.isdir("/software/public/databases")):
                    # Use SLURM for execution
                    with open('batch.sh', 'w') as rsh:
                        rsh.write(
                            f'#!/bin/bash\n\n'
                            f'#SBATCH --ntasks=48\n'
                            f'#SBATCH --job-name="blast_db"\n'
                            f'#SBATCH --export=NONE\n\n'
                        )
                        # Add system-specific module loads
                        if self.hpc_system == "ames":
                            rsh.write(
                                # f'module load vsnp-2.0.3-gcc-9.2.0-ypwbztb\n'
                                f'module load blast-plus-2.12.0-gcc-9.2.0-pj4bk\n\n'
                            )
                        elif self.hpc_system == "scomp":
                            rsh.write(
                                f'module load ncbi-blast/2.17.0\n\n'
                            )
                        
                        rsh.write(
                            f'blastn -query {FASTA} -db {blast_db} -word_size 11 -out {self.sample_name}_blast_out.txt '
                            f'-outfmt "{format}" -num_alignments {num_alignment} -num_threads {self.cpus}'
                        )
                    os.system(f'sbatch -W ./batch.sh')
                else:
                    # Direct execution without SLURM
                    os.system(f'blastn -query {FASTA} -db {blast_db} -word_size 11 -out {self.sample_name}_blast_out.txt '
                            f'-outfmt "{format}" -num_alignments {num_alignment} -num_threads {self.cpus}')
            else:
                # Handle large files
                file_size = os.path.getsize(FASTA)
                if 10000 <= file_size <= 100000:
                    max_file_size = int(file_size/10)
                    print(f'Splitting to: {int(file_size/10):,} size files')
                elif file_size > 100000:
                    max_file_size = int(file_size/20)
                    print(f'Splitting to: {int(file_size/20):,} size files')
                else:
                    max_file_size = file_size
                
                # Split files
                records = SeqIO.parse(FASTA, "fasta")
                total_size = 0
                subset_list = []
                grp_list_of_list = []
                for record in records:
                    if max_file_size > total_size:
                        total_size = len(record) + total_size
                        subset_list.append(record)
                    else:
                        subset_list.append(record)
                        grp_list_of_list.append(subset_list)
                        total_size = 0
                        subset_list = []
                grp_list_of_list.append(subset_list)
                print(f'Split to {len(grp_list_of_list)} files')
                
                # Write split files
                count = 0
                self.split_file_list = []  # Reset the list before populating
                for each_list in grp_list_of_list:
                    count += 1
                    outfile = f'group_for_blast_{count}.fasta'
                    self.split_file_list.append(outfile)
                    SeqIO.write(each_list, outfile, "fasta")
                
                if is_slurm and (os.path.isdir("/project") or os.path.isdir("/software/public/databases")):
                    # SLURM batch processing
                    with open('batch.sh', 'w') as rsh:
                        rsh.write(
                            f'#!/bin/bash\n\n'
                            f'#SBATCH --ntasks=48\n'
                            f'#SBATCH --job-name="blst splt"\n'
                            f'#SBATCH --output=blastout-%A_%a.out\n'
                            f'#SBATCH --array=0-20\n'
                            f'#SBATCH --export=NONE\n'
                            f'#SBATCH --wait\n'
                            f'#SBATCH --nice\n\n'
                        )
                        # Add system-specific module loads
                        if self.hpc_system == "ames":
                            rsh.write(f'module load blast-plus-2.12.0-gcc-9.2.0-pj4bk\n\n')
                        elif self.hpc_system == "scomp":
                            rsh.write(f'module load ncbi-blast/2.17.0\n\n')
                        
                        rsh.write(
                            f'ls group_for_blast_*.fasta |cut -d. -f1 > jobs\n'
                            f'names=($(cat jobs))\n'
                            f'echo ${{names[${{SLURM_ARRAY_TASK_ID}}]}}\n'
                            f'blastn -query ${{names[${{SLURM_ARRAY_TASK_ID}}]}}.fasta -db {blast_db} '
                            f'-out ${{names[${{SLURM_ARRAY_TASK_ID}}]}}_blastout.txt -outfmt "{format}" '
                            f'-num_alignments {num_alignment} -num_threads {self.cpus}\n'
                            f'rm jobs\n'
                        )
                    os.system('sbatch -W ./batch.sh')
                else:
                    # Sequential processing without SLURM
                    total_files = len(self.split_file_list)
                    for i, split_file in enumerate(self.split_file_list, 1):
                        print(f'  BLAST chunk {i}/{total_files}: {split_file}', flush=True)
                        os.system(f'blastn -query {split_file} -db {blast_db} '
                                f'-out {split_file}_blastout.txt -outfmt "{format}" '
                                f'-num_alignments {num_alignment} -num_threads {self.cpus}')
                        print(f'  BLAST chunk {i}/{total_files} complete', flush=True)

                # Concatenate results
                concatenation = f'{self.sample_name}_blast_out.txt'
                with open(concatenation, 'wb') as outfile:
                    for filename in glob.glob('*_blastout.txt'):
                        if filename == concatenation:
                            continue
                        with open(filename, 'rb') as readfile:
                            shutil.copyfileobj(readfile, outfile)
                        os.remove(filename)
                
                # Cleanup
                for each in glob.glob('blastout-*.out'):
                    os.remove(each)
                for each in glob.glob('group_for_blast_*.fasta'):
                    os.remove(each)
        
        self.blastout_file = blastout_file
        
        # Process BLAST results
        try:
            blast_dict = defaultdict(list)
            with open(blastout_file, 'r') as blast_file:
                for line in blast_file:
                    line = line.rstrip()
                    line = line.split('\t')
                    blast_dict.setdefault(line[0], []).append(line[1:])
            
            top_hit_acc = []
            descriptions = {}
            with open(f'{sample_name}_blast_all.txt', 'w') as all_blast:
                for item, value in blast_dict.items():
                    print(f'{item}', file=all_blast)
                    for val in value:
                        print(f'\t{val[0]} {val[1]} {val[2]} {val[3]}', file=all_blast)
                    print(f'', file=all_blast)

            if search:
                node_list = []
                for item, value in blast_dict.items():
                    for val in value:
                        for v in val:
                            if search.lower() in v.lower():
                                node_list.append(item.lower())
                node_list = set(node_list)
                found_record = []
                for seq_record in SeqIO.parse(FASTA, "fasta"):
                    if seq_record.description.lower() in node_list:
                        found_record.append(seq_record)
                term = search.lower()
                term = re.sub('[\/.!@#$%^&*()+,"= ]', '_', term)
                SeqIO.write(found_record, f'{sample_name}_search_{term}.fasta', 'fasta')

            acc_frequency = []
            top_hit_acc_norm = []
            norm_dict = {}
            descriptions = {}
            for header, description in blast_dict.items():
                # Get accession frequencies
                acc_frequency.append(description[0][0])  # top hit, 1st item in hit
                acc_count = Counter()
                for acc in acc_frequency:
                    acc_count[acc] += 1
                # Get FASTA sizes per accessions
                fasta_dict = SeqIO.to_dict(SeqIO.parse(FASTA, "fasta"))
                seq_length = len(fasta_dict[header])
                top_hit_acc_norm.append((description[0][0], seq_length))
                acc_size_collection = defaultdict(list)
                for acc, size in top_hit_acc_norm:  # collect sizes by accession
                    acc_size_collection[acc].append(size)
                for acc, sizes in acc_size_collection.items():  # add sizes by accession
                    norm_dict[acc] = sum(sizes)
                # Get accession descriptions
                descriptions[description[0][0]] = description[0][3]
            sorted_norm_dict = {k: v for k, v in sorted(norm_dict.items(), key=lambda item: item[1])}

            for value in blast_dict.values():
                top_hit_acc.append(value[0][0])  # top hit, 1st item in hit
                descriptions[value[0][0]] = value[0][3]
            cnt = Counter()
            for acc in top_hit_acc:
                cnt[acc] += 1
            self.top_hit_acc = top_hit_acc

            summary_dict = {}
            summary_list = []  # list of tuples (nt rep, contigs, description list)
            blast_summary_file = f'{sample_name}_blast_summary.txt'
            self.blast_summary_file = blast_summary_file
            with open(blast_summary_file, 'w') as summary_blast:
                for acc, count in sorted_norm_dict.items():
                    summary_dict[f'{acc} {descriptions[acc]}'] = f'{count}'
                    summary_list.append((f'{count:,}', f'{acc_count[acc]:,}', descriptions[acc]))
                    print(f'{count:,}\t{acc_count[acc]:,}\t{acc} {descriptions[acc]}', file=summary_blast)
                    print(f'{bcolors.YELLOW}{count:,}{bcolors.ENDC} nt\t{bcolors.RED}{acc_count[acc]:,}{bcolors.ENDC} contigs\t'
                        f'{bcolors.BLUE}{int(round(count/acc_count[acc])):,}{bcolors.ENDC} nt mean length\t'
                        f'{bcolors.WHITE}{acc} {descriptions[acc]}{bcolors.ENDC}')
            
            # Get highest hit on most nucleotide identifying as a single accession
            try:
                highest_hit_accession = max(sorted_norm_dict.items(), key=operator.itemgetter(1))[0]
                self.highest_hit_description_list = descriptions[highest_hit_accession]
            except ValueError:
                print(f"{bcolors.RED}Warning: No BLAST hits found{bcolors.ENDC}")
                self.highest_hit_description_list = None
            
            self.summary_dict = summary_dict
            self.summary_list = summary_list

        except FileNotFoundError:
            print(f"{bcolors.RED}Error: BLAST output file not found. BLAST may have failed.{bcolors.ENDC}")
            self.highest_hit_description_list = None
            self.summary_dict = {}
            self.summary_list = []

    def detect_hpc_system(self):
        """
        Detect which HPC system we're running on based on available directories
        Returns: 'ames', 'scomp', or 'unknown'
        """
        if os.path.isdir("/project/bioinformatic_databases"):
            return "ames"
        elif os.path.isdir("/software/public/databases"):
            return "scomp"
        else:
            print(f"{bcolors.YELLOW}Warning: Could not detect HPC system. Using default settings.{bcolors.ENDC}")
            return "unknown"

    def latex(self, tex):
        blast_string_int = {k: int(v) for k, v in self.summary_dict.items()}
        blast_sorted = sorted(blast_string_int.items(), key=operator.itemgetter(1), reverse=True)
        self.blast_sorted = blast_sorted
        table_breaks = [0, 55, 110, 165, 220, 275, 330, 385, 440, 495, 550, 605, 660, 715, 770, 825, 880, 935, 990, 1045, 1100, 1155, 1210]
        count = 0
        basename = os.path.basename(self.blast_db)
        # basename = basename.replace("_", "\_")
        blast_banner = Banner(f'BLAST {basename} - Assembly Identification')
        for break_start, break_end in zip(table_breaks, table_breaks[1:]):
            #split table if needed:
            if len(blast_sorted) > break_start:
                if count > 0:
                    blast_banner = Banner(f'BLAST {basename} - Assembly Identification - continued')
                print(r'\begin{table}[H]', file=tex)
                print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
                print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
                print(r'\end{adjustbox}', file=tex)
                print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
                print(r'\begin{tabular}{ l | p{1.3cm} | l }', file=tex)
                print(f'nt base count & contigs & Description \\\\', file=tex)
                print(r'\hline', file=tex)
                # for each_row in self.summary_list[-1]:
                #     flip_list.append(each_row)
                for each_row in self.summary_list[::-1][break_start:(break_end - 1)]:
                    description = f'{each_row[2]}'
                    description = description.replace("_", "\_")[:108]
                    print(f'{each_row[0]} & {each_row[1]} & {description} \\\\', file=tex)
                print(r'\hline', file=tex)
                # Close tabular before the adjustbox that wraps it (LIFO); drop stray \\.
                print(r'\end{tabular}', file=tex)
                print(r'\end{adjustbox}', file=tex)
                print(r'\vspace{0.1 mm}', file=tex)
                basename_slashed = basename.replace("_", "\_")
                print(r'\begin{flushleft}Results provided by: \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi}{BLAST ' + f'{basename_slashed}' + r' database}\end{flushleft}', file=tex)
                print(r'\end{table}', file=tex)
                count += 1

    def excel(self, excel_dict):
        basename = os.path.basename(self.blast_db)
        try:
            excel_dict[f'Top BLAST {basename} Hit - 1'] = f'{self.summary_list[-1][0]} {basename}, {self.summary_list[-1][1]} contigs of {self.summary_list[-1][2][0]} {self.summary_list[-1][2][3]}'
        except IndexError:
            excel_dict[f'Top BLAST {basename} Hit - 1'] = 'BLAST failed - No Results Output From BLAST search'
        try:
            excel_dict[f'Top BLAST {basename} Hit - 2'] = f'{self.summary_list[-2][0]} {basename} {self.summary_list[-2][1]} contigs of {self.summary_list[-2][2][0]} {self.summary_list[-2][2][3]}'
        except IndexError:
            pass
        try:
            excel_dict[f'Top BLAST {basename} Hit - 3'] = f'{self.summary_list[-3][0]} {basename} {self.summary_list[-3][1]} contigs of {self.summary_list[-3][2][0]} {self.summary_list[-3][2][3]}'
        except IndexError:
            pass


class Spades_Stats:
    '''

    '''
    def __init__(self, fasta_in,):
        base_name = os.path.basename(fasta_in)
        self.sample_name = re.sub('[._].*', '', base_name)
        records = list(SeqIO.parse(fasta_in, "fasta"))
        cov_length_list=[]
        contig_count = 0
        coverage_list=[]
        length_list=[]
        small_contigs=[]
        greater_one_kb=[]
        mid_size = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
            except IndexError:
                coverage_value = 1
            try:
                coverage_value = int(float(coverage_value))
            except ValueError:
                coverage_value = 1
            cov_length_list.append({'name': rec.description, 'cov': coverage_value, 'length': len(rec)})
            coverage_list.append(coverage_value)
            length_list.append(len(rec))
            contig_count += 1
            if len(rec) <= 300:
                small_contigs.append(len(rec))
            elif len(rec) >= 1000:
                greater_one_kb.append(len(rec))
            else:
                mid_size.append(len(rec))
        total_contig_lengths = int(sum(length_list))

        #Calculate mean coverage
        normalized_list = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
            except IndexError:
                coverage_value = 1
            try:
                coverage_value = int(float(coverage_value))
            except ValueError:
                coverage_value = 1
            normalized_list.append((len(rec) / total_contig_lengths) * coverage_value)
        spades_mean_coverage = sum(normalized_list)

        #N50 calculation
        all_len = sorted(length_list, reverse=True)
        csum = np.cumsum(all_len)
        n2 = int(sum(length_list)/2)
        csumn2 = min(csum[csum >= n2])
        ind = np.where(csum == csumn2)
        self.n50 = all_len[int(ind[0])] # n50 smallest size contig which, along with the larger contigs, contain half of sequence of a particular genome
        self.l50 = int(ind[0][0]) + 1 # l50 smallest number of contigs whose length sum makes up half of genome
        self.cov_length_list = cov_length_list
        self.contig_count = int(contig_count)
        self.longest_contig = int(max(length_list))
        self.total_contig_lengths = total_contig_lengths
        self.spades_mean_coverage = spades_mean_coverage
        self.small_contigs_count = len(small_contigs)
        self.greater_one_kb_count = len(greater_one_kb)
        self.mid_size = len(mid_size)
        self.spades_version = os.popen("spades.py --version").readlines()[0]
        
    def print_by_coverage(self,):
        for each_dict in sorted(self.cov_length_list, key=itemgetter('cov')):
            print(f'{each_dict["cov"]:,}X  {each_dict["length"]:,}  {each_dict["name"]}')

    def print_by_length(self,):
        for each_dict in sorted(self.cov_length_list, key=itemgetter('length')):
            print(f'{each_dict["length"]:,}  {each_dict["cov"]:,}X  {each_dict["name"]}')

    def write_stats(self, fq=None, build_latex=None, build_excel=None, message=None):

        fastq_coverage_title = 'Mean Read Depth: read count * read size / total assembly length'
        spades_coverage_title = 'Largest k Value Mean Depth: SPAdes reporting'

        print(
            f'\n'
            f'Contig count: {bcolors.YELLOW}{self.contig_count:,}{bcolors.ENDC} \n'
            f'Contig length counts <|301-999bp|>: {bcolors.RED}{self.small_contigs_count:,}{bcolors.ENDC}|{bcolors.BLUE}{self.mid_size:,}{bcolors.ENDC}|{bcolors.GREEN}{self.greater_one_kb_count:,}{bcolors.ENDC} \n'
            f'Longest contig: {bcolors.GREEN}{self.longest_contig:,}{bcolors.ENDC} \n'
            f'Total length: {bcolors.BLUE}{self.total_contig_lengths:,}{bcolors.ENDC} \n'
            f'N50: {bcolors.UNDERLINE}{self.n50:,}{bcolors.ENDC} \n'
            f'{spades_coverage_title}: {bcolors.YELLOW}{self.spades_mean_coverage:,.1f}X{bcolors.ENDC}\n'
            )

        if fq: #calculating from FASTQ reads more accurate than spades reporting
            fastq_mean_coverage = ((fq.read1.total_read_count * fq.read1.length_mean) * 2) / self.total_contig_lengths
            mean_coverage_latex_title = 'Mean Read Depth'
            self.mean_coverage = fastq_mean_coverage  #default to read coverage on latex document
            f'{fastq_coverage_title}: {bcolors.YELLOW}{fastq_mean_coverage:,.1f}X{bcolors.ENDC} \n'
        else:
            mean_coverage_latex_title = 'Largest k Value Mean Depth'
            self.mean_coverage = self.spades_mean_coverage

        if build_excel is None: #just write out a default excel file.
            #stats to excel
            ts = time.time()
            st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
            sample_name = self.sample_name
            df = pd.DataFrame(index=[sample_name], \
                columns=[ \
                    'Read 1', 'R1 File Size', 'R1 Total Reads', 'R1 Mean Length', 'R1 Mean Quality', 'R1 Passing Q30', \
                    'Read 2', 'R2 File Size', 'R2 Total Reads', 'R2 Mean Length', 'R2 Mean Quality', 'R2 Passing Q30', \
                    'Assembly Contig Count', '<300bp Count', '301-999bp Count', '>1kb Count', 'Total Length', 'Longest Contig', 'N50', fastq_coverage_title, spades_coverage_title,])
            try:
                df.at[sample_name, 'Read 1'] = f'{fq.read1.fastq}'
                df.at[sample_name, 'R1 File Size'] = f'{fq.read1.file_size}'
                df.at[sample_name, 'R1 Total Reads'] = f'{fq.read1.total_read_count:,}'
                df.at[sample_name, 'R1 Mean Length'] = f'{fq.read1.length_mean:.1f}'
                df.at[sample_name, 'R1 Mean Quality'] = f'{fq.read1.read_average:.1f}'
                df.at[sample_name, 'R1 Passing Q30'] = f'{fq.read1.reads_gt_q30/fq.read1.sampling_size:0.1%}'
                all_reads = fq.read1.total_read_count
            except AttributeError:
                df.at[sample_name, 'Read 1'] = 'NA'
                df.at[sample_name, 'R1 File Size'] = 'NA'
                df.at[sample_name, 'R1 Total Reads'] = 'NA'
                df.at[sample_name, 'R1 Mean Length'] = 'NA'
                df.at[sample_name, 'R1 Mean Quality'] = 'NA'
                df.at[sample_name, 'R1 Passing Q30'] = 'NA'
            try: 
                df.at[sample_name, 'Read 2'] = f'{fq.read2.fastq}'
                df.at[sample_name, 'R2 File Size'] = f'{fq.read2.file_size}'
                df.at[sample_name, 'R2 Total Reads'] = f'{fq.read2.total_read_count:,}'
                df.at[sample_name, 'R2 Mean Length'] = f'{fq.read2.length_mean:.1f}'
                df.at[sample_name, 'R2 Mean Quality'] = f'{fq.read2.read_average:.1f}'
                df.at[sample_name, 'R2 Passing Q30'] = f'{fq.read2.reads_gt_q30/fq.read2.sampling_size:0.1%}'
                all_reads = all_reads + fq.read2.total_read_count
            except AttributeError:
                df.at[sample_name, 'Read 2'] = 'NA'
                df.at[sample_name, 'R2 File Size'] = 'NA'
                df.at[sample_name, 'R2 Total Reads'] = 'NA'
                df.at[sample_name, 'R2 Mean Length'] = 'NA'
                df.at[sample_name, 'R2 Mean Quality'] = 'NA'
                df.at[sample_name, 'R2 Passing Q30'] = 'NA'
            df.at[sample_name, 'Assembly Contig Count'] = f'{self.contig_count:,}'
            df.at[sample_name, '<300bp Count'] = f'{self.small_contigs_count:,}'
            df.at[sample_name, '301-999bp Count'] = f'{self.mid_size:,}'
            df.at[sample_name, '>1kb Count'] = f'{self.greater_one_kb_count:,}'
            df.at[sample_name, 'Total Length'] = f'{self.total_contig_lengths:,}'
            df.at[sample_name, 'Longest Contig'] = f'{self.longest_contig:,}'
            df.at[sample_name, 'N50'] = f'{self.n50:,}'
            if fq: #calculating from FASTQ reads more accurate than spades reporting
                df.at[sample_name, fastq_coverage_title] = f'{fastq_mean_coverage:,.1f}X'
                df.at[sample_name, spades_coverage_title] = f'{self.spades_mean_coverage:,.1f}X'
            else:  #when no FASTQ info just use spades reporting
                df.at[sample_name, fastq_coverage_title] = 'NA'
                df.at[sample_name, spades_coverage_title] = f'{self.spades_mean_coverage:,.1f}X'
            df.index.name = 'sample'
            df.to_excel(f'{sample_name}_{st}_stats.xlsx')
            self.self_excel = f'{os.getcwd()}/{sample_name}_{st}_stats.xlsx'

        if build_latex:
            tex = build_latex
            print(r'\begin{table}[H]', file=tex)
            print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
            assembly_banner = Banner("Assembly Metrics")
            print('\includegraphics[scale=1]{' + assembly_banner.banner + '}', file=tex)
            print(r'\end{adjustbox}', file=tex)
            print(r'\centering', file=tex)
            print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
            print(r'\begin{tabular}{l|l|l|l|l|l}', file=tex)
            print(r'Scaffolds & Total Length & Longest Scaffold & Scaffolds \textgreater 1K nt & N50 & ' + mean_coverage_latex_title +  r' \\', file=tex)
            print(r'\hline', file=tex)
            print(f'{self.contig_count:,} & {self.total_contig_lengths:,} & {int(self.longest_contig):,} & {self.greater_one_kb_count:,} & {self.n50:,} & {self.mean_coverage:,.1f}X \\\\', file=tex)
            print(r'\hline', file=tex)
            print(r'\end{tabular}', file=tex)
            print(r'\end{adjustbox}', file=tex)
            print(r'\begin{flushleft}Results provided by: \href{http://cab.spbu.ru/software/spades/}{SPAdes}\end{flushleft}', file=tex)
            print(r'\end{table}', file=tex)

        if isinstance(build_excel, pd.DataFrame): #cannot use df as true value, must test with isinstance
            df = build_excel
            df.at[df.index[0], 'Assembly Contig Count'] = f'{self.contig_count:,}'
            df.at[df.index[0], '<300bp Count'] = f'{self.small_contigs_count:,}'
            df.at[df.index[0], '301-999bp Count'] = f'{self.mid_size:,}'
            df.at[df.index[0], '>1kb Count'] = f'{self.greater_one_kb_count:,}'
            df.at[df.index[0], 'Total Length'] = f'{self.total_contig_lengths:,}'
            df.at[df.index[0], 'Longest Contig'] = f'{self.longest_contig:,}'
            df.at[df.index[0], 'N50'] = f'{self.n50:,}'
            df.at[df.index[0], spades_coverage_title] = f'{self.spades_mean_coverage:,.1f}X'
            if fq:
                fastq_mean_coverage = ((fq.read1.total_read_count * fq.read1.length_mean) * 2) / self.total_contig_lengths
                df.at[df.index[0], fastq_coverage_title] = f'{fastq_mean_coverage:,.1f}X'


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    To search after BLAST has been completed use -b option along with -s options.
    If a fasta file is not provide the _blast_out.txt and search must be provide (-b -s).

    blast_fasta_and_search.py --> $ blast_fasta_and_search.py -f <FASTA in file> -s "search term"
                                    blast_fasta_and_search.py -f *fasta -x ref_viruses_rep_genomes

    blast_report_split_sample.py can be used for large assembly file.
    blast_report_split_sample.py -f *.fasta
    blast_fasta_and_search.py -f *.fasta -b *_blast_out.txt -s "search term"
    Search term is not case sensitive

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=False, help='REQUIRED: In file to be processed')
    parser.add_argument('-s', '--search', action='store', dest='search', required=False, help='Search Term: provide a term to select on')
    parser.add_argument('-b', '--blast_out', action='store', dest='blast_out', required=False, help='Provide _blast_out.txt file to skip the BLAST command portion if it has alread been done')
    parser.add_argument('-n', '--num_alignment', action='store', dest='num_alignment', default=3, required=False, help='Number of alignments to return')
    parser.add_argument('-t', '--format', action='store', dest='format', default="6 qseqid sacc bitscore pident stitle", required=False, help='BLAST output format')
    parser.add_argument('-x', '--blast_db', action='store', dest='blast_db', default="nt", required=False, help='BLAST output format')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')

    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)

    #Main script
    blast = Blast_Fasta(FASTA=args.fasta, search=args.search, blast_out=args.blast_out, format=args.format, num_alignment=args.num_alignment, blast_db=args.blast_db)

    #Latex report
    latex_report = Latex_Report(blast.sample_name)
    blast.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(blast.sample_name)
    blast.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        safe_move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# blast_fasta_and_search.py - Created January 2021 by Tod Stuber