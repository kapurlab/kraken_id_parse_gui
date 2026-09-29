#!/usr/bin/env python

import os
import re
import shutil
import glob
import subprocess
import operator
from collections import defaultdict
from collections import Counter
import argparse
import textwrap
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from file_setup import Setup, bcolors, Excel_Stats, safe_move


# DB prefixes already pulled into the page cache by this process, so the two or
# three Blast_Fasta rounds of one pipeline run pay the disk cost only once.
_PREWARMED_DBS = set()


def _prewarm_blast_db(blast_db):
    """Sequentially read a local BLAST DB's volumes into the OS page cache.

    blastn memory-maps its database and touches it in an essentially random
    order. Against a cold cache on rotational storage that is a 4 KB seek per
    page fault — a single query can sit in disk-wait for an hour while every
    thread idles, which reads as "BLAST stopped working". One sequential pass
    beforehand moves the same bytes at streaming speed, after which blastn is
    CPU-bound and its -num_threads actually help.

    Only does anything when the DB is local files and comfortably fits in the
    memory the kernel could actually give to cache (MemAvailable); otherwise it
    says why and lets blastn run as before. Linux-only by design — /proc/meminfo
    is the honest availability signal, and the macOS installs run small DBs.
    """
    if blast_db in _PREWARMED_DBS:
        return
    volume_re = re.compile(r'\.(\d+\.)?[np][a-z]{2}$')
    volumes = [f for f in glob.glob(glob.escape(blast_db) + '*')
               if volume_re.search(f) and os.path.isfile(f)]
    if not volumes:
        return  # not a local DB path ("nt" remote or a bad path) — blastn will say so
    total = sum(os.path.getsize(f) for f in volumes)
    try:
        with open('/proc/meminfo') as fh:
            mem_available = next(int(line.split()[1]) * 1024
                                 for line in fh if line.startswith('MemAvailable:'))
    except (OSError, StopIteration, ValueError, IndexError):
        return  # no /proc/meminfo (macOS) — keep prior behavior
    if total > 0.6 * mem_available:
        print(f'{bcolors.YELLOW}BLAST DB is {total / 1e9:.1f} GB but only '
              f'{mem_available / 1e9:.1f} GB RAM is available — skipping cache '
              f'prewarm; BLAST may be disk-bound.{bcolors.ENDC}', flush=True)
        return
    print(f'{bcolors.BLUE}Prewarming BLAST DB into RAM cache '
          f'({total / 1e9:.1f} GB, {len(volumes)} files)...{bcolors.ENDC}', flush=True)
    started = datetime.now()
    for f in sorted(volumes):
        try:
            with open(f, 'rb') as fh:
                while fh.read(32 * 1024 * 1024):
                    pass
        except OSError:
            continue
    elapsed = (datetime.now() - started).total_seconds()
    print(f'{bcolors.BLUE}BLAST DB cached in {elapsed:.0f}s.{bcolors.ENDC}', flush=True)
    _PREWARMED_DBS.add(blast_db)


# --- BLAST tuning ----------------------------------------------------------- #
# Assemblies often contain thousands of short SPAdes fragments; BLASTing them
# all is the main cause of slow/stalled runs and adds only noise. Keep only
# informative-length contigs as the query (with a longest-few fallback so we
# never search an empty file).
BLAST_MIN_CONTIG_LEN = int(os.environ.get('BLAST_MIN_CONTIG_LEN', '300'))
BLAST_FALLBACK_KEEP = 50
# 'dc-megablast' (discontiguous megablast) keeps cross-species sensitivity for
# viral/segmented targets while being far faster than blastn with a tiny word
# size. Override with the BLAST_TASK env var if a run needs a different mode.
BLAST_TASK = os.environ.get('BLAST_TASK', 'dc-megablast')
# Queries >= 50 kb (the split path: bacterial assemblies and their consensus,
# large DNA viruses) hit their own species, where megablast (the blastn default
# used here before dc-megablast) is the right tool; dc-megablast turned a single
# MTB consensus search from minutes into over an hour. BLAST_TASK_LARGE (or
# BLAST_TASK) overrides it.
BLAST_TASK_LARGE = os.environ.get('BLAST_TASK_LARGE') or os.environ.get('BLAST_TASK') or 'megablast'
# Hard ceiling so a pathological query/database can never hang the pipeline.
BLAST_TIMEOUT_SECONDS = int(os.environ.get('BLAST_TIMEOUT_SECONDS', '3600'))


def _prepare_blast_query(fasta, sample_name, min_len=BLAST_MIN_CONTIG_LEN,
                         fallback_keep=BLAST_FALLBACK_KEEP):
    """Return a query FASTA with short SPAdes micro-contigs removed.

    Keeps records >= min_len; if that removes everything, keeps the longest
    `fallback_keep` so we never BLAST nothing. Returns the original path
    unchanged when no filtering is needed (e.g. a small consensus FASTA).
    """
    try:
        records = list(SeqIO.parse(fasta, "fasta"))
    except Exception:
        return fasta
    if not records:
        return fasta
    kept = [r for r in records if len(r.seq) >= min_len]
    if not kept:
        kept = sorted(records, key=lambda r: len(r.seq), reverse=True)[:fallback_keep]
    if len(kept) == len(records):
        return fasta
    out = f'{sample_name}_blast_query.fasta'
    SeqIO.write(kept, out, "fasta")
    print(f'{bcolors.BLUE}BLAST query reduced to {len(kept)}/{len(records)} '
          f'contigs >= {min_len} bp (skipped {len(records) - len(kept)} short '
          f'fragments to speed up search){bcolors.ENDC}')
    return out


def _run_blast_cmd(cmd):
    """Run a blastn command with a timeout so it can never stall indefinitely.

    Returns True on success. On timeout or failure it prints a warning and
    returns False; downstream parsing already tolerates an empty/partial
    result file, so the pipeline degrades gracefully instead of hanging.
    """
    try:
        subprocess.run(cmd, shell=True, check=True, timeout=BLAST_TIMEOUT_SECONDS)
        return True
    except subprocess.TimeoutExpired:
        print(f'{bcolors.RED}BLAST timed out after {BLAST_TIMEOUT_SECONDS}s; '
              f'continuing with partial/empty results.{bcolors.ENDC}')
    except subprocess.CalledProcessError as e:
        print(f'{bcolors.RED}BLAST exited with code {e.returncode}; '
              f'continuing with partial/empty results.{bcolors.ENDC}')
    return False


# Columns appended to the BLAST outfmt when chimeric contigs are split, so each
# hit carries the query span it explains.
CHIMERA_COLUMNS = 'qstart qend qlen'


def _uncovered_parts(blastout_file, query_fasta, min_len=BLAST_MIN_CONTIG_LEN):
    """Return SeqRecords for contig stretches (>= min_len bp) that the contig's
    top BLAST hit does not cover.

    Segmented viruses (BTV/EHD, ISAV, ...) share conserved segment termini, so
    SPAdes can fuse two segments into one contig (e.g. BTV VP5 + VP1). BLAST then
    reports only the top hit per contig and the fused segment is never seen. The
    uncovered stretch is returned so it can be searched on its own. Requires
    CHIMERA_COLUMNS as the last three outfmt columns.
    """
    top = {}   # qseqid -> [top subject, [(start, end), ...], qlen]
    with open(blastout_file) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            try:
                qstart, qend, qlen = int(f[-3]), int(f[-2]), int(f[-1])
            except ValueError:
                continue
            hit = top.setdefault(f[0], [f[1], [], qlen])
            if hit[0] == f[1]:
                hit[1].append((min(qstart, qend), max(qstart, qend)))
    if not top:
        return []
    parts = []
    for rec in SeqIO.parse(query_fasta, 'fasta'):
        if rec.id not in top:
            continue
        _sacc, spans, qlen = top[rec.id]
        gaps, pos = [], 1
        for start, end in sorted(spans):
            if start > pos:
                gaps.append((pos, start - 1))
            pos = max(pos, end + 1)
        if pos <= qlen:
            gaps.append((pos, qlen))
        n = 0
        for start, end in gaps:
            seq = rec.seq[start - 1:end]
            if len(seq) < min_len or seq.upper().count('N') > len(seq) / 2:
                continue
            n += 1
            parts.append(SeqRecord(seq, id=f'{rec.id}_part{n}',
                                   description=f'{rec.id}:{start}-{end}'))
    return parts


class Blast_Fasta(Setup, bcolors):
    def __init__(self, FASTA=None, search=None, blast_out=None, format="6 qseqid sacc bitscore pident stitle", num_alignment=3, blast_db="nt", split_chimeras=False):
        """split_chimeras: after the search, re-BLAST any contig stretch its top
        hit does not explain (fused segments of a segmented virus) and add those
        hits to the results. Only used for the de novo assembly search."""
        Setup.__init__(self, FASTA=FASTA)
        sample_name = self.sample_name
        self.blast_db = blast_db
        blastout_file = f'{sample_name}_blast_out.txt'
        self.split_file_list = []  # Initialize split_file_list here
        split_chimeras = split_chimeras and not blast_out
        if split_chimeras:
            format = f'{format} {CHIMERA_COLUMNS}'
        part_lengths = {}   # chimera part id -> length; parent id -> length explained by its top hit

        # Detect HPC system
        self.hpc_system = self.detect_hpc_system()
        print(f'{bcolors.BLUE}Detected HPC system: {self.hpc_system}{bcolors.ENDC}')

        if blast_out:
            blastout_file = blast_out
        else:
            # Drop short micro-contigs before searching (biggest speed win).
            FASTA = _prepare_blast_query(FASTA, self.sample_name)
            fasta_size = sum([len(seq_record.seq) for seq_record in SeqIO.parse(FASTA, "fasta")])
            
            # Check if we're in a SLURM environment
            is_slurm = shutil.which('sbatch') is not None
            use_slurm = is_slurm and (os.path.isdir("/project") or os.path.isdir("/software/public/databases"))
            if not use_slurm:
                # blastn runs on THIS host — make sure the DB is served from RAM,
                # not one 4 KB disk seek at a time (see _prewarm_blast_db).
                _prewarm_blast_db(blast_db)

            if fasta_size < 50000:
                print(f'\n{bcolors.YELLOW}Running BLAST...{bcolors.ENDC}\n')
                if use_slurm:
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
                            f'blastn -query {FASTA} -db {blast_db} -task {BLAST_TASK} -out {self.sample_name}_blast_out.txt '
                            f'-outfmt "{format}" -num_alignments {num_alignment} -num_threads {self.cpus}'
                        )
                    os.system(f'sbatch -W ./batch.sh')
                else:
                    # Direct execution without SLURM (timeout-guarded).
                    _run_blast_cmd(
                        f'blastn -query {FASTA} -db {blast_db} -task {BLAST_TASK} -out {self.sample_name}_blast_out.txt '
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
                
                if use_slurm:
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
                            f'blastn -query ${{names[${{SLURM_ARRAY_TASK_ID}}]}}.fasta -db {blast_db} -task {BLAST_TASK_LARGE} '
                            f'-out ${{names[${{SLURM_ARRAY_TASK_ID}}]}}_blastout.txt -outfmt "{format}" '
                            f'-num_alignments {num_alignment} -num_threads {self.cpus}\n'
                            f'rm jobs\n'
                        )
                    os.system('sbatch -W ./batch.sh')
                else:
                    # Sequential processing without SLURM (timeout-guarded).
                    total_files = len(self.split_file_list)
                    for i, split_file in enumerate(self.split_file_list, 1):
                        print(f'  BLAST chunk {i}/{total_files}: {split_file}', flush=True)
                        _run_blast_cmd(
                            f'blastn -query {split_file} -db {blast_db} -task {BLAST_TASK_LARGE} '
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

        # If BLAST timed out or failed, its output file may be missing. Create an
        # empty one so parsing yields empty (not crashing) results and still sets
        # every attribute the pipeline reads downstream (e.g. blast_summary_file).
        if not os.path.exists(blastout_file):
            print(f'{bcolors.RED}No BLAST output produced; proceeding with empty '
                  f'results.{bcolors.ENDC}')
            open(blastout_file, 'a').close()

        # Fused segments: search each contig stretch its top hit leaves unexplained
        # and append those hits, so e.g. a VP5 segment fused onto VP1 is identified.
        if split_chimeras:
            parts = _uncovered_parts(blastout_file, FASTA)
            if parts:
                parts_fasta = f'{sample_name}_blast_chimera_parts.fasta'
                parts_out = f'{sample_name}_blast_chimera_parts_out.txt'
                SeqIO.write(parts, parts_fasta, 'fasta')
                print(f'{bcolors.YELLOW}{len(parts)} contig region(s) not explained by the contig\'s '
                      f'top BLAST hit (likely fused segments); searching them separately: '
                      f'{", ".join(p.description for p in parts)}{bcolors.ENDC}')
                _run_blast_cmd(
                    f'blastn -query {parts_fasta} -db {blast_db} -task {BLAST_TASK} -out {parts_out} '
                    f'-outfmt "{format}" -num_alignments {num_alignment} -num_threads {self.cpus}')
                if os.path.exists(parts_out):
                    with open(parts_out) as extra, open(blastout_file, 'a') as out:
                        out.write(extra.read())
                    os.remove(parts_out)
                os.remove(parts_fasta)
                # Credit each stretch to its own hit: the part gets its length and
                # the parent contig keeps only the length its top hit explains.
                parent_len = {rec.id: len(rec.seq) for rec in SeqIO.parse(FASTA, 'fasta')}
                for p in parts:
                    parent = p.description.rsplit(':', 1)[0]
                    part_lengths[p.id] = len(p.seq)
                    part_lengths[parent] = part_lengths.get(parent, parent_len[parent]) - len(p.seq)

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
                term = re.sub(r'[/.!@#$%^&*()+,"= ]', '_', term)
                SeqIO.write(found_record, f'{sample_name}_search_{term}.fasta', 'fasta')

            acc_frequency = []
            top_hit_acc_norm = []
            norm_dict = {}
            descriptions = {}
            # Contig lengths (read once, not once per contig); fused-segment parts
            # and their trimmed parents come from part_lengths.
            seq_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(FASTA, "fasta")}
            seq_lengths.update(part_lengths)
            for header, description in blast_dict.items():
                # Get accession frequencies
                acc_frequency.append(description[0][0])  # top hit, 1st item in hit
                acc_count = Counter()
                for acc in acc_frequency:
                    acc_count[acc] += 1
                # Get FASTA sizes per accessions
                seq_length = seq_lengths[header]
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
                    summary_list.append((f'{count:,}', f'{acc_count[acc]:,}', descriptions[acc], acc))
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

    def excel(self, excel_dict):
        basename = os.path.basename(self.blast_db)
        if not self.summary_list:
            excel_dict[f'Top BLAST {basename} Hit - 1'] = 'BLAST failed - No Results Output From BLAST search'
            return
        # summary_list is sorted ascending by nt count, so the top hits are last.
        # Each entry is (nt count, contigs, description, accession); the
        # description is a plain string (indexing it used to emit single letters).
        for rank, (count, contigs, description, acc) in enumerate(reversed(self.summary_list[-3:]), 1):
            excel_dict[f'Top BLAST {basename} Hit - {rank}'] = f'{count} {basename}, {contigs} contigs of {acc} {description}'


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