#!/usr/bin/env python3
"""
Coverage Graph Generator
Creates PNG coverage graphs from FASTA reference and FASTQ reads
"""

__version__ = "2.0.0"

import os
import sys
import re
import glob
import argparse
import textwrap
import subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO

import coverage_plot
from file_setup import resolve_cpus

class CoverageGraphGenerator:
    """Generate coverage graphs from FASTA references and FASTQ reads"""

    def __init__(self, fasta_file, fastq_r1, fastq_r2=None, debug=False):
        # Convert relative paths to absolute paths to avoid BWA/samtools path issues.
        # fastq_r2 is None for single-end (e.g. Nanopore / single FASTQ) runs.
        self.fasta_file = os.path.abspath(fasta_file)
        self.fastq_r1 = os.path.abspath(fastq_r1)
        self.fastq_r2 = os.path.abspath(fastq_r2) if fastq_r2 else None
        self.debug = debug
        self.cpus = resolve_cpus()
        self.output_files = []
        # Set to a human-readable reason if variant calling could not run (e.g.
        # freebayes missing / broken). Distinguishes "SNP calling failed" from the
        # very different "calling ran and found 0 SNPs" everywhere downstream.
        self.snp_call_error = None

        # Variables to store full paths of PNG files
        self.png_file_paths = {}  # Dictionary to store PNG file paths
        self.png_file_descriptions = {}  # Dictionary to store sequence descriptions for each PNG
        self.png_file_variables = []  # List of variable names

        # Validate input files
        self._validate_inputs()

        # Get sample name from FASTQ file
        self.sample_name = re.sub('[._].*', '', os.path.basename(fastq_r1))

        # Create clean display name without _# suffixes for graph titles
        self.display_name = re.sub(r'_\d+$', '', self.sample_name)

        if self.debug:
            print(f"Initialized CoverageGraphGenerator:")
            print(f"  FASTA: {fasta_file}")
            print(f"  FASTQ R1: {fastq_r1}")
            print(f"  FASTQ R2: {fastq_r2}")
            print(f"  Sample name (for files): {self.sample_name}")
            print(f"  Display name (for graphs): {self.display_name}")
            print(f"  Sample name: {self.sample_name}")
            print(f"  CPUs: {self.cpus}")

    def _validate_inputs(self):
        """Validate that input files exist and are readable"""
        for file_path, file_type in [(self.fasta_file, "FASTA"),
                                    (self.fastq_r1, "FASTQ R1"),
                                    (self.fastq_r2, "FASTQ R2")]:
            if file_path is None:  # R2 is optional (single-end run)
                continue
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"{file_type} file not found: {file_path}")
            if not os.path.isfile(file_path):
                raise ValueError(f"{file_type} path is not a file: {file_path}")

    def _process_reference(self):
        """Process reference FASTA file for alignment"""
        if self.debug:
            print(f"Processing reference: {self.fasta_file}")

        # Clean up sequence IDs and ensure uniqueness
        record_list = []
        used_ids = set()
        counter = 1
        id_to_description = {}  # Map final cleaned IDs to original descriptions

        with open(self.fasta_file) as handle:
            records = SeqIO.parse(handle, "fasta")
            for rec in records:
                original_description = rec.description  # Store original full description

                # Clean the ID (replace all problematic characters with underscores)
                clean_id = rec.id
                # Replace spaces, pipes, parentheses, brackets, and other shell/bioinformatics problematic chars
                problematic_chars = [' ', ':', '|', '(', ')', '[', ']', ';', '&', '$', '*', '?', '/', '\\', '"', "'", '#', '%', '<', '>', '=', '+', '!', '@', '^', '~', '`']
                for char in problematic_chars:
                    clean_id = clean_id.replace(char, '_')

                # Remove any double underscores created by consecutive special characters
                while '__' in clean_id:
                    clean_id = clean_id.replace('__', '_')

                # Remove leading/trailing underscores
                clean_id = clean_id.strip('_')

                # Ensure uniqueness by adding counter if needed
                if clean_id in used_ids:
                    unique_id = f"{clean_id}_{counter}"
                    while unique_id in used_ids:
                        counter += 1
                        unique_id = f"{clean_id}_{counter}"
                    rec.id = unique_id
                    counter += 1
                else:
                    rec.id = clean_id

                # Map final ID to original description
                id_to_description[rec.id] = original_description

                used_ids.add(rec.id)
                record_list.append(rec)

                if self.debug:
                    print(f"  Processed sequence: {rec.id} -> {original_description[:60]}...")

        # Store the mapping for later use
        self._id_to_description = id_to_description

        # Create temporary reference file with cleaned IDs in current directory
        temp_reference = f'{os.path.basename(self.fasta_file)}.temp.fasta'
        SeqIO.write(record_list, temp_reference, "fasta")

        if self.debug:
            print(f"Created temporary reference: {temp_reference}")

        return temp_reference

    def _align_reads(self, temp_reference):
        """Align reads to reference using BWA"""
        if self.debug:
            print("Starting read alignment with BWA")

        # File names for alignment pipeline (keep output files relative to current directory)
        sam_file = f"{self.sample_name}.sam"
        bam_file = f"{self.sample_name}.bam"
        sorted_bam = f"{self.sample_name}.sorted.bam"

        try:
            # Index reference
            print("Indexing reference...")
            if self.debug:
                print(f"  Command: bwa index {temp_reference}")
                print(f"  Reference file exists: {os.path.exists(temp_reference)}")
                if os.path.exists(temp_reference):
                    print(f"  Reference file size: {os.path.getsize(temp_reference)} bytes")

            result = subprocess.run(["bwa", "index", temp_reference],
                                  capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"BWA index failed:\n  Command: bwa index {temp_reference}\n  Return code: {result.returncode}\n  Stderr: {result.stderr}\n  Stdout: {result.stdout}")

            if self.debug:
                print(f"  BWA index completed successfully")

            # Align reads (R2 only when the run is paired-end)
            bwa_cmd = ["bwa", "mem", "-t", str(self.cpus), temp_reference, self.fastq_r1]
            if self.fastq_r2:
                bwa_cmd.append(self.fastq_r2)
            print("Aligning reads to reference...")
            if self.debug:
                print(f"  Command: {' '.join(bwa_cmd)}")
                print(f"  FASTQ R1 exists: {os.path.exists(self.fastq_r1)}")
                if os.path.exists(self.fastq_r1):
                    print(f"  FASTQ R1 size: {os.path.getsize(self.fastq_r1)} bytes")
                if self.fastq_r2:
                    print(f"  FASTQ R2 exists: {os.path.exists(self.fastq_r2)}")
                    if os.path.exists(self.fastq_r2):
                        print(f"  FASTQ R2 size: {os.path.getsize(self.fastq_r2)} bytes")

            with open(sam_file, 'w') as sam_out:
                result = subprocess.run(bwa_cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                raise RuntimeError(f"BWA mem alignment failed:\n  Command: {' '.join(bwa_cmd)}\n  Return code: {result.returncode}\n  Stderr: {result.stderr}")

            # Check SAM file was created properly
            if not os.path.exists(sam_file):
                raise RuntimeError(f"SAM file was not created: {sam_file}")

            sam_size = os.path.getsize(sam_file)
            if sam_size == 0:
                raise RuntimeError(f"SAM file is empty: {sam_file}")

            if self.debug:
                print(f"  BWA mem alignment completed successfully")
                print(f"  SAM file created: {sam_file}")
                print(f"  SAM file size: {sam_size} bytes")

            # Convert SAM to BAM
            print("Converting SAM to BAM...")
            if self.debug:
                print(f"  Command: samtools view -bS {sam_file}")

            with open(bam_file, 'wb') as bam_out:
                result = subprocess.run(["samtools", "view", "-bS", sam_file],
                                      stdout=bam_out, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                raise RuntimeError(f"samtools view (SAM to BAM conversion) failed:\n  Command: samtools view -bS {sam_file}\n  Return code: {result.returncode}\n  Stderr: {result.stderr}\n  SAM file size: {sam_size} bytes")

            # Sort BAM file
            print("Sorting BAM file...")
            subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam],
                         check=True)

            # Index BAM file
            print("Indexing BAM file...")
            subprocess.run(["samtools", "index", sorted_bam], check=True)

            # Print alignment statistics
            mapped_reads = subprocess.check_output(
                ["samtools", "view", "-c", "-F", "4", sorted_bam]
            )
            print(f"Mapped reads: {mapped_reads.decode().strip()}")

            if self.debug:
                bam_stats = subprocess.check_output(["samtools", "flagstat", sorted_bam])
                print(f"BAM statistics:\n{bam_stats.decode()}")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Alignment failed: {e}")

        # Clean up intermediate files
        if not self.debug:
            for f in [sam_file, bam_file]:
                if os.path.exists(f):
                    os.remove(f)

        return sorted_bam

    def _organize_alignment_files(self, sorted_bam, vcf_file, temp_reference):
        """Organize BAM, VCF, and index files into a dedicated folder"""
        alignment_dir = "coverage_graph_alignment"

        if self.debug:
            print(f"Organizing alignment files into {alignment_dir}/")

        # Create alignment directory
        os.makedirs(alignment_dir, exist_ok=True)

        # List of files to move
        files_to_move = []

        # BAM file and its index
        if os.path.exists(sorted_bam):
            files_to_move.append(sorted_bam)
            bam_index = sorted_bam + '.bai'
            if os.path.exists(bam_index):
                files_to_move.append(bam_index)

        # VCF file
        if vcf_file and os.path.exists(vcf_file):
            files_to_move.append(vcf_file)

        # Reference index files
        fai_file = temp_reference + '.fai'
        if os.path.exists(fai_file):
            files_to_move.append(fai_file)

        # BWA index files
        for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            index_file = temp_reference + ext
            if os.path.exists(index_file):
                files_to_move.append(index_file)

        # Move files to alignment directory
        moved_files = []
        for file_path in files_to_move:
            if os.path.exists(file_path):
                dest_path = os.path.join(alignment_dir, os.path.basename(file_path))
                try:
                    # Use shutil.move to handle cross-filesystem moves
                    import shutil
                    shutil.move(file_path, dest_path)
                    moved_files.append(dest_path)
                    if self.debug:
                        print(f"Moved: {file_path} → {dest_path}")
                except Exception as e:
                    print(f"Warning: Could not move {file_path}: {e}")

        print(f"Organized {len(moved_files)} alignment files into {alignment_dir}/")

        # Update paths for return
        new_sorted_bam = os.path.join(alignment_dir, os.path.basename(sorted_bam)) if sorted_bam in files_to_move else sorted_bam
        new_vcf_file = os.path.join(alignment_dir, os.path.basename(vcf_file)) if vcf_file and vcf_file in files_to_move else vcf_file

        return new_sorted_bam, new_vcf_file, alignment_dir

    def _generate_vcf(self, sorted_bam, temp_reference):
        """Call variants (reads vs the reference) with freebayes for the report's
        SNP track. freebayes is the pipeline's variant caller everywhere else, so
        the report has no separate/duplicate caller dependency (no bcftools)."""
        if self.debug:
            print("Generating VCF file for SNP analysis (freebayes)...")

        vcf_file = f"{self.sample_name}_variants.vcf"

        try:
            # freebayes needs a .fai index on the reference; create it up front.
            subprocess.run(["samtools", "faidx", temp_reference],
                           check=True, stderr=subprocess.PIPE)

            print("Calling variants with freebayes...")
            # Same caller/quality gating idea as the main reference-guided step:
            # min mapping quality 30, min base quality 20, spec-compliant output.
            fb_cmd = [
                "freebayes", "-f", temp_reference,
                "-m", "30", "-q", "20", "--strict-vcf",
                sorted_bam,
            ]
            with open(vcf_file, "w") as out:
                proc = subprocess.run(fb_cmd, stdout=out, stderr=subprocess.PIPE)
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, 'freebayes',
                                                    stderr=proc.stderr)

            if self.debug:
                print(f"VCF file created: {vcf_file}")

        except subprocess.CalledProcessError as e:
            self.snp_call_error = self._describe_caller_error(
                'freebayes',
                e.stderr.decode('utf-8', 'replace') if getattr(e, 'stderr', None) else str(e))
            print(f"Warning: variant calling (freebayes) failed: {self.snp_call_error}")
            print("SNP calling will be reported as unavailable in the report")
            return None
        except FileNotFoundError:
            self.snp_call_error = ("freebayes not found on PATH — install it "
                                   "(conda install -c bioconda -c conda-forge freebayes)")
            print(f"Warning: {self.snp_call_error}")
            return None

        return vcf_file

    @staticmethod
    def _describe_caller_error(tool, stderr_text):
        """Turn raw variant-caller stderr into a short, actionable one-liner."""
        text = (stderr_text or '').strip()
        low = text.lower()
        # Pull out the most useful single line (the shared-library / error line).
        line = next((l.strip() for l in text.splitlines()
                     if 'error' in l.lower() or '.so' in l.lower()),
                    text.splitlines()[0] if text.splitlines() else text)
        if 'error while loading shared libraries' in low or '.so' in low:
            return (f"{tool} cannot start — {line}. Its conda environment is broken "
                    "(a shared library is missing). Reinstall it: "
                    f"conda install -c bioconda -c conda-forge {tool}")
        return f"{tool} failed: {line[:300]}"

    def _add_zero_coverage_to_vcf(self, sorted_bam, vcf_file, temp_reference):
        """Add zero coverage positions to VCF file"""
        if not vcf_file or not os.path.exists(vcf_file):
            return vcf_file

        if self.debug:
            print("Adding zero coverage positions to VCF...")

        try:
            # Get coverage data using samtools depth
            coverage_dict = {}
            result = subprocess.run(['samtools', 'depth', sorted_bam],
                                  capture_output=True, text=True, check=True)

            for line in result.stdout.splitlines():
                parts = line.split('\t')
                if len(parts) >= 3:
                    chrom, position, depth = parts[0], parts[1], parts[2]
                    coverage_dict[chrom + "-" + position] = int(depth)

            # Create complete position set (all positions from 1 to sequence length)
            zero_positions = []
            for record in SeqIO.parse(temp_reference, "fasta"):
                chrom = record.id
                total_len = len(record.seq)

                for pos in range(1, total_len + 1):
                    key = f"{chrom}-{pos}"
                    if key not in coverage_dict:  # Position has zero coverage
                        zero_positions.append({
                            'CHROM': chrom,
                            'POS': pos,
                            'ID': '.',
                            'REF': 'N',
                            'ALT': '.',
                            'QUAL': '.',
                            'FILTER': '.',
                            'INFO': '.',
                            'FORMAT': 'GT',
                            'Sample': './.'
                        })

            if not zero_positions:
                if self.debug:
                    print("No zero coverage positions found")
                return vcf_file

            print(f"Found {len(zero_positions)} zero coverage positions")

            # Read existing VCF data
            vcf_data = []
            header_lines = []

            with open(vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        header_lines.append(line.strip())
                    else:
                        parts = line.strip().split('\t')
                        if len(parts) >= 10:
                            vcf_data.append({
                                'CHROM': parts[0],
                                'POS': int(parts[1]),
                                'ID': parts[2],
                                'REF': parts[3],
                                'ALT': parts[4],
                                'QUAL': parts[5],
                                'FILTER': parts[6],
                                'INFO': parts[7],
                                'FORMAT': parts[8],
                                'Sample': parts[9]
                            })

            # Combine and sort all data
            all_data = vcf_data + zero_positions
            all_data.sort(key=lambda x: (x['CHROM'], x['POS']))

            # Write enhanced VCF file
            zero_coverage_vcf = f"{self.sample_name}_zc.vcf"
            with open(zero_coverage_vcf, 'w') as f:
                # Write header
                for header in header_lines:
                    f.write(header + '\n')

                # Write all data
                for entry in all_data:
                    f.write(f"{entry['CHROM']}\t{entry['POS']}\t{entry['ID']}\t"
                           f"{entry['REF']}\t{entry['ALT']}\t{entry['QUAL']}\t"
                           f"{entry['FILTER']}\t{entry['INFO']}\t{entry['FORMAT']}\t"
                           f"{entry['Sample']}\n")

            if self.debug:
                print(f"Created zero coverage VCF: {zero_coverage_vcf}")

            return zero_coverage_vcf

        except subprocess.CalledProcessError as e:
            print(f"Warning: Error getting coverage data: {e}")
            return vcf_file
        except Exception as e:
            print(f"Warning: Error adding zero coverage positions: {e}")
            return vcf_file

    def _parse_vcf_file(self, vcf_file):
        """Parse VCF file to extract SNPs and no-coverage regions"""
        if not vcf_file or not os.path.exists(vcf_file):
            return {}, {}

        snps = defaultdict(list)
        no_coverage_regions = defaultdict(list)

        if self.debug:
            print(f"Parsing VCF file: {vcf_file}")

        snp_count = 0
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 10:
                        continue

                    chrom, pos, id_val, ref, alt, qual, filter_val, info, format_val, sample = parts[:10]
                    pos = int(pos)

                    # Check for real SNPs (not indels, with quality score)
                    if (alt != '.' and ref != 'N' and qual != '.' and
                        len(ref) == 1 and len(alt) == 1):
                        try:
                            qual_score = float(qual)
                            if qual_score > coverage_plot.GOOD_SNP_MIN_QUAL:  # keep only high-confidence SNPs
                                snps[chrom].append({
                                    'pos': pos,
                                    'ref': ref,
                                    'alt': alt,
                                    'qual': qual_score
                                })
                                snp_count += 1
                        except ValueError:
                            pass

                    # Check for no-coverage regions (REF=N, ALT=., GT=./.)
                    elif ref == 'N' and alt == '.' and './.' in sample:
                        # Group consecutive no-coverage positions into regions
                        if no_coverage_regions[chrom] and no_coverage_regions[chrom][-1][1] == pos - 1:
                            # Extend the last region
                            no_coverage_regions[chrom][-1] = (no_coverage_regions[chrom][-1][0], pos)
                        else:
                            # Start a new region
                            no_coverage_regions[chrom].append((pos, pos))

            if self.debug:
                print(f"Found {snp_count} SNPs across {len(snps)} sequences")
                total_no_cov_regions = sum(len(regions) for regions in no_coverage_regions.values())
                if total_no_cov_regions > 0:
                    print(f"Found {total_no_cov_regions} zero coverage regions")

        except Exception as e:
            print(f"Warning: Error parsing VCF file: {e}")
            return {}, {}

        return snps, no_coverage_regions

    def _get_reference_info(self, temp_reference):
        """Get reference sequence information"""
        ref_info = {}

        # Get info from temporary reference
        for record in SeqIO.parse(temp_reference, "fasta"):
            ref_info[record.id] = len(record.seq)

        # Use the mapping created during reference processing
        original_headers = self._id_to_description.copy()

        if self.debug:
            print(f"Reference info mapping:")
            for ref_id in ref_info.keys():
                print(f"  {ref_id} -> {original_headers.get(ref_id, 'NOT FOUND')[:80]}...")

        return ref_info, original_headers

    def _calculate_alignment_stats(self, sorted_bam, ref_info, original_headers):
        """Calculate alignment statistics for each reference sequence"""
        if self.debug:
            print("Calculating alignment statistics...")

        stats = {}
        bam = pysam.AlignmentFile(sorted_bam, "rb")

        for ref_id, ref_len in ref_info.items():
            # Initialize coverage array
            coverage_array = np.zeros(ref_len)

            try:
                # Calculate coverage for this reference
                for pileup_column in bam.pileup(ref_id):
                    if pileup_column.pos < ref_len:
                        coverage_array[pileup_column.pos] = pileup_column.n
            except Exception as e:
                if self.debug:
                    print(f"Warning: Error processing reference {ref_id}: {str(e)}")
                continue

            # Calculate statistics
            covered_bases = np.count_nonzero(coverage_array)
            mean_coverage = np.mean(coverage_array) if ref_len > 0 else 0.0
            percent_covered = (covered_bases / ref_len) * 100 if ref_len > 0 else 0

            stats[ref_id] = {
                'header': original_headers.get(ref_id, ref_id),
                'length': ref_len,
                'mean_coverage': mean_coverage,
                'percent_covered': percent_covered,
                'coverage_array': coverage_array
            }

            if self.debug:
                print(f"Stats for {ref_id}:")
                print(f"  Length: {ref_len}")
                print(f"  Mean coverage: {mean_coverage:.1f}X")
                print(f"  Percent covered: {percent_covered:.1f}%")

        bam.close()
        return stats

    def _create_individual_sequence_graphs(self, alignment_stats, snps=None, no_coverage_regions=None):
        """Render one coverage+variant figure per reference and record its path
        and description for the report.

        Plotting is delegated to coverage_plot.static_png — the single coverage
        renderer shared with the report — which downsamples the depth line and
        switches to a binned SNP-density track on large genomes, so big
        references never look cluttered with hundreds of individual SNP markers.
        """
        snps = snps or {}
        no_coverage_regions = no_coverage_regions or {}
        individual_graphs = []

        # Coverage graphs are produced for ALL references, even low-coverage ones.
        for ref_id, stats in alignment_stats.items():
            if self.debug:
                print(f"Creating coverage graph for {ref_id}: "
                      f"{stats['mean_coverage']:.1f}X coverage, "
                      f"{stats['percent_covered']:.1f}% covered")

            output_png = coverage_plot.static_png(
                ref_id, stats, snps.get(ref_id, []),
                no_coverage_regions.get(ref_id, []), out_dir='.')

            safe_ref_id = re.sub(r'[^\w\-_.]', '_', ref_id)
            variable_name = f"png_{safe_ref_id}_path"
            self.png_file_paths[variable_name] = output_png
            self.png_file_descriptions[variable_name] = stats['header']  # original FASTA header
            self.png_file_variables.append(variable_name)
            setattr(self, variable_name, output_png)

            individual_graphs.append(output_png)
            self.output_files.append(output_png)
            if self.debug:
                print(f"Created coverage graph: {output_png}")

        return individual_graphs

    def generate_coverage_graphs(self):
        """Main method to generate all coverage graphs with SNP analysis"""
        print(f"\n=== Coverage Graph Generation Started ===")
        print(f"Sample: {self.sample_name}")
        print(f"Reference: {os.path.basename(self.fasta_file)}")
        print("=" * 50)

        temp_reference = None
        sorted_bam = None
        vcf_file = None

        try:
            # Process reference
            temp_reference = self._process_reference()

            # Get reference information
            ref_info, original_headers = self._get_reference_info(temp_reference)
            print(f"Found {len(ref_info)} reference sequences")

            # Align reads
            sorted_bam = self._align_reads(temp_reference)

            # Generate VCF for SNP analysis
            vcf_file = self._generate_vcf(sorted_bam, temp_reference)

            # Add zero coverage positions to VCF
            if vcf_file:
                vcf_file = self._add_zero_coverage_to_vcf(sorted_bam, vcf_file, temp_reference)

            # Organize alignment files into dedicated folder
            organized_bam, organized_vcf, alignment_dir = self._organize_alignment_files(sorted_bam, vcf_file, temp_reference)

            # Parse VCF to extract SNPs and coverage gaps (use organized VCF path)
            snps, no_coverage_regions = self._parse_vcf_file(organized_vcf)

            # Calculate alignment statistics (use organized BAM path)
            alignment_stats = self._calculate_alignment_stats(organized_bam, ref_info, original_headers)

            # Record the high-quality SNP count on each reference's stats. This is
            # the authoritative count the report displays (Alignment Statistics
            # table + graph title): it lives in the SAME dict as length/coverage,
            # so it can never fall out of sync with a separate SNP structure.
            # If variant calling could NOT run, store None (= "not called") rather
            # than 0, so the report never presents a failure as a real zero.
            for _rid, _stats in alignment_stats.items():
                _stats['snp_count'] = None if self.snp_call_error else len(snps.get(_rid, []))

            # Store for convenience function access (consumed by the HTML/PDF report)
            self._last_alignment_stats = alignment_stats
            self._last_snps = snps
            self._last_no_coverage = no_coverage_regions
            self._last_snp_call_error = self.snp_call_error

            # Create individual sequence graphs with SNPs
            individual_graphs = self._create_individual_sequence_graphs(alignment_stats, snps, no_coverage_regions)

            # Print summary
            print("\n=== Generation Complete ===")
            print(f"Total output files: {len(self.output_files)}")
            for output_file in self.output_files:
                size_kb = os.path.getsize(output_file) / 1024
                print(f"  - {output_file} ({size_kb:.1f} KB)")

            # Print PNG file variables
            if self.png_file_variables:
                print(f"\n=== PNG File Variables ===")
                print(f"Variable names list: {self.png_file_variables}")
                print("Full paths:")
                for var_name in self.png_file_variables:
                    full_path = self.png_file_paths[var_name]
                    print(f"  {var_name}: {full_path}")

            # Print alignment and SNP summary
            print("\n=== Alignment and SNP Summary ===")
            total_snps = 0
            total_zero_coverage_regions = 0
            for ref_id, stats in alignment_stats.items():
                ref_snps = len(snps.get(ref_id, []))
                ref_zero_regions = len(no_coverage_regions.get(ref_id, []))
                total_snps += ref_snps
                total_zero_coverage_regions += ref_zero_regions
                print(f"{ref_id}: {stats['mean_coverage']:.1f}X coverage, "
                     f"{stats['percent_covered']:.1f}% genome covered, "
                     f"{ref_snps} SNPs, {ref_zero_regions} zero coverage regions")

            if total_snps > 0:
                print(f"\nTotal SNPs found: {total_snps}")
                # Show nucleotide distribution
                nuc_counts = defaultdict(int)
                if snps:
                    for seq_snps in snps.values():
                        for snp in seq_snps:
                            nuc_counts[snp['alt']] += 1
                if nuc_counts:
                    print("SNP nucleotide distribution:", dict(nuc_counts))
            else:
                print("\nNo SNPs detected (or freebayes may not be installed)")

            if total_zero_coverage_regions > 0:
                print(f"Total zero coverage regions: {total_zero_coverage_regions}")

            # Print alignment folder info
            if 'alignment_dir' in locals():
                print(f"\nAlignment files organized in: {alignment_dir}/")
                alignment_files = os.listdir(alignment_dir) if os.path.exists(alignment_dir) else []
                for af in alignment_files:
                    af_path = os.path.join(alignment_dir, af)
                    if os.path.isfile(af_path):
                        size_mb = os.path.getsize(af_path) / (1024 * 1024)
                        print(f"  - {af} ({size_mb:.1f} MB)")

        except Exception as e:
            print(f"Error during coverage graph generation: {e}")
            raise
        finally:
            # Cleanup only temporary files (not organized alignment files)
            self._cleanup_files(temp_reference)

    def _cleanup_files(self, temp_reference):
        """Clean up temporary files only (organized alignment files are preserved)"""
        if not self.debug:
            # Only clean up truly temporary files that weren't organized
            cleanup_patterns = [
                '*.sam',  # SAM files are always temporary
            ]

            # Remove temporary reference file
            files_to_remove = []
            if temp_reference and os.path.exists(temp_reference):
                files_to_remove.append(temp_reference)

            # Remove pattern-based temporary files (only in current directory, not in alignment folder)
            for pattern in cleanup_patterns:
                for f in glob.glob(pattern):
                    # Only remove if not in alignment directory
                    if not f.startswith('coverage_graph_alignment/'):
                        try:
                            os.remove(f)
                            if self.debug:
                                print(f"Removed temporary file: {f}")
                        except OSError:
                            pass

            # Remove specific temporary files
            for f in files_to_remove:
                try:
                    if os.path.exists(f):
                        os.remove(f)
                        if self.debug:
                            print(f"Removed temporary file: {f}")
                except OSError:
                    pass

    def get_png_file_paths(self):
        """
        Get dictionary of PNG file paths with variable names as keys

        Returns:
            dict: Dictionary with variable names as keys and full file paths as values
        """
        return self.png_file_paths.copy()

    def get_png_variable_names(self):
        """
        Get list of variable names for PNG file paths

        Returns:
            list: List of variable names that store PNG file paths
        """
        return self.png_file_variables.copy()


# Convenience function for importing this script
def generate_coverage_graphs_with_snps(fasta_file, fastq_r1, fastq_r2, output_dir=None, debug=False):
    """
    Convenience function for generating coverage graphs with SNP visualization

    Args:
        fasta_file (str): Path to reference FASTA file
        fastq_r1 (str): Path to R1 FASTQ file
        fastq_r2 (str): Path to R2 FASTQ file
        output_dir (str, optional): Output directory. Defaults to current directory.
        debug (bool, optional): Enable debug mode. Defaults to False.

    Returns:
        dict: Dictionary containing output information
            - 'output_files': List of generated PNG files
            - 'alignment_dir': Path to alignment files directory
            - 'sample_name': Sample name used for files
            - 'alignment_stats': Dictionary of alignment statistics
            - 'snp_stats': Dictionary of SNP statistics

    Example:
        >>> from coverage_graph_generator import generate_coverage_graphs_with_snps
        >>> result = generate_coverage_graphs_with_snps(
        ...     fasta_file="reference.fasta",
        ...     fastq_r1="sample_R1.fastq.gz",
        ...     fastq_r2="sample_R2.fastq.gz",
        ...     output_dir="results/",
        ...     debug=True
        ... )
        >>> print(f"Generated {len(result['output_files'])} files")
        >>> print(f"Found {result['snp_stats']['total_snps']} SNPs")
    """
    # Change to output directory if specified
    original_dir = None
    if output_dir:
        original_dir = os.getcwd()
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

    try:
        # Create generator and run analysis
        generator = CoverageGraphGenerator(
            fasta_file=fasta_file,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            debug=debug
        )

        generator.generate_coverage_graphs()

        # Collect SNP statistics
        snp_stats = {
            'total_snps': 0,
            'nucleotide_distribution': {},
            'sequences_with_snps': 0
        }

        # Calculate SNP stats if available
        if hasattr(generator, '_last_snps') and generator._last_snps:
            nuc_counts = defaultdict(int)
            total_snps = 0
            for seq_snps in generator._last_snps.values():
                if seq_snps:
                    snp_stats['sequences_with_snps'] += 1
                for snp in seq_snps:
                    nuc_counts[snp['alt']] += 1
                    total_snps += 1

            snp_stats['total_snps'] = total_snps
            snp_stats['nucleotide_distribution'] = dict(nuc_counts)

        return {
            'output_files': generator.output_files,
            'alignment_dir': 'coverage_graph_alignment',
            'sample_name': generator.sample_name,
            'alignment_stats': getattr(generator, '_last_alignment_stats', {}),
            'snp_stats': snp_stats,
            'png_file_paths': generator.get_png_file_paths(),
            'png_variable_names': generator.get_png_variable_names()
        }

    finally:
        # Return to original directory if we changed it
        if original_dir:
            os.chdir(original_dir)


def main():
    parser = argparse.ArgumentParser(
        prog='Coverage Graph Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Coverage Graph Generator with SNP Visualization
        ===============================================

        Generate PNG coverage graphs from FASTA reference sequences and FASTQ reads
        with SNP (Single Nucleotide Polymorphism) overlay visualization.

        This script:
        1. Aligns FASTQ reads to FASTA reference sequences using BWA
        2. Calls variants using freebayes to identify SNPs
        3. Identifies zero coverage positions and adds them to VCF
        4. Calculates coverage statistics for each reference sequence
        5. Creates individual coverage graphs as PNG files (one per reference sequence)
        6. Overlays SNP positions and zero coverage regions with visualization

        Output files:
        - {sample_name}_{sequence_id}_coverage_snps.png: Individual sequence plots with SNPs
        - {sample_name}_zc.vcf: Enhanced VCF with zero coverage positions

        SNP and Coverage Visualization:
        - SNPs shown as colored triangular markers on coverage curves
        - Colors: A=green, T=red, G=orange, C=blue
        - Individual plots include variant track showing zero coverage regions (gray rectangles)
        - Only high-confidence SNPs (QUAL > coverage_plot.GOOD_SNP_MIN_QUAL) are kept
        - Large genomes show a binned SNP-density track; short genomes show individual markers

        Requirements:
        - BWA (Burrows-Wheeler Aligner)
        - SAMtools
        - freebayes (for variant calling)
        - Python packages: matplotlib, pandas, numpy, pysam, biopython
        ''')
    )

    # Required arguments
    parser.add_argument('-f', '--fasta',
                       required=True,
                       help='Reference FASTA file containing sequences to align against')
    parser.add_argument('-r1', '--fastq1',
                       required=True,
                       help='Input R1 FASTQ file (can be gzipped)')
    parser.add_argument('-r2', '--fastq2',
                       required=False, default=None,
                       help='Input R2 FASTQ file (can be gzipped); omit for single-end reads')

    # Optional arguments
    parser.add_argument('-d', '--debug',
                       action='store_true',
                       help='Enable debug mode (keep intermediate files, verbose output)')
    parser.add_argument('-o', '--output-dir',
                       help='Output directory (default: current directory)')
    parser.add_argument('-v', '--version',
                       action='version',
                       version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    # Print configuration
    print(f"\nCoverage Graph Generator v{__version__}")
    print(f"Reference FASTA: {args.fasta}")
    print(f"FASTQ R1: {args.fastq1}")
    print(f"FASTQ R2: {args.fastq2}")
    print(f"Debug mode: {args.debug}")

    # Change to output directory if specified
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        original_dir = os.getcwd()
        os.chdir(args.output_dir)
        print(f"Output directory: {args.output_dir}")
    else:
        original_dir = None
        print(f"Output directory: {os.getcwd()}")

    try:
        # Create generator and run analysis
        generator = CoverageGraphGenerator(
            fasta_file=args.fasta,
            fastq_r1=args.fastq1,
            fastq_r2=args.fastq2,
            debug=args.debug
        )

        generator.generate_coverage_graphs()

        print(f"\n✓ Coverage graph generation completed successfully!")
        print(f"Generated {len(generator.output_files)} PNG files")

    except KeyboardInterrupt:
        print("\n⚠️ Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)
    finally:
        # Return to original directory if we changed it
        if original_dir:
            os.chdir(original_dir)


if __name__ == "__main__":
    main()