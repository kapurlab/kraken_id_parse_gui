#!/usr/bin/env python3

"""
BTV Reference Database Builder
==============================
Builds improved VP2 and VP5 reference databases for BTV serotyping
using validated samples with known serotypes.
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
import argparse
import textwrap

def extract_segment_sequences(fasta_file, sample_id, serotype):
    """Extract VP2 and VP5 sequences from a reference-guided FASTA file."""
    vp2_seq = None
    vp5_seq = None

    print(f"Processing {sample_id} (expected serotype: {serotype})")

    for record in SeqIO.parse(fasta_file, 'fasta'):
        header_lower = f"{record.id} {record.description}".lower()
        print(f"  Found sequence: {record.id}")
        print(f"    Description: {record.description}")

        # Check for VP2 (segment 2)
        if any(pattern in header_lower for pattern in [
            'segment 2', 'segment_2', '_segment2', 'vp2', 'viral protein 2'
        ]):
            print(f"    → Identified as VP2 (segment 2)")
            vp2_seq = record

        # Check for VP5 (segment 6)
        elif any(pattern in header_lower for pattern in [
            'segment 6', 'segment_6', '_segment6', 'vp5', 'viral protein 5'
        ]):
            print(f"    → Identified as VP5 (segment 6)")
            vp5_seq = record
        else:
            print(f"    → Not VP2 or VP5, skipping")

    return vp2_seq, vp5_seq

def create_reference_header(original_record, sample_id, serotype, protein_type):
    """Create a properly formatted reference database header."""
    # Format: >SampleID|Serotype|Description
    new_header = f"{sample_id}|{serotype}|Bluetongue virus {serotype[-2:]} sample {sample_id} {protein_type} gene, complete cds"

    return new_header

def build_reference_databases(input_dir, serotype_mapping, output_dir):
    """Build VP2 and VP5 reference databases from validated samples."""

    vp2_sequences = []
    vp5_sequences = []

    print("=" * 60)
    print("BTV REFERENCE DATABASE BUILDER")
    print("=" * 60)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Known serotypes: {len(serotype_mapping)} samples")
    print()

    # Process each sample
    for sample_id, serotype in serotype_mapping.items():
        print(f"Looking for sample: {sample_id}")

        # Find matching FASTA files
        fasta_files = []
        for pattern in [f"{sample_id}*.fasta", f"{sample_id}*reference_guided.fasta"]:
            matches = list(Path(input_dir).glob(pattern))
            fasta_files.extend(matches)

        if not fasta_files:
            print(f"  ⚠️  No FASTA files found for {sample_id}")
            continue

        # Use the first matching file
        fasta_file = fasta_files[0]
        print(f"  Using file: {fasta_file.name}")

        # Extract VP2 and VP5 sequences
        vp2_seq, vp5_seq = extract_segment_sequences(fasta_file, sample_id, serotype)

        # Process VP2 if found
        if vp2_seq:
            # Check sequence quality
            n_count = str(vp2_seq.seq).upper().count('N')
            seq_length = len(vp2_seq.seq)
            n_percent = (n_count / seq_length) * 100

            if n_percent > 50:
                print(f"    ⚠️  VP2 sequence has {n_percent:.1f}% N's - may be low quality")
            else:
                print(f"    ✅ VP2 sequence quality good ({n_percent:.1f}% N's)")

            # Create new record with proper header
            new_header = create_reference_header(vp2_seq, sample_id, serotype, "VP2")
            new_record = vp2_seq
            new_record.id = new_header
            new_record.description = ""
            vp2_sequences.append(new_record)
        else:
            print(f"    ❌ No VP2 sequence found")

        # Process VP5 if found
        if vp5_seq:
            # Check sequence quality
            n_count = str(vp5_seq.seq).upper().count('N')
            seq_length = len(vp5_seq.seq)
            n_percent = (n_count / seq_length) * 100

            if n_percent > 50:
                print(f"    ⚠️  VP5 sequence has {n_percent:.1f}% N's - may be low quality")
            else:
                print(f"    ✅ VP5 sequence quality good ({n_percent:.1f}% N's)")

            # Create new record with proper header
            new_header = create_reference_header(vp5_seq, sample_id, serotype, "VP5")
            new_record = vp5_seq
            new_record.id = new_header
            new_record.description = ""
            vp5_sequences.append(new_record)
        else:
            print(f"    ❌ No VP5 sequence found")

        print()

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Write VP2 reference database
    if vp2_sequences:
        vp2_output = os.path.join(output_dir, "BTV_VP2_references_improved.fasta")
        SeqIO.write(vp2_sequences, vp2_output, "fasta")
        print(f"✅ Created VP2 database: {vp2_output}")
        print(f"   Contains {len(vp2_sequences)} sequences")
    else:
        print("❌ No VP2 sequences found - cannot create VP2 database")

    # Write VP5 reference database
    if vp5_sequences:
        vp5_output = os.path.join(output_dir, "BTV_VP5_references_improved.fasta")
        SeqIO.write(vp5_sequences, vp5_output, "fasta")
        print(f"✅ Created VP5 database: {vp5_output}")
        print(f"   Contains {len(vp5_sequences)} sequences")
    else:
        print("❌ No VP5 sequences found - cannot create VP5 database")

    # Summary
    print()
    print("=" * 60)
    print("DATABASE CREATION SUMMARY")
    print("=" * 60)
    print(f"VP2 sequences: {len(vp2_sequences)}")
    print(f"VP5 sequences: {len(vp5_sequences)}")

    if vp2_sequences and vp5_sequences:
        print("✅ Both databases created successfully!")
        print()
        print("Next steps:")
        print("1. Review the generated databases")
        print("2. Replace the existing reference files:")
        print(f"   cp {vp2_output} reference_sequences/BTV_VP2_references.fasta")
        print(f"   cp {vp5_output} reference_sequences/BTV_VP5_references.fasta")
        print("3. Test serotyping with known samples")
    else:
        print("❌ Database creation incomplete - check input files")

def main():
    parser = argparse.ArgumentParser(
        prog='BTV Reference Database Builder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Build improved BTV serotype reference databases from validated samples.

        Extracts VP2 (segment 2) and VP5 (segment 6) sequences from reference-guided
        FASTA files with known serotypes and creates properly formatted reference
        databases for improved serotype prediction.
        '''))

    parser.add_argument('-i', '--input-dir',
                       required=True,
                       help='Directory containing reference-guided FASTA files')

    parser.add_argument('-o', '--output-dir',
                       default='new_reference_databases',
                       help='Output directory for new reference databases')

    parser.add_argument('-m', '--mapping-file',
                       help='Tab-separated file with sample_id -> serotype mapping')

    args = parser.parse_args()

    # No built-in mapping: pass the lab's own sample_id -> serotype list with
    # -m (tab-separated, one "SAMPLE-ID<TAB>BTV-17" per line).
    serotype_mapping = {}

    # Load mapping from file if provided
    if args.mapping_file:
        print(f"Loading serotype mapping from: {args.mapping_file}")
        serotype_mapping = {}
        with open(args.mapping_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        serotype_mapping[parts[0]] = parts[1]

    if not serotype_mapping:
        print("Error: No serotype mapping provided", file=sys.stderr)
        sys.exit(1)

    # Validate input directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    # Build the databases
    build_reference_databases(args.input_dir, serotype_mapping, args.output_dir)

if __name__ == "__main__":
    main()