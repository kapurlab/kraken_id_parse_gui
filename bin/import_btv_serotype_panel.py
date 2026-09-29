#!/usr/bin/env python3

"""
BTV Serotype Panel Importer
===========================
Merge a serotype-labeled VP2 (segment 2) reference FASTA into the curated BTV
serotyping panel used by ``bin/btv_serotyping.py``.

The lab's ``BTV_segment2_serotypes.txt`` carries one or more VP2 references per
serotype with headers like ``>BTV_segment2_serotype17``. This script:

  * skips malformed headers with no serotype number (e.g. ``>BTV_segment2_serotype``),
  * de-duplicates by sequence (the source file repeats serotype blocks),
  * drops length outliers that are not a single VP2 segment (a sanity guard that
    removes e.g. two concatenated VP2s),
  * rewrites each kept record with the panel header convention
    ``seg2ref_BTV-N_k|BTV-N|Bluetongue virus N segment 2 (VP2) serotype reference``
    (so ``panel_serotypes()`` / ``_serotype_from()`` in btv_serotyping.py pick it up),
  * backs up and APPENDS the new entries to ``reference_sequences/BTV_VP2_references.fasta``
    (existing lab consensuses are kept), and
  * regenerates ``reference_sequences/BTV_panel_manifest.tsv`` from the panel files.

Idempotent: sequences already present in the panel are not re-added, so the
script can be re-run safely.

Usage:
  python bin/import_btv_serotype_panel.py --serotype-fasta /path/to/BTV_segment2_serotypes.txt
  # options: --vp2-panel, --vp5-panel, --min-len, --max-len, --dry-run, --keep-raw-copy
"""

import argparse
import os
import re
import shutil
import sys
import textwrap
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# A single BTV VP2 (segment 2) CDS is ~2.7-3.1 kb. Anything well outside this is
# not a clean single-segment reference (e.g. two VP2s concatenated) -> skip it.
DEFAULT_MIN_LEN = 2400
DEFAULT_MAX_LEN = 3300

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_VP2 = REPO_ROOT / 'reference_sequences' / 'BTV_VP2_references.fasta'
DEFAULT_VP5 = REPO_ROOT / 'reference_sequences' / 'BTV_VP5_references.fasta'
DEFAULT_MANIFEST = REPO_ROOT / 'reference_sequences' / 'BTV_panel_manifest.tsv'

SERO_RE = re.compile(r'serotype\s*(\d{1,2})', re.IGNORECASE)
PANEL_SERO_RE = re.compile(r'\|(BTV-\d+)\|')


def norm_seq(seq):
    """Normalized (uppercase, whitespace-stripped) sequence for de-duplication."""
    return re.sub(r'\s+', '', str(seq)).upper()


def serotype_of(header):
    """Return 'BTV-N' parsed from a source header, or None if absent."""
    m = SERO_RE.search(header)
    return f"BTV-{int(m.group(1))}" if m else None


def existing_panel_seqs(panel_path):
    """Set of normalized sequences already in the panel (for idempotent merge)."""
    seqs = set()
    if os.path.exists(panel_path):
        for rec in SeqIO.parse(panel_path, 'fasta'):
            seqs.add(norm_seq(rec.seq))
    return seqs


def panel_counts(panel_path):
    """serotype -> count of records in a panel FASTA (by |BTV-N| header tag)."""
    counts = defaultdict(int)
    if os.path.exists(panel_path):
        with open(panel_path) as fh:
            for line in fh:
                if line.startswith('>'):
                    m = PANEL_SERO_RE.search(line)
                    if m:
                        counts[m.group(1)] += 1
    return counts


def load_new_records(serotype_fasta, min_len, max_len):
    """Parse + clean the serotype-labeled source FASTA.

    Returns (kept, skipped) where kept is a list of (serotype, seq_str) with
    unique sequences per serotype, and skipped is a list of (header, reason).
    """
    kept_keys = set()          # (serotype, normseq) already taken
    kept = []
    skipped = []
    for rec in SeqIO.parse(serotype_fasta, 'fasta'):
        header = f"{rec.id} {rec.description}".strip()
        sero = serotype_of(header)
        seq = norm_seq(rec.seq)
        if not seq:
            skipped.append((header, 'empty sequence'))
            continue
        if sero is None:
            skipped.append((header, 'no serotype number in header'))
            continue
        if not (min_len <= len(seq) <= max_len):
            skipped.append((header, f'length {len(seq)} outside VP2 range '
                                    f'{min_len}-{max_len} bp'))
            continue
        key = (sero, seq)
        if key in kept_keys:
            skipped.append((header, 'duplicate sequence'))
            continue
        kept_keys.add(key)
        kept.append((sero, str(rec.seq).replace(' ', '')))
    return kept, skipped


def regenerate_manifest(vp2_path, vp5_path, manifest_path):
    """Rewrite the manifest from the current VP2/VP5 panel files."""
    rows = []
    callable_counts = defaultdict(lambda: {'VP2': 0, 'VP5': 0})
    for segment, path in (('VP2', vp2_path), ('VP5', vp5_path)):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, 'fasta'):
            header = rec.id  # ID|BTV-N|... -> ID is the first pipe field
            parts = header.split('|')
            sample_id = parts[0]
            m = PANEL_SERO_RE.search(f"|{'|'.join(parts[1:])}|")
            sero = m.group(1) if m else (parts[1] if len(parts) > 1 else 'BTV-?')
            seq = str(rec.seq)
            length = len(seq)
            n_pct = (seq.upper().count('N') / length * 100) if length else 0.0
            rows.append((sero, sample_id, segment, length, n_pct))
            if sero in callable_counts or True:
                callable_counts[sero][segment] += 1

    def sero_key(s):
        m = re.search(r'(\d+)', s)
        return int(m.group(1)) if m else 999
    rows.sort(key=lambda r: (sero_key(r[0]), r[2], r[1]))

    callable_summary = ', '.join(
        f"{s} (VP2 n={callable_counts[s]['VP2']}, VP5 n={callable_counts[s]['VP5']})"
        for s in sorted(callable_counts, key=sero_key)
    )

    with open(manifest_path, 'w') as fh:
        fh.write("# BTV serotyping curated reference panel manifest\n")
        fh.write(f"# Regenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"# Serotypes callable: {callable_summary}\n")
        fh.write("serotype\tsample_id\tsegment\tlength\tn_percent\n")
        for sero, sample_id, segment, length, n_pct in rows:
            fh.write(f"{sero}\t{sample_id}\t{segment}\t{length}\t{n_pct:.1f}\n")
    return callable_summary


def main():
    parser = argparse.ArgumentParser(
        prog='import_btv_serotype_panel.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__))
    parser.add_argument('-s', '--serotype-fasta', required=True,
                        help='Serotype-labeled VP2 FASTA (e.g. BTV_segment2_serotypes.txt)')
    parser.add_argument('--vp2-panel', default=str(DEFAULT_VP2),
                        help='VP2 panel FASTA to append to (default: reference_sequences/BTV_VP2_references.fasta)')
    parser.add_argument('--vp5-panel', default=str(DEFAULT_VP5),
                        help='VP5 panel FASTA (read-only; used for manifest)')
    parser.add_argument('--manifest', default=str(DEFAULT_MANIFEST),
                        help='Manifest TSV to regenerate')
    parser.add_argument('--min-len', type=int, default=DEFAULT_MIN_LEN)
    parser.add_argument('--max-len', type=int, default=DEFAULT_MAX_LEN)
    parser.add_argument('--keep-raw-copy', action='store_true',
                        help='Also copy the raw serotype FASTA next to the panel for provenance')
    parser.add_argument('--dry-run', action='store_true',
                        help='Report what would change without writing')
    args = parser.parse_args()

    if not os.path.exists(args.serotype_fasta):
        print(f"Error: serotype FASTA not found: {args.serotype_fasta}", file=sys.stderr)
        sys.exit(1)

    print("=" * 64)
    print("BTV SEROTYPE PANEL IMPORTER")
    print("=" * 64)
    print(f"Source serotype FASTA: {args.serotype_fasta}")
    print(f"VP2 panel:             {args.vp2_panel}")
    print(f"VP2 length window:     {args.min_len}-{args.max_len} bp\n")

    kept, skipped = load_new_records(args.serotype_fasta, args.min_len, args.max_len)

    if skipped:
        print(f"Skipped {len(skipped)} record(s):")
        for header, reason in skipped:
            print(f"  - [{reason}] {header[:70]}")
        print()

    # Drop sequences already in the panel (idempotent re-runs).
    already = existing_panel_seqs(args.vp2_panel)
    to_add = [(sero, seq) for (sero, seq) in kept if norm_seq(seq) not in already]
    n_dup_existing = len(kept) - len(to_add)
    if n_dup_existing:
        print(f"{n_dup_existing} cleaned record(s) already present in the panel - skipping those.\n")

    # Build new SeqRecords with panel-convention headers.
    per_sero_idx = defaultdict(int)
    new_records = []
    for sero, seq in sorted(to_add, key=lambda x: (int(x[0].split('-')[1]),)):
        per_sero_idx[sero] += 1
        k = per_sero_idx[sero]
        n = sero.split('-')[1]
        seq_id = f"seg2ref_{sero}_{k}|{sero}|Bluetongue virus {n} segment 2 (VP2) serotype reference"
        rec = SeqRecord(Seq(seq), id=seq_id, description="")
        new_records.append(rec)

    added_by_sero = defaultdict(int)
    for sero, _ in to_add:
        added_by_sero[sero] += 1
    print(f"New VP2 references to add: {len(new_records)} across "
          f"{len(added_by_sero)} serotype(s)")
    for sero in sorted(added_by_sero, key=lambda s: int(s.split('-')[1])):
        print(f"  {sero}: +{added_by_sero[sero]}")
    print()

    if args.dry_run:
        print("[dry-run] no files written.")
        return 0

    # Back up and append.
    if new_records:
        if os.path.exists(args.vp2_panel):
            backup = args.vp2_panel + '.bak'
            shutil.copy2(args.vp2_panel, backup)
            print(f"Backed up existing VP2 panel -> {backup}")
        with open(args.vp2_panel, 'a') as fh:
            SeqIO.write(new_records, fh, 'fasta')
        print(f"Appended {len(new_records)} record(s) to {args.vp2_panel}")
    else:
        print("Nothing new to append (panel already covers these sequences).")

    if args.keep_raw_copy:
        dest = os.path.join(os.path.dirname(args.vp2_panel),
                            os.path.basename(args.serotype_fasta))
        if os.path.abspath(dest) != os.path.abspath(args.serotype_fasta):
            shutil.copy2(args.serotype_fasta, dest)
            print(f"Copied raw serotype FASTA -> {dest} (provenance)")

    summary = regenerate_manifest(args.vp2_panel, args.vp5_panel, args.manifest)
    print(f"\nRegenerated manifest: {args.manifest}")
    print(f"Serotypes callable: {summary}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
