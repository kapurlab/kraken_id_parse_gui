#!/usr/bin/env python3

"""
Add a validated sample to the curated BTV serotyping reference panel.
=====================================================================
The BTV serotyper (bin/btv_serotyping.py) calls serotypes ONLY against the
curated panel in reference_sequences/ (BTV_VP2_references.fasta and
BTV_VP5_references.fasta). This script is the controlled, documented way to
grow that panel as new validated samples become available.

It is intentionally conservative:
  * Serotype is defined by VP2 (segment 2); VP2 is required, VP5 is optional.
  * The claimed serotype is cross-checked against the EXISTING panel (real
    blastn). If the new VP2 looks like a DIFFERENT serotype already in the
    panel (possible mislabel), or the sample is already present, the add is
    REFUSED unless --force.
  * Panel FASTAs are backed up before any change.
  * A manifest (reference_sequences/BTV_panel_manifest.tsv) documents exactly
    which serotypes/samples are in the panel, and the script prints the current
    coverage plus suggestions every run.

NOTE: GenBank is NOT used to confirm serotype — public BTV serotype labels are
unreliable. Only the lab's curated panel is trusted.

Usage:
  add_btv_reference.py --fasta SAMPLE_reference_guided.fasta \\
        --sample-id SAMPLE-ID --serotype BTV-13
  add_btv_reference.py --vp2 my_vp2.fasta [--vp5 my_vp5.fasta] \\
        --sample-id ID --serotype BTV-13 [--force]
"""

import os
import re
import sys
import shutil
import argparse
import subprocess
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Shared, canonical segment resolver (VP2->2, VP5->6, segment/seg N, ...)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orbivirus_specific import resolve_segment_number

# --- tunables -----------------------------------------------------------------
# Expected nucleotide length ranges (warn, not fatal, if outside)
EXPECTED_LEN = {'VP2': (2700, 3100), 'VP5': (1500, 1850)}
MAX_N_PERCENT_WARN = 5.0     # warn above this
MAX_N_PERCENT_FATAL = 50.0   # refuse above this (unless --force)
# Conflict detection: a hit to a DIFFERENT serotype at/above these is a conflict
CONFLICT_PIDENT = 90.0
CONFLICT_QCOVS = 70.0
NEAR_DUP_PIDENT = 99.5       # near-identical to an existing same-serotype entry

PROTEIN_SEGMENT = {'VP2': 2, 'VP5': 6}


class bcolors:
    GREEN = '\033[92m'; YELLOW = '\033[93m'; RED = '\033[91m'
    BLUE = '\033[94m'; BOLD = '\033[1m'; ENDC = '\033[0m'


def info(msg):  print(f"{bcolors.BLUE}{msg}{bcolors.ENDC}")
def good(msg):  print(f"{bcolors.GREEN}{msg}{bcolors.ENDC}")
def warn(msg):  print(f"{bcolors.YELLOW}WARNING: {msg}{bcolors.ENDC}")
def err(msg):   print(f"{bcolors.RED}ERROR: {msg}{bcolors.ENDC}", file=sys.stderr)


def normalize_serotype(value):
    """Accept 'BTV-13', 'btv13', '13' -> 'BTV-13'; return None if unparseable."""
    if value is None:
        return None
    m = re.search(r'(\d{1,2})', str(value))
    if not m:
        return None
    n = int(m.group(1))
    if not (1 <= n <= 29):
        return None
    return f"BTV-{n}"


def seq_stats(record):
    s = str(record.seq).upper()
    length = len(s)
    n_pct = (s.count('N') / length * 100) if length else 100.0
    return length, n_pct


def extract_segments(args):
    """Return {'VP2': record, 'VP5': record} from the provided inputs."""
    found = {}
    if args.fasta:
        for record in SeqIO.parse(args.fasta, 'fasta'):
            seg = resolve_segment_number(f"{record.id} {record.description}")
            if seg == 2 and 'VP2' not in found:
                found['VP2'] = record
            elif seg == 6 and 'VP5' not in found:
                found['VP5'] = record
    if args.vp2:
        found['VP2'] = next(SeqIO.parse(args.vp2, 'fasta'))
    if args.vp5:
        found['VP5'] = next(SeqIO.parse(args.vp5, 'fasta'))
    return found


def load_panel(refs_dir, protein):
    path = os.path.join(refs_dir, f"BTV_{protein}_references.fasta")
    records = list(SeqIO.parse(path, 'fasta')) if os.path.exists(path) else []
    return path, records


def parse_header_serotype(record):
    m = re.search(r'\|(BTV-\d+)\|', record.id)
    return m.group(1) if m else None


def parse_header_sample(record):
    return record.id.split('|')[0]


def run_blast_subject(query_record, subject_path, tmp_dir):
    """blastn one query record against a panel FASTA (subject). Returns sorted hits."""
    if not subject_path or not os.path.exists(subject_path) or os.path.getsize(subject_path) == 0:
        return []
    qpath = os.path.join(tmp_dir, 'query.fasta')
    SeqIO.write([query_record], qpath, 'fasta')
    outfmt = '6 sseqid pident qcovs bitscore'
    cmd = ['blastn', '-query', qpath, '-subject', subject_path,
           '-outfmt', outfmt, '-max_target_seqs', '10', '-task', 'megablast']
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except FileNotFoundError:
        warn("blastn not found on PATH - skipping panel cross-check (cannot verify serotype)")
        return None  # signal "could not check"
    except subprocess.TimeoutExpired:
        warn("blastn timed out - skipping panel cross-check")
        return None
    if res.returncode != 0:
        warn(f"blastn failed during cross-check: {res.stderr.strip()[:200]}")
        return None
    hits = []
    for line in res.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 4:
            continue
        sseqid = f[0]
        m = re.search(r'\|(BTV-\d+)\|', sseqid)
        hits.append({
            'sseqid': sseqid,
            'serotype': m.group(1) if m else 'Unknown',
            'sample': sseqid.split('|')[0],
            'pident': float(f[1]), 'qcovs': float(f[2]), 'bitscore': float(f[3]),
        })
    hits.sort(key=lambda h: h['bitscore'], reverse=True)
    return hits


def make_header(sample_id, serotype, protein):
    n = serotype.split('-')[1]
    return f"{sample_id}|{serotype}|Bluetongue virus {n} sample {sample_id} {protein} gene, complete cds"


def coverage_summary(refs_dir):
    """Return {serotype: {'VP2': count, 'VP5': count}} from the current panel."""
    cov = {}
    for protein in ('VP2', 'VP5'):
        _, records = load_panel(refs_dir, protein)
        for r in records:
            st = parse_header_serotype(r)
            if not st:
                continue
            cov.setdefault(st, {'VP2': 0, 'VP5': 0})[protein] += 1
    return cov


def write_manifest(refs_dir):
    """(Re)build a human-readable manifest from the current panel contents."""
    manifest = os.path.join(refs_dir, 'BTV_panel_manifest.tsv')
    rows = []
    for protein in ('VP2', 'VP5'):
        _, records = load_panel(refs_dir, protein)
        for r in records:
            length, n_pct = seq_stats(r)
            rows.append({
                'serotype': parse_header_serotype(r) or 'Unknown',
                'sample_id': parse_header_sample(r),
                'segment': protein,
                'length': length,
                'n_percent': f"{n_pct:.1f}",
            })
    rows.sort(key=lambda x: (int(x['serotype'].split('-')[1]) if x['serotype'].startswith('BTV-') else 99,
                             x['sample_id'], x['segment']))
    stamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(manifest, 'w') as fh:
        fh.write(f"# BTV serotyping curated reference panel manifest\n")
        fh.write(f"# Regenerated: {stamp}\n")
        cov = coverage_summary(refs_dir)
        covered = ', '.join(f"{s} (VP2 n={c['VP2']}, VP5 n={c['VP5']})"
                            for s, c in sorted(cov.items(), key=lambda kv: int(kv[0].split('-')[1])))
        fh.write(f"# Serotypes callable: {covered or 'none'}\n")
        fh.write("serotype\tsample_id\tsegment\tlength\tn_percent\n")
        for r in rows:
            fh.write(f"{r['serotype']}\t{r['sample_id']}\t{r['segment']}\t{r['length']}\t{r['n_percent']}\n")
    return manifest


def append_record(refs_dir, protein, sample_id, serotype, record):
    path = os.path.join(refs_dir, f"BTV_{protein}_references.fasta")
    new = SeqRecord(record.seq, id=make_header(sample_id, serotype, protein), description="")
    with open(path, 'a') as fh:
        SeqIO.write([new], fh, 'fasta')


def backup_panel(refs_dir):
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backups = []
    for protein in ('VP2', 'VP5'):
        path = os.path.join(refs_dir, f"BTV_{protein}_references.fasta")
        if os.path.exists(path):
            bak = f"{path}.{stamp}.bak"
            shutil.copy2(path, bak)
            backups.append(bak)
    return backups


def main():
    parser = argparse.ArgumentParser(
        prog='add_btv_reference.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    parser.add_argument('--fasta', help='Reference-guided consensus FASTA (VP2/VP5 auto-extracted)')
    parser.add_argument('--vp2', help='Explicit VP2 (segment 2) FASTA (overrides --fasta VP2)')
    parser.add_argument('--vp5', help='Explicit VP5 (segment 6) FASTA (overrides --fasta VP5)')
    parser.add_argument('--sample-id', required=True, help='Sample identifier (e.g. SAMPLE-ID)')
    parser.add_argument('--serotype', required=True, help='Validated serotype, e.g. BTV-13')
    parser.add_argument('--force', action='store_true',
                        help='Proceed despite conflicts/duplicates/quality warnings')
    repo_refs = str(Path(__file__).resolve().parent.parent / 'reference_sequences')
    parser.add_argument('--refs-dir', default=repo_refs,
                        help=f'Curated panel directory (default: {repo_refs})')
    args = parser.parse_args()

    print(f"\n{bcolors.BOLD}=== Add BTV reference to curated panel ==={bcolors.ENDC}")

    # 1. Validate serotype
    serotype = normalize_serotype(args.serotype)
    if not serotype:
        err(f"Invalid --serotype '{args.serotype}'. Use BTV-N (1-29).")
        sys.exit(2)
    if not (args.fasta or args.vp2):
        err("Provide --fasta (consensus) or --vp2 (explicit VP2 FASTA).")
        sys.exit(2)
    if not os.path.isdir(args.refs_dir):
        err(f"Reference panel directory not found: {args.refs_dir}")
        sys.exit(2)

    info(f"Sample: {args.sample_id}   Claimed serotype: {serotype}")
    info(f"Panel:  {args.refs_dir}")

    # 2. Extract segments
    segs = extract_segments(args)
    if 'VP2' not in segs:
        err("VP2 (segment 2) not found in the input. VP2 is required to add a "
            "serotype reference (serotype is defined by VP2).")
        if args.fasta:
            info("Tip: check the consensus headers contain 'VP2'/'segment 2'.")
        sys.exit(1)
    if 'VP5' not in segs:
        warn("VP5 (segment 6) not found - adding VP2 only. VP2-based calls still "
             "work; consider adding VP5 later for confirmation.")

    # 3. Quality checks
    fatal_quality = False
    for protein, rec in segs.items():
        length, n_pct = seq_stats(rec)
        lo, hi = EXPECTED_LEN[protein]
        flag = ''
        if not (lo <= length <= hi):
            flag += f" [length {length} outside expected {lo}-{hi}]"
        if n_pct > MAX_N_PERCENT_WARN:
            flag += f" [{n_pct:.1f}% N]"
        line = f"  {protein}: length={length}, N={n_pct:.1f}%{flag}"
        if n_pct >= MAX_N_PERCENT_FATAL:
            err(line + "  <-- too many ambiguous bases")
            fatal_quality = True
        elif flag:
            warn(line)
        else:
            good(line)
    if fatal_quality and not args.force:
        err("Refusing due to sequence quality. Re-sequence, or override with --force.")
        sys.exit(1)

    # 4. Duplicate + conflict checks against the EXISTING panel (VP2)
    vp2_path, vp2_records = load_panel(args.refs_dir, 'VP2')
    existing_samples = {parse_header_sample(r) for r in vp2_records}
    covered_before = set(coverage_summary(args.refs_dir).keys())

    conflicts = []
    if args.sample_id in existing_samples:
        conflicts.append(f"sample '{args.sample_id}' is already in the panel (duplicate)")

    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        hits = run_blast_subject(segs['VP2'], vp2_path, tmp)
    if hits is None:
        warn("Could not run the panel cross-check (blastn unavailable). "
             "Serotype label will NOT be verified against the panel.")
    elif hits:
        top = hits[0]
        info(f"Closest existing panel VP2: {top['serotype']} ({top['sample']}) "
             f"{top['pident']:.1f}% id, {top['qcovs']:.0f}% cov")
        if (top['serotype'] not in ('Unknown', serotype)
                and top['pident'] >= CONFLICT_PIDENT and top['qcovs'] >= CONFLICT_QCOVS):
            conflicts.append(
                f"new VP2 best-matches an EXISTING different serotype {top['serotype']} "
                f"({top['sample']}, {top['pident']:.1f}% id) but you labeled it {serotype} "
                f"- possible mislabel")
        same = [h for h in hits if h['serotype'] == serotype]
        if same and same[0]['pident'] >= NEAR_DUP_PIDENT:
            warn(f"near-duplicate of existing {serotype} entry {same[0]['sample']} "
                 f"({same[0]['pident']:.1f}% id) - adding is OK but adds little new information.")
    else:
        info("No comparable hit in the existing panel (novel serotype is expected here).")

    if conflicts:
        err("Conflicts detected:")
        for c in conflicts:
            print(f"  - {c}", file=sys.stderr)
        if not args.force:
            err("Refusing to modify the panel. Re-check the label/sample, or override "
                "with --force if you are certain.")
            sys.exit(1)
        warn("Proceeding despite conflicts because --force was given.")

    # 5. Backup, then append
    backups = backup_panel(args.refs_dir)
    if backups:
        info("Backed up existing panel files:")
        for b in backups:
            print(f"  {os.path.basename(b)}")

    added = []
    for protein in ('VP2', 'VP5'):
        if protein in segs:
            append_record(args.refs_dir, protein, args.sample_id, serotype, segs[protein])
            added.append(protein)
    good(f"Appended {', '.join(added)} for {args.sample_id} ({serotype}).")

    # 6. Rebuild manifest + report coverage
    manifest = write_manifest(args.refs_dir)
    cov = coverage_summary(args.refs_dir)
    print(f"\n{bcolors.BOLD}Panel now covers:{bcolors.ENDC}")
    for st in sorted(cov, key=lambda s: int(s.split('-')[1])):
        c = cov[st]
        print(f"  {st}: VP2 n={c['VP2']}, VP5 n={c['VP5']}")
    info(f"Manifest updated: {manifest}")

    # 7. Suggestions
    print(f"\n{bcolors.BOLD}Notes:{bcolors.ENDC}")
    if serotype not in covered_before:
        good(f"  + {serotype} is NEW - the serotyper can now call {serotype}.")
    if 'VP5' not in segs:
        print("  - VP5 not added; add it later with --vp5 for VP2/VP5 confirmation.")
    print("  - Review the .bak files if you need to roll back this change.")
    print("  - Re-run a known sample of this serotype to confirm it now calls correctly.")


if __name__ == '__main__':
    main()
