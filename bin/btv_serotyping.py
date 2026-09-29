#!/usr/bin/env python3

"""
BTV Serotyping Module
=====================
Determines BTV serotype by BLASTing VP2 (segment 2) and VP5 (segment 6)
consensus sequences against curated reference databases.
"""

import os
import re
import csv
import subprocess
import logging
import sys
from pathlib import Path
from datetime import datetime

# Ensure the bin/ directory is importable when run standalone
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orbivirus_specific import resolve_segment_number, orbivirus_species

logger = logging.getLogger('kraken_pipeline')


class BTVSerotyping:
    """Run BTV serotyping analysis using BLAST against VP2/VP5 reference databases."""

    # Minimum identity/coverage for a confident curated-panel serotype call.
    # A best hit below either threshold yields "additional reference data needed"
    # rather than a guess. Tunable as the curated panel grows.
    PANEL_MIN_PIDENT = 90.0
    PANEL_MIN_QCOVS = 70.0

    def __init__(self, consensus_fasta, output_dir='btv_serotype', debug=False):
        self.consensus_fasta = consensus_fasta
        self.output_dir = output_dir
        self.debug = debug
        self.predictions = []
        self.all_hits = {}  # protein_type -> list of all hits for detailed results
        self.interpretation = ''
        self.tsv_file = None

        # Curated BTV serotype reference panel (the ONLY reference used for calls).
        # These files live in the repo and grow via bin/add_btv_reference.py.
        script_dir = Path(__file__).resolve().parent.parent
        ref_dir = script_dir / 'reference_sequences'
        self.vp2_db = str(ref_dir / 'BTV_VP2_references.fasta')
        self.vp5_db = str(ref_dir / 'BTV_VP5_references.fasta')

    def panel_serotypes(self):
        """Return the sorted list of serotypes the curated panel can currently
        call (parsed from VP2/VP5 reference headers, e.g. 'ACC|BTV-17|...')."""
        seen = set()
        for db in (self.vp2_db, self.vp5_db):
            if not db or not os.path.exists(db):
                continue
            try:
                with open(db) as fh:
                    for line in fh:
                        if line.startswith('>'):
                            m = re.search(r'\|(BTV-\d+)\|', line)
                            if m:
                                seen.add(m.group(1))
            except Exception:
                continue
        return sorted(seen, key=lambda s: int(s.split('-')[1]))

    def _extract_segment_fastas(self):
        """Extract VP2 (segment 2) from the consensus FASTA. Serotype is called
        from segment 2 only; other segments (incl. VP5/segment 6) are ignored."""
        from Bio import SeqIO

        segment_files = {}
        logger.debug(f"Parsing consensus FASTA: {self.consensus_fasta}")

        for record in SeqIO.parse(self.consensus_fasta, 'fasta'):
            desc = f"{record.id} {record.description}"
            logger.debug(f"Processing sequence: {record.id}")
            logger.debug(f"Full description: {desc}")

            # Use the canonical resolver (shared with the coverage-graph segment
            # status) so VP2->segment 2 detection is consistent.
            seg_num = resolve_segment_number(desc)
            if seg_num is not None:
                logger.debug(f"Resolved segment number {seg_num} for {record.id}")

            if seg_num == 2 and orbivirus_species(record.description) == 'EHD':
                # A mixed BTV + EHD sample has a segment 2 for each virus; the EHD
                # one must not be serotyped against the BTV panel (it used to
                # overwrite the BTV VP2 query and yield "Undetermined").
                logger.info(f"Skipping EHD segment 2 (not a BTV serotyping target): {record.id}")
            elif seg_num == 2:  # VP2 -- the only segment used for serotyping
                outpath = os.path.join(self.output_dir, 'vp2_query.fasta')
                SeqIO.write([record], outpath, 'fasta')
                segment_files['VP2'] = outpath
                logger.info(f"Found VP2 (segment 2): {record.id}")
            elif seg_num is not None:
                logger.debug(f"Segment {seg_num} not used for serotyping (segment 2 only)")
            else:
                logger.debug(f"Could not detect segment number for: {record.id}")

        logger.info(f"Extracted {len(segment_files)} segment file(s) for serotyping")
        return segment_files

    @staticmethod
    def _serotype_from(sseqid, stitle):
        """Parse the BTV serotype from a hit's subject id and/or title.

        Handles both the curated panel header format (``ACC|BTV-17|...``) and
        free-text GenBank titles (``Bluetongue virus 17 ... VP2``)."""
        for text in (sseqid or '', stitle or ''):
            m = re.search(r'BTV[-\s]?(\d{1,2})\b', text, re.IGNORECASE)
            if m:
                return f"BTV-{int(m.group(1))}"
        # GenBank titles: "Bluetongue virus 17 ..." / "...serotype 17..."
        m = re.search(r'[Bb]luetongue\s+virus\s+(?:serotype\s+)?(\d{1,2})\b', stitle or '')
        if m:
            return f"BTV-{int(m.group(1))}"
        m = re.search(r'serotype\s+(\d{1,2})\b', stitle or '', re.IGNORECASE)
        if m:
            return f"BTV-{int(m.group(1))}"
        return 'Unknown'

    def _run_blastn(self, query_fasta, protein_type, subject, max_hits=5):
        """Run real blastn of the consensus query against the curated panel
        FASTA (via ``-subject`` — no makeblastdb needed). Returns a list of hit
        dicts sorted by bitscore (best first)."""
        import subprocess
        if not subject or not os.path.exists(subject):
            logger.error(f"Curated panel reference not found: {subject}")
            return []
        outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle'
        cmd = ['blastn', '-query', query_fasta, '-subject', subject, '-outfmt', outfmt,
               '-max_target_seqs', str(max_hits), '-task', 'megablast']

        logger.info(f"Running blastn for {protein_type} against curated panel...")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        except FileNotFoundError:
            logger.error("blastn not found on PATH - cannot run BTV serotyping BLAST")
            return []
        except subprocess.TimeoutExpired:
            logger.error(f"blastn timed out for {protein_type} against curated panel")
            return []
        if result.returncode != 0:
            logger.error(f"blastn failed for {protein_type} (panel): {result.stderr.strip()[:300]}")
            return []

        hits = []
        for line in result.stdout.splitlines():
            f = line.split('\t')
            if len(f) < 14:
                continue
            sseqid, stitle = f[1], f[13]
            try:
                hit = {
                    'query': f[0], 'subject': sseqid,
                    'pident': float(f[2]), 'length': int(f[3]),
                    'mismatch': int(f[4]), 'gapopen': int(f[5]),
                    'qstart': int(f[6]), 'qend': int(f[7]),
                    'sstart': int(f[8]), 'send': int(f[9]),
                    'evalue': float(f[10]), 'bitscore': float(f[11]),
                    'qcovs': float(f[12]), 'stitle': stitle,
                    'accession': sseqid.split('|')[0],
                    'serotype': self._serotype_from(sseqid, stitle),
                    'source': 'panel',
                }
            except ValueError:
                continue
            hits.append(hit)

        hits.sort(key=lambda x: x['bitscore'], reverse=True)

        # Debug dump
        out_file = os.path.join(self.output_dir, f'{protein_type}_panel_results.txt')
        try:
            with open(out_file, 'w') as fh:
                fh.write(f"# blastn results for {protein_type} vs curated panel\n")
                fh.write("query\tsubject\tpident\tlength\tevalue\tbitscore\tqcovs\tserotype\tstitle\n")
                for h in hits[:10]:
                    fh.write(f"{h['query']}\t{h['subject']}\t{h['pident']:.1f}\t{h['length']}\t"
                             f"{h['evalue']:.2e}\t{h['bitscore']:.1f}\t{h['qcovs']:.0f}\t{h['serotype']}\t{h['stitle']}\n")
        except Exception:
            pass

        logger.info(f"  {protein_type}/panel: {len(hits)} hits"
                    + (f"; top={hits[0]['serotype']} {hits[0]['pident']:.1f}% id, {hits[0]['qcovs']:.0f}% cov" if hits else ""))
        return hits[:max_hits]

    def _determine_serotype(self, hits, protein_type):
        """Determine best serotype call from BLAST hits."""
        if not hits:
            logger.warning(f"No BLAST hits found for {protein_type}")
            return {
                'Protein': protein_type,
                'Serotype': 'No hit',
                'Top Hit Accession': '-',
                'Percent Identity': '-',
                'Query Coverage': '-',
                'Bitscore': '-',
                'E-value': '-',
            }

        hits.sort(key=lambda x: x['bitscore'], reverse=True)
        best = hits[0]

        # Log detailed information about top hits for debugging
        logger.info(f"Top 3 {protein_type} hits:")
        for i, hit in enumerate(hits[:3]):
            logger.info(f"  {i+1}. {hit['serotype']} ({hit['accession']}): "
                       f"{hit['pident']:.1f}% ID, {hit['qcovs']:.0f}% Cov, "
                       f"Score: {hit['bitscore']:.1f}")

        # Check for potential issues
        if best['pident'] < 80:
            logger.warning(f"Low identity for {protein_type}: {best['pident']:.1f}%")
        if best['qcovs'] < 70:
            logger.warning(f"Low coverage for {protein_type}: {best['qcovs']:.0f}%")

        # Check for close secondary hits of a DIFFERENT serotype (ambiguity)
        if len(hits) > 1:
            second_best = hits[1]
            score_diff = best['bitscore'] - second_best['bitscore']
            if score_diff < 10 and second_best['serotype'] != best['serotype']:
                logger.warning(f"Close secondary hit for {protein_type}: "
                             f"{second_best['serotype']} (score diff: {score_diff:.1f})")

        # Confidence gate: only assign a serotype when the best curated-panel hit
        # is strong enough; otherwise report that more reference data is needed.
        confident = best['pident'] >= self.PANEL_MIN_PIDENT and best['qcovs'] >= self.PANEL_MIN_QCOVS
        serotype = best['serotype'] if confident else 'Undetermined'
        if not confident:
            logger.warning(
                f"{protein_type} best panel match {best['serotype']} "
                f"({best['pident']:.1f}% id, {best['qcovs']:.0f}% cov) is below the "
                f"confidence threshold ({self.PANEL_MIN_PIDENT:.0f}% id / {self.PANEL_MIN_QCOVS:.0f}% cov) "
                f"- additional reference data needed")

        return {
            'Protein': protein_type,
            'Serotype': serotype,
            'Top Hit Accession': best['accession'],
            'Percent Identity': f"{best['pident']:.1f}%",
            'Query Coverage': f"{best['qcovs']:.0f}%",
            'Bitscore': f"{best['bitscore']:.1f}",
            'E-value': f"{best['evalue']:.2e}",
            'Source': best.get('source', 'panel') or 'panel',
            'Best Panel Match': best['serotype'],  # closest serotype even if sub-threshold
            'Confident': confident,
        }

    # Serotype values that do not represent a confident serotype call
    _INVALID_SEROTYPES = ('No hit', 'Unknown', 'Segments not found', 'Undetermined', '-', '')

    def _serotype_for_protein(self, protein_type):
        """Return the confident serotype call for a given protein (VP2/VP5), or None."""
        for pred in self.predictions:
            if pred.get('Protein') == protein_type:
                sero = pred.get('Serotype')
                if sero not in self._INVALID_SEROTYPES:
                    return sero
        return None

    def _prediction_for_protein(self, protein_type):
        """Return the raw prediction dict for a protein, or None."""
        for pred in self.predictions:
            if pred.get('Protein') == protein_type:
                return pred
        return None

    def _coverage_note(self):
        covered = self.panel_serotypes()
        return ('current reference panel covers: ' + ', '.join(covered)) if covered \
            else 'reference panel is empty'

    def _determine_interpretation(self):
        """Determine interpretation. Serotype is determined SOLELY by VP2
        (segment 2) matched against the curated panel. VP5 (segment 6), if
        recovered, is reported only as supporting context and never affects the
        serotype call. Calls are made only against the curated panel."""
        vp2 = self._serotype_for_protein('VP2')
        vp2_pred = self._prediction_for_protein('VP2')

        if vp2 is None:
            # run() records a 'Segments not found' placeholder when no segment 2
            # was recovered, so treat that the same as no VP2 prediction at all.
            if vp2_pred is None or vp2_pred.get('Serotype') == 'Segments not found':
                # VP2 segment was not recovered from the consensus at all
                self.interpretation = ('Serotype cannot be determined: Segment 2 (VP2) '
                                       'is required and was not found')
            else:
                # VP2 present but no confident panel match -> need more reference data
                best = vp2_pred.get('Best Panel Match', '-')
                pid = vp2_pred.get('Percent Identity', '-')
                cov = vp2_pred.get('Query Coverage', '-')
                self.interpretation = (
                    'Serotype undetermined - additional reference data needed. '
                    f'Best curated-panel match for VP2 was {best} ({pid} identity, {cov} coverage), '
                    f'below the confidence threshold '
                    f'({self.PANEL_MIN_PIDENT:.0f}% identity / {self.PANEL_MIN_QCOVS:.0f}% coverage); '
                    f'{self._coverage_note()}.')
        else:
            self.interpretation = f'Serotype assigned from Segment 2 (VP2) only: {vp2}'

    def _write_detailed_results(self):
        """Write btv_serotyping_detailed_results.tsv with all BLAST hits."""
        detailed_file = os.path.join(self.output_dir, 'btv_serotyping_detailed_results.tsv')
        fieldnames = ['Protein', 'Rank', 'Serotype', 'Accession', 'Percent Identity',
                     'Alignment Length', 'Mismatches', 'Gap Opens',
                     'Query Start', 'Query End', 'Subject Start', 'Subject End',
                     'E-value', 'Bitscore', 'Query Coverage']
        with open(detailed_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for protein_type, hits in self.all_hits.items():
                for rank, hit in enumerate(sorted(hits, key=lambda x: x['bitscore'], reverse=True), 1):
                    writer.writerow({
                        'Protein': protein_type,
                        'Rank': rank,
                        'Serotype': hit['serotype'],
                        'Accession': hit['accession'],
                        'Percent Identity': f"{hit['pident']:.1f}",
                        'Alignment Length': hit['length'],
                        'Mismatches': hit['mismatch'],
                        'Gap Opens': hit['gapopen'],
                        'Query Start': hit['qstart'],
                        'Query End': hit['qend'],
                        'Subject Start': hit['sstart'],
                        'Subject End': hit['send'],
                        'E-value': f"{hit['evalue']:.2e}",
                        'Bitscore': f"{hit['bitscore']:.1f}",
                        'Query Coverage': f"{hit['qcovs']:.0f}",
                    })

    def _write_summary(self):
        """Write btv_serotyping_summary.txt human-readable summary."""
        summary_file = os.path.join(self.output_dir, 'btv_serotyping_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("BTV Serotyping Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input: {self.consensus_fasta}\n\n")

            consensus = self.get_consensus_serotype()
            f.write(f"Consensus Serotype: {consensus}\n")
            f.write(f"Interpretation: {self.interpretation}\n\n")

            f.write("Per-Protein Results:\n")
            f.write("-" * 50 + "\n")
            for pred in self.predictions:
                f.write(f"  {pred['Protein']}:\n")
                f.write(f"    Serotype:       {pred['Serotype']}\n")
                f.write(f"    Top Hit:        {pred['Top Hit Accession']}\n")
                f.write(f"    Identity:       {pred['Percent Identity']}\n")
                f.write(f"    Query Coverage: {pred['Query Coverage']}\n")
                f.write(f"    Bitscore:       {pred['Bitscore']}\n")
                f.write(f"    E-value:        {pred['E-value']}\n\n")

            f.write("Reference Database (segment 2 / VP2 only):\n")
            f.write(f"  VP2: {self.vp2_db}\n")

    def run(self):
        """Execute serotyping analysis. Returns list of prediction dicts."""
        os.makedirs(self.output_dir, exist_ok=True)

        logger.info("Running BTV serotyping analysis...")

        segment_files = self._extract_segment_fastas()

        if not segment_files:
            logger.warning("No VP2 (segment 2) found in consensus FASTA")
            self.predictions = [{
                'Protein': 'VP2',
                'Serotype': 'Segments not found',
                'Top Hit Accession': '-',
                'Percent Identity': '-',
                'Query Coverage': '-',
                'Bitscore': '-',
                'E-value': '-',
            }]
        else:
            db_map = {'VP2': self.vp2_db}
            for protein_type, query_file in segment_files.items():
                panel_db = db_map[protein_type]
                hits = self._run_blastn(query_file, protein_type, subject=panel_db)
                self.all_hits[protein_type] = hits
                prediction = self._determine_serotype(hits, protein_type)
                self.predictions.append(prediction)
                logger.info(f"  {protein_type} serotype call: {prediction['Serotype']} "
                          f"({prediction['Percent Identity']} identity, "
                          f"{prediction['Query Coverage']} coverage)")

        # Determine interpretation
        self._determine_interpretation()

        # Write predictions TSV
        self.tsv_file = os.path.join(self.output_dir, 'btv_serotyping_predictions.tsv')
        if self.predictions:
            fieldnames = ['Protein', 'Serotype', 'Top Hit Accession', 'Percent Identity',
                         'Query Coverage', 'Bitscore', 'E-value', 'Source']
            with open(self.tsv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                        extrasaction='ignore')
                writer.writeheader()
                writer.writerows(self.predictions)

        # Write detailed results and summary
        self._write_detailed_results()
        self._write_summary()

        # Clean up intermediate files unless debug
        if not self.debug:
            for fname in ['vp2_query.fasta', 'vp5_query.fasta',
                         'VP2_blast_results.txt', 'VP5_blast_results.txt']:
                fpath = os.path.join(self.output_dir, fname)
                if os.path.exists(fpath):
                    os.remove(fpath)

        logger.info(f"Serotyping results written to: {self.output_dir}")
        return self.predictions

    def debug_consensus_analysis(self):
        """Debug method to analyze consensus FASTA and report findings."""
        from Bio import SeqIO

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        debug_file = os.path.join(self.output_dir, 'consensus_debug_analysis.txt')

        with open(debug_file, 'w') as f:
            f.write("CONSENSUS FASTA DEBUG ANALYSIS\n")
            f.write("=" * 70 + "\n")
            f.write(f"Input file: {self.consensus_fasta}\n")
            f.write(f"Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            if not os.path.exists(self.consensus_fasta):
                f.write(f"ERROR: Consensus FASTA file not found: {self.consensus_fasta}\n")
                return

            f.write("SEQUENCES FOUND:\n")
            f.write("-" * 70 + "\n")

            total_sequences = 0
            segments_found = {}

            for record in SeqIO.parse(self.consensus_fasta, 'fasta'):
                total_sequences += 1
                full_header = f"{record.id} {record.description}"
                f.write(f"\nSequence {total_sequences}:\n")
                f.write(f"  ID: {record.id}\n")
                f.write(f"  Description: {record.description}\n")
                f.write(f"  Full header: {full_header}\n")
                f.write(f"  Length: {len(record.seq)} bp\n")

                # Test all segment detection patterns
                desc_lower = full_header.lower()
                patterns_tested = [
                    (r'segment\s*(\d+)', "segment N"),
                    (r'_segment(\d+)', "_segmentN"),
                    (r'seg(\d+)', "segN"),
                    (r's(\d+)', "sN"),
                    (r'segment_(\d+)', "segment_N"),
                    (r'\.(\d+)\.', ".N."),
                ]

                segment_detected = False
                f.write("  Pattern matching results:\n")
                for pattern, description in patterns_tested:
                    match = re.search(pattern, desc_lower)
                    if match:
                        seg_num = int(match.group(1))
                        segments_found[seg_num] = record.id
                        f.write(f"    ✓ {description}: Found segment {seg_num}\n")
                        segment_detected = True

                        # Check if this is VP2 or VP5
                        if seg_num == 2:
                            f.write(f"      → This is VP2 (needed for serotyping)\n")
                        elif seg_num == 6:
                            f.write(f"      → This is VP5 (needed for serotyping)\n")
                        else:
                            f.write(f"      → Segment {seg_num} not used for serotyping\n")
                    else:
                        f.write(f"    ✗ {description}: No match\n")

                if not segment_detected:
                    f.write("    ⚠️  NO SEGMENT NUMBER DETECTED\n")
                    f.write("    → This sequence will be ignored for serotyping\n")

                    # Check for common issues
                    if 'guided by' in desc_lower:
                        f.write("    → Detected 'guided by' - this looks like a renamed consensus\n")
                    if any(virus in desc_lower for virus in ['bluetongue', 'btv', 'epizootic']):
                        f.write("    → Contains BTV/Orbivirus keywords - should have segment info\n")

            f.write(f"\nSUMMARY:\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total sequences: {total_sequences}\n")
            f.write(f"Segments detected: {len(segments_found)}\n")

            if segments_found:
                f.write("Segment mapping:\n")
                for seg_num in sorted(segments_found.keys()):
                    f.write(f"  Segment {seg_num}: {segments_found[seg_num]}\n")
            else:
                f.write("⚠️  NO SEGMENTS DETECTED!\n")

            # Check for required segments
            vp2_found = 2 in segments_found
            vp5_found = 6 in segments_found
            f.write(f"\nSerotyping requirements:\n")
            f.write(f"  VP2 (segment 2): {'✓ Found' if vp2_found else '✗ Missing'}\n")
            f.write(f"  VP5 (segment 6): {'✓ Found' if vp5_found else '✗ Missing'}\n")

            if not vp2_found and not vp5_found:
                f.write("\n❌ CRITICAL: No VP2 or VP5 segments found - serotyping will fail\n")
                f.write("RECOMMENDATIONS:\n")
                f.write("1. Check that consensus sequences have proper segment labels\n")
                f.write("2. Verify the reference-guided assembly worked correctly\n")
                f.write("3. Check if segments 2 and 6 have sufficient coverage\n")
            elif not vp2_found:
                f.write("\n⚠️  WARNING: VP2 (segment 2) missing - limited serotyping\n")
            elif not vp5_found:
                f.write("\n⚠️  WARNING: VP5 (segment 6) missing - limited serotyping\n")
            else:
                f.write("\n✓ Both VP2 and VP5 found - full serotyping possible\n")

        logger.info(f"Debug analysis written to: {debug_file}")
        return debug_file

    def get_consensus_serotype(self):
        """Return the predicted serotype. Serotype is defined by VP2 (segment 2)
        matched against the curated panel; without a confident VP2 panel call no
        serotype is reported (VP5 alone is never sufficient)."""
        vp2 = self._serotype_for_protein('VP2')
        if vp2 is None:
            pred = self._prediction_for_protein('VP2')
            if pred is None or pred.get('Serotype') == 'Segments not found':
                # Nothing to match: more panel references would not help here.
                return 'Undetermined - segment 2 (VP2) not found'
            return 'Undetermined - additional reference data needed'
        return vp2

    def excel(self, excel_dict):
        """Add serotyping results to the Excel summary dictionary."""
        consensus = self.get_consensus_serotype()
        excel_dict['BTV Serotype'] = consensus
        excel_dict['Serotype Interpretation'] = self.interpretation

        for pred in self.predictions:
            protein = pred['Protein']
            excel_dict[f'{protein} Serotype'] = pred['Serotype']
            excel_dict[f'{protein} Identity'] = pred['Percent Identity']
            excel_dict[f'{protein} Coverage'] = pred['Query Coverage']
            excel_dict[f'{protein} Top Hit'] = pred['Top Hit Accession']


def main():
    """Standalone BTV serotyping script."""
    import argparse
    import textwrap
    import sys

    parser = argparse.ArgumentParser(
        prog='BTV Serotyping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        ---------------------------------------------------------
        BTV Serotyping Analysis

        Determines BTV serotype by BLASTing VP2 (segment 2) and
        VP5 (segment 6) consensus sequences against curated
        reference databases.

        Input: Consensus FASTA file containing multiple segments
        Output: Serotype predictions and detailed results
        ---------------------------------------------------------
        '''))

    # Required arguments
    parser.add_argument('-i', '--input',
                       required=True,
                       help='Input consensus FASTA file')

    # Optional arguments
    parser.add_argument('-o', '--output',
                       default='btv_serotype',
                       help='Output directory (default: btv_serotype)')
    parser.add_argument('-d', '--debug',
                       action='store_true',
                       help='Enable debug mode and keep intermediate files')
    parser.add_argument('-v', '--version',
                       action='version',
                       version='%(prog)s 1.0')

    args = parser.parse_args()

    # Set up logging
    import logging
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Check reference databases
    script_dir = Path(__file__).resolve().parent.parent
    ref_dir = script_dir / 'reference_sequences'
    vp2_db = str(ref_dir / 'BTV_VP2_references.fasta')
    vp5_db = str(ref_dir / 'BTV_VP5_references.fasta')

    missing_dbs = []
    if not os.path.exists(vp2_db):
        missing_dbs.append(f"VP2 database: {vp2_db}")
    if not os.path.exists(vp5_db):
        missing_dbs.append(f"VP5 database: {vp5_db}")

    if missing_dbs:
        print("Error: Missing reference databases:", file=sys.stderr)
        for db in missing_dbs:
            print(f"  {db}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{os.path.basename(__file__)} Settings:")
    print(f"  Input FASTA: {args.input}")
    print(f"  Output directory: {args.output}")
    print(f"  Debug mode: {args.debug}")
    print(f"  VP2 database: {vp2_db}")
    print(f"  VP5 database: {vp5_db}\n")

    try:
        # Run serotyping analysis
        serotyper = BTVSerotyping(
            consensus_fasta=args.input,
            output_dir=args.output,
            debug=args.debug
        )

        # Run debug analysis if requested
        if args.debug:
            debug_file = serotyper.debug_consensus_analysis()
            print(f"Debug analysis written to: {debug_file}")

        predictions = serotyper.run()
        consensus_serotype = serotyper.get_consensus_serotype()

        # Print results to console
        print("=" * 60)
        print("BTV SEROTYPING RESULTS")
        print("=" * 60)
        print(f"Consensus Serotype: {consensus_serotype}")
        print(f"Interpretation: {serotyper.interpretation}")
        print()

        print("Per-Protein Results:")
        print("-" * 60)
        for pred in predictions:
            print(f"{pred['Protein']}:")
            print(f"  Serotype:       {pred['Serotype']}")
            print(f"  Top Hit:        {pred['Top Hit Accession']}")
            print(f"  Identity:       {pred['Percent Identity']}")
            print(f"  Query Coverage: {pred['Query Coverage']}")
            print(f"  Bitscore:       {pred['Bitscore']}")
            print(f"  E-value:        {pred['E-value']}")
            print()

        print("Output Files:")
        print(f"  Summary: {args.output}/btv_serotyping_summary.txt")
        print(f"  Predictions: {args.output}/btv_serotyping_predictions.tsv")
        print(f"  Detailed: {args.output}/btv_serotyping_detailed_results.tsv")
        print()

        return 0

    except Exception as e:
        print(f"Error during serotyping analysis: {e}", file=sys.stderr)
        if args.debug:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
