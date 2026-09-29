#!/usr/bin/env python3
"""Shared engine and declarative profiles for segmented / multi-chromosome
organisms.

Background
----------
Historically each segmented organism (Orbivirus, Isavirus, Apicomplexa) had its
own hand-copied ``*_specific.py`` with an identical ``order_segments_by_coverage
/ concatenate_fasta_files / run`` trio that differed ONLY in a small slot table
and a segment count. The copies drifted (only Orbivirus got the strict-first +
cautious-tentative-fallback fix), so oddly-labeled references were silently
dropped or mis-assigned for the others.

This module removes the duplication:

  * ``SegmentedEngine`` holds the one, correct implementation of segment
    ordering, downloading, and concatenation (lifted verbatim from the fixed
    Orbivirus code, parametrized by an ``OrganismProfile``).

  * ``OrganismProfile`` is the per-organism DATA: how many slots, how to resolve
    a slot number from a reference description (strict + optional cautious
    fallback), and cosmetic labels.

Adding a new segmented virus or multi-chromosome organism is then a matter of
declaring one ``OrganismProfile`` — no new orchestration code.

Design note (safety): Orbivirus keeps its own, unchanged resolver functions
(``orbivirus_specific.resolve_segment_number`` / ``fallback_segment_number``) and
its original display table, so routing it through this engine is a pure code
move with byte-identical output (see bin/test_segment_profiles.py). Only the
previously-stale organisms change behavior — they gain the fix.
"""
import os
import re
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

from download_fasta_by_acc import Downloader


# --------------------------------------------------------------------------- #
# Generic resolver factories (used by organisms that do NOT need Orbivirus'
# bespoke serotype-aware resolver — i.e. Isavirus, Apicomplexa, and future ones).
# --------------------------------------------------------------------------- #

# Numeric contexts that are never a slot number: isolate / strain / clone / type
# identifiers, plus GenBank accession.version tokens. Stripping these before
# counting bare numbers keeps the cautious fallback from mistaking, e.g., an
# isolate year or an accession version for a segment/chromosome number.
_GENERIC_NON_SLOT_CONTEXT = re.compile(
    r'(?:ISOLATE|STRAIN|CLONE|TYPE|SEROTYPE)[\s:\-/]*\d+')


def make_strict_resolver(index_keywords, name_to_slot, n_slots):
    """Build a strict (high-confidence) slot resolver.

    index_keywords: words that precede an explicit slot index, e.g.
        ['SEGMENT', 'SEG'] for viruses or ['CHROMOSOME', 'CHR'] for genomes.
    name_to_slot:   exact tokens mapped to a slot, e.g. {'PB2': 1, 'NP': 3} or
        {'MITOCHONDRION': 5, 'APICOPLAST': 6}. Matched on word boundaries.
    n_slots:        maximum valid slot number.

    Returns a function(text) -> int | None. Explicit "KEYWORD N" notation wins;
    otherwise a name token; otherwise None.
    """
    kw_patterns = [re.compile(kw + r'[\s_]*(\d{1,2})') for kw in index_keywords]
    name_items = list(name_to_slot.items())

    def resolve(text):
        if not text:
            return None
        t = str(text).upper()
        for pat in kw_patterns:
            m = pat.search(t)
            if m:
                n = int(m.group(1))
                if 1 <= n <= n_slots:
                    return n
        for token, slot in name_items:
            if re.search(r'\b' + re.escape(token) + r'\b', t):
                return slot
        return None

    return resolve


def make_fallback_resolver(n_slots, extra_context=None):
    """Build a cautious, low-confidence fallback slot resolver.

    Returns a slot number (1..n_slots) only when the description contains exactly
    ONE unambiguous standalone number in range after removing accession.version
    tokens and obvious non-slot numeric contexts (isolate/strain/clone/type, plus
    any organism-specific ``extra_context`` regex). Zero or conflicting
    candidates -> None (refuses to guess). Callers must flag results 'tentative'.
    """
    def resolve(text):
        if not text:
            return None
        t = str(text).upper()
        # Drop GenBank accession.version (e.g. "AB000077.1", "NC_012345.2") and
        # any remaining ".N" version suffix so version digits are never counted.
        t = re.sub(r'\b[A-Z]{1,2}_?\d{3,}(?:\.\d+)?', ' ', t)
        t = re.sub(r'\.\d+', ' ', t)
        t = _GENERIC_NON_SLOT_CONTEXT.sub(' ', t)
        if extra_context is not None:
            t = extra_context.sub(' ', t)
        candidates = set()
        for m in re.findall(r'\b(\d{1,2})\b', t):
            n = int(m)
            if 1 <= n <= n_slots:
                candidates.add(n)
        if len(candidates) == 1:
            return next(iter(candidates))
        return None

    return resolve


# --------------------------------------------------------------------------- #
# Organism profile (declarative data)
# --------------------------------------------------------------------------- #

@dataclass
class OrganismProfile:
    """Declarative description of how to handle one organism.

    name:            value compared against --taxon.
    kind:            'segmented' (multi-slot: segments or chromosomes) or
                     'single' (one reference; handled by the generic path).
    n_segments:      number of slots (segments / chromosomes) for 'segmented'.
    strict_resolver: fn(description) -> slot int | None (high confidence).
    fallback_resolver: fn(description) -> slot int | None (tentative) or None to
                     disable the cautious fallback for this organism.
    display_order:   {display_token: slot} used ONLY for the cosmetic segment
                     name appended to concatenated headers.
    unit_noun:       'segment' or 'sequence' — wording in progress messages.
    empty_placeholder: FASTA record written when nothing downloads.
    match_keywords:  LOWERCASE substrings identifying this organism's BLAST hits
                     (matched against a lowercased FASTA header), e.g.
                     ['infectious salmon anemia virus', 'isav']. Used by the
                     generic dispatcher to pick which hits to process.
    species_split:   optional list of (match_substrings, label) groups used to
                     split BLAST hits into sub-species (e.g. BTV vs EHD). None
                     means treat all hits as one group.
    serotyper:       optional callable/class for serotyping (e.g. BTV). Consumed
                     by the dispatcher, not the engine.
    """
    name: str
    kind: str = 'segmented'
    n_segments: int = 10
    strict_resolver: Optional[Callable] = None
    fallback_resolver: Optional[Callable] = None
    display_order: Dict[str, int] = field(default_factory=dict)
    unit_noun: str = 'segment'
    empty_placeholder: str = '>No_Segments_Available\nNNNNNNNNNN\n'
    match_keywords: list = field(default_factory=list)
    species_split: Optional[list] = None
    serotyper: Optional[Callable] = None


# --------------------------------------------------------------------------- #
# Shared engine (one correct implementation for all segmented organisms)
# --------------------------------------------------------------------------- #

class SegmentedEngine:
    """Segment ordering + download + concatenation, driven by an OrganismProfile.

    This is the single copy of what used to be duplicated across every
    ``*_specific.py``. The logic matches the fixed Orbivirus implementation:
    strict (high-confidence) assignment first, then a cautious tentative
    fallback that only fills MISSING slots from references the strict pass could
    not place.
    """

    def __init__(self, profile: OrganismProfile):
        self.profile = profile

    def order_segments_by_coverage(self, virus_dict):
        """Return {seq_id: seq_info} with 'segment_number', 'confidence'
        ('high'|'tentative') and 'segment_name' set; one (highest-coverage)
        reference per slot."""
        p = self.profile

        # Tier 1: strict, high-confidence.
        strict = defaultdict(list)
        unassigned = []
        for seq_id, seq_info in virus_dict.items():
            slot = p.strict_resolver(seq_info['header']) if p.strict_resolver else None
            if slot:
                strict[slot].append((seq_id, seq_info))
            else:
                unassigned.append((seq_id, seq_info))

        # Tier 2: cautious fallback — only fill MISSING slots, only from
        # references the strict pass could not place, only unambiguous matches.
        fallback = defaultdict(list)
        missing = set(range(1, p.n_segments + 1)) - set(strict.keys())
        if missing and p.fallback_resolver:
            for seq_id, seq_info in unassigned:
                slot = p.fallback_resolver(seq_info['header'])
                if slot and slot in missing:
                    fallback[slot].append((seq_id, seq_info))

        ordered_dict = {}
        for slot in range(1, p.n_segments + 1):
            if slot in strict:
                candidates, confidence = strict[slot], 'high'
            elif slot in fallback:
                candidates, confidence = fallback[slot], 'tentative'
            else:
                continue

            best = max(candidates, key=lambda x: x[1]['percent_covered'])
            seq_info = best[1].copy()
            seq_info['segment_number'] = slot
            seq_info['confidence'] = confidence

            header = best[1]['header'].upper()
            segment_name = f"Segment {slot}"
            for seg_name, seg_num in p.display_order.items():
                if seg_num == slot and seg_name in header:
                    segment_name = seg_name
                    break
            seq_info['segment_name'] = segment_name

            if confidence == 'tentative':
                print(f"  ⚠️  Tentative (low-confidence) {p.unit_noun} assignment: "
                      f"{best[0]} → segment {slot} (from '{best[1]['header']}')")
            ordered_dict[best[0]] = seq_info

        return ordered_dict

    def concatenate_fasta_files(self, ordered_dict, output_file="concatenated_specific.fasta"):
        """Concatenate downloaded per-slot FASTA files, appending segment info to
        each header (identical format to the original ``*_specific.py``)."""
        with open(output_file, 'w') as outfile:
            for seq_id, seq_info in ordered_dict.items():
                input_file = f"{seq_id}.fasta"
                try:
                    with open(input_file, 'r') as infile:
                        content = infile.read().strip()
                        parts = content.split('\n', 1)
                        if len(parts) != 2:
                            print(f"Warning: Unexpected format in {input_file}")
                            continue
                        header, sequence = parts
                        segment_info = f"segment_{seq_info['segment_number']}"
                        if 'segment_name' in seq_info:
                            segment_info = f"{segment_info}_{seq_info['segment_name']}"
                        outfile.write(f"{header} {segment_info}\n")
                        outfile.write(f"{sequence}\n")
                except FileNotFoundError:
                    print(f"Warning: Could not find file for {seq_id}")
                    continue
        print(f"Concatenation complete. Output written to {output_file}")

    def run(self, alignment_stats=None):
        """Order, download (with retries + rate limiting), and concatenate.
        Returns the concatenated {seq_id: seq_info} dict (incl. per-slot
        'confidence') so callers can carry the high/tentative flag downstream."""
        p = self.profile
        dict_ordered = self.order_segments_by_coverage(alignment_stats)

        print(f"Downloading {len(dict_ordered)} {p.name} {p.unit_noun}s...")
        downloaded_count = 0
        failed_downloads = []

        for i, (seq_id, seq_inf) in enumerate(dict_ordered.items(), 1):
            print(f"Downloading {p.unit_noun} {i}/{len(dict_ordered)}: {seq_id}")
            if i > 1:
                print("  Waiting 2 seconds (NCBI rate limiting)...")
                time.sleep(2)

            downloader = Downloader(seq_id)
            max_retries = 3
            success = False
            for attempt in range(max_retries):
                try:
                    downloader.fasta()
                    success = True
                    downloaded_count += 1
                    print(f"  ✓ Successfully downloaded {seq_id}")
                    break
                except Exception as e:
                    print(f"  ⚠ Attempt {attempt + 1} failed for {seq_id}: {e}")
                    if attempt < max_retries - 1:
                        wait_time = (attempt + 1) * 5
                        print(f"    Waiting {wait_time} seconds before retry...")
                        time.sleep(wait_time)

            if not success:
                failed_downloads.append(seq_id)
                print(f"  ✗ Failed to download {seq_id} after {max_retries} attempts")

        print(f"\nDownload summary:")
        print(f"  Successfully downloaded: {downloaded_count}/{len(dict_ordered)}")
        if failed_downloads:
            print(f"  Failed downloads: {len(failed_downloads)}")
            print(f"  Failed accessions: {', '.join(failed_downloads)}")

        if downloaded_count > 0:
            print(f"Proceeding with concatenation of successfully downloaded {p.unit_noun}s...")
            successful_dict = {k: v for k, v in dict_ordered.items() if k not in failed_downloads}
            self.concatenate_fasta_files(successful_dict)
            for seq_id, seq_inf in successful_dict.items():
                try:
                    os.remove(f'{seq_id}.fasta')
                except FileNotFoundError:
                    pass
            return successful_dict
        else:
            print(f"Warning: No {p.unit_noun}s were successfully downloaded")
            with open("concatenated_specific.fasta", 'w') as f:
                f.write(p.empty_placeholder)
            return {}
