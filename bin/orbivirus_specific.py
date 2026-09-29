import os
import sys
import re
import time
from pathlib import Path
from collections import defaultdict

# Add package directory to Python path
# home = str(Path.home())
# sys.path.append(f'{home}/git/gitlab/vsnp3/bin')
from download_fasta_by_acc import Downloader
from organism_profiles import OrganismProfile, SegmentedEngine


# Canonical Orbivirus (BTV/EHD) protein-code -> genome segment number mapping.
# Serotype is determined by Segment 2 (VP2); Segment 6 is VP5.
BTV_PROTEIN_TO_SEGMENT = {
    'VP1': 1, 'VP2': 2, 'VP3': 3, 'VP4': 4, 'NS1': 5,
    'VP5': 6, 'VP7': 7, 'NS2': 8, 'VP6': 9, 'NS3': 10,
}


# Whole-word Orbivirus abbreviations, optionally with a serotype suffix
# ("BTV", "BTV-8", "BTV8", "EHD", "EHDV", "EHDV-2"). The word boundaries keep them
# from firing on unrelated text that merely contains the letters.
_BTV_ABBREV = re.compile(r'\bbtv(?:[-_]?\d{1,2})?\b')
_EHD_ABBREV = re.compile(r'\behdv?(?:[-_]?\d{1,2})?\b')


def orbivirus_species(description):
    """Return 'BTV', 'EHD' or None for a BLAST/FASTA description.

    The full species names (what GenBank descriptions use) are checked first, so
    an abbreviation that happens to appear elsewhere in a header (e.g. inside a
    sample name) can never override the real species; whole-word abbreviations
    are only the fallback. The ICTV binomial species names (Orbivirus
    caerulinguae = BTV, Orbivirus ruminantium = EHDV) are recognized too, since
    newer GenBank titles may use them (as ISAV titles already use "Isavirus
    salaris"). Shared by the pipeline (species split, report sections) and BTV
    serotyping so they never disagree."""
    d = str(description).lower()
    if 'bluetongue' in d or 'orbivirus caerulinguae' in d:
        return 'BTV'
    if 'epizootic hemorrhagic disease' in d or 'orbivirus ruminantium' in d:
        return 'EHD'
    if _BTV_ABBREV.search(d):
        return 'BTV'
    if _EHD_ABBREV.search(d):
        return 'EHD'
    return None


def resolve_segment_number(text):
    """Resolve an Orbivirus genome segment number (1-10) from a FASTA header,
    reference id, or description.

    This is the single source of truth for segment detection, shared by the
    coverage-graph segment-status banner and serotyping so they never disagree.
    It recognizes both explicit "segment N" / "seg N" notation and the protein
    codes (e.g. VP2 -> segment 2, VP5 -> segment 6). Returns None if no segment
    can be determined.
    """
    if not text:
        return None
    t = str(text).upper()

    # Explicit segment notation is most authoritative (handles "segment2",
    # "segment 2", "segment_2", "seg2", "seg 2"), including the size-class
    # form "segment S5" / "segment L2" used by some GenBank records (the L/M/S
    # letter is the segment's size class; the number is still the segment).
    m = re.search(r'SEGMENT[\s_]*[LMS]?(\d{1,2})', t)
    if not m:
        m = re.search(r'\bSEG[\s_]*[LMS]?(\d{1,2})\b', t)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 10:
            return n

    # Fall back to protein codes (VP1-VP7, NS1-NS3).
    for code, seg in BTV_PROTEIN_TO_SEGMENT.items():
        if re.search(r'\b' + code + r'\b', t):
            return seg

    return None


# Numeric contexts in Orbivirus descriptions that are NOT genome-segment
# numbers and must be ignored by the cautious fallback: serotype
# ("Bluetongue virus 11", "BTV-8", "serotype 4", "type 2", "genotype 3") and
# isolate / strain / clone identifiers. Stripping these before counting bare
# numbers keeps the fallback from mistaking a serotype (BTV-1..BTV-10) for a
# segment. Over-stripping is safe: it only makes the fallback decline to guess.
_NON_SEGMENT_CONTEXT = re.compile(
    r'(?:BLUETONGUE\s+VIRUS'
    r'|EPIZOOTIC\s+HEMORRHAGIC\s+DISEASE(?:\s+VIRUS)?'
    r'|BTV|EHDV|EHD|SEROTYPE|TYPE|ISOLATE|STRAIN|CLONE)'
    r'[\s:\-/]*\d+')


def fallback_segment_number(text):
    """Cautious, low-confidence fallback used ONLY to fill a segment the strict
    resolver (resolve_segment_number) could not place.

    Returns a segment number (1-10) only when the description contains exactly
    ONE unambiguous standalone number in that range, after removing obvious
    non-segment numeric contexts (serotype, isolate, strain, clone, type). If
    there are zero or conflicting candidates it returns None — deliberately
    refusing to guess. Callers must flag any result as 'tentative'.
    """
    if not text:
        return None
    t = str(text).upper()
    # Remove GenBank accession.version tokens (e.g. "AY636071.1", "NC_012345.2")
    # and any remaining ".N" version suffix so the version digit is never
    # counted as a segment number.
    t = re.sub(r'\b[A-Z]{1,2}_?\d{3,}(?:\.\d+)?', ' ', t)
    t = re.sub(r'\.\d+', ' ', t)
    t = _NON_SEGMENT_CONTEXT.sub(' ', t)
    candidates = set()
    for m in re.findall(r'\b(\d{1,2})\b', t):
        n = int(m)
        if 1 <= n <= 10:
            candidates.add(n)
    if len(candidates) == 1:
        return next(iter(candidates))
    return None


# Cosmetic display tokens for the segment name appended to concatenated headers.
# (Preserved verbatim from the original order_segments_by_coverage so the
# reference-guided FASTA descriptions are byte-identical to before.)
ORBIVIRUS_DISPLAY_ORDER = {
    '(VP1)': 1, '(VP2)': 2, '(VP3)': 3, '(VP4)': 4, '(NS1)': 5,
    '(VP5)': 6, '(VP7)': 7, '(NS2)': 8, '(VP6)': 9, '(NS3)': 10,
    'SEGMENT 1 ': 1, 'SEGMENT 2': 2, 'SEGMENT 3': 3, 'SEGMENT 4': 4,
    'SEGMENT 5': 5, 'SEGMENT 6': 6, 'SEGMENT 7': 7, 'SEGMENT 8': 8,
    'SEGMENT 9': 9, 'SEGMENT 10': 10,
    'SEG 1 ': 1, 'SEG 2': 2, 'SEG 3': 3, 'SEG 4': 4, 'SEG 5': 5,
    'SEG 6': 6, 'SEG 7': 7, 'SEG 8': 8, 'SEG 9': 9, 'SEG 10': 10,
}


# Declarative profile for Orbivirus (BTV/EHD): 10 segments, Orbivirus' own
# bespoke serotype-aware resolvers (unchanged), split into BTV vs EHD species.
ORBIVIRUS_PROFILE = OrganismProfile(
    name='Orbivirus',
    kind='segmented',
    n_segments=10,
    strict_resolver=resolve_segment_number,
    fallback_resolver=fallback_segment_number,
    display_order=ORBIVIRUS_DISPLAY_ORDER,
    unit_noun='segment',
    empty_placeholder='>No_Segments_Available\nNNNNNNNNNN\n',
    species_split=[
        (('bluetongue virus',), 'BTV'),
        (('epizootic hemorrhagic disease',), 'EHD'),
    ],
)


class Orbivirus_Specific(SegmentedEngine):
    """Thin subclass kept for backward compatibility with existing call sites
    (``Orbivirus_Specific().run(...)`` / ``.order_segments_by_coverage(...)``).
    All behavior lives in the shared SegmentedEngine, driven by ORBIVIRUS_PROFILE.
    """

    def __init__(self):
        super().__init__(ORBIVIRUS_PROFILE)
