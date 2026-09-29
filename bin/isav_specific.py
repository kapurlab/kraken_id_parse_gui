"""Isavirus salaris (infectious salmon anemia virus) profile.

Now a thin declaration over the shared SegmentedEngine. Previously this module
carried a hand-copied ordering/download/concatenation implementation with the
old brittle matcher; it now uses the same strict-first + cautious-tentative
fallback as every other segmented organism (see organism_profiles.py).

NOTE: ported from the original slot table; validate against real ISAV
references when convenient (no ISAV sample was available during the refactor).
"""
from organism_profiles import (
    OrganismProfile, SegmentedEngine, make_strict_resolver, make_fallback_resolver,
)


# Influenza-style protein code -> ISAV segment (8 segments). NS1 and NEP are
# both encoded on segment 7.
ISAV_PROTEIN_TO_SEGMENT = {
    'PB2': 1, 'PB1': 2, 'NP': 3, 'PA': 4, 'F': 5,
    'HE': 6, 'NS1': 7, 'NEP': 7, 'M1': 8,
}

# Cosmetic display tokens for the concatenated-header segment name (verbatim
# from the original module).
ISAV_DISPLAY_ORDER = {
    '(PB2)': 1, '(PB1)': 2, '(NP)': 3, '(PA)': 4, '(F)': 5,
    '(HE)': 6, '(NS1)': 7, '(NEP)': 7, '(M1)': 8,
    'SEGMENT 1': 1, 'SEGMENT 2': 2, 'SEGMENT 3': 3, 'SEGMENT 4': 4,
    'SEGMENT 5': 5, 'SEGMENT 6': 6, 'SEGMENT 7': 7, 'SEGMENT 8': 8,
    'SEG 1': 1, 'SEG 2': 2, 'SEG 3': 3, 'SEG 4': 4, 'SEG 5': 5,
    'SEG 6': 6, 'SEG 7': 7, 'SEG 8': 8,
}

ISAV_PROFILE = OrganismProfile(
    name='Isavirus salaris',
    kind='segmented',
    n_segments=8,
    strict_resolver=make_strict_resolver(['SEGMENT', 'SEG'], ISAV_PROTEIN_TO_SEGMENT, 8),
    fallback_resolver=make_fallback_resolver(8),
    display_order=ISAV_DISPLAY_ORDER,
    unit_noun='segment',
    empty_placeholder='>No_Segments_Available\nNNNNNNNNNN\n',
    match_keywords=['infectious salmon anemia virus', 'isavirus', 'isav'],
)


class ISAV_Specific(SegmentedEngine):
    """Backward-compatible wrapper; all behavior lives in SegmentedEngine."""

    def __init__(self):
        super().__init__(ISAV_PROFILE)
