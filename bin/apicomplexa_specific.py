"""Apicomplexa (e.g. Theileria equi) profile — a multi-chromosome organism.

Demonstrates that the shared SegmentedEngine handles multi-chromosome genomes
just as it handles segmented viruses: chromosomes and organellar genomes are
simply the "slots". Now a thin declaration; previously a hand-copied module with
the old brittle matcher.

NOTE: ported from the original slot table; validate against real references when
convenient (no Apicomplexa sample was available during the refactor).
"""
from organism_profiles import (
    OrganismProfile, SegmentedEngine, make_strict_resolver, make_fallback_resolver,
)


# Organellar genome names -> slot. Nuclear chromosomes are matched via the
# "CHROMOSOME N" index notation (see strict_resolver below), not this table.
APICOMPLEXA_NAME_TO_SLOT = {
    'MITOCHONDRION': 5,
    'APICOPLAST': 6,
}

# Cosmetic display tokens for the concatenated-header name (verbatim from the
# original module).
APICOMPLEXA_DISPLAY_ORDER = {
    'CHROMOSOME 1': 1, 'CHROMOSOME 2': 2, 'CHROMOSOME 3': 3, 'CHROMOSOME 4': 4,
    'MITOCHONDRION': 5, 'APICOPLAST': 6,
}

APICOMPLEXA_PROFILE = OrganismProfile(
    name='Apicomplexa',
    kind='segmented',
    n_segments=6,
    strict_resolver=make_strict_resolver(
        ['CHROMOSOME', 'CHR', 'SEGMENT', 'SEG'], APICOMPLEXA_NAME_TO_SLOT, 6),
    fallback_resolver=make_fallback_resolver(6),
    display_order=APICOMPLEXA_DISPLAY_ORDER,
    unit_noun='sequence',
    empty_placeholder='>No_Sequences_Available\nNNNNNNNNNN\n',
    match_keywords=['theileria equi'],
)


class Apicomplexa(SegmentedEngine):
    """Backward-compatible wrapper; all behavior lives in SegmentedEngine."""

    def __init__(self):
        super().__init__(APICOMPLEXA_PROFILE)
