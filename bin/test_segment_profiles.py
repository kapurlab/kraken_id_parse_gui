#!/usr/bin/env python3
"""Characterization + regression tests for per-organism segment resolution.

Purpose
-------
This is the safety net for the organism-profile refactor. It pins:

  1. ORBIVIRUS BEHAVIOR MUST NOT CHANGE. The refactor moves Orbivirus onto a
     shared engine; these cases assert the (accession -> segment, confidence)
     result is byte-identical to the pre-refactor behavior.

  2. ISAV / APICOMPLEXA GET THE FIX. These organisms previously used the old
     brittle matcher; after the refactor they must use the same strict-first +
     cautious-tentative-fallback logic as Orbivirus. Their cases assert the
     intended post-fix behavior.

Run:  python bin/test_segment_profiles.py
Exits non-zero on any failure so it can gate the refactor.
"""
import sys


def _summarize(ordered):
    """Reduce an order_segments_by_coverage result to {accession: (segment, confidence)}."""
    return {acc: (info['segment_number'], info.get('confidence', 'high'))
            for acc, info in ordered.items()}


# --- Representative inputs -------------------------------------------------- #
# Each value is a virus_dict as passed to order_segments_by_coverage:
#   {accession: {'header': <FASTA description>, 'percent_covered': <float>}}

ORBIVIRUS_INPUT = {
    'KM580469': {'header': 'KM580469.1 Bluetongue virus 11 isolate USA2012/XX 000000-12 VP1 protein (VP1) gene, complete cds', 'percent_covered': 99},
    'AY636071': {'header': 'AY636071.1 Bluetongue virus strain BT00L2-000 VP2 (L2) mRNA, complete cds', 'percent_covered': 98},
    'KM580445': {'header': 'KM580445.1 Bluetongue virus 11 isolate USA2011/XX 11-00000-18 VP3 protein (VP3) gene, complete cds', 'percent_covered': 97},
    'KM580419': {'header': 'KM580419.1 Bluetongue virus 11 isolate USA2011/XX 11-00000-3 VP4 protein (VP4) gene, complete cds', 'percent_covered': 96},
    'KX164083': {'header': 'KX164083.1 Bluetongue virus 12 isolate USA2008/XX 000000 segment 5, complete sequence', 'percent_covered': 90},
    'KM580487': {'header': 'KM580487.1 Bluetongue virus 11 isolate USA2013/XX 13-000001 VP7 protein (VP7) gene, complete cds', 'percent_covered': 92},
    'KM580452': {'header': 'KM580452.1 Bluetongue virus 11 isolate USA2011/XX 11-00000-18 NS2 protein (NS2) gene, complete cds', 'percent_covered': 91},
    'KX164107': {'header': 'KX164107.1 Bluetongue virus 13 isolate USA2013/XX 13-000002 segment 9, complete sequence', 'percent_covered': 89},
    # Oddly-labeled reference the strict resolver can't place; a single unambiguous
    # number (6) survives after serotype/isolate stripping -> tentative segment 6.
    'AY123456': {'header': 'AY123456.1 Bluetongue virus 12 isolate XYZ complete sequence 6', 'percent_covered': 80},
    # Only a serotype number (8) present -> stripped -> no candidate -> dropped.
    'KX999999': {'header': 'KX999999.1 Bluetongue virus 8 isolate ABC complete genome', 'percent_covered': 70},
}

# Golden result: Orbivirus behavior that must be preserved by the refactor.
ORBIVIRUS_GOLDEN = {
    'KM580469': (1, 'high'), 'AY636071': (2, 'high'), 'KM580445': (3, 'high'),
    'KM580419': (4, 'high'), 'KX164083': (5, 'high'), 'KM580487': (7, 'high'),
    'KM580452': (8, 'high'), 'KX164107': (9, 'high'),
    'AY123456': (6, 'tentative'),
    # KX999999 dropped (not in result)
}

ISAV_INPUT = {
    'AB000001': {'header': 'AB000001.1 Infectious salmon anemia virus segment 1 polymerase PB2 protein (PB2) gene', 'percent_covered': 99},
    'AB000003': {'header': 'AB000003.1 Infectious salmon anemia virus nucleoprotein (NP) gene', 'percent_covered': 98},
    'AB000005': {'header': 'AB000005.1 Infectious salmon anemia virus segment 5, complete sequence', 'percent_covered': 95},
    # Fusion protein bare code (F) -> segment 5 by table; but seg5 already taken by
    # AB000005 with higher coverage, so this is a lower-coverage duplicate.
    'AB000055': {'header': 'AB000055.1 Infectious salmon anemia virus fusion protein (F) gene', 'percent_covered': 50},
    # Oddly labeled, single unambiguous number 7 -> tentative segment 7.
    'AB000077': {'header': 'AB000077.1 Infectious salmon anemia virus isolate NOR complete sequence 7', 'percent_covered': 60},
}

APICOMPLEXA_INPUT = {
    'CP000001': {'header': 'CP000001.1 Theileria equi strain WA chromosome 1, complete sequence', 'percent_covered': 99},
    'CP000002': {'header': 'CP000002.1 Theileria equi strain WA chromosome 2, complete sequence', 'percent_covered': 98},
    'CP000005': {'header': 'CP000005.1 Theileria equi mitochondrion, complete genome', 'percent_covered': 88},
    'CP000006': {'header': 'CP000006.1 Theileria equi apicoplast, complete genome', 'percent_covered': 80},
}


def run_current():
    """Capture behavior of the CURRENT modules (pre-refactor baseline)."""
    from orbivirus_specific import Orbivirus_Specific
    from isav_specific import ISAV_Specific
    from apicomplexa_specific import Apicomplexa
    print('ORBIVIRUS   :', _summarize(Orbivirus_Specific().order_segments_by_coverage(ORBIVIRUS_INPUT)))
    print('ISAV        :', _summarize(ISAV_Specific().order_segments_by_coverage(ISAV_INPUT)))
    print('APICOMPLEXA :', _summarize(Apicomplexa().order_segments_by_coverage(APICOMPLEXA_INPUT)))


def _check(label, got, expected):
    ok = got == expected
    print(f'[{"PASS" if ok else "FAIL"}] {label}')
    if not ok:
        print(f'    expected: {expected}')
        print(f'    got     : {got}')
    return ok


def run_asserts():
    """Assert post-refactor behavior for all three organisms."""
    from orbivirus_specific import Orbivirus_Specific
    from isav_specific import ISAV_Specific
    from apicomplexa_specific import Apicomplexa

    orb = _summarize(Orbivirus_Specific().order_segments_by_coverage(ORBIVIRUS_INPUT))
    isav = _summarize(ISAV_Specific().order_segments_by_coverage(ISAV_INPUT))
    api = _summarize(Apicomplexa().order_segments_by_coverage(APICOMPLEXA_INPUT))

    isav_expected = {
        'AB000001': (1, 'high'),   # segment 1 / PB2
        'AB000003': (3, 'high'),   # NP
        'AB000005': (5, 'high'),   # segment 5 (higher coverage than the (F) dup)
        # Previously mis-merged to seg1 via the ".1" version artifact and lost;
        # the fix places it correctly as a tentative segment 7.
        'AB000077': (7, 'tentative'),
    }
    api_expected = {
        'CP000001': (1, 'high'), 'CP000002': (2, 'high'),
        'CP000005': (5, 'high'), 'CP000006': (6, 'high'),
    }

    from orbivirus_specific import resolve_segment_number, orbivirus_species
    # GenBank size-class notation ("segment S5") must resolve like "segment 5".
    size_class = {t: resolve_segment_number(t) for t in (
        'MT013328.1 Epizootic hemorrhagic disease virus strain JC00C000 segment S5',
        'MT013321.1 Epizootic hemorrhagic disease virus strain JC00C001 segment S8',
        'Bluetongue virus segment L2', 'complete genome segments 1-10')}
    size_class_expected = {
        'MT013328.1 Epizootic hemorrhagic disease virus strain JC00C000 segment S5': 5,
        'MT013321.1 Epizootic hemorrhagic disease virus strain JC00C001 segment S8': 8,
        'Bluetongue virus segment L2': 2, 'complete genome segments 1-10': None}
    # Species split used by the report sections and BTV serotyping (mixed samples).
    species = {t: orbivirus_species(t) for t in (
        'S1_segment2 guided by KM580445.1 Bluetongue virus 11 VP3',
        'S1_segment2 guided by MH845320.1 Epizootic hemorrhagic disease virus 2 segment 2',
        'EHD-2024 guided by AY636071.1 Bluetongue virus strain BT00L2-000 VP2',
        'BTV-8 isolate X', 'EHDV-2 isolate Y', 'Whispovirus xiabaidian',
        'Orbivirus caerulinguae isolate Z segment 2', 'Orbivirus ruminantium isolate W VP7')}
    species_expected = {
        'S1_segment2 guided by KM580445.1 Bluetongue virus 11 VP3': 'BTV',
        'S1_segment2 guided by MH845320.1 Epizootic hemorrhagic disease virus 2 segment 2': 'EHD',
        'EHD-2024 guided by AY636071.1 Bluetongue virus strain BT00L2-000 VP2': 'BTV',
        'BTV-8 isolate X': 'BTV', 'EHDV-2 isolate Y': 'EHD', 'Whispovirus xiabaidian': None,
        'Orbivirus caerulinguae isolate Z segment 2': 'BTV',
        'Orbivirus ruminantium isolate W VP7': 'EHD'}

    results = [
        _check('Orbivirus behavior preserved (golden)', orb, ORBIVIRUS_GOLDEN),
        _check('ISAV strict + tentative fallback', isav, isav_expected),
        _check('Apicomplexa chromosomes/organelles', api, api_expected),
        _check('Segment size-class notation (segment S5 / L2)', size_class, size_class_expected),
        _check('BTV/EHD species split', species, species_expected),
    ]
    return all(results)


if __name__ == '__main__':
    if '--capture' in sys.argv:
        run_current()
    else:
        sys.exit(0 if run_asserts() else 1)
