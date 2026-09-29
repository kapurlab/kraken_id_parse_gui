#!/usr/bin/env python3
"""The Orbivirus (BTV/EHD) analysis is beta: untested, and its results cannot
be trusted.

One wording, used everywhere a BTV/EHD result reaches a person — both reports
(HTML and PDF), the stats workbook, the run and summary logs, the standalone
serotyping output, and the GUI's run form and Results pane — so nobody meets a
BTV/EHD result without the caveat beside it. Remove it (and its callers) only
once the analysis has been validated.
"""
import re

TITLE = 'BETA — Orbivirus (BTV/EHD) analysis is untested'
TEXT = ('The Bluetongue virus / Epizootic hemorrhagic disease virus analysis — the '
        'BTV/EHD split, segment assignments and BTV serotype — is still beta: it is '
        'untested and needs further testing. Its results cannot be trusted. Confirm '
        'them by a validated method before acting on them.')
SHORT = 'BETA: untested, needs further testing; results cannot be trusted'
# Beside the serotype call itself, in the report's BTV Serotyping section.
SEROTYPE = 'This serotype call is untested and needs further testing; it cannot be trusted.'
# One line for the top margin of every PDF page, so a lone printed page says it.
PAGE = 'BETA — Orbivirus (BTV/EHD) analysis untested; results cannot be trusted'
# The Results pane shows a one-line reason beside its REVIEW chip (48 chars max).
CHIP = 'BETA: Orbivirus (BTV/EHD) results untested'

# A taxon that selects the BTV/EHD workflow: the Orbivirus genus, either species
# by common or ICTV name, or their abbreviations.
_TAXON = re.compile(
    r'\borbivirus\b|\bbluetongue\b|\bepizootic\s+ha?emorrhagic\b|\bBTV\b|\bEHDV?\b',
    re.IGNORECASE)


def applies_to_taxon(taxon):
    """True when a run's target taxon selects the BTV/EHD (Orbivirus) workflow."""
    return bool(taxon) and bool(_TAXON.search(str(taxon)))


def notice():
    """The notice as stored in run_manifest.json."""
    return {'title': TITLE, 'text': TEXT, 'short': SHORT}


def notice_for_manifest(manifest):
    """The notice for a run_manifest.json, or None.

    A manifest written since the notice existed records it (the pipeline also
    flags BTV/EHD results found under a broader taxon such as Viruses). Older
    Orbivirus runs are recognized by their target taxon, so they are not shown
    as trustworthy just because they predate it."""
    manifest = manifest or {}
    if manifest.get('beta_notice'):
        return manifest['beta_notice']
    params = manifest.get('parameters') or {}
    target = ((manifest.get('legacy_report') or {}).get('target') or {})
    if applies_to_taxon(params.get('taxon')) or applies_to_taxon(target.get('taxon')):
        return notice()
    return None
