#!/usr/bin/env python3
"""Registry of organism profiles keyed by --taxon value.

This is the single place to register an organism for the segmented /
multi-chromosome workflow. To add a new segmented virus or multi-chromosome
organism:

  1. Declare an OrganismProfile in its own module (see isav_specific.py /
     apicomplexa_specific.py for minimal examples), giving it match_keywords so
     the dispatcher can pick its BLAST hits.
  2. Import it here and add it to PROFILES.

That's it — the generic dispatcher in kraken_id_parse.py picks it up
automatically; no new orchestration code is required.

Note: Orbivirus is registered here for completeness, but it is still handled by
a dedicated branch in kraken_id_parse.py because it needs BTV/EHD species
splitting and BTV serotyping. Its profile carries that metadata (species_split)
so it can be migrated onto the generic handler once validated end-to-end.
"""
from orbivirus_specific import ORBIVIRUS_PROFILE
from isav_specific import ISAV_PROFILE
from apicomplexa_specific import APICOMPLEXA_PROFILE

PROFILES = {
    p.name: p for p in (
        ORBIVIRUS_PROFILE,
        ISAV_PROFILE,
        APICOMPLEXA_PROFILE,
    )
}

# Taxa handled by a dedicated (bespoke) branch rather than the generic handler.
BESPOKE_TAXA = {'Orbivirus'}
