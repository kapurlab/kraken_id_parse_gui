# `database/` — self-contained Orbivirus (BTV/EHD) reference data

This folder makes the Orbivirus workflow reproducible from inside the repo. It holds the
**source reference sequences** and the **documentation/scripts** needed to (re)build the
databases the pipeline uses for Bluetongue virus (BTV) and Epizootic hemorrhagic disease
virus (EHDV).

```
database/
  source_sequences/        # COMMITTED FASTA inputs (the build inputs live in git)
    BTV_sequence.fasta      # BTV genome segments S1–S10 (serotype/segment-labeled)
    EHD_sequence.fasta      # EHDV genome segments S1–S10 (GenBank + lab consensuses)
  kraken_orbivirus/        # custom Orbivirus-enriched Kraken2 DB (built here, binaries git-ignored)
    README.md               # what it contains, how it is built, how to add genomes
    .gitignore              # keeps the heavy *.k2d binaries out of git
  blast_orbivirus/         # BLAST side: how to add these sequences to nt_viruses
    README.md
```

## Why the layout is split this way

| Artifact | Size | In git? | How it's reproduced |
|---|---|---|---|
| `source_sequences/*.fasta` | small | **yes** | committed directly |
| Kraken2 DB (`kraken_orbivirus/*.k2d`) | ~680 MB | no (git-ignored) | `bin/build_orbivirus_kraken_db.sh` rebuilds it here from the committed FASTAs + NCBI downloads |
| BLAST `nt_viruses` | ~70 GB | no (external) | stays external; `bin/add_orbivirus_to_blast.sh` adds these sequences to it — see `blast_orbivirus/README.md` |

The Kraken DB binaries are **not** committed (no Git LFS; ~680 MB would bloat history). Instead the
repo carries everything needed to rebuild an identical DB locally into `kraken_orbivirus/`. The
BLAST `nt_viruses` database is far too large to live in the repo, so the repo carries the source
FASTAs plus a documented procedure to add them to whatever `nt_viruses` a host already has.

## Quick start

```bash
conda activate kraken_id_parse_env

# (Re)build the custom Orbivirus Kraken DB into database/kraken_orbivirus/
bin/build_orbivirus_kraken_db.sh                 # see database/kraken_orbivirus/README.md

# Add the BTV/EHD reference sequences to a local nt_viruses BLAST DB
bin/add_orbivirus_to_blast.sh                    # see database/blast_orbivirus/README.md
```

The BTV **serotyping** panel (VP2/VP5, serotypes 1–27) lives separately under the repo's
top-level `reference_sequences/` and is consumed directly by `bin/btv_serotyping.py` (no
built database needed). See that script and `reference_sequences/BTV_panel_manifest.tsv`.
