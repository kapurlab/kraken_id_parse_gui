# Custom Orbivirus-enriched Kraken2 database

This directory holds the **custom Orbivirus (BTV/EHD) Kraken2 database** used by the
`*_orbivirus_custom` presets. The heavy binaries (`*.k2d`, Bracken `*.kmer_distrib`,
`_build_tmp/`, `taxonomy/`, `library/`) are **git-ignored** (see `.gitignore`) — they are
rebuilt locally from committed inputs by `bin/build_orbivirus_kraken_db.sh`.

## Why this DB exists

The pipeline extracts the **entire `Orbivirus` genus clade** (`extract_kraken_reads.py
--include-children`) before assembly/serotyping, so the number of reads that reach VP2
assembly is a direct function of how many reads Kraken classifies as Orbivirus — i.e. of
**database content**. The stock `k2_standard_08gb` is downsampled to an 8 GB cap and built
only from sparse RefSeq orbivirus diversity, so it badly under-extracts BTV/EHD. This DB is
**uncapped** and orbivirus-rich, recovering far more BTV/EHD reads and segments.

## What it contains

| Source | Taxid tagging | Notes |
|---|---|---|
| Entire NCBI **Orbivirus** genus (GenBank + RefSeq) | species-level (from `datasets` metadata) | the enrichment that drives recall |
| **RefSeq viral** release background | coarse `Viruses` (`10239`); orbivirus accessions excluded first | non-orbivirus context without coarsening orbivirus LCA |
| `database/source_sequences/BTV_sequence.fasta` | Bluetongue virus `40051` | committed lab BTV segment panel (S1–S10) |
| `database/source_sequences/EHD_sequence.fasta` | EHDV `40054` | committed lab EHDV segment panel (S1–S10) |
| Optional `--panel-dir` lab consensuses | `40051` (or `--panel-taxid`) | reference-guided assemblies, if supplied |

## How it was built

All downloads go over **HTTPS** (NCBI rsync port 873 and FTP are blocked on the build host);
`kraken2-build` is used only for `--add-to-library` / `--build`. Bracken distributions are
built **before** `kraken2-build --clean` (clean removes `library/` + `taxonomy/`, which
`bracken-build` needs). The pipeline calls `bracken ... -r 250`, so
`database250mers.kmer_distrib` is **required**.

```bash
conda activate kraken_id_parse
# Builds INTO this directory by default (database/kraken_orbivirus):
bin/build_orbivirus_kraken_db.sh

# Optional: also inject reference-guided lab consensuses
bin/build_orbivirus_kraken_db.sh --panel-dir /path/to/btv_fastas
```

The build re-downloads the Orbivirus genus (~22k seqs) + RefSeq viral over HTTPS, so it takes
minutes to a few hours depending on the network. It is idempotent: re-running reuses cached
downloads in `_build_tmp/`.

### Runtime footprint
Only `*.k2d` (~680 MB) + `database250mers.kmer_distrib` are needed at runtime; `_build_tmp/`
is a removable download cache.

## How to add additional genomes later

1. **A handful of curated sequences (BTV/EHD lab references):** append them to the matching
   committed FASTA and rebuild:
   ```bash
   cat new_btv_segments.fasta >> database/source_sequences/BTV_sequence.fasta   # BTV -> taxid 40051
   cat new_ehd_segments.fasta >> database/source_sequences/EHD_sequence.fasta   # EHD -> taxid 40054
   bin/build_orbivirus_kraken_db.sh
   ```
   Headers need no `|kraken:taxid|` — the build script tags BTV with `40051` and EHD with
   `40054` automatically. For a different organism, add a FASTA and pass it with the right
   taxid via `--btv-ref` / `--ehd-ref`, or inject it as a `--panel-dir` with `--panel-taxid`.

2. **A whole new taxon clade from NCBI:** extend the `fetch_and_tag <taxon>` calls in
   `bin/build_orbivirus_kraken_db.sh` (the Orbivirus genus is fetched this way) and rebuild.

After any change, sanity-check:
```bash
ls database/kraken_orbivirus/hash.k2d database/kraken_orbivirus/database250mers.kmer_distrib
bin/build_orbivirus_kraken_db.sh --test-r1 R1.fastq.gz --test-r2 R2.fastq.gz   # classifies + greps Orbivirus
```

## Transfer to the HPC

The built DB is self-contained. Copy it (excluding the cache) and point an HPC preset's
`kraken_db` at the destination:
```bash
rsync -av --exclude _build_tmp database/kraken_orbivirus/ user@hpc:/dest/k2_orbi_virus/
# then set kraken_db: /dest/k2_orbi_virus in internal/scomp_kraken_configs.yaml (orbivirus_custom)
```
