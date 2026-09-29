# Kraken ID and Parse

## Running (quick start)

Run from a directory containing **one sample's** paired FASTQs (`*_R1*.fastq.gz` / `*_R2*.fastq.gz`).

Each run writes an interactive **HTML** report (`*_report.html`, with zoomable coverage
graphs) and a **PDF** rendered from it (`*_report.pdf`), plus `*_stats.xlsx`. No LaTeX/TeX
install is required — the PDF is produced by WeasyPrint.

There are two ways to run, both driving the **same** Python pipeline: a single sample by
hand, or through **SLURM** (one job per sample). Pick whichever fits.

**Single sample, direct invocation** — `cd` into the sample's working directory (the one
containing its paired FASTQs), activate the environment, then run:
```bash
conda activate kraken_id_parse

python ${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 SAMPLE_R1.fastq.gz \
  -r2 SAMPLE_R2.fastq.gz \
  --taxon "Orbivirus" \
  --kraken_db /path/to/kraken2_db \
  --blast_db /path/to/blast_db \
  --logo /path/to/logo.png
```
`--taxon` is whatever you're targeting (e.g. `"Mycobacterium tuberculosis complex"`, `"Orbivirus"`,
`"Viruses"`) — see [Overview](#overview) for how taxon-specific handling works. `--logo` is
optional. Outputs land in the working directory: `SAMPLE_<stamp>_report.html`,
`SAMPLE_<stamp>_report.pdf`, `SAMPLE_<stamp>_stats.xlsx`, `SAMPLE_<stamp>_summary.log`.

**Local (macOS dev) — Orbivirus/BTV, custom DB + serotyping (preset-based, this repo's dev setup):**
```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate kraken_id_parse && \
export BLASTDB=/Users/todstuber/databases/blast && \
python /Users/todstuber/git/temp/kraken_id_parse/bin/run_with_config.py --preset macos_dev_orbivirus_custom
```

**HPC (SLURM):**
```bash
sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus
# many samples (one per subdirectory):
packagefastqs.sh; d=$(pwd); for f in */; do (cd "$f" && sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus); done
```
Edit the 3 settings at the top of `internal/kraken_id_parse.slurm` (`CONDA_PATH`, `CONDA_ENV`, `REPO_ROOT`) once for your HPC. Presets are in `internal/kraken_configs.yaml` (local) / `internal/scomp_kraken_configs.yaml` (HPC).

**BTV serotyping only (standalone, segment 2 / VP2, serotypes 1–27) — BETA, untested; results cannot be trusted:**
```bash
python bin/btv_serotyping.py -i SAMPLE_reference_guided.fasta -o ./sero_out
```

**First time / Kraken DB missing?** Build the custom Orbivirus DB once (details: [`database/`](database/README.md)):
```bash
bin/build_orbivirus_kraken_db.sh      # builds into database/kraken_orbivirus/
```

---

## Repository Setup
Before using the pipeline, set up the repository root path for easier command execution:

```bash
# Set repository root location for convenience
export REPO_ROOT="${HOME}/git/gitlab/kraken_id_parse"
```

This variable will be used throughout this documentation to reference the repository location.

## Overview
This pipeline performs taxonomic read filtering, assembly, BLAST analysis, and coverage visualization of WGS data. It consists of several interconnected scripts that can be run individually or as a complete workflow using the wrapper script.

Reads can be extracted by providing a taxonomical name using the `--taxon` option. Additionally, species-specific functions are available that can further parse difficult to distinguish organisms. For example, an Orbivirus function is used when the `--taxon` search "Orbivirus" is called. This will further distinguish reads as Bluetongue Virus or Epizootic Hemorrhagic Disease.

> **BETA — the Orbivirus (BTV/EHD) analysis is untested.** The BTV/EHD split, segment
> assignments and BTV serotype need further testing, and their results cannot be trusted:
> confirm them by a validated method before acting on them. Every report, the stats workbook
> and the GUI say so for each Orbivirus run (`bin/orbivirus_beta.py` holds the wording).

Upon completion, an interactive HTML report, a matching PDF, and an Excel stats workbook are
created to summarize the run (see [Running (quick start)](#running-quick-start) above for exact
filenames).

## Workflow
1. Kraken
2. Reads parsed on taxon
3. Parsed reads assembled
4. Assembly identified with BLAST
   - For segmented organisms (Orbivirus, ISAV, Apicomplexa), any contig stretch its top hit does
     not explain (e.g. two segments the assembler fused into one contig) is BLASTed on its own
5. Download FASTAs from BLAST findings
6. Coverage graph on downloaded FASTAs
7. If available continue with taxon specific workflow:
   - If Orbivirus, split the BLAST hits into Bluetongue virus and Epizootic hemorrhagic disease virus
     (common or ICTV binomial names, e.g. *Orbivirus caerulinguae*)
     - Order by segment and make coverage graph
8. Reference guided assembly using the original FASTQ files
9. BLAST reference guided assembly
10. Alignment using top BLAST results
11. Merged VCF to the top BLAST results as the final consensus
12. Make final coverage graph

# Installation

In the suite, `bdtools install kraken_id_parse_gui` builds this tool's environment from
[`conda_setup/environment.yml`](./conda_setup/environment.yml) (environment `kraken_id_parse`);
standalone, follow [conda_setup/SETUP_INSTRUCTIONS.md](./conda_setup/SETUP_INSTRUCTIONS.md).
That one environment holds Kraken2, SPAdes, BLAST, BWA, samtools, freebayes, the HTML/PDF
report stack and the GUI's web backend. Bracken is optional: where it cannot share that
environment (the osx-64 Bracken 2.x builds link an old zlib), it can live in its own small
environment, [`conda_setup/environment.bracken.yml`](./conda_setup/environment.bracken.yml) — see
[`conda_setup/conda_setup.md`](./conda_setup/conda_setup.md). The pipeline finds a sibling `bracken`
environment, `$BRACKEN_BIN`, or `bracken` on PATH, and skips the Bracken abundance chart (with a
warning) when there is none.

Quick setup:
```bash
# Create environment (use conda, not mamba on macOS)
CONDA_OVERRIDE_OSX=11.0 conda env create -f conda_setup/environment.yml
# Optional: Bracken in its own environment (named "bracken")
conda env create -f conda_setup/environment.bracken.yml

# Activate environment
conda activate kraken_id_parse
```

On Apple Silicon Macs, bioconda has no native `osx-arm64` builds for several required tools
(Kraken2, Bracken, SPAdes); build both environments under Rosetta instead:
```bash
CONDA_SUBDIR=osx-64 mamba env create -f conda_setup/environment.yml
CONDA_SUBDIR=osx-64 mamba env create -f conda_setup/environment.bracken.yml
```

Then activate the **main** environment before running the pipeline (the `bracken` environment is
found automatically and never activated directly):
```bash
mamba activate kraken_id_parse
```

Full setup details (matplotlib backend config, verifying the install, etc.): [`conda_setup/conda_setup.md`](./conda_setup/conda_setup.md)

## Prerequisites

### Kraken Database
1. Download Prebuilt Kraken Database from [Prebuilt databases](https://benlangmead.github.io/aws-indexes/k2)

Example:
```bash
cd ${HOME}
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
```

2. Extract the database:
```bash
mkdir k2_standard_08gb
tar -xzf k2_standard_08gb_*.tar.gz -C k2_standard_08gb
```

Note: The large Standard or core_nt Database Collections are preferred, but their size makes downloading and memory usage limiting.

3. Create alias to link taxonomy file to conda environment:
```bash
rm -rf ${HOME}/miniconda3/envs/kraken_id_parse/opt/krona/taxonomy
ln -s ${HOME}/k2_standard_08gb ${HOME}/miniconda3/envs/kraken_id_parse/opt/krona/taxonomy
cd ${HOME}/k2_standard_08gb
ktUpdateTaxonomy.sh
```

#### Build a custom Orbivirus-enriched Kraken DB (optional, recommended for BTV/EHD)

The small prebuilt `k2_standard_08gb` is downsampled to an 8 GB minimizer cap and built only from
RefSeq, so it classifies far fewer Orbivirus reads than the large `nt_core` DB — leaving BTV/EHD
segments (notably VP2) under-assembled. `bin/build_orbivirus_kraken_db.sh` builds an **uncapped**
DB = RefSeq **viral** library + the **entire NCBI *Orbivirus* genus** + (optionally) your lab's
validated consensus sequences, which recovers many more orbivirus reads with no pipeline changes.

```bash
conda activate kraken_id_parse            # needs ncbi-datasets-cli + blast (dustmasker)
bin/build_orbivirus_kraken_db.sh \
  --db /Users/todstuber/databases/k2_orbi_virus \
  --panel-dir /path/to/reference_guided_consensuses \   # optional; injects lab strains
  --test-r1 SAMPLE_R1.fastq.gz --test-r2 SAMPLE_R2.fastq.gz   # optional sanity check
```

The script downloads taxonomy, fetches the entire NCBI *Orbivirus* genus (via `datasets`, tagged
with species-level taxids) plus the RefSeq viral release as background context (via `curl`, tagged
at the "Viruses" level with orbivirus accessions excluded so orbivirus recall is never coarsened),
builds the index, generates the **Bracken** read-length distributions the pipeline needs
(`database250mers.kmer_distrib` and others), cleans up, links Krona's taxonomy, and reports the DB
size. (HTTPS is used throughout because NCBI rsync/FTP are blocked on many networks; `kraken2-build`
is invoked only for `--add-to-library`/`--build`. If the viral background download fails, the script
still builds an Orbivirus-only DB — the primary goal.) The
resulting directory is self-contained (`*.k2d` + `taxonomy/`) and can be `rsync`'d to the HPC.
Then run with the matching preset:

```bash
python bin/run_with_config.py --preset macos_dev_orbivirus_custom   # kraken_db -> k2_orbi_virus
```
On the HPC, set `orbivirus_custom`'s `kraken_db` in `internal/scomp_kraken_configs.yaml` to the
transferred path.

### BLAST Database
1. Create directory for BLAST databases:
```bash
cd ${HOME}
mkdir blast_databases
cd blast_databases
```

2. View available databases:
```bash
update_blastdb.pl --showall
```

3. Download and decompress database:
```bash
update_blastdb.pl --source aws ref_prok_rep_genomes
```

### Logo
The reports can include your organization's logo (`--logo`, or `logo:` in a preset), as a .png or .svg file.

**USDA reports:** use the official USDA signature lockup (USDA symbol + Department name) from the
[USDA Style Guide → Logo](https://www.usda.gov/about-usda/policies-and-links/digital/usda-style-guide/logo)
("Download logo and lockup files"). Per DR 1430-002 the USDA logo is the only identifier for USDA
agencies and programs, and lockups must not be recreated, so supply the Office of Communications
artwork as-is. The report follows the USDA visual standards:

- The logo is shown unmodified at the top left of a solid white Signature Iso-Bar (≥ 8% of the page
  height, min 0.75" in the PDF) with clear space around it; it is never stretched, recolored, boxed,
  shadowed or placed on a gradient.
- Colors are the official USDA Dark Blue `#002D72` (PMS 288) and Dark Green `#005440` (PMS 343), with
  U.S. Web Design System color tokens and usda.gov neutrals for everything else.
- Type is Public Sans / Source Sans Pro (usda.gov and USWDS), falling back to Helvetica/Arial. Install
  [Public Sans](https://github.com/uswds/public-sans) on the host for the exact usda.gov typeface.
- Accessibility (Section 508): all text has ≥ 4.5:1 contrast, charts mark SNP bases by shape as well as
  color, images carry alt text, and the PDF is tagged.

# Running the Pipeline

There are two main ways to run the pipeline:

## Option 1: Running with the Python Script

### Direct Method with Manual Parameters

You can run the pipeline directly using the `kraken_id_parse.py` script with all parameters specified:

```bash
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 *_R1*fastq.gz \
  -r2 *_R2*fastq.gz \
  --taxon "Mycobacterium tuberculosis complex" \
  --kraken_db ${HOME}/k2_standard_08gb \
  --blast_db ${HOME}/blast_databases/ref_prok_rep_genomes \
  --logo ${HOME}/logo.png
```

### Using Predefined Configurations

You can use predefined configurations with the `run_with_config.py` script:

```bash
# Activate the environment first
conda activate kraken_id_parse

# Run with a preset
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb

# Run with a custom config file
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb --config /path/to/custom/config.yaml
```

### Overriding Preset Parameters

You can use a preset while overriding specific parameters such as the taxon:

```bash
# Use the mtb preset but search for a different species
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb --override taxon="Mycobacterium bovis"

# Multiple overrides can be specified
python ${REPO_ROOT}/bin/run_with_config.py --preset orbivirus --override taxon="Bluetongue virus" --override logo="/path/to/custom/logo.png"
```

These same override options can be used with the SLURM script:

```bash
# Override taxon when running with SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb --override taxon="Mycobacterium bovis"

# Multiple overrides with SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus" --override logo="/path/to/custom/logo.png"
```

This is particularly useful when you want to use the database paths and other settings from a preset, but need to search for a different organism.

## Option 2: Using Predefined Configurations with SLURM

For systems with SLURM job scheduler, you can use the provided SLURM script with preset configurations:

1. First, define your presets in `internal/kraken_configs.yaml`:
   ```yaml
   presets:
     mtb:
       taxon: "Mycobacterium tuberculosis complex"
       kraken_db: "/path/to/kraken/database"
       blast_db: "/path/to/blast/database"
       logo: "/path/to/logo.png"
     
     orbivirus:
       taxon: "Orbivirus"
       kraken_db: "/path/to/kraken/database"
       blast_db: "/path/to/blast/database"
   ```

2. Run with a preset using the SLURM script:
   ```bash
   # Single run
   sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb
   
   # Process multiple directories
   currentdir=`pwd`
   for f in */; do 
     cd $f
     sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb
     cd $currentdir
   done
   ```

3. Run with a preset while overriding parameters:
   ```bash
   # Single run with custom taxon
   sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb --override taxon="Mycobacterium bovis"
   
   # Multiple parameters can be overridden
   sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb --override taxon="Mycobacterium bovis" --override logo="/path/to/custom/logo.png"
   
   # Process multiple directories with custom taxon
   currentdir=`pwd`
   for f in */; do 
     cd $f
     sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus"
     cd $currentdir
   done
   ```

The SLURM script will automatically:
- Load the appropriate conda environment
- Find your FASTQ files based on patterns in the config
- Run the analysis with the specified preset parameters and any overrides

# Testing

To Download files, if needed, download a [SRA ToolKit compiled binary package](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).  I have found using a compiled package from GitHub works better then installing from `conda`.

## Mycobacterium Tuberculosis Test

### Download Test Files
```bash
sra_number="SRR28623786"
wget -O "${sra_number}.fastq.gz" "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra_number}/${sra_number}"
fasterq-dump -S ${sra_number}.fastq.gz
```

### Prepare Files
```bash
rm ${sra_number}.fastq.gz
mv ${sra_number}.fastq.gz_1.fastq ${sra_number}_R1.fastq
mv ${sra_number}.fastq.gz_2.fastq ${sra_number}_R2.fastq
pigz *fastq
```

### Run Test
```bash
# Direct method
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 *_R1*fastq.gz \
  -r2 *_R2*fastq.gz \
  --taxon "Mycobacterium tuberculosis complex" \
  --kraken_db ${HOME}/k2_standard_08gb \
  --blast_db ${HOME}/blast_databases/ref_prok_rep_genomes \
  --logo ${HOME}/logo.png

# Or using preset configuration with Python script
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb

# With taxon override
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb --override taxon="Mycobacterium bovis"

# With multiple overrides
python ${REPO_ROOT}/bin/run_with_config.py --preset mtb --override taxon="Mycobacterium bovis" --override logo="${HOME}/custom_logo.png"

# Or using preset configuration with SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb

# With taxon override using SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb --override taxon="Mycobacterium bovis"

# With multiple overrides using SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm mtb --override taxon="Mycobacterium bovis" --override logo="${HOME}/custom_logo.png"
```

## Orbivirus Test

### Download Test Files
```bash
sra_number="SRR9598511"
wget -O "${sra_number}.fastq.gz" "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra_number}/${sra_number}"
fasterq-dump -S ${sra_number}.fastq.gz
```

### Prepare Files
```bash
rm ${sra_number}.fastq.gz
mv ${sra_number}.fastq.gz_1.fastq ${sra_number}_R1.fastq
mv ${sra_number}.fastq.gz_2.fastq ${sra_number}_R2.fastq
pigz *fastq
```

### Run Test
```bash
# Direct method
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 *_R1*fastq.gz \
  -r2 *_R2*fastq.gz \
  --taxon Orbivirus \
  --kraken_db ${HOME}/k2_standard_08gb \
  --blast_db ${HOME}/blast_databases/ref_prok_rep_genomes \
  --logo ${HOME}/logo.png

# Or using preset configuration with Python script
python ${REPO_ROOT}/bin/run_with_config.py --preset orbivirus

# With taxon override
python ${REPO_ROOT}/bin/run_with_config.py --preset orbivirus --override taxon="Bluetongue virus"

# With multiple overrides
python ${REPO_ROOT}/bin/run_with_config.py --preset orbivirus --override taxon="Bluetongue virus" --override logo="${HOME}/custom_logo.png"

# Or using preset configuration with SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm orbivirus

# With taxon override using SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus"

# With multiple overrides using SLURM
sbatch ${REPO_ROOT}/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus" --override logo="${HOME}/custom_logo.png"
```

## Orbivirus (BTV/EHD) databases and serotyping

> **BETA — the Orbivirus (BTV/EHD) analysis is untested.** The BTV/EHD split, segment
> assignments and BTV serotype need further testing, and their results cannot be trusted:
> confirm them by a validated method before acting on them. Every report, the stats workbook
> and the GUI say so for each Orbivirus run (`bin/orbivirus_beta.py` holds the wording).

The Orbivirus workflow is self-contained under [`database/`](database/README.md):

- **Custom Orbivirus Kraken DB** — `bin/build_orbivirus_kraken_db.sh` builds an uncapped,
  orbivirus-rich Kraken2 DB into `database/kraken_orbivirus/` (RefSeq viral + the entire NCBI
  Orbivirus genus + lab consensuses + the committed `database/source_sequences/BTV_sequence.fasta`
  and `EHD_sequence.fasta`). The heavy `.k2d` binaries are git-ignored and rebuilt on demand;
  only the source FASTAs + build script are committed. See
  [`database/kraken_orbivirus/README.md`](database/kraken_orbivirus/README.md). Use it via the
  `macos_dev_orbivirus_custom` preset (or the HPC `orbivirus_custom` stub after transfer).
- **BLAST references** — `bin/add_orbivirus_to_blast.sh` adds the BTV/EHD sequences to the
  external `nt_viruses` non-destructively (small supplement DB + union alias; `nt_viruses` is
  too large to commit). See [`database/blast_orbivirus/README.md`](database/blast_orbivirus/README.md).

### BTV serotyping (standalone or in-pipeline)

`bin/btv_serotyping.py` calls the BTV serotype by BLASTing VP2 (segment 2 — the serotype-defining
segment; VP5 is not used for the call) against the curated panel in `reference_sequences/` —
**no built database needed**. The panel
covers **serotypes BTV-1 … BTV-27** (VP2). It runs standalone on any consensus FASTA:

```bash
python bin/btv_serotyping.py -i SAMPLE_reference_guided.fasta -o ./sero_out
```

Grow the panel with validated samples via `bin/add_btv_reference.py`, or merge a
serotype-labeled VP2 FASTA with `bin/import_btv_serotype_panel.py`
(see `reference_sequences/BTV_panel_manifest.tsv` for current coverage).

**Mixed BTV + EHD samples:** each virus gets its own report section, and the consensus IDs carry
the species (`SAMPLE_BTV_segment2`, `SAMPLE_EHD_segment2`) so the reference-guided FASTA has no
duplicate IDs; single-virus samples keep `SAMPLE_segmentN`. Serotyping uses only the BTV segment 2,
and the Excel segment columns are reported per virus.

Segment assignment is regression-tested (Orbivirus, ISAV, Apicomplexa):
```bash
python bin/test_segment_profiles.py
```

## Customization for Different Systems

To port the **direct** and **SLURM** runners to another system:

1. Keep the repository structure the same.
2. Update `CONDA_PATH`, `CONDA_ENV`, and `REPO_ROOT` at the top of the SLURM script
   (`internal/kraken_id_parse.slurm` or `internal/scomp_kraken_id_parse_slurm.sh`).
3. Create system-specific presets (DB paths, taxon, logo) in
   `internal/kraken_configs.yaml` — these are shared by the direct and SLURM runners.
4. Use environment variables (e.g. `${HOME}`) for paths that change between systems.

**Packing samples onto one node under SLURM:** the pipeline auto-detects cores (leaving 2 free),
but you can pin threads so several jobs share a node — pass `--cpus-per-task` (honored
automatically) or add `--override threads=<n>`:
```bash
sbatch --cpus-per-task=8 internal/kraken_id_parse.slurm orbivirus
```

If needed, you can manually set REPO_ROOT in the SLURM script:
```bash
REPO_ROOT="${HOME}/git/gitlab/kraken_id_parse"
```

**Optional environment variables:**

| Variable | Default | Effect |
|---|---|---|
| `NCBI_EMAIL`, `NCBI_API_KEY` | unset | Identify NCBI downloads (raises the E-utilities rate limit, fewer refused requests) |
| `KIP_THREADS` | auto (cores − 2) | Worker threads (same as `--threads`) |
| `BLAST_TASK` | `dc-megablast` | BLAST task for queries under 50 kb (viral/segmented assemblies and consensus) |
| `BLAST_TASK_LARGE` | `megablast` (or `BLAST_TASK` if that is set) | BLAST task for larger queries (bacterial assemblies/consensus, large DNA viruses) |
| `BLAST_MIN_CONTIG_LEN` | `300` | Shorter de novo contigs are not BLASTed |
| `BLAST_TIMEOUT_SECONDS` | `3600` | Limit for each BLAST command |