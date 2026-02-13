# Kraken ID and Parse

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

Upon completion, Excel and PDF files are created to provide summaries and reporting.

## Workflow
1. Kraken2 classification and Bracken abundance estimation
2. Reads parsed on taxon
3. **Platform auto-detection** — reads > 701 bp → ONT mode; otherwise Illumina
4. Parsed reads assembled with SPAdes (Illumina) or Flye (ONT)
5. Assembly identified with BLAST
6. Download FASTAs from BLAST findings (with **local reference cache**)
7. Coverage graph on downloaded FASTAs (BWA for Illumina, minimap2 for ONT)
8. If available continue with taxon specific workflow:
   - If Orbivirus, check BLAST coverage stats as "bluetongue virus" or "epizootic hemorrhagic disease"
     - Order by segment and make coverage graph
9. Reference guided assembly using the original FASTQ files
10. BLAST reference guided assembly
11. Alignment using top BLAST results
12. Merged VCF to the top BLAST results as the final consensus
13. Make final coverage graph
14. Generate PDF report with platform-aware methodology appendix

## Output Files

After a successful run, the output directory contains:

| File | Description |
|------|-------------|
| `*_report.pdf` | PDF report — FASTQ quality, Kraken/Bracken pie chart, classification tables, assembly stats, BLAST identification, coverage graphs, methodology appendix |
| `*_stats.xlsx` | Excel spreadsheet with summary statistics |
| `*_denovo.fasta` | *De novo* SPAdes assembly of extracted reads |
| `*_reference_guided.fasta` | Consensus sequence from reference-guided assembly |
| `*_blast_summary.txt` | BLAST results summary (de novo assembly) |
| `consensus_blast_summary.txt` | BLAST results summary (consensus assembly) |
| `CAUTION_SITES.xlsx` | Ambiguous/heterozygous sites flagged during consensus calling |
| `kraken/` | Kraken2 report, Bracken species estimates, Krona interactive HTML visualization |

Intermediate files (section banners, coverage graph PDFs, LaTeX source, pie chart PNG, raw BLAST output) are automatically cleaned up after the PDF report is compiled.

Use `--keep-extracted-reads` to preserve the taxon-filtered FASTQ files (default: removed to save space).

# Installation

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/kapurlab/kraken_id_parse_gui.git
cd kraken_id_parse_gui

# 2. Create conda environment
mamba env create -f conda_setup/environment.yml
mamba activate kraken_id_parse

# 3. Install LaTeX (required for PDF report generation — see below)

# 4. Download Kraken and BLAST databases (see Prerequisites below)

# 5. Run
export REPO_ROOT="$(pwd)"
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz \
  --taxon Orbivirus \
  --kraken_db /path/to/kraken_db \
  --blast_db /path/to/blast_db
```

## Conda Environment

Follow detailed instructions at [conda_setup/conda_setup.md](./conda_setup/conda_setup.md)

```bash
mamba env create -f conda_setup/environment.yml
mamba activate kraken_id_parse
```

## LaTeX (pdflatex)

The pipeline requires `pdflatex` to compile the PDF report. This is **not** included in the conda environment and must be installed separately.

**macOS:**
```bash
# Full MacTeX (recommended, ~3.5 GB)
brew install --cask mactex-no-gui

# Or TinyTeX (lighter, ~150 MB — then install required packages)
curl -sL "https://yihui.org/tinytex/install-bin-unix.sh" | sh
tlmgr install adjustbox collectbox fancyhdr grfext float xcolor helvet psnfss
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

**Verify installation:**
```bash
which pdflatex
pdflatex --version
```

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
update_blastdb.pl --decompress ref_prok_rep_genomes
```

### Logo
The reports can include your organization's logo. Provide an image as a .png file for enhanced report presentation.

# Running the Pipeline

## Command-Line Options

```
kraken_id_parse.py [-h] [-r1 FASTQ_R1] [-r2 FASTQ_R2] [-l LOGO]
                   -t TAXON -k KRAKEN_DB [-b BLAST_DB]
                   [--database-root DIR] [--reference-cache DIR]
                   [--platform {illumina,ont,auto}]
                   [-s SPECIFIC] [--keep-extracted-reads] [-d] [-v]

Required:
  -r1, --read1          R1 FASTQ file (or single read)
  -t,  --taxon          Target taxon name (e.g., "Orbivirus", "Flaviviridae")
  -k,  --kraken_db      Path to Kraken2 database directory

Optional:
  -r2, --read2          R2 FASTQ file (paired-end)
  -b,  --blast_db       BLAST database path (overrides auto-resolution)
  --database-root       Root dir containing BLAST databases (enables auto-resolution)
  --reference-cache     Directory to cache downloaded NCBI reference FASTAs
  --platform            Sequencing platform: illumina, ont, or auto (default: auto)
  -l,  --logo           Organization logo PNG for report header
  -s,  --specific       Custom taxon-specific function name
  --keep-extracted-reads  Keep taxon-filtered FASTQ files (default: remove)
  -d,  --debug          Keep all intermediate/temp files
  -v,  --version        Show version
```

## Automatic BLAST Database Resolution

When `--database-root` is provided (or set in `~/.kraken_id_parse.yaml`), the pipeline automatically selects the appropriate BLAST database based on the taxon name:

| Taxon Category | BLAST Database |
|----------------|----------------|
| Viral families & genera (Orbivirus, Flaviviridae, etc.) | `nt_viruses` |
| Bacterial (Mycobacterium, Brucellaceae, Leptospirales) | `ref_prok_rep_genomes` |
| Protozoan/Eukaryotic (Apicomplexa) | `nt` |
| Unknown taxa with "virus"/"viridae" in name | `nt_viruses` |
| Fallback | `nt` |

### Setup

Create `~/.kraken_id_parse.yaml`:
```yaml
database_root: "/path/to/blast/databases"
reference_cache: "/path/to/reference_cache"
```

Then simply run:
```bash
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz \
  --taxon Orbivirus \
  --kraken_db /path/to/kraken_db
# BLAST database auto-resolved to: /path/to/blast/databases/nt_viruses
```

Use `-b` to override auto-resolution for a specific run.

## Platform Auto-Detection (Illumina / ONT)

The pipeline automatically detects the sequencing platform by sampling the first 100 reads from the R1 FASTQ file:

- **Reads > 701 bp** → ONT mode (minimap2 + Flye)
- **Reads ≤ 701 bp** → Illumina mode (BWA + SPAdes)

| Component | Illumina | ONT |
|-----------|----------|-----|
| Assembler | SPAdes | Flye (`--nano-raw`) |
| Aligner | BWA-MEM | minimap2 (`-x map-ont`) |
| VCF QUAL | As-is | +100 normalization |

Override with `--platform illumina` or `--platform ont` to skip auto-detection.

## Local Reference Cache

When `--reference-cache` is provided (or set in `~/.kraken_id_parse.yaml`), downloaded NCBI reference FASTAs are cached locally. Subsequent runs reuse cached files instead of re-downloading, which significantly speeds up repeated analyses on the same taxa.

```bash
# First run — downloads and caches references
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 reads_R1.fastq.gz -r2 reads_R2.fastq.gz \
  --taxon Orbivirus --kraken_db /path/to/db \
  --reference-cache /path/to/cache

# Second run — prints "Cache hit" and skips downloads
${REPO_ROOT}/bin/kraken_id_parse.py \
  -r1 reads_R1.fastq.gz -r2 reads_R2.fastq.gz \
  --taxon Orbivirus --kraken_db /path/to/db \
  --reference-cache /path/to/cache
```

Or set it permanently in `~/.kraken_id_parse.yaml`:
```yaml
reference_cache: "/path/to/reference_cache"
```

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

## Customization for Different Systems

To customize for a different system:

1. Keep the repository structure the same
2. Update the CONDA_PATH and CONDA_ENV in the SLURM script
3. Create system-specific presets in the kraken_configs.yaml file
4. Use environment variables for paths that change between systems

If needed, you can manually set REPO_ROOT in the SLURM script:
```bash
REPO_ROOT="${HOME}/git/gitlab/kraken_id_parse"
```

# Changelog

## v0.2.0

### ONT / MinION Long-Read Support
- **Platform auto-detection** — reads > 701 bp automatically trigger ONT mode
- **Flye assembler** for ONT de novo assembly (`--nano-raw`)
- **minimap2 aligner** for ONT read alignment (`-x map-ont`)
- **VCF QUAL normalization** — +100 to all QUAL scores for ONT reads (compensates for lower nanopore base quality)
- **`--platform` CLI option** — override auto-detection with `illumina`, `ont`, or `auto`
- **Platform-aware PDF report** — methodology appendix dynamically references correct tools

### Local Reference Cache
- **`--reference-cache` CLI option** — cache downloaded NCBI reference FASTAs locally
- **`~/.kraken_id_parse.yaml` support** — set `reference_cache` for persistent configuration
- Cache hit/miss logging — "Cache hit: {accession}" when a cached reference is reused
- Zero-byte or missing cached files treated as cache miss

### BLAST Database Auto-Resolution
- **`--database-root` CLI option** — automatically selects the appropriate BLAST database based on taxon name
- **`~/.kraken_id_parse.yaml` support** — set `database_root` for persistent configuration
- Taxon-to-database mapping: viral → `nt_viruses`, bacterial → `ref_prok_rep_genomes`, fallback → `nt`

### PDF Report Enhancements
- Methodology appendix with detailed pipeline description
- FASTQ quality statistics section
- Kraken/Bracken pie chart visualization
- Automatic cleanup of intermediate LaTeX and image files
- `--keep-extracted-reads` option to preserve taxon-filtered FASTQs

### Bug Fixes
- Fixed `coverage_graph.py` crash when R2 is None (single-end reads)
- Fixed LaTeX report rendering — all sections now compile correctly
- Skip missing logo path to avoid LaTeX failures

### Infrastructure
- Preset configuration system with `run_with_config.py` and `--override` support
- SLURM job script with preset integration
- Preflight tool checks in GUI
- `minimap2` and `flye` added to conda environment

## v0.1.0
- Initial release
- Kraken2 classification and Bracken abundance estimation
- SPAdes de novo assembly
- BLAST identification
- BWA coverage graphs
- Orbivirus-specific segment ordering
- Reference-guided assembly with consensus VCF
- Excel and PDF report generation