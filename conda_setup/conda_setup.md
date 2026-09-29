# Conda/Mamba Setup Instructions

> In this repository `environment.yml` is the GUI's single environment, `kraken_id_parse`
> (the pipeline, Bracken where bioconda builds it, and the web backend) — the one
> `bdtools install kraken_id_parse_gui` builds. `environment.bracken.yml` is the fallback for
> osx-64, where Bracken 2.x cannot share it. The rest of this page is the pipeline's own setup
> guide.

## Prerequisites

We recommend `mamba` over `conda` — it resolves this environment (Kraken2, SPAdes, WeasyPrint,
Plotly, etc.) considerably faster and more reliably.

### Installing Mamba
If you don't have Mamba installed, install it using Conda:
```bash
conda install -n base -c conda-forge mamba
```

## Two environments

The pipeline uses **two** conda environments, both defined in this directory:

| Environment | File | Contents |
|---|---|---|
| `kraken_id_parse` (main) | `environment.yml` | Python 3.11, Kraken2, SPAdes, BLAST, BWA, samtools, freebayes, and the HTML/PDF report stack (Jinja2, WeasyPrint, Plotly). This is the environment you activate to run the pipeline. |
| `bracken` | `environment.bracken.yml` | Bracken only. |

**Why Bracken is split out:** on macOS (`osx-64`), the only available Bracken 2.x builds link an
old `zlib` (<1.3.0) that cannot coexist in the same environment as the modern scientific stack
(numpy/WeasyPrint pull `zlib` 1.3.2). Isolating Bracken into its own environment avoids that
conflict. You never call Bracken directly — `bin/kraken.py` auto-detects
`<conda base>/envs/bracken/bin/bracken` next to the main environment (override with the
`BRACKEN_BIN` environment variable if your Bracken env lives somewhere else).

## Installation steps

### 1. Create both environments
```bash
cd conda_setup
mamba env create -f environment.yml
mamba env create -f environment.bracken.yml
```

#### macOS note
Bioconda does not publish native `osx-arm64` builds for several packages this pipeline needs
(Kraken2, Bracken, SPAdes). On Apple Silicon, create both environments under the Intel (`osx-64`)
subdirectory via Rosetta instead:
```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml
CONDA_SUBDIR=osx-64 mamba env create -f environment.bracken.yml
```

### 2. Configure Matplotlib (macOS only)
```bash
mkdir -p ~/.matplotlib
echo "backend: Agg" > ~/.matplotlib/matplotlibrc
```

### 3. Activate the main environment
```bash
mamba activate kraken_id_parse
```
(The `bracken` environment is never activated directly — the pipeline finds it automatically.)

### 4. Verify the installation
```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
python -c "import weasyprint, plotly, jinja2; print('report stack OK')"
kraken2 --version
mamba run -n bracken bracken --help >/dev/null && echo "bracken OK"
```

## Building a custom Kraken database
`bin/build_orbivirus_kraken_db.sh` builds an Orbivirus-enriched Kraken2 DB (see the main
[README](../README.md), section "Build a custom Orbivirus-enriched Kraken DB"). It uses
`datasets`/`dataformat` from the `ncbi-datasets-cli` package, `curl`, `seqkit`, and `dustmasker`
(from `blast`) — all included in `environment.yml`. All NCBI downloads go over **HTTPS** (NCBI
rsync on port 873 and FTP are blocked on many networks), and `kraken2-build` is used only for
`--add-to-library`/`--build`, so kraken2's bundled rsync/FTP downloaders are never invoked.

## Environment file locations
Both `environment.yml` and `environment.bracken.yml` live in this (`conda_setup/`) directory.
