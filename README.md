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
1. Kraken
2. Reads parsed on taxon
3. Parsed reads assembled
4. Assembly identified with BLAST
5. Download FASTAs from BLAST findings
6. Coverage graph on downloaded FASTAs
7. If available continue with taxon specific workflow:
   - If Orbivirus, check BLAST coverage stats as "bluetongue virus" or "epizootic hemorrhagic disease"
     - Order by segment and make coverage graph
8. Reference guided assembly using the original FASTQ files
9. BLAST reference guided assembly
10. Alignment using top BLAST results
11. Merged VCF to the top BLAST results as the final consensus
12. Make final coverage graph

# Installation

Follow instructions at [conda_setup/SETUP_INSTRUCTIONS.md](./conda_setup/SETUP_INSTRUCTIONS.md)

Quick setup:
```bash
# Create environment (use conda, not mamba on macOS)
CONDA_OVERRIDE_OSX=11.0 conda env create -f conda_setup/environment.yml

# Activate environment
conda activate kraken_id_parse
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
update_blastdb.pl --source aws ref_prok_rep_genomes
```

### Logo
The reports can include your organization's logo. Provide an image as a .png file for enhanced report presentation.

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