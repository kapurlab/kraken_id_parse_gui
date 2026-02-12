#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --job-name="parse_id"
#SBATCH --partition=prod-compute-mem
#SBATCH --nice
#SBATCH --export=NONE

# Usage:
# --partition=prod-compute,prod-compute-mem
# Single Run:
# sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm lepto
# sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus                             # use preset name
# sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus"  # override taxon
# sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm isav
#
# Multiple runs (one per directory):
#
# packagefastqs.sh; currentdir=`pwd`; for f in */; do echo $f; cd $f; sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm lepto; cd $currentdir; done
# packagefastqs.sh; currentdir=`pwd`; for f in */; do echo $f; cd $f; sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus; cd $currentdir; done
# packagefastqs.sh; currentdir=`pwd`; for f in */; do echo $f; cd $f; sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus --override taxon="Bluetongue virus"; cd $currentdir; done
# packagefastqs.sh; currentdir=`pwd`; for f in */; do echo $f; cd $f; sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm isav; cd $currentdir; done

# ===== CONFIGURATION (modify as needed) =====
# Initialize conda
echo "Initializing miniconda..."
. /project/shared/miniconda3/etc/profile.d/conda.sh
echo "Activating kraken_parse_and_id"
conda activate kraken_parse_and_id
# ===== CONFIGURATION (modify as needed) =====

# Path to conda environment
CONDA_PATH="/project/shared/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="kraken_parse_and_id"

# Path to repository root - THE MOST IMPORTANT SETTING
# IMPORTANT: When running as a Slurm job, automatic path detection often fails
# Either set this manually to your absolute repo path:
REPO_ROOT="${HOME}/git/gitlab/kraken_id_parse/"

# Or uncomment and modify this line if you need a different path:
# REPO_ROOT="/path/to/your/kraken_id_parse"

# ===== END CONFIGURATION =====

# Check if a preset was provided
if [ -z "$1" ]; then
    echo "Error: No preset specified"
    echo "Usage: sbatch $0 <preset_name> [--override key=value [--override key=value ...]]"
    echo "Available presets:"
    grep "^  [a-z]" "${REPO_ROOT}/internal/kraken_configs.yaml" | cut -d':' -f1 | sed 's/^  /    /'
    exit 1
fi

# Extract the preset name and shift arguments
PRESET=$1
shift  # Remove preset from arguments

# Initialize conda if available
if [ -f "$CONDA_PATH" ]; then
    echo "Initializing conda from $CONDA_PATH..."
    . "$CONDA_PATH"
    echo "Activating $CONDA_ENV environment"
    conda activate $CONDA_ENV
else
    echo "Warning: Conda initialization file not found at $CONDA_PATH"
    echo "Proceeding without conda activation"
fi

# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Working Directory: $PWD"
echo "Repository Root: $REPO_ROOT"
echo "Preset: $PRESET"
echo "Additional arguments: $@"
date

# Run the Python script with the configuration
echo "Running: python ${REPO_ROOT}/bin/run_with_config.py --preset $PRESET $@"

# Check if the script exists
if [ ! -f "${REPO_ROOT}/bin/run_with_config.py" ]; then
    echo "Error: Cannot find run_with_config.py at ${REPO_ROOT}/bin/run_with_config.py"
    echo "Please check that REPO_ROOT is set correctly in this script."
    echo "Current REPO_ROOT: ${REPO_ROOT}"
    exit 1
fi

# Run with the preset and any additional arguments
python "${REPO_ROOT}/bin/run_with_config.py" --preset "$PRESET" "$@" 

# Print completion
echo "Job finished"
date