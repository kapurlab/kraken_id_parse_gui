#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --job-name="parse_id"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --nice
#SBATCH --export=NONE

# Single Run
# sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus # preset name
#
# Many
# packagefastqs.sh
# currentdir=`pwd`; for f in *; do echo $f; cd ./$f; sbatch ${HOME}/git/gitlab/kraken_id_parse/internal/kraken_id_parse.slurm orbivirus; cd $currentdir; done


# Set path to internal directory
INTERNAL_DIR="${HOME}/git/gitlab/kraken_id_parse/internal"

# Check if a preset was provided
if [ -z "$1" ]; then
    echo "Error: No preset specified"
    echo "Usage: sbatch $0 <preset_name>"
    echo "Available presets:"
    grep "^  [a-z]" "${INTERNAL_DIR}/kraken_configs.yaml" | cut -d':' -f1 | sed 's/^  /    /'
    exit 1
fi

PRESET=$1

# Initialize conda
echo "Initializing miniconda..."
. /project/shared/miniconda3/etc/profile.d/conda.sh
echo "Activating kraken_parse_and_id"
conda activate kraken_parse_and_id

# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Working Directory: $PWD"
date

# Run the Python script with the configuration
python ${INTERNAL_DIR}/run_with_config.py --preset $PRESET

# Print completion
echo "Job finished"
date