# Conda/Mamba Setup Instructions

## Prerequisites

Note: We use `mamba` instead of `conda` as the Bracken tool installation has been more reliable with Mamba. However, one may be successful with `conda` also.

### Installing Mamba
If you don't have Mamba installed, install it using Conda:
```bash
conda install -n base -c conda-forge mamba
```

## Installation Steps

### 1. Create Environment
Create environment using the provided `environment.yml` file:
```bash
mamba env create -f environment.yml
```

### 2. Configure Matplotlib (macOS Only)
Create matplotlib configuration for macOS systems:
```bash
mkdir -p ~/.matplotlib
echo "backend: Agg" > ~/.matplotlib/matplotlibrc
```

### 3. Activate Environment
```bash
mamba activate kraken_id_parse
```

### 4. Verify Installation
Test the matplotlib backend configuration:
```bash
python -c "import matplotlib; print(matplotlib.get_backend())"
```

## Environment File Location
The `environment.yml` file can be found in the `conda_setup` directory of this repository.