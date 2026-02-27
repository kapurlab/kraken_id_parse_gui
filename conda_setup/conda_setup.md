# Conda Environment Setup Instructions

## Cross-Platform Installation

The environment creation has been updated to work reliably across all platforms (macOS, Linux, Windows).

### Automated Installation (Recommended)

Use the provided installation script that handles platform detection and setup automatically:

```bash
bash conda_setup/install_environment.sh
conda activate kraken_id_parse
```

This script will:
- Detect your operating system
- Set appropriate environment variables for macOS
- Try conda first, fallback to mamba if needed
- Configure matplotlib backend for macOS
- Provide clear success/failure feedback

### Manual Installation

### Method 1: Using Conda (Recommended)

```bash
# For macOS with recent versions, set OSX override to avoid compatibility issues
export CONDA_OVERRIDE_OSX=12.0  # macOS only

# Create the environment
conda env create -f conda_setup/environment.yml

# Activate the environment
conda activate kraken_id_parse
```

### Method 2: Using Mamba (Alternative)

If you prefer mamba or encounter issues with conda:

```bash
# Install mamba if not already available
conda install -n base -c conda-forge mamba

# For macOS with recent versions, set OSX override
export CONDA_OVERRIDE_OSX=12.0  # macOS only

# Create the environment with mamba
mamba env create -f conda_setup/environment.yml

# Activate the environment
conda activate kraken_id_parse
```

## Platform-Specific Notes

### macOS
- **macOS 11.0+**: The `CONDA_OVERRIDE_OSX=12.0` environment variable is required to resolve package dependencies
- **matplotlib backend**: Configure for headless operation:
  ```bash
  mkdir -p ~/.matplotlib
  echo "backend: Agg" > ~/.matplotlib/matplotlibrc
  ```

### Linux
- No special configuration required
- Standard conda/mamba installation should work without issues

### Windows
- **Recommended**: Use Windows Subsystem for Linux (WSL) + conda for best compatibility
- **Alternative**: Native conda on Windows (some bioinformatics tools may have limitations)
- **Automated installer**: Works with WSL, Cygwin, or MINGW environments
- **Note**: Many bioinformatics tools are primarily designed for Unix-like systems

## Troubleshooting

### If Environment Creation Fails:
1. **Update conda**: `conda update -n base -c conda-forge conda`
2. **Clear cache**: `conda clean --all`
3. **Try conda instead of mamba** or vice versa
4. **For macOS**: Ensure you set `export CONDA_OVERRIDE_OSX=12.0`

### Package Compatibility Issues:
The environment file has been updated to use modern, compatible package versions:
- Python 3.10 (instead of 3.8)
- Updated numpy, pandas, and other core packages
- `pyvcf3` instead of the deprecated `pyvcf`
- Compatible bracken version (1.0.0)

### Verification
After successful installation, verify key packages:
```bash
conda activate kraken_id_parse
python -c "import numpy; import pandas; import biopython; print('All packages imported successfully!')"
```

### Optional: Install LaTeX via Conda
To avoid using system package managers like brew or apt, you can install LaTeX through conda:
```bash
conda activate kraken_id_parse
conda install -c conda-forge texlive-core pdflatex
```

This provides a complete LaTeX installation for PDF report generation without requiring brew or system packages.

## Environment File Location
The updated `environment.yml` file is in the `conda_setup` directory of this repository and has been tested on macOS, Linux, and Windows systems.