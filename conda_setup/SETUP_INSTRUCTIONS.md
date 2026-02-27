# Kraken ID Parse GUI - Environment Setup

## Quick Start

### Prerequisites
- Conda or Miniconda installed
- macOS, Linux, or Windows

### Basic Environment Creation

```bash
# Navigate to the project directory
cd /path/to/kraken_id_parse_gui

# Create the environment (use conda, not mamba on macOS due to codesigning issues)
CONDA_OVERRIDE_OSX=11.0 conda env create -f conda_setup/environment.yml

# Activate the environment
conda activate kraken_id_parse

# Verify installation
python -c "import numpy, pandas, matplotlib, Bio, allel, vcf; print('Environment ready!')"

# Set up Krona taxonomy (required for taxonomic visualizations)
rm -rf ${HOME}/miniconda3/envs/kraken_id_parse/opt/krona/taxonomy
ln -s ${HOME}/k2_standard_08gb ${HOME}/miniconda3/envs/kraken_id_parse/opt/krona/taxonomy
```

**Note**: The krona taxonomy setup assumes you have a kraken database at `${HOME}/k2_standard_08gb`. Adjust the path as needed for your kraken database location.

## macOS Specific Notes

- **Use `conda` instead of `mamba`** to avoid codesigning issues
- Set `CONDA_OVERRIDE_OSX=11.0` environment variable if you encounter macOS version detection issues
- If you get permission errors, you may need to run with appropriate permissions

## Extended Bioinformatics Tools

After the basic environment is working, you can add additional tools:

```bash
conda activate kraken_id_parse

# Install additional bioinformatics tools
CONDA_OVERRIDE_OSX=11.0 conda install -c conda-forge -c bioconda \
    blast bwa spades seqkit krona bracken picard vcflib freebayes parallel pigz
```

## Testing the Installation

Test that key components work:

```bash
conda activate kraken_id_parse

# Test Python packages
python -c "import numpy, pandas, matplotlib, Bio, allel, vcf; print('✅ Python packages OK')"

# Test bioinformatics tools
kraken2 --version
samtools --version
```

## Troubleshooting

### Common Issues

1. **Codesigning errors on macOS**: Use `conda` instead of `mamba`
2. **OSX version not detected**: Set `CONDA_OVERRIDE_OSX=11.0`
3. **Package not found**: Some packages may not be available for your platform
4. **Environment creation fails**: Try creating with minimal packages first, then add more

### Minimal Environment Creation

If the full environment fails, create a minimal one first:

```bash
# Create minimal environment
conda create -n kraken_id_parse python=3.10 pip numpy pandas matplotlib biopython pytest

# Activate and add tools incrementally
conda activate kraken_id_parse
conda install -c bioconda kraken2 krakentools samtools
pip install PyVCF3
```

### Environment Management

```bash
# List environments
conda env list

# Remove environment if needed
conda env remove -n kraken_id_parse

# Update environment
conda env update -f conda_setup/environment.yml

# Export current environment
conda env export -n kraken_id_parse > my_environment.yml
```

## Package Information

### Core Analysis Packages
- **numpy, pandas**: Data manipulation and analysis
- **matplotlib, seaborn**: Data visualization
- **scipy, scikit-learn**: Scientific computing and machine learning
- **biopython**: Bioinformatics toolkit

### Bioinformatics Tools
- **kraken2**: Taxonomic classification of metagenomic sequences
- **krakentools**: Utilities for working with kraken2 output
- **samtools**: SAM/BAM file manipulation
- **scikit-allel**: Genetic variation analysis
- **PyVCF3**: VCF file handling

### Development Tools
- **pytest**: Testing framework
- **h5py**: HDF5 file support
- **openpyxl**: Excel file handling

## Support

If you encounter issues:
1. Check the troubleshooting section above
2. Verify your conda/mamba installation
3. Try creating a minimal environment first
4. Check package availability for your platform