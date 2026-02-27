#!/bin/bash

# Cross-platform conda environment installation script for kraken_id_parse
# This script handles the platform-specific setup automatically

set -e  # Exit on any error

echo "🧬 Installing Kraken ID Parse Environment"
echo "=========================================="

# Detect operating system
OS="$(uname -s)"
case "${OS}" in
    Linux*)     MACHINE=Linux;;
    Darwin*)    MACHINE=Mac;;
    CYGWIN*)    MACHINE=Cygwin;;
    MINGW*)     MACHINE=MinGw;;
    *)          MACHINE="UNKNOWN:${OS}"
esac

echo "Detected OS: ${MACHINE}"

# Set macOS-specific environment variable
if [ "${MACHINE}" = "Mac" ]; then
    echo "Setting macOS compatibility override..."
    export CONDA_OVERRIDE_OSX=12.0
fi

# Check if environment already exists
if conda env list | grep -q "kraken_id_parse"; then
    echo "⚠️  Environment 'kraken_id_parse' already exists."
    echo "To reinstall, first remove it: conda env remove -n kraken_id_parse -y"
    exit 1
fi

# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"

# Check if environment file exists
if [ ! -f "${ENV_FILE}" ]; then
    echo "❌ Environment file not found: ${ENV_FILE}"
    exit 1
fi

# Try conda first, fallback to mamba
echo "📦 Creating conda environment..."
if command -v conda &> /dev/null; then
    if conda env create -f "${ENV_FILE}"; then
        echo "✅ Environment created successfully with conda"
    else
        echo "❌ Failed to create environment with conda"

        # Try mamba as fallback
        if command -v mamba &> /dev/null; then
            echo "🔄 Trying with mamba as fallback..."
            if mamba env create -f "${ENV_FILE}"; then
                echo "✅ Environment created successfully with mamba"
            else
                echo "❌ Failed to create environment with both conda and mamba"
                exit 1
            fi
        else
            echo "❌ Neither conda nor mamba succeeded. Please check the installation logs above."
            exit 1
        fi
    fi
else
    echo "❌ conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

# Configure matplotlib for macOS
if [ "${MACHINE}" = "Mac" ]; then
    echo "🎨 Configuring matplotlib for macOS..."
    mkdir -p ~/.matplotlib
    echo "backend: Agg" > ~/.matplotlib/matplotlibrc
    echo "✅ matplotlib configured for headless operation"
fi

# Optional LaTeX installation
echo ""
echo "📄 LaTeX Installation (Optional)"
echo "The pipeline requires pdflatex for PDF reports."
echo ""
read -p "Install LaTeX via conda? (y/N): " install_latex
case "$install_latex" in
    [Yy]* )
        echo "📦 Installing LaTeX packages..."
        if conda install -n kraken_id_parse -c conda-forge texlive-core pdflatex -y; then
            echo "✅ LaTeX installed successfully"
        else
            echo "⚠️ LaTeX installation failed. You can install it manually later:"
            echo "  conda activate kraken_id_parse"
            echo "  conda install -c conda-forge texlive-core pdflatex"
        fi
        ;;
    * )
        echo "⏭️  Skipping LaTeX installation"
        echo "You can install LaTeX later with:"
        echo "  conda activate kraken_id_parse"
        echo "  conda install -c conda-forge texlive-core pdflatex"
        echo "Or use system packages (brew, apt, etc.)"
        ;;
esac

echo ""
echo "🎉 Installation complete!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate kraken_id_parse"
echo ""
echo "To verify the installation:"
echo "  conda activate kraken_id_parse"
echo "  python -c \"import numpy; import pandas; import biopython; print('All packages working!')\""
if command -v pdflatex &> /dev/null; then
    echo "  pdflatex --version  # Verify LaTeX installation"
fi