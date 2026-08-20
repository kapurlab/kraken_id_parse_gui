#!/usr/bin/env python3

"""
Core setup and utility classes for bioinformatics pipeline tools.
Provides base classes for file handling, reporting, and formatting.
"""

import os
import shutil
import re
from typing import Optional, Dict, List, Union, Tuple
import pandas as pd
import multiprocessing
from datetime import datetime
from pathlib import Path

def safe_move(src, dst):
    """shutil.move that silently overwrites an existing destination file or directory."""
    src_path = Path(src)
    dst_path = Path(dst)
    actual_dst = dst_path / src_path.name if dst_path.is_dir() else dst_path
    if actual_dst.exists():
        shutil.rmtree(actual_dst) if actual_dst.is_dir() else actual_dst.unlink()
    shutil.move(str(src), str(dst))


try:
    import svgwrite
    from cairosvg import svg2png
    from PIL import Image
    import numpy as np
    import colorsys
    HAS_SVG_SUPPORT = True
except ImportError:
    HAS_SVG_SUPPORT = False


def apply_mpl_style(style_candidates: Optional[List[str]] = None):
    """Apply a Matplotlib style with fallback across Matplotlib versions."""
    import matplotlib

    if not os.environ.get("MPLBACKEND"):
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    candidates = style_candidates or [
        "seaborn-v0_8-colorblind",
        "seaborn-colorblind",
        "seaborn-v0_8",
        "seaborn",
        "ggplot",
    ]
    for style in candidates:
        try:
            plt.style.use(style)
            break
        except OSError:
            continue
    return plt


class bcolors:
    """ANSI color codes for terminal output formatting"""
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

class Setup:
    def __init__(self, 
                 FASTA: Optional[str] = None,
                 FASTQ_R1: Optional[str] = None,
                 FASTQ_R2: Optional[str] = None,
                 debug: bool = False):
        """
        Initialize Setup class with input files.
        """
        self.cwd = os.getcwd()
        self.debug = debug
        self.paired = bool(FASTQ_R2)
        
        # Initialize storage for FASTQ files
        self.FASTQ_list = []
        self.FASTQ_dict = {}
        self.FASTQ_R1 = None
        self.FASTQ_R2 = None
        self.FASTA = None

        # Set database paths
        self._set_database_paths()  # Add this call here

        # Process FASTQ files if provided
        if FASTQ_R1:
            self.FASTQ_R1 = self._setup_fastq(FASTQ_R1, 'R1')
            if FASTQ_R2:
                self.FASTQ_R2 = self._setup_fastq(FASTQ_R2, 'R2')
        
        # Process FASTA if provided
        if FASTA:
            self.FASTA = self._setup_fasta(FASTA)
        
        # Set up sample name
        if FASTQ_R1:
            self.sample_name = re.sub('[_.].*', '', os.path.basename(FASTQ_R1))
        elif FASTA:
            self.sample_name = re.sub('[_.].*', '', os.path.basename(FASTA))
        
        # Set up additional attributes
        self.startTime = datetime.now()
        self.cpus = max(1, multiprocessing.cpu_count() - 2)
        # A dashboard session (or Slurm allocation) declares a core budget for
        # every tool it launches. Without this cap a single kraken2/bwa/blastn
        # run asks for every core on a shared box, and a room of concurrent
        # users turns into thread-thrash. Absent the variables, standalone
        # behavior is unchanged.
        for _var in ("BDTOOLS_SESSION_CORES", "SLURM_CPUS_PER_TASK"):
            _val = os.environ.get(_var, "").strip()
            if _val.isdigit() and int(_val) > 0:
                self.cpus = max(1, min(self.cpus, int(_val)))
                break
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        
        if self.debug:
            print(f"Debug - Setup initialized with:")
            print(f"  FASTQ_R1: {self.FASTQ_R1}")
            print(f"  FASTQ_R2: {self.FASTQ_R2}")
            print(f"  FASTA: {self.FASTA}")
            print(f"  FASTQ_list: {self.FASTQ_list}")
            print(f"  FASTQ_dict: {self.FASTQ_dict}")

    def _setup_fastq(self, fastq_path: str, read_type: str) -> str:
        """Process and setup FASTQ file"""
        if not fastq_path:
            return None
            
        try:
            local_path = os.path.abspath(fastq_path)
            if not os.path.exists(local_path):
                local_path = os.path.join(self.cwd, os.path.basename(fastq_path))
                if not os.path.exists(local_path):
                    shutil.copy2(fastq_path, local_path)
            
            # Add to FASTQ_list and FASTQ_dict
            self.FASTQ_list.append(local_path)
            self.FASTQ_dict[f'FASTQ_{read_type}'] = local_path
            
            return local_path
            
        except (TypeError, OSError) as e:
            if self.debug:
                print(f"Error setting up FASTQ {read_type}: {e}")
            return None

    def _setup_fasta(self, fasta_path: str) -> str:
        """Process and setup FASTA file"""
        if not fasta_path:
            return None
            
        try:
            local_path = os.path.abspath(fasta_path)
            if not os.path.exists(local_path):
                local_path = os.path.join(self.cwd, os.path.basename(fasta_path))
                if not os.path.exists(local_path):
                    shutil.copy2(fasta_path, local_path)
            
            return local_path
            
        except (TypeError, OSError) as e:
            if self.debug:
                print(f"Error setting up FASTA: {e}")
            return None

    def _set_database_paths(self) -> None:
        """Set paths to various databases and resources"""
        pass

    def print_time(self) -> None:
        """Print total runtime since initialization"""
        print(f'\n\nruntime: {datetime.now() - self.startTime}\n')

    def print_run_time(self, tool: str) -> None:
        """Print start time for a specific tool"""
        print(f'{bcolors.RED}\n{tool} running... {bcolors.ENDC}')
        now = datetime.now()
        print(f'{bcolors.WHITE}{now.strftime("%Y-%m-%d %H:%M:%S")}{bcolors.ENDC}')

# The LaTeX report machinery (analyze_logo_color, Banner, Latex_Report and
# its tectonic/pdflatex compile step) lived here until 2026-08. Reporting is
# HTML-first now: run_manifest.json + report.html, with report.pdf produced
# by WeasyPrint from that HTML when wanted (see reporting.py).


class Excel_Stats:
    """Generate Excel statistics reports"""

    def __init__(self, sample_name: str):
        """
        Initialize Excel stats report.
        
        Args:
            sample_name: Name of the sample for the report
        """
        self.sample_name = sample_name
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.excel_filename = f'{sample_name}_{self.date_stamp}_stats.xlsx'
        
        self.excel_dict = {
            'sample': sample_name,
            'date': self.date_stamp
        }

    def post_excel(self) -> None:
        """Save the Excel report to file"""
        df = pd.DataFrame.from_dict(self.excel_dict, orient='index').T
        df = df.set_index('sample')
        df.to_excel(self.excel_filename)
