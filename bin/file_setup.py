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
    """
    Apply a matplotlib style with graceful fallback across versions.
    Returns the pyplot module for convenience.
    """
    import matplotlib
    # Force non-GUI backend on macOS to avoid AppKit font crashes
    if not os.environ.get("MPLBACKEND"):
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # local import after backend selection
    candidates = style_candidates or [
        'seaborn-v0_8-colorblind',
        'seaborn-colorblind',
        'seaborn-v0_8',
        'seaborn',
        'ggplot',
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

def analyze_logo_color(logo_path: str) -> str:
    """
    Analyze the logo image to extract a dominant color for use in banners.
    
    Args:
        logo_path: Path to the logo image file
        
    Returns:
        str: RGB color code in format "r, g, b"
    """
    # Default color if analysis fails (original shade of blue)
    default_color = "56, 68, 117"
    
    try:
        # Check if PIL is available
        if 'Image' not in globals() or not os.path.exists(logo_path):
            return default_color
            
        # Open the logo and convert to RGB mode
        img = Image.open(logo_path).convert('RGB')
        
        # Resize for faster processing
        img = img.resize((100, 100))
        
        # Convert to numpy array
        img_array = np.array(img)
        
        # Skip transparent or white background pixels
        # Create a mask for pixels that are not white/transparent
        mask = ~np.all(img_array > 240, axis=2)
        
        # If we have valid pixels, use them for color calculation
        if np.any(mask):
            # Extract only non-white pixels
            valid_pixels = img_array[mask]
            
            # Calculate mean color from valid pixels
            avg_color = np.mean(valid_pixels, axis=0).astype(int)
        else:
            # Fallback to using all pixels if no valid pixels found
            avg_color = np.mean(img_array, axis=(0, 1)).astype(int)
            
        r, g, b = avg_color
        
        # Convert RGB to HSV for better color manipulation
        h, s, v = colorsys.rgb_to_hsv(r/255, g/255, b/255)
        
        # Enhance saturation for more vibrant color
        s = min(1.0, s * 1.3)
        
        # Ensure value (brightness) is appropriate for white text
        # Lower brightness for better contrast with white text
        v = min(0.85, v)
        
        # Convert back to RGB
        enhanced_r, enhanced_g, enhanced_b = colorsys.hsv_to_rgb(h, s, v)
        enhanced_color = (int(enhanced_r * 255), int(enhanced_g * 255), int(enhanced_b * 255))
        
        # Check if the color is too light for white text
        brightness = (enhanced_r * 299 + enhanced_g * 587 + enhanced_b * 114) / 1000
        if brightness > 0.7:  # If too bright
            # Darken the color more dramatically
            v = max(0.3, v - 0.4)
            enhanced_r, enhanced_g, enhanced_b = colorsys.hsv_to_rgb(h, s, v)
            enhanced_color = (int(enhanced_r * 255), int(enhanced_g * 255), int(enhanced_b * 255))
        
        return f"{enhanced_color[0]}, {enhanced_color[1]}, {enhanced_color[2]}"
        
    except Exception as e:
        print(f"Error analyzing logo color: {e}")
        return default_color

class Banner:
    """Generate banner images for reports with modern design elements"""

    def __init__(self, title: str, hexcode: str = "56, 68, 117"):
        """
        Create a professional banner image with specified title and color.
        
        Args:
            title: Text to display on banner
            hexcode: RGB color code for banner background
        
        Raises:
            ImportError: If svgwrite or cairosvg are not available
        """
        if not HAS_SVG_SUPPORT:
            raise ImportError("svgwrite, cairosvg, and PIL are required for banner generation")
            
        width = 2600
        height = 90
        
        # Parse the RGB values for gradient creation
        r, g, b = map(int, hexcode.split(','))
        # Create slightly lighter shade for gradient
        light_r = min(r + 20, 255)
        light_g = min(g + 20, 255)
        light_b = min(b + 20, 255)
        
        # Create SVG
        svgimg = svgwrite.Drawing(size=(width, height))
        
        # Create a basic rounded rectangle with gradient fill
        gradient = svgimg.linearGradient(
            start=(0, 0), end=(0, height), 
            id="banner_gradient",
            gradientUnits="userSpaceOnUse"
        )
        gradient.add_stop_color(offset='0%', color=f'rgb({light_r}, {light_g}, {light_b})')
        gradient.add_stop_color(offset='100%', color=f'rgb({r}, {g}, {b})')
        svgimg.defs.add(gradient)
        
        # Create a basic rounded rectangle
        svgimg.add(svgimg.rect([0, 0], [width, height], 
                              rx=10, ry=10,  # Rounded corners
                              fill="url(#banner_gradient)", 
                              stroke="none"))
        
        # Add a thin highlight line at the top for a polished look
        svgimg.add(svgimg.rect([2, 2], [width-4, 3], 
                              rx=8, ry=3,
                              fill="white",
                              stroke="none",
                              opacity=0.3))
        
        # Add a subtle shadow at the bottom
        svgimg.add(svgimg.rect([2, height-5], [width-4, 3], 
                              rx=0, ry=0,
                              fill="black",
                              stroke="none",
                              opacity=0.2))
        
        # Add title text
        svgimg.add(svgimg.text(title, 
                              insert=(30, 60), 
                              fill='white',
                              font_size='50px',
                              font_weight='bold'))
        
        # Save and convert
        temp_svg = f'temp_{title.replace(" ", "_")}.svg'
        svgimg.saveas(temp_svg)
        
        with open(temp_svg, 'r') as content_file:
            content = content_file.read()
            
        output_name = f'{title.replace(" ", "_")}-banner.png'
        svg2png(bytestring=content, write_to=output_name)
        
        # Cleanup
        os.remove(temp_svg)
        self.banner = str(Path(os.getcwd()) / output_name)
        
class Latex_Report:
    """Generate LaTeX reports with consistent formatting"""

    def __init__(self, sample_name: str, report_description: Optional[str] = None, logo: Optional[str] = None):
        """
        Initialize LaTeX report with standard formatting.
        
        Args:
            sample_name: Name of the sample for the report
            report_description: Optional description to include in header
            logo: Optional path to logo file
        """
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.tex_file = f'{sample_name}_{date_stamp}_report.tex'
        self.logo = logo
        
        # Extract color from logo if provided, otherwise use default
        self.banner_color = analyze_logo_color(logo) if logo else "56, 68, 117"
        
        self.tex = open(self.tex_file, 'w')
        self._write_preamble()
        self._write_header(sample_name, report_description, logo)
        
    def _write_preamble(self) -> None:
        """Write LaTeX preamble with package imports and styling"""
        preamble = [
            r'\documentclass{article}',
            r'\usepackage[margin=1in]{geometry}',
            r'\usepackage{adjustbox}',
            r'\usepackage{float}',
            r'\usepackage{graphicx}',
            r'\usepackage{fancyhdr}',
            r'\usepackage{hyperref}',
            r'\usepackage{longtable}',
            r'\usepackage[scaled]{helvet}',
            r'\renewcommand\familydefault{\sfdefault}',
            r'\usepackage[T1]{fontenc}',
            r'\usepackage{xcolor}'
        ]
        print('\n'.join(preamble), file=self.tex)
        
    def _write_header(self, sample_name: str, report_description: Optional[str], logo: Optional[str]):
        """Write report header with logo and title"""
        current_date = datetime.now().strftime('%B %d, %Y')
        
        # Setup the page style without logo for all pages
        header = [
            r'\pagestyle{fancy}',
            r'\fancyhf{}', # Clear all header/footer fields
            r'\rhead{\textbf{\large ' + current_date + r'}}',
            r'\renewcommand{\headrulewidth}{0.4pt}',
            r'\begin{document}'
        ]
        
        # Add title page content with logo
        if logo:
            header.extend([
                r'\thispagestyle{empty}', # No headers on first page
                r'\vspace*{-0.5cm}', # Move up to minimize space at top
                r'\noindent\begin{minipage}[t]{0.5\textwidth}',
                r'\includegraphics[width=1.0\textwidth,keepaspectratio]{' + f'{logo}' + r'}',
                r'\end{minipage}',
                r'\begin{minipage}[t]{0.5\textwidth}',
                r'\begin{flushright}',
                r'\textbf{\large ' + current_date + r'}',
                r'\end{flushright}',
                r'\end{minipage}',
                r'\vspace{0.2cm}',
                r'\hrule', # Add horizontal line
                r'\vspace{0.2cm}'
            ])
        else:
            header.extend([
                r'\thispagestyle{empty}',
                r'\vspace*{-0.5cm}',
                r'\begin{flushright}',
                r'\textbf{\large ' + current_date + r'}',
                r'\end{flushright}',
                r'\vspace{0.2cm}',
                r'\hrule',
                r'\vspace{0.2cm}'
            ])
        
        # Add sample name
        header.append(r'\noindent\textbf{\large{\fontfamily{\sfdefault}\selectfont Sample: ' + sample_name + r'}}')
        
        # Add report description if provided
        if report_description:
            header.append(r'\vspace{0.2cm}')
            header.append(r'\noindent\textbf{\large{' + report_description + r'}}')
        
        # Add some space before content starts
        header.append(r'\vspace{0.2cm}')
        
        print('\n'.join(header), file=self.tex)

    def add_section(self, title: str) -> None:
        """Add a new section with a banner to the report
        
        Args:
            title: Section title
        """
        # Create banner with title and matching color from logo analysis
        banner = Banner(title, hexcode=self.banner_color)
        
        # Add the banner to the LaTeX document
        print(r'\begin{figure}[H]', file=self.tex)
        print(r'\centering', file=self.tex)
        print(r'\includegraphics[width=\textwidth]{' + banner.banner + '}', file=self.tex)
        print(r'\end{figure}', file=self.tex)

    def add_table_section(self, title: str) -> None:
        """Add a section with a smaller banner suitable for tables
        
        Args:
            title: Section title
        """
        # Create banner with title and matching color from logo analysis
        banner = Banner(title, hexcode=self.banner_color)
        
        # Add the banner to the LaTeX document, but with smaller size to fit tables
        print(r'\begin{center}', file=self.tex)
        print(r'\includegraphics[width=\textwidth]{' + banner.banner + '}', file=self.tex)
        print(r'\end{center}', file=self.tex)

    def latex_ending(self) -> None:
        """Finalize and compile the LaTeX document"""
        print(r'\end{document}', file=self.tex)
        self.tex.close()

        # Run pdflatex twice for proper rendering
        log_file = self.tex_file.replace('.tex', '_pdflatex.log')
        for i in range(2):
            ret = os.system(f'pdflatex -interaction=nonstopmode {self.tex_file} > {log_file} 2>&1')
            if ret != 0:
                print(f"  pdflatex pass {i+1} returned exit code {ret}")
        # Print any errors/warnings from the log
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    if any(kw in line for kw in ['!', 'Error', 'Warning', 'Missing']):
                        print(f"  pdflatex: {line.rstrip()}")
        except FileNotFoundError:
            pass

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
