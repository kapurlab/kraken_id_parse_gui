#!/usr/bin/env python3

"""
Core setup and utility classes for bioinformatics pipeline tools.
Provides base classes for file handling, reporting, and formatting.
"""

import os
import sys
import shutil
import re
from typing import Optional, Dict, List, Union, Tuple
import pandas as pd
import multiprocessing
from datetime import datetime
from pathlib import Path

try:
    from PIL import Image
    import numpy as np
    import colorsys
except ImportError:
    pass

def resolve_cpus(explicit: Optional[Union[int, str]] = None) -> int:
    """Resolve the worker thread count for the CPU-heavy steps (kraken2, BLAST,
    SPAdes, bwa).

    Priority: explicit value > ``KIP_THREADS`` env var > ``SLURM_CPUS_PER_TASK``
    > CPU affinity (cgroup-aware, Linux) > ``multiprocessing.cpu_count()``.

    An explicit allocation (``--threads N`` / ``--override threads=N``, or
    ``sbatch --cpus-per-task``) is honored EXACTLY so many samples can be
    packed onto one node. The auto-detect path (a direct/Mac run with nothing
    set) keeps the historical behavior of leaving 2 cores free so an interactive
    machine stays responsive.
    """
    # Explicit request (--threads / KIP_THREADS) wins and is honored exactly,
    # even if 1.
    for val in (explicit, os.environ.get('KIP_THREADS')):
        try:
            n = int(val)
        except (TypeError, ValueError):
            continue
        if n > 0:
            return n
    # SLURM_CPUS_PER_TASK is honored only when it reflects a real per-task
    # allocation (>1). Some clusters export "1" alongside a whole-node --ntasks
    # request; treating that as authoritative would silently single-thread the
    # existing production SLURM runs, so we ignore it and fall through to affinity.
    #
    # A bdtools dashboard session declares a core budget for every tool it
    # launches (BDTOOLS_SESSION_CORES). Without this cap a single kraken2/bwa/
    # blastn run asks for every core on a shared box, and a room of concurrent
    # users turns into thread-thrash. Only an explicit request (above) is exempt;
    # absent the variable, behavior is unchanged.
    try:
        budget = int(os.environ.get('BDTOOLS_SESSION_CORES', '').strip())
    except ValueError:
        budget = 0

    def _capped(n):
        return max(1, min(n, budget)) if budget > 0 else n

    try:
        n = int(os.environ.get('SLURM_CPUS_PER_TASK'))
        if n > 1:
            return _capped(n)
    except (TypeError, ValueError):
        pass
    try:
        avail = len(os.sched_getaffinity(0))  # honors cgroup/Slurm binding
    except AttributeError:
        avail = multiprocessing.cpu_count()   # macOS has no sched_getaffinity
    return _capped(max(1, avail - 2))


def color_enabled(stream=None) -> bool:
    """
    Decide whether to emit ANSI colors. True only for an interactive terminal.
    When stdout is redirected to a file (e.g. a SLURM .out log) colors are OFF so
    the log stays clean instead of full of raw escape codes. Honors the
    conventional NO_COLOR (force off) and FORCE_COLOR (force on) env vars.
    """
    if os.environ.get('NO_COLOR') is not None:
        return False
    if os.environ.get('FORCE_COLOR') is not None:
        return True
    stream = stream if stream is not None else sys.stdout
    try:
        return stream.isatty() and os.environ.get('TERM', '') != 'dumb'
    except Exception:
        return False


class bcolors:
    """ANSI color codes, blanked automatically when output is not a terminal
    (so SLURM/redirected logs are not littered with escape sequences)."""
    _enabled = color_enabled()
    if _enabled:
        PURPLE = '\033[95m'
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        WHITE = '\033[37m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        ENDC = '\033[0m'
    else:
        PURPLE = BLUE = GREEN = YELLOW = RED = WHITE = BOLD = UNDERLINE = ENDC = ''


class UI:
    """
    TTY-aware console formatting for pipeline status lines.

    On an interactive terminal these render with color, bold, and box-drawing
    rules; when redirected to a file (SLURM logs) they degrade to clean ASCII
    with no escape codes, so both surfaces look intentional. Each helper returns
    a single line, which keeps them tidy even when passed through the logger
    (whose 'timestamp - INFO -' prefix would break a multi-line box).
    """
    color = bcolors._enabled

    @staticmethod
    def banner(msg: str) -> str:
        if UI.color:
            return f"\n{bcolors.BOLD}{bcolors.GREEN}{'═' * 60}{bcolors.ENDC}\n" \
                   f"{bcolors.BOLD}{bcolors.GREEN}  {msg}{bcolors.ENDC}\n" \
                   f"{bcolors.BOLD}{bcolors.GREEN}{'═' * 60}{bcolors.ENDC}"
        bar = '=' * 60
        return f"\n{bar}\n  {msg}\n{bar}"

    @staticmethod
    def step(n, total, msg: str) -> str:
        if UI.color:
            return f"{bcolors.BOLD}{bcolors.BLUE}━━▶ Step {n}/{total} · {msg}{bcolors.ENDC}"
        return f"==> Step {n}/{total}: {msg}"

    @staticmethod
    def ok(msg: str) -> str:
        return f"{bcolors.GREEN}✓{bcolors.ENDC} {msg}" if UI.color else f"[OK] {msg}"

    @staticmethod
    def warn(msg: str) -> str:
        return f"{bcolors.YELLOW}⚠ {msg}{bcolors.ENDC}" if UI.color else f"[WARN] {msg}"

    @staticmethod
    def error(msg: str) -> str:
        return f"{bcolors.RED}✗ {msg}{bcolors.ENDC}" if UI.color else f"[ERROR] {msg}"


def move_overwrite(src: str, dst_dir: str) -> str:
    """
    Move ``src`` into ``dst_dir``, replacing any existing entry of the same name.

    Re-run safe: plain ``shutil.move`` raises "Destination path ... already exists"
    when the pipeline is run a second time in the same folder. Returns the final path.
    """
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, os.path.basename(src.rstrip('/')))
    if os.path.exists(dst):
        shutil.rmtree(dst) if os.path.isdir(dst) else os.remove(dst)
    shutil.move(src, dst)
    return dst


def safe_move(src, dst):
    """shutil.move that silently overwrites an existing destination file or directory."""
    src_path = Path(src)
    dst_path = Path(dst)
    actual_dst = dst_path / src_path.name if dst_path.is_dir() else dst_path
    if actual_dst.exists():
        shutil.rmtree(actual_dst) if actual_dst.is_dir() else actual_dst.unlink()
    shutil.move(str(src), str(dst))


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
        self.cpus = resolve_cpus()
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


class SummaryLog:
    """Condensed run summary log that records tool versions, section timings,
    and any errors encountered during a pipeline run."""

    def __init__(self, sample_name: str, script_version: str, args=None):
        self.sample_name = sample_name
        self.script_version = script_version
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.log_file = f'{sample_name}_{self.date_stamp}_summary.log'
        self._versions: Dict[str, str] = {}
        self._sections: list = []  # list of (name, start, end, status)
        self._errors: list = []    # list of (timestamp, message)
        self._current_section: Optional[tuple] = None
        self._args = args
        self.git_branch: str = 'unknown'
        self.git_commit: str = 'unknown'

    # -- tool versions --------------------------------------------------------

    def collect_versions(self):
        """Probe external tools and record their version strings."""
        import subprocess
        import shutil as _shutil

        self._collect_git_info()

        self._versions['Python'] = self._run_version('python3', ['python3', '--version'])
        self._versions['Biopython'] = self._python_module_version('Bio')
        self._versions['pandas'] = self._python_module_version('pandas')

        tool_cmds = [
            ('kraken2',    ['kraken2', '--version']),
            ('SPAdes',     ['spades.py', '--version']),
            ('seqkit',     ['seqkit', 'version']),
            ('BWA',        ['bwa']),
            ('samtools',   ['samtools', '--version']),
            ('blastn',     ['blastn', '-version']),
            ('Picard',     ['picard', 'MarkDuplicates', '--version']),
            ('freebayes',  ['freebayes', '--version']),
            ('pigz',       ['pigz', '--version']),
        ]

        for name, cmd in tool_cmds:
            if _shutil.which(cmd[0]) is None:
                continue
            self._versions[name] = self._run_version(name, cmd)

        # bracken — has no simple --version flag; extract from conda list
        if _shutil.which('bracken'):
            ver = self._run_version('bracken', ['conda', 'list', 'bracken'])
            # Parse conda list output: "bracken  2.9  ..."
            for line in ver.splitlines() if '\n' in ver else [ver]:
                if line.startswith('bracken'):
                    parts = line.split()
                    if len(parts) >= 2:
                        self._versions['bracken'] = f'bracken {parts[1]}'
                        break
            else:
                self._versions['bracken'] = 'installed (version unknown)'

    def _collect_git_info(self):
        """Capture the git branch and short commit of the running code (for QA)."""
        import subprocess
        repo_dir = os.path.dirname(os.path.abspath(__file__))
        try:
            self.git_branch = subprocess.run(
                ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                cwd=repo_dir, capture_output=True, text=True, timeout=10,
            ).stdout.strip() or 'unknown'
            self.git_commit = subprocess.run(
                ['git', 'rev-parse', '--short', 'HEAD'],
                cwd=repo_dir, capture_output=True, text=True, timeout=10,
            ).stdout.strip() or 'unknown'
        except Exception:
            self.git_branch = 'unknown'
            self.git_commit = 'unknown'
        if self.git_branch == 'unknown' or self.git_commit == 'unknown':
            branch, commit = self._git_info_from_files(repo_dir)
            if branch and self.git_branch == 'unknown':
                self.git_branch = branch
            if commit and self.git_commit == 'unknown':
                self.git_commit = commit

    @staticmethod
    def _git_info_from_files(start_dir: str):
        """Read (branch, short commit) straight from the .git directory. Used when
        the git command cannot run - e.g. an osx-64 (Rosetta) conda env on Apple
        Silicon, where /usr/bin/git fails to load the arm64-only Command Line
        Tools. Returns (None, None) outside a repository."""
        d = os.path.abspath(start_dir)
        while not os.path.exists(os.path.join(d, '.git')):
            parent = os.path.dirname(d)
            if parent == d:
                return None, None
            d = parent
        git_dir = os.path.join(d, '.git')
        try:
            if os.path.isfile(git_dir):   # worktree / submodule: "gitdir: <path>"
                with open(git_dir) as fh:
                    git_dir = os.path.join(d, fh.read().split('gitdir:', 1)[1].strip())
            with open(os.path.join(git_dir, 'HEAD')) as fh:
                head = fh.read().strip()
            if not head.startswith('ref: '):
                return 'HEAD', head[:7]      # detached HEAD
            ref = head[5:]
            branch = ref[len('refs/heads/'):] if ref.startswith('refs/heads/') else ref
            common = os.path.join(git_dir, 'commondir')   # worktrees keep refs here
            if os.path.exists(common):
                with open(common) as fh:
                    git_dir = os.path.normpath(os.path.join(git_dir, fh.read().strip()))
            commit = None
            ref_path = os.path.join(git_dir, ref)
            if os.path.exists(ref_path):
                with open(ref_path) as fh:
                    commit = fh.read().strip()
            else:
                packed = os.path.join(git_dir, 'packed-refs')
                if os.path.exists(packed):
                    with open(packed) as fh:
                        for line in fh:
                            parts = line.split()
                            if len(parts) == 2 and parts[1] == ref:
                                commit = parts[0]
            return branch, commit[:7] if commit else None
        except (OSError, IndexError):
            return None, None

    @staticmethod
    def _run_version(name: str, cmd: list) -> str:
        import subprocess
        import re as _re
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            output = (result.stdout + result.stderr).strip()
            # Prefer lines containing a version number (e.g. 1.21, v4.2.0)
            for line in output.splitlines():
                line = line.strip()
                if line and _re.search(r'[vV]?\d+\.\d+', line):
                    return line
            # Fallback to first non-empty line
            for line in output.splitlines():
                line = line.strip()
                if line:
                    return line
            return 'unknown'
        except Exception:
            return 'not found'

    @staticmethod
    def _python_module_version(module_name: str) -> str:
        try:
            mod = __import__(module_name)
            return f'{module_name} {mod.__version__}'
        except Exception:
            return 'not installed'

    # -- section timing -------------------------------------------------------

    def start_section(self, name: str):
        """Mark the start of a pipeline section."""
        if self._current_section is not None:
            # Auto-close the previous section
            self.end_section()
        self._current_section = (name, datetime.now())

    def end_section(self, status: str = 'OK'):
        """Mark the end of the current pipeline section."""
        if self._current_section is None:
            return
        name, start = self._current_section
        self._sections.append((name, start, datetime.now(), status))
        self._current_section = None

    # -- error capture --------------------------------------------------------

    def log_error(self, message: str):
        """Record an error for the summary."""
        self._errors.append((datetime.now().strftime('%H:%M:%S'), message))

    # -- write the summary file -----------------------------------------------

    def write(self):
        """Write the condensed summary log to disk."""
        lines = []
        sep = '-' * 70

        # Header
        lines.append(sep)
        lines.append(f'PIPELINE SUMMARY — {self.sample_name}')
        lines.append(sep)
        lines.append(f'Date            : {self.date_stamp}')
        lines.append(f'Script version  : {self.script_version}')
        lines.append(f'Git branch      : {self.git_branch}')
        lines.append(f'Git commit      : {self.git_commit}')
        if self._args:
            lines.append(f'Target taxon    : {getattr(self._args, "taxon", "N/A")}')
            lines.append(f'BLAST database  : {getattr(self._args, "blast_db", "N/A")}')
            lines.append(f'Kraken database : {getattr(self._args, "kraken_db", "N/A")}')
            lines.append(f'Fast mode       : {getattr(self._args, "fast_mode", False)}')
            lines.append(f'Debug mode      : {getattr(self._args, "debug", False)}')
        lines.append('')

        # Tool versions
        lines.append(sep)
        lines.append('TOOL VERSIONS')
        lines.append(sep)
        for tool, ver in self._versions.items():
            lines.append(f'  {tool:<16}: {ver}')
        lines.append('')

        # Section timings
        lines.append(sep)
        lines.append(f'{"SECTION":<42} {"TIME":>10}  {"STATUS"}')
        lines.append(sep)
        total_seconds = 0
        for name, start, end, status in self._sections:
            elapsed = end - start
            total_seconds += elapsed.total_seconds()
            elapsed_str = str(elapsed).split('.')[0]  # drop microseconds
            status_marker = status if status == 'OK' else f'** {status} **'
            lines.append(f'  {name:<40} {elapsed_str:>10}  {status_marker}')
        lines.append(sep)
        # Total
        m, s = divmod(int(total_seconds), 60)
        h, m = divmod(m, 60)
        lines.append(f'  {"TOTAL":<40} {h:02d}:{m:02d}:{s:02d}')
        lines.append('')

        # Errors
        if self._errors:
            lines.append(sep)
            lines.append(f'ERRORS / WARNINGS  ({len(self._errors)})')
            lines.append(sep)
            for ts, msg in self._errors:
                lines.append(f'  [{ts}] {msg}')
        else:
            lines.append(sep)
            lines.append('No errors recorded.')
        lines.append(sep)
        lines.append('')

        with open(self.log_file, 'w') as f:
            f.write('\n'.join(lines))