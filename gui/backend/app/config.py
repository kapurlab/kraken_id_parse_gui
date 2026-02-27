import json
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CONFIG_PATH = DATA_DIR / "config.json"

HOME_DIR = Path.home()

def _detect_conda_env_path() -> str:
    """Auto-detect kraken_id_parse conda environment path"""
    try:
        # Try to get conda base path
        result = subprocess.run(
            ["conda", "info", "--base"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            conda_base = Path(result.stdout.strip())
            kraken_env = conda_base / "envs" / "kraken_id_parse"
            if kraken_env.exists():
                return str(kraken_env)
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
        pass

    # Try common conda locations
    common_conda_bases = [
        HOME_DIR / "miniconda3",
        HOME_DIR / "anaconda3",
        HOME_DIR / "mambaforge",
        Path("/opt/conda"),
        Path("/usr/local/anaconda3"),
        Path("/usr/local/miniconda3")
    ]

    for conda_base in common_conda_bases:
        kraken_env = conda_base / "envs" / "kraken_id_parse"
        if kraken_env.exists():
            return str(kraken_env)

    return ""

def _detect_repo_path() -> str:
    """Auto-detect repository path (go up 3 levels from backend/app)"""
    # From gui/backend/app -> gui -> repo_root
    repo_root = BASE_DIR.parent.parent
    if (repo_root / "bin" / "kraken_id_parse.py").exists():
        return str(repo_root)
    return ""

def _get_default_projects_path() -> str:
    """Get default projects path in user's home directory"""
    projects_path = HOME_DIR / "kraken_projects"
    projects_path.mkdir(exist_ok=True)
    return str(projects_path)

DEFAULTS: Dict[str, Any] = {
    "kraken_repo_path": "",
    "projects_root": "",
    "conda_env_path": "",
    "default_preset": ""
}

def _normalize_cfg(cfg: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(cfg)
    for key in ("kraken_repo_path", "projects_root", "conda_env_path", "default_preset"):
        if isinstance(out.get(key), str):
            out[key] = out[key].strip()
    return out

def _auto_detect_missing_paths(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Auto-detect any missing or invalid paths"""
    out = dict(cfg)

    # Auto-detect conda environment if missing or invalid
    conda_env = out.get("conda_env_path", "")
    if not conda_env or not Path(conda_env).exists():
        detected_conda = _detect_conda_env_path()
        if detected_conda:
            out["conda_env_path"] = detected_conda
            print(f"Auto-detected conda environment: {detected_conda}")

    # Auto-detect repository path if missing or invalid
    repo_path = out.get("kraken_repo_path", "")
    if not repo_path or not Path(repo_path).exists():
        detected_repo = _detect_repo_path()
        if detected_repo:
            out["kraken_repo_path"] = detected_repo
            print(f"Auto-detected repository path: {detected_repo}")

    # Set default projects path if missing
    projects_root = out.get("projects_root", "")
    if not projects_root:
        out["projects_root"] = _get_default_projects_path()
        print(f"Set default projects path: {out['projects_root']}")

    return out

def load_config() -> Dict[str, Any]:
    DATA_DIR.mkdir(exist_ok=True)
    if not CONFIG_PATH.exists():
        save_config(DEFAULTS)

    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    cfg = _normalize_cfg(cfg)

    # Auto-detect missing paths
    cfg = _auto_detect_missing_paths(cfg)

    # Persist auto-detected values
    save_config(cfg)
    return cfg

def save_config(cfg: Dict[str, Any]) -> None:
    DATA_DIR.mkdir(exist_ok=True)
    cfg = _normalize_cfg(cfg)
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, sort_keys=True)