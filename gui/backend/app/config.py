import json
from pathlib import Path
from typing import Any, Dict

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CONFIG_PATH = DATA_DIR / "config.json"

HOME_DIR = Path.home()
DEFAULTS: Dict[str, Any] = {
    "kraken_repo_path": str(HOME_DIR / "kraken" / "pipeline"),
    "projects_root": str(HOME_DIR / "kraken" / "projects"),
    "conda_env_path": "",
    "default_preset": ""
}


def _normalize_cfg(cfg: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(cfg)
    for key in ("kraken_repo_path", "projects_root", "conda_env_path", "default_preset"):
        if isinstance(out.get(key), str):
            out[key] = out[key].strip()
    return out


def load_config() -> Dict[str, Any]:
    DATA_DIR.mkdir(exist_ok=True)
    if not CONFIG_PATH.exists():
        save_config(DEFAULTS)
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    cfg = _normalize_cfg(cfg)
    # Light validation for common dev machine paths
    kraken_repo = Path(cfg.get("kraken_repo_path", ""))
    projects_root = Path(cfg.get("projects_root", ""))
    if str(kraken_repo).startswith(str(HOME_DIR)) and not kraken_repo.exists():
        cfg["kraken_repo_path"] = DEFAULTS["kraken_repo_path"]
    if str(projects_root).startswith(str(HOME_DIR)) and not projects_root.exists():
        cfg["projects_root"] = DEFAULTS["projects_root"]
    # Persist normalized values if they differ
    save_config(cfg)
    return cfg


def save_config(cfg: Dict[str, Any]) -> None:
    DATA_DIR.mkdir(exist_ok=True)
    cfg = _normalize_cfg(cfg)
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, sort_keys=True)
