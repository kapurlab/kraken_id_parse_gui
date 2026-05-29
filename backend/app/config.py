import json
import os
from pathlib import Path
from typing import Any, Dict

def _user_config_dir() -> Path:
    xdg = os.environ.get("XDG_CONFIG_HOME", "").strip()
    if xdg:
        return Path(xdg) / "kraken_id_parse_gui"
    return Path.home() / ".config" / "kraken_id_parse_gui"


DATA_DIR = _user_config_dir()
CONFIG_PATH = DATA_DIR / "config.json"

_SHARED_PROJECTS_ROOT = Path("/srv/kapurlab/projects")
_DEFAULT_SHARED_PROJECTS_ROOT = (
    str(_SHARED_PROJECTS_ROOT) if _SHARED_PROJECTS_ROOT.is_dir() else ""
)

DEFAULTS: Dict[str, Any] = {
    "projects_root": str(Path.home() / "projects"),
    "shared_projects_root": _DEFAULT_SHARED_PROJECTS_ROOT,
    "kraken_db": "/srv/kapurlab/databases/kraken2/k2_standard_08gb",
    "blast_db": "/srv/kapurlab/databases/blast/ref_prok_rep_genomes",
}


def load_config() -> Dict[str, Any]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not CONFIG_PATH.exists():
        save_config(DEFAULTS)
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    for k, v in DEFAULTS.items():
        cfg.setdefault(k, v)
    return cfg


def save_config(cfg: Dict[str, Any]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, sort_keys=True)
