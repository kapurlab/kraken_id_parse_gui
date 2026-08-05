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

# Where this deployment keeps its shared data. The launcher exports both values
# (bdtools bin/lib/site_paths.py resolves them from the site's configuration), so
# the tool ASKS the deployment instead of containing a path.
#
# These used to be absolute literals from the site this tool was written on. That
# is correct on exactly one machine and fails silently everywhere else: the
# directory simply isn't there, so a fresh user's Settings page looked configured
# while pointing at something that could never exist — and the only remedy was for
# every user to retype both paths, or for an admin to fake up the original site's
# directory layout.
_DB_ROOT = os.environ.get("BDTOOLS_DB_ROOT", "").strip()
_PROJECTS_ROOT = os.environ.get("BDTOOLS_SHARED_PROJECTS_ROOT", "").strip()

# Conventional layout under the database root, matching `bdtools setup-databases`.
_DB_SUBPATHS = {
    "kraken_db": "kraken2/k2_standard_08gb",
    "blast_db": "blast/ref_prok_rep_genomes",
}


def _default_db(key: str) -> str:
    """A default only when we can point at something real.

    An empty string means "not configured", which the GUI surfaces as a prompt to
    set it. That is a better first experience than a plausible-looking path that
    was never going to resolve. BLAST databases are a file *prefix* rather than a
    directory, so test the parent.
    """
    if not _DB_ROOT:
        return ""
    candidate = Path(_DB_ROOT) / _DB_SUBPATHS[key]
    if candidate.is_dir() or candidate.parent.is_dir():
        return str(candidate)
    return ""


_DEFAULT_SHARED_PROJECTS_ROOT = (
    _PROJECTS_ROOT if _PROJECTS_ROOT and Path(_PROJECTS_ROOT).is_dir() else ""
)

DEFAULTS: Dict[str, Any] = {
    "projects_root": str(Path.home() / "projects"),
    "shared_projects_root": _DEFAULT_SHARED_PROJECTS_ROOT,
    "saved_project_roots": [],
    "kraken_db": _default_db("kraken_db"),
    "blast_db": _default_db("blast_db"),
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
