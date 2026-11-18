"""
Kraken ID Parse GUI — FastAPI backend.

Serves the React SPA from frontend/dist/ and provides:
  /api/projects        — list shared + personal projects (FASTQ browser)
  /api/projects/{n}/samples — list FASTQ pairs in project/download/
  /api/config          — get/set user config (DB paths)
  /api/run             — start a kraken_id_parse.py run
  /api/jobs            — list running/completed jobs
  /api/jobs/{id}       — job detail
  /api/jobs/{id}/log   — SSE stream of the job log

All URLs served from / (uvicorn is behind OOD rnode proxy — relative paths only).
"""

import asyncio
import glob
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import aiofiles
from fastapi import FastAPI, File, HTTPException, Query, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .config import DEFAULTS, load_config, save_config
from .jobs import JobManager
from .request_safety import install_request_safety
from .sra import (
    SRAExpansionError,
    build_download_script,
    expand_accessions_with_mapping,
    is_valid_accession,
    write_crosswalk_tsv,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parent.parent          # /srv/kapurlab/tools/kraken_id_parse_gui
_BIN_DIR = _REPO_ROOT / "bin"
_FRONTEND_DIST = _REPO_ROOT / "frontend" / "dist"

# Shared project root — resolved by config.py from the launcher's env
# (BDTOOLS_SHARED_PROJECTS_ROOT), "" when this deployment has none. Keep a
# definitely-nonexistent placeholder when unset so the listing loops skip it
# without special cases.
_SHARED_PROJECTS = (
    Path(DEFAULTS["shared_projects_root"])
    if DEFAULTS.get("shared_projects_root")
    else Path.home() / ".local" / "share" / "kraken_id_parse_gui" / "no-shared-projects"
)

# Jobs log directory (inside repo so it survives across sessions)
_JOBS_DIR = _REPO_ROOT / "backend" / "jobs"

# Target names (the "Target taxon" dropdown, shared with the vSNP GUI): the
# list that ships in config/taxa.yaml plus the names users save, which live in
# the per-user config folder, never in this checkout. See target_names.py.
from . import target_names  # noqa: E402

try:
    _MIGRATED = target_names.migrate_checkout_additions(_REPO_ROOT)
    if _MIGRATED.get("carried"):
        logger.info("Moved %d saved target name(s) out of config/taxa.yaml into %s",
                    len(_MIGRATED["carried"]), target_names.saved_path())
except Exception as exc:  # never let a carry-over stop the app starting
    logger.warning("Target-name carry-over skipped: %s", exc)

# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------
def _resolve_app_version() -> str:
    """Version of the deployed checkout — the exact string the Diagnostic
    Tools Dashboard shows for this tool (`git describe --tags --always`,
    the same command bdtools runs). Resolved once at startup; empty when
    git or the .git dir is unavailable, in which case the frontend falls
    back to its built-in constant."""
    try:
        out = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "describe", "--tags", "--always"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


APP_VERSION = _resolve_app_version()


# ---------------------------------------------------------------------------
# App & job manager
# ---------------------------------------------------------------------------
app = FastAPI(title="Kraken ID Parse GUI")
install_request_safety(app)

job_manager = JobManager(_JOBS_DIR)


# ---------------------------------------------------------------------------
# Helpers — project listing
# ---------------------------------------------------------------------------
_SCOPE_SHARED = "shared"
_SCOPE_PERSONAL = "personal"


def _safe_mtime(p: Path) -> float:
    try:
        return p.stat().st_mtime if p.is_dir() else 0
    except PermissionError:
        return 0


def _count_project_reads(download_dir: Path, step1_dir: Path) -> int:
    """Count input read files (*.fastq.gz) across download/ and step1/.

    Native projects keep reads in download/; vSNP/Roar-imported projects keep
    them in step1/<sample>/ (and may symlink them into download/). Count the
    union, deduped by resolved path (so a download/ symlink pointing at a step1
    read isn't double-counted), skipping *_unmapped_* (the unmapped-read subset
    vSNP3 emits alongside the real reads — not an input read set). step1 is
    globbed one level deep (reads live directly under each sample dir)."""
    seen: set = set()
    candidates = []
    if download_dir.is_dir():
        candidates += download_dir.rglob("*.fastq.gz")
    if step1_dir.is_dir():
        candidates += step1_dir.glob("*/*.fastq.gz")
    for f in candidates:
        if "_unmapped_" in f.name:
            continue
        try:
            key = f.resolve()
        except OSError:
            key = f
        seen.add(key)
    return len(seen)


def _list_projects_from_root(root: Path, scope: str) -> List[Dict]:
    if not root.is_dir():
        return []
    projects = []
    try:
        entries = sorted(root.iterdir(), key=_safe_mtime, reverse=True)
    except PermissionError:
        return []
    for p in entries:
        try:
            if not p.is_dir() or p.name.startswith("."):
                continue
        except PermissionError:
            continue
        download_dir = p / "download"
        try:
            fastq_count = _count_project_reads(download_dir, p / "step1")
        except PermissionError:
            fastq_count = -1  # signals "no access" to frontend
        kraken_runs = []
        kraken_dir = p / "kraken"
        try:
            if kraken_dir.is_dir():
                kraken_runs = [d.name for d in sorted(kraken_dir.iterdir()) if d.is_dir()]
        except PermissionError:
            pass
        projects.append({
            "name": p.name,
            "path": str(p),
            "scope": scope,
            "fastq_count": fastq_count,
            "kraken_runs": kraken_runs,
        })
    return projects


def _get_project_dir(name: str) -> Optional[Path]:
    """Find a project dir in shared then personal roots."""
    if "/" in name or name.startswith("."):
        return None
    cfg = load_config()
    for root in [_SHARED_PROJECTS, Path(cfg.get("projects_root", ""))]:
        candidate = root / name
        if candidate.is_dir():
            return candidate
    return None


# ---------------------------------------------------------------------------
# Project creation.
#
# A project created here uses the SAME on-disk skeleton vSNP GUI creates, so a
# project made in Kraken ID Parse is immediately usable in vSNP (and vice
# versa) — both tools share /srv/kapurlab/projects and per-user ~/projects and
# list whatever is on disk. We add the kraken/ subdir up front so the sample
# browser and results endpoints have a stable layout.
# ---------------------------------------------------------------------------
_PROJECT_NAME_OK_CHARSET = re.compile(r"^[A-Za-z0-9._-]+$")


def _normalize_project_name(name: str) -> str:
    """Filesystem-safe project dir name. Mirrors vsnp_gui's rules so a name
    accepted in one tool is accepted in the other: spaces auto-convert to
    underscores (downstream CLIs don't quote paths), other unsafe chars are
    rejected with a clear message."""
    if not isinstance(name, str):
        raise ValueError("Project name must be a string")
    cleaned = re.sub(r"\s+", "_", name.strip())
    if not cleaned:
        raise ValueError("Project name is empty")
    if cleaned.startswith("."):
        raise ValueError("Project name cannot start with '.'")
    if len(cleaned) > 100:
        raise ValueError("Project name too long (max 100 characters)")
    if not _PROJECT_NAME_OK_CHARSET.match(cleaned):
        bad = sorted(set(ch for ch in cleaned if not re.match(r"[A-Za-z0-9._-]", ch)))
        raise ValueError(
            f"Project name contains unsupported characters: {''.join(bad)!r}. "
            "Only letters, digits, _ - . are allowed (spaces become underscores)."
        )
    return cleaned


def _ensure_project_dirs(project_dir: Path) -> None:
    (project_dir / "download").mkdir(parents=True, exist_ok=True)
    (project_dir / "kraken").mkdir(parents=True, exist_ok=True)
    # vSNP-compatible layout so the project is shared cleanly between tools.
    (project_dir / "step1").mkdir(parents=True, exist_ok=True)
    (project_dir / "step2" / "vcf_source").mkdir(parents=True, exist_ok=True)
    (project_dir / f"{project_dir.name}_VCFs").mkdir(parents=True, exist_ok=True)


def _create_project(name: str, scope: str) -> Path:
    """Create a project under the requested scope ('personal' or 'shared').

    Personal projects land under the user's configured projects_root (defaults
    to ~/projects). Shared projects land under /srv/kapurlab/projects, which is
    only writable where the lab has granted it — a PermissionError there is
    surfaced as a clear 400 rather than a stack trace.
    """
    name = _normalize_project_name(name)
    cfg = load_config()
    if scope == _SCOPE_SHARED:
        root = _SHARED_PROJECTS
    else:
        root = Path(cfg.get("projects_root", "") or (Path.home() / "projects"))
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"Cannot create projects root {root}: {exc}")
    project_dir = root / name
    if project_dir.exists():
        raise ValueError(f"Project already exists: {name}")
    try:
        _ensure_project_dirs(project_dir)
    except PermissionError:
        raise ValueError(
            f"No permission to create a project under {root}. "
            "Shared projects require lab write access; create it as a personal "
            "project instead."
        )
    meta = {"name": name, "created_at": _now_iso(), "status": "created"}
    try:
        with open(project_dir / "project.json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, sort_keys=True)
    except OSError:
        pass
    return project_dir


def _now_iso() -> str:
    from datetime import datetime
    return datetime.now().isoformat(timespec="seconds")


# Matches _R1/_R2 (with optional _001 etc.) or _1/_2 immediately before .fastq.gz
_READ_TAG_RE = re.compile(r'(?:_R([12])(?:_\d+)?|_([12]))\.fastq\.gz$', re.IGNORECASE)


def _strip_read_tag(filename: str):
    """Return (base, read_num) where read_num is '1', '2', or None."""
    m = _READ_TAG_RE.search(filename)
    if m:
        tag = m.group(1) or m.group(2)
        return filename[:m.start()], tag
    return filename[:-len(".fastq.gz")], None


def _list_fastq_pairs(download_dir: Path) -> List[Dict]:
    """Return samples as {sample, paired, r1, r1_name, r2, r2_name} dicts.

    Handles both Illumina (_R1/_R2) and SRA (_1/_2) naming conventions.
    Files with no read suffix are treated as single-end.
    """
    try:
        all_fq = sorted(download_dir.glob("*.fastq.gz"))
    except PermissionError:
        return []
    # A link whose source is gone (moved, deleted, unmounted) is not a sample:
    # its stat() below raised, and one such file failed the whole listing. The
    # Inputs pane still shows it, marked "link broken", so it can be removed.
    def _readable(f: Path) -> bool:
        try:
            f.stat()
            return True
        except OSError:
            return False
    all_fq = [f for f in all_fq if _readable(f)]

    # Preserve insertion order so samples appear sorted
    groups: Dict[str, Dict] = {}
    for fq in all_fq:
        base, tag = _strip_read_tag(fq.name)
        if base not in groups:
            groups[base] = {"r1": None, "r2": None, "extras": []}
        g = groups[base]
        if tag == "1":
            g["r1"] = fq
        elif tag == "2":
            g["r2"] = fq
        else:
            g["extras"].append(fq)

    pairs = []
    for base, g in groups.items():
        r1, r2 = g["r1"], g["r2"]
        if r1 or r2:
            eff_r1 = r1 or r2
            eff_r2 = r2 if r1 else None
            pairs.append({
                "sample": base,
                "paired": bool(r1 and r2),
                "r1": str(eff_r1), "r1_name": eff_r1.name,
                "r1_size": eff_r1.stat().st_size,
                "r2": str(eff_r2) if eff_r2 else None,
                "r2_name": eff_r2.name if eff_r2 else None,
                "r2_size": eff_r2.stat().st_size if eff_r2 else None,
            })
        for fq in g["extras"]:
            pairs.append({
                "sample": fq.name[:-len(".fastq.gz")],
                "paired": False,
                "r1": str(fq), "r1_name": fq.name,
                "r1_size": fq.stat().st_size,
                "r2": None, "r2_name": None,
                "r2_size": None,
            })

    return pairs


def _list_project_samples(project_dir: Path) -> List[Dict]:
    """Enumerate a project's samples from BOTH download/ and step1/.

    download/ is the canonical FASTQ location (SRA/upload land here) and wins on
    dedup. But the vSNP GUI also stages samples as step1/<sample>/ subfolders —
    e.g. folders dragged/placed into a project outside the download flow — whose
    FASTQs never appear in download/. vSNP lists those; without this, the Kraken
    GUI wouldn't, so a sample visible in one tool was invisible in the other.
    step1 FASTQs are symlinks that resolve to real files, so the r1/r2 paths
    returned here are directly runnable via /api/run.
    """
    download_dir = project_dir / "download"
    pairs = _list_fastq_pairs(download_dir) if download_dir.is_dir() else []
    seen = {p["sample"] for p in pairs}
    # An SRA download into a named folder lands in download/<folder>/.
    if download_dir.is_dir():
        for sub in _download_subdirs(download_dir):
            for p in _list_fastq_pairs(sub):
                if p["sample"] in seen:
                    continue
                p["source"] = f"download/{sub.name}"
                seen.add(p["sample"])
                pairs.append(p)

    step1_dir = project_dir / "step1"
    if step1_dir.is_dir():
        for sub in sorted(step1_dir.iterdir()):
            if not sub.is_dir() or sub.name.startswith(("_", ".")):
                continue
            for p in _list_fastq_pairs(sub):
                if p["sample"] in seen:
                    continue
                p["source"] = "step1"
                seen.add(p["sample"])
                pairs.append(p)

    return pairs


# ---------------------------------------------------------------------------
# API routes
# ---------------------------------------------------------------------------

@app.get("/api/projects")
def api_list_projects():
    cfg = load_config()
    projects = _list_projects_from_root(_SHARED_PROJECTS, _SCOPE_SHARED)
    personal_root = Path(cfg.get("projects_root", ""))
    if personal_root != _SHARED_PROJECTS:
        personal = _list_projects_from_root(personal_root, _SCOPE_PERSONAL)
        seen = {p["name"] for p in projects}
        projects += [p for p in personal if p["name"] not in seen]
    return JSONResponse(projects)


class ProjectCreate(BaseModel):
    name: str
    scope: Optional[str] = None   # "personal" (default) | "shared"


@app.post("/api/projects")
def api_create_project(payload: ProjectCreate):
    scope = (payload.scope or _SCOPE_PERSONAL).strip() or _SCOPE_PERSONAL
    if scope not in (_SCOPE_PERSONAL, _SCOPE_SHARED):
        raise HTTPException(400, f"Invalid scope: {scope!r}")
    try:
        project_dir = _create_project(payload.name, scope)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    return JSONResponse({"name": project_dir.name, "path": str(project_dir), "scope": scope})


# ---------------------------------------------------------------------------
# Loading samples into a project — import (link), upload (drag & drop), and
# SRA download. Mirrors vSNP GUI so a project can be populated from within
# either tool. All three land FASTQs in <project>/download/, which the sample
# browser and the run endpoints already read.
# ---------------------------------------------------------------------------
def _writable_project_dir(name: str) -> Path:
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    (project_dir / "download").mkdir(parents=True, exist_ok=True)
    return project_dir


def _download_subdirs(download_dir: Path) -> List[Path]:
    """download/'s visible sub-folders (an SRA download into a named folder)."""
    try:
        return sorted(d for d in download_dir.iterdir() if d.is_dir() and not d.name.startswith((".", "_")))
    except OSError:
        return []


@app.get("/api/projects/{name}/inputs")
def api_project_inputs(name: str):
    """List files in <project>/download/ and its sub-folders (one level).

    `name` is the path under download/ ("sub/file" for a sub-folder), which is
    also what DELETE takes. `dangling` marks a link whose source is gone.
    """
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    download_dir = project_dir / "download"
    files: List[Dict] = []
    total = 0

    def _add(p: Path, rel: str) -> None:
        nonlocal total
        if p.name.startswith("."):
            return
        if p.is_symlink() and not p.exists():
            files.append({"name": rel, "size": 0, "mtime": 0, "dangling": True})
            return
        if not p.is_file():
            return
        try:
            st = p.stat()
        except OSError:
            return
        files.append({"name": rel, "size": st.st_size, "mtime": st.st_mtime})
        total += st.st_size

    if download_dir.is_dir():
        for p in sorted(download_dir.iterdir()):
            _add(p, p.name)
        for sub in _download_subdirs(download_dir):
            try:
                children = sorted(sub.iterdir())
            except OSError:
                continue
            for p in children:
                _add(p, f"{sub.name}/{p.name}")
    return JSONResponse({"files": files, "total_bytes": total, "count": len(files)})


@app.delete("/api/projects/{name}/inputs/{filename:path}")
def api_project_input_delete(name: str, filename: str):
    """Delete one file (or link) from download/, or from one of its sub-folders."""
    parts = [x for x in (filename or "").split("/") if x]
    if (not parts or len(parts) > 2 or "\\" in filename
            or any(x.startswith(".") or x == ".." for x in parts)):
        raise HTTPException(400, "Invalid filename")
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    download_dir = (project_dir / "download")
    target = download_dir.joinpath(*parts)
    # The link itself must sit inside download/ (its source may be anywhere).
    try:
        target.parent.resolve().relative_to(download_dir.resolve())
    except (ValueError, OSError):
        raise HTTPException(400, "Invalid filename")
    if not target.is_file() and not target.is_symlink():
        raise HTTPException(404, f"File not found: {filename}")
    target.unlink()
    return JSONResponse({"deleted": filename})


@app.post("/api/projects/{name}/upload")
async def api_project_upload(name: str, files: List[UploadFile] = File(...)):
    """Save drag-and-dropped / chosen FASTQ files into <project>/download/."""
    project_dir = _writable_project_dir(name)
    download_dir = project_dir / "download"
    saved = 0
    for f in files:
        if not f.filename:
            continue
        target = download_dir / Path(f.filename).name
        async with aiofiles.open(target, "wb") as out:
            while True:
                chunk = await f.read(1024 * 1024)
                if not chunk:
                    break
                await out.write(chunk)
        saved += 1
    return JSONResponse({"uploaded": saved})


class LinkLocalRequest(BaseModel):
    path: Optional[str] = None
    paths: Optional[List[str]] = None


@app.post("/api/projects/{name}/link-local")
def api_project_link_local(name: str, payload: LinkLocalRequest):
    """Symlink server-side FASTQs into download/ — folders and/or single files.

    Lets users bring in reads that already live on the shared filesystem
    without copying gigabytes around. Same contract as the vSNP GUI's
    link-local: a folder contributes every *.fastq.gz directly inside it; one
    bad path does not throw away the rest of a multi-select. Returns what
    happened to each: linked, already there, not found, not a .fastq.gz.
    """
    project_dir = _writable_project_dir(name)
    raw_paths = [str(payload.path or "")] + list(payload.paths or [])
    raw_paths = [r.strip() for r in raw_paths if r and r.strip()]
    if not raw_paths:
        raise HTTPException(400, "No input path given")
    download_dir = project_dir / "download"
    candidates: List[Path] = []
    missing: List[str] = []
    skipped: List[str] = []
    for raw in raw_paths:
        src = Path(raw).expanduser()
        if not src.exists():
            missing.append(str(src))
            continue
        if src.is_file():
            # Only .fastq.gz: the sample list and the pipeline pair reads by
            # that suffix, so a .fq.gz would link and then never show up.
            if src.name.endswith(".fastq.gz"):
                candidates.append(src)
            else:
                skipped.append(src.name)
        else:
            candidates.extend(sorted(src.glob("*.fastq.gz")))
            skipped.extend(sorted(f.name for f in src.glob("*.fq.gz")))
    if missing and not candidates:
        raise HTTPException(400, f"Input path not found: {missing[0]}")
    count = 0
    already: List[str] = []
    for f in sorted(set(candidates)):
        target = download_dir / f.name
        # is_symlink too: a dangling link still occupies the name.
        if target.exists() or target.is_symlink():
            already.append(f.name)
            continue
        target.symlink_to(f.resolve())
        count += 1
    return JSONResponse({
        "linked": count,
        "already_present": already,
        "missing": missing,
        "skipped_not_fastq": skipped,
    })


class SraRequest(BaseModel):
    accessions: List[str]
    folder: Optional[str] = None


@app.post("/api/projects/{name}/sra/download")
def api_project_sra_download(name: str, payload: SraRequest):
    """Resolve SRA accessions and kick off a background download into
    download/. Uses curl/ENA + (if present) fasterq-dump — the script skips
    whatever tool isn't installed, so it works with just curl/wget."""
    project_dir = _writable_project_dir(name)
    # Refuse a mistyped accession with a message that says so (expansion
    # failures are reported as NCBI trouble, which would mislead here).
    bad = [a for a in (payload.accessions or []) if str(a).strip() and not is_valid_accession(str(a))]
    if bad:
        shown = ", ".join(repr(b) for b in bad[:5])
        more = f" (+{len(bad) - 5} more)" if len(bad) > 5 else ""
        raise HTTPException(400, f"Not an SRA/ENA accession: {shown}{more}. "
                                 "Accessions are letters, digits, dot, underscore and hyphen.")
    try:
        expanded, mapping = expand_accessions_with_mapping(payload.accessions, strict=True)
    except SRAExpansionError as e:
        raise HTTPException(
            502,
            f"Could not resolve SRA accessions via NCBI eutils: {e}. "
            "This is usually NCBI rate-limiting; wait ~30 s and retry.",
        )
    download_root = project_dir / "download"
    if payload.folder:
        download_root = download_root / Path(payload.folder).name
    download_root.mkdir(parents=True, exist_ok=True)
    try:
        write_crosswalk_tsv(download_root, mapping)
    except OSError as e:
        logger.warning("Failed to write sra_crosswalk.tsv: %s", e)
    script = build_download_script(
        download_root,
        expanded,
        allow_insecure_https=False,
        # A sample the project already holds (reads kept in step1/<acc>/ by a
        # vSNP run) is listed here too, so it is skipped, not fetched again.
        step1_dir=project_dir / "step1",
        # Always top-level download/, whatever the sub-folder, so the report
        # endpoint has one place to read.
        report_path=project_dir / "download" / "sra_download_report.tsv",
    )
    script_path = download_root / "download_sra.sh"
    script_path.write_text(script, encoding="utf-8")
    script_path.chmod(0o755)
    env = {"PATH": os.environ.get("PATH", "")}
    job_id = job_manager.start_job(
        name=f"sra_download — {name}",
        command=["bash", str(script_path)],
        cwd=download_root,
        env=env,
    )
    return JSONResponse({"job_id": job_id})


@app.get("/api/projects/{name}/sra-download-report")
def api_project_sra_download_report(name: str):
    """The last SRA download's outcome per accession, in buckets.

    Written by the download script as download/sra_download_report.tsv. Empty
    buckets (not 404) when no download has run, so the UI treats both alike.
    """
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    report = project_dir / "download" / "sra_download_report.tsv"
    buckets: Dict[str, List[str]] = {"downloaded": [], "already_in_step1": [], "failed": []}
    if report.is_file():
        try:
            lines = report.read_text(encoding="utf-8").splitlines()
        except OSError:
            lines = []
        for line in lines:
            if not line or line.startswith("#"):
                continue
            outcome, _, acc = line.partition("\t")
            acc = acc.strip()
            if acc and outcome in buckets:
                buckets[outcome].append(acc)
    return JSONResponse(buckets)


@app.get("/api/projects/{name}/sra-crosswalk")
def api_project_sra_crosswalk(name: str):
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    crosswalk = project_dir / "download" / "sra_crosswalk.tsv"
    if not crosswalk.is_file():
        raise HTTPException(404, "No SRA crosswalk for this project")
    return FileResponse(crosswalk, media_type="text/plain")


@app.get("/api/projects/{name}/samples")
def api_project_samples(name: str):
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    return JSONResponse(_list_project_samples(project_dir))


# ---------------------------------------------------------------------------
# Per-sample Kraken results (decoupled from a single job).
#
# Results are read straight from <project>/kraken/<sample>/ on disk so any
# previously-run sample's outputs can be revisited — not just the last job.
# ---------------------------------------------------------------------------
def _collect_result_files(run_dir: Path, include_all: bool) -> List[Dict]:
    """List result files under a kraken run dir, categorized + sorted."""
    files: List[Dict] = []
    if not run_dir.is_dir():
        return files
    for p in sorted(run_dir.rglob("*")):
        if not p.is_file() or p.name.endswith(".log"):
            continue
        rel = str(p.relative_to(run_dir))
        category = _result_category(rel)
        if not include_all and category is None:
            continue
        stat = p.stat()
        files.append({
            "name": rel,
            "path": str(p),
            "label": _result_label(rel, category),
            "size": stat.st_size,
            # _dedupe_primary_results compares mtimes to keep only the newest
            # file in singleton categories (krona/stats/coverage). Without this
            # key it raised KeyError once a sample had 2+ krona files (reruns),
            # 500-ing the per-sample results endpoint so the GUI showed
            # "No Kraken results yet" even though the files were on disk.
            "mtime": stat.st_mtime,
            "openable": _can_open_inline(rel),
            "category": category,
        })
    if not include_all:
        files = _dedupe_primary_results(files)

    def sort_key(f):
        category = f.get("category")
        if category in _CATEGORY_ORDER:
            return (_CATEGORY_ORDER[category], _TIMESTAMP_RE.sub("", f["name"]))
        return (50, f["name"])

    files.sort(key=sort_key)
    for f in files:
        f.pop("mtime", None)
        if include_all and f.get("category") is None:
            f["label"] = f["name"]
    return files


def _sample_run_status(run_dir: Path) -> str:
    """Status for a sample: 'running' if a live job owns its dir, else 'done'
    if the dir holds output, else 'none'."""
    run_dir_str = str(run_dir)
    for job in job_manager.list_jobs():
        if job.get("cwd") == run_dir_str and job.get("status") == "running":
            return "running"
    try:
        if run_dir.is_dir() and any(p.is_file() for p in run_dir.rglob("*")):
            return "done"
    except PermissionError:
        pass
    return "none"


@app.get("/api/projects/{name}/samples/{sample}/kraken-results")
def api_sample_kraken_results(name: str, sample: str, all: int = Query(0)):
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    run_dir = project_dir / "kraken" / sample
    return JSONResponse({
        "project": name,
        "sample": sample,
        "present": run_dir.is_dir(),
        "status": _sample_run_status(run_dir),
        "run_dir": str(run_dir),
        "files": _collect_result_files(run_dir, bool(all)),
    })


# ---------------------------------------------------------------------------
# Cross-tool visibility — surface vSNP results for a sample.
#
# vSNP GUI and this tool share /srv/kapurlab/projects, so vSNP writes its
# per-sample alignment + VCF into <project>/step1/<sample>/ inside the same
# project dir we browse. These endpoints let a user looking at a sample in
# Kraken ID Parse see (and download) what vSNP produced for it. Read-only.
# ---------------------------------------------------------------------------
def _resolve_vsnp_sample_dir(step1_dir: Path, sample: str) -> Optional[Path]:
    """Resolve a sample name to its vSNP step1 subdirectory.

    Exact match first, then the same prefix fallback vSNP itself uses
    (e.g. sample ``13-1941-6`` matches dir ``13-1941-6_S4_L001``).
    """
    exact = step1_dir / sample
    if exact.is_dir():
        return exact
    try:
        candidates = sorted(
            d for d in step1_dir.iterdir()
            if d.is_dir() and d.name.startswith(f"{sample}_")
        )
    except (OSError, PermissionError):
        return None
    return candidates[0] if candidates else None


def _sample_in_step2_run(run_dir: Path, sample: str) -> bool:
    """True if the sample appears as a leaf in any group FASTA of a step2 run.

    vSNP step2 writes per-group alignment FASTAs whose headers are the sample
    VCF leaf names (e.g. ``>hCoV-19-...-EPI_zc.vcf``). Scanning header lines is
    robust across both the timestamped ``runs/`` layout and the legacy flat
    layout, and only reads ``>`` lines so it stays cheap.
    """
    try:
        fastas = list(run_dir.rglob("*.fasta"))
    except (OSError, PermissionError):
        return False
    for fa in fastas:
        try:
            with fa.open("r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    leaf = line[1:].strip()
                    for suf in ("_zc.vcf.gz", "_zc.vcf"):
                        if leaf.endswith(suf):
                            leaf = leaf[: -len(suf)]
                            break
                    if leaf == sample or leaf.startswith(f"{sample}_"):
                        return True
        except (OSError, PermissionError):
            continue
    return False


def _find_latest_step2_for_sample(project_dir: Path, sample: str) -> Dict:
    """Find the most-recent step2 comparison the sample took part in.

    step2 runs accumulate as timestamped folders and each report covers many
    samples, so rather than dumping every run we surface only the latest one
    that actually contains this sample. Returns {"present": False} if none.
    """
    step2_dir = project_dir / "step2"
    if not step2_dir.is_dir():
        return {"present": False}
    runs_dir = step2_dir / "runs"
    candidates: List = []
    if runs_dir.is_dir():
        # Run ids are timestamped, so reverse-sorted name == newest first.
        try:
            candidates = [(d, d.name) for d in sorted(runs_dir.iterdir(), reverse=True) if d.is_dir()]
        except (OSError, PermissionError):
            candidates = []
    else:
        # Legacy flat layout: groups live directly under step2/.
        candidates = [(step2_dir, "legacy")]
    for run_dir, run_id in candidates:
        if not _sample_in_step2_run(run_dir, sample):
            continue
        html = sorted(run_dir.glob("*.html"), key=_safe_mtime)
        report = html[-1] if html else None
        started_at = None
        meta = run_dir / "run_metadata.json"
        if meta.is_file():
            try:
                started_at = json.loads(meta.read_text(encoding="utf-8")).get("started_at")
            except (json.JSONDecodeError, OSError):
                pass
        try:
            groups = [
                g.name for g in sorted(run_dir.iterdir())
                if g.is_dir() and g.name not in ("vcf_source", "runs", "_provenance")
                and not g.name.startswith(".")
            ]
        except (OSError, PermissionError):
            groups = []
        return {
            "present": True,
            "run_id": run_id,
            "started_at": started_at,
            "report_name": report.name if report else None,
            "report_path": str(report) if report else None,
            "groups": groups,
        }
    return {"present": False}


@app.get("/api/projects/{name}/vsnp/samples/{sample}/files")
def api_vsnp_sample_files(name: str, sample: str):
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    step1_dir = project_dir / "step1"
    sample_dir = _resolve_vsnp_sample_dir(step1_dir, sample) if step1_dir.is_dir() else None
    files: List[Dict] = []
    sample_dir_str = ""
    if sample_dir:
        base = sample_dir.resolve()
        sample_dir_str = str(base)
        for p in sorted(base.rglob("*")):
            if not p.is_file() or p.name.startswith(".~lock"):
                continue
            try:
                rel = p.relative_to(base).as_posix()
                st = p.stat()
            except (OSError, ValueError):
                continue
            files.append({
                "name": p.name,
                "relpath": rel,
                "path": str(p),
                "size": st.st_size,
                "openable": _can_open_inline(p.name),
                "type": p.suffix.lstrip(".").lower() or "file",
            })
    return JSONResponse({
        "project": name,
        "sample": sample,
        "step1_present": bool(sample_dir),
        "step1_dir": sample_dir_str,
        "files": files,
        "step2": _find_latest_step2_for_sample(project_dir, sample),
    })


@app.get("/api/projects/{name}/file")
def api_project_file(name: str, path: str = Query(...), inline: int = 0):
    """Serve a file from anywhere inside a project dir (cross-tool downloads).

    Used to download/open vSNP outputs surfaced in the results panel. The
    path-traversal guard keeps `target` inside the resolved project dir.
    """
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    root = project_dir.resolve()
    target = Path(path).resolve()
    if root != target and root not in target.parents:
        raise HTTPException(403, "Path outside project directory")
    if not target.is_file():
        raise HTTPException(404, f"File not found: {path}")
    media_type = _media_type_for(target.name)
    want_inline = bool(inline) and _can_open_inline(target.name)
    if want_inline and target.name == "report.html":
        healed = _report_with_working_links(target)
        if healed is not None:
            return HTMLResponse(healed, headers={"Cache-Control": "no-cache"})
    disposition = "inline" if want_inline else "attachment"
    headers = {"Content-Disposition": f'{disposition}; filename="{target.name}"'}
    return FileResponse(target, media_type=media_type, headers=headers)


def _report_with_working_links(report: Path) -> Optional[str]:
    """A report.html written before its Krona link worked when served, rendered
    again (in memory, nothing written) from the run's manifest; None otherwise.

    v0.2.29 reports linked the Krona chart by a path relative to the run folder,
    which this endpoint's ?path= URL turned into "Not Found". Reports written
    since carry the fix (a data-run-file link and the script that rewrites
    it) and are served exactly as written.
    """
    manifest_path = report.parent / "run_manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        if b"data-run-file" in report.read_bytes():
            return None
        if str(_BIN_DIR) not in sys.path:
            sys.path.insert(0, str(_BIN_DIR))
        from reporting.html_renderer import render_html_text
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return render_html_text(manifest, asset_mode="embed", base_dir=report.parent)
    except Exception as exc:
        logger.warning("Could not re-render %s: %s", report, exc)
        return None


def _home_short(p: Path) -> str:
    """A path as people read it: ~/… when it is under the home folder."""
    try:
        return "~/" + str(Path(p).relative_to(Path.home()))
    except ValueError:
        return str(p)


def _taxa_payload(**extra) -> Dict[str, Any]:
    entries = target_names.entries(_REPO_ROOT)
    return {
        "taxa": [e["name"] for e in entries],
        "entries": entries,
        "saved_file": _home_short(target_names.saved_path()),
        **extra,
    }


@app.get("/api/taxa")
def api_get_taxa():
    """Every target name: the shipped list, then the names users saved."""
    return JSONResponse(_taxa_payload())


class TaxonPayload(BaseModel):
    name: str


@app.post("/api/taxa")
def api_add_taxon(payload: TaxonPayload):
    """Save a new target name (no duplicates, case-insensitive).

    Returns the whole list so the caller refreshes in one round trip.
    """
    try:
        res = target_names.add(_REPO_ROOT, payload.name)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    except OSError as exc:
        raise HTTPException(500, f"Could not save to {target_names.saved_path()}: {exc}")
    return JSONResponse(_taxa_payload(name=res["name"], added=res["added"]))


@app.post("/api/taxa/remove")
def api_remove_taxon(payload: TaxonPayload):
    """Remove a saved target name. Names that ship with the tool stay."""
    try:
        res = target_names.remove(_REPO_ROOT, payload.name)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    except OSError as exc:
        raise HTTPException(500, f"Could not save to {target_names.saved_path()}: {exc}")
    return JSONResponse(_taxa_payload(removed=res["removed"]))


def _project_report_names(project_dir: Path, min_percent: float) -> List[Dict[str, Any]]:
    """Taxon names from every Kraken report in a project, best evidence first.

    One entry per name with the highest clade % seen, its rank, and how many
    samples reported it. These are the names exactly as Kraken spells them —
    the spelling the read parser needs.
    """
    best: Dict[str, Dict[str, Any]] = {}
    kdir = project_dir / "kraken"
    if not kdir.is_dir():
        return []
    try:
        sample_dirs = [d for d in kdir.iterdir() if d.is_dir() and not d.name.startswith(".")]
    except OSError:
        return []
    for sdir in sample_dirs:
        reports = sorted((sdir / "kraken").glob("*_reportkraken.txt")) if (sdir / "kraken").is_dir() else []
        if not reports:
            continue
        report = max(reports, key=lambda r: r.stat().st_mtime)
        seen_here = set()
        for row in target_names.report_names(report, min_percent):
            key = row["name"].lower()
            if key in seen_here:
                continue
            seen_here.add(key)
            cur = best.get(key)
            if cur is None:
                best[key] = {**row, "samples": 1, "best_sample": sdir.name}
            else:
                cur["samples"] += 1
                if row["percent"] > cur["percent"]:
                    cur.update(percent=row["percent"], reads=row["reads"], best_sample=sdir.name)
    return sorted(best.values(), key=lambda r: (-r["percent"], r["name"].lower()))


@app.get("/api/projects/{name}/report-taxa")
def api_project_report_taxa(name: str, min_percent: float = Query(0.01)):
    """Names Kraken reported in this project's runs, for picking a target name."""
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    rows = _project_report_names(project_dir, max(0.0, min_percent))
    return JSONResponse({"names": rows[:2000], "total": len(rows)})


@app.get("/api/config")
def api_get_config():
    cfg = load_config()
    # Deployed checkout's version (git describe) — what the dashboard shows.
    # Injected into the response only, never written back to the config file
    # (and POST /api/config's ConfigPayload has no such field, so a round-trip
    # can't persist it either).
    cfg["app_version"] = APP_VERSION
    return JSONResponse(cfg)


def _db_entry(p: str) -> Dict[str, Any]:
    """One dropdown entry for a database path: name, size (sum of *.k2d) and a
    ``missing`` flag when the path no longer holds a Kraken2 DB (moved DB,
    unmounted share) — shown rather than hidden so the user can remove it."""
    d = Path(p)
    entry: Dict[str, Any] = {"name": d.name, "path": p, "size_bytes": 0}
    try:
        if (d / "hash.k2d").is_file():
            entry["size_bytes"] = sum(f.stat().st_size for f in d.glob("*.k2d"))
        else:
            entry["missing"] = True
    except OSError:
        entry["missing"] = True
    return entry


@app.get("/api/kraken-dbs")
def api_kraken_dbs():
    """Kraken2 databases for the settings/run dropdowns: exactly the paths the
    user has added (``saved_kraken_dbs``) plus the currently active one.

    This list is managed explicitly in Settings — add and remove, like the
    saved project locations. Nothing is discovered or auto-included: a sibling
    database sitting next to a configured one stays out of the dropdowns until
    the user adds it. (Earlier versions scanned the parent directories of every
    known DB and the site database root, which dumped every database on the
    system into the list.)
    """
    cfg = load_config()
    current = (cfg.get("kraken_db") or "").strip()
    saved = [str(p).strip() for p in (cfg.get("saved_kraken_dbs") or []) if str(p).strip()]
    dbs: List[Dict[str, Any]] = []
    seen: set = set()
    for p in (*saved, current):
        if not p:
            continue
        try:
            rp = str(Path(p).resolve())
        except OSError:
            rp = p
        if rp in seen:
            continue
        seen.add(rp)
        dbs.append(_db_entry(p))
    return JSONResponse({"databases": dbs, "current": current, "saved": saved})


class DbPathPayload(BaseModel):
    path: str


@app.post("/api/kraken-dbs/add")
def api_kraken_dbs_add(payload: DbPathPayload):
    """Explicitly add one database location to the saved list. Validates the
    path up front (same files kraken2 itself needs) so a typo is a clear error
    now instead of a failed run later."""
    p = (payload.path or "").strip().rstrip("/")
    if not p:
        raise HTTPException(400, "A database path is required.")
    d = Path(p)
    if not d.is_dir():
        raise HTTPException(400, f"Not a directory: {p}")
    missing = [n for n in ("hash.k2d", "opts.k2d", "taxo.k2d") if not (d / n).is_file()]
    if missing:
        raise HTTPException(
            400,
            f"Not a Kraken2 database (missing {', '.join(missing)}). "
            "Point at the directory that contains hash.k2d, opts.k2d and taxo.k2d.",
        )
    cfg = load_config()
    saved = [str(x).strip() for x in (cfg.get("saved_kraken_dbs") or []) if str(x).strip()]
    if p not in saved:
        saved.append(p)
        cfg["saved_kraken_dbs"] = saved
        save_config(cfg)
    return JSONResponse({"ok": True, "saved": saved})


@app.post("/api/kraken-dbs/remove")
def api_kraken_dbs_remove(payload: DbPathPayload):
    """Remove one database location from the saved list (nothing on disk is
    touched). If it was the active database, the active choice is cleared too —
    a removed database shouldn't stay silently selected."""
    p = (payload.path or "").strip()
    cfg = load_config()
    saved = [str(x).strip() for x in (cfg.get("saved_kraken_dbs") or []) if str(x).strip()]
    cfg["saved_kraken_dbs"] = [x for x in saved if x != p]
    cleared = False
    if (cfg.get("kraken_db") or "").strip() == p:
        cfg["kraken_db"] = ""
        cleared = True
    save_config(cfg)
    return JSONResponse({"ok": True, "saved": cfg["saved_kraken_dbs"], "cleared_active": cleared})


def _blast_db_entry(p: str) -> Dict[str, Any]:
    """One dropdown entry for a BLAST database.

    A BLAST database is a file PREFIX, not a directory: `/db/blast/nt` is
    backed by nt.nhr / nt.nin / nt.nsq (or nt.00.nhr … plus an nt.nal alias
    for a volume set). A bare name with no path separator (`nt`, `nr`) is
    NCBI's remote database, which needs nothing on disk.
    """
    entry: Dict[str, Any] = {"name": Path(p).name, "path": p, "size_bytes": 0,
                             "remote": "/" not in p}
    if entry["remote"]:
        return entry
    d = Path(p)
    try:
        vols = [f for f in d.parent.glob(d.name + ".*")
                if f.is_file() and re.search(r"\.(\d+\.)?[np][a-z]{2}$", f.name)]
        if vols:
            entry["size_bytes"] = sum(f.stat().st_size for f in vols)
        else:
            entry["missing"] = True
    except OSError:
        entry["missing"] = True
    return entry


@app.get("/api/blast-dbs")
def api_blast_dbs():
    """BLAST databases for the settings/run dropdowns: the paths the user has
    added (``saved_blast_dbs``) plus the active one. Managed explicitly, like
    the Kraken2 list — nothing is discovered on disk."""
    cfg = load_config()
    current = (cfg.get("blast_db") or "").strip()
    saved = [str(p).strip() for p in (cfg.get("saved_blast_dbs") or []) if str(p).strip()]
    dbs: List[Dict[str, Any]] = []
    seen: set = set()
    for p in (*saved, current):
        if not p or p in seen:
            continue
        seen.add(p)
        dbs.append(_blast_db_entry(p))
    return JSONResponse({"databases": dbs, "current": current, "saved": saved})


def _validate_blast_db(p: str) -> None:
    """Raise HTTPException unless `p` looks like a usable BLAST database."""
    if "/" not in p:
        return  # bare name — NCBI remote (nt, nr); nothing to check on disk
    d = Path(p)
    if d.is_dir():
        raise HTTPException(
            400,
            f"{p} is a directory. A BLAST database is a file PREFIX — include the "
            "database name itself (e.g. .../blast/ref_prok_rep_genomes).",
        )
    try:
        vols = [f for f in d.parent.glob(d.name + ".*")
                if f.is_file() and re.search(r"\.(\d+\.)?[np][a-z]{2}$", f.name)]
    except OSError:
        vols = []
    if not vols:
        raise HTTPException(
            400,
            f"No BLAST database files found for prefix {p} (expected {d.name}.nhr / "
            f".nin / .nsq, or a {d.name}.nal volume alias beside it).",
        )


@app.post("/api/blast-dbs/add")
def api_blast_dbs_add(payload: DbPathPayload):
    """Add one BLAST database to the saved list, validating it up front so a
    typo is an error here rather than a failed run an hour into a pipeline."""
    p = (payload.path or "").strip().rstrip("/")
    if not p:
        raise HTTPException(400, "A database path or name is required.")
    _validate_blast_db(p)
    cfg = load_config()
    saved = [str(x).strip() for x in (cfg.get("saved_blast_dbs") or []) if str(x).strip()]
    if p not in saved:
        saved.append(p)
        cfg["saved_blast_dbs"] = saved
        save_config(cfg)
    return JSONResponse({"ok": True, "saved": saved})


@app.post("/api/blast-dbs/remove")
def api_blast_dbs_remove(payload: DbPathPayload):
    """Remove one BLAST database from the saved list (nothing on disk is
    touched); clear the active choice if that is what was removed."""
    p = (payload.path or "").strip()
    cfg = load_config()
    saved = [str(x).strip() for x in (cfg.get("saved_blast_dbs") or []) if str(x).strip()]
    cfg["saved_blast_dbs"] = [x for x in saved if x != p]
    cleared = False
    if (cfg.get("blast_db") or "").strip() == p:
        cfg["blast_db"] = ""
        cleared = True
    save_config(cfg)
    return JSONResponse({"ok": True, "saved": cfg["saved_blast_dbs"], "cleared_active": cleared})


def _remember_blast_db(cfg: Dict[str, Any], db_path: str) -> bool:
    """Remember a BLAST database actually run with, so it is a dropdown pick
    next time. Returns True when the config changed."""
    db_path = (db_path or "").strip()
    if not db_path:
        return False
    if "/" in db_path:
        d = Path(db_path)
        try:
            if not any(f.is_file() and re.search(r"\.(\d+\.)?[np][a-z]{2}$", f.name)
                       for f in d.parent.glob(d.name + ".*")):
                return False
        except OSError:
            return False
    saved = [str(p).strip() for p in (cfg.get("saved_blast_dbs") or []) if str(p).strip()]
    if db_path in saved:
        return False
    saved.append(db_path)
    cfg["saved_blast_dbs"] = saved
    return True


def _remember_kraken_db(cfg: Dict[str, Any], db_path: str) -> bool:
    """Add a Kraken DB path to the remembered list (saved_kraken_dbs) so
    switching databases later is a pick, not a re-type. Returns True if the
    config changed. Only real DBs are remembered — a typo'd path should fail
    its run, not haunt the dropdown forever."""
    db_path = (db_path or "").strip()
    if not db_path or not (Path(db_path) / "hash.k2d").is_file():
        return False
    saved = [str(p).strip() for p in (cfg.get("saved_kraken_dbs") or []) if str(p).strip()]
    if db_path in saved:
        return False
    saved.append(db_path)
    cfg["saved_kraken_dbs"] = saved
    return True


class ConfigPayload(BaseModel):
    kraken_db: Optional[str] = None
    blast_db: Optional[str] = None
    projects_root: Optional[str] = None
    shared_projects_root: Optional[str] = None
    saved_project_roots: Optional[List[str]] = None
    saved_kraken_dbs: Optional[List[str]] = None
    saved_blast_dbs: Optional[List[str]] = None


@app.post("/api/config")
def api_save_config(payload: ConfigPayload):
    cfg = load_config()
    updates = payload.model_dump(exclude_none=True)
    cfg.update(updates)
    # Normalize the curated Projects-root bookmarks (dedupe, drop blanks).
    roots = cfg.get("saved_project_roots") or []
    if isinstance(roots, list):
        seen, cleaned = set(), []
        for r in roots:
            r = (r or "").strip()
            if r and r not in seen:
                seen.add(r); cleaned.append(r)
        cfg["saved_project_roots"] = cleaned
    # Same for the remembered Kraken DBs, and auto-remember the active one so
    # every DB ever configured here stays one dropdown pick away.
    dbs = cfg.get("saved_kraken_dbs")
    if isinstance(dbs, list):
        seen, cleaned = set(), []
        for r in dbs:
            r = (str(r) or "").strip()
            if r and r not in seen:
                seen.add(r); cleaned.append(r)
        cfg["saved_kraken_dbs"] = cleaned
    _remember_kraken_db(cfg, cfg.get("kraken_db", ""))
    save_config(cfg)
    return JSONResponse({"ok": True})


def _disk_usage_for(p: Path) -> Optional[Dict[str, int]]:
    try:
        u = shutil.disk_usage(p)
    except OSError:
        return None
    return {"total": u.total, "free": u.free, "used": u.used}


@app.get("/api/browse-dirs")
def api_browse_dirs(path: str = "", include_files: bool = False, exts: str = ""):
    """List sub-directories of `path` for the project-root folder picker.

    Runs as the OOD session user, so the OS filesystem permissions are the only
    limit on what can be browsed. Defaults to the user's home when no path is
    given. Returns the resolved path, its parent (null at the filesystem root),
    and the immediate sub-directories.
    """
    try:
        p = (Path(path).expanduser() if path.strip() else Path.home()).resolve()
    except (OSError, RuntimeError):
        raise HTTPException(400, "Invalid path")
    if not p.is_dir():
        raise HTTPException(400, f"Not a directory: {p}")
    # include_files=1 (the "Browse server files" picker) also lists the
    # folder's files, narrowed to ?exts= suffixes, so reads already on the
    # server are linked in place instead of uploaded.
    suffixes = tuple(x.strip().lower() for x in exts.split(",") if x.strip())
    entries: List[Dict[str, str]] = []
    files: List[Dict[str, Any]] = []
    try:
        for child in sorted(p.iterdir(), key=lambda c: c.name.lower()):
            if child.name.startswith("."):
                continue
            try:
                if child.is_dir():
                    entries.append({"name": child.name, "path": str(child)})
                elif include_files and child.is_file():
                    if suffixes and not child.name.lower().endswith(suffixes):
                        continue
                    try:
                        size = child.stat().st_size
                    except OSError:
                        size = 0
                    files.append({"name": child.name, "path": str(child), "size": size})
            except OSError:
                continue
    except PermissionError:
        raise HTTPException(403, f"Permission denied: {p}")
    parent = str(p.parent) if p.parent != p else None
    return JSONResponse({
        "path": str(p), "parent": parent, "entries": entries, "files": files,
        "disk": _disk_usage_for(p),
        "writable": os.access(p, os.W_OK | os.X_OK),
    })


class RunPayload(BaseModel):
    project: str
    r1: str           # absolute path to R1 FASTQ
    r2: Optional[str] = None
    taxon: Optional[str] = None
    kraken_db: Optional[str] = None
    blast_db: Optional[str] = None
    kraken_only: bool = False   # run only Kraken2 + Krona graph; skip parsing/assembly/BLAST
    no_blast: bool = False       # run Kraken2 + read parsing only; skip assembly/BLAST/coverage


@app.post("/api/run")
def api_run(payload: RunPayload):
    cfg = load_config()
    kraken_db = payload.kraken_db or cfg.get("kraken_db", "")
    blast_db = payload.blast_db or cfg.get("blast_db", "nt")
    taxon = (payload.taxon or "").strip()
    # Remember any DB actually run with, so it's a dropdown pick next time —
    # including ad-hoc per-run overrides typed into the run form.
    remembered = _remember_kraken_db(cfg, kraken_db)
    # Only full runs actually BLAST; don't remember a DB a skipped step never used.
    if not payload.kraken_only and not payload.no_blast:
        remembered = _remember_blast_db(cfg, blast_db) or remembered
    if remembered:
        save_config(cfg)

    # Taxon drives the read-parsing/assembly/BLAST identification, which
    # kraken-only mode skips — so it's only required for runs that parse.
    if not payload.kraken_only and not taxon:
        raise HTTPException(400, "A target taxon is required for a full run.")
    if payload.kraken_only and not kraken_db:
        raise HTTPException(400, "Kraken-only mode requires a Kraken DB path.")
    if payload.no_blast and not kraken_db:
        raise HTTPException(400, "Parse-only (skip BLAST) mode requires a Kraken DB path.")

    # Validate paths
    r1 = Path(payload.r1)
    if not r1.exists():
        raise HTTPException(400, f"R1 file not found: {payload.r1}")

    project_dir = _get_project_dir(payload.project)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {payload.project}")

    # Derive sample name — strip _R1/_R2 or _1/_2 read tags
    sample_name, _ = _strip_read_tag(r1.name)

    # Output directory: <project>/kraken/<sample_name>/
    run_dir = project_dir / "kraken" / sample_name

    # Refuse to start a second pipeline in the same output directory: two runs
    # sharing a working directory race on the same temp/output folders and can
    # delete each other's CWD mid-run (manifests as "getcwd() failed" and
    # spurious SPAdes/Excel failures).
    for existing in job_manager.list_jobs():
        if existing.get("status") == "running" and existing.get("cwd") == str(run_dir):
            raise HTTPException(
                409,
                f"A run is already in progress for {sample_name} "
                f"(job {existing['id'][:8]}). Wait for it to finish before re-running.",
            )

    run_dir.mkdir(parents=True, exist_ok=True)

    # Build argv without a shell so project/sample/taxon/db values cannot be
    # interpreted as shell syntax. A Kraken-only run goes through a GUI-side
    # wrapper that runs the same pipeline and then writes report.html/.pdf —
    # the pipeline itself stops at the Krona chart in that mode and writes no
    # report (bin/reporting/kraken_only_report.py).
    script = (_BIN_DIR / "reporting" / "kraken_only_report.py") if payload.kraken_only \
        else (_BIN_DIR / "kraken_id_parse.py")
    # -u: unbuffered stdout/stderr so the parent script's progress prints stream
    # to the log in real time and in order (otherwise Python block-buffers when
    # stdout is a pipe and the log appears to freeze between subprocess outputs).
    command = [sys.executable, "-u", str(script), "-r1", str(r1)]
    if payload.r2:
        r2 = Path(payload.r2)
        if r2.exists():
            command.extend(["-r2", str(r2)])
        else:
            raise HTTPException(400, f"R2 file not found: {payload.r2}")
    if taxon:
        command.extend(["-t", taxon])
    if kraken_db:
        command.extend(["-k", kraken_db])
    if payload.kraken_only:
        pass  # the wrapper adds --kraken-only itself
    elif payload.no_blast:
        command.append("--no-blast")
    elif blast_db:
        command.extend(["-b", blast_db])

    env = {
        "PYTHONPATH": str(_BIN_DIR),
        "PATH": os.environ.get("PATH", ""),
        "PYTHONUNBUFFERED": "1",  # belt-and-suspenders with `python -u`
    }

    if payload.kraken_only:
        job_label = "Kraken-only (Krona)"
    elif payload.no_blast:
        job_label = f"{taxon} (parse only)" if taxon else "parse only"
    else:
        job_label = taxon or "run"
    job_name = f"{payload.project}/{sample_name} — {job_label}"
    job_id = job_manager.start_job(name=job_name, command=command, cwd=run_dir, env=env)
    return JSONResponse({"job_id": job_id, "run_dir": str(run_dir)})


@app.get("/api/jobs")
def api_list_jobs():
    return JSONResponse(job_manager.list_jobs())


@app.get("/api/jobs/{job_id}")
def api_get_job(job_id: str):
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")
    return JSONResponse(job)


@app.get("/api/jobs/{job_id}/log")
async def api_job_log(job_id: str, request: Request):
    """SSE stream of the job's log file. Tails from beginning, closes when job finishes."""
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")

    log_path = Path(job["log_path"])

    _ansi_re = re.compile(r'\x1b\[[0-9;]*[mGKHFABCDJsur]')

    async def event_stream():
        position = 0
        while True:
            if await request.is_disconnected():
                break
            current_job = job_manager.get_job(job_id)
            if log_path.exists():
                async with aiofiles.open(log_path, "r", encoding="utf-8", errors="replace") as f:
                    await f.seek(position)
                    chunk = await f.read(4096)
                    if chunk:
                        lines = chunk.splitlines(keepends=True)
                        for line in lines:
                            clean = _ansi_re.sub("", line.rstrip())
                            if clean:
                                yield f"data: {clean}\n\n"
                        position += len(chunk.encode("utf-8"))
            if current_job and current_job["status"] in ("succeeded", "failed"):
                yield "data: [DONE]\n\n"
                break
            await asyncio.sleep(0.5)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


# ---------------------------------------------------------------------------
# Plain (non-streaming) log + status in ONE GET — polled by the UI.
#
# The SSE route above stays for anything still listening, but the UI no longer
# uses it. Behind an Open OnDemand /rnode Apache reverse proxy a held-open SSE
# connection is hazardous: the proxy hands the SSE buffer back as the body of
# concurrent sibling GETs, so a status poll arrives as log text, JSON parsing
# fails, and a run that succeeded is reported "failed". The same proxy truncates
# buffered responses at roughly 43.5 KB. A plain GET returning BOTH the recorded
# status (from the real exit code) and the current log text, tail-truncated well
# under that ceiling, drives status and live-ish logs from one proxy-safe poll
# loop. Pattern from mhc_gui, which polled from the start and never had the bug;
# the 30000-char tail is the value field-tested behind that proxy (ICAR-NIVEDI).
# ---------------------------------------------------------------------------
_LOGTEXT_MAX_CHARS = 30000
_LOGTEXT_ANSI_RE = re.compile(r'\x1b\[[0-9;]*[mGKHFABCDJsur]')


@app.get("/api/jobs/{job_id}/logtext")
def api_job_logtext(job_id: str):
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")
    text = ""
    lp = job.get("log_path")
    try:
        if lp and Path(lp).is_file():
            text = Path(lp).read_text(encoding="utf-8", errors="replace")
    except OSError:
        pass
    text = _LOGTEXT_ANSI_RE.sub("", text)
    truncated = len(text) > _LOGTEXT_MAX_CHARS
    if truncated:
        text = "...(earlier log truncated)...\n" + text[-_LOGTEXT_MAX_CHARS:]
    return JSONResponse({
        "status": job.get("status"),
        "exit_code": job.get("exit_code"),
        "log": text,
        "truncated": truncated,
    })


# File extensions a browser can render in a tab (open inline); everything else
# is sent as a download. Maps extension -> MIME type.
_INLINE_MEDIA = {
    ".pdf": "application/pdf",
    ".html": "text/html",
    ".htm": "text/html",
    ".txt": "text/plain",
    ".log": "text/plain",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".svg": "image/svg+xml",
    ".csv": "text/plain",
}
_DOWNLOAD_MEDIA = {
    ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    ".xls": "application/vnd.ms-excel",
    ".vcf": "text/plain",
    ".fasta": "text/plain",
    ".fa": "text/plain",
    ".gz": "application/gzip",
}
_TIMESTAMP_RE = re.compile(r"_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}")


def _can_open_inline(name: str) -> bool:
    return Path(name).suffix.lower() in _INLINE_MEDIA


def _media_type_for(name: str) -> str:
    ext = Path(name).suffix.lower()
    return _INLINE_MEDIA.get(ext) or _DOWNLOAD_MEDIA.get(ext) or "application/octet-stream"


def _result_category(rel: str) -> Optional[str]:
    """Return the primary-results category for a relative run output path.

    The run directory intentionally keeps intermediate files for audit/debugging,
    but the GUI should show the small set users normally open or download.
    """
    path = Path(rel)
    name = path.name
    parts = path.parts

    if any(part.startswith(".") for part in parts):
        return None
    if parts and parts[0] == "previews":
        return None
    if name.endswith("-banner.png") or name == "bracken_pie.png":
        return None
    if name.endswith(".fastq.gz"):
        return None
    if name in {"run_manifest.json", "combined_references.fasta", "CAUTION_SITES.xlsx"}:
        return None

    if rel == "report.pdf":
        return "report_pdf"
    if rel == "report.html":
        return "report_html"
    # The pipeline's own report, <sample>_<stamp>_report.html/.pdf. It is the
    # one that carries the organism-specific sections — BTV serotyping, segment
    # status, tentative segments — so it is listed, not left under "all files".
    if len(parts) == 1 and name.endswith("_report.html"):
        return "pipeline_report_html"
    if len(parts) == 1 and name.endswith("_report.pdf"):
        return "pipeline_report_pdf"
    if name.endswith("_stats.xlsx"):
        return "stats"
    if name.endswith("_krona.html"):
        return "krona"
    if name.endswith("-bracken.xlsx"):
        return "bracken"
    if name.endswith("-coverage_interactive.html"):
        return "coverage_interactive"
    if name.endswith("_blast_summary.txt") or name == "consensus_blast_summary.txt":
        return "blast_summary"
    if name.endswith("_denovo.fasta"):
        return "denovo_fasta"
    if name.endswith("_reference_guided.fasta"):
        return "reference_fasta"
    if name.endswith("-coverage_graph.pdf"):
        return "coverage_pdf"
    if name.endswith("-coverage_stats.json"):
        return None
    if parts and parts[0] == "kraken":
        return None
    return None


_CATEGORY_ORDER = {
    # HTML first: it is the primary report (it carries the interactive
    # Coverage & Variants charts); the PDF is its print rendering.
    "report_html": 0,
    "report_pdf": 1,
    "pipeline_report_html": 1.3,
    "pipeline_report_pdf": 1.6,
    "stats": 2,
    "krona": 3,
    "blast_summary": 4,
    "denovo_fasta": 5,
    "reference_fasta": 6,
    "coverage_interactive": 6.5,
    "coverage_pdf": 7,
    "bracken": 8,
    "log": 99,
}


def _result_label(rel: str, category: Optional[str]) -> str:
    if category == "coverage_interactive":
        return "Coverage & Variants (interactive)"
    if category == "report_pdf":
        return "Report PDF"
    if category == "report_html":
        return "Report HTML"
    if category == "pipeline_report_html":
        return "Pipeline report HTML"
    if category == "pipeline_report_pdf":
        return "Pipeline report PDF"
    if category == "stats":
        return "Run statistics workbook"
    if category == "krona":
        return "Krona taxonomy report"
    if category == "blast_summary":
        return "BLAST summary" if Path(rel).name != "consensus_blast_summary.txt" else "Consensus BLAST summary"
    if category == "denovo_fasta":
        return "De novo assembly FASTA"
    if category == "reference_fasta":
        return "Reference-guided FASTA"
    if category == "coverage_pdf":
        return "Coverage graph PDF"
    if category == "bracken":
        return "Bracken results workbook"
    if category == "log":
        return "Pipeline log"
    return rel


def _dedupe_primary_results(files: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Keep only latest timestamped/rerun artifacts in each primary category."""
    latest_by_category: Dict[str, Dict[str, Any]] = {}
    passthrough: List[Dict[str, Any]] = []

    singleton_categories = {"stats", "krona", "coverage_pdf",
                            "pipeline_report_html", "pipeline_report_pdf"}
    for file in files:
        category = file.get("category")
        if category in singleton_categories:
            current = latest_by_category.get(category)
            if current is None or file["mtime"] > current["mtime"]:
                latest_by_category[category] = file
        else:
            passthrough.append(file)

    return passthrough + list(latest_by_category.values())


def _render_api_report_html(run_dir: Path, report_html: Path) -> Optional[str]:
    """Render report.html with API-relative asset URLs for GUI viewing."""
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        if str(_BIN_DIR) not in sys.path:
            sys.path.insert(0, str(_BIN_DIR))
        from reporting.html_renderer import render_html_text

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        # base_dir: the run folder as it is now, so a run whose folder has
        # moved since still finds its Kraken report and Krona chart.
        return render_html_text(manifest, asset_mode="api", base_dir=run_dir)
    except Exception as exc:
        logger.warning("Could not render API report HTML from %s: %s", manifest_path, exc)
        try:
            return report_html.read_text(encoding="utf-8")
        except OSError:
            return None


@app.get("/api/jobs/{job_id}/results")
def api_job_results(job_id: str, all: int = Query(0)):
    """List output files in the job's run directory, plus the pipeline log.

    Each entry carries `openable` (can the browser render it inline) so the UI
    can show an Open link; everything is downloadable via /file.
    """
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")

    files = []
    cwd = job.get("cwd")
    if cwd and Path(cwd).is_dir():
        run_dir = Path(cwd)
        for p in sorted(run_dir.rglob("*")):
            if p.is_file() and not p.name.endswith(".log"):
                rel = str(p.relative_to(run_dir))
                category = _result_category(rel)
                if not all and category is None:
                    continue
                files.append(
                    {
                        "name": rel,
                        "label": _result_label(rel, category),
                        "size": p.stat().st_size,
                        "mtime": p.stat().st_mtime,
                        "openable": _can_open_inline(rel),
                        "category": category,
                    }
                )

    # Surface the pipeline log itself (lives outside run_dir, in the jobs dir).
    log_path = Path(job.get("log_path", ""))
    if log_path.is_file():
        files.append(
            {
                "name": "pipeline_log.txt",
                "label": "Pipeline log",
                "size": log_path.stat().st_size,
                "mtime": log_path.stat().st_mtime,
                "openable": True,
                "category": "log",
                "is_log": True,
            }
        )

    if not all:
        files = _dedupe_primary_results(files)

    # Sort: primary outputs in workflow order, log last.
    def sort_key(f):
        if f.get("is_log"):
            return (_CATEGORY_ORDER["log"], f["name"])
        category = f.get("category")
        if category in _CATEGORY_ORDER:
            return (_CATEGORY_ORDER[category], _TIMESTAMP_RE.sub("", f["name"]))
        return (50, f["name"])

    files.sort(key=sort_key)
    for file in files:
        file.pop("mtime", None)
        if all and file.get("category") is None:
            file["label"] = file["name"]
    return JSONResponse(files)


@app.get("/api/jobs/{job_id}/file")
def api_job_file(job_id: str, path: str = Query(...), inline: int = 0):
    """Serve a single result file. `inline=1` renders in the browser (PDF/HTML/
    text/images); otherwise it downloads as an attachment.

    `path` is relative to the job's run directory, except the sentinel
    "pipeline_log.txt" which maps to the job's log file.
    """
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")

    if path == "pipeline_log.txt":
        target = Path(job.get("log_path", ""))
        display_name = f"{job_id[:8]}_pipeline_log.txt"
    else:
        cwd = job.get("cwd")
        if not cwd:
            raise HTTPException(404, "No run directory for job")
        run_dir = Path(cwd).resolve()
        target = (run_dir / path).resolve()
        # Path-traversal guard: target must stay inside the run directory.
        if run_dir != target and run_dir not in target.parents:
            raise HTTPException(403, "Path outside run directory")
        display_name = target.name

    if not target.is_file():
        raise HTTPException(404, f"File not found: {path}")

    media_type = _media_type_for(target.name)
    want_inline = bool(inline) and _can_open_inline(target.name)
    if path == "report.html" and want_inline and target.name == "report.html":
        html = _render_api_report_html(run_dir, target) if path != "pipeline_log.txt" else None
        if html is not None:
            return HTMLResponse(
                html,
                headers={"Cache-Control": "no-cache, no-store, must-revalidate"},
            )

    disposition = "inline" if want_inline else "attachment"
    headers = {"Content-Disposition": f'{disposition}; filename="{display_name}"'}
    return FileResponse(target, media_type=media_type, headers=headers)


# ---------------------------------------------------------------------------
# Static frontend — must be last (catches everything not matched above)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# The Results pane — every completed sample in one searchable table.
#
# Modelled on vSNP's Step 1 Results. The per-sample endpoints elsewhere in
# this file answer "tell me about THIS one", which is why results were only ever
# visible by expanding a row in the Projects tree; there was no way to see, sort,
# search or export everything a project had produced. This answers that.
# ---------------------------------------------------------------------------
_RP_SUBDIR = "kraken"

def _rp_finished_at(run_dir: Path) -> str:
    """When this run finished, as an ISO string ("" if unknown).

    Prefer the pipeline's own record over filesystem mtimes: a later re-read or
    an rsync can touch files long after the analysis actually ran."""
    try:
        data = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        for key in ("pipeline_finished_at", "finished_at_utc", "finished_at", "timestamp"):
            val = str(data.get(key) or "").strip()
            if val:
                return val
    except (OSError, ValueError, AttributeError):
        pass
    try:
        newest = max((p.stat().st_mtime for p in run_dir.rglob("*") if p.is_file()), default=0)
        if newest:
            from datetime import datetime, timezone
            return datetime.fromtimestamp(newest, tz=timezone.utc).isoformat()
    except (OSError, ValueError):
        pass
    return ""


def _rp_flags(run_dir: Path) -> Dict[str, Any]:
    """pass / review / fail plus reasons.

    The run's exit status outranks any QC verdict: qc.json grades the INPUT, so a
    run whose analysis step exited non-zero — producing nothing — would otherwise
    still be reported as "pass". A Results pane that says PASS for a failed run is
    worse than no pane at all."""
    level, reasons = "pass", []
    try:
        qc = json.loads((run_dir / "qc.json").read_text(encoding="utf-8"))
        verdict = str(qc.get("verdict") or "").strip().lower()
        if verdict in ("fail", "failed"):
            level = "fail"
        elif verdict in ("review", "warn", "warning"):
            level = "review"
        notes = qc.get("notes") or qc.get("reasons") or []
        if isinstance(notes, str):
            notes = [notes]
        reasons = [str(n) for n in notes if str(n).strip()]
    except (OSError, ValueError, AttributeError):
        pass
    try:
        man = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        rc = man.get("return_code")
        if rc not in (None, 0):
            level = "fail"
            lines = [l.strip() for l in str(man.get("stderr_tail") or "").splitlines() if l.strip()]
            detail = ""
            for i, l in enumerate(lines):
                if l.startswith("*** ERROR"):
                    detail = next((x for x in lines[i + 1:]
                                   if not x.startswith(("Command line:", "Running:"))), "")
                    break
            if not detail:
                detail = next((l for l in lines if "error" in l.lower()), "")
            if not detail:
                detail = next((l for l in reversed(lines)
                               if not l.startswith(("Command line:", "Running:"))), "")
            reasons.insert(0, "the analysis exited %s%s" % (rc, (" — " + detail[:120]) if detail else ""))
    except (OSError, ValueError, AttributeError):
        pass
    # The Orbivirus (BTV/EHD) analysis is beta — untested, its results not to be
    # trusted — so its runs are never graded PASS. Runs made before the pipeline
    # recorded the notice are recognized by their target taxon.
    try:
        man = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        if str(_BIN_DIR) not in sys.path:
            sys.path.insert(0, str(_BIN_DIR))
        import orbivirus_beta
        if orbivirus_beta.notice_for_manifest(man):
            if level == "pass":
                level = "review"
            reasons.append(orbivirus_beta.CHIP)
    except (OSError, ValueError, AttributeError, ImportError):
        pass
    return {"level": level, "reasons": reasons}


def _rp_json(run_dir: Path, name: str) -> Dict[str, Any]:
    """Read one of the run's small result JSONs; {} when absent or malformed."""
    try:
        data = json.loads((run_dir / name).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


_rp_status = _sample_run_status


def _rp_kraken_top(report: Path) -> List[tuple]:
    """Species-level calls from a Kraken report, best first, one per lineage.

    Candidates are species (S), the ranks below species (S1, S2…) and the
    species groups between genus and species (G1, G2… — e.g. "Mycobacterium
    tuberculosis complex"). Kraken often assigns most reads to such a group
    rather than to one species, so looking at S rows alone called a 95% MTBC
    sample "M. orygis 2.9%". A candidate nested inside a better one (or holding
    it) is the same answer at another depth, so the runner-up is the best
    candidate from a DIFFERENT lineage. Returns [(pct, name), …].
    """
    cands = []           # (pct, name, path_of_names)
    stack: List[tuple] = []  # (indent, name)
    with open(report, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            raw = parts[-1]
            name = raw.strip()
            indent = len(raw) - len(raw.lstrip(" "))
            while stack and stack[-1][0] >= indent:
                stack.pop()
            stack.append((indent, name))
            rank = parts[-3].strip()
            if not (rank == "S" or re.match(r"^S\d+$", rank) or re.match(r"^G\d+$", rank)):
                continue
            try:
                pct = float(parts[0])
            except ValueError:
                continue
            cands.append((pct, name, tuple(n for _, n in stack)))
    cands.sort(key=lambda c: -c[0])
    picked: List[tuple] = []
    for pct, name, path in cands:
        if any(name in other_path or other_name in path for _, other_name, other_path in picked):
            continue
        picked.append((pct, name, path))
        if len(picked) == 2:
            break
    return [(pct, name) for pct, name, _ in picked]


def _rp_metrics(run_dir: Path) -> Dict[str, Any]:
    """What Kraken found (top call and runner-up) plus what the run was for.

    Parsed from the Kraken report rather than a summary file because this
    pipeline does not write one; see _rp_kraken_top for how the calls are
    chosen and _rp_run_facts for the rest."""
    empty = {"top_taxon": None, "top_pct": None, "runner_up": None}
    try:
        reports = sorted(run_dir.glob("kraken/*_reportkraken.txt")) or sorted(run_dir.rglob("*_reportkraken.txt"))
        top2 = _rp_kraken_top(max(reports, key=lambda r: r.stat().st_mtime)) if reports else []
    except OSError:
        top2 = []
    top = top2[0] if top2 else None
    second = top2[1] if len(top2) > 1 else None
    def _pct(x: float) -> str:
        return "%.1f%%" % x if x >= 1 else "%.2f%%" % x if x >= 0.01 else "<0.01%"
    out = dict(empty)
    if top:
        out.update(top_taxon=top[1], top_pct=_pct(top[0]))
    if second:
        out["runner_up"] = "%s (%s)" % (second[1], _pct(second[0]))
    out.update(_rp_run_facts(run_dir))
    return out


_RP_MODE_LABEL = {"kraken_only": "Krona only", "parse": "Parsed reads", "full": "Full ID"}


def _rp_last_command(run_dir: Path) -> str:
    """The pipeline command line of the latest run in this folder ("" if unknown).

    The GUI appends one block per run to .provenance/kraken_id_parse_gui_run_cmd.txt;
    the last non-comment line is the command that produced what is there now."""
    try:
        text = (run_dir / ".provenance" / "kraken_id_parse_gui_run_cmd.txt").read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    lines = [l.strip() for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")]
    return lines[-1] if lines else ""


def _rp_short_subject(description: str) -> str:
    """'NC_000962 Mycobacterium tuberculosis H37Rv, complete genome' -> the organism."""
    d = re.sub(r"^\S+\s+", "", (description or "").strip(), count=1) if re.match(r"^[A-Z]{1,4}_?[A-Z]*\d", description or "") else (description or "").strip()
    d = re.sub(r",\s*(complete genome|complete sequence|whole genome shotgun sequence|chromosome.*)$", "", d, flags=re.I)
    return d.strip()


def _rp_run_facts(run_dir: Path) -> Dict[str, Any]:
    """What was asked of the run and what it answered, for the Results table.

    mode            Krona only / Parsed reads / Full ID
    target          the target name the reads were parsed for
    target_pct      share of reads extracted for it ("95.3%")
    target_reads    reads extracted ("559,228 of 586,737")
    identification  the strongest BLAST match's organism (Full ID runs)
    ident_share     that match's share of the assembled bases ("92%")
    """
    facts: Dict[str, Any] = {"mode": None, "target": None, "target_pct": None,
                             "target_reads": None, "identification": None, "ident_share": None}
    cmd = _rp_last_command(run_dir)
    man = _rp_json(run_dir, "run_manifest.json")
    legacy = man.get("legacy_report") if isinstance(man.get("legacy_report"), dict) else {}
    target = legacy.get("target") if isinstance(legacy.get("target"), dict) else {}
    blast = [r for r in (legacy.get("blast_rows") or []) if isinstance(r, dict)]
    # What the report says wins: the manifest's run_mode (recorded since
    # 2026-11), else what the manifest's contents show. A later Kraken-only run
    # in the same folder keeps a richer report (kraken_only_report.py), so the
    # row keeps describing it. No manifest: the recorded command, then the files.
    if man.get("run_mode") in _RP_MODE_LABEL:
        mode = man["run_mode"]
    elif man:
        if blast or any(run_dir.glob("*_denovo.fasta")) or any(run_dir.glob("*_blast_summary.txt")):
            mode = "full"
        elif target or "--no-blast" in cmd:
            mode = "parse"
        else:
            mode = "kraken_only" if ("--kraken-only" in cmd or "kraken_only_report.py" in cmd) else "parse"
    elif "--kraken-only" in cmd or "kraken_only_report.py" in cmd:
        mode = "kraken_only"
    elif "--no-blast" in cmd:
        mode = "parse"
    elif cmd:
        mode = "full"
    else:
        mode = "kraken_only" if any(run_dir.glob("kraken/*_krona.html")) else None
    facts["mode"] = _RP_MODE_LABEL.get(mode) if mode else None
    if mode == "kraken_only":
        return facts
    taxon = target.get("taxon") or (man.get("parameters") or {}).get("taxon")
    facts["target"] = taxon or None
    rate = target.get("extraction_rate")
    try:
        facts["target_pct"] = "%.1f%%" % float(str(rate).rstrip("%")) if rate not in (None, "") else None
    except ValueError:
        facts["target_pct"] = None
    ex, tot = target.get("extracted_reads"), target.get("total_input_reads")
    try:
        if ex not in (None, ""):
            facts["target_reads"] = f"{int(ex):,}" + (f" of {int(tot):,}" if tot not in (None, "") else "")
    except (TypeError, ValueError):
        pass
    if blast:
        def _nt(r):
            try:
                return int(str(r.get("nt_base_count") or "0").replace(",", ""))
            except ValueError:
                return 0
        best = max(blast, key=_nt)
        total = sum(_nt(r) for r in blast)
        facts["identification"] = _rp_short_subject(str(best.get("description") or "")) or None
        if total:
            facts["ident_share"] = "%.0f%%" % (100.0 * _nt(best) / total)
    return facts


# No cross-tool probes: the entry the shared scaffold shipped with pointed at
# this tool's OWN Krona chart through a route that only exists in the vSNP
# backend (/api/projects/{p}/kraken/samples/{s}/krona), so every sample with a
# Krona grew a dead "from kraken" link in the Files cell. The Krona chart is a
# first-class file of this tool's runs — it is in `files` (and the pane's
# quick-open column) already.
_RP_CROSS_PROBES: List = []


def _rp_cross_tool(project: str, project_dir: Path, sample: str) -> List[Dict]:
    """Sibling tools' outputs for the same sample.

    Every tool builds the same project skeleton, so a sample analysed here may
    also have a Kraken run. Those files live outside this tool's run dir, which is
    exactly why they never appeared in a results view."""
    out: List[Dict] = []
    for tool, kind, label, subdir, glob, href in _RP_CROSS_PROBES:
        base = project_dir / subdir
        if not base.is_dir():
            continue
        d = base / sample
        if not d.is_dir():
            try:
                cands = sorted(x for x in base.iterdir()
                               if x.is_dir() and x.name.startswith(sample + "_"))
            except OSError:
                cands = []
            d = cands[0] if cands else None
        if d is None:
            continue
        try:
            if not any(d.rglob(glob)):
                continue
        except OSError:
            continue
        out.append({"tool": tool, "kind": kind, "label": label,
                    "href": href.format(p=project, s=sample)})
    return out


def _rp_rows(name: str, include_all: bool = False) -> List[Dict]:
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, "Project not found: %s" % name)
    root = project_dir / _RP_SUBDIR
    rows: List[Dict] = []
    if not root.is_dir():
        return rows
    try:
        entries = sorted(p for p in root.iterdir() if p.is_dir())
    except (OSError, PermissionError):
        return rows
    for d in entries:
        rows.append({
            "sample": d.name,
            "status": _rp_status(d),
            "run_date": _rp_finished_at(d),
            "run_dir": str(d),
            "flags": _rp_flags(d),
            "metrics": _rp_metrics(d),
            "files": _collect_result_files(d, include_all),
            "cross_tool": _rp_cross_tool(name, project_dir, d.name),
        })
    rows.sort(key=lambda r: (r["run_date"] or "", r["sample"]), reverse=True)
    return rows


def _rp_filter(rows: List[Dict], start: str, end: str, q: str) -> List[Dict]:
    """Apply the pane's filters server-side too, so an export matches the view."""
    ql = (q or "").strip().lower()
    out = []
    for r in rows:
        if ql and ql not in r["sample"].lower():
            continue
        day = (r.get("run_date") or "")[:10]
        if start and (not day or day < start):
            continue
        if end and (not day or day > end):
            continue
        out.append(r)
    return out


@app.get("/api/projects/{name}/results")
def api_project_results(name: str, all: int = Query(0)):
    return JSONResponse({"project": name, "tool": _RP_SUBDIR,
                         "rows": _rp_rows(name, include_all=bool(all))})


_RP_EXPORT_COLUMNS = [("sample", "Sample"), ("status", "Status"),
                      ("run_date", "Run date"), ("qc", "QC"),
                      ("qc_reasons", "QC notes")] + [("mode", "Run mode"), ("target", "Target name"),
                                                    ("target_pct", "Target reads %"), ("target_reads", "Target reads"),
                                                    ("identification", "Identification (BLAST)"),
                                                    ("ident_share", "Share of assembled bases"),
                                                    ("top_taxon", "Top species (Kraken)"), ("top_pct", "Percent"),
                                                    ("runner_up", "Runner-up")] + \
                     [("run_dir", "Run directory")]


def _rp_export_records(name, start, end, q):
    for r in _rp_filter(_rp_rows(name), start, end, q):
        m = r.get("metrics") or {}
        rec = {"sample": r["sample"], "status": r["status"],
               "run_date": (r.get("run_date") or "")[:19],
               "qc": (r.get("flags") or {}).get("level", ""),
               "qc_reasons": "; ".join((r.get("flags") or {}).get("reasons", [])),
               "run_dir": r.get("run_dir", "")}
        for key, _label in _RP_EXPORT_COLUMNS:
            rec.setdefault(key, m.get(key))
        yield rec


@app.get("/api/projects/{name}/results.csv")
def api_project_results_csv(name: str, start: str = Query(""), end: str = Query(""),
                            q: str = Query("")):
    import csv
    import io
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=[k for k, _ in _RP_EXPORT_COLUMNS], extrasaction="ignore")
    w.writerow({k: label for k, label in _RP_EXPORT_COLUMNS})
    for rec in _rp_export_records(name, start, end, q):
        w.writerow(rec)
    return Response(content=buf.getvalue(), media_type="text/csv",
                    headers={"Content-Disposition":
                             'attachment; filename="%s_kraken_results.csv"' % name})


@app.get("/api/projects/{name}/results.xlsx")
def api_project_results_xlsx(name: str, start: str = Query(""), end: str = Query(""),
                             q: str = Query("")):
    try:
        from openpyxl import Workbook
    except ImportError:
        raise HTTPException(501, "Excel export needs openpyxl in this tool's environment.")
    import io
    wb = Workbook()
    ws = wb.active
    ws.title = "Results"
    ws.append([label for _, label in _RP_EXPORT_COLUMNS])
    for rec in _rp_export_records(name, start, end, q):
        ws.append([rec.get(k) for k, _ in _RP_EXPORT_COLUMNS])
    for i, (_k, label) in enumerate(_RP_EXPORT_COLUMNS, start=1):
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = max(12, len(label) + 2)
    ws.freeze_panes = "A2"
    buf = io.BytesIO()
    wb.save(buf)
    return Response(content=buf.getvalue(),
                    media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    headers={"Content-Disposition":
                             'attachment; filename="%s_kraken_results.xlsx"' % name})


if _FRONTEND_DIST.is_dir():
    # Serve index.html explicitly with no-cache so the browser (and the OOD
    # rnode proxy) always revalidate it. The hashed assets it points to are
    # immutable and safe to cache; a stale index.html, however, would keep
    # referencing old asset hashes after a rebuild and silently break styling.
    _INDEX_HTML = _FRONTEND_DIST / "index.html"

    @app.get("/")
    def index():
        return FileResponse(
            _INDEX_HTML,
            headers={"Cache-Control": "no-cache, no-store, must-revalidate"},
        )

    app.mount("/", StaticFiles(directory=str(_FRONTEND_DIST), html=True), name="static")
else:
    @app.get("/")
    def root():
        return JSONResponse(
            {"error": "Frontend not built. Run: cd frontend && npm run build"},
            status_code=503,
        )
