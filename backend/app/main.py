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
import logging
import os
import re
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import aiofiles
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .config import load_config, save_config
from .jobs import JobManager

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parent.parent          # /srv/kapurlab/tools/kraken_id_parse_gui
_BIN_DIR = _REPO_ROOT / "bin"
_FRONTEND_DIST = _REPO_ROOT / "frontend" / "dist"

# Shared project root
_SHARED_PROJECTS = Path("/srv/kapurlab/projects")

# Jobs log directory (inside repo so it survives across sessions)
_JOBS_DIR = _REPO_ROOT / "backend" / "jobs"

# ---------------------------------------------------------------------------
# App & job manager
# ---------------------------------------------------------------------------
app = FastAPI(title="Kraken ID Parse GUI")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

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
            fastq_count = len(list(download_dir.rglob("*.fastq.gz"))) if download_dir.is_dir() else 0
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


@app.get("/api/projects/{name}/samples")
def api_project_samples(name: str):
    project_dir = _get_project_dir(name)
    if project_dir is None:
        raise HTTPException(404, f"Project not found: {name}")
    download_dir = project_dir / "download"
    if not download_dir.is_dir():
        return JSONResponse([])
    return JSONResponse(_list_fastq_pairs(download_dir))


@app.get("/api/config")
def api_get_config():
    return JSONResponse(load_config())


class ConfigPayload(BaseModel):
    kraken_db: Optional[str] = None
    blast_db: Optional[str] = None
    projects_root: Optional[str] = None
    shared_projects_root: Optional[str] = None


@app.post("/api/config")
def api_save_config(payload: ConfigPayload):
    cfg = load_config()
    updates = payload.model_dump(exclude_none=True)
    cfg.update(updates)
    save_config(cfg)
    return JSONResponse({"ok": True})


class RunPayload(BaseModel):
    project: str
    r1: str           # absolute path to R1 FASTQ
    r2: Optional[str] = None
    taxon: str
    kraken_db: Optional[str] = None
    blast_db: Optional[str] = None


@app.post("/api/run")
def api_run(payload: RunPayload):
    cfg = load_config()
    kraken_db = payload.kraken_db or cfg.get("kraken_db", "")
    blast_db = payload.blast_db or cfg.get("blast_db", "nt")

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

    # Build command
    script = _BIN_DIR / "kraken_id_parse.py"
    # -u: unbuffered stdout/stderr so the parent script's progress prints stream
    # to the log in real time and in order (otherwise Python block-buffers when
    # stdout is a pipe and the log appears to freeze between subprocess outputs).
    cmd_parts = [f'python -u "{script}"', f'-r1 "{payload.r1}"']
    if payload.r2:
        r2 = Path(payload.r2)
        if r2.exists():
            cmd_parts.append(f'-r2 "{payload.r2}"')
    cmd_parts.append(f'-t "{payload.taxon}"')
    if kraken_db:
        cmd_parts.append(f'-k "{kraken_db}"')
    if blast_db:
        cmd_parts.append(f'-b "{blast_db}"')

    command = " ".join(cmd_parts)
    env = {
        "PYTHONPATH": str(_BIN_DIR),
        "PATH": os.environ.get("PATH", ""),
        "PYTHONUNBUFFERED": "1",  # belt-and-suspenders with `python -u`
    }

    job_name = f"{payload.project}/{sample_name} — {payload.taxon}"
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
# Files surfaced first in the Results pane (the ones users actually want).
_PRIORITY_SUFFIXES = ("_report.pdf", "_stats.xlsx", "_krona.html")


def _can_open_inline(name: str) -> bool:
    return Path(name).suffix.lower() in _INLINE_MEDIA


def _media_type_for(name: str) -> str:
    ext = Path(name).suffix.lower()
    return _INLINE_MEDIA.get(ext) or _DOWNLOAD_MEDIA.get(ext) or "application/octet-stream"


@app.get("/api/jobs/{job_id}/results")
def api_job_results(job_id: str):
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
                files.append(
                    {"name": rel, "size": p.stat().st_size, "openable": _can_open_inline(rel)}
                )

    # Surface the pipeline log itself (lives outside run_dir, in the jobs dir).
    log_path = Path(job.get("log_path", ""))
    if log_path.is_file():
        files.append(
            {
                "name": "pipeline_log.txt",
                "size": log_path.stat().st_size,
                "openable": True,
                "is_log": True,
            }
        )

    # Sort: priority files first (PDF, Excel, Krona), then the rest, log last.
    def sort_key(f):
        if f.get("is_log"):
            return (2, f["name"])
        if any(f["name"].endswith(s) for s in _PRIORITY_SUFFIXES):
            return (0, f["name"])
        return (1, f["name"])

    files.sort(key=sort_key)
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
    disposition = "inline" if want_inline else "attachment"
    headers = {"Content-Disposition": f'{disposition}; filename="{display_name}"'}
    return FileResponse(target, media_type=media_type, headers=headers)


# ---------------------------------------------------------------------------
# Static frontend — must be last (catches everything not matched above)
# ---------------------------------------------------------------------------
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
