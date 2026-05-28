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
from fastapi import FastAPI, HTTPException, Request
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


def _list_projects_from_root(root: Path, scope: str) -> List[Dict]:
    if not root.is_dir():
        return []
    projects = []
    for p in sorted(root.iterdir(), key=lambda d: d.stat().st_mtime if d.is_dir() else 0, reverse=True):
        if not p.is_dir() or p.name.startswith("."):
            continue
        download_dir = p / "download"
        fastq_count = len(list(download_dir.rglob("*.fastq.gz"))) if download_dir.is_dir() else 0
        kraken_runs = []
        kraken_dir = p / "kraken"
        if kraken_dir.is_dir():
            kraken_runs = [d.name for d in sorted(kraken_dir.iterdir()) if d.is_dir()]
        projects.append({
            "name": p.name,
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


def _list_fastq_pairs(download_dir: Path) -> List[Dict]:
    """Return paired FASTQ files as {sample, r1, r2} dicts."""
    r1_files = sorted(download_dir.glob("*_R1*.fastq.gz"))
    pairs = []
    for r1 in r1_files:
        sample = re.sub(r"_R1.*", "", r1.name)
        # Find matching R2 in same directory
        r2_candidates = list(download_dir.glob(f"*{sample}*_R2*.fastq.gz"))
        r2 = r2_candidates[0] if r2_candidates else None
        pairs.append({
            "sample": sample,
            "r1": str(r1),
            "r1_name": r1.name,
            "r2": str(r2) if r2 else None,
            "r2_name": r2.name if r2 else None,
        })
    # Also catch singletons (no _R1 pattern) not already included
    all_fq = set(download_dir.glob("*.fastq.gz"))
    covered = {Path(p["r1"]) for p in pairs} | {Path(p["r2"]) for p in pairs if p["r2"]}
    for fq in sorted(all_fq - covered):
        pairs.append({
            "sample": re.sub(r"\.fastq\.gz$", "", fq.name),
            "r1": str(fq),
            "r1_name": fq.name,
            "r2": None,
            "r2_name": None,
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

    # Derive sample name (strip everything from _R1 onward, or .fastq.gz)
    sample_name = re.sub(r"[._]R1.*", "", r1.name)
    sample_name = re.sub(r"\.fastq\.gz$", "", sample_name)

    # Output directory: <project>/kraken/<sample_name>/
    run_dir = project_dir / "kraken" / sample_name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    script = _BIN_DIR / "kraken_id_parse.py"
    cmd_parts = [f'python "{script}"', f'-r1 "{payload.r1}"']
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
                            yield f"data: {line.rstrip()}\n\n"
                        position += len(chunk.encode("utf-8"))
            if current_job and current_job["status"] in ("succeeded", "failed"):
                yield "data: [DONE]\n\n"
                break
            await asyncio.sleep(0.5)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


@app.get("/api/jobs/{job_id}/results")
def api_job_results(job_id: str):
    """List output files in the job's run directory."""
    job = job_manager.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Job not found")
    cwd = job.get("cwd")
    if not cwd or not Path(cwd).is_dir():
        return JSONResponse([])
    run_dir = Path(cwd)
    files = []
    for p in sorted(run_dir.rglob("*")):
        if p.is_file() and not p.name.endswith(".log"):
            rel = p.relative_to(run_dir)
            size = p.stat().st_size
            files.append({"name": str(rel), "size": size})
    return JSONResponse(files)


# ---------------------------------------------------------------------------
# Static frontend — must be last (catches everything not matched above)
# ---------------------------------------------------------------------------
if _FRONTEND_DIST.is_dir():
    app.mount("/", StaticFiles(directory=str(_FRONTEND_DIST), html=True), name="static")
else:
    @app.get("/")
    def root():
        return JSONResponse(
            {"error": "Frontend not built. Run: cd frontend && npm run build"},
            status_code=503,
        )
