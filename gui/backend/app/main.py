from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from pathlib import Path
from typing import Dict, List, Optional
import json
import os
import shutil
import time
import logging
import yaml
import subprocess
import sys

from app.config import load_config, save_config
from app.jobs import JobManager
from app.projects import create_project, list_projects, ensure_project_dirs, archive_project, delete_project, update_project_meta
from app.sra import expand_accessions, build_download_script

app = FastAPI(title="Kraken ID Parse GUI API")
logger = logging.getLogger("uvicorn.error")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

cfg = load_config()
projects_root = Path(cfg["projects_root"])
projects_root.mkdir(parents=True, exist_ok=True)
job_manager = JobManager(Path(cfg["projects_root"]) / ".jobs")


class ConfigPayload(BaseModel):
    kraken_repo_path: str
    projects_root: str
    conda_env_path: Optional[str] = ""
    default_preset: Optional[str] = ""


class ProjectPayload(BaseModel):
    name: str = Field(..., min_length=1)


class LinkLocalPayload(BaseModel):
    path: str


class RunPayload(BaseModel):
    preset: str
    overrides: Dict[str, str] = Field(default_factory=dict)


@app.get("/api/health")
def health():
    return {"ok": True}


@app.get("/api/config")
def get_config():
    return load_config()


@app.post("/api/config")
def set_config(payload: ConfigPayload):
    cfg = payload.model_dump()
    save_config(cfg)
    return cfg


def _kraken_paths(cfg: Dict[str, str]) -> Dict[str, Path]:
    repo = Path(cfg.get("kraken_repo_path", "")).expanduser()
    return {
        "repo": repo,
        "run_with_config": repo / "bin" / "run_with_config.py",
        "configs": repo / "internal" / "kraken_configs.yaml"
    }


def _load_presets(cfg: Dict[str, str]) -> Dict:
    paths = _kraken_paths(cfg)
    config_path = paths["configs"]
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with open(config_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data.get("presets", {})


@app.get("/api/presets")
def list_presets():
    cfg = load_config()
    try:
        presets = _load_presets(cfg)
    except FileNotFoundError:
        return []
    out = []
    for name, payload in presets.items():
        out.append({
            "name": name,
            "taxon": payload.get("taxon", ""),
            "kraken_db": payload.get("kraken_db", ""),
            "blast_db": payload.get("blast_db", "")
        })
    out.sort(key=lambda x: x["name"])
    return out


def _check_tool(tool_name: str, conda_env_bin: str = "") -> dict:
    """Check if a command-line tool is available.

    Returns dict with 'found' (bool), 'path' (str or None), 'version' (str or None).
    Searches conda env bin first, then system PATH.
    """
    import shutil as _shutil

    search_path = os.environ.get("PATH", "")
    if conda_env_bin:
        search_path = f"{conda_env_bin}:{search_path}"

    tool_path = _shutil.which(tool_name, path=search_path)
    version = None
    if tool_path:
        try:
            result = subprocess.run(
                [tool_path, "--version"],
                capture_output=True, text=True, timeout=10
            )
            version = (result.stdout.strip() or result.stderr.strip())[:120]
        except Exception:
            version = "(version check failed)"
    return {"found": tool_path is not None, "path": tool_path, "version": version}


@app.get("/api/preflight")
def preflight():
    cfg = load_config()
    paths = _kraken_paths(cfg)
    issues = []
    checks = {
        "kraken_repo_path": paths["repo"].is_dir(),
        "run_with_config": paths["run_with_config"].exists(),
        "kraken_configs": paths["configs"].exists(),
        "projects_root": Path(cfg.get("projects_root", "")).expanduser().is_dir(),
        "conda_env_path": True,
    }
    conda_env = Path(cfg.get("conda_env_path", "")).expanduser()
    conda_env_bin = ""
    if str(conda_env).strip():
        checks["conda_env_path"] = conda_env.is_dir()
        conda_env_bin = str(conda_env / "bin")
    if not checks["kraken_repo_path"]:
        issues.append(f"Missing kraken repo path: {paths['repo']}")
    if not checks["run_with_config"]:
        issues.append(f"Missing run_with_config.py: {paths['run_with_config']}")
    if not checks["kraken_configs"]:
        issues.append(f"Missing kraken_configs.yaml: {paths['configs']}")
    if not checks["projects_root"]:
        issues.append(f"Missing projects root: {cfg.get('projects_root', '')}")
    if not checks["conda_env_path"]:
        issues.append(f"Missing conda env path: {cfg.get('conda_env_path', '')}")

    # --- Preflight tool checks ---
    required_tools = ["kraken2", "blastn", "spades.py", "bwa", "samtools", "pdflatex"]
    optional_tools = ["bracken", "picard", "freebayes", "seqkit"]
    tool_checks = {}
    for tool in required_tools:
        info = _check_tool(tool, conda_env_bin)
        tool_checks[tool] = info
        if not info["found"]:
            issues.append(f"Required tool not found: {tool}")
    for tool in optional_tools:
        tool_checks[tool] = _check_tool(tool, conda_env_bin)

    return {
        "ok": not issues,
        "issues": issues,
        "checks": checks,
        "tools": tool_checks,
    }


@app.get("/api/projects")
def get_projects():
    cfg = load_config()
    return list_projects(Path(cfg["projects_root"]))


@app.post("/api/projects")
def create_project_api(payload: ProjectPayload):
    cfg = load_config()
    projects_root = Path(cfg["projects_root"])
    projects_root.mkdir(parents=True, exist_ok=True)
    project_dir = create_project(projects_root, payload.name)
    return {"name": project_dir.name}


@app.post("/api/projects/{project}/archive")
def archive_project_api(project: str):
    cfg = load_config()
    archive_project(Path(cfg["projects_root"]), project)
    return {"archived": True}


@app.delete("/api/projects/{project}")
def delete_project_api(project: str):
    """Archive the project (move to projects_archive/) instead of permanent delete."""
    cfg = load_config()
    try:
        archive_path = archive_project(Path(cfg["projects_root"]), project)
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    return {"archived": True, "archive_path": str(archive_path)}


@app.post("/api/projects/{project}/link-local")
def link_local(project: str, payload: LinkLocalPayload):
    cfg = load_config()
    project_dir = Path(cfg["projects_root"]) / project
    if not project_dir.exists():
        raise HTTPException(status_code=404, detail="Project not found")
    ensure_project_dirs(project_dir)
    source = Path(payload.path).expanduser()
    if not source.exists():
        raise HTTPException(status_code=400, detail=f"Path not found: {source}")
    download_dir = project_dir / "download"
    fastqs = list(source.rglob("*.fastq.gz"))
    if not fastqs:
        raise HTTPException(status_code=400, detail="No .fastq.gz files found")
    linked = 0
    for f in fastqs:
        # Skip hidden or empty files and anything under runs/kraken subfolders
        try:
            rel_parts = f.relative_to(source).parts
        except ValueError:
            rel_parts = f.parts
        if f.name.startswith("."):
            continue
        if "runs" in rel_parts or "kraken" in rel_parts:
            continue
        try:
            if f.stat().st_size == 0:
                continue
        except OSError:
            continue
        target = download_dir / f.name
        if target.exists():
            continue
        try:
            target.symlink_to(f)
            linked += 1
        except OSError:
            shutil.copy2(f, target)
            linked += 1
    update_project_meta(project_dir, {"last_input_path": str(source)})
    return {"linked": linked}


def _build_env(cfg: Dict[str, str]) -> Dict[str, str]:
    conda_env = cfg.get("conda_env_path", "").strip()
    if not conda_env:
        return {"MPLBACKEND": "Agg"}
    env_bin = str(Path(conda_env) / "bin")
    jvm_home = Path(conda_env) / "lib" / "jvm"
    current_path = os.environ.get("PATH", "")
    env = {"PATH": f"{env_bin}:{current_path}", "MPLBACKEND": "Agg"}
    if jvm_home.exists():
        env["JAVA_HOME"] = str(jvm_home)
    return env


def _collect_outputs(run_dir: Path) -> List[Dict[str, str]]:
    exts = (".pdf", ".xlsx", ".xls", ".csv", ".png", ".html", ".txt")
    outputs = []
    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path.suffix.lower() in exts:
            outputs.append({
                "name": path.name,
                "path": str(path),
                "relative": str(path.relative_to(run_dir))
            })
    return outputs


@app.post("/api/projects/{project}/run")
def run_project(project: str, payload: RunPayload):
    cfg = load_config()
    project_dir = Path(cfg["projects_root"]) / project
    if not project_dir.exists():
        raise HTTPException(status_code=404, detail="Project not found")
    ensure_project_dirs(project_dir)

    presets = _load_presets(cfg)
    if payload.preset not in presets:
        raise HTTPException(status_code=400, detail=f"Unknown preset: {payload.preset}")

    download_dir = project_dir / "download"
    fastqs = list(download_dir.rglob("*.fastq.gz"))
    if not fastqs:
        raise HTTPException(status_code=400, detail="No FASTQ files in download directory")

    run_id = time.strftime("%Y%m%d_%H%M%S")
    run_dir = project_dir / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # Link FASTQs into run folder for pattern matching
    for f in fastqs:
        target = run_dir / f.name
        if not target.exists():
            target.symlink_to(f)

    repo_paths = _kraken_paths(cfg)
    run_with_config = repo_paths["run_with_config"]
    if not run_with_config.exists():
        raise HTTPException(status_code=400, detail=f"Missing run_with_config.py at {run_with_config}")

    clean_overrides: Dict[str, str] = {}
    for key, value in payload.overrides.items():
        if value is None:
            continue
        if isinstance(value, str):
            v = value.strip()
            if v.startswith(f"{key}="):
                v = v.split("=", 1)[1].strip()
            if not v:
                continue
            clean_overrides[key] = v
        else:
            clean_overrides[key] = str(value)

    # Ensure run_with_config finds FASTQs even if run_dir is empty
    clean_overrides.setdefault("r1_pattern", str(download_dir / "*_R1*.fastq.gz"))
    clean_overrides.setdefault("r2_pattern", str(download_dir / "*_R2*.fastq.gz"))

    cmd = ["python", str(run_with_config), "--preset", payload.preset]
    for key, value in clean_overrides.items():
        cmd.extend(["--override", f"{key}={value}"])
    cmd_str = " ".join(cmd)

    job_id = job_manager.start_job(
        name=f"kraken:{project}:{run_id}",
        command=cmd_str,
        cwd=run_dir,
        env=_build_env(cfg)
    )
    (run_dir / ".job_id").write_text(job_id, encoding="utf-8")
    (run_dir / "run_meta.json").write_text(json.dumps({
        "run_id": run_id,
        "preset": payload.preset,
        "overrides": clean_overrides
    }, indent=2), encoding="utf-8")

    update_project_meta(project_dir, {"last_run_id": run_id, "last_preset": payload.preset})
    return {"job_id": job_id, "run_id": run_id}


@app.get("/api/projects/{project}/runs")
def list_runs(project: str):
    cfg = load_config()
    project_dir = Path(cfg["projects_root"]) / project
    runs_dir = project_dir / "runs"
    if not runs_dir.exists():
        return []
    runs = []
    for run_dir in sorted([d for d in runs_dir.iterdir() if d.is_dir()], reverse=True):
        job_id_path = run_dir / ".job_id"
        job_id = job_id_path.read_text(encoding="utf-8").strip() if job_id_path.exists() else ""
        job = job_manager.get_job(job_id) if job_id else None
        status = None
        if job:
            status = job.get("status")
        else:
            log_path = Path(cfg["projects_root"]) / ".jobs" / f"{job_id}.log"
            if log_path.exists():
                text = log_path.read_text(errors="ignore")
                finished = ("# finished_at_utc:" in text) or ("# exit_code:" in text)
                has_error = any(token in text for token in ("Traceback", "Exception", "Error:", "ERROR"))
                if finished:
                    status = "failed" if has_error else "succeeded"
                else:
                    status = "failed" if has_error else "running"
            else:
                status = "unknown"
        runs.append({
            "run_id": run_dir.name,
            "job_id": job_id,
            "status": status or "unknown",
            "outputs": len(_collect_outputs(run_dir))
        })
    return runs


@app.get("/api/projects/{project}/runs/{run_id}/outputs")
def run_outputs(project: str, run_id: str):
    cfg = load_config()
    run_dir = Path(cfg["projects_root"]) / project / "runs" / run_id
    if not run_dir.exists():
        raise HTTPException(status_code=404, detail="Run not found")
    return _collect_outputs(run_dir)


class OpenPayload(BaseModel):
    path: str


@app.post("/api/projects/{project}/runs/{run_id}/open")
def open_output(project: str, run_id: str, payload: OpenPayload):
    cfg = load_config()
    run_dir = Path(cfg["projects_root"]) / project / "runs" / run_id
    if not run_dir.exists():
        raise HTTPException(status_code=404, detail="Run not found")
    target = (run_dir / payload.path).resolve()
    if not str(target).startswith(str(run_dir.resolve())):
        raise HTTPException(status_code=400, detail="Invalid path")
    if not target.exists():
        raise HTTPException(status_code=404, detail="File not found")
    if sys.platform == "darwin":
        subprocess.run(["open", str(target)], check=False)
    elif sys.platform.startswith("linux"):
        subprocess.run(["xdg-open", str(target)], check=False)
    else:
        raise HTTPException(status_code=400, detail="Open not supported on this OS")
    return {"opened": True}


@app.delete("/api/projects/{project}/runs/{run_id}")
def delete_run(project: str, run_id: str):
    cfg = load_config()
    run_dir = Path(cfg["projects_root"]) / project / "runs" / run_id
    if not run_dir.exists():
        raise HTTPException(status_code=404, detail="Run not found")
    job_id_path = run_dir / ".job_id"
    if job_id_path.exists():
        job_id = job_id_path.read_text(encoding="utf-8").strip()
        job = job_manager.get_job(job_id) if job_id else None
        if job and job.get("status") == "running":
            raise HTTPException(status_code=409, detail="Run is still running")
    shutil.rmtree(run_dir)
    return {"deleted": True}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@app.get("/api/jobs/{job_id}/events")
def job_events(job_id: str):
    job = job_manager.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    log_path = Path(job["log_path"])

    def event_stream():
        last_size = 0
        while True:
            if log_path.exists():
                size = log_path.stat().st_size
                if size > last_size:
                    with open(log_path, "r", encoding="utf-8") as f:
                        f.seek(last_size)
                        chunk = f.read()
                        for line in chunk.splitlines():
                            yield f"data: {line}\n\n"
                    last_size = size
            job = job_manager.get_job(job_id)
            if job and job["status"] in {"succeeded", "failed"}:
                yield f"data: [job:{job['status']}]\n\n"
                break
            time.sleep(0.5)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


# ---------------------------------------------------------------------------
# SRA Download
# ---------------------------------------------------------------------------

class SRAPayload(BaseModel):
    accessions: str = Field(..., min_length=1, description="Comma or space separated SRA accession(s)")


@app.post("/api/projects/{project}/sra-download")
def sra_download(project: str, payload: SRAPayload):
    """Download SRA FASTQ files into a project's download directory."""
    cfg = load_config()
    project_dir = Path(cfg["projects_root"]) / project
    if not project_dir.exists():
        raise HTTPException(status_code=404, detail="Project not found")
    ensure_project_dirs(project_dir)

    # Parse accessions (accept comma, space, newline separated)
    raw = payload.accessions.replace(",", " ").replace("\n", " ")
    raw_list = [a.strip() for a in raw.split() if a.strip()]
    if not raw_list:
        raise HTTPException(status_code=400, detail="No accessions provided")

    # Expand study/experiment accessions to run accessions
    expanded = expand_accessions(raw_list)
    if not expanded:
        raise HTTPException(status_code=400, detail="Could not resolve any run accessions")

    download_dir = project_dir / "download"
    download_dir.mkdir(parents=True, exist_ok=True)

    # Build and write bash script
    script_content = build_download_script(download_dir, expanded, allow_insecure_https=False)
    script_path = download_dir / "sra_download.sh"
    script_path.write_text(script_content, encoding="utf-8")
    script_path.chmod(0o755)

    # Run as a tracked job
    env = _build_env(cfg)
    job_id = job_manager.start_job(
        name=f"sra-download:{project}:{','.join(expanded[:3])}",
        command=f"bash {script_path}",
        cwd=download_dir,
        env=env,
    )

    return {"job_id": job_id, "accessions": expanded}
