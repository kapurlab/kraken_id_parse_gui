import json
import logging
import os
import subprocess
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


class JobManager:
    def __init__(self, jobs_dir: Path):
        self.jobs_dir = jobs_dir
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._jobs: Dict[str, Dict] = {}
        self._restore_jobs()

    def _restore_jobs(self) -> None:
        """Re-attach to any pipeline subprocess that outlived a previous uvicorn instance."""
        for state_file in self.jobs_dir.glob("*.json"):
            try:
                job = json.loads(state_file.read_text())
            except (json.JSONDecodeError, OSError):
                continue
            if job.get("status") != "running":
                continue
            job_id = job["id"]
            pid_path = self.jobs_dir / f"{job_id}.pid"
            alive = False
            pid = None
            if pid_path.exists():
                try:
                    pid = int(pid_path.read_text().strip())
                    os.kill(pid, 0)
                    alive = True
                except (ValueError, OSError):
                    pass
            if alive:
                logger.info("Restored running job %s (pid %s)", job_id, pid)
                with self._lock:
                    self._jobs[job_id] = job
                t = threading.Thread(
                    target=self._watch_pid,
                    args=(job_id, pid, Path(job["log_path"]), pid_path, state_file),
                    daemon=True,
                )
                t.start()
            else:
                pid_path.unlink(missing_ok=True)
                job["status"] = "failed"
                job["exit_code"] = -1
                try:
                    state_file.write_text(json.dumps(job))
                except OSError:
                    pass

    def _watch_pid(
        self, job_id: str, pid: int, log_path: Path, pid_path: Path, state_path: Path
    ) -> None:
        """Poll until a detached subprocess exits, then update job state."""
        while True:
            try:
                os.kill(pid, 0)
                time.sleep(2)
            except OSError:
                break
        pid_path.unlink(missing_ok=True)
        finished_at = datetime.now(timezone.utc)
        # Determine success by checking if the pipeline wrote its completion marker
        status = "failed"
        try:
            if "# finished_at_utc:" in log_path.read_text(encoding="utf-8", errors="replace"):
                status = "succeeded"
        except OSError:
            pass
        exit_code = 0 if status == "succeeded" else -1
        with self._lock:
            job = self._jobs.get(job_id)
            if job:
                job["status"] = status
                job["exit_code"] = exit_code
                job["finished_at"] = finished_at.isoformat()
        try:
            data = json.loads(state_path.read_text())
            data["status"] = status
            data["exit_code"] = exit_code
            data["finished_at"] = finished_at.isoformat()
            state_path.write_text(json.dumps(data))
        except OSError:
            pass

    def start_job(
        self,
        name: str,
        command: str,
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
    ) -> str:
        job_id = uuid.uuid4().hex
        log_path = self.jobs_dir / f"{job_id}.log"
        pid_path = self.jobs_dir / f"{job_id}.pid"
        state_path = self.jobs_dir / f"{job_id}.json"
        started_at = datetime.now(timezone.utc)
        job = {
            "id": job_id,
            "name": name,
            "command": command,
            "cwd": str(cwd) if cwd else None,
            "status": "running",
            "exit_code": None,
            "log_path": str(log_path),
            "started_at": started_at.isoformat(),
            "finished_at": None,
            "duration_seconds": None,
        }
        state_path.write_text(json.dumps(job))
        with self._lock:
            self._jobs[job_id] = job
        thread = threading.Thread(
            target=self._run,
            args=(job_id, command, cwd, env, log_path, pid_path, state_path),
            daemon=True,
        )
        thread.start()
        return job_id

    def get_job(self, job_id: str) -> Optional[Dict]:
        with self._lock:
            return dict(self._jobs[job_id]) if job_id in self._jobs else None

    def list_jobs(self) -> list:
        with self._lock:
            return [dict(j) for j in self._jobs.values()]

    def _run(
        self,
        job_id: str,
        command: str,
        cwd: Optional[Path],
        env: Optional[Dict[str, str]],
        log_path: Path,
        pid_path: Path,
        state_path: Path,
    ) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w", encoding="utf-8") as log:
            started_at = datetime.now(timezone.utc)
            log.write(f"# started_at_utc: {started_at.isoformat()}\n")
            log.write(f"$ {command}\n\n")
            log.flush()
            process = subprocess.Popen(
                command,
                cwd=str(cwd) if cwd else None,
                env={**os.environ, **(env or {})},
                stdout=log,
                stderr=subprocess.STDOUT,
                shell=True,
                text=True,
                start_new_session=True,  # detach from uvicorn's process group
                stdin=subprocess.DEVNULL,
            )
            pid_path.write_text(str(process.pid))
            exit_code = process.wait()
            finished_at = datetime.now(timezone.utc)
            duration = (finished_at - started_at).total_seconds()
            log.write(f"\n# finished_at_utc: {finished_at.isoformat()}\n")
            log.write(f"# duration_seconds: {duration:.2f}\n")
            log.flush()
        pid_path.unlink(missing_ok=True)
        status = "succeeded" if exit_code == 0 else "failed"
        with self._lock:
            job = self._jobs.get(job_id)
            if job:
                job["exit_code"] = exit_code
                job["status"] = status
                job["finished_at"] = finished_at.isoformat()
                job["duration_seconds"] = duration
        try:
            data = json.loads(state_path.read_text())
            data["exit_code"] = exit_code
            data["status"] = status
            data["finished_at"] = finished_at.isoformat()
            data["duration_seconds"] = duration
            state_path.write_text(json.dumps(data))
        except OSError:
            pass
