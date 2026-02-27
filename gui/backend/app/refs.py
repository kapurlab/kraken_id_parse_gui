from pathlib import Path
from typing import Dict, List


def list_references(repo_path: Path) -> List[Dict]:
    refs: List[Dict] = []
    ref_paths = _load_reference_paths(repo_path)
    for parent in ref_paths:
        if not parent.exists():
            continue
        for child in sorted(parent.iterdir()):
            if child.is_dir():
                refs.append({"name": child.name, "path": str(child)})
    return refs


def _load_reference_paths(repo_path: Path) -> List[Path]:
    deps_file = repo_path / "dependencies" / "reference_options_paths.txt"
    if deps_file.exists():
        paths = []
        for line in deps_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                paths.append(Path(line))
        return paths
    fallback = repo_path / "dependencies"
    return [fallback] if fallback.exists() else []


def get_reference_paths(repo_path: Path) -> List[str]:
    deps_file = repo_path / "dependencies" / "reference_options_paths.txt"
    if not deps_file.exists():
        return []
    paths = []
    for line in deps_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            paths.append(line)
    return paths


def add_reference_path(repo_path: Path, new_path: str) -> List[str]:
    deps_dir = repo_path / "dependencies"
    deps_dir.mkdir(parents=True, exist_ok=True)
    deps_file = deps_dir / "reference_options_paths.txt"
    existing = get_reference_paths(repo_path)
    resolved = str(Path(new_path).resolve())
    if resolved not in [str(Path(p).resolve()) for p in existing]:
        existing.append(resolved)
        deps_file.write_text("\n".join(existing) + "\n", encoding="utf-8")
    return existing


def remove_reference_path(repo_path: Path, remove_path: str) -> List[str]:
    deps_dir = repo_path / "dependencies"
    deps_file = deps_dir / "reference_options_paths.txt"
    if not deps_file.exists():
        return []
    existing = get_reference_paths(repo_path)
    resolved_remove = str(Path(remove_path).resolve())
    updated = [p for p in existing if str(Path(p).resolve()) != resolved_remove]
    deps_file.write_text("\n".join(updated) + "\n", encoding="utf-8")
    return updated
