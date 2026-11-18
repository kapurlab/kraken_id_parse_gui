"""Results rows say what the run found; the Inputs endpoints match the vSNP GUI's.

Results (per sample):
  - the top Kraken call counts species groups (G1), so a 95% M. tuberculosis
    complex sample is not reported as "M. orygis 2.9%"; the runner-up comes from
    a different lineage;
  - the run mode comes from the manifest's run_mode, then the recorded command
    (a Kraken-only run is recorded as the kraken_only_report.py wrapper);
  - target %, reads and the BLAST identification are read from the manifest.
Inputs:
  - link-local takes many paths and says what happened to each;
  - download/ sub-folders (an SRA download into a named folder) are listed and
    their samples runnable; a broken link is listed as such and does not break
    the sample list.

Synthetic data only. Run directly:  <env>/bin/python backend/app/test_results_and_inputs.py
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
TMP = Path(tempfile.mkdtemp(prefix="kgui_rows_"))
os.environ["XDG_CONFIG_HOME"] = str(TMP / "xdg")

from app import main as M  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    print(f"  {'OK ' if cond else 'FAIL'} {label}")
    if not cond:
        FAILED += 1


REPORT = (
    "  2.00\t200\t200\tU\t0\tunclassified\n"
    " 98.00\t9800\t0\tR\t1\troot\n"
    " 96.80\t9680\t10\tG\t1763\t              Mycobacterium\n"
    " 95.40\t9540\t9000\tG1\t77643\t                Mycobacterium tuberculosis complex\n"
    "  2.92\t292\t292\tS\t1305738\t                  Mycobacterium orygis\n"
    "  1.92\t192\t192\tS\t1773\t                  Mycobacterium tuberculosis\n"
    "  0.80\t80\t0\tG\t9605\t      Homo\n"
    "  0.80\t80\t80\tS\t9606\t        Homo sapiens\n"
)


def run_dir(name: str, manifest=None, cmd: str = "") -> Path:
    d = TMP / "proj" / "kraken" / name
    (d / "kraken").mkdir(parents=True)
    (d / "kraken" / f"{name}_reportkraken.txt").write_text(REPORT, encoding="utf-8")
    if manifest is not None:
        (d / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    if cmd:
        (d / ".provenance").mkdir()
        (d / ".provenance" / "kraken_id_parse_gui_run_cmd.txt").write_text(
            f"# run cmd\n# ---\n{cmd}\n", encoding="utf-8")
    return d


def main() -> int:
    print("Results rows")
    full = run_dir("S_full", {
        "legacy_report": {
            "target": {"taxon": "Mycobacterium tuberculosis complex", "extraction_rate": 95.31,
                       "extracted_reads": 559228, "total_input_reads": 586737},
            "blast_rows": [
                {"nt_base_count": "128", "description": "NZ_X1 Other bug strain 9, whole genome shotgun sequence"},
                {"nt_base_count": "3,923,730", "description": "NC_000962 Mycobacterium tuberculosis H37Rv, complete genome"},
            ]}}, cmd="python kraken_id_parse.py -r1 a -t 'Mycobacterium tuberculosis complex' -k db -b blastdb")
    m = M._rp_metrics(full)
    check("top call is the species group, not its biggest species",
          m["top_taxon"] == "Mycobacterium tuberculosis complex" and m["top_pct"] == "95.4%")
    check("runner-up is from another lineage", m["runner_up"] == "Homo sapiens (0.80%)")
    check("full run: mode, target and identification",
          m["mode"] == "Full ID" and m["target_pct"] == "95.3%"
          and m["target_reads"] == "559,228 of 586,737"
          and m["identification"] == "Mycobacterium tuberculosis H37Rv" and m["ident_share"] == "100%")

    konly = run_dir("S_konly", {"run_mode": "kraken_only", "legacy_report": {}},
                    cmd="python bin/reporting/kraken_only_report.py -r1 a -k db")
    m = M._rp_metrics(konly)
    check("Kraken-only (wrapper command, manifest run_mode) reads as Krona only, no target",
          m["mode"] == "Krona only" and m["target"] is None and m["identification"] is None)
    check("...from the command alone, too",
          M._rp_run_facts(run_dir("S_konly2", None, cmd="python x/kraken_only_report.py -r1 a -k db"))["mode"] == "Krona only")
    parse = run_dir("S_parse", {"legacy_report": {"target": {"taxon": "T", "extraction_rate": "64.05"}}},
                    cmd="python kraken_id_parse.py -r1 a -t T -k db --no-blast")
    check("parse run", M._rp_run_facts(parse)["mode"] == "Parsed reads"
          and M._rp_run_facts(parse)["target_pct"] == "64.0%")
    after = run_dir("S_after", {"run_mode": "full", "legacy_report": {"blast_rows": [
        {"nt_base_count": "10", "description": "NC_1 Thing"}]}},
        cmd="python bin/reporting/kraken_only_report.py -r1 a -k db")
    check("a Krona-only re-run after a full ID still shows the full ID the report holds",
          M._rp_run_facts(after)["mode"] == "Full ID")

    print("Inputs")
    proj = TMP / "projects" / "P"
    (proj / "download" / "batch1").mkdir(parents=True)
    src = TMP / "elsewhere"
    src.mkdir()
    for n in ("A_R1.fastq.gz", "A_R2.fastq.gz", "B.fq.gz", "notes.txt"):
        (src / n).write_bytes(b"x" * 10)
    real = M._get_project_dir
    M._get_project_dir = lambda name: proj if name == "P" else None
    try:
        r = json.loads(M.api_project_link_local("P", M.LinkLocalRequest(
            paths=[str(src), str(src / "notes.txt"), str(TMP / "nope.fastq.gz")])).body)
        check("a folder links its .fastq.gz; the rest is accounted for",
              r["linked"] == 2 and r["skipped_not_fastq"] == ["B.fq.gz", "notes.txt"]
              and r["missing"] == [str(TMP / "nope.fastq.gz")])
        r = json.loads(M.api_project_link_local("P", M.LinkLocalRequest(path=str(src / "A_R1.fastq.gz"))).body)
        check("linking again reports it as already there", r["linked"] == 0 and r["already_present"] == ["A_R1.fastq.gz"])
        (proj / "download" / "batch1" / "C_1.fastq.gz").write_bytes(b"y" * 5)
        (proj / "download" / "batch1" / "C_2.fastq.gz").symlink_to(TMP / "gone.fastq.gz")
        files = json.loads(M.api_project_inputs("P").body)["files"]
        names = {f["name"]: f for f in files}
        check("sub-folder files are listed as folder/name",
              "batch1/C_1.fastq.gz" in names and names["batch1/C_2.fastq.gz"].get("dangling") is True)
        samples = json.loads(M.api_project_samples("P").body)
        by = {s["sample"]: s for s in samples}
        check("the sample list survives a broken link and includes the sub-folder sample",
              "A" in by and by["A"]["paired"] and "C" in by and by["C"]["source"] == "download/batch1")
        M.api_project_input_delete("P", "batch1/C_2.fastq.gz")
        check("a sub-folder file can be removed", not (proj / "download" / "batch1" / "C_2.fastq.gz").is_symlink())
        try:
            M.api_project_input_delete("P", "../project.json")
            check("a path out of download/ is refused", False)
        except M.HTTPException as exc:
            check("a path out of download/ is refused", exc.status_code == 400)
    finally:
        M._get_project_dir = real

    print("PASS" if not FAILED else f"{FAILED} FAILURE(S)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
