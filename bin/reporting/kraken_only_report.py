"""Run a Kraken-only classification, then write its report.

kraken_id_parse.py's --kraken-only mode stops right after the Krona chart
(``sys.exit(0)`` before the structured report is set up), so a Kraken-only run
— the GUI's default, quickest mode — left no report.html, report.pdf or
run_manifest.json: only the chart. This wrapper runs the pipeline unchanged
and then builds the same structured report every other mode gets, from what
the run left in the working directory.

It lives here, in the GUI-only bin/reporting package, rather than as an edit to
kraken_id_parse.py, which is vendored from the USDA pipeline and merged on
every sync.

    python kraken_only_report.py -r1 R1 [-r2 R2] -k KRAKEN_DB [--no-pdf]
    python kraken_only_report.py --report-only -r1 R1 [-r2 R2] -k KRAKEN_DB

Run from the sample's output folder (the GUI's cwd for the job). The exit code
is the pipeline's; a report that cannot be written is a warning, never a
failed run.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import subprocess
import sys
from pathlib import Path

_BIN = Path(__file__).resolve().parents[1]
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))

from reporting import build_run_manifest, render_html_report, render_pdf_report, write_manifest  # noqa: E402


def _strip_read_tag(name: str) -> str:
    import re
    m = re.search(r"(?:_R([12])(?:_\d+)?|_([12]))\.fastq\.gz$", name, re.IGNORECASE)
    if m:
        return name[:m.start()]
    return name[:-len(".fastq.gz")] if name.endswith(".fastq.gz") else name


def _richer_report_present(output_dir: Path) -> str:
    """The run mode of a parse/full report already in this folder, else "".

    A Krona-only re-run of a sample that already had a full identification
    must not replace that report with a thinner one: the identification is
    still the most useful thing in the folder, and its files are still there.
    """
    import json
    from reporting.summary import infer_run_mode
    try:
        man = json.loads((output_dir / "run_manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    mode = infer_run_mode(man) if isinstance(man, dict) else ""
    return mode if mode in ("parse", "full") else ""


def write_report(r1: str, r2: str | None, kraken_db: str, started: _dt.datetime | None,
                 no_pdf: bool = False, output_dir: Path | None = None) -> Path | None:
    output_dir = Path(output_dir or Path.cwd())
    kept = _richer_report_present(output_dir)
    if kept:
        print(f"This folder already holds a report from a {'full identification' if kept == 'full' else 'parse'} "
              f"run; it is kept. The new Krona chart is in kraken/.", flush=True)
        return None
    sample = _strip_read_tag(Path(r1).name)
    kwargs = dict(
        sample_id=sample,
        status="completed",
        warnings=[],
        inputs={"r1": r1, "r2": r2},
        parameters={"taxon": None, "kraken_db": kraken_db, "blast_db": None},
        output_dir=output_dir,
        started_at=started,
        run_mode="kraken_only",
    )
    manifest = build_run_manifest(**kwargs)
    write_manifest(manifest, output_dir)
    html_path = render_html_report(manifest, output_dir)
    if not no_pdf:
        pdf_path, pdf_warning = render_pdf_report(html_path, output_dir)
        if pdf_warning:
            print(f"WARNING: {pdf_warning}")
            kwargs["warnings"] = [pdf_warning]
        # Rebuild so the artifact list includes the PDF (or says why not).
        manifest = build_run_manifest(**kwargs)
        write_manifest(manifest, output_dir)
        html_path = render_html_report(manifest, output_dir)
    return html_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-r1", required=True)
    ap.add_argument("-r2")
    ap.add_argument("-k", "--kraken_db", required=True)
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--report-only", action="store_true",
                    help="Do not run the pipeline; write the report for a finished Kraken-only run.")
    args = ap.parse_args(argv)

    started = _dt.datetime.now().astimezone()
    rc = 0
    if not args.report_only:
        cmd = [sys.executable, "-u", str(_BIN / "kraken_id_parse.py"), "-r1", args.r1]
        if args.r2:
            cmd += ["-r2", args.r2]
        cmd += ["-k", args.kraken_db, "--kraken-only"]
        if args.no_pdf:
            cmd.append("--no-pdf")
        sys.stdout.flush()
        rc = subprocess.call(cmd, env=os.environ.copy())
        if rc != 0:
            return rc
    print("\nWriting the Kraken-only report (report.html, report.pdf)…", flush=True)
    try:
        html = write_report(args.r1, args.r2, args.kraken_db, started, no_pdf=args.no_pdf)
        if html:
            print(f"HTML report generated: {html}", flush=True)
    except Exception as exc:  # the classification itself succeeded
        print(f"WARNING: the Kraken-only report could not be written: {exc}", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
