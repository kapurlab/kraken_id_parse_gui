"""The Results pane lists both reports a run writes.

A full run writes the suite's report (report.html / report.pdf, rendered from
run_manifest.json) and the pipeline's own <sample>_<stamp>_report.html/.pdf,
which carries the organism-specific sections (BTV serotyping, segment status).
The pipeline's report used to be reachable only under "all files". The suite
report stays the primary one: it is what the Open column links to.

Run directly:  <conda>/bin/python backend/app/test_result_categories.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import main as M  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    if cond:
        print(f"  OK  {label}")
    else:
        FAILED += 1
        print(f"  FAIL {label}")


def main() -> int:
    cat = M._result_category
    check("suite report stays report_html", cat("report.html") == "report_html")
    check("suite PDF stays report_pdf", cat("report.pdf") == "report_pdf")
    check("pipeline HTML report is listed",
          cat("S1_2026-09-29_13-17-49_report.html") == "pipeline_report_html")
    check("pipeline PDF report is listed",
          cat("S1_2026-09-29_13-17-49_report.pdf") == "pipeline_report_pdf")
    check("a report-named file in a subfolder is not",
          cat("btv_serotype/S1_report.html") is None)
    check("the interactive coverage chart keeps its own category",
          cat("S1-abc-coverage_interactive.html") == "coverage_interactive")

    order = M._CATEGORY_ORDER
    check("the suite report sorts first",
          order["report_html"] < order["report_pdf"] < order["pipeline_report_html"]
          < order["pipeline_report_pdf"] < order["stats"])
    check("labels say which report is which",
          M._result_label("x_report.html", "pipeline_report_html") == "Pipeline report HTML"
          and M._result_label("report.html", "report_html") == "Report HTML")

    reruns = [
        {"name": "S1_2026-09-01_10-00-00_report.html", "category": "pipeline_report_html", "mtime": 1},
        {"name": "S1_2026-09-29_13-17-49_report.html", "category": "pipeline_report_html", "mtime": 2},
        {"name": "report.html", "category": "report_html", "mtime": 2},
    ]
    kept = [f["name"] for f in M._dedupe_primary_results(reruns)]
    check("a rerun shows only the newest pipeline report",
          "S1_2026-09-29_13-17-49_report.html" in kept
          and "S1_2026-09-01_10-00-00_report.html" not in kept
          and "report.html" in kept)

    print("PASS" if not FAILED else f"{FAILED} FAILED")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
