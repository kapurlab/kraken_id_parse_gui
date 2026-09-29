"""The Orbivirus (BTV/EHD) analysis is beta, and every output says so.

It is untested and its results cannot be trusted, so a BTV/EHD result must
never reach anyone without that caveat beside it: not in either report (HTML or
PDF), not in the stats workbook, and not as a PASS in the Results pane — also
for Orbivirus runs made before the notice existed.

Run directly:  <conda>/bin/python backend/app/test_orbivirus_beta.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from app import main as M  # noqa: E402
import orbivirus_beta as B  # noqa: E402
from reporting.html_renderer import render_html_text  # noqa: E402
from reporting.manifest import build_run_manifest  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    if cond:
        print(f"  OK  {label}")
    else:
        FAILED += 1
        print(f"  FAIL {label}")


def manifest(taxon, beta=None, return_code=None):
    out = Path(tempfile.mkdtemp(prefix="kip_run_"))
    man = build_run_manifest(sample_id="S1", status="completed", inputs={"r1": "S1_R1.fastq.gz"},
                             parameters={"taxon": taxon, "kraken_db": "db", "blast_db": "nt"},
                             output_dir=out, beta_notice=beta)
    if return_code is not None:
        man["return_code"] = return_code
    (out / "run_manifest.json").write_text(json.dumps(man), encoding="utf-8")
    return out, man


def main() -> int:
    for t in ("Orbivirus", "Bluetongue virus", "Epizootic hemorrhagic disease virus",
              "Orbivirus caerulinguae", "BTV", "EHDV"):
        check(f"taxon {t!r} is beta", B.applies_to_taxon(t))
    for t in ("Mycobacterium tuberculosis complex", "Viruses", "Isavirus salaris", "", None):
        check(f"taxon {t!r} is not", not B.applies_to_taxon(t))
    check("the Results chip reason fits its 48 characters", len(B.CHIP) <= 48)

    # --- the suite report --------------------------------------------------
    _, orbi = manifest("Orbivirus", beta=B.notice())
    _, old_orbi = manifest("Orbivirus")                    # written before the notice existed
    _, viruses = manifest("Viruses", beta=B.notice())      # BTV/EHD found under a broader taxon
    _, mtb = manifest("Mycobacterium tuberculosis complex")
    check("the manifest records the notice", orbi.get("beta_notice") == B.notice())
    for label, man, want in (("an Orbivirus run", orbi, True),
                             ("an older Orbivirus run", old_orbi, True),
                             ("BTV/EHD under Viruses", viruses, True),
                             ("an MTB run", mtb, False)):
        html = render_html_text(man)
        check(f"report.html banner for {label}: {want}", (B.TITLE in html) is want)

    # --- the Results pane ----------------------------------------------------
    d, _ = manifest("Orbivirus", beta=B.notice())
    f = M._rp_flags(d)
    check("an Orbivirus run is REVIEW, not PASS", f["level"] == "review")
    check("...and says why", B.CHIP in f["reasons"])
    d, _ = manifest("Orbivirus")
    check("an older Orbivirus run is REVIEW too", M._rp_flags(d)["level"] == "review")
    d, _ = manifest("Mycobacterium tuberculosis complex")
    check("an MTB run still passes", M._rp_flags(d) == {"level": "pass", "reasons": []})
    d, _ = manifest("Orbivirus", beta=B.notice(), return_code=1)
    f = M._rp_flags(d)
    check("a failed Orbivirus run stays FAIL, failure first",
          f["level"] == "fail" and f["reasons"][0].startswith("the analysis exited")
          and B.CHIP in f["reasons"])

    # --- the pipeline's own report (HTML + PDF) ---------------------------------
    from report_html import HtmlReport
    rep = HtmlReport("S1", out_dir=tempfile.mkdtemp())
    rep.add_overview([("Target Taxon", "Orbivirus")])
    rep.add_beta_notice(B.TITLE, B.TEXT, B.PAGE)
    html, pdf = rep._render("html"), rep._render("pdf")
    check("the notice sits above the run overview",
          B.TITLE in html and html.index(B.TITLE) < html.index("Run Overview"))
    check("every PDF page carries the beta line", "@top-left" in pdf and B.PAGE in pdf)
    plain = HtmlReport("S2", out_dir=tempfile.mkdtemp())
    plain.add_overview([("Target Taxon", "Mycobacterium tuberculosis complex")])
    check("a non-Orbivirus report has neither", B.TITLE not in plain._render("html")
          and "@top-left" not in plain._render("pdf"))
    rep.add_serotype("BTV-17", "Serotype assigned from Segment 2 (VP2) only: BTV-17", [],
                     caveat=B.SEROTYPE)
    check("the serotyping section carries the caveat", B.SEROTYPE in rep._render("html"))

    # --- the stats workbook ------------------------------------------------------
    from btv_serotyping import BTVSerotyping

    class Fake:
        interpretation = "Serotype assigned from Segment 2 (VP2) only: BTV-17"
        predictions = []

        def get_consensus_serotype(self):
            return "BTV-17"

    row = {}
    BTVSerotyping.excel(Fake(), row)
    check("the workbook's serotype interpretation carries the caveat",
          row["Serotype Interpretation"].endswith(f"({B.SHORT})") and row["BTV Serotype"] == "BTV-17")

    print("PASS" if not FAILED else f"{FAILED} FAILED")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
