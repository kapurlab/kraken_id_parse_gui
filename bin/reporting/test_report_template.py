"""report.html: the summary comes first, in every run mode.

Synthetic manifests only (this repository is public): no real sample, host or
path appears here. Covers the three run modes (Kraken + Krona only, Kraken +
parse reads, Kraken + parse + BLAST), the Needs attention rules, the Orbivirus
beta notice, numeric header alignment, the screen/print split of the coverage
figures, and the Kraken overview read from a Kraken2 report.

Run directly:  <env>/bin/python bin/reporting/test_report_template.py
"""
from __future__ import annotations

import base64
import copy
import re
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))

from reporting import summary as S  # noqa: E402
from reporting.html_renderer import _css_string, render_html_text  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    if cond:
        print(f"  OK  {label}")
    else:
        FAILED += 1
        print(f"  FAIL {label}")


KRAKEN_REPORT = "\n".join([
    "  4.00\t400\t400\tU\t0\tunclassified",
    " 96.00\t9600\t10\tR\t1\troot",
    " 90.00\t9000\t5\tD\t2\t  Bacteria",
    " 89.00\t8900\t20\tG\t1000\t    Examplebacter",
    " 85.00\t8500\t7000\tG1\t1001\t      Examplebacter alpha complex",
    " 10.00\t1000\t1000\tS\t1002\t        Examplebacter alpha",
    "  5.00\t500\t500\tS\t1003\t        Examplebacter beta",
    "  6.00\t600\t0\tD\t2759\t  Eukaryota",
    "  6.00\t600\t600\tS\t9999\t    Hostus exemplaris",
]) + "\n"

OVERVIEW = S.kraken_overview(None)  # filled in main()

PNG_1PX = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8/5+hHgAHggJ/PchI7wAAAABJRU5ErkJggg==")


def manifest(mode: str, **over):
    """A run manifest shaped like bin/reporting/manifest.py writes, per mode."""
    m = {
        "schema_version": "0.1",
        "report_date": "January 2, 2026",
        "sample_id": "SAMPLE-A",
        "run_id": "SAMPLE-A",
        "status": "completed",
        "warnings": [],
        "inputs": {"r1": "/data/SAMPLE-A_R1.fastq.gz", "r2": "/data/SAMPLE-A_R2.fastq.gz"},
        "input_sha256": {"r1": "a" * 64, "r2": "b" * 64},
        "parameters": {"taxon": "Examplebacter alpha complex", "kraken_db": "/db/kraken2/k2_example",
                       "blast_db": "/db/blast/example_refs"},
        "software": {"pipeline_git_commit": "abc1234", "python": "3.10.0", "kraken2": "Kraken version 2.1"},
        "environment": {"hostname": "example-host", "cwd": "/nonexistent/run"},
        "timing": {"started_at": "2026-01-02T10:00:00", "finished_at": "2026-01-02T11:00:00"},
        "artifacts": [
            {"name": "SAMPLE-A_2026_krona.html", "path": "kraken/SAMPLE-A_2026_krona.html", "kind": "krona"},
        ],
        "legacy_report": {
            "kraken_overview": copy.deepcopy(OVERVIEW),
            "krona": [{"name": "SAMPLE-A_2026_krona.html", "path": "kraken/SAMPLE-A_2026_krona.html", "kind": "krona"}],
            "fastq": [
                {"read": "R1", "filename": "SAMPLE-A_R1.fastq.gz", "file_size": "10.0 MB", "read_count": "10,000",
                 "q30": "92.0%", "mean_score": "33.1", "avg_length": "150.0"},
                {"read": "R2", "filename": "SAMPLE-A_R2.fastq.gz", "file_size": "10.0 MB", "read_count": "10,000",
                 "q30": "90.0%", "mean_score": "32.0", "avg_length": "150.0"},
            ],
            "target": {"taxon": "Examplebacter alpha complex", "taxon_id": "1001", "total_input_reads": 10000,
                       "extracted_reads": 8500, "extraction_rate": 85.0},
            "target_taxonomy": [
                {"name": "Examplebacter alpha complex", "level": "G1", "reads": "8500", "percent": "85.00"},
                {"name": "Examplebacter alpha", "level": "S", "reads": "1000", "percent": "10.00"},
            ],
            "assembly": {}, "blast_rows": [], "coverage_figures": [], "coverage_interactive": "",
            "coverage_stats": [],
        },
    }
    if mode == "kraken_only":
        m["parameters"]["taxon"] = None
        lr = m["legacy_report"]
        lr["fastq"] = [{"read": "R1", "filename": None}, {"read": "R2", "filename": None}]
        lr["target"] = {}
        lr["target_taxonomy"] = []
    if mode == "parse":
        m["artifacts"].append({"name": "SAMPLE-A_Examplebacter_alpha_complex_R1.fastq.gz",
                               "path": "SAMPLE-A_Examplebacter_alpha_complex_R1.fastq.gz", "kind": "compressed"})
    if mode == "full":
        lr = m["legacy_report"]
        lr["assembly"] = {"contig_count": "40", "contig_bins": "3|7|30", "longest_contig": "120,000",
                          "total_length": "2,000,000", "n50": "60,000", "mean_coverage": "50.0X"}
        # File order is weakest-to-strongest, as the pipeline writes it.
        lr["blast_rows"] = [
            {"nt_base_count": "1,000", "contigs": "1", "description": "NZ_XX000001 Otherus minor strain 1, contig"},
            {"nt_base_count": "1,900,000", "contigs": "38", "description": "NC_999999 Examplebacter alpha strain T, complete genome"},
        ]
        lr["coverage_stats"] = [{"header": "SAMPLE-A_guided_by_NC_999999", "reference": "SAMPLE-A_guided",
                                 "length": 2000000, "mean_coverage": 48.0, "percent_covered": 98.0}]
        lr["coverage_figures"] = [{"name": "cov.pdf", "path": "cov.pdf", "preview_path": "cov.png"}]
        lr["coverage_interactive"] = '<div class="plotly-fragment">chart</div>'
        m["artifacts"].append({"name": "SAMPLE-A_blast_summary.txt", "path": "SAMPLE-A_blast_summary.txt", "kind": "text"})
    for key, value in over.items():
        m[key] = value
    return m


def verdict_title(html: str) -> str:
    match = re.search(r'<div class="verdict-title">(.*?)</div>', html, re.S)
    return match.group(1).strip() if match else ""


def main() -> int:
    global OVERVIEW
    tmp = Path(tempfile.mkdtemp(prefix="kip_report_"))
    (tmp / "kraken").mkdir()
    (tmp / "kraken" / "SAMPLE-A_reportkraken.txt").write_text(KRAKEN_REPORT, encoding="utf-8")
    (tmp / "bracken_pie.png").write_bytes(PNG_1PX)

    print("Kraken overview")
    OVERVIEW = S.kraken_overview(S.find_kraken_report(tmp))
    check("classified / unclassified shares", OVERVIEW["unclassified_pct"] == 4.0 and OVERVIEW["classified_pct"] == 96.0)
    check("total reads = unclassified + root", OVERVIEW["total_reads"] == 10000)
    names = [d["name"] for d in OVERVIEW["domains"]]
    check("domains largest first, then unclassified", names == ["Bacteria", "Eukaryota", "Unclassified"])
    top = OVERVIEW["top_taxa"]
    check("a species group outranks its 10% species",
          top[0]["name"] == "Examplebacter alpha complex" and top[0]["rank_label"] == "Species group")
    child = next(t for t in top if t["name"] == "Examplebacter alpha")
    check("a species nests under its listed group", child["parent"] == "Examplebacter alpha complex")
    host = next(t for t in top if t["name"] == "Hostus exemplaris")
    check("a taxon outside any listed group has no parent", host["parent"] is None)
    check("no report -> empty overview", S.kraken_overview(None) == {})

    print("Rank labels and numbers")
    check("rank codes read as words",
          [S.rank_label(c) for c in ("G1", "S", "S1", "S2", "F", "D")]
          == ["Species group", "Species", "Subspecies / strain", "Strain (level 2)", "Family", "Domain"])
    check("numbers from recorded text", S.num("4,285,951 bp") == 4285951 and S.num("93.11%") == 93.11
          and S.num("68.2X") == 68.2 and S.num("No Reads Found") is None)
    check("CSS strings are escaped", _css_string('a"b\\c') == 'a\\"b\\\\c')

    print("Run modes")
    for mode, label in (("kraken_only", "Kraken + Krona only"), ("parse", "Kraken + parse reads"),
                        ("full", "Kraken + parse + BLAST")):
        m = manifest(mode)
        check(f"{mode}: inferred from what the run produced", S.infer_run_mode(m) == mode)
        try:
            html = render_html_text(m, asset_mode="embed", base_dir=tmp)
        except Exception as exc:  # noqa: BLE001
            check(f"{mode}: renders ({exc})", False)
            continue
        check(f"{mode}: renders, mode named in the header", f"<dd>{label}</dd>" in html)
        check(f"{mode}: summary before the first numbered section",
              0 < html.find('class="verdict') < html.find('<h2 class="sec">'))
        check(f"{mode}: Needs attention present", "Needs attention" in html)
        check(f"{mode}: a clean run says nothing is flagged", "Nothing flagged" in html and "No issues flagged" in html)
    bare = manifest("full")
    bare["legacy_report"]["kraken_overview"] = {}
    try:
        html = render_html_text(bare, asset_mode="api")
        check("a run with no Kraken report still renders", "No Kraken report was found" in html)
    except Exception as exc:  # noqa: BLE001
        check(f"a run with no Kraken report still renders ({exc})", False)
    check("a recorded run_mode wins over inference", S.infer_run_mode({**manifest("full"), "run_mode": "parse"}) == "parse")

    print("Summary content")
    full = render_html_text(manifest("full"), asset_mode="embed", base_dir=tmp)
    check("full: headline is the strongest BLAST subject, not the file's first row",
          verdict_title(full) == "Examplebacter alpha strain T")
    check("full: share of assembled bases shown", "99.9% of assembled bases" in full)
    check("full: contig bins as three labelled numbers",
          "300 bp or shorter" in full and "<b>30</b>" in full and "3|7|30" not in full)
    konly = render_html_text(manifest("kraken_only"), asset_mode="embed", base_dir=tmp)
    check("kraken_only: headline is the most abundant taxon", verdict_title(konly) == "Examplebacter alpha complex")
    check("kraken_only: no target / assembly sections",
          "Target taxon reads" not in konly and "De novo assembly" not in konly)
    parse = render_html_text(manifest("parse"), asset_mode="embed", base_dir=tmp)
    check("parse: headline is the extracted target", "Target reads extracted" in parse
          and verdict_title(parse) == "Examplebacter alpha complex")
    check("parse: parsed read files listed", "SAMPLE-A_Examplebacter_alpha_complex_R1.fastq.gz" in parse)
    check("parse: no BLAST section", "Assembly identification" not in parse)

    print("Layout rules")
    css = full[: full.find("</style>")]
    left = css.find("table.grid thead th {")
    right = css.find("table.grid thead th.num, table.grid td.num { text-align: right")
    check("numeric headers right-align (rule after, and more specific than, the left default)", 0 < left < right)
    check("numeric headers carry the class", '<th class="num">Reads</th>' in full)
    check("with an interactive chart, the static coverage figures are print-only",
          'class="print-only coverage-static"' in full and ".print-only { display: none; }" in css
          and re.search(r"@media print \{.*?\.print-only \{ display: block; \}", css, re.S))
    no_inter = manifest("full")
    no_inter["legacy_report"]["coverage_interactive"] = ""
    html = render_html_text(no_inter, asset_mode="embed", base_dir=tmp)
    check("without one, they show on screen", 'class="coverage-static"' in html and "print-only coverage-static" not in html)
    check("coverage table shows without the interactive chart", "Reference covered</th>" in html)

    odd = manifest("full", sample_id='S<1>&"2"')
    odd["legacy_report"]["blast_rows"][1]["description"] = "NC_999999 Examplebacter <i>alpha</i> & co, complete"
    html = render_html_text(odd, asset_mode="api")
    check("manifest text is HTML-escaped", "<i>alpha</i>" not in html and "&lt;i&gt;alpha&lt;/i&gt; &amp; co" in html
          and "<dd>S&lt;1&gt;&amp;&#34;2&#34;</dd>" in html)
    check("...and CSS-escaped, not HTML-escaped, in the page header",
          'content: "Kraken ID Parse \\00B7  S1&\\"2\\""' in html)
    check("the Plotly fragment is kept as markup", '<div class="plotly-fragment">chart</div>' in html)

    print("Asset modes")
    with_pie = manifest("full")
    with_pie["artifacts"].append({"name": "bracken_pie.png", "path": "bracken_pie.png", "kind": "figure"})
    emb = render_html_text(with_pie, asset_mode="embed", base_dir=tmp)
    api = render_html_text(with_pie, asset_mode="api")
    check("embed: images are data URIs", 'src="data:image/png;base64,' in emb)
    check("api: images are served by the GUI", 'src="./file?path=bracken_pie.png&amp;inline=1"' in api)
    check("api: Krona link served by the GUI", 'href="./file?path=kraken/SAMPLE-A_2026_krona.html&amp;inline=1"' in api)
    check("embed: Krona link stays relative to the run folder", 'href="kraken/SAMPLE-A_2026_krona.html"' in emb)
    # Served by a GUI as <endpoint>?path=<run>/report.html, the relative link
    # resolved to a route that does not exist; the page rewrites it (the live
    # behaviour is checked in a browser, this pins that the pieces are there).
    check("embed: the Krona link names its run file for the rewrite",
          'data-run-file="kraken/SAMPLE-A_2026_krona.html"' in emb)
    check("embed: the rewrite script is in the page",
          'querySelectorAll("a[data-run-file]")' in emb and 'params.get("path")' in emb)

    print("Needs attention")
    bad = manifest("full", status="completed_with_warnings",
                   warnings=["Target taxon not found by Kraken: Examplebacter alpha complex"])
    lr = bad["legacy_report"]
    lr["fastq"][1]["q30"] = "61.0%"
    lr["target"]["extracted_reads"] = "No Reads Found"
    lr["blast_rows"] = []
    lr["assembly"] = {}
    lr["coverage_stats"] = []
    lr["coverage_figures"] = []
    lr["coverage_interactive"] = ""
    sm = S.build_summary(bad, beta_notice=None)
    titles = [a["title"] for a in sm["attention"]]
    check("status not completed is flagged first", titles[0].startswith("Run completed with warnings"))
    check("low Q30 names the read and the rule", any("R2 quality: 61.0%" in t for t in titles))
    check("taxon-not-found warning listed", any("not found" in t for t in titles))
    check("no duplicate 'no BLAST' flag when the target was not found", not any("No BLAST" in t for t in titles))
    check("headline says the target was not found", sm["headline"]["title"].endswith("not found"))
    html = render_html_text(bad, asset_mode="embed", base_dir=tmp)
    check("flagged run: amber verdict and item count", 'verdict lvl-warn' in html and "to check</span>" in html)

    low = manifest("full")
    low["legacy_report"]["coverage_stats"][0].update(percent_covered=70.0, mean_coverage=4.0)
    low["legacy_report"]["blast_rows"][0]["nt_base_count"] = "1,800,000"
    low["legacy_report"]["blast_rows"].insert(0, {"nt_base_count": "500,000", "contigs": "4",
                                                 "description": "NZ_XX000002 Thirdus example, contig"})
    titles = [a["title"] for a in S.build_summary(low)["attention"]]
    check("thin coverage flagged (breadth and depth)",
          any("70.0% of" in t for t in titles) and any("Mean depth 4.0X" in t for t in titles))
    check("a split BLAST call is flagged", any("Top BLAST match holds only" in t for t in titles))

    few = manifest("parse")
    few["legacy_report"]["target"].update(extracted_reads=40, extraction_rate=0.4)
    sm = S.build_summary(few)
    check("under 100 target reads is a failure", sm["level"] == "fail"
          and any(a["level"] == "fail" and "Only 40" in a["title"] for a in sm["attention"]))

    mismatch = manifest("full", status="completed_with_warnings",
                        warnings=["Only 12 reads matched target taxon 'X' in database 'k2'. This usually means..."])
    sm = S.build_summary(mismatch)
    check("taxon/database mismatch is a failure, said once",
          [a["title"] for a in sm["attention"] if a["level"] == "fail"] == ["Possible taxon / database mismatch"])

    print("Orbivirus beta")
    orbi = manifest("full")
    orbi["parameters"]["taxon"] = "Orbivirus"
    orbi["legacy_report"]["target"]["taxon"] = "Orbivirus"
    import orbivirus_beta
    orbi["warnings"] = [f"{orbivirus_beta.TITLE}: {orbivirus_beta.TEXT}"]
    html = render_html_text(orbi, asset_mode="embed", base_dir=tmp)
    check("beta banner with the shared wording", orbivirus_beta.TITLE in html and orbivirus_beta.TEXT in html)
    check("beta run is never graded clean", "Review before use" in html and "No issues flagged" not in html)
    check("every PDF page carries the beta line", f'@top-right {{ content: "{orbivirus_beta.PAGE}"' in html)
    sm = S.build_summary(orbi, beta_notice=orbivirus_beta.notice())
    check("the beta warning is not listed twice",
          sum(orbivirus_beta.TITLE in a["title"] for a in sm["attention"]) == 1)

    print("Kraken overview nesting")
    rep = tmp / "nest_reportkraken.txt"
    rep.write_text(
        "  2.00\t200\t200\tU\t0\tunclassified\n"
        " 98.00\t9800\t0\tR\t1\troot\n"
        " 96.00\t9600\t0\tG\t10\t      Genus a\n"
        " 95.00\t9500\t9000\tG1\t11\t        Group a\n"
        "  4.00\t400\t400\tS\t12\t          Species a1\n"
        "  0.23\t23\t23\tS\t13\t          Species a2\n"
        "  0.24\t24\t24\tS\t20\t      Species b\n",
        encoding="utf-8")
    names = [(t["name"], t["parent"]) for t in S.kraken_overview(rep)["top_taxa"]]
    check("a species is listed directly under its own group, not under the row before it",
          names == [("Group a", None), ("Species a1", "Group a"), ("Species a2", "Group a"),
                    ("Species b", None)])

    print("PASS" if not FAILED else f"{FAILED} FAILURE(S)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
