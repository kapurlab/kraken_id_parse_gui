"""The at-a-glance top of report.html: what was run, what the sample is, and
what needs a person's attention.

The report used to open on a paragraph of method and reach the answer — the
BLAST identification, or for a Kraken-only run nothing at all — several screens
down, under tables of input statistics. A diagnostician reads the top of a
report first and often only that, so the answer, the evidence behind it and the
reasons to doubt it now come first, computed here from the run manifest (and the
run's Kraken report when an older manifest did not record an overview).

Everything here is derived, never invented: a condition is flagged only from a
number the run recorded, and each flag says which number and which rule.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

# --- Flag rules -------------------------------------------------------------
# Each is stated in the flag's own text, so the reader sees the rule applied.
Q30_WARN_PCT = 75.0          # Illumina's MiSeq 2x250 spec is >= 75% of bases at Q30
MIN_INPUT_READS = 10_000     # fewer reads than this cannot support much of anything
MIN_TARGET_READS = 100       # kraken_id_parse.MIN_TARGET_READS: no assembly below this
LOW_EXTRACTION_PCT = 1.0     # under 1% of reads on target: a trace, not an isolate
UNCLASSIFIED_WARN_PCT = 50.0 # most reads unknown to the database
COVERAGE_BREADTH_WARN = 90.0 # % of the reference covered by reads
COVERAGE_DEPTH_WARN = 10.0   # mean depth (X)
BLAST_SPLIT_PCT = 50.0       # top BLAST subject holds under half the assembled bases

MODE_LABELS = {
    "kraken_only": "Kraken + Krona only",
    "parse": "Kraken + parse reads",
    "full": "Kraken + parse + BLAST",
}

_RANKS = {
    "U": "Unclassified", "R": "Root", "D": "Domain", "K": "Kingdom",
    "P": "Phylum", "C": "Class", "O": "Order", "F": "Family", "G": "Genus",
    "S": "Species",
}


def rank_label(code: Any) -> str:
    """Kraken2 rank code -> words. Numbered codes are ranks Kraken places
    between the named ones: G1 sits between genus and species (a species group
    or complex), S1 below species (subspecies, serovar, strain), S2 below that."""
    code = str(code or "").strip()
    if not code:
        return ""
    base, depth = code[0], code[1:]
    name = _RANKS.get(base)
    if not name:
        return code
    if not depth:
        return name
    if base == "G" and depth == "1":
        return "Species group"
    if base == "S":
        return "Subspecies / strain" if depth == "1" else f"Strain (level {depth})"
    return f"Below {name.lower()}"


def num(value: Any) -> Optional[float]:
    """A number from what the pipeline records: 1234, "1,234", "93.11%",
    "68.2X", "4,285,951 bp". None when there is none."""
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    match = re.search(r"-?\d[\d,]*(?:\.\d+)?", str(value))
    if not match:
        return None
    try:
        return float(match.group(0).replace(",", ""))
    except ValueError:
        return None


# --- Kraken report ------------------------------------------------------------
def find_kraken_report(run_dir: Optional[Path]) -> Optional[Path]:
    if not run_dir:
        return None
    try:
        run_dir = Path(run_dir)
        reports = sorted((run_dir / "kraken").glob("*_reportkraken.txt"),
                         key=lambda p: p.stat().st_mtime)
    except OSError:
        return None
    return reports[-1] if reports else None


def kraken_overview(report_path: Optional[Path], top_n: int = 10) -> Dict[str, Any]:
    """What is in the sample, from a Kraken2 report: the classified share, the
    split across domains, and the most abundant species-level taxa.

    "Species-level" is rank S plus G1 (species group/complex), because Kraken
    resolves many reads only to a complex: in an MTBC isolate ~95% of reads stop
    at "Mycobacterium tuberculosis complex" (G1) and the best rank-S row holds
    3%. Listing rank S alone would name a 3% species as the sample's top taxon.
    Rows are ordered by clade reads; `parent` names the listed row a taxon sits
    under, so the table can indent a species beneath its complex."""
    out: Dict[str, Any] = {}
    if not report_path:
        return out
    rows = []
    try:
        with open(report_path, "r", encoding="utf-8", errors="replace") as fh:
            for order, line in enumerate(fh):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                try:
                    pct = float(parts[0])
                    clade = int(parts[1])
                except ValueError:
                    continue
                raw_name = parts[5]
                rows.append({
                    "pct": pct, "reads": clade, "rank": parts[3].strip(),
                    "taxid": parts[4].strip(), "name": raw_name.strip(),
                    "indent": len(raw_name) - len(raw_name.lstrip(" ")),
                    "order": order,
                })
    except OSError:
        return out
    if not rows:
        return out

    unclassified = next((r for r in rows if r["rank"] == "U"), None)
    root = next((r for r in rows if r["rank"] == "R"), None)
    total = (unclassified["reads"] if unclassified else 0) + (root["reads"] if root else 0)
    out["total_reads"] = total or None
    out["unclassified_pct"] = unclassified["pct"] if unclassified else 0.0
    out["classified_pct"] = round(100.0 - out["unclassified_pct"], 2)

    domains = [{"name": r["name"], "pct": r["pct"], "reads": r["reads"]}
               for r in rows if r["rank"] == "D"]
    domains.sort(key=lambda d: d["reads"], reverse=True)
    other = round(out["classified_pct"] - sum(d["pct"] for d in domains), 2)
    if other >= 0.01:
        domains.append({"name": "Other classified", "pct": other, "reads": None})
    if unclassified:
        domains.append({"name": "Unclassified", "pct": unclassified["pct"],
                        "reads": unclassified["reads"]})
    out["domains"] = domains

    # Parent of each row = the nearest earlier row with a smaller indent.
    stack: List[Dict[str, Any]] = []
    for r in rows:
        while stack and stack[-1]["indent"] >= r["indent"]:
            stack.pop()
        r["_ancestors"] = [a["name"] for a in stack]
        stack.append(r)
    species = [r for r in rows if r["rank"] in ("S", "G1") and r["reads"] > 0]
    species.sort(key=lambda r: r["reads"], reverse=True)
    top = species[:top_n]
    listed = {r["name"] for r in top}
    entries = [{
        "name": r["name"], "rank": r["rank"], "rank_label": rank_label(r["rank"]),
        "reads": r["reads"], "pct": r["pct"], "taxid": r["taxid"],
        "parent": next((a for a in reversed(r["_ancestors"]) if a in listed), None),
    } for r in top]
    # Each listed taxon directly under its listed parent (both tables draw a
    # child indented beneath the row above it), groups ordered by reads. In
    # plain read order a 23-read species could land under an unrelated 24-read
    # one — M. orygis drawn as a child of Homo sapiens.
    children: Dict[Optional[str], List[Dict[str, Any]]] = {}
    for e in entries:
        children.setdefault(e["parent"], []).append(e)
    ordered: List[Dict[str, Any]] = []

    def _walk(parent: Optional[str]) -> None:
        for e in children.get(parent, []):   # already in read order
            ordered.append(e)
            _walk(e["name"])

    _walk(None)
    out["top_taxa"] = ordered
    return out


# --- Run mode -----------------------------------------------------------------
def infer_run_mode(manifest: Dict[str, Any]) -> str:
    """kraken_only / parse / full. A manifest written since run_mode existed
    says so; older ones are read from what the run produced."""
    manifest = manifest or {}
    recorded = manifest.get("run_mode") or (manifest.get("parameters") or {}).get("mode")
    if recorded in MODE_LABELS:
        return recorded
    legacy = manifest.get("legacy_report") or {}
    target = legacy.get("target") or {}
    params = manifest.get("parameters") or {}
    names = [str(a.get("name", "")) for a in manifest.get("artifacts") or []]
    if (legacy.get("blast_rows") or (legacy.get("assembly") or {}).get("contig_count")
            or legacy.get("coverage_stats") or legacy.get("coverage_figures")
            or any(n.endswith(("_blast_summary.txt", "_denovo.fasta")) for n in names)):
        return "full"
    if not (params.get("taxon") or target.get("taxon")):
        return "kraken_only"
    # A full run that stopped at assembly says so in its warnings; a parse run
    # ends cleanly with nothing past the reads.
    text = " ".join(str(w) for w in manifest.get("warnings") or []).lower()
    if "assembly" in text or "blast" in text:
        return "full"
    return "parse"


# --- BLAST --------------------------------------------------------------------
def _split_subject(description: str) -> Dict[str, str]:
    """'NC_000962 Mycobacterium tuberculosis H37Rv, complete genome' ->
    accession + organism (the text before the first comma)."""
    text = str(description or "").strip()
    accession, _, rest = text.partition(" ")
    if not re.match(r"^[A-Z]{1,6}_?[A-Z0-9]*\d[\w.]*$", accession):
        accession, rest = "", text
    organism = rest.split(",")[0].strip() or rest.strip()
    return {"accession": accession, "organism": organism, "description": text,
            "subject": rest.strip()}


def blast_table(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """BLAST summary rows, strongest first, each with its share of the
    assembled bases. The summary file is written weakest-to-strongest, so file
    order buries the identification at the bottom."""
    parsed = []
    for row in rows or []:
        bases = num(row.get("nt_base_count"))
        parsed.append({**row, **_split_subject(row.get("description")),
                       "bases": bases, "contig_n": num(row.get("contigs"))})
    parsed.sort(key=lambda r: (r["bases"] or 0), reverse=True)
    total = sum(r["bases"] or 0 for r in parsed)
    for r in parsed:
        r["share_pct"] = (100.0 * r["bases"] / total) if (total and r["bases"]) else None
    return parsed


# --- Attention ----------------------------------------------------------------
def _flag(items, level, title, detail=""):
    items.append({"level": level, "title": title, "detail": detail})


def _split_warning(text: str) -> Dict[str, str]:
    """A long pipeline warning as a short title plus its explanation."""
    text = str(text).strip()
    if len(text) > 90 and ". " in text:
        head, _, rest = text.partition(". ")
        return {"title": head + ".", "detail": rest}
    return {"title": text, "detail": ""}


def _fmt_int(v: Optional[float]) -> str:
    return f"{int(round(v)):,}" if v is not None else "?"


def build_summary(manifest: Dict[str, Any], kraken: Optional[Dict[str, Any]] = None,
                  beta_notice: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    manifest = manifest or {}
    legacy = manifest.get("legacy_report") or {}
    target = legacy.get("target") or {}
    params = manifest.get("parameters") or {}
    # Every key present, so the template can test a value without tripping on
    # a run that recorded no Kraken report at all.
    kraken = {"classified_pct": None, "unclassified_pct": None, "total_reads": None,
              "domains": [], "top_taxa": [],
              **(kraken or legacy.get("kraken_overview") or {})}
    mode = infer_run_mode(manifest)
    status = str(manifest.get("status") or "completed")
    taxon = target.get("taxon") or params.get("taxon")
    warnings = [str(w) for w in manifest.get("warnings") or [] if str(w).strip()]
    attention: List[Dict[str, str]] = []

    # Status first: it outranks every other verdict in the report.
    status_level = "ok"
    if status == "completed_with_warnings":
        status_level = "warn"
    elif status not in ("completed", "succeeded", "success"):
        status_level = "fail"
    if status_level != "ok":
        _flag(attention, status_level, "Run " + status.replace("_", " "),
              "The run did not finish cleanly; read the warnings below before using any result.")

    if beta_notice:
        # The banner above the summary carries the full wording; one line here.
        _flag(attention, "fail", beta_notice.get("title", "Beta analysis"),
              beta_notice.get("short") or beta_notice.get("text", ""))

    beta_title = (beta_notice or {}).get("title", "")
    not_found = mismatch = assembly_failed = False
    for w in warnings:
        if beta_title and w.startswith(beta_title):
            continue  # the beta flag above already says it, in full
        low = w.lower()
        if "taxon not found" in low:
            not_found = True
        if "matched target taxon" in low:
            # kraken_id_parse._taxon_db_mismatch_hint: under MIN_TARGET_READS
            # matched, the run stopped before assembly.
            mismatch = True
            _flag(attention, "fail", "Possible taxon / database mismatch", w)
            continue
        if "assembly did not complete" in low:
            assembly_failed = True
        parts = _split_warning(w)
        _flag(attention, "warn" if status_level != "fail" else "fail", parts["title"], parts["detail"])

    # Input reads.
    total_reads = num(target.get("total_input_reads")) or kraken.get("total_reads")
    if total_reads is not None and total_reads < MIN_INPUT_READS:
        _flag(attention, "warn", f"Only {_fmt_int(total_reads)} input reads",
              f"Below {MIN_INPUT_READS:,}: too few reads to support a confident result.")
    for row in legacy.get("fastq") or []:
        q30 = num(row.get("q30"))
        if row.get("filename") and q30 is not None and q30 < Q30_WARN_PCT:
            _flag(attention, "warn", f"{row.get('read') or 'Read'} quality: {q30:.1f}% of bases at Q30",
                  f"{row.get('filename')} is below {Q30_WARN_PCT:.0f}% at Q30; expect more sequencing errors.")

    # Classification.
    unclassified = kraken.get("unclassified_pct")
    if unclassified is not None and unclassified >= UNCLASSIFIED_WARN_PCT:
        db = str(params.get("kraken_db") or "the Kraken database").rstrip("/").split("/")[-1]
        _flag(attention, "warn", f"{unclassified:.1f}% of reads unclassified",
              f"Most reads matched nothing in {db}; the organism may be missing from it.")

    # Target extraction.
    extracted = num(target.get("extracted_reads"))
    rate = num(target.get("extraction_rate"))
    if mode != "kraken_only" and taxon:
        if not_found or mismatch or (extracted is not None and extracted < MIN_TARGET_READS):
            if not (not_found or mismatch):
                _flag(attention, "fail", f"Only {_fmt_int(extracted)} {taxon} reads",
                      f"Below {MIN_TARGET_READS}: too few to assemble or identify.")
        elif rate is not None and rate < LOW_EXTRACTION_PCT:
            _flag(attention, "warn", f"{rate:.2f}% of reads assigned to {taxon}",
                  f"Under {LOW_EXTRACTION_PCT:g}% of the sample: a trace of the target, not an isolate.")

    # Identification.
    blast = blast_table(legacy.get("blast_rows"))
    top_blast = blast[0] if blast else None
    if mode == "full":
        if not blast and not (not_found or mismatch or assembly_failed):
            _flag(attention, "warn", "No BLAST identification",
                  "No assembled contig produced a BLAST summary row; there is no species call.")
        elif top_blast and top_blast.get("share_pct") is not None and top_blast["share_pct"] < BLAST_SPLIT_PCT:
            _flag(attention, "warn",
                  f"Top BLAST match holds only {top_blast['share_pct']:.0f}% of assembled bases",
                  "The assembly matched several references; read the full BLAST table.")

    # Coverage.
    for row in legacy.get("coverage_stats") or []:
        breadth = num(row.get("percent_covered"))
        depth = num(row.get("mean_coverage"))
        ref = row.get("reference") or row.get("header") or "reference"
        if breadth is not None and breadth < COVERAGE_BREADTH_WARN:
            _flag(attention, "warn", f"{breadth:.1f}% of {ref} covered",
                  f"Under {COVERAGE_BREADTH_WARN:.0f}% of the reference has reads; gaps weaken the consensus.")
        if depth is not None and depth < COVERAGE_DEPTH_WARN:
            _flag(attention, "warn", f"Mean depth {depth:.1f}X on {ref}",
                  f"Under {COVERAGE_DEPTH_WARN:.0f}X; variant calls there are unreliable.")

    # The answer.
    top_taxa = kraken.get("top_taxa") or []
    coverage = (legacy.get("coverage_stats") or [None])[0]
    if mode == "full" and top_blast:
        headline = {
            "eyebrow": "Identification (BLAST)",
            "title": top_blast["organism"] or top_blast["description"],
            "detail": " · ".join(x for x in [
                top_blast.get("accession"),
                (f"{top_blast['share_pct']:.1f}% of assembled bases" if top_blast.get("share_pct") is not None else ""),
                (f"{_fmt_int(top_blast['contig_n'])} contigs" if top_blast.get("contig_n") is not None else ""),
            ] if x),
        }
    elif mode != "kraken_only" and taxon and not not_found and extracted is not None:
        headline = {
            "eyebrow": "Target reads extracted",
            "title": taxon,
            "detail": f"{_fmt_int(extracted)} reads" + (f" · {rate:.2f}% of the sample" if rate is not None else ""),
        }
    elif mode != "kraken_only" and taxon and not_found:
        headline = {
            "eyebrow": "Target taxon",
            "title": f"{taxon} not found",
            "detail": "Kraken assigned no reads to it; the most abundant taxa are below.",
        }
    elif top_taxa:
        t = top_taxa[0]
        headline = {
            "eyebrow": "Most abundant taxon (Kraken)",
            "title": t["name"],
            "detail": f"{t['pct']:.2f}% of reads · {t['rank_label'].lower()}",
        }
    else:
        headline = {"eyebrow": "Result", "title": "No classification recorded",
                    "detail": "This run recorded no Kraken report."}

    worst = "ok"
    for item in attention:
        if item["level"] == "fail":
            worst = "fail"
            break
        if item["level"] == "warn":
            worst = "warn"

    return {
        "mode": mode,
        "mode_label": MODE_LABELS[mode],
        "status": status,
        "status_label": status.replace("_", " ").capitalize(),
        "status_level": status_level,
        "level": worst,
        "taxon": taxon,
        "headline": headline,
        "total_reads": total_reads,
        "extracted": extracted,
        "extraction_rate": rate,
        "kraken": kraken,
        "blast": blast,
        "top_blast": top_blast,
        "coverage": coverage,
        "attention": attention,
    }
