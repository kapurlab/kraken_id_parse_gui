#!/usr/bin/env python
"""
Interactive Coverage & Variants charts for the HTML report (Plotly).

The report has always carried a static coverage PNG: you could see that depth
dipped somewhere, but not where, not how far, and not what the reads actually
said there. This builds the interactive version — one zoomable chart per
reference, with:

  * per-position depth (binned for genome-scale references, see BIN_TARGET),
  * a diamond at every variant position,
  * directly beneath each diamond, a stacked bar coloured by the nucleotides
    observed at that position, in proportion — so a 50/50 G/A call reads as a
    half-yellow, half-green bar (A green, C blue, G yellow, T red),
  * zero-coverage stretches shaded, and the 100X depth threshold drawn in.

Allele counts come from pysam's count_coverage(), which applies the base
quality filter in C over the whole contig — fast enough to run on every
reference we chart, and the reason the titles can honestly say Q>20.

Everything is best-effort: without plotly (or on any failure) the caller keeps
the static PNG it already had.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Nucleotide colours — the conventional scheme shared with the IRMA report, so
# a mixed call looks the same in every tool in the suite.
NT_COLORS = {"A": "#3AA655", "C": "#3B6FD4", "G": "#E8B93B",
             "T": "#D8453E", "-": "#8B96A0", "N": "#B9C0C6"}
NT_ORDER = ["A", "C", "G", "T"]

TEAL_DARK = "#2F6F6C"
INK = "#1F2A2E"
DANGER = "#C46A6A"

MIN_BASE_QUALITY = 20      # the "Q>20" the chart titles quote
DEPTH_THRESHOLD = 100      # the 100X line the tool's coverage graphs have always drawn
MIN_VARIANT_DEPTH = 10     # below this a "variant" is noise, not a call
MIN_MINOR_FRACTION = 0.10  # a mixed call has to be at least this mixed to plot
BIN_TARGET = 6000          # plotted points per reference (a 4.4 Mb genome would
                           # otherwise put 4.4 M points in the browser)
MAX_REFERENCES = 6         # charts per report, best-covered first
MAX_VARIANTS_PER_REF = 3000


def _bin_depths(depths: List[float], target: int = BIN_TARGET):
    """(x, y) for plotting: exact positions when short, else per-bin MINIMUM.

    The minimum, not the mean: this chart's job is to show where coverage
    fails, and averaging a 200 bp dropout into a 700 bp bin hides exactly the
    thing the reader is looking for."""
    n = len(depths)
    if n <= target:
        return list(range(1, n + 1)), list(depths)
    step = (n + target - 1) // target
    xs, ys = [], []
    for start in range(0, n, step):
        chunk = depths[start:start + step]
        if not chunk:
            continue
        xs.append(start + 1)
        ys.append(min(chunk))
    return xs, ys


def _zero_runs(depths: List[float]) -> List[Tuple[int, int]]:
    runs: List[Tuple[int, int]] = []
    start = None
    for i, d in enumerate(depths, start=1):
        if d <= 0 and start is None:
            start = i
        elif d > 0 and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(depths)))
    return runs


def allele_stem_traces(go, variants, depths, line_width: float = 3.0):
    """Stacked allele "stems" drawn ON the coverage curve, one trace per base.

    Each variant gets a vertical line at its position running from 0 up to the
    read depth there — the height the diamond sits at — split into coloured
    lengths by the nucleotides observed. A 50/50 G/A call is a stem that is
    half yellow, half green.

    Scatter rather than Bar: scatter shares the layer the coverage curve is in,
    so the stems sit ON TOP of the area fill instead of behind it, and the
    width is in PIXELS — a hairline stays a hairline whether the reference is
    1 kb or 4 Mb, and zooming does not fatten it. One trace per base keeps the
    legend to four entries no matter how many variants there are.
    """
    by_base: Dict[str, Dict[str, List]] = {}
    for v in variants:
        pos = v["position"]
        depth = depths[pos - 1] if 0 < pos <= len(depths) else 0
        if depth <= 0:
            continue
        bottom = 0.0
        for allele, frac, count in v["composition"]:
            top = bottom + frac * depth
            key = allele if allele in NT_COLORS else "other"
            entry = by_base.setdefault(key, {"x": [], "y": [], "text": []})
            # None breaks the polyline so segments stay separate strokes.
            entry["x"] += [pos, pos, None]
            entry["y"] += [bottom, top, None]
            label = f"{allele}: {frac * 100:.1f}%" + (f" ({count:,} reads)" if count else "")
            entry["text"] += [label, label, None]
            bottom = top
    traces = []
    ordered = [k for k in NT_ORDER if k in by_base] + \
              [k for k in by_base if k not in NT_ORDER]
    for key in ordered:
        e = by_base[key]
        traces.append(go.Scatter(
            x=e["x"], y=e["y"], mode="lines", name=key,
            line=dict(color=NT_COLORS.get(key, "#8B96A0"), width=line_width),
            text=e["text"], hovertemplate="pos %{x:,} — %{text}<extra></extra>",
            connectgaps=False, showlegend=True,
        ))
    return traces


def call_variants(bam_path: str, ref_id: str, ref_seq: str,
                  depths: Optional[List[float]] = None) -> List[Dict[str, Any]]:
    """Variant positions on one reference, with their full allele composition.

    A position is reported when it has usable depth and either the majority
    base disagrees with the reference, or a second base is present above
    MIN_MINOR_FRACTION (the mixed call this chart exists to show)."""
    import pysam
    out: List[Dict[str, Any]] = []
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            counts = bam.count_coverage(
                ref_id, quality_threshold=MIN_BASE_QUALITY, read_callback="all")
    except (ValueError, OSError):
        return out
    a_arr, c_arr, g_arr, t_arr = counts
    n = min(len(a_arr), len(ref_seq))
    for i in range(n):
        comp = ((a_arr[i], "A"), (c_arr[i], "C"), (g_arr[i], "G"), (t_arr[i], "T"))
        total = comp[0][0] + comp[1][0] + comp[2][0] + comp[3][0]
        if total < MIN_VARIANT_DEPTH:
            continue
        ordered = sorted(comp, reverse=True)
        top_count, top_base = ordered[0]
        second_count, _second_base = ordered[1]
        ref_base = ref_seq[i].upper()
        minor_frac = second_count / total
        if top_base == ref_base and minor_frac < MIN_MINOR_FRACTION:
            continue
        out.append({
            "position": i + 1,
            "ref": ref_base,
            "consensus": top_base,
            "depth": total,
            "minor_fraction": minor_frac,
            "composition": [(b, cnt / total, cnt) for cnt, b in ordered if cnt > 0],
        })
        if len(out) >= MAX_VARIANTS_PER_REF:
            break
    return out


def _figure(sample: str, ref_id: str, header: str, depths: List[float],
            variants: List[Dict[str, Any]], guided_by: str = ""):
    """One reference's Coverage & Variants figure.

    A single panel: the depth curve, a diamond at every variant, and under each
    diamond a stem down to zero coloured by the nucleotides observed there in
    proportion.
    """
    try:
        import plotly.graph_objects as go
    except Exception:  # noqa: BLE001
        return None
    total_len = len(depths)
    if not total_len:
        return None
    covered = sum(1 for d in depths if d > 0)
    pct_cov = covered / total_len * 100
    mean = sum(depths) / total_len

    fig = go.Figure()
    xs, ys = _bin_depths(depths)
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="lines", line=dict(color=TEAL_DARK, width=1),
        fill="tozeroy", fillcolor="rgba(76,140,138,0.30)",
        name="depth", showlegend=False,
        hovertemplate="pos %{x:,}: %{y:.0f}×<extra></extra>",
    ))
    for a, b in _zero_runs(depths):
        if (b - a) < max(total_len * 0.0005, 1):
            continue  # single-base dropouts would just stipple the chart
        if (b - a) > total_len * 0.02:
            fig.add_vrect(x0=a, x1=b, fillcolor="rgba(196,106,106,0.20)", line_width=0,
                          annotation_text="no coverage", annotation_position="top left",
                          annotation_font_size=10, annotation_font_color=DANGER)
        else:
            fig.add_vrect(x0=a, x1=b, fillcolor="rgba(196,106,106,0.20)", line_width=0)
    fig.add_hline(y=DEPTH_THRESHOLD, line_dash="dash", line_color=DANGER,
                  annotation_text=f"{DEPTH_THRESHOLD}×", annotation_position="top left",
                  annotation_font_color=DANGER)

    if variants:
        # Stems first so the diamonds stay on top of them.
        for tr in allele_stem_traces(go, variants, depths):
            fig.add_trace(tr)
        vx = [v["position"] for v in variants]
        vy = [depths[v["position"] - 1] if v["position"] <= total_len else 0 for v in variants]
        vtext = [
            f"{v['ref']}→{v['consensus']} · depth {v['depth']:,}× · "
            + ", ".join(f"{b} {f * 100:.0f}%" for b, f, _c in v["composition"])
            for v in variants
        ]
        fig.add_trace(go.Scatter(
            x=vx, y=vy, mode="markers", name="variants", showlegend=False,
            marker=dict(color=DANGER, size=9, symbol="diamond",
                        line=dict(color="#7E3B3B", width=1)),
            text=vtext, hovertemplate="pos %{x:,}: %{text}<extra>SNP</extra>",
        ))

    guided = f"<br><sup>{guided_by}</sup>" if guided_by else ""
    fig.update_layout(
        title=dict(text=(f"{sample} — {ref_id} — Coverage &amp; Variants "
                         f"({len(variants)} SNPs, Q&gt;{MIN_BASE_QUALITY})"
                         f"<br><sup>Mean {mean:.1f}× · {total_len:,} bp · "
                         f"{pct_cov:.1f}% covered</sup>{guided}"),
                   font=dict(size=15, color=INK)),
        template="plotly_white",
        height=380,
        margin=dict(l=70, r=24, t=96, b=48),
        xaxis_title="Position (bp)",
        yaxis_title="Coverage depth",
        showlegend=bool(variants),
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0,
                    font=dict(size=11), title=dict(text="SNP alleles  ", side="left")),
    )
    return fig


def build_charts(sample: str, sorted_bam: str, ref_seqs: Dict[str, str],
                 depth_by_ref: Dict[str, List[float]],
                 headers: Dict[str, str], outpath: Path, log=print) -> Optional[str]:
    """Write one HTML fragment holding a chart per covered reference.

    Returns the path written, or None when nothing could be built — and says
    WHY on the way out. A silent None is indistinguishable from "this run had
    nothing to chart", which is how a missing dependency hides for months.
    """
    try:
        import plotly  # noqa: F401
    except Exception:  # noqa: BLE001
        log("  interactive coverage chart skipped: plotly is not installed in "
            "this environment (the static coverage figure is unaffected).")
        return None
    ranked = sorted(
        ((r, d) for r, d in depth_by_ref.items() if d and sum(d) > 0),
        key=lambda kv: sum(kv[1]) / max(len(kv[1]), 1), reverse=True,
    )[:MAX_REFERENCES]
    if not ranked:
        log("  interactive coverage chart skipped: no reference had any coverage.")
        return None
    parts: List[str] = []
    included_js = False
    for ref_id, depths in ranked:
        try:
            variants = call_variants(sorted_bam, ref_id, ref_seqs.get(ref_id, ""), depths)
        except Exception:  # noqa: BLE001
            variants = []
        header = headers.get(ref_id, ref_id)
        guided = ""
        # The pipeline underscores the first four spaces of each consensus
        # header (clean_fasta_headers), so the marker arrives as "guided_by".
        if "guided by" in header or "guided_by" in header:
            # Restore only the marker's own underscores. A blanket replace would
            # also rewrite NC_000962.3 into "NC 000962.3" — corrupting the very
            # accession the caption exists to name.
            guided = ("Coverage analysis with SNPs for: "
                      + header.replace("_guided_by_", " guided by "))
        fig = _figure(sample, ref_id, header, depths, variants, guided)
        if fig is None:
            continue
        parts.append(fig.to_html(full_html=False,
                                 include_plotlyjs=(not included_js),
                                 config={"responsive": True, "displaylogo": False}))
        included_js = True
    if not parts:
        log("  interactive coverage chart skipped: no figure could be built.")
        return None
    Path(outpath).write_text("\n".join(parts), encoding="utf-8")
    return str(outpath)
