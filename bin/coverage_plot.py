#!/usr/bin/env python3
"""
coverage_plot.py
================
Render per-reference coverage + variant figures for the HTML/PDF report.

Two renderers from the SAME data (coverage depth array + SNP list):

  * interactive_div()  -> a self-contained Plotly <div> for the HTML report.
    Zoom/pan, hover for SNP detail, and a legend toggle for the SNP track so a
    crowded genome can be de-cluttered on demand.

  * static_png()       -> a matplotlib PNG for the WeasyPrint PDF (static). For
    references longer than DENSITY_THRESHOLD bp the individual SNP markers are
    replaced by a binned SNP-density track, which stays legible when a large
    genome carries hundreds/thousands of SNPs.

Inputs mirror what CoverageGraphGenerator already computes:
  coverage_array : 1-D numpy array, 0-based, depth per position
  snps           : list of {'pos'(1-based), 'ref', 'alt', 'qual'}
  no_coverage    : list of (start, end) 1-based inclusive zero-coverage regions
"""

import os
import re
import numpy as np

# Only "good" (high-confidence) SNPs are counted and drawn. freebayes QUAL is the
# Phred-scaled confidence that a variant is real (QUAL 20 ≈ 99% probability), and
# it scales with the read depth supporting the call — so at modest coverage even
# solid SNPs sit well below the hundreds. We therefore use the SAME gate the
# pipeline already trusts to build its consensus: freebayes QUAL > 20
# (alignment_vcf.py: `vcffilter -f "QUAL > 20"`). This keeps genuine variants at
# any coverage while dropping low-quality noise, and the count stays consistent
# with the variants the pipeline itself accepts. Shared by the VCF parser (count)
# and the renderers (markers/title). Tune here if a data type needs a different cut.
GOOD_SNP_MIN_QUAL = 20

# References longer than this get a SNP-density track (not individual markers)
# in the STATIC (PDF) figure, and individual markers are hidden by default in the
# interactive figure — a large genome carries too many SNPs to draw as individual
# markers in the space available. Short (e.g. viral) genomes show individual SNPs.
DENSITY_THRESHOLD = 100_000
# Cap points drawn in the interactive line so the HTML stays a sane size on
# multi-Mb genomes; a min/max envelope preserves peaks/troughs for zooming.
INTERACTIVE_MAX_POINTS = 30_000
# Same idea for the static (PDF) line: plotting millions of raw points overflows
# matplotlib's Agg renderer, so draw a downsampled envelope instead.
STATIC_MAX_POINTS = 12_000
# Render the interactive line with SVG (Scatter) below this many plotted points,
# and WebGL (Scattergl) only above it. A multi-segment viral report draws many
# small graphs on one page; each WebGL graph consumes a scarce browser WebGL
# context, so too many of them makes the browser drop the oldest and blank the
# first graphs. SVG has no such limit; WebGL is reserved for the few very large
# single-genome graphs where its speed matters.
WEBGL_MIN_POINTS = 8_000

# ---- USDA / USWDS palette ---------------------------------------------------- #
# Colors follow the usda.gov Design and Brand guidance: the official USDA logo
# colors plus U.S. Web Design System color tokens (no custom colors). Every mark
# keeps >= 3:1 contrast on the white plot background (WCAG 1.4.11), and SNP alt
# bases are told apart by marker SHAPE as well as color, per the USDA guidance
# "don't use color alone to convey meaning, particularly for charts or graphs".
USDA_BLUE = '#002D72'      # USDA Dark Blue (PMS 288) - coverage line
USDA_GREEN = '#005440'     # USDA Dark Green (PMS 343) - SNP-density bars
INK = '#0B2437'            # usda.gov text color
MUTED = '#565C65'          # USWDS base-dark
GRID = '#DFE1E2'           # USWDS base-lighter
REF_LINE = '#D83933'       # USWDS secondary - 100x reference line
NO_COVERAGE = '#A9AEB1'    # USWDS base-light - zero-coverage bands
NOT_CALLED = '#936F38'     # USWDS warning-darker - "SNP calling unavailable" note
# alt base -> (color, matplotlib marker, plotly marker symbol)
_NUC_STYLE = {
    'A': ('#008817', 'o', 'circle'),          # USWDS success-dark
    'C': ('#005EA2', 's', 'square'),          # USWDS primary
    'G': ('#C05600', 'D', 'diamond'),         # USWDS accent-warm-dark
    'T': ('#B50909', '^', 'triangle-up'),     # USWDS red-60v
    'N': ('#565C65', 'X', 'x'),               # USWDS base-dark (N / other)
}
_NUC_ORDER = 'ACGTN'
# Report typefaces: Public Sans / Source Sans Pro (usda.gov and USWDS), falling
# back to Helvetica/Arial (the USDA signage guide's preferred typefaces).
PLOTLY_FONT = '"Public Sans", "Source Sans Pro", "Source Sans 3", Helvetica, Arial, sans-serif'
_PREFERRED_FONTS = ('Public Sans', 'Source Sans Pro', 'Source Sans 3',
                    'Helvetica', 'Arial', 'Liberation Sans')


def _nuc_key(nt):
    k = str(nt).upper()[:1]
    return k if k in _NUC_STYLE else 'N'


def nucleotide_color(nt):
    return _NUC_STYLE[_nuc_key(nt)][0]


def mpl_report_font():
    """First installed report typeface, resolved once so matplotlib never logs a
    'font family not found' warning for every figure."""
    global _MPL_FONT
    if _MPL_FONT is None:
        from matplotlib import font_manager
        installed = {f.name for f in font_manager.fontManager.ttflist}
        _MPL_FONT = next((f for f in _PREFERRED_FONTS if f in installed), 'DejaVu Sans')
    return _MPL_FONT


_MPL_FONT = None


def _snps_by_base(snps, length):
    """Group drawable SNPs (1 <= pos <= length) by alt base, in A C G T N order."""
    groups = {}
    for s in snps:
        if 1 <= s['pos'] <= length:
            groups.setdefault(_nuc_key(s.get('alt')), []).append(s)
    return [(k, groups[k]) for k in _NUC_ORDER if k in groups]


def _safe_id(ref_id):
    return re.sub(r'[^\w\-.]', '_', str(ref_id))


def _envelope_downsample(coverage_array, max_points):
    """Downsample to <=max_points while preserving shape via per-bin min & max.
    Returns (positions, depths) with each bin contributing its min then max
    point so peaks and dropouts survive. Small arrays are returned as-is."""
    n = len(coverage_array)
    if n <= max_points:
        pos = np.arange(1, n + 1)
        return pos, coverage_array
    bins = max_points // 2
    edges = np.linspace(0, n, bins + 1, dtype=int)
    xs, ys = [], []
    for i in range(bins):
        lo, hi = edges[i], edges[i + 1]
        if hi <= lo:
            continue
        seg = coverage_array[lo:hi]
        lo_idx = lo + int(np.argmin(seg))
        hi_idx = lo + int(np.argmax(seg))
        for idx in sorted((lo_idx, hi_idx)):
            xs.append(idx + 1)          # 1-based
            ys.append(float(coverage_array[idx]))
    return np.array(xs), np.array(ys)


_NOT_CALLED = object()  # sentinel: variant calling did not run


def _snp_count(stats, snps):
    """Reported SNP count. Prefer the per-reference count stored on the stats
    dict by the generator (authoritative, VCF-derived, high-quality only). A
    stored value of None means variant calling did not run — return the
    _NOT_CALLED sentinel so callers can say 'not called' instead of '0'. If the
    key is absent entirely, fall back to the length of the drawn SNP list."""
    if 'snp_count' not in stats:
        return len(snps or [])
    n = stats['snp_count']
    return _NOT_CALLED if n is None else int(n)


def _snp_label(snp_count):
    """Text fragment for the figure title."""
    if snp_count is _NOT_CALLED:
        return "SNPs not called"
    return f"{snp_count:,} SNPs, Q>{GOOD_SNP_MIN_QUAL}"


def _fig_title(ref_id, stats, snp_count):
    return (f"{ref_id} — Coverage & Variants ({_snp_label(snp_count)})<br>"
            f"<sub>Mean {stats['mean_coverage']:.1f}× · "
            f"{stats['length']:,} bp · {stats['percent_covered']:.1f}% covered</sub>")


def interactive_div(ref_id, stats, snps=None, no_coverage=None, include_plotlyjs=False):
    """Return an interactive Plotly coverage figure as an HTML <div> string."""
    import plotly.graph_objects as go
    from plotly.offline import plot as _offline_plot

    coverage_array = np.asarray(stats['coverage_array'], dtype=float)
    snps = snps or []
    no_coverage = no_coverage or []

    fig = go.Figure()

    # Coverage depth (downsampled envelope for large genomes). Use SVG for
    # small graphs so many-segment reports don't exhaust browser WebGL contexts;
    # reserve WebGL for the few very large single-genome graphs.
    xs, ys = _envelope_downsample(coverage_array, INTERACTIVE_MAX_POINTS)
    Trace = go.Scattergl if len(xs) > WEBGL_MIN_POINTS else go.Scatter
    fig.add_trace(Trace(
        x=xs, y=ys, mode='lines', name='Coverage',
        line=dict(color=USDA_BLUE, width=1.2),
        fill='tozeroy', fillcolor='rgba(0,45,114,0.12)',
        hovertemplate='pos %{x:,}<br>depth %{y:.0f}×<extra></extra>'))

    # 100× reference line.
    fig.add_hline(y=100, line=dict(color=REF_LINE, width=1, dash='dash'),
                  annotation_text='100× reference', annotation_position='top left',
                  annotation_font=dict(color=INK, size=11))

    # SNP markers with rich hover: one legend entry per alt base, each with its
    # own marker shape AND color so the base can be read without color vision.
    # On large genomes a marker per SNP clutters the zoomed-out view, so they
    # start hidden and the reader clicks the legend to reveal them (then zoom/pan).
    dense = len(coverage_array) > DENSITY_THRESHOLD
    for base, group in _snps_by_base(snps, len(coverage_array)):
        color, _mpl_marker, symbol = _NUC_STYLE[base]
        fig.add_trace(Trace(
            x=[s['pos'] for s in group], y=[coverage_array[s['pos'] - 1] for s in group],
            mode='markers', name=f'SNP alt {base} ({len(group):,})',
            visible='legendonly' if dense else True,
            marker=dict(color=color, size=8, symbol=symbol,
                        line=dict(color='black', width=0.6)),
            text=[f"pos {s['pos']:,}<br>{s['ref']}→{s['alt']}<br>Q {s['qual']:.0f}" for s in group],
            hovertemplate='%{text}<extra></extra>'))

    # Zero-coverage regions as light bands (plus a legend key for them).
    for start, end in no_coverage:
        fig.add_vrect(x0=start, x1=end, fillcolor='rgba(169,174,177,0.30)',
                      line_width=0)
    if no_coverage:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name='No coverage',
                                 marker=dict(symbol='square', size=11, color='rgba(169,174,177,0.6)'),
                                 hoverinfo='skip'))

    fig.update_layout(
        title=dict(text=_fig_title(ref_id, stats, _snp_count(stats, snps)),
                   font=dict(size=14, color=INK)),
        font=dict(family=PLOTLY_FONT, color=INK),
        xaxis_title='Position (bp)', yaxis_title='Coverage depth',
        template='plotly_white', height=430,
        margin=dict(l=60, r=20, t=60, b=45),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
        hovermode='closest')
    fig.update_xaxes(rangeslider=dict(visible=False), showspikes=True,
                     spikemode='across', spikethickness=1, gridcolor=GRID)
    fig.update_yaxes(gridcolor=GRID)

    return _offline_plot(fig, include_plotlyjs=include_plotlyjs, output_type='div',
                         config={'displaylogo': False, 'responsive': True,
                                 'modeBarButtonsToRemove': ['lasso2d', 'select2d']})


def static_png(ref_id, stats, snps=None, no_coverage=None, out_dir='.', dpi=150):
    """Render the static coverage PNG for the PDF. Large references (> threshold)
    get a SNP-density track instead of individual markers. Returns the path."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    coverage_array = np.asarray(stats['coverage_array'], dtype=float)
    snps = snps or []
    no_coverage = no_coverage or []
    length = len(coverage_array)
    dense = length > DENSITY_THRESHOLD

    style = {'font.family': mpl_report_font(), 'text.color': INK, 'axes.labelcolor': INK,
             'axes.edgecolor': MUTED, 'xtick.color': INK, 'ytick.color': INK,
             'axes.titlecolor': INK, 'grid.color': GRID}
    with plt.rc_context(style):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 1]})

        # Downsample the depth line for very large genomes. Plotting millions of
        # jagged points overflows the Agg renderer ("Exceeded cell block limit")
        # and bloats the PDF; a min/max envelope preserves peaks and dropouts.
        px, py = _envelope_downsample(coverage_array, STATIC_MAX_POINTS)
        ax1.plot(px, py, color=USDA_BLUE, linewidth=0.7)
        ax1.fill_between(px, py, alpha=0.15, color=USDA_BLUE)
        ax1.axhline(100, color=REF_LINE, linestyle='--', linewidth=0.9, label='100× reference')
        ax1.set_ylabel('Coverage depth')
        ax1.grid(True, alpha=0.8)
        snp_count = _snp_count(stats, snps)
        not_called = snp_count is _NOT_CALLED
        title = f'{ref_id} — Coverage & Variants ({_snp_label(snp_count)})'
        ax1.set_title(f'{title}\nMean {stats["mean_coverage"]:.1f}×, '
                      f'{stats["length"]:,} bp, {stats["percent_covered"]:.1f}% covered',
                      fontsize=10)

        # zero-coverage bands
        for i, (start, end) in enumerate(no_coverage):
            ax1.axvspan(start, end, color=NO_COVERAGE, alpha=0.35,
                        label='No coverage' if i == 0 else None)

        if not_called:
            # Variant calling did not run — say so plainly instead of an empty track
            # that reads as "0 SNPs".
            ax2.text(0.5, 0.5, 'SNP calling unavailable — variant caller did not run',
                     transform=ax2.transAxes, ha='center', va='center',
                     fontsize=9, color=NOT_CALLED, style='italic')
            ax2.set_yticks([])
            ax2.set_ylabel('Variants')
        elif dense and snps:
            # SNP-density track: binned SNP counts across the genome.
            nbins = min(300, max(50, length // 20_000))
            counts, edges = np.histogram([s['pos'] for s in snps], bins=nbins, range=(1, length))
            centers = (edges[:-1] + edges[1:]) / 2
            width = (edges[1] - edges[0])
            ax2.bar(centers, counts, width=width, color=USDA_GREEN, align='center')
            ax2.set_ylabel('SNP density')
            ax2.set_ylim(0, max(1, counts.max()) * 1.15)
            ax2.text(0.995, 0.9, f'{snp_count:,} SNPs in {nbins} bins',
                     transform=ax2.transAxes, ha='right', va='top', fontsize=7,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85))
        else:
            # Individual SNP markers (no per-SNP text labels — those caused clutter),
            # one shape + color per alt base so it reads without color vision.
            ax2.set_ylim(-0.6, 1.2)
            ax2.set_ylabel('Variants')
            rank = {id(s): i for i, s in enumerate(snps)}  # alternate track heights
            for base, group in _snps_by_base(snps, length):
                color, marker, _symbol = _NUC_STYLE[base]
                pos = np.array([s['pos'] for s in group])
                label = f'SNP alt {base} ({len(group):,})'
                ax1.scatter(pos, coverage_array[pos - 1], color=color, s=30, label=label,
                            marker=marker, zorder=5, edgecolor='black', linewidth=0.4)
                heights = [0.5 if rank[id(s)] % 2 == 0 else 0.8 for s in group]
                ax2.scatter(pos, heights, color=color, s=22,
                            marker=marker, zorder=5, edgecolor='black', linewidth=0.4)
            if not snps:
                # Calling ran and found nothing: say so rather than leave a blank track.
                ax2.text(0.5, 0.5, f'No high-confidence SNPs (QUAL > {GOOD_SNP_MIN_QUAL})',
                         transform=ax2.transAxes, ha='center', va='center',
                         fontsize=9, color=MUTED, style='italic')
                ax2.set_yticks([])

        ax2.set_xlabel('Position (bp)')
        # Plain position numbers (200,000 rather than a 1e6 axis offset).
        ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _pos: f'{int(x):,}'))
        ax2.grid(True, alpha=0.8)
        # One key below the plots (100x line, no-coverage bands, SNP bases), so it
        # can never cover data points.
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=max(1, len(handles)),
                   fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
        plt.tight_layout(rect=(0, 0.06, 1, 1))

        out_png = os.path.join(out_dir, f'{_safe_id(ref_id)}_coverage.png')
        plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
    return os.path.abspath(out_png)


# --------------------------------------------------------------------------- #
# Standalone self-test — troubleshoot rendering without a pipeline run:
#   python coverage_plot.py --selftest
# Renders a small viral segment (SVG path, individual SNP markers) and a large
# genome (WebGL path, binned SNP-density track + downsampled line) to both an
# interactive HTML and static PNGs, so every rendering branch is exercised.
# --------------------------------------------------------------------------- #
def _selftest(out_dir='coverage_plot_selftest'):
    os.makedirs(out_dir, exist_ok=True)
    # (ref_id, length, mean_depth, num_snps, zero-coverage regions)
    cases = [
        ('demo_segment2_VP2', 2_926, 500, 8, [(700, 780)]),        # small -> SVG
        ('demo_genome_large', 1_200_000, 40, 600, [(300_000, 302_000)]),  # large -> WebGL/density
    ]
    parts = []
    for i, (rid, length, mean, nsnp, drop) in enumerate(cases):
        rng = np.random.default_rng(i)
        arr = np.clip(np.abs(rng.normal(mean, mean * 0.3, length)), 0, None)
        for lo, hi in drop:
            arr[lo:hi] = 0
        stats = {'coverage_array': arr, 'length': length,
                 'mean_coverage': float(arr.mean()),
                 'percent_covered': 100.0 * float((arr > 0).mean()), 'header': rid}
        pos = sorted(rng.choice(np.arange(1, length + 1), size=min(nsnp, length), replace=False))
        # cycle the alt base so every marker shape/color is exercised
        snps = [{'pos': int(p), 'ref': 'N', 'alt': 'ACGT'[k % 4], 'qual': 500.0}
                for k, p in enumerate(pos)]
        stats['snp_count'] = len(snps)
        div = interactive_div(rid, stats, snps, drop, include_plotlyjs=(i == 0))
        png = static_png(rid, stats, snps, drop, out_dir=out_dir)
        trace = 'WebGL (Scattergl)' if length > WEBGL_MIN_POINTS else 'SVG (Scatter)'
        print(f'  {rid}: {length:,} bp -> {trace}; static PNG {png}')
        parts.append(f'<h3>{rid} <small>({trace})</small></h3>{div}'
                     f'<img src="{os.path.basename(png)}" style="max-width:760px">')
    html_path = os.path.join(out_dir, 'coverage_plot_selftest.html')
    with open(html_path, 'w') as fh:
        fh.write('<!doctype html><meta charset="utf-8"><body style="font-family:sans-serif">'
                 '<h2>coverage_plot self-test</h2>' + ''.join(parts) + '</body>')
    print(f'  interactive demo -> {html_path}')
    return html_path


def main():
    import argparse
    ap = argparse.ArgumentParser(description='Coverage graph renderer (report helper). '
                                 'Run with --selftest to render demo graphs for troubleshooting.')
    ap.add_argument('--selftest', action='store_true',
                    help='Render demo coverage graphs (small + large) to verify rendering')
    ap.add_argument('-o', '--out', default='coverage_plot_selftest', help='Output directory')
    args = ap.parse_args()
    if args.selftest:
        print('coverage_plot self-test:')
        _selftest(args.out)
    else:
        ap.print_help()


if __name__ == '__main__':
    main()
