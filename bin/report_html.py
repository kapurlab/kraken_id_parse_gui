#!/usr/bin/env python3
"""
report_html.py
==============
HTML-first report builder for the kraken_id_parse pipeline.

Produces a professional, self-contained HTML report (with interactive Plotly
coverage graphs) and renders a matching PDF from the same content via WeasyPrint
(no TeX/LaTeX toolchain required).

Usage sketch:
    r = HtmlReport(sample_name, logo=logo_path, out_dir='.')
    r.add_fastq_quality(r1, r2)
    r.add_pie(pie_png_path, 'FASTQ Read Identification')
    r.add_kv('Kraken Classification', [('Target Taxon', taxon), ...])
    r.add_table('Kraken Detailed Classification',
                ['Taxonomic Name','Level','Reads','Percent'], rows)
    r.add_assembly(assembler)
    r.add_table('BLAST nt_viruses - Assembly Identification',
                ['nt base count','contigs','Description'], rows, footer='...')
    r.add_serotype(consensus, interpretation, predictions)     # optional
    r.add_alignment_stats(stats_rows)                          # Table 1
    r.add_coverage('Coverage Graphs', items)                   # interactive+static
    html_path = r.write_html()
    pdf_path = r.write_pdf()
"""

import base64
import html as _html
import os
from datetime import datetime

from jinja2 import Environment, BaseLoader


# --------------------------------------------------------------------------- #
# Styling
# --------------------------------------------------------------------------- #
# Follows the USDA visual standards (usda.gov USDA Style Guide > Logo, Digital
# Strategy > Design and Brand plays, Facility Signage Guide, DR 1430-002):
#   * The USDA logo / signature lockup is the only identifier, shown unmodified at
#     the top left on a solid white Signature Iso-Bar (>= 8% of the page height,
#     min 0.75 in) with clear space >= the width of the logo's letter "A" - never
#     distorted, recolored, boxed/framed, shadowed, or placed on a gradient.
#   * Colors: the official USDA Dark Blue #002D72 (PMS 288) and Dark Green
#     #005440 (PMS 343); everything else is a U.S. Web Design System color token
#     or usda.gov's own neutrals (ink #0B2437, tan #F7F1E9, gold #CEA467 - gold is
#     decorative only, never text). All text meets >= 4.5:1 contrast (Section 508).
#   * Type: Public Sans / Source Sans Pro (usda.gov, USWDS) with Helvetica/Arial
#     (preferred by the USDA signage guide) as fallbacks.
BASE_CSS = """
:root {
  --brand: {{ brand_rgb }};   /* section headings: USDA Dark Blue unless overridden */
  --usda-blue: #002D72;       /* USDA Dark Blue, PMS 288 */
  --usda-green: #005440;      /* USDA Dark Green, PMS 343 */
  --usda-gold: #CEA467;       /* usda.gov accent - decorative rules only */
  --ink: #0B2437;             /* usda.gov text */
  --muted: #565C65;           /* USWDS base-dark */
  --line: #DFE1E2;            /* USWDS base-lighter */
  --zebra: #F0F0F0;           /* USWDS base-lightest */
  --tan: #F7F1E9;             /* usda.gov warm background */
  --on-dark: #DBE3EF;         /* secondary text on USDA Dark Blue (10:1) */
}
* { box-sizing: border-box; }
body {
  font-family: "Public Sans Web", "Public Sans", "Source Sans Pro", "Source Sans 3",
               Helvetica, Arial, sans-serif;
  color: var(--ink); margin: 0; font-size: 13px; line-height: 1.5;
  background: var(--zebra);
  -webkit-font-smoothing: antialiased;
}
.page { max-width: 980px; margin: 0 auto; background: #fff; }

/* ---- Signature Iso-Bar: solid white band holding only the USDA lockup ---- */
.iso-bar {
  display: flex; align-items: center; background: #fff;
  min-height: 92px; padding: 20px 34px;   /* clear space >= the logo's letter "A" */
}
.iso-bar .lockup {
  display: block; height: 52px; width: auto; max-width: 70%;  /* proportions preserved */
}

/* ---- title band (below the Iso-Bar) ------------------------------------ */
.report-head {
  display: flex; align-items: flex-end; justify-content: space-between; gap: 18px;
  padding: 18px 34px 16px; background: var(--usda-blue); color: #fff;
  border-bottom: 4px solid var(--usda-gold);
}
.report-head .titles { min-width: 0; }
.report-head .kicker {
  text-transform: uppercase; letter-spacing: .12em; font-size: 11px;
  font-weight: 700; margin: 0 0 4px; color: var(--on-dark);
}
.report-head h1.sample {
  font-size: 26px; font-weight: 700; margin: 0; line-height: 1.15;
  overflow-wrap: anywhere;
}
.report-head .meta { text-align: right; flex: none; }
.report-head .meta .date { font-weight: 700; font-size: 14px; }
.report-head .meta .tool { font-size: 11px; margin-top: 3px; color: var(--on-dark); }

.content { padding: 22px 34px 40px; }

/* ---- sections ------------------------------------------------------------ */
.section {
  margin-top: 18px; page-break-inside: avoid;
  background: #fff; border: 1px solid var(--line); border-radius: 4px;
  overflow: hidden;
}
.section > h2 {
  margin: 0; font-size: 13px; font-weight: 700; letter-spacing: .05em;
  text-transform: uppercase; color: #fff; padding: 8px 16px;
  background: rgb(var(--brand));
}
.section > .sec-body { padding: 14px 16px 16px; }

table.data { border-collapse: collapse; width: 100%; font-size: 12px; }
table.data th, table.data td {
  border: 1px solid var(--line); padding: 6px 10px; text-align: left; vertical-align: top;
}
table.data thead th {
  background: var(--line); font-weight: 700; color: var(--ink);
  border-bottom: 2px solid var(--usda-blue);
}
table.data td.num, table.data th.num { text-align: right; font-variant-numeric: tabular-nums; }
table.data tbody tr:nth-child(even) { background: var(--zebra); }
/* compact: tighter rows for long tables (Kraken classification, BLAST hits) */
table.data.compact { font-size: 11px; }
table.data.compact th { padding: 4px 10px; }
table.data.compact td { padding: 2px 10px; line-height: 1.3; }
.kv th { font-weight: 700; width: 34%; color: var(--ink); background: var(--zebra); }
.note { color: var(--muted); font-size: 11px; margin: 8px 2px 0; }
.figure { text-align: center; margin: 6px 0 2px; page-break-inside: avoid; }
.figure img { max-width: 100%; height: auto; }
.figure .cap { color: var(--muted); font-size: 11px; margin-top: 4px; }
.cov-item { margin-top: 14px; page-break-inside: avoid; }
.cov-item:first-child { margin-top: 0; }
.cov-item .cap {
  font-weight: 700; font-size: 12px; margin: 0 0 6px; color: var(--ink);
  padding-bottom: 5px; border-bottom: 1px solid var(--line);
}
.status { font-size: 12px; color: var(--ink); font-weight: 600; margin: 0 0 10px; }
.serotype-call { text-align: center; margin: 2px 0 14px; }
.serotype-call .label {
  display: block; font-size: 11px; letter-spacing: .1em; text-transform: uppercase;
  color: var(--muted); margin-bottom: 6px; font-weight: 700;
}
.serotype-call .badge {
  display: inline-block; font-size: 22px; font-weight: 700; color: #fff;
  padding: 8px 26px; border-radius: 4px; background: var(--usda-green);
}
.interp {
  background: var(--tan); border-left: 4px solid var(--usda-gold);
  padding: 9px 13px; font-size: 12px; margin-top: 10px;
}
/* notices use the USWDS state palettes; the level is also spelled out in text */
.notice {
  border-left: 5px solid var(--notice-border, #E5A000); background: var(--notice-bg, #FAF3D1);
  padding: 10px 14px; font-size: 13px; white-space: pre-wrap; color: var(--ink);
}
.notice .level {
  display: block; font-size: 11px; font-weight: 700; letter-spacing: .08em;
  text-transform: uppercase; margin-bottom: 2px; white-space: normal;
}

/* ---- overview stat strip ------------------------------------------------ */
.overview { display: flex; flex-wrap: wrap; gap: 10px; padding: 14px 16px; }
.stat {
  flex: 1 1 130px; background: var(--tan); border-top: 3px solid var(--usda-green);
  padding: 9px 12px;
}
.stat .stat-label {
  font-size: 10.5px; letter-spacing: .07em; text-transform: uppercase;
  color: var(--muted); font-weight: 700; margin-bottom: 3px;
}
.stat .stat-value {
  font-size: 15px; font-weight: 700; color: var(--usda-blue);
  overflow-wrap: break-word; word-break: break-word; hyphens: auto;
}

@media print {
  body { background: #fff; font-size: 10.5px; }
  .page { max-width: none; }
  /* Iso-Bar >= 8% of an 11 in page and >= 0.75 in; lockup clear space kept */
  .iso-bar { min-height: 0.9in; padding: 0.18in 0.05in; }
  .iso-bar .lockup { height: 0.55in; }
  .report-head { padding: 12px 14px 10px; }
  .content { padding: 0 2px 8px; }
  .section { margin-top: 12px; }
  a[href]:after { content: ""; }
}
@page {
  size: letter; margin: 16mm 12mm 16mm 12mm;
  font-family: "Public Sans", "Source Sans Pro", "Source Sans 3", Helvetica, Arial, sans-serif;
  @top-right { content: "{{ date_long }}"; font-size: 9px; color: #565C65; }
  @bottom-left { content: "{{ sample_css }} \\2014  Pathogen Identification Report"; font-size: 9px; color: #565C65; }
  @bottom-right { content: "Page " counter(page) " of " counter(pages); font-size: 9px; color: #565C65; }
}
/* first page: the Iso-Bar sits at the top of the sheet and the date is in the title band */
@page :first {
  margin-top: 8mm;
  @top-right { content: none; }
}
"""

SHELL = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="generator" content="Kraken ID &amp; Parse">
<title>{{ sample }} &ndash; Pathogen Identification Report</title>
{% if mode == 'html' and plotly_js %}<script>{{ plotly_js }}</script>{% endif %}
<style>{{ css }}</style>
</head>
<body>
<div class="page">
  {% if logo_uri %}<header class="iso-bar"><img class="lockup" src="{{ logo_uri }}" alt="{{ logo_alt }}"></header>{% endif %}
  <div class="report-head">
    <div class="titles">
      <p class="kicker">Pathogen Identification Report &middot; Kraken ID &amp; Parse</p>
      <h1 class="sample">{{ sample }}</h1>
    </div>
    <div class="meta"><div class="date">{{ date_long }}</div><div class="tool">Automated taxonomic ID, assembly &amp; consensus</div></div>
  </div>
  <main class="content">
    {{ body }}
  </main>
</div>
</body>
</html>
"""


USDA_DARK_BLUE_RGB = '0,45,114'   # #002D72, PMS 288 (official USDA logo color)


def _rgb(color):
    """Accept 'r,g,b' or 'r, g, b' -> 'r,g,b' string usable in rgb()."""
    try:
        parts = [int(float(x)) for x in str(color).replace(' ', '').split(',')[:3]]
        if len(parts) == 3:
            return ','.join(str(p) for p in parts)
    except Exception:
        pass
    return USDA_DARK_BLUE_RGB


def _default_logo_alt(logo):
    """Meaningful alt text for the header logo (Section 508)."""
    if logo and 'usda' in os.path.basename(logo).lower():
        return 'U.S. Department of Agriculture (USDA) logo'
    return 'Organization logo'


def _css_string(text):
    """Escape text for use inside a CSS content: "..." string."""
    return str(text).replace('\\', '\\\\').replace('"', '\\"')


def _file_uri(path):
    return 'file://' + os.path.abspath(path)


def _data_uri(path):
    ext = os.path.splitext(path)[1].lower().lstrip('.') or 'png'
    mime = {'jpg': 'jpeg', 'svg': 'svg+xml'}.get(ext, ext)
    with open(path, 'rb') as fh:
        b64 = base64.b64encode(fh.read()).decode('ascii')
    return f'data:image/{mime};base64,{b64}'


class HtmlReport:
    def __init__(self, sample_name, logo=None, out_dir='.', brand_color=USDA_DARK_BLUE_RGB,
                 logo_alt=None):
        """logo: for USDA reports pass the official USDA signature lockup artwork
        supplied by the Office of Communications (PNG or SVG). It is placed
        unmodified at the top left of a white Iso-Bar, as the USDA style guide
        requires. brand_color only tints the section headings and defaults to the
        official USDA Dark Blue."""
        self.sample = sample_name
        self.logo = logo if (logo and os.path.exists(logo)) else None
        self.logo_alt = logo_alt or _default_logo_alt(self.logo)
        self.out_dir = out_dir
        self.brand = _rgb(brand_color)
        now = datetime.now()
        self.date_long = now.strftime('%B %d, %Y')
        self.stamp = now.strftime('%Y-%m-%d_%H-%M-%S')
        self.sections = []            # list of html-string blocks (order preserved)
        self._cov_slots = []          # coverage items, swapped per render mode
        # autoescape OFF by design: we inject raw HTML/CSS/JS blocks and escape
        # all dynamic *data* ourselves via self._esc().
        self._env = Environment(loader=BaseLoader(), autoescape=False)

    # ---- low-level helpers ------------------------------------------------- #
    @staticmethod
    def _esc(v):
        return _html.escape('' if v is None else str(v))

    def _section(self, title, inner):
        self.sections.append(
            f'<div class="section"><h2>{self._esc(title)}</h2>'
            f'<div class="sec-body">{inner}</div></div>')

    # ---- section builders -------------------------------------------------- #
    def add_overview(self, rows):
        """A glanceable strip of stat cards, inserted as the very first section
        regardless of call order (run summary belongs right under the header).
        rows = [(label, value), ...]."""
        cards = ''.join(
            f'<div class="stat"><div class="stat-label">{self._esc(k)}</div>'
            f'<div class="stat-value">{self._esc(v)}</div></div>' for k, v in rows)
        block = (f'<div class="section"><h2>Run Overview</h2>'
                 f'<div class="overview">{cards}</div></div>')
        self.sections.insert(0, block)

    # USWDS state palettes: (label, border, background). The level is written out
    # in the notice as well, so it never relies on color alone.
    _NOTICE_STYLES = {
        'info': ('Note', '#009EC1', '#E7F6F8'),       # info-dark / info-lighter
        'warn': ('Warning', '#E5A000', '#FAF3D1'),    # warning-dark / warning-lighter
        'error': ('Error', '#B50909', '#F4E3DB'),     # error-dark / error-lighter
    }

    def add_message(self, title, text, level='warn'):
        """A prominent notice block (e.g. target taxon not found / assembly failed)."""
        label, border, bg = self._NOTICE_STYLES.get(level, self._NOTICE_STYLES['warn'])
        inner = (f'<div class="notice" style="--notice-border:{border};--notice-bg:{bg};">'
                 f'<span class="level">{label}</span>{self._esc(text)}</div>')
        self._section(title, inner)

    def add_kv(self, title, rows):
        """Two-column key/value table. rows = [(key, value), ...]."""
        body = ['<table class="data kv"><tbody>']
        for k, v in rows:
            body.append(f'<tr><th scope="row">{self._esc(k)}</th><td>{self._esc(v)}</td></tr>')
        body.append('</tbody></table>')
        self._section(title, ''.join(body))

    def add_table(self, title, headers, rows, footer=None, num_cols=None, compact=False):
        """Generic table. num_cols = set/list of column indices to right-align.
        compact=True tightens row height for long tables (still legible)."""
        num_cols = set(num_cols or [])
        cls = 'data compact' if compact else 'data'
        out = [f'<table class="{cls}"><thead><tr>']
        for i, h in enumerate(headers):
            cls = ' class="num"' if i in num_cols else ''
            out.append(f'<th{cls}>{self._esc(h)}</th>')
        out.append('</tr></thead><tbody>')
        for row in rows:
            out.append('<tr>')
            for i, cell in enumerate(row):
                cls = ' class="num"' if i in num_cols else ''
                out.append(f'<td{cls}>{self._esc(cell)}</td>')
            out.append('</tr>')
        out.append('</tbody></table>')
        if footer:
            out.append(f'<div class="note">{self._esc(footer)}</div>')
        self._section(title, ''.join(out))

    def add_fastq_quality(self, r1, r2=None):
        """r1/r2 expose: file_name, file_size, passQ30, read_quality_average, avg_len."""
        def col(c):
            return [os.path.basename(c.file_name), c.file_size, f'{c.passQ30}%',
                    f'{float(c.read_quality_average):.1f}', c.avg_len]
        labels = ['Filename', 'File Size', 'Q30 Passing', 'Mean Read Score', 'Average Read Length']
        c1 = col(r1)
        c2 = col(r2) if r2 else None
        headers = ['Metric', 'R1'] + (['R2'] if c2 else [])
        rows = []
        for i, lab in enumerate(labels):
            row = [lab, c1[i]] + ([c2[i]] if c2 else [])
            rows.append(row)
        self.add_table('FASTQ Quality', headers, rows)

    def add_pie(self, img_path, caption=None):
        """caption is optional — the chart PNG already carries its own title,
        so a caption is only added when the caller wants extra context."""
        if not img_path or not os.path.exists(img_path):
            return
        uri = _data_uri(img_path)
        cap = f'<div class="cap">{self._esc(caption)}</div>' if caption else ''
        inner = (f'<div class="figure"><img src="{uri}" '
                 f'alt="Pie chart of read identifications by taxon (Kraken/Bracken)">{cap}</div>')
        self._section('FASTQ Identifications', inner)

    def add_assembly(self, a):
        headers = ['Contig count', 'Contigs <301 | 301-999 | >999bp', 'Longest contig',
                   'Total length', 'N50', a.coverage_title]
        rows = [[f'{a.contig_count:,}',
                 f'{a.small_contigs_count:,} | {a.mid_size:,} | {a.greater_one_kb_count:,}',
                 f'{a.longest_contig:,}', f'{a.total_contig_lengths:,}', f'{a.n50:,}',
                 f'{a.mean_coverage:,.1f}X']]
        self.add_table('Assembly', headers, rows, num_cols={0, 2, 3, 4, 5})

    def add_serotype(self, consensus, interpretation, predictions, tentative_segments=None):
        tentative_segments = tentative_segments or set()
        call = self._esc(consensus)
        headers = ['Segment', 'Protein', 'Serotype', 'Top Hit', '% Identity',
                   'Query Cov', 'Bitscore', 'E-value']
        seg = {'VP2': '2', 'VP5': '6'}
        rows = []
        any_tentative = False
        for p in predictions:
            seg_str = seg.get(p.get('Protein'), '-')
            if seg_str != '-' and int(seg_str) in tentative_segments:
                seg_str = f'{seg_str} (tentative)'
                any_tentative = True
            rows.append([seg_str, p.get('Protein', '-'),
                         p.get('Serotype', '-'), p.get('Top Hit Accession', '-'),
                         p.get('Percent Identity', '-'), p.get('Query Coverage', '-'),
                         p.get('Bitscore', '-'), p.get('E-value', '-')])
        tbl = ['<div class="serotype-call"><span class="label">Predicted Serotype</span>'
               f'<span class="badge">{call}</span></div>']
        tbl.append('<table class="data"><thead><tr>')
        for i, h in enumerate(headers):
            cls = ' class="num"' if i >= 4 else ''
            tbl.append(f'<th{cls}>{self._esc(h)}</th>')
        tbl.append('</tr></thead><tbody>')
        for row in rows:
            tbl.append('<tr>')
            for i, c in enumerate(row):
                cls = ' class="num"' if i >= 4 else ''
                tbl.append(f'<td{cls}>{self._esc(c)}</td>')
            tbl.append('</tr>')
        tbl.append('</tbody></table>')
        if any_tentative:
            tbl.append('<div class="note">(tentative) = segment assigned by the '
                       'low-confidence fallback (loose match); verify manually.</div>')
        tbl.append(f'<div class="interp"><b>Interpretation:</b> {self._esc(interpretation)}</div>')
        self._section('BTV Serotyping', ''.join(tbl))

    def add_alignment_stats(self, rows, tentative_note=False):
        """rows = [(reference, length, mean_cov, pct_covered[, snp_count]), ...].
        The optional 5th element is the high-quality SNP count vs the reference.
        tentative_note appends a footnote explaining the '(tentative)' badge when
        any reference's segment was assigned by the low-confidence fallback."""
        try:
            from coverage_plot import GOOD_SNP_MIN_QUAL as min_qual
        except Exception:
            min_qual = 20
        headers = ['Reference', 'Length', 'Mean Coverage', '% Genome Covered', f'SNPs (Q>{min_qual})']
        disp = [[r[0], f'{int(r[1]):,}', f'{float(r[2]):.1f}X', f'{float(r[3]):.1f}%',
                 f'{int(r[4]):,}' if len(r) > 4 and r[4] is not None else '—'] for r in rows]
        footer = ('Table 1: Alignment statistics for each reference sequence. '
                  f'SNPs = high-confidence variants (freebayes QUAL > {min_qual}) vs the reference.')
        if tentative_note:
            footer += (' (tentative) = segment assigned by the low-confidence fallback '
                       '(loose match); verify manually.')
        self.add_table('Alignment Statistics', headers, disp, num_cols={1, 2, 3, 4},
                       footer=footer)

    def add_coverage(self, title, items, status=None):
        """items = [{'caption':.., 'div':interactive_html, 'img':static_png_path}]."""
        inner = []
        if status:
            inner.append(f'<div class="status">{self._esc(status)}</div>')
        for it in items:
            inner.append('<div class="cov-item">')
            if it.get('caption'):
                inner.append(f'<div class="cap">{self._esc(it["caption"])}</div>')
            # Placeholder token swapped per render mode (interactive div vs static img).
            slot = len(self._cov_slots)
            self._cov_slots.append(it)
            inner.append(f'@@COVSLOT{slot}@@')
            inner.append('</div>')
        self._section(title, ''.join(inner))

    # ---- rendering --------------------------------------------------------- #
    def _render(self, mode):
        css = self._env.from_string(BASE_CSS).render(
            brand_rgb=self.brand, date_long=self.date_long,
            sample_css=_css_string(self.sample))
        body = '\n'.join(self.sections)
        # Swap coverage slots for the chosen mode.
        for slot, it in enumerate(self._cov_slots or []):
            token = f'@@COVSLOT{slot}@@'
            if mode == 'html' and it.get('div'):
                repl = f'<div class="figure">{it["div"]}</div>'
            elif it.get('img') and os.path.exists(it['img']):
                uri = _data_uri(it['img']) if mode == 'html' else _file_uri(it['img'])
                alt = self._esc(f'Coverage depth and SNP graph. {it.get("caption") or ""}'.strip())
                repl = f'<div class="figure"><img src="{uri}" alt="{alt}"></div>'
            else:
                repl = ''
            body = body.replace(token, repl)

        plotly_js = ''
        if mode == 'html' and self._cov_slots:
            try:
                from plotly.offline import get_plotlyjs
                plotly_js = get_plotlyjs()
            except Exception:
                plotly_js = ''

        logo_uri = None
        if self.logo:
            logo_uri = _data_uri(self.logo) if mode == 'html' else _file_uri(self.logo)

        return self._env.from_string(SHELL).render(
            css=css, body=body, sample=self._esc(self.sample), date_long=self.date_long,
            logo_uri=logo_uri, logo_alt=self._esc(self.logo_alt), mode=mode,
            plotly_js=plotly_js)

    def write_html(self):
        if self._cov_slots is None:
            self._cov_slots = []
        path = os.path.join(self.out_dir, f'{self.sample}_{self.stamp}_report.html')
        with open(path, 'w') as fh:
            fh.write(self._render('html'))
        return os.path.abspath(path)

    def write_pdf(self):
        if self._cov_slots is None:
            self._cov_slots = []
        from weasyprint import HTML
        html_str = self._render('pdf')
        path = os.path.join(self.out_dir, f'{self.sample}_{self.stamp}_report.pdf')
        doc = HTML(string=html_str, base_url=os.path.abspath(self.out_dir))
        try:
            # Tagged PDF (headings, tables, alt text) for screen readers / Section 508.
            doc.write_pdf(path, pdf_tags=True)
        except TypeError:
            doc.write_pdf(path)   # older WeasyPrint without the pdf_tags option
        return os.path.abspath(path)


# --------------------------------------------------------------------------- #
# Standalone self-test — troubleshoot report layout/CSS/PDF without a pipeline
# run:   python report_html.py --selftest [--logo path.png]
# Builds a full demo report (every section: overview, FASTQ, pie-less/pie,
# compact tables, assembly, serotype badge, alignment table, coverage graphs)
# from synthetic data and writes both the interactive HTML and the PDF.
# --------------------------------------------------------------------------- #
def _selftest(out_dir='report_html_selftest', logo=None):
    import numpy as np
    from types import SimpleNamespace
    import coverage_plot
    os.makedirs(out_dir, exist_ok=True)

    def fq(name, size, q30, mq, ln):
        return SimpleNamespace(file_name=name, file_size=size, passQ30=q30,
                               read_quality_average=mq, avg_len=ln)

    r = HtmlReport('DEMO-0001', logo=logo, out_dir=out_dir)
    r.add_fastq_quality(fq('DEMO_R1.fastq.gz', '9.1 MB', '98.7', 34.3, '145'),
                        fq('DEMO_R2.fastq.gz', '10.1 MB', '97.8', 31.0, '145'))
    r.add_overview([('Target Taxon', 'Orbivirus'), ('Total Reads', '128,204'),
                    ('Kraken Database', 'k2_orbi_virus'), ('Run Date', r.date_long)])
    r.add_kv('Kraken Classification', [('Target Taxon', 'Orbivirus'), ('Taxon ID', '10892'),
                                       ('Extracted Reads', '127,980'), ('Extraction Rate', '99.83%')])
    r.add_table('Kraken Detailed Classification',
                ['Taxonomic Name', 'Level', 'Reads', 'Percent'],
                [['Orbivirus', 'G', '127,980', '99.83'], ['Bluetongue virus', 'S', '126,101', '98.36'],
                 ['unclassified', 'U', '224', '0.17']], num_cols={2, 3}, compact=True)
    r.add_assembly(SimpleNamespace(contig_count=14, small_contigs_count=2, mid_size=3,
                                   greater_one_kb_count=9, longest_contig=3944,
                                   total_contig_lengths=19204, n50=2772, mean_coverage=512.4,
                                   coverage_title='Mean Coverage'))
    r.add_table('BLAST nt_viruses - Assembly Identification',
                ['nt base count', 'contigs', 'Description'],
                [['3,944', '1', 'Bluetongue virus segment 2 (VP2), complete cds'],
                 ['2,772', '1', 'Bluetongue virus segment 3 (VP3), complete cds']],
                footer='Results provided by: BLAST nt_viruses database', num_cols={0, 1}, compact=True)
    r.add_serotype('BTV-17', 'VP2 matches BTV-17 in the curated panel with high identity.',
                   [{'Protein': 'VP2', 'Serotype': 'BTV-17', 'Top Hit Accession': 'seg2ref_BTV-17_1',
                     'Percent Identity': '98.9', 'Query Coverage': '99.4', 'Bitscore': '5312', 'E-value': '0.0'}])
    # Coverage graphs (small segments -> exercises SVG path + static PNG).
    items, rows = [], []
    for i, (name, length, mean) in enumerate([('DEMO segment 1 (VP1)', 3944, 480),
                                              ('DEMO segment 2 (VP2)', 2926, 510)]):
        rng = np.random.default_rng(i)
        arr = np.clip(np.abs(rng.normal(mean, mean * 0.3, length)), 0, None)
        st = {'coverage_array': arr, 'length': length, 'mean_coverage': float(arr.mean()),
              'percent_covered': 100.0 * float((arr > 0).mean()), 'header': name}
        snps = [{'pos': int(p), 'ref': 'A', 'alt': 'G', 'qual': 55.0}
                for p in rng.choice(np.arange(1, length), size=6, replace=False)]
        items.append({'caption': f'Coverage analysis with SNPs for: {name}',
                      'div': coverage_plot.interactive_div(name, st, snps, [], include_plotlyjs=(i == 0)),
                      'img': coverage_plot.static_png(name, st, snps, [], out_dir=out_dir)})
        rows.append((name, length, float(arr.mean()), 100.0 * float((arr > 0).mean())))
    r.add_alignment_stats(rows)
    r.add_coverage('Bluetongue Virus Coverage Graphs', items)

    html_path = r.write_html()
    print(f'  HTML -> {html_path}')
    try:
        print(f'  PDF  -> {r.write_pdf()}')
    except Exception as e:
        print(f'  PDF skipped ({e})')
    return html_path


def main():
    import argparse
    ap = argparse.ArgumentParser(description='HTML/PDF report builder (pipeline helper). '
                                 'Run with --selftest to render a demo report for troubleshooting.')
    ap.add_argument('--selftest', action='store_true', help='Render a full demo report (HTML + PDF)')
    ap.add_argument('-o', '--out', default='report_html_selftest', help='Output directory')
    ap.add_argument('--logo', default=None, help='Optional logo image for the header')
    args = ap.parse_args()
    if args.selftest:
        print('report_html self-test:')
        _selftest(args.out, logo=args.logo)
    else:
        ap.print_help()


if __name__ == '__main__':
    main()
