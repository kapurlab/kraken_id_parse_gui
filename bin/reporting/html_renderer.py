import base64
import mimetypes
import re
from pathlib import Path
from typing import Any, Dict, Optional
from urllib.parse import quote


def _format_value(value: Any, decimals: int = None, suffix: str = "") -> str:
    if value is None or value == "":
        return "not recorded"
    if isinstance(value, str):
        raw = value.strip()
        if not raw:
            return "not recorded"
        if raw.endswith("%"):
            suffix = suffix or "%"
            raw = raw[:-1].strip()
        try:
            number = float(raw.replace(",", ""))
        except ValueError:
            return value
    else:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)

    if decimals is not None:
        return f"{number:,.{decimals}f}{suffix}"
    if number.is_integer():
        return f"{int(number):,}{suffix}"
    return f"{number:,.1f}{suffix}"


def _asset_url(path: Any, mode: str = "filesystem", base_dir: Optional[Path] = None) -> str:
    path_text = str(path or "")
    if not path_text:
        return path_text
    if mode == "api":
        return f"./file?path={quote(path_text)}&inline=1"
    if mode == "embed":
        # Inline the asset as a base64 data URI so the report is fully
        # self-contained — it renders identically no matter which GUI (or
        # neither) serves it, how it's proxied, or where it's copied to. This
        # is what keeps image links from breaking across the vSNP / Kraken
        # tools, which serve the same file from different URL roots.
        p = Path(path_text)
        if not p.is_absolute() and base_dir is not None:
            p = Path(base_dir) / path_text
        try:
            data = p.read_bytes()
        except OSError:
            return path_text
        mime = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        return f"data:{mime};base64,{base64.b64encode(data).decode('ascii')}"
    return path_text


def _file_url(path: Any, mode: str = "filesystem") -> str:
    """A link to another file of the run (the Krona chart). Never a data URI:
    browsers refuse to navigate to one, so in the self-contained report the
    link stays relative to the run folder the report sits in."""
    path_text = str(path or "")
    if not path_text:
        return path_text
    if mode == "api":
        return f"./file?path={quote(path_text)}&inline=1"
    return quote(path_text)


def _css_string(value: Any) -> str:
    """Text for a CSS string literal (the PDF's running page header).

    Escaped for CSS here; the template marks the result |safe, because inside
    <style> the browser does not decode HTML entities, so autoescaping it
    would print "&#34;". (No markupsafe import: it is not a declared
    dependency of this env, only jinja2's.)"""
    text = re.sub(r"[\x00-\x1f]", " ", str(value or ""))
    text = text.replace("\\", "\\\\").replace('"', '\\"').replace("<", "").replace(">", "")
    return text


def _run_dir(manifest: Dict[str, Any], base_dir: Optional[Path]) -> Optional[Path]:
    """The run folder, for manifests too old to carry the Kraken overview:
    the caller's base_dir, else the folder the run recorded as its own."""
    for cand in (base_dir, (manifest.get("environment") or {}).get("cwd")):
        if cand:
            try:
                if Path(cand).is_dir():
                    return Path(cand)
            except OSError:
                continue
    return None


def render_html_text(
    manifest: Dict[str, Any],
    asset_mode: str = "filesystem",
    base_dir: Optional[Path] = None,
) -> str:
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    try:
        from .summary import build_summary, find_kraken_report, kraken_overview, rank_label
    except ImportError:  # imported as a top-level module (bin/reporting on sys.path)
        from summary import build_summary, find_kraken_report, kraken_overview, rank_label

    template_dir = Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        # The template is report.html.j2: matched on "html" alone, autoescape
        # was never on for it, so a BLAST title or sample name holding "<" or
        # "&" went into the page as markup. Escape everything; the one
        # pre-rendered HTML fragment (the Plotly chart) is marked |safe.
        autoescape=select_autoescape(["html", "xml", "j2"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["fmt"] = _format_value
    env.filters["rank"] = rank_label
    env.filters["css_str"] = _css_string
    env.globals["asset_url"] = lambda path: _asset_url(path, asset_mode, base_dir)
    env.globals["file_url"] = lambda path: _file_url(path, asset_mode)
    template = env.get_template("report.html.j2")
    try:
        from orbivirus_beta import PAGE as beta_page, notice_for_manifest
        beta_notice = notice_for_manifest(manifest)
    except ImportError:
        beta_notice, beta_page = manifest.get("beta_notice"), None
    kraken = (manifest.get("legacy_report") or {}).get("kraken_overview")
    if not kraken:
        kraken = kraken_overview(find_kraken_report(_run_dir(manifest, base_dir)))
    summary = build_summary(manifest, kraken=kraken, beta_notice=beta_notice)
    return template.render(report=manifest, beta_notice=beta_notice,
                           beta_page=(beta_page or (beta_notice or {}).get("title")) if beta_notice else None,
                           summary=summary)


def render_html_report(manifest: Dict[str, Any], output_dir: Path) -> Path:
    # Write a SELF-CONTAINED report.html (images embedded as data URIs). A
    # self-contained file is the robust choice here: the same report is opened
    # by both the Kraken GUI and the vSNP GUI, which serve it from different
    # URL roots, plus it's downloaded and emailed — relative/filesystem asset
    # paths break in all of those, embedded assets break in none.
    html = render_html_text(manifest, asset_mode="embed", base_dir=output_dir)
    out = output_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    return out
