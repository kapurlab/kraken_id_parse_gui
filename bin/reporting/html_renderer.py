import base64
import mimetypes
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


def render_html_text(
    manifest: Dict[str, Any],
    asset_mode: str = "filesystem",
    base_dir: Optional[Path] = None,
) -> str:
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    template_dir = Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["fmt"] = _format_value
    env.globals["asset_url"] = lambda path: _asset_url(path, asset_mode, base_dir)
    template = env.get_template("report.html.j2")
    return template.render(report=manifest)


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
