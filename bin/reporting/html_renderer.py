from pathlib import Path
from typing import Any, Dict
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


def _asset_url(path: Any, mode: str = "filesystem") -> str:
    path_text = str(path or "")
    if mode == "api":
        return f"./file?path={quote(path_text)}&inline=1"
    return path_text


def render_html_text(manifest: Dict[str, Any], asset_mode: str = "filesystem") -> str:
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    template_dir = Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["fmt"] = _format_value
    env.globals["asset_url"] = lambda path: _asset_url(path, asset_mode)
    template = env.get_template("report.html.j2")
    return template.render(report=manifest)


def render_html_report(manifest: Dict[str, Any], output_dir: Path) -> Path:
    html = render_html_text(manifest, asset_mode="filesystem")
    out = output_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    return out
