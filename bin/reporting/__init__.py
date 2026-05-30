"""Structured report generation for Kraken ID Parse."""

from .html_renderer import render_html_report
from .manifest import build_run_manifest, write_manifest
from .pdf_renderer import render_pdf_report

__all__ = ["build_run_manifest", "render_html_report", "render_pdf_report", "write_manifest"]
