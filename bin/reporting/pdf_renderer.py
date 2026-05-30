import asyncio
import base64
import json
import os
import shutil
import subprocess
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Optional, Tuple


def _truthy(value: str) -> bool:
    return value.strip().lower() not in {"0", "false", "no", "off"}


def _find_browser_executable() -> Optional[str]:
    env_path = os.environ.get("KRAKEN_REPORT_BROWSER", "").strip()
    if env_path and Path(env_path).exists():
        return env_path

    names = ["chromium", "chromium-browser", "google-chrome", "google-chrome-stable", "microsoft-edge"]
    for name in names:
        found = shutil.which(name)
        if found:
            return found

    mac_paths = [
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    ]
    for path in mac_paths:
        if Path(path).exists():
            return path
    return None


def _render_with_playwright(html_path: Path, pdf_path: Path) -> Optional[str]:
    try:
        from playwright.sync_api import Error as PlaywrightError
        from playwright.sync_api import sync_playwright
    except ImportError:
        return "Playwright Python package is not installed"

    try:
        with sync_playwright() as p:
            browser = None
            try:
                browser = p.chromium.launch()
            except PlaywrightError:
                browser = p.chromium.launch(channel="chrome")
            page = browser.new_page(viewport={"width": 1280, "height": 1600})
            page.goto(html_path.resolve().as_uri(), wait_until="networkidle")
            page.pdf(
                path=str(pdf_path),
                format="Letter",
                print_background=True,
                margin={"top": "0.5in", "right": "0.45in", "bottom": "0.5in", "left": "0.45in"},
            )
            browser.close()
    except Exception as exc:
        return f"Playwright PDF export failed: {exc}"
    return None


async def _cdp_send(websocket, method: str, params: Optional[dict] = None, msg_id: int = 1) -> dict:
    await websocket.send(json.dumps({"id": msg_id, "method": method, "params": params or {}}))
    while True:
        message = json.loads(await websocket.recv())
        if message.get("id") == msg_id:
            return message


def _read_json_url(url: str, timeout: float = 1.0) -> Optional[object]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception:
        return None


def _render_with_chrome_cdp(html_path: Path, pdf_path: Path) -> Optional[str]:
    try:
        import websockets
    except ImportError:
        return "websockets Python package is not installed"

    browser = _find_browser_executable()
    if not browser:
        return "No Chromium/Chrome/Edge executable found for PDF export"

    with tempfile.TemporaryDirectory(prefix="kraken-report-chrome-") as user_data_dir:
        cmd = [
            browser,
            "--headless=new",
            "--disable-gpu",
            "--no-sandbox",
            "--remote-debugging-port=0",
            f"--user-data-dir={user_data_dir}",
            "about:blank",
        ]
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except OSError as exc:
            return f"Chrome DevTools PDF export failed to launch: {exc}"

        port = None
        stderr_lines = []
        try:
            deadline = time.time() + 10
            while time.time() < deadline and port is None:
                line = proc.stderr.readline() if proc.stderr else ""
                if line:
                    stderr_lines.append(line.strip())
                    marker = "DevTools listening on ws://"
                    if marker in line:
                        host_port = line.split(marker, 1)[1].split("/", 1)[0]
                        port = int(host_port.rsplit(":", 1)[1])
                elif proc.poll() is not None:
                    break

            if port is None:
                msg = "\n".join(stderr_lines[-4:])
                return f"Chrome DevTools PDF export did not start: {msg}"

            targets = None
            deadline = time.time() + 10
            while time.time() < deadline:
                targets = _read_json_url(f"http://127.0.0.1:{port}/json/list")
                if isinstance(targets, list) and targets:
                    break
                time.sleep(0.1)
            if not isinstance(targets, list) or not targets:
                return "Chrome DevTools PDF export could not find a page target"

            page_target = next((target for target in targets if target.get("type") == "page"), targets[0])
            ws_url = page_target.get("webSocketDebuggerUrl")
            if not ws_url:
                return "Chrome DevTools PDF export page target had no websocket URL"

            async def print_pdf() -> Optional[str]:
                async with websockets.connect(ws_url, max_size=32 * 1024 * 1024) as websocket:
                    msg_id = 1
                    await _cdp_send(websocket, "Page.enable", msg_id=msg_id)
                    msg_id += 1
                    await _cdp_send(
                        websocket,
                        "Page.navigate",
                        {"url": html_path.resolve().as_uri()},
                        msg_id=msg_id,
                    )
                    msg_id += 1
                    while True:
                        event = json.loads(await websocket.recv())
                        if event.get("method") == "Page.loadEventFired":
                            break
                    result = await _cdp_send(
                        websocket,
                        "Page.printToPDF",
                        {
                            "displayHeaderFooter": False,
                            "printBackground": True,
                            "paperWidth": 8.5,
                            "paperHeight": 11,
                            "marginTop": 0.5,
                            "marginBottom": 0.5,
                            "marginLeft": 0.5,
                            "marginRight": 0.5,
                            "preferCSSPageSize": True,
                        },
                        msg_id=msg_id,
                    )
                    if "error" in result:
                        return str(result["error"])
                    data = result.get("result", {}).get("data")
                    if not data:
                        return "Chrome DevTools PDF export returned no PDF data"
                    pdf_path.write_bytes(base64.b64decode(data))
                    return None

            warning = asyncio.run(print_pdf())
            if warning:
                return warning
            if not pdf_path.exists() or pdf_path.stat().st_size == 0:
                return "Chrome DevTools PDF export completed but did not produce a non-empty PDF"
        finally:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=5)
    return None


def _render_with_browser_cli(html_path: Path, pdf_path: Path) -> Optional[str]:
    browser = _find_browser_executable()
    if not browser:
        return "No Chromium/Chrome/Edge executable found for PDF export"

    cmd = [
        browser,
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--print-to-pdf-no-header",
        f"--print-to-pdf={pdf_path}",
        html_path.resolve().as_uri(),
    ]
    try:
        result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"Browser PDF export failed to launch: {exc}"
    if result.returncode != 0:
        msg = (result.stderr or result.stdout or "").strip()
        return f"Browser PDF export failed with exit {result.returncode}: {msg}"
    if not pdf_path.exists() or pdf_path.stat().st_size == 0:
        return "Browser PDF export completed but did not produce a non-empty PDF"
    return None


def render_pdf_report(html_path: Path, output_dir: Path) -> Tuple[Optional[Path], Optional[str]]:
    """Render report.html to report.pdf when enabled.

    Returns (pdf_path, warning). PDF export is intentionally non-fatal because
    report.html is the canonical OOD report artifact.
    """
    if not _truthy(os.environ.get("REPORT_EXPORT_PDF", "1")):
        return None, "PDF export skipped because REPORT_EXPORT_PDF is disabled"

    pdf_path = output_dir / "report.pdf"
    if pdf_path.exists():
        pdf_path.unlink()

    warning = _render_with_playwright(html_path, pdf_path)
    if warning and "not installed" in warning:
        warning = _render_with_chrome_cdp(html_path, pdf_path)
    if warning and ("not installed" in warning or "DevTools" in warning):
        warning = _render_with_browser_cli(html_path, pdf_path)
    if warning:
        return None, warning
    return pdf_path, None
