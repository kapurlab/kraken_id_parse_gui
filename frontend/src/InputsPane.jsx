import { useEffect, useMemo, useRef, useState } from "react";
import Elapsed from "./Elapsed";
import { askOnPage } from "./pageAsk.jsx";
import { BusyButton, DiskSpace, formatBytes } from "./PaneBits";

// The Inputs pane: every way reads get into a project's download/ folder.
// It offers what the vSNP GUI's Inputs pane offers, laid out the same way, so
// a project can be filled from either tool and look the same in both:
//   - Browse server files…  (pick FASTQs or a folder already on the server)
//   - Pick a folder on the server…  (another project's download/, or a path)
//   - upload / drag & drop, with progress and Cancel
//   - SRA download into an optional sub-folder, with the per-accession outcome
//     of the last download and the accession crosswalk
//   - "Ready to run": what is in download/, one row per sample
// vSNP's Quarantine is not here: it holds samples taken out of vSNP Step 1,
// which this tool does not have.

const READ_TAG_RE = /(?:_R([12])(?:_\d+)?|_([12]))\.fastq\.gz$/i;
const EMPTY_REPORT = { downloaded: [], already_in_step1: [], failed: [] };
const api = (project, rest) => `./api/projects/${encodeURIComponent(project)}${rest}`;

function parseAccessions(text) {
  return (text || "").split(/[\s,]+/).map((s) => s.trim()).filter(Boolean);
}

// download/ files grouped into samples the way the backend pairs them: R1 + R2
// (or _1 + _2) under one name, anything else on its own.
function groupInputs(files) {
  const groups = new Map();
  for (const f of files) {
    const base = f.name.split("/").pop();
    if (!base.endsWith(".fastq.gz")) continue;
    const m = base.match(READ_TAG_RE);
    const dir = f.name.includes("/") ? f.name.slice(0, f.name.lastIndexOf("/") + 1) : "";
    const sample = m ? base.slice(0, m.index) : base.slice(0, -".fastq.gz".length);
    const key = dir + sample;
    if (!groups.has(key)) groups.set(key, { sample, folder: dir.replace(/\/$/, ""), files: [] });
    groups.get(key).files.push(f);
  }
  return [...groups.values()].map((g) => ({
    ...g,
    isPair: g.files.length === 2,
    totalSize: g.files.reduce((n, f) => n + (f.size || 0), 0),
    dangling: g.files.some((f) => f.dangling),
  }));
}

export default function InputsPane({
  project, projects, onSelectProject, onInputsChanged, onSraJob, sraBlocked,
  refreshKey, paneClass, expandButton,
}) {
  const [inputs, setInputs] = useState({ files: [], total_bytes: 0, count: 0 });
  const [inputsSince, setInputsSince] = useState(0);   // > 0 while a listing is in flight
  const [localPath, setLocalPath] = useState("");
  const [pickerMode, setPickerMode] = useState("dropdown");   // dropdown | custom
  const [linkMsg, setLinkMsg] = useState(null);       // {text, since?, error?}
  const [upload, setUpload] = useState(null);         // {text, pct?, since?, error?, active}
  const uploadXhrRef = useRef(null);
  const uploadInputRef = useRef(null);
  const [dragOver, setDragOver] = useState(false);
  const [sraText, setSraText] = useState("");
  const [sraFolder, setSraFolder] = useState("");
  const [sraStatus, setSraStatus] = useState(null);   // {text, since?, error?, active}
  const [sraReport, setSraReport] = useState(EMPTY_REPORT);
  const [browser, setBrowser] = useState({ open: false });
  const projectRef = useRef(project);
  projectRef.current = project;

  const proj = projects.find((p) => p.name === project);
  const krakenRuns = useMemo(() => new Set(proj?.kraken_runs || []), [proj]);
  const groups = useMemo(() => groupInputs(inputs.files), [inputs.files]);
  const readyBytes = groups.reduce((n, g) => n + g.totalSize, 0);
  const hasCrosswalk = inputs.files.some((f) => f.name === "sra_crosswalk.tsv");

  async function loadInputs(p = project) {
    if (!p) { setInputs({ files: [], total_bytes: 0, count: 0 }); return; }
    setInputsSince(Date.now());
    try {
      const r = await fetch(api(p, "/inputs"));
      const d = r.ok ? await r.json() : { files: [], total_bytes: 0, count: 0 };
      if (projectRef.current === p) setInputs(d);
    } catch {
      if (projectRef.current === p) setInputs({ files: [], total_bytes: 0, count: 0 });
    } finally {
      if (projectRef.current === p) setInputsSince(0);
    }
  }

  async function loadSraReport(p = project) {
    if (!p) { setSraReport(EMPTY_REPORT); return; }
    try {
      const r = await fetch(api(p, "/sra-download-report"));
      const d = r.ok ? await r.json() : EMPTY_REPORT;
      if (projectRef.current === p) setSraReport({ ...EMPTY_REPORT, ...d });
    } catch {
      if (projectRef.current === p) setSraReport(EMPTY_REPORT);
    }
  }

  useEffect(() => {
    setLinkMsg(null);
    setLocalPath("");
    setPickerMode("dropdown");
    if (!upload?.active) setUpload(null);
    if (!sraStatus?.active) setSraStatus(null);
    loadInputs(project);
    loadSraReport(project);
  }, [project]);

  useEffect(() => { if (refreshKey) loadInputs(project); }, [refreshKey]);

  async function afterChange() {
    await loadInputs(project);
    onInputsChanged?.(project);
  }

  // --- Link reads already on the server -------------------------------------
  async function linkPaths(paths) {
    if (!project || !paths.length) return null;
    setLinkMsg({ text: `Linking ${paths.length} item${paths.length === 1 ? "" : "s"}…`, since: Date.now() });
    let res;
    try {
      res = await fetch(api(project, "/link-local"), {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ paths }),
      });
    } catch (e) {
      setLinkMsg({ text: `Link failed: ${e.message || "network error"}`, error: true });
      return null;
    }
    const data = await res.json().catch(() => ({}));
    if (!res.ok) {
      setLinkMsg({ text: data.detail || `Link failed (HTTP ${res.status})`, error: true });
      return null;
    }
    const n = data.linked || 0;
    setLinkMsg({
      text: [
        `${n} file${n === 1 ? "" : "s"} linked into download/`,
        data.already_present?.length ? `${data.already_present.length} already there` : null,
        data.skipped_not_fastq?.length ? `${data.skipped_not_fastq.length} skipped (not .fastq.gz)` : null,
        data.missing?.length ? `${data.missing.length} not found` : null,
      ].filter(Boolean).join(" · "),
    });
    await afterChange();
    return data;
  }

  // --- Server file browser (multi-select FASTQs, or a whole folder) ---------
  function browseTo(path) {
    const qs = new URLSearchParams({ path: path || "", include_files: "1", exts: ".fastq.gz" });
    setBrowser((s) => ({ ...s, loading: true, since: Date.now(), error: "" }));
    fetch(`./api/browse-dirs?${qs.toString()}`)
      .then((r) => (r.ok ? r.json() : r.json().then((e) => { throw new Error(e.detail || "Cannot open folder"); })))
      .then((d) => setBrowser((s) => ({
        ...s, loading: false, path: d.path, parent: d.parent,
        entries: d.entries || [], files: d.files || [], disk: d.disk || null, selected: [],
      })))
      .catch((err) => setBrowser((s) => ({ ...s, loading: false, error: err.message, files: [], entries: [] })));
  }
  function openBrowser() {
    setBrowser({ open: true, path: "", parent: null, entries: [], files: [], selected: [], loading: true, error: "" });
    browseTo(localPath || (proj?.path ? `${proj.path}/download` : ""));
  }
  function toggleBrowserFile(path) {
    setBrowser((s) => {
      const sel = new Set(s.selected || []);
      if (sel.has(path)) sel.delete(path); else sel.add(path);
      return { ...s, selected: [...sel] };
    });
  }
  async function linkBrowserSelection() {
    const chosen = (browser.selected || []).length ? browser.selected : (browser.files || []).map((f) => f.path);
    if (!chosen.length) {
      setBrowser((s) => ({ ...s, error: "No .fastq.gz files here to link." }));
      return;
    }
    const result = await linkPaths(chosen);
    if (result) {
      setLocalPath(browser.path);
      setBrowser({ open: false });
    }
  }

  // --- Upload / drag & drop -------------------------------------------------
  async function uploadFiles(fileList) {
    if (!project) return;
    const all = Array.from(fileList || []);
    const files = all.filter((f) => f.name.endsWith(".fastq.gz"));
    const skipped = all.length - files.length;
    if (!files.length) {
      setUpload({ text: all.length ? "Only .fastq.gz files can be added." : "No files selected.", error: true });
      return;
    }
    const total = files.reduce((n, f) => n + f.size, 0);
    const fd = new FormData();
    files.forEach((f) => fd.append("files", f));
    const since = Date.now();
    const label = `${files.length} file${files.length === 1 ? "" : "s"}`;
    const skipNote = skipped ? ` · ${skipped} skipped (not .fastq.gz)` : "";
    setUpload({ text: `Uploading ${label} (0 / ${formatBytes(total)})`, pct: 0, since, active: true });
    await new Promise((resolve) => {
      const xhr = new XMLHttpRequest();
      uploadXhrRef.current = xhr;
      xhr.upload.addEventListener("progress", (e) => {
        if (!e.lengthComputable) return;
        const secs = (Date.now() - since) / 1000;
        const rate = secs > 0 ? `${formatBytes(e.loaded / secs)}/s` : "";
        setUpload({
          text: `Uploading ${label}: ${formatBytes(e.loaded)} of ${formatBytes(total)}${rate ? ` · ${rate}` : ""}`,
          pct: (e.loaded / e.total) * 100, since, active: true,
        });
      });
      xhr.addEventListener("load", () => {
        const secs = Math.round((Date.now() - since) / 1000);
        if (xhr.status >= 200 && xhr.status < 300) {
          let n = files.length;
          try { n = JSON.parse(xhr.responseText).uploaded ?? n; } catch { /* keep count */ }
          setUpload({ text: `Uploaded ${n} file${n === 1 ? "" : "s"} (${formatBytes(total)}) in ${secs}s${skipNote}` });
        } else {
          let detail = `Upload failed (HTTP ${xhr.status})`;
          try { detail = JSON.parse(xhr.responseText).detail || detail; } catch { /* plain body */ }
          setUpload({ text: detail, error: true });
        }
        resolve();
      });
      xhr.addEventListener("error", () => { setUpload({ text: "Upload failed (network error)", error: true }); resolve(); });
      xhr.addEventListener("abort", () => {
        const secs = Math.round((Date.now() - since) / 1000);
        setUpload({ text: `Upload cancelled after ${secs}s. A partly copied file may be listed under Ready to run — remove it with ×.`, error: true });
        resolve();
      });
      xhr.open("POST", api(project, "/upload"));
      xhr.send(fd);
    });
    uploadXhrRef.current = null;
    await afterChange();
  }
  function cancelUpload() {
    uploadXhrRef.current?.abort();
    uploadXhrRef.current = null;
  }

  // --- SRA ------------------------------------------------------------------
  async function sraDownload() {
    const accessions = parseAccessions(sraText);
    if (!project || !accessions.length) return;
    const p = project;
    const since = Date.now();
    const n = accessions.length;
    setSraStatus({ text: `Resolving ${n} accession${n === 1 ? "" : "s"} with NCBI…`, since, active: true });
    let res;
    try {
      res = await fetch(api(p, "/sra/download"), {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ accessions, folder: sraFolder.trim() || null }),
      });
    } catch (e) {
      setSraStatus({ text: `Download failed: ${e.message || "network error"}`, error: true });
      return;
    }
    const data = await res.json().catch(() => ({}));
    if (!res.ok) {
      setSraStatus({ text: `Download failed: ${data.detail || `HTTP ${res.status}`}`, error: true });
      return;
    }
    setSraText("");
    setSraStatus({ text: `Downloading ${n} accession${n === 1 ? "" : "s"} — progress is in the Pipeline Log`, since, active: true });
    const outcome = await onSraJob(data.job_id, p);
    if (projectRef.current !== p) return;
    setSraStatus(outcome === "failed"
      ? { text: "Download finished with errors — see Last download and the Pipeline Log.", error: true }
      : { text: "Download finished." });
    await Promise.all([loadSraReport(p), afterChange()]);
  }

  async function deleteGroup(g) {
    const label = g.files.length > 1 ? `${g.files.map((f) => f.name).join(" + ")}` : g.files[0].name;
    if (!(await askOnPage({ okLabel: "Remove", message: `Remove ${label} from this project's download/ folder? A linked file only loses its link; an uploaded or downloaded file is deleted.` }))) return;
    for (const f of g.files) {
      const r = await fetch(api(project, `/inputs/${f.name.split("/").map(encodeURIComponent).join("/")}`), { method: "DELETE" }).catch(() => null);
      if (!r || !r.ok) {
        const e = r ? await r.json().catch(() => ({})) : {};
        askOnPage({ notice: true, message: `Could not remove ${f.name}: ${e.detail || "network error"}` });
        break;
      }
    }
    await afterChange();
  }

  const otherDownloads = projects.filter((p) => p.name !== project && p.path);
  const sraCount = parseAccessions(sraText).length;
  const reportShown = sraReport.downloaded.length + sraReport.already_in_step1.length + sraReport.failed.length > 0;

  return (
    <section className={paneClass("panel inputs-pane", "inputs")}>
      <div className="pane-head">
        <h2>Inputs</h2>
        {projects.length > 0 && (
          <label className="pane-project" title="The project these reads go into">
            <span className="muted">Project</span>
            <select value={project} onChange={(e) => onSelectProject(e.target.value)}>
              {projects.map((p) => <option key={p.name} value={p.name}>{p.name}</option>)}
            </select>
          </label>
        )}
        {expandButton("inputs", "the Inputs pane")}
      </div>

      {!project ? (
        <div className="empty-msg">Create a project first (top of the Projects pane), then add reads to it here.</div>
      ) : (
        <div className="input-columns">
          {/* ── LEFT: get reads in ───────────────────────────────────── */}
          <div className="input-column">
            <h3>Add local FASTQ</h3>
            <div className="stack-6">
              <button type="button" onClick={openBrowser}
                title="Browse the server and pick FASTQ files, or a whole folder, to link into this project">
                Browse server files…
              </button>
              <select
                value={pickerMode === "custom" ? "__custom__" : ""}
                onChange={async (e) => {
                  const v = e.target.value;
                  if (v === "__custom__") { setPickerMode("custom"); return; }
                  if (!v) return;
                  setLocalPath(v);
                  await linkPaths([v]);
                }}
              >
                <option value="">Pick a folder on the server…</option>
                {otherDownloads.map((p) => (
                  <option key={p.path} value={`${p.path}/download`}>
                    {p.name}/download{p.scope === "shared" ? " (shared)" : ""}
                  </option>
                ))}
                <option value="__custom__">Custom path…</option>
              </select>
              {pickerMode === "custom" && (
                <div className="row tight">
                  <input
                    autoFocus
                    value={localPath}
                    onChange={(e) => setLocalPath(e.target.value)}
                    placeholder="/path/to/folder  or  /path/to/sample_R1.fastq.gz"
                    onKeyDown={(e) => { if (e.key === "Enter" && localPath.trim()) linkPaths([localPath.trim()]); }}
                  />
                  <BusyButton className="ghost action" busyLabel="Linking…"
                    disabled={!localPath.trim()} onClick={() => linkPaths([localPath.trim()])}>Link</BusyButton>
                  <button type="button" className="ghost" title="Back to the list"
                    onClick={() => { setPickerMode("dropdown"); setLocalPath(""); }}>↶</button>
                </div>
              )}
              <div className="form-hint">Linked, not copied: every .fastq.gz in the folder (or the files you pick) is symlinked into download/.</div>
              {linkMsg && (
                <div className={`note inline-note${linkMsg.error ? " error" : ""}`}>
                  <span>{linkMsg.text}{linkMsg.since ? <> <Elapsed since={linkMsg.since} /></> : null}</span>
                  {!linkMsg.since && <button type="button" className="ghost small" onClick={() => setLinkMsg(null)}>Dismiss</button>}
                </div>
              )}
            </div>

            <div className="block">
              <div className="muted sub-label">…or upload / drag &amp; drop</div>
              <div
                className={`dropzone${dragOver ? " drag-over" : ""}`}
                onDragOver={(e) => { e.preventDefault(); setDragOver(true); }}
                onDragLeave={() => setDragOver(false)}
                onDrop={(e) => { e.preventDefault(); setDragOver(false); uploadFiles(e.dataTransfer.files); }}
              >
                <input ref={uploadInputRef} type="file" multiple accept=".gz,application/gzip" style={{ display: "none" }}
                  onChange={(e) => { const f = Array.from(e.target.files); e.target.value = ""; uploadFiles(f); }} />
                <button type="button" disabled={!!upload?.active} onClick={() => uploadInputRef.current?.click()}>Choose Files</button>
                <span className="drop-hint">Or drop FASTQ.GZ files here</span>
              </div>
              {upload && (
                <div className={`note inline-note${upload.error ? " error" : ""}`}>
                  <span style={{ flex: 1 }}>
                    {upload.text}{upload.since && upload.active ? <> <Elapsed since={upload.since} /></> : null}
                    {upload.active && <span className="progress"><span style={{ width: `${upload.pct || 0}%` }} /></span>}
                  </span>
                  {upload.active
                    ? <button type="button" className="ghost small" onClick={cancelUpload} title="Stop the upload">Cancel</button>
                    : <button type="button" className="ghost small" onClick={() => setUpload(null)}>Dismiss</button>}
                </div>
              )}
            </div>

            <h3>SRA Download</h3>
            <div className="stack-6">
              <textarea
                rows={5}
                placeholder={"SRR/ERR/DRR, SRX/SRS or PRJNA accessions\n(one per line)"}
                value={sraText}
                onChange={(e) => setSraText(e.target.value)}
              />
              <input
                placeholder="Optional sub-folder (e.g. 2026-02-01_batch1)"
                value={sraFolder}
                onChange={(e) => setSraFolder(e.target.value)}
              />
              <BusyButton
                onClick={sraDownload}
                busyLabel="Starting…"
                disabled={!sraCount || !!sraStatus?.active || sraBlocked}
                title={sraBlocked ? "Wait for the run in progress to finish — the download shares its log" : ""}
              >
                Download{sraCount ? ` (${sraCount})` : ""}
              </BusyButton>
              {sraStatus && (
                <div className={`note inline-note${sraStatus.error ? " error" : ""}`}>
                  <span>
                    {sraStatus.active && <span className="pulse-dot" />}
                    {sraStatus.text}{sraStatus.active && sraStatus.since ? <> <Elapsed since={sraStatus.since} /></> : null}
                  </span>
                  {!sraStatus.active && <button type="button" className="ghost small" onClick={() => setSraStatus(null)}>Dismiss</button>}
                </div>
              )}
              {reportShown && (
                <div className="note sra-report">
                  <div className="sra-report-head">Last download</div>
                  <div className="muted">
                    {sraReport.downloaded.length} downloaded
                    {sraReport.already_in_step1.length ? ` · ${sraReport.already_in_step1.length} already in this project (skipped)` : ""}
                    {sraReport.failed.length ? ` · ${sraReport.failed.length} failed` : ""}
                  </div>
                  {sraReport.already_in_step1.length > 0 && (
                    <details>
                      <summary>Already in this project — {sraReport.already_in_step1.length} skipped</summary>
                      <pre>{sraReport.already_in_step1.join("\n")}</pre>
                    </details>
                  )}
                  {sraReport.failed.length > 0 && (
                    <details open>
                      <summary className="danger-text">Failed — {sraReport.failed.length}</summary>
                      <pre>{sraReport.failed.join("\n")}</pre>
                    </details>
                  )}
                </div>
              )}
              {hasCrosswalk && (
                <div className="form-hint">
                  <a href={api(project, "/sra-crosswalk")} target="_blank" rel="noopener noreferrer">View crosswalk</a>
                  {" "}— the sample/study accessions you entered → the run accessions downloaded.
                </div>
              )}
            </div>
          </div>

          {/* ── RIGHT: what is in download/ ──────────────────────────── */}
          <div className="input-column">
            <h3 className="ready-head">
              <span>
                Ready to run
                {groups.length > 0 && (
                  <span className="muted ready-count">
                    {groups.length} sample{groups.length === 1 ? "" : "s"} · {formatBytes(readyBytes)}
                  </span>
                )}
              </span>
              <button type="button" className="ghost small" onClick={() => loadInputs()} disabled={!!inputsSince}
                title="Read download/ again">
                {inputsSince ? <Elapsed since={inputsSince} /> : "Refresh"}
              </button>
            </h3>
            {groups.length === 0 ? (
              <div className="empty-msg">
                {inputsSince ? <>Reading download/… <Elapsed since={inputsSince} /></> : "No reads yet. Add some on the left."}
              </div>
            ) : (
              <div className="ready-list">
                {groups.map((g) => {
                  const done = krakenRuns.has(g.sample);
                  const tip = g.files.map((f) => `${f.name} (${formatBytes(f.size)})`).join("\n");
                  return (
                    <div key={g.files.map((f) => f.name).join("|")} className={`ready-row${g.dangling ? " dangling" : ""}`}>
                      <span className="ready-name" title={tip}>
                        {g.sample}
                        {g.folder && <span className="muted ready-folder"> · {g.folder}/</span>}
                      </span>
                      <span className={`read-badge ${g.isPair ? "badge-pe" : "badge-se"}`}
                        title={g.isPair ? "Paired-end: R1 + R2" : "Single file"}>{g.isPair ? "PE" : "SE"}</span>
                      {g.dangling
                        ? <span className="state-chip warn" title="The file this link points to is gone">link broken</span>
                        : done
                          ? <span className="state-chip ok" title="This sample has Kraken results in this project">✓ results</span>
                          : <span className="state-chip" title="Not run through Kraken yet">not run</span>}
                      <span className="muted ready-size">{formatBytes(g.totalSize)}</span>
                      <button type="button" className="chip-remove" onClick={() => deleteGroup(g)}
                        title={g.isPair ? `Remove both ${g.files.map((f) => f.name).join(" + ")}` : `Remove ${g.files[0].name}`}>×</button>
                    </div>
                  );
                })}
              </div>
            )}
            <div className="form-hint">
              Everything here is in <code>{project}/download/</code> — the same folder the vSNP GUI reads,
              so reads added in either tool show up in both. Select samples to run in the Projects pane.
            </div>
          </div>
        </div>
      )}

      {browser.open && (
        <div className="modal-backdrop" onClick={() => setBrowser({ open: false })}>
          <div className="modal" onClick={(e) => e.stopPropagation()}>
            <div className="modal-head">
              <span>Select FASTQ files on the server</span>
              <DiskSpace disk={browser.disk} />
            </div>
            <div className="modal-bar">
              <button type="button" className="ghost" disabled={!browser.parent || browser.loading} onClick={() => browseTo(browser.parent)}>↑ Up</button>
              <input
                value={browser.path || ""}
                onChange={(e) => setBrowser((s) => ({ ...s, path: e.target.value }))}
                onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); browseTo(browser.path); } }}
              />
              <button type="button" className="ghost" onClick={() => browseTo(browser.path)}>Go</button>
            </div>
            <div className="modal-body">
              {browser.loading ? (
                <div className="muted pad">Loading… <Elapsed since={browser.since} /></div>
              ) : browser.error ? (
                <div className="pad danger-text">{browser.error}</div>
              ) : !(browser.entries || []).length && !(browser.files || []).length ? (
                <div className="muted pad">No sub-folders or .fastq.gz files here.</div>
              ) : (
                <>
                  {(browser.entries || []).map((e) => (
                    <div key={e.path} className="browse-row" onClick={() => browseTo(e.path)}>
                      <span>📁</span><span>{e.name}</span>
                    </div>
                  ))}
                  {(browser.files || []).map((f) => (
                    <label key={f.path} className="browse-row">
                      <input type="checkbox" checked={(browser.selected || []).includes(f.path)} onChange={() => toggleBrowserFile(f.path)} />
                      <span>🧬</span>
                      <span style={{ flex: 1 }}>{f.name}</span>
                      <span className="muted small">{formatBytes(f.size)}</span>
                    </label>
                  ))}
                </>
              )}
            </div>
            {(browser.files || []).length > 0 && (
              <div className="modal-bar small">
                <button type="button" className="ghost small"
                  onClick={() => setBrowser((s) => ({ ...s, selected: (s.files || []).map((f) => f.path) }))}>
                  Select all {browser.files.length}
                </button>
                <button type="button" className="ghost small" disabled={!(browser.selected || []).length}
                  onClick={() => setBrowser((s) => ({ ...s, selected: [] }))}>Clear</button>
                <span className="muted">{(browser.selected || []).length} selected</span>
              </div>
            )}
            <div className="modal-foot">
              <button type="button" className="ghost" onClick={() => setBrowser({ open: false })}>Cancel</button>
              <BusyButton busyLabel="Linking…" onClick={linkBrowserSelection}
                disabled={browser.loading || !browser.path || !(browser.files || []).length}
                onError={(err) => setBrowser((s) => ({ ...s, error: err?.message || String(err) }))}
                title="Symlink the selection into this project's download/ folder">
                {(browser.selected || []).length
                  ? `Link ${browser.selected.length} selected`
                  : `Link all ${(browser.files || []).length || ""} FASTQs in this folder`}
              </BusyButton>
            </div>
          </div>
        </div>
      )}
    </section>
  );
}
