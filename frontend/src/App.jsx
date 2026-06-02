import { useState, useEffect, useRef } from "react";
import "./App.css";

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const APP_VERSION = "0.2.0";

const TAXON_PRESETS = [
  "Mycobacterium tuberculosis complex",
  "Mycobacterium bovis",
  "Orbivirus",
  "Apicomplexa",
  "Isavirus salaris",
];

function fileIcon(name) {
  if (name.endsWith(".json")) return "📁";
  if (name.endsWith(".xlsx")) return "📊";
  if (name.endsWith(".pdf")) return "📄";
  if (name.endsWith(".png")) return "🖼";
  if (name.endsWith(".fasta") || name.endsWith(".fa")) return "🧬";
  if (name.endsWith(".vcf")) return "🔬";
  if (name.endsWith(".txt")) return "📝";
  if (name.endsWith(".html")) return "🌐";
  return "📁";
}

function fmtSize(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / 1024 / 1024).toFixed(1)} MB`;
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------
export default function App() {
  const [projects, setProjects] = useState([]);
  const [projectsLoading, setProjectsLoading] = useState(true);
  const [expanded, setExpanded] = useState({});          // project name → bool
  const [samples, setSamples] = useState({});            // project name → [sample]
  const [checkedKeys, setCheckedKeys] = useState({});    // key → {project, ...sample}  (batch selection)
  const [openResults, setOpenResults] = useState({});    // key → bool (inline results expanded)
  const [sampleResults, setSampleResults] = useState({}); // key → {loading, status, present, files}
  const [vsnpResults, setVsnpResults] = useState({});    // key → {loading, step1_present, files, step2} (cross-tool)
  const [activeRun, setActiveRun] = useState(null);      // {project, sample} currently running
  const [queueInfo, setQueueInfo] = useState({ total: 0, done: 0 }); // batch progress
  const [taxon, setTaxon] = useState("");
  const [krakenOnly, setKrakenOnly] = useState(false);   // Kraken2 + Krona only, no read parsing
  const [krakenDb, setKrakenDb] = useState("");
  const [blastDb, setBlastDb] = useState("nt");
  const [running, setRunning] = useState(false);
  const [jobId, setJobId] = useState(null);
  const [jobStatus, setJobStatus] = useState("idle"); // idle | running | succeeded | failed
  const [logLines, setLogLines] = useState([]);
  const [settingsDraft, setSettingsDraft] = useState({});
  const [currentStep, setCurrentStep] = useState("");

  // Section visibility (collapsible flow, adapted from latex gui)
  const [showSettings, setShowSettings] = useState(false);
  const [showProjects, setShowProjects] = useState(true);
  const [showRun, setShowRun] = useState(true);
  const [showLogs, setShowLogs] = useState(true);

  const logRef = useRef(null);
  const eventSourceRef = useRef(null);

  // Load config & projects on mount; reconnect to any pipeline still running
  useEffect(() => {
    fetch("./api/config")
      .then((r) => r.json())
      .then((cfg) => {
        setKrakenDb(cfg.kraken_db || "");
        setBlastDb(cfg.blast_db || "nt");
        setSettingsDraft(cfg);
      })
      .catch(() => {});
    loadProjects();
    fetch("./api/jobs")
      .then((r) => r.json())
      .then((jobs) => {
        const live = jobs.find((j) => j.status === "running");
        if (live) {
          setJobId(live.id);
          setJobStatus("running");
          setRunning(true);
          // Reconstruct the running sample from the job name ("project/sample — taxon")
          let samp = null;
          const m = (live.name || "").match(/^(.*?)\/(.*?) — /);
          if (m) {
            samp = { project: m[1], sample: m[2] };
            setActiveRun(samp);
          }
          streamLogUntilDone(live.id, samp, () => {});
        }
      })
      .catch(() => {});
  }, []);

  function loadProjects() {
    setProjectsLoading(true);
    fetch("./api/projects")
      .then((r) => r.json())
      .then((data) => {
        setProjects(data);
        setProjectsLoading(false);
      })
      .catch(() => setProjectsLoading(false));
  }

  // Auto-scroll to the bottom on new lines, but only if the user is already
  // near the bottom — otherwise leave their scroll position alone so they can
  // read back through the log while the pipeline is still running.
  useEffect(() => {
    const el = logRef.current;
    if (!el) return;
    const nearBottom =
      el.scrollHeight - el.scrollTop - el.clientHeight < 80;
    if (nearBottom) {
      el.scrollTop = el.scrollHeight;
    }
  }, [logLines]);

  function toggleProject(name) {
    const isExpanded = expanded[name];
    setExpanded((e) => ({ ...e, [name]: !isExpanded }));
    if (!isExpanded && !samples[name]) {
      fetch(`./api/projects/${encodeURIComponent(name)}/samples`)
        .then((r) => r.json())
        .then((data) => setSamples((s) => ({ ...s, [name]: data })))
        .catch(() => setSamples((s) => ({ ...s, [name]: [] })));
    }
  }

  // --- Sample selection / results (per-sample, decoupled from a single job) ---
  const sampleKey = (project, s) => `${project}::${s.sample}`;
  const isActive = (project, s) =>
    activeRun && activeRun.project === project && activeRun.sample === s.sample;

  function toggleChecked(project, s) {
    const key = sampleKey(project, s);
    setCheckedKeys((m) => {
      const next = { ...m };
      if (next[key]) delete next[key];
      else next[key] = { project, ...s };
      return next;
    });
  }

  function loadSampleResults(project, s) {
    const key = sampleKey(project, s);
    setSampleResults((m) => ({ ...m, [key]: { ...(m[key] || {}), loading: true } }));
    fetch(`./api/projects/${encodeURIComponent(project)}/samples/${encodeURIComponent(s.sample)}/kraken-results`)
      .then((r) => r.json())
      .then((data) => setSampleResults((m) => ({ ...m, [key]: { loading: false, ...data } })))
      .catch(() => setSampleResults((m) => ({ ...m, [key]: { loading: false, present: false, status: "none", files: [] } })));
  }

  // Cross-tool: vSNP results for the same sample (step1 files + latest step2).
  function loadVsnpResults(project, s) {
    const key = sampleKey(project, s);
    setVsnpResults((m) => ({ ...m, [key]: { ...(m[key] || {}), loading: true } }));
    fetch(`./api/projects/${encodeURIComponent(project)}/vsnp/samples/${encodeURIComponent(s.sample)}/files`)
      .then((r) => r.json())
      .then((data) => setVsnpResults((m) => ({ ...m, [key]: { loading: false, ...data } })))
      .catch(() => setVsnpResults((m) => ({ ...m, [key]: { loading: false, step1_present: false, files: [], step2: { present: false } } })));
  }

  function toggleResults(project, s) {
    const key = sampleKey(project, s);
    const willOpen = !openResults[key];
    setOpenResults((m) => ({ ...m, [key]: willOpen }));
    if (willOpen && !sampleResults[key]) loadSampleResults(project, s);
    if (willOpen && !vsnpResults[key]) loadVsnpResults(project, s);
  }

  // Run one or more samples back-to-back (sequential — avoids overloading the
  // box with concurrent heavy pipelines, and keeps a single coherent live log).
  async function runSamples(list) {
    if (running || !list.length) return;
    if (!krakenOnly && !taxon.trim()) return;
    setShowLogs(true);
    setQueueInfo({ total: list.length, done: 0 });
    for (let i = 0; i < list.length; i++) {
      await runOne(list[i]);
      setQueueInfo({ total: list.length, done: i + 1 });
    }
    setActiveRun(null);
  }

  function runSelected() {
    runSamples(Object.values(checkedKeys));
  }

  function runOne(samp) {
    return new Promise((resolve) => {
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
        eventSourceRef.current = null;
      }
      setRunning(true);
      setActiveRun({ project: samp.project, sample: samp.sample });
      setJobStatus("running");
      setLogLines([]);
      setCurrentStep("");
      // Mark this sample as running in its inline panel immediately.
      const key = sampleKey(samp.project, samp);
      setSampleResults((m) => ({ ...m, [key]: { ...(m[key] || {}), status: "running" } }));
      setOpenResults((m) => ({ ...m, [key]: true }));

      fetch("./api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          project: samp.project,
          r1: samp.r1,
          r2: samp.r2 || null,
          taxon: taxon.trim(),
          kraken_db: krakenDb.trim() || null,
          blast_db: blastDb.trim() || null,
          kraken_only: krakenOnly,
        }),
      })
        .then((r) => (r.ok ? r.json() : r.json().then((e) => { throw new Error(e.detail || "Run failed"); })))
        .then(({ job_id }) => {
          setJobId(job_id);
          streamLogUntilDone(job_id, samp, resolve);
        })
        .catch((err) => {
          setLogLines((prev) => [...prev, `ERROR: ${err.message}`]);
          setRunning(false);
          setJobStatus("failed");
          resolve();
        });
    });
  }

  function streamLogUntilDone(id, samp, done) {
    const es = new EventSource(`./api/jobs/${id}/log`);
    eventSourceRef.current = es;
    es.onmessage = (evt) => {
      const data = evt.data;
      if (data === "[DONE]") {
        es.close();
        setRunning(false);
        fetch(`./api/jobs/${id}`)
          .then((r) => r.json())
          .then((job) => {
            setJobStatus(job.status);
            setCurrentStep("");
            if (samp) loadSampleResults(samp.project, samp); // refresh this sample's inline results
            loadProjects();                                  // refresh kraken_runs badges
          })
          .catch(() => {})
          .finally(() => done());
      } else {
        setLogLines((prev) => [...prev, data]);
        if (/Step \d+:/i.test(data) ||
            /Starting bioinformatics/i.test(data) ||
            /Generating analysis reports/i.test(data) ||
            /Pipeline completed/i.test(data) ||
            /Downloading.*reference/i.test(data)) {
          setCurrentStep(data.trim().replace(/^#+\s*/, ""));
        }
      }
    };
    es.onerror = () => {
      es.close();
      setRunning(false);
      setJobStatus("failed");
      done();
    };
  }

  function saveSettings() {
    fetch("./api/config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        kraken_db: settingsDraft.kraken_db,
        blast_db: settingsDraft.blast_db,
        projects_root: settingsDraft.projects_root,
      }),
    })
      .then((r) => r.json())
      .then(() => {
        setKrakenDb(settingsDraft.kraken_db || "");
        setBlastDb(settingsDraft.blast_db || "nt");
        loadProjects();
      })
      .catch(() => {});
  }

  const logLineClass = (line) => {
    if (line.startsWith("$ ")) return "log-line cmd";
    if (line.startsWith("ERROR") || line.startsWith("error")) return "log-line error";
    if (line === "[DONE]") return "log-line done";
    return "log-line";
  };

  const statusText = {
    idle: "idle",
    running: "running",
    succeeded: "succeeded",
    failed: "failed",
  }[jobStatus];

  return (
    <div className="app">
      {/* ── Header ─────────────────────────────────────────────── */}
      <header className="app-header">
        <div className="app-brand">
          <img className="app-logo" src="./kraken_icon.png" alt="Kraken icon" />
          <div>
            <h1>
              Kraken ID Parse <span className="version-tag">v{APP_VERSION}</span>
            </h1>
            <p>Classify and isolate reads for species-level identification and contamination screening</p>
          </div>
        </div>
        <div className="status-pill">
          <span className="dot" data-state={jobStatus} />
          <span>{statusText}</span>
        </div>
      </header>

      <main className="layout">
        {/* ── Alert banner ─────────────────────────────────────── */}
        {!krakenDb && (
          <div className="alert-banner">
            <strong>Setup needed:</strong> No Kraken2 database configured. Open{" "}
            <button
              className="ghost action"
              style={{ padding: "2px 8px", fontSize: 12 }}
              onClick={() => setShowSettings(true)}
            >
              Settings
            </button>{" "}
            to set the database path before running.
          </div>
        )}

        {/* ── Status strip ─────────────────────────────────────── */}
        <section className="status-strip">
          <div className="status-item">
            <span className="status-label">Selected</span>
            <span className="status-value">
              {Object.keys(checkedKeys).length
                ? `${Object.keys(checkedKeys).length} sample${Object.keys(checkedKeys).length > 1 ? "s" : ""}`
                : "—"}
            </span>
          </div>
          <div className="status-item">
            <span className="status-label">Running</span>
            <span className="status-value">
              {activeRun ? activeRun.sample : "—"}
              {queueInfo.total > 1 ? ` (${queueInfo.done}/${queueInfo.total})` : ""}
            </span>
          </div>
          <div className="status-item">
            <span className="status-label">Target Taxon</span>
            <span className="status-value">{taxon.trim() || "—"}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Job</span>
            <span className="status-value cap">
              {jobStatus === "running" ? <><span className="pulse-dot" />running</> : statusText}
            </span>
          </div>
        </section>

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Settings (collapsed by default)                */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Settings</h2>
          <button className="ghost" onClick={() => {
            if (!showSettings) {
              fetch("./api/config").then((r) => r.json()).then(setSettingsDraft).catch(() => {});
            }
            setShowSettings(!showSettings);
          }}>
            {showSettings ? "Hide" : "Show"}
          </button>
        </div>
        {showSettings && (
          <div className="row-grid row-grid-single">
            <section className="panel">
              <div className="form-section">
                <label className="form-label">Kraken2 database path</label>
                <input
                  placeholder="/srv/kapurlab/databases/kraken2/k2_standard"
                  value={settingsDraft.kraken_db || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, kraken_db: e.target.value }))}
                />
                <div className="form-hint">Directory containing hash.k2d, opts.k2d, taxo.k2d</div>
              </div>
              <div className="form-section">
                <label className="form-label">BLAST database path or name</label>
                <input
                  placeholder="nt"
                  value={settingsDraft.blast_db || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, blast_db: e.target.value }))}
                />
                <div className="form-hint">Use "nt" for NCBI remote, or an absolute path to a local BLAST db</div>
              </div>
              <div className="form-section">
                <label className="form-label">Personal projects root</label>
                <input
                  value={settingsDraft.projects_root || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, projects_root: e.target.value }))}
                />
                <div className="form-hint">Shared projects at /srv/kapurlab/projects/ are always visible</div>
              </div>
              <div style={{ display: "flex", justifyContent: "flex-end", gap: 8 }}>
                <button onClick={saveSettings}>Save</button>
              </div>
            </section>
          </div>
        )}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Projects & Samples                             */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Projects &amp; Samples</h2>
          <button className="ghost" onClick={() => setShowProjects(!showProjects)}>
            {showProjects ? "Hide" : "Show"}
          </button>
        </div>
        {showProjects && (
          <div className="row-grid row-grid-split">
            {/* LEFT — project / sample browser */}
            <section className="panel">
              <div className="panel-header">
                <h2>Projects</h2>
                <div className="panel-actions">
                  <button className="ghost action" onClick={loadProjects}>↻ Refresh</button>
                </div>
              </div>
              <div className="list project-list">
                {projectsLoading && <div className="loading-text">Loading projects…</div>}
                {!projectsLoading && projects.length === 0 && (
                  <div className="note">No projects found. Check Settings for the projects path.</div>
                )}
                {projects.map((proj) => (
                  <div
                    key={proj.name}
                    className={`list-item ${activeRun?.project === proj.name ? "active" : ""}`}
                  >
                    <div className="item-top" onClick={() => toggleProject(proj.name)}>
                      <span className="expand-icon">{expanded[proj.name] ? "▾" : "▸"}</span>
                      <div className="list-title" title={proj.name}>{proj.name}</div>
                      <span className={`scope-badge scope-${proj.scope}`}>{proj.scope}</span>
                    </div>
                    {proj.path && <div className="list-path" title={proj.path}>{proj.path}</div>}
                    <div className="list-meta">
                      {proj.fastq_count} FASTQ
                      {proj.kraken_runs?.length > 0 &&
                        ` · ${proj.kraken_runs.length} Kraken run${proj.kraken_runs.length > 1 ? "s" : ""}`}
                    </div>
                    {expanded[proj.name] && (
                      <div className="sample-list">
                        {!samples[proj.name] && <div className="loading-text">Loading samples…</div>}
                        {samples[proj.name]?.length === 0 && (
                          <div className="empty-msg" style={{ paddingLeft: 4 }}>No FASTQ files in download/</div>
                        )}
                        {samples[proj.name]?.map((s) => {
                          const key = sampleKey(proj.name, s);
                          const res = sampleResults[key];
                          const vres = vsnpResults[key];
                          const hasRun = proj.kraken_runs?.includes(s.sample);
                          const status = res?.status || (hasRun ? "done" : "none");
                          const checked = !!checkedKeys[key];
                          const open = !!openResults[key];
                          const statusLabel =
                            status === "running" ? "● running" : status === "done" ? "✓ results" : "not run";
                          return (
                          <div
                            key={s.r1}
                            className={`sample-item ${isActive(proj.name, s) ? "active" : ""}`}
                          >
                            <div className="sample-name-row" style={{ display: "flex", alignItems: "center", gap: 8 }}>
                              <input
                                type="checkbox"
                                checked={checked}
                                onChange={() => toggleChecked(proj.name, s)}
                                title="Select for batch run"
                              />
                              <div
                                className="sample-name"
                                title={`${s.sample} — click to show results`}
                                style={{ flex: 1, cursor: "pointer" }}
                                onClick={() => toggleResults(proj.name, s)}
                              >
                                {s.sample}
                              </div>
                              <span className={`read-badge ${s.paired ? "badge-pe" : "badge-se"}`}>
                                {s.paired ? "PE" : "SE"}
                              </span>
                              <span
                                className={`run-status run-status-${status}`}
                                title={`Run status: ${status}`}
                                style={{ fontSize: 11, whiteSpace: "nowrap" }}
                              >
                                {statusLabel}
                              </span>
                              <button
                                className="ghost"
                                style={{ fontSize: 11 }}
                                onClick={() => toggleResults(proj.name, s)}
                                title="Show/hide results for this sample"
                              >
                                {open ? "▾" : "▸"}
                              </button>
                            </div>
                            <div className="sample-files">
                              {s.paired ? (
                                <>
                                  <div className="sample-file-row">
                                    <span className="file-label">R1</span>
                                    <span className="file-name" title={s.r1_name}>{s.r1_name}</span>
                                    <span className="file-size">{fmtSize(s.r1_size)}</span>
                                  </div>
                                  <div className="sample-file-row">
                                    <span className="file-label">R2</span>
                                    <span className="file-name" title={s.r2_name}>{s.r2_name}</span>
                                    <span className="file-size">{fmtSize(s.r2_size)}</span>
                                  </div>
                                </>
                              ) : (
                                <div className="sample-file-row">
                                  <span className="file-name" title={s.r1_name}>{s.r1_name}</span>
                                  <span className="file-size">{fmtSize(s.r1_size)}</span>
                                </div>
                              )}
                            </div>
                            {open && (
                              <div className="sample-results-inline" style={{ marginTop: 6, paddingLeft: 22 }}>
                                <div style={{ display: "flex", gap: 8, marginBottom: 4 }}>
                                  <button
                                    className="ghost action"
                                    disabled={running || (!krakenOnly && !taxon.trim())}
                                    onClick={() => runSamples([{ project: proj.name, ...s }])}
                                    title={!krakenOnly && !taxon.trim() ? "Enter a target taxon first (or tick Kraken only)" : ""}
                                  >
                                    {status === "done" ? "↻ Re-run this sample" : "▶ Run this sample"}
                                  </button>
                                  <button className="ghost action" onClick={() => loadSampleResults(proj.name, s)}>
                                    ↻ Refresh
                                  </button>
                                </div>
                                {res?.loading ? (
                                  <div className="loading-text">Loading results…</div>
                                ) : !res || !res.present || (res.files || []).length === 0 ? (
                                  <div className="empty-msg" style={{ paddingLeft: 0 }}>
                                    {status === "running"
                                      ? "Running… results will appear here when finished."
                                      : "No Kraken results yet for this sample."}
                                  </div>
                                ) : (
                                  <div className="results-list">
                                    {res.files.map((f) => {
                                      const base = `./api/projects/${encodeURIComponent(proj.name)}/file?path=${encodeURIComponent(f.path)}`;
                                      return (
                                        <div key={f.name} className="results-item">
                                          <span className="result-icon">{fileIcon(f.name)}</span>
                                          {f.openable ? (
                                            <a className="result-name result-link" href={`${base}&inline=1`}
                                               target="_blank" rel="noopener noreferrer" title={`Open ${f.name}`}>
                                              {f.label || f.name}
                                            </a>
                                          ) : (
                                            <a className="result-name result-link" href={`${base}&inline=0`}
                                               title={`Download ${f.name}`}>
                                              {f.label || f.name}
                                            </a>
                                          )}
                                          <span className="result-size">{fmtSize(f.size)}</span>
                                          <a className="result-download" href={`${base}&inline=0`} title={`Download ${f.name}`}>⬇</a>
                                        </div>
                                      );
                                    })}
                                  </div>
                                )}

                                {/* Cross-tool: vSNP results for this sample */}
                                <div className="vsnp-cross-tool" style={{ marginTop: 10, borderTop: "1px solid var(--border, #e2e2e2)", paddingTop: 8 }}>
                                  <div style={{ fontWeight: 600, fontSize: 12, marginBottom: 4 }}>vSNP results</div>
                                  {vres?.loading ? (
                                    <div className="loading-text">Loading vSNP results…</div>
                                  ) : !vres || !vres.step1_present ? (
                                    <div className="empty-msg" style={{ paddingLeft: 0 }}>No vSNP run for this sample yet.</div>
                                  ) : (
                                    <>
                                      <div className="results-list">
                                        {(vres.files || []).map((f) => {
                                          const vbase = `./api/projects/${encodeURIComponent(proj.name)}/file?path=${encodeURIComponent(f.path)}`;
                                          return (
                                            <div key={f.relpath} className="results-item">
                                              <span className="result-icon">{fileIcon(f.name)}</span>
                                              {f.openable ? (
                                                <a className="result-name result-link" href={`${vbase}&inline=1`}
                                                   target="_blank" rel="noopener noreferrer" title={`Open ${f.name}`}>
                                                  {f.relpath}
                                                </a>
                                              ) : (
                                                <a className="result-name result-link" href={`${vbase}&inline=0`}
                                                   title={`Download ${f.name}`}>
                                                  {f.relpath}
                                                </a>
                                              )}
                                              <span className="result-size">{fmtSize(f.size)}</span>
                                              <a className="result-download" href={`${vbase}&inline=0`} title={`Download ${f.name}`}>⬇</a>
                                            </div>
                                          );
                                        })}
                                      </div>
                                      {vres.step2 && vres.step2.present && (
                                        <div className="vsnp-step2" style={{ marginTop: 6 }}>
                                          {vres.step2.report_path ? (
                                            <a className="result-name result-link"
                                               href={`./api/projects/${encodeURIComponent(proj.name)}/file?path=${encodeURIComponent(vres.step2.report_path)}&inline=1`}
                                               target="_blank" rel="noopener noreferrer"
                                               title="Open the latest SNP comparison report this sample appears in">
                                              📊 Latest SNP comparison{vres.step2.started_at ? ` (${vres.step2.started_at})` : ""}
                                            </a>
                                          ) : (
                                            <span className="muted">
                                              In latest SNP comparison{vres.step2.started_at ? ` (${vres.step2.started_at})` : ""}
                                            </span>
                                          )}
                                          {vres.step2.groups && vres.step2.groups.length > 0 && (
                                            <span className="muted" style={{ marginLeft: 6 }}>
                                              — group{vres.step2.groups.length > 1 ? "s" : ""}: {vres.step2.groups.join(", ")}
                                            </span>
                                          )}
                                        </div>
                                      )}
                                    </>
                                  )}
                                </div>
                              </div>
                            )}
                          </div>
                          );
                        })}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </section>

            {/* RIGHT — samples selected for a batch run */}
            <section className="panel">
              <div className="panel-header">
                <h2>Selected for run</h2>
                {Object.keys(checkedKeys).length > 0 && (
                  <button className="ghost action" onClick={() => setCheckedKeys({})}>Clear</button>
                )}
              </div>
              {Object.keys(checkedKeys).length === 0 ? (
                <div className="empty-msg">
                  Check one or more samples on the left, then run them as a batch from “Run Kraken” below.
                  Click a sample’s name to view its results inline.
                </div>
              ) : (
                <div className="selection-box">
                  <div className="sel-title">{Object.keys(checkedKeys).length} sample(s) queued</div>
                  {Object.entries(checkedKeys).map(([key, samp]) => (
                    <div key={key} className="sel-row" style={{ display: "flex", alignItems: "center", gap: 8 }}>
                      <span className="sel-name" style={{ flex: 1 }}>{samp.sample}</span>
                      <span className="muted" style={{ fontSize: 11 }}>{samp.project}</span>
                      <button className="ghost" style={{ fontSize: 11 }}
                              onClick={() => toggleChecked(samp.project, samp)} title="Remove from batch">✕</button>
                    </div>
                  ))}
                </div>
              )}
            </section>
          </div>
        )}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Run Kraken — configure + results               */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Run Kraken</h2>
          <button className="ghost" onClick={() => setShowRun(!showRun)}>
            {showRun ? "Hide" : "Show"}
          </button>
        </div>
        {showRun && (
          <div className="row-grid row-grid-split">
            {/* LEFT — configure & run */}
            <section className="panel">
              <h2>Configure &amp; Run</h2>

              <div className="form-section">
                <label className="checkbox-label" style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
                  <input
                    type="checkbox"
                    checked={krakenOnly}
                    onChange={(e) => setKrakenOnly(e.target.checked)}
                    disabled={running}
                  />
                  <span>Kraken only (Krona graph, no read parsing)</span>
                </label>
                <div className="note" style={{ marginTop: 4 }}>
                  Runs Kraken2 and produces the Krona graph only — skips read parsing, assembly, and BLAST. No target taxon needed.
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">Target Taxon</label>
                <input
                  placeholder='e.g. "Mycobacterium tuberculosis complex"'
                  value={taxon}
                  onChange={(e) => setTaxon(e.target.value)}
                  disabled={running || krakenOnly}
                />
                <div className="taxon-presets">
                  {TAXON_PRESETS.map((p) => (
                    <button
                      key={p}
                      className={`preset-btn ${taxon === p ? "active" : ""}`}
                      onClick={() => setTaxon(p)}
                      disabled={running || krakenOnly}
                    >
                      {p}
                    </button>
                  ))}
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">
                  Kraken2 DB path
                  {!krakenDb && <span style={{ color: "var(--danger)", marginLeft: 6, fontSize: 11 }}>⚠ not configured</span>}
                </label>
                <input
                  placeholder="/srv/kapurlab/databases/kraken2/k2_standard"
                  value={krakenDb}
                  onChange={(e) => setKrakenDb(e.target.value)}
                  disabled={running}
                />
              </div>

              <div className="form-section">
                <label className="form-label">BLAST DB path (or name)</label>
                <input
                  placeholder="nt  or  /srv/kapurlab/databases/blast/nt"
                  value={blastDb}
                  onChange={(e) => setBlastDb(e.target.value)}
                  disabled={running || krakenOnly}
                />
              </div>

              <button
                className="run-btn"
                onClick={runSelected}
                disabled={running || Object.keys(checkedKeys).length === 0 || (!krakenOnly && !taxon.trim())}
              >
                {running
                  ? `Running… ${queueInfo.total > 1 ? `(${queueInfo.done}/${queueInfo.total})` : ""}`
                  : `▶ Run selected${Object.keys(checkedKeys).length ? ` (${Object.keys(checkedKeys).length})` : ""}`}
              </button>
              {Object.keys(checkedKeys).length === 0 && (
                <div className="note">Check one or more samples on the left to enable the run. (Or use “Run this sample” under any sample.)</div>
              )}
              {!krakenOnly && !taxon.trim() && Object.keys(checkedKeys).length > 0 && (
                <div className="note">Enter a target taxon above (or tick “Kraken only”) to enable the run.</div>
              )}
            </section>

            {/* RIGHT — current run status (per-sample results live inline at left) */}
            <section className="panel">
              <div className="panel-header">
                <h2>Current run</h2>
                {jobId && <span className="muted" style={{ fontSize: 12 }}>job {jobId.slice(0, 8)}</span>}
              </div>
              {activeRun ? (
                <div className="selection-box">
                  <div className="sel-title">
                    {jobStatus === "running" ? "Running" : jobStatus === "succeeded" ? "Done" : jobStatus}
                    {queueInfo.total > 1 ? ` — ${queueInfo.done}/${queueInfo.total} in batch` : ""}
                  </div>
                  <div><span className="sel-name">{activeRun.sample}</span></div>
                  <div style={{ marginTop: 2 }}>
                    <span className="muted">Project:</span> <strong>{activeRun.project}</strong>
                  </div>
                  {currentStep && <div className="muted" style={{ marginTop: 4 }}>{currentStep}</div>}
                  <div className="note" style={{ marginTop: 8 }}>
                    Output files appear inline under each sample on the left (click a sample’s name to expand).
                  </div>
                </div>
              ) : (
                <div className="empty-msg">
                  No active run. Results for any sample are shown inline under that sample on the left.
                </div>
              )}
            </section>
          </div>
        )}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Pipeline Log                                   */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Pipeline Log</h2>
          <button className="ghost" onClick={() => setShowLogs(!showLogs)}>
            {showLogs ? "Hide" : "Show"}
          </button>
        </div>
        {showLogs && (
          <div className="row-grid row-grid-single">
            <section className="panel">
              <div className="log-meta">
                <span className="dot" data-state={jobStatus} />
                <span style={{ fontWeight: 600 }}>
                  {jobStatus === "idle" && "Idle"}
                  {jobStatus === "running" && "Running"}
                  {jobStatus === "succeeded" && "Done"}
                  {jobStatus === "failed" && "Failed"}
                </span>
                {jobStatus === "running" && currentStep && (
                  <span className="log-step" title={currentStep}>— {currentStep}</span>
                )}
              </div>
              <div className="log" ref={logRef}>
                {logLines.length === 0 ? (
                  <span className="log-placeholder">
                    {jobStatus === "idle"
                      ? "Select a sample and click Run to start."
                      : "Waiting for output…"}
                  </span>
                ) : (
                  logLines.map((line, i) => (
                    <div key={i} className={logLineClass(line)}>{line}</div>
                  ))
                )}
              </div>
            </section>
          </div>
        )}
      </main>
    </div>
  );
}
