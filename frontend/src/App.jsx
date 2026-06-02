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
  const [selectedSample, setSelectedSample] = useState(null); // {project, sample, r1, r2}
  const [taxon, setTaxon] = useState("");
  const [krakenDb, setKrakenDb] = useState("");
  const [blastDb, setBlastDb] = useState("nt");
  const [running, setRunning] = useState(false);
  const [jobId, setJobId] = useState(null);
  const [jobStatus, setJobStatus] = useState("idle"); // idle | running | succeeded | failed
  const [logLines, setLogLines] = useState([]);
  const [results, setResults] = useState([]);
  const [vsnpResults, setVsnpResults] = useState(null);  // null=loading; {step1_present, files, step2}
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
          streamLog(live.id);
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

  function selectSample(project, sample) {
    setSelectedSample({ project, ...sample });
    // Cross-tool: pull whatever vSNP produced for this sample (step1 files +
    // the latest step2 comparison it appears in). Both tools share the same
    // project dir, so these already exist on disk if vSNP was run.
    setVsnpResults(null);
    fetch(`./api/projects/${encodeURIComponent(project)}/vsnp/samples/${encodeURIComponent(sample.sample)}/files`)
      .then((r) => r.json())
      .then(setVsnpResults)
      .catch(() => setVsnpResults({ step1_present: false, files: [], step2: { present: false } }));
  }

  function startRun() {
    if (!selectedSample || !taxon.trim()) return;
    // Close any existing SSE stream
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
      eventSourceRef.current = null;
    }
    setRunning(true);
    setJobStatus("running");
    setLogLines([]);
    setResults([]);
    setCurrentStep("");
    setShowLogs(true);

    fetch("./api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        project: selectedSample.project,
        r1: selectedSample.r1,
        r2: selectedSample.r2 || null,
        taxon: taxon.trim(),
        kraken_db: krakenDb.trim() || null,
        blast_db: blastDb.trim() || null,
      }),
    })
      .then((r) => {
        if (!r.ok) return r.json().then((e) => { throw new Error(e.detail || "Run failed"); });
        return r.json();
      })
      .then(({ job_id }) => {
        setJobId(job_id);
        streamLog(job_id);
      })
      .catch((err) => {
        setLogLines([`ERROR: ${err.message}`]);
        setRunning(false);
        setJobStatus("failed");
      });
  }

  function streamLog(id) {
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
            if (job.status === "succeeded") loadResults(id);
          });
      } else {
        setLogLines((prev) => [...prev, data]);
        // Extract current pipeline step for the header
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
    };
  }

  function loadResults(id) {
    fetch(`./api/jobs/${id}/results`)
      .then((r) => r.json())
      .then(setResults)
      .catch(() => {});
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
            <span className="status-label">Project</span>
            <span className="status-value">{selectedSample?.project || "—"}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Sample</span>
            <span className="status-value">{selectedSample?.sample || "—"}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Reads</span>
            <span className="status-value">
              {selectedSample ? (selectedSample.paired ? "Paired-end" : "Single-end") : "—"}
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
                    className={`list-item ${selectedSample?.project === proj.name ? "active" : ""}`}
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
                        {samples[proj.name]?.map((s) => (
                          <div
                            key={s.r1}
                            className={`sample-item ${selectedSample?.r1 === s.r1 ? "active" : ""}`}
                            onClick={() => selectSample(proj.name, s)}
                          >
                            <div className="sample-name-row">
                              <div className="sample-name" title={s.sample}>{s.sample}</div>
                              <span className={`read-badge ${s.paired ? "badge-pe" : "badge-se"}`}>
                                {s.paired ? "PE" : "SE"}
                              </span>
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
                          </div>
                        ))}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </section>

            {/* RIGHT — selected sample details */}
            <section className="panel">
              <h2>Selected Sample</h2>
              {selectedSample ? (
                <div className="selection-box">
                  <div className="sel-title">Ready to run</div>
                  <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                    <span className="sel-name">{selectedSample.sample}</span>
                    <span className={`read-badge ${selectedSample.paired ? "badge-pe" : "badge-se"}`}>
                      {selectedSample.paired ? "Paired-end" : "Single-end"}
                    </span>
                  </div>
                  <div className="sel-row">
                    {selectedSample.paired && <span className="file-label">R1</span>}
                    <span className="file-name" title={selectedSample.r1_name}>{selectedSample.r1_name}</span>
                  </div>
                  {selectedSample.r2_name && (
                    <div className="sel-row">
                      <span className="file-label">R2</span>
                      <span className="file-name" title={selectedSample.r2_name}>{selectedSample.r2_name}</span>
                    </div>
                  )}
                  <div style={{ marginTop: 2 }}>
                    <span className="muted">Project:</span> <strong>{selectedSample.project}</strong>
                  </div>
                </div>
              ) : (
                <div className="empty-msg">
                  Select a sample from a project on the left to configure and run Kraken classification.
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
                <label className="form-label">Target Taxon</label>
                <input
                  placeholder='e.g. "Mycobacterium tuberculosis complex"'
                  value={taxon}
                  onChange={(e) => setTaxon(e.target.value)}
                  disabled={running}
                />
                <div className="taxon-presets">
                  {TAXON_PRESETS.map((p) => (
                    <button
                      key={p}
                      className={`preset-btn ${taxon === p ? "active" : ""}`}
                      onClick={() => setTaxon(p)}
                      disabled={running}
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
                  disabled={running}
                />
              </div>

              <button
                className="run-btn"
                onClick={startRun}
                disabled={running || !selectedSample || !taxon.trim()}
              >
                {running ? "Running…" : "▶ Run Kraken ID Parse"}
              </button>
              {!selectedSample && (
                <div className="note">Select a sample first to enable the run.</div>
              )}
            </section>

            {/* RIGHT — results */}
            <section className="panel">
              <div className="panel-header">
                <h2>Results</h2>
                {jobId && <span className="muted" style={{ fontSize: 12 }}>job {jobId.slice(0, 8)}</span>}
              </div>
              {results.length === 0 ? (
                <div className="empty-msg">
                  {jobStatus === "succeeded" ? "No output files found." : "Run a sample to see results here."}
                </div>
              ) : (
                <div className="results-list">
                  {results.map((f) => {
                    const base = `./api/jobs/${jobId}/file?path=${encodeURIComponent(f.name)}`;
                    const displayName = f.label || f.name;
                    return (
                      <div key={f.name} className="results-item">
                        <span className="result-icon">{fileIcon(f.name)}</span>
                        {f.openable ? (
                          <a className="result-name result-link" href={`${base}&inline=1`}
                             target="_blank" rel="noopener noreferrer" title={`Open ${f.name}`}>
                            {displayName}
                          </a>
                        ) : (
                          <a className="result-name result-link" href={`${base}&inline=0`}
                             title={`Download ${f.name}`}>
                            {displayName}
                          </a>
                        )}
                        <span className="result-size">{fmtSize(f.size)}</span>
                        <a className="result-download" href={`${base}&inline=0`} title={`Download ${f.name}`}>⬇</a>
                      </div>
                    );
                  })}
                </div>
              )}

              {/* Cross-tool: vSNP results for the selected sample */}
              {selectedSample && (
                <div className="vsnp-cross-tool" style={{ marginTop: 18, borderTop: "1px solid var(--border, #e2e2e2)", paddingTop: 12 }}>
                  <h3 style={{ marginBottom: 6 }}>vSNP results — {selectedSample.sample}</h3>
                  {vsnpResults === null ? (
                    <div className="muted">Loading…</div>
                  ) : (
                    <>
                      {!vsnpResults.step1_present ? (
                        <div className="muted">No vSNP run for this sample yet.</div>
                      ) : (
                        <div className="results-list">
                          {vsnpResults.files.map((f) => {
                            const vbase = `./api/projects/${encodeURIComponent(selectedSample.project)}/file?path=${encodeURIComponent(f.path)}`;
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
                      )}
                      {vsnpResults.step2 && vsnpResults.step2.present && (
                        <div className="vsnp-step2" style={{ marginTop: 8 }}>
                          {vsnpResults.step2.report_path ? (
                            <a className="result-name result-link"
                               href={`./api/projects/${encodeURIComponent(selectedSample.project)}/file?path=${encodeURIComponent(vsnpResults.step2.report_path)}&inline=1`}
                               target="_blank" rel="noopener noreferrer"
                               title="Open the latest SNP comparison report this sample appears in">
                              📊 Latest SNP comparison{vsnpResults.step2.started_at ? ` (${vsnpResults.step2.started_at})` : ""}
                            </a>
                          ) : (
                            <span className="muted">
                              In latest SNP comparison{vsnpResults.step2.started_at ? ` (${vsnpResults.step2.started_at})` : ""}
                            </span>
                          )}
                          {vsnpResults.step2.groups && vsnpResults.step2.groups.length > 0 && (
                            <span className="muted" style={{ marginLeft: 6 }}>
                              — group{vsnpResults.step2.groups.length > 1 ? "s" : ""}: {vsnpResults.step2.groups.join(", ")}
                            </span>
                          )}
                        </div>
                      )}
                    </>
                  )}
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
