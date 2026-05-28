import { useState, useEffect, useRef, useCallback } from "react";

// ---------------------------------------------------------------------------
// Styles (injected as a <style> tag so there's no separate CSS file to build)
// ---------------------------------------------------------------------------
const CSS = `
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'IBM Plex Sans', system-ui, sans-serif; font-size: 14px;
  background: #f4f6f8; color: #1a2536; height: 100vh; overflow: hidden; }
#root { height: 100vh; display: flex; flex-direction: column; }

/* Header */
.header { background: #1c3754; color: #fff; padding: 0 16px;
  display: flex; align-items: center; gap: 12px; height: 48px; flex-shrink: 0; }
.header h1 { font-size: 16px; font-weight: 600; letter-spacing: 0.02em; }
.header-badge { background: #b85a3e; color: #fff; font-size: 11px;
  padding: 2px 8px; border-radius: 10px; font-weight: 600; }
.header-spacer { flex: 1; }
.header-settings { background: none; border: 1px solid rgba(255,255,255,0.3);
  color: #fff; padding: 4px 12px; border-radius: 4px; cursor: pointer; font-size: 13px; }
.header-settings:hover { background: rgba(255,255,255,0.1); }

/* Layout */
.layout { display: flex; flex: 1; overflow: hidden; gap: 0; }

/* Left panel — project/sample browser */
.panel-left { width: 280px; flex-shrink: 0; background: #fff;
  border-right: 1px solid #dde2e8; display: flex; flex-direction: column; overflow: hidden; }
.panel-title { padding: 10px 14px; font-size: 12px; font-weight: 700;
  text-transform: uppercase; letter-spacing: 0.06em; color: #64748b;
  border-bottom: 1px solid #e8ecf0; background: #f8fafc; flex-shrink: 0; }
.project-list { flex: 1; overflow-y: auto; }
.project-item { border-bottom: 1px solid #e8ecf0; }
.project-header { padding: 10px 14px; cursor: pointer; display: flex;
  align-items: center; gap: 8px; transition: background 0.1s; }
.project-header:hover { background: #f0f4f8; }
.project-header.selected { background: #e8f0fe; }
.project-name { font-weight: 600; font-size: 13px; flex: 1; min-width: 0;
  overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.project-meta { font-size: 11px; color: #64748b; }
.scope-badge { font-size: 10px; padding: 1px 6px; border-radius: 8px;
  font-weight: 600; flex-shrink: 0; }
.scope-shared { background: #dbeafe; color: #1d4ed8; }
.scope-personal { background: #f0fdf4; color: #166534; }
.expand-icon { font-size: 10px; color: #94a3b8; flex-shrink: 0; }
.sample-list { background: #f8fafc; }
.sample-item { padding: 8px 14px 8px 28px; cursor: pointer; border-top: 1px solid #e8ecf0;
  transition: background 0.1s; }
.sample-item:hover { background: #eff6ff; }
.sample-item.active { background: #dbeafe; border-left: 3px solid #2563eb; padding-left: 25px; }
.sample-name { font-size: 12px; font-weight: 500; overflow: hidden;
  text-overflow: ellipsis; white-space: nowrap; }
.sample-files { font-size: 11px; color: #94a3b8; margin-top: 2px; }
.no-projects { padding: 24px 14px; color: #64748b; font-size: 13px; text-align: center; }
.loading-text { padding: 16px; color: #64748b; font-size: 13px; }

/* Center panel — run form */
.panel-center { flex: 1; display: flex; flex-direction: column; overflow: hidden; }
.run-form { padding: 16px; border-bottom: 1px solid #dde2e8; background: #fff;
  flex-shrink: 0; overflow-y: auto; max-height: 60%; }
.form-section { margin-bottom: 16px; }
.form-label { font-size: 12px; font-weight: 600; color: #374151; margin-bottom: 4px; display: block; }
.form-input { width: 100%; padding: 7px 10px; border: 1px solid #cbd5e1;
  border-radius: 4px; font-size: 13px; outline: none; transition: border-color 0.15s; }
.form-input:focus { border-color: #2563eb; box-shadow: 0 0 0 2px rgba(37,99,235,0.1); }
.form-input:disabled { background: #f1f5f9; color: #94a3b8; cursor: not-allowed; }
.taxon-presets { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 6px; }
.preset-btn { padding: 3px 10px; font-size: 11px; border: 1px solid #cbd5e1;
  border-radius: 12px; cursor: pointer; background: #f8fafc; transition: all 0.1s; }
.preset-btn:hover { background: #dbeafe; border-color: #93c5fd; }
.preset-btn.active { background: #2563eb; color: #fff; border-color: #2563eb; }
.selected-sample-box { background: #f0f7ff; border: 1px solid #bfdbfe;
  border-radius: 6px; padding: 10px 12px; margin-bottom: 16px; }
.selected-sample-title { font-size: 11px; font-weight: 700; color: #1d4ed8;
  text-transform: uppercase; letter-spacing: 0.06em; margin-bottom: 4px; }
.selected-sample-name { font-size: 14px; font-weight: 600; }
.selected-sample-path { font-size: 11px; color: #64748b; margin-top: 2px;
  overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.no-sample-msg { font-size: 13px; color: #94a3b8; font-style: italic; margin-bottom: 16px; }
.run-btn { width: 100%; padding: 10px; background: #1c3754; color: #fff;
  border: none; border-radius: 6px; font-size: 14px; font-weight: 600; cursor: pointer;
  transition: background 0.15s; }
.run-btn:hover:not(:disabled) { background: #2563eb; }
.run-btn:disabled { background: #94a3b8; cursor: not-allowed; }

/* Log panel */
.log-panel { flex: 1; display: flex; flex-direction: column; overflow: hidden;
  background: #1a2536; }
.log-header { padding: 8px 14px; display: flex; align-items: center; gap: 8px;
  background: #0f172a; flex-shrink: 0; }
.log-title { font-size: 12px; font-weight: 600; color: #94a3b8;
  text-transform: uppercase; letter-spacing: 0.06em; }
.status-dot { width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }
.status-running { background: #f59e0b; animation: pulse 1s infinite; }
.status-succeeded { background: #22c55e; }
.status-failed { background: #ef4444; }
.status-idle { background: #475569; }
@keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.4} }
.log-content { flex: 1; overflow-y: auto; padding: 12px 14px; font-family: 'IBM Plex Mono', monospace;
  font-size: 12px; line-height: 1.6; color: #94a3b8; }
.log-line { white-space: pre-wrap; word-break: break-all; }
.log-line.cmd { color: #60a5fa; }
.log-line.done { color: #22c55e; font-weight: 600; }
.log-line.error { color: #f87171; }
.log-placeholder { color: #475569; font-style: italic; }

/* Right panel — results */
.panel-right { width: 260px; flex-shrink: 0; background: #fff;
  border-left: 1px solid #dde2e8; display: flex; flex-direction: column; overflow: hidden; }
.results-list { flex: 1; overflow-y: auto; }
.result-file { padding: 8px 14px; border-bottom: 1px solid #e8ecf0;
  display: flex; align-items: flex-start; gap: 8px; font-size: 12px; }
.result-icon { flex-shrink: 0; font-size: 14px; }
.result-name { word-break: break-all; flex: 1; color: #1a2536; }
.result-size { font-size: 11px; color: #94a3b8; flex-shrink: 0; }
.no-results { padding: 24px 14px; color: #94a3b8; font-size: 12px; text-align: center; font-style: italic; }

/* Settings modal */
.modal-overlay { position: fixed; inset: 0; background: rgba(0,0,0,0.4);
  display: flex; align-items: center; justify-content: center; z-index: 100; }
.modal { background: #fff; border-radius: 8px; padding: 24px; width: 500px;
  max-width: 90vw; box-shadow: 0 20px 60px rgba(0,0,0,0.2); }
.modal h2 { font-size: 16px; font-weight: 700; margin-bottom: 20px; }
.modal-actions { display: flex; gap: 10px; justify-content: flex-end; margin-top: 20px; }
.btn-cancel { padding: 8px 16px; background: #f1f5f9; border: 1px solid #cbd5e1;
  border-radius: 4px; cursor: pointer; font-size: 13px; }
.btn-save { padding: 8px 16px; background: #1c3754; color: #fff;
  border: none; border-radius: 4px; cursor: pointer; font-size: 13px; font-weight: 600; }
.btn-save:hover { background: #2563eb; }
.form-hint { font-size: 11px; color: #64748b; margin-top: 4px; }
.refresh-btn { padding: 3px 8px; font-size: 11px; background: none;
  border: 1px solid #cbd5e1; border-radius: 4px; cursor: pointer; color: #64748b; margin-left: auto; }
.refresh-btn:hover { background: #f0f4f8; }
`;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const TAXON_PRESETS = [
  "Mycobacterium tuberculosis complex",
  "Mycobacterium bovis",
  "Orbivirus",
  "Apicomplexa",
  "Isavirus salaris",
];

function fileIcon(name) {
  if (name.endsWith(".xlsx")) return "📊";
  if (name.endsWith(".pdf")) return "📄";
  if (name.endsWith(".png") || name.endsWith(".pdf")) return "🖼";
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
  const [showSettings, setShowSettings] = useState(false);
  const [settingsDraft, setSettingsDraft] = useState({});
  const logRef = useRef(null);
  const eventSourceRef = useRef(null);

  // Load config & projects on mount
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

  // Scroll log to bottom when new lines arrive
  useEffect(() => {
    if (logRef.current) {
      logRef.current.scrollTop = logRef.current.scrollHeight;
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
            if (job.status === "succeeded") loadResults(id);
          });
      } else {
        setLogLines((prev) => [...prev, data]);
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
        setShowSettings(false);
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

  return (
    <>
      <style>{CSS}</style>
      <div id="root" style={{ display: "flex", flexDirection: "column", height: "100vh" }}>
        {/* Header */}
        <header className="header">
          <span style={{ fontSize: 20 }}>🦠</span>
          <h1>Kraken ID Parse</h1>
          <span className="header-badge">contamination screen</span>
          <div className="header-spacer" />
          <button className="header-settings" onClick={() => {
            fetch("./api/config").then(r => r.json()).then(cfg => setSettingsDraft(cfg));
            setShowSettings(true);
          }}>
            ⚙ Settings
          </button>
        </header>

        {/* Main layout */}
        <div className="layout">
          {/* LEFT — project browser */}
          <aside className="panel-left">
            <div className="panel-title" style={{ display: "flex", alignItems: "center" }}>
              Projects
              <button className="refresh-btn" onClick={loadProjects}>↻ refresh</button>
            </div>
            <div className="project-list">
              {projectsLoading && <div className="loading-text">Loading projects…</div>}
              {!projectsLoading && projects.length === 0 && (
                <div className="no-projects">
                  No projects found.<br />
                  <span style={{ fontSize: 11 }}>Check Settings for the projects path.</span>
                </div>
              )}
              {projects.map((proj) => (
                <div key={proj.name} className="project-item">
                  <div
                    className={`project-header ${selectedSample?.project === proj.name ? "selected" : ""}`}
                    onClick={() => toggleProject(proj.name)}
                  >
                    <span className="expand-icon">{expanded[proj.name] ? "▾" : "▸"}</span>
                    <span className="project-name" title={proj.name}>{proj.name}</span>
                    <span className={`scope-badge scope-${proj.scope}`}>{proj.scope}</span>
                  </div>
                  <div className="project-meta" style={{ paddingLeft: 32, paddingBottom: expanded[proj.name] ? 0 : 6, fontSize: 11, color: "#94a3b8" }}>
                    {proj.fastq_count} FASTQ
                    {proj.kraken_runs?.length > 0 && ` · ${proj.kraken_runs.length} Kraken run${proj.kraken_runs.length > 1 ? "s" : ""}`}
                  </div>
                  {expanded[proj.name] && (
                    <div className="sample-list">
                      {!samples[proj.name] && <div className="loading-text">Loading samples…</div>}
                      {samples[proj.name]?.length === 0 && (
                        <div style={{ padding: "8px 14px 8px 28px", fontSize: 12, color: "#94a3b8" }}>
                          No FASTQ files in download/
                        </div>
                      )}
                      {samples[proj.name]?.map((s) => (
                        <div
                          key={s.sample}
                          className={`sample-item ${selectedSample?.r1 === s.r1 ? "active" : ""}`}
                          onClick={() => selectSample(proj.name, s)}
                        >
                          <div className="sample-name" title={s.sample}>{s.sample}</div>
                          <div className="sample-files">
                            R1 {s.r1_name}
                            {s.r2_name && <><br />R2 {s.r2_name}</>}
                          </div>
                        </div>
                      ))}
                    </div>
                  )}
                </div>
              ))}
            </div>
          </aside>

          {/* CENTER — run form + log */}
          <main className="panel-center">
            {/* Run form */}
            <div className="run-form">
              {selectedSample ? (
                <div className="selected-sample-box">
                  <div className="selected-sample-title">Selected Sample</div>
                  <div className="selected-sample-name">{selectedSample.sample}</div>
                  <div className="selected-sample-path">{selectedSample.r1_name}</div>
                  {selectedSample.r2_name && (
                    <div className="selected-sample-path">{selectedSample.r2_name}</div>
                  )}
                  <div style={{ fontSize: 11, color: "#1d4ed8", marginTop: 4 }}>
                    Project: {selectedSample.project}
                  </div>
                </div>
              ) : (
                <div className="no-sample-msg">
                  ← Select a sample from a project to run Kraken classification
                </div>
              )}

              <div className="form-section">
                <label className="form-label">Target Taxon</label>
                <input
                  className="form-input"
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
                  {!krakenDb && <span style={{ color: "#b85a3e", marginLeft: 6, fontSize: 11 }}>⚠ not configured</span>}
                </label>
                <input
                  className="form-input"
                  placeholder="/srv/kapurlab/databases/kraken2/k2_standard"
                  value={krakenDb}
                  onChange={(e) => setKrakenDb(e.target.value)}
                  disabled={running}
                />
              </div>

              <div className="form-section">
                <label className="form-label">BLAST DB path (or name)</label>
                <input
                  className="form-input"
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
            </div>

            {/* Log stream */}
            <div className="log-panel">
              <div className="log-header">
                <div className={`status-dot status-${jobStatus}`} />
                <span className="log-title">
                  {jobStatus === "idle" && "Log output"}
                  {jobStatus === "running" && "Running…"}
                  {jobStatus === "succeeded" && "Completed"}
                  {jobStatus === "failed" && "Failed"}
                </span>
                {jobId && (
                  <span style={{ fontSize: 11, color: "#475569", marginLeft: 8 }}>
                    job {jobId.slice(0, 8)}
                  </span>
                )}
              </div>
              <div className="log-content" ref={logRef}>
                {logLines.length === 0 && (
                  <span className="log-placeholder">
                    {jobStatus === "idle"
                      ? "Select a sample and click Run to start."
                      : "Waiting for output…"}
                  </span>
                )}
                {logLines.map((line, i) => (
                  <div key={i} className={logLineClass(line)}>{line}</div>
                ))}
              </div>
            </div>
          </main>

          {/* RIGHT — results */}
          <aside className="panel-right">
            <div className="panel-title">Results</div>
            <div className="results-list">
              {results.length === 0 ? (
                <div className="no-results">
                  {jobStatus === "succeeded"
                    ? "No output files found."
                    : "Run a sample to see results here."}
                </div>
              ) : (
                results.map((f) => (
                  <div key={f.name} className="result-file">
                    <span className="result-icon">{fileIcon(f.name)}</span>
                    <span className="result-name">{f.name}</span>
                    <span className="result-size">{fmtSize(f.size)}</span>
                  </div>
                ))
              )}
            </div>
          </aside>
        </div>

        {/* Settings modal */}
        {showSettings && (
          <div className="modal-overlay" onClick={() => setShowSettings(false)}>
            <div className="modal" onClick={(e) => e.stopPropagation()}>
              <h2>⚙ Settings</h2>

              <div className="form-section">
                <label className="form-label">Kraken2 database path</label>
                <input
                  className="form-input"
                  placeholder="/srv/kapurlab/databases/kraken2/k2_standard"
                  value={settingsDraft.kraken_db || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, kraken_db: e.target.value }))}
                />
                <div className="form-hint">
                  Directory containing hash.k2d, opts.k2d, taxo.k2d
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">BLAST database path or name</label>
                <input
                  className="form-input"
                  placeholder="nt"
                  value={settingsDraft.blast_db || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, blast_db: e.target.value }))}
                />
                <div className="form-hint">
                  Use "nt" for NCBI remote, or an absolute path to a local BLAST db
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">Personal projects root</label>
                <input
                  className="form-input"
                  value={settingsDraft.projects_root || ""}
                  onChange={(e) => setSettingsDraft((d) => ({ ...d, projects_root: e.target.value }))}
                />
                <div className="form-hint">
                  Shared projects at /srv/kapurlab/projects/ are always visible
                </div>
              </div>

              <div className="modal-actions">
                <button className="btn-cancel" onClick={() => setShowSettings(false)}>Cancel</button>
                <button className="btn-save" onClick={saveSettings}>Save</button>
              </div>
            </div>
          </div>
        )}
      </div>
    </>
  );
}
