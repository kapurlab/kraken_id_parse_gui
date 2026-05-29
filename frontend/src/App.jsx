import { useState, useEffect, useRef } from "react";
import "./App.css";


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
  const [currentStep, setCurrentStep] = useState("");
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
      <div className="app-root">
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
                    <div className="project-name-block">
                      <div className="project-name" title={proj.name}>{proj.name}</div>
                      {proj.path && (
                        <div className="project-path" title={proj.path}>{proj.path}</div>
                      )}
                    </div>
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
                                  <span style={{ fontSize: 10, color: "#64748b", flexShrink: 0 }}>{fmtSize(s.r1_size)}</span>
                                </div>
                                <div className="sample-file-row">
                                  <span className="file-label">R2</span>
                                  <span className="file-name" title={s.r2_name}>{s.r2_name}</span>
                                  <span style={{ fontSize: 10, color: "#64748b", flexShrink: 0 }}>{fmtSize(s.r2_size)}</span>
                                </div>
                              </>
                            ) : (
                              <div className="sample-file-row">
                                <span className="file-name" title={s.r1_name}>{s.r1_name}</span>
                                <span style={{ fontSize: 10, color: "#64748b", flexShrink: 0 }}>{fmtSize(s.r1_size)}</span>
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
          </aside>

          {/* CENTER — run form + log */}
          <main className="panel-center">
            {/* Run form */}
            <div className="run-form">
              {selectedSample ? (
                <div className="selected-sample-box">
                  <div className="selected-sample-title">Selected Sample</div>
                  <div style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 6 }}>
                    <div className="selected-sample-name">{selectedSample.sample}</div>
                    <span className={`read-badge ${selectedSample.paired ? "badge-pe" : "badge-se"}`}>
                      {selectedSample.paired ? "Paired-end" : "Single-end"}
                    </span>
                  </div>
                  <div style={{ fontSize: 11, color: "#475569" }}>
                    <div className="sample-file-row" style={{ marginBottom: 1 }}>
                      {selectedSample.paired && <span className="file-label" style={{ color: "#93c5fd" }}>R1</span>}
                      <span style={{ overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}
                            title={selectedSample.r1_name}>{selectedSample.r1_name}</span>
                    </div>
                    {selectedSample.r2_name && (
                      <div className="sample-file-row">
                        <span className="file-label" style={{ color: "#93c5fd" }}>R2</span>
                        <span style={{ overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}
                              title={selectedSample.r2_name}>{selectedSample.r2_name}</span>
                      </div>
                    )}
                  </div>
                  <div style={{ fontSize: 11, color: "#1d4ed8", marginTop: 5 }}>
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
                  {jobStatus === "idle" && "Log"}
                  {jobStatus === "running" && "Running"}
                  {jobStatus === "succeeded" && "Done"}
                  {jobStatus === "failed" && "Failed"}
                </span>
                {jobStatus === "running" && currentStep && (
                  <span className="log-step" title={currentStep}>— {currentStep}</span>
                )}
                {jobId && (
                  <span style={{ fontSize: 11, color: "#475569", marginLeft: "auto", flexShrink: 0 }}>
                    job {jobId.slice(0, 8)}
                  </span>
                )}
              </div>
              <div className="log-body">
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
                results.map((f) => {
                  const base = `./api/jobs/${jobId}/file?path=${encodeURIComponent(f.name)}`;
                  return (
                    <div key={f.name} className="result-file">
                      <span className="result-icon">{fileIcon(f.name)}</span>
                      {f.openable ? (
                        <a className="result-name result-link" href={`${base}&inline=1`}
                           target="_blank" rel="noopener noreferrer" title={`Open ${f.name}`}>
                          {f.name}
                        </a>
                      ) : (
                        <a className="result-name result-link" href={`${base}&inline=0`}
                           title={`Download ${f.name}`}>
                          {f.name}
                        </a>
                      )}
                      <span className="result-size">{fmtSize(f.size)}</span>
                      <a className="result-download" href={`${base}&inline=0`}
                         title={`Download ${f.name}`}>⬇</a>
                    </div>
                  );
                })
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
