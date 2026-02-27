import { useEffect, useMemo, useRef, useState } from "react";
import "./styles.css";

const API_BASE =
  (typeof window !== "undefined" && window?.kraken?.apiBase) ||
  import.meta.env.VITE_API_URL ||
  "http://localhost:8000";
const APP_VERSION = "0.2.0";

export default function App() {
  /* ── Section visibility ────────────────────────────────────── */
  const [showSetup, setShowSetup] = useState(false);
  const [showProjects, setShowProjects] = useState(true);
  const [showRun, setShowRun] = useState(true);
  const [showLogs, setShowLogs] = useState(true);

  /* ── Settings ──────────────────────────────────────────────── */
  const [settings, setSettings] = useState({
    kraken_repo_path: "",
    projects_root: "",
    conda_env_path: "",
    default_preset: ""
  });
  const [settingsReady, setSettingsReady] = useState(false);
  const [settingsMsg, setSettingsMsg] = useState("");
  const [settingsChecks, setSettingsChecks] = useState({});
  const [toolChecks, setToolChecks] = useState({});

  /* ── Projects ──────────────────────────────────────────────── */
  const [projects, setProjects] = useState([]);
  const [selectedProject, setSelectedProject] = useState("");
  const [newProjectName, setNewProjectName] = useState("");

  /* ── Presets & overrides ───────────────────────────────────── */
  const [presets, setPresets] = useState([]);
  const [preset, setPreset] = useState("");
  const [overrides, setOverrides] = useState({
    taxon: "",
    kraken_db: "",
    blast_db: "",
    logo: ""
  });

  /* ── SRA download ──────────────────────────────────────────── */
  const [sraAccessions, setSraAccessions] = useState("");

  /* ── Job tracking ──────────────────────────────────────────── */
  const [jobId, setJobId] = useState("");
  const [jobStatus, setJobStatus] = useState("idle");
  const [logText, setLogText] = useState("");
  const [runs, setRuns] = useState([]);
  const [outputs, setOutputs] = useState([]);
  const [selectedRunId, setSelectedRunId] = useState("");

  const logRef = useRef(null);
  const canPickPath = Boolean(window?.kraken?.selectPath);

  /* ── Derived: selected preset details ──────────────────────── */
  const selectedPreset = useMemo(
    () => presets.find((p) => p.name === preset),
    [presets, preset]
  );

  /* ── Derived: effective DB values (preset + override) ──────── */
  const effectiveTaxon = overrides.taxon.trim() || selectedPreset?.taxon || "";
  const effectiveKrakenDb = overrides.kraken_db.trim() || selectedPreset?.kraken_db || "";
  const effectiveBlastDb = overrides.blast_db.trim() || selectedPreset?.blast_db || "";

  /* ── Derived: projects list — last 5, newest first ─────────── */
  const visibleProjects = useMemo(() => {
    const sorted = [...projects].sort((a, b) => {
      // Sort by name descending (newest projects are usually later)
      // If there's a created_at field, prefer it
      if (a.created_at && b.created_at) return b.created_at.localeCompare(a.created_at);
      return b.name.localeCompare(a.name);
    });
    return sorted;
  }, [projects]);

  /* ── Lifecycle ─────────────────────────────────────────────── */
  useEffect(() => {
    refreshAll();
  }, []);

  useEffect(() => {
    if (!jobId) return;
    const source = new EventSource(`${API_BASE}/api/jobs/${jobId}/events`);
    source.onmessage = (event) => {
      if (event.data.startsWith("[job:")) {
        const status = event.data.replace("[job:", "").replace("]", "");
        setJobStatus(status);
        source.close();
        refreshRuns();
        refreshProjects(selectedProject);
        return;
      }
      setLogText((prev) => (prev ? `${prev}\n${event.data}` : event.data));
    };
    source.onerror = () => {
      source.close();
    };
    return () => source.close();
  }, [jobId]);

  useEffect(() => {
    if (logRef.current) {
      logRef.current.scrollTop = logRef.current.scrollHeight;
    }
  }, [logText]);

  /* ── Auto-refresh project stats while job is running ─────── */
  useEffect(() => {
    if (jobStatus !== "running" || !selectedProject) return;
    const interval = setInterval(() => {
      refreshProjects(selectedProject);
    }, 10000);
    return () => clearInterval(interval);
  }, [jobStatus, selectedProject]);

  /* ── API helpers ───────────────────────────────────────────── */
  async function refreshAll() {
    await loadConfig();
    await loadPresets();
    await refreshProjects();
  }

  async function loadConfig() {
    const res = await fetch(`${API_BASE}/api/config`);
    const cfg = await res.json();
    setSettings(cfg);
    setPreset(cfg.default_preset || "");
    const preflight = await fetch(`${API_BASE}/api/preflight`);
    const pf = await preflight.json();
    setSettingsReady(pf.ok);
    setSettingsMsg(pf.issues?.join(" | ") || "");
    setSettingsChecks(pf.checks || {});
    setToolChecks(pf.tools || {});
  }

  async function saveSettings() {
    const res = await fetch(`${API_BASE}/api/config`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(settings)
    });
    if (res.ok) {
      await loadConfig();
    }
  }

  async function loadPresets() {
    try {
      const res = await fetch(`${API_BASE}/api/presets`);
      if (!res.ok) return;
      const data = await res.json();
      setPresets(data);
      if (!preset && data.length) {
        setPreset(data[0].name);
      }
    } catch {
      setPresets([]);
    }
  }

  async function refreshProjects(nextSelected = selectedProject) {
    const res = await fetch(`${API_BASE}/api/projects`);
    if (!res.ok) return;
    const data = await res.json();
    setProjects(data);
    if (nextSelected) {
      setSelectedProject(nextSelected);
      await refreshRuns(nextSelected);
    }
  }

  async function createProject() {
    const name = newProjectName.trim();
    if (!name) return;
    const res = await fetch(`${API_BASE}/api/projects`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name })
    });
    if (res.ok) {
      setNewProjectName("");
      await refreshProjects(name);
    }
  }

  async function deleteProject(name) {
    if (!window.confirm(`Archive project "${name}"?\n\nThe project folder will be moved to projects_archive/ and can be restored later.`)) return;
    const res = await fetch(`${API_BASE}/api/projects/${name}`, { method: "DELETE" });
    if (res.ok) {
      await refreshProjects("");
      setSelectedProject("");
      setRuns([]);
      setOutputs([]);
    }
  }

  async function linkLocalFolder() {
    if (!selectedProject) return;
    let picked = "";
    if (canPickPath) {
      picked = await window.kraken.selectPath({
        kind: "folder",
        title: "Select FASTQ folder",
        defaultPath: settings.projects_root || undefined
      });
    } else {
      picked = window.prompt("Enter FASTQ folder path:") || "";
    }
    if (!picked) return;
    await fetch(`${API_BASE}/api/projects/${selectedProject}/link-local`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path: picked })
    });
    await refreshProjects(selectedProject);
  }

  async function downloadSRA() {
    if (!selectedProject || !sraAccessions.trim()) return;
    setLogText("");
    setJobStatus("running");
    const res = await fetch(`${API_BASE}/api/projects/${selectedProject}/sra-download`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ accessions: sraAccessions.trim() })
    });
    if (res.ok) {
      const data = await res.json();
      setJobId(data.job_id || "");
      setSraAccessions("");
    } else {
      setJobStatus("failed");
      const msg = await res.json();
      setLogText(msg.detail || "SRA download failed");
    }
  }

  async function runPreset() {
    if (!selectedProject || !preset) return;
    setLogText("");
    setJobStatus("running");
    const cleanOverrides = Object.fromEntries(
      Object.entries(overrides)
        .map(([key, value]) => {
          const trimmed = (value || "").trim();
          if (!trimmed) return [key, ""];
          if (trimmed.includes("=")) {
            const [maybeKey, ...rest] = trimmed.split("=");
            if (maybeKey === key) {
              return [key, rest.join("=").trim()];
            }
          }
          return [key, trimmed];
        })
        .filter(([, v]) => v)
    );
    const res = await fetch(`${API_BASE}/api/projects/${selectedProject}/run`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ preset, overrides: cleanOverrides })
    });
    if (res.ok) {
      const data = await res.json();
      setJobId(data.job_id || "");
      await refreshRuns(selectedProject);
    } else {
      setJobStatus("failed");
      const msg = await res.json();
      setLogText(msg.detail || "Failed to start run");
    }
  }

  async function refreshRuns(projectName = selectedProject) {
    if (!projectName) return;
    const res = await fetch(`${API_BASE}/api/projects/${projectName}/runs`);
    if (!res.ok) return;
    const data = await res.json();
    setRuns(data);
    if (data.length) {
      await loadOutputs(projectName, data[0].run_id);
    }
  }

  async function deleteRun(runId) {
    if (!selectedProject || !runId) return;
    if (!window.confirm(`Delete run ${runId}?`)) return;
    const res = await fetch(`${API_BASE}/api/projects/${selectedProject}/runs/${runId}`, {
      method: "DELETE"
    });
    if (res.ok) {
      await refreshRuns(selectedProject);
      setOutputs([]);
    }
  }

  async function loadOutputs(projectName, runId) {
    if (!projectName || !runId) return;
    const res = await fetch(`${API_BASE}/api/projects/${projectName}/runs/${runId}/outputs`);
    if (!res.ok) return;
    const data = await res.json();
    setOutputs(data);
    setSelectedRunId(runId);
  }

  async function openOutput(runId, relativePath) {
    if (!selectedProject || !runId || !relativePath) return;
    await fetch(`${API_BASE}/api/projects/${selectedProject}/runs/${runId}/open`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path: relativePath })
    });
  }

  const selected = useMemo(
    () => projects.find((p) => p.name === selectedProject),
    [projects, selectedProject]
  );

  /* ── File icon helper ──────────────────────────────────────── */
  function fileIcon(name) {
    if (name.endsWith(".pdf")) return "📄";
    if (name.endsWith(".xlsx") || name.endsWith(".xls")) return "📊";
    if (name.endsWith(".html")) return "🌐";
    if (name.endsWith(".fasta")) return "🧬";
    if (name.endsWith(".txt")) return "📝";
    if (name.endsWith(".png") || name.endsWith(".jpg")) return "🖼️";
    return "📁";
  }

  /* ================================================================== */
  /* RENDER                                                              */
  /* ================================================================== */
  return (
    <div className="app">
      {/* ── Header ─────────────────────────────────────────────── */}
      <header className="app-header">
        <div className="app-brand">
          <img className="app-logo" src="/kraken_icon.png" alt="Kraken icon" />
          <div>
            <h1>
              Kraken ID Parse <span className="version-tag">v{APP_VERSION}</span>
            </h1>
            <p>Classify and isolate reads for species-level assembly and analyses</p>
          </div>
        </div>
        <div className="status-pill">
          <span className="dot" data-state={jobStatus}></span>
          <span>{jobStatus}</span>
        </div>
      </header>

      <main className="layout">
        {/* ── Alert banner ───────────────────────────────────────── */}
        {!settingsReady ? (
          <div className="panel alert-banner">
            <strong>Setup required:</strong> Configure paths in Settings, then Save.
            {settingsMsg ? <div className="note warning" style={{ marginTop: "8px" }}>{settingsMsg}</div> : null}
          </div>
        ) : null}

        {/* ── Status Strip ───────────────────────────────────────── */}
        <section className="status-strip">
          <div className="status-item">
            <span className="status-label">Project</span>
            <span className="status-value">{selected?.name || "—"}</span>
          </div>
          <div className="status-item">
            <span className="status-label">FASTQ</span>
            <span className="status-value">{selected?.fastq_count ?? 0}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Runs</span>
            <span className="status-value">{selected?.run_count ?? 0}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Reports</span>
            <span className="status-value">{selected?.report_count ?? 0}</span>
          </div>
          <div className="status-item">
            <span className="status-label">Job</span>
            <span className="status-value">
              {jobStatus === "running" ? <><span className="pulse-dot" />{jobStatus}</> : jobStatus}
            </span>
          </div>
        </section>

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Settings (collapsed by default)                */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Settings</h2>
          <button className="ghost" onClick={() => setShowSetup(!showSetup)}>
            {showSetup ? "Hide" : "Show"}
          </button>
        </div>
        {showSetup ? (
          <div className="row-grid" style={{ gridTemplateColumns: "1fr" }}>
            <section className="panel">
              <div className="input-columns" style={{ gridTemplateColumns: "1fr" }}>
                <div className="input-column">
                  {/* Kraken repo path */}
                  <div className="settings-row">
                    <label className="label">Kraken repo path</label>
                    <input
                      placeholder="/path/to/kraken_id_parse_gui"
                      value={settings.kraken_repo_path}
                      onChange={(e) => setSettings({ ...settings, kraken_repo_path: e.target.value })}
                    />
                    <span style={{ display: "inline-flex", alignItems: "center", gap: "6px" }}>
                      {canPickPath ? (
                        <button className="ghost action" onClick={() =>
                          window.kraken.selectPath({ kind: "folder", title: "Select Kraken repo", defaultPath: settings.kraken_repo_path || undefined })
                            .then((v) => v && setSettings({ ...settings, kraken_repo_path: v }))
                        }>Choose</button>
                      ) : null}
                      {settingsChecks.kraken_repo_path === true
                        ? <span style={{ color: "var(--success)", fontWeight: 700 }}>✓</span>
                        : <span style={{ color: "var(--danger)", fontWeight: 700 }}>✕</span>}
                    </span>
                  </div>
                  {/* Projects root */}
                  <div className="settings-row">
                    <label className="label">Projects root</label>
                    <input
                      placeholder="/path/to/projects"
                      value={settings.projects_root}
                      onChange={(e) => setSettings({ ...settings, projects_root: e.target.value })}
                    />
                    <span style={{ display: "inline-flex", alignItems: "center", gap: "6px" }}>
                      {canPickPath ? (
                        <button className="ghost action" onClick={() =>
                          window.kraken.selectPath({ kind: "folder", title: "Select projects root", defaultPath: settings.projects_root || undefined })
                            .then((v) => v && setSettings({ ...settings, projects_root: v }))
                        }>Choose</button>
                      ) : null}
                      {settingsChecks.projects_root === true
                        ? <span style={{ color: "var(--success)", fontWeight: 700 }}>✓</span>
                        : <span style={{ color: "var(--danger)", fontWeight: 700 }}>✕</span>}
                    </span>
                  </div>
                  {/* Conda env */}
                  <div className="settings-row">
                    <label className="label">Conda env path</label>
                    <input
                      placeholder="/path/to/conda/env"
                      value={settings.conda_env_path}
                      onChange={(e) => setSettings({ ...settings, conda_env_path: e.target.value })}
                    />
                    <span style={{ display: "inline-flex", alignItems: "center", gap: "6px" }}>
                      {canPickPath ? (
                        <button className="ghost action" onClick={() =>
                          window.kraken.selectPath({ kind: "folder", title: "Select conda env", defaultPath: settings.conda_env_path || undefined })
                            .then((v) => v && setSettings({ ...settings, conda_env_path: v }))
                        }>Choose</button>
                      ) : null}
                      {settingsChecks.conda_env_path === true
                        ? <span style={{ color: "var(--success)", fontWeight: 700 }}>✓</span>
                        : <span style={{ color: "var(--danger)", fontWeight: 700 }}>✕</span>}
                    </span>
                  </div>
                </div>
              </div>
              <div style={{ display: "flex", justifyContent: "flex-end", marginTop: "12px" }}>
                <button onClick={saveSettings}>Save</button>
              </div>
              {/* Pipeline tool badges */}
              {Object.keys(toolChecks).length > 0 && (
                <div style={{ marginTop: "16px", borderTop: "1px dashed var(--border)", paddingTop: "12px" }}>
                  <label className="label" style={{ marginBottom: "8px", display: "block" }}>Pipeline Tools</label>
                  <div style={{ display: "flex", flexWrap: "wrap", gap: "6px" }}>
                    {Object.entries(toolChecks).map(([name, info]) => (
                      <span
                        key={name}
                        title={info.found ? `${info.path}\n${info.version || ""}` : `${name} not found — install in conda env`}
                        className={`badge ${info.found ? "complete" : "error"}`}
                        style={{ cursor: "help" }}
                      >
                        {info.found ? "✓" : "✕"} {name}
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </section>
          </div>
        ) : null}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Projects & Inputs                              */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Projects &amp; Inputs</h2>
          <button className="ghost" onClick={() => setShowProjects(!showProjects)}>
            {showProjects ? "Hide" : "Show"}
          </button>
        </div>
        {showProjects ? (
          <div className="row-grid">
            {/* LEFT: Projects list */}
            <section className="panel">
              <h2>Projects</h2>
              <div className="row">
                <input
                  placeholder="New project name"
                  value={newProjectName}
                  onChange={(e) => setNewProjectName(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && createProject()}
                />
                <button onClick={createProject}>Create</button>
              </div>
              <div className="list" style={{ maxHeight: "260px", minHeight: 0 }}>
                {visibleProjects.map((p) => (
                  <div
                    key={p.name}
                    className={`list-item ${p.name === selectedProject ? "active" : ""}`}
                    role="button"
                    tabIndex={0}
                    onClick={() => {
                      setSelectedProject(p.name);
                      refreshRuns(p.name);
                    }}
                  >
                    <div className="list-details">
                      <div className="list-title">{p.name}</div>
                      <div className="list-meta">
                        FASTQ: {p.fastq_count} &nbsp;|&nbsp; Runs: {p.run_count} &nbsp;|&nbsp; Reports: {p.report_count}
                      </div>
                    </div>
                    <div className="list-actions">
                      <button
                        className="ghost-btn danger"
                        title="Archive project"
                        onClick={(e) => {
                          e.stopPropagation();
                          deleteProject(p.name);
                        }}
                      >
                        Archive
                      </button>
                    </div>
                  </div>
                ))}
                {!projects.length && <div className="note">No projects yet.</div>}
              </div>
            </section>

            {/* RIGHT: Inputs — FASTQ linking + SRA download */}
            <section className="panel">
              <h2>Inputs</h2>
              <button onClick={linkLocalFolder} disabled={!selectedProject || !settingsReady}>
                Link FASTQ Folder
              </button>
              {selected?.last_input_path ? (
                <div className="note">Linked: <strong>{selected.last_input_path}</strong></div>
              ) : null}

              <div style={{ marginTop: "16px", borderTop: "1px dashed var(--border)", paddingTop: "12px" }}>
                <label className="label">SRA Download</label>
                <div className="row">
                  <input
                    placeholder="SRR28623786, SRR9598511, ..."
                    value={sraAccessions}
                    onChange={(e) => setSraAccessions(e.target.value)}
                    onKeyDown={(e) => e.key === "Enter" && downloadSRA()}
                    disabled={!selectedProject || !settingsReady}
                    style={{ flex: 1 }}
                  />
                  <button
                    onClick={downloadSRA}
                    disabled={!selectedProject || !settingsReady || !sraAccessions.trim()}
                  >
                    Download
                  </button>
                </div>
                <div className="note">Enter SRA accession(s) — comma, space, or newline separated.</div>
              </div>
              {!selectedProject || !settingsReady ? (
                <div className="note warning" style={{ marginTop: "12px" }}>
                  Select a project and complete Settings first.
                </div>
              ) : null}
            </section>
          </div>
        ) : null}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Run Kraken — LEFT/RIGHT SPLIT                  */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Run Kraken</h2>
          <button className="ghost" onClick={() => setShowRun(!showRun)}>
            {showRun ? "Hide" : "Show"}
          </button>
        </div>
        {showRun ? (
          <div className="row-grid row-grid-split">
            {/* LEFT: Preset selection + Run controls */}
            <section className="panel">
              <h2>Configure &amp; Run</h2>

              {/* Preset selector */}
              <label className="label">Preset</label>
              <select value={preset} onChange={(e) => setPreset(e.target.value)}>
                {presets.map((p) => (
                  <option key={p.name} value={p.name}>{p.name}</option>
                ))}
              </select>

              {/* Dynamic DB info from selected preset */}
              {selectedPreset ? (
                <div className="selection-box" style={{ marginTop: "10px" }}>
                  <div><span className="muted">Taxon:</span> <strong>{effectiveTaxon || "—"}</strong></div>
                  <div><span className="muted">Kraken DB:</span> <strong>{effectiveKrakenDb ? effectiveKrakenDb.split("/").pop() : "—"}</strong></div>
                  <div><span className="muted">BLAST DB:</span> <strong>{effectiveBlastDb ? effectiveBlastDb.split("/").pop() : "auto-resolved"}</strong></div>
                </div>
              ) : null}

              {/* Overrides */}
              <details className="step2-options-panel" style={{ marginTop: "12px" }}>
                <summary style={{ cursor: "pointer", fontWeight: 600, fontSize: "13px" }}>
                  Overrides
                  {Object.values(overrides).filter((v) => v.trim()).length > 0
                    ? ` (${Object.values(overrides).filter((v) => v.trim()).length} active)`
                    : " (optional)"}
                </summary>
                <div style={{ display: "grid", gap: "8px", marginTop: "8px" }}>
                  <div>
                    <label className="label">Taxon</label>
                    <input
                      placeholder={selectedPreset?.taxon || "e.g. Orbivirus"}
                      value={overrides.taxon}
                      onChange={(e) => setOverrides({ ...overrides, taxon: e.target.value })}
                    />
                  </div>
                  <div>
                    <label className="label">Kraken DB</label>
                    <div className="row" style={{ marginBottom: 0 }}>
                      <input
                        placeholder={selectedPreset?.kraken_db ? selectedPreset.kraken_db.split("/").pop() : "auto"}
                        value={overrides.kraken_db}
                        onChange={(e) => setOverrides({ ...overrides, kraken_db: e.target.value })}
                        style={{ flex: 1 }}
                      />
                      {canPickPath ? (
                        <button className="ghost action" style={{ whiteSpace: "nowrap" }} onClick={() =>
                          window.kraken.selectPath({
                            kind: "folder",
                            title: "Select Kraken DB folder",
                            defaultPath: selectedPreset?.kraken_db
                              ? selectedPreset.kraken_db.split("/").slice(0, -1).join("/")
                              : settings.projects_root || undefined
                          }).then((v) => v && setOverrides({ ...overrides, kraken_db: v }))
                        }>Choose</button>
                      ) : null}
                    </div>
                  </div>
                  <div>
                    <label className="label">BLAST DB</label>
                    <div className="row" style={{ marginBottom: 0 }}>
                      <input
                        placeholder={selectedPreset?.blast_db ? selectedPreset.blast_db.split("/").pop() : "auto-resolved from taxon"}
                        value={overrides.blast_db}
                        onChange={(e) => setOverrides({ ...overrides, blast_db: e.target.value })}
                        style={{ flex: 1 }}
                      />
                      {canPickPath ? (
                        <button className="ghost action" style={{ whiteSpace: "nowrap" }} onClick={() =>
                          window.kraken.selectPath({
                            kind: "folder",
                            title: "Select BLAST DB folder",
                            defaultPath: selectedPreset?.blast_db
                              ? selectedPreset.blast_db.split("/").slice(0, -1).join("/")
                              : settings.projects_root || undefined
                          }).then((v) => v && setOverrides({ ...overrides, blast_db: v }))
                        }>Choose</button>
                      ) : null}
                    </div>
                  </div>
                  <div>
                    <label className="label">Logo</label>
                    <div className="row" style={{ marginBottom: 0 }}>
                      <input
                        placeholder="path to logo PNG"
                        value={overrides.logo}
                        onChange={(e) => setOverrides({ ...overrides, logo: e.target.value })}
                        style={{ flex: 1 }}
                      />
                      {canPickPath ? (
                        <button className="ghost action" style={{ whiteSpace: "nowrap" }} onClick={() =>
                          window.kraken.selectPath({
                            kind: "file",
                            title: "Select logo image",
                            defaultPath: settings.kraken_repo_path || undefined
                          }).then((v) => v && setOverrides({ ...overrides, logo: v }))
                        }>Choose</button>
                      ) : null}
                    </div>
                  </div>
                </div>
              </details>

              <div className="step1-actions" style={{ marginTop: "12px" }}>
                <button
                  onClick={runPreset}
                  disabled={!selectedProject || !settingsReady || !preset}
                  style={{ width: "100%" }}
                >
                  {jobStatus === "running" ? "Running…" : "Run Pipeline"}
                </button>
              </div>

              {/* Recent runs list */}
              <div style={{ marginTop: "16px", borderTop: "1px dashed var(--border)", paddingTop: "12px" }}>
                <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
                  <label className="label" style={{ marginBottom: 0 }}>Recent Runs</label>
                  <button className="ghost action" style={{ fontSize: "11px", padding: "4px 8px" }}
                    onClick={() => refreshRuns(selectedProject)} disabled={!selectedProject}>
                    Refresh
                  </button>
                </div>
                <div className="list" style={{ maxHeight: "200px", marginTop: "8px" }}>
                  {runs.length ? (
                    runs.map((r) => (
                      <div
                        key={r.run_id}
                        className={`list-item ${r.run_id === selectedRunId ? "active" : ""}`}
                        onClick={() => loadOutputs(selectedProject, r.run_id)}
                      >
                        <div className="list-details">
                          <div className="list-title">{r.run_id}</div>
                          <div className="list-meta">
                            <span className={`badge ${r.status}`}>{r.status}</span>
                            &nbsp; {r.outputs} output{r.outputs !== 1 ? "s" : ""}
                          </div>
                        </div>
                        <div className="list-actions">
                          <button className="ghost-btn danger" onClick={(e) => { e.stopPropagation(); deleteRun(r.run_id); }}>✕</button>
                        </div>
                      </div>
                    ))
                  ) : (
                    <div className="note">No runs yet.</div>
                  )}
                </div>
              </div>
            </section>

            {/* RIGHT: Outputs + Log */}
            <section className="panel results-panel" style={{ minHeight: "360px" }}>
              <div className="panel-header">
                <h2>Outputs</h2>
                <div className="panel-actions">
                  <button
                    className="ghost action"
                    onClick={() => openOutput(selectedRunId || runs[0]?.run_id, ".")}
                    disabled={!selectedRunId && !(runs[0]?.run_id)}
                  >
                    Open Folder
                  </button>
                </div>
              </div>
              <div className="results-list" style={{ maxHeight: "480px", minHeight: "120px" }}>
                {outputs.length ? (
                  outputs.map((o) => (
                    <div key={o.path} className="results-item">
                      <div className="results-main">
                        <div className="results-name">{fileIcon(o.name)} {o.name}</div>
                        <div className="results-path">{o.relative}</div>
                      </div>
                      <div className="results-actions">
                        <button
                          className="ghost-btn"
                          onClick={() => openOutput(selectedRunId || runs[0]?.run_id, o.relative)}
                        >
                          Open
                        </button>
                      </div>
                    </div>
                  ))
                ) : (
                  <div className="note">
                    {selectedRunId ? "No output files found." : "Select a run to view outputs."}
                  </div>
                )}
              </div>

            </section>
          </div>
        ) : null}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Logs (optional expanded view)                  */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Pipeline Log</h2>
          <button className="ghost" onClick={() => setShowLogs(!showLogs)}>
            {showLogs ? "Hide" : "Show"}
          </button>
        </div>
        {showLogs ? (
          <div className="row-grid" style={{ gridTemplateColumns: "1fr" }}>
            <section className="panel" style={{ minHeight: "unset" }}>
              <div style={{ marginTop: 0 }}>
                <div ref={logRef} style={{ maxHeight: "480px", overflow: "auto" }}>
                  <pre style={{
                    background: "#151b21", color: "#f3f4f6", padding: "12px",
                    borderRadius: "8px", fontSize: "12px", lineHeight: "1.4",
                    whiteSpace: "pre-wrap", wordBreak: "break-word", margin: 0,
                    minHeight: "100px"
                  }}>
                    {logText || "No log yet."}
                  </pre>
                </div>
              </div>
            </section>
          </div>
        ) : null}
      </main>
    </div>
  );
}
