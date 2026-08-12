import { useEffect, useRef, useState } from "react";
import "./App.css";
import ThemeToggle from "./ThemeToggle";
import CitationFooter from "./Citations";
import ResultsPane from "./ResultsPane";
import { useResults } from "./useResults";

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const APP_VERSION = "0.2.6";

// Taxon presets are loaded at runtime from the shared config/taxa.yaml via
// /api/taxa. This fallback is only used if that fetch fails.
const TAXON_PRESETS_FALLBACK = [
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
// The tool-specific columns of the shared Results table. Everything else
// about the pane is identical across the suite.
const RESULT_COLUMNS = [
  { key: "top_taxon", label: "Top taxon" },
  { key: "top_pct", label: "%", align: "right" },
  { key: "runner_up", label: "Runner-up" },
];

export default function App() {
  const [projects, setProjects] = useState([]);
  const [projectsLoading, setProjectsLoading] = useState(true);
  const [newProjectName, setNewProjectName] = useState("");
  const [creatingProject, setCreatingProject] = useState(false);
  // Sample-loading state, keyed by project name so multiple expanded projects
  // don't clobber each other.
  const [activeProject, setActiveProject] = useState(""); // project the Inputs pane targets
  /* Every completed sample for the active project. Refreshed when a run
     finishes rather than polled, matching the rest of the suite. */
  const results = useResults(activeProject);
  const [showResultsPane, setShowResultsPane] = useState(true);
  const [addPath, setAddPath] = useState({});       // proj -> import path string
  const [sraText, setSraText] = useState({});       // proj -> SRA accessions string
  const [addStatus, setAddStatus] = useState({});   // proj -> status message
  const [inputsByProj, setInputsByProj] = useState({}); // proj -> {files,count,total_bytes}
  const uploadProjRef = useRef("");                 // which project the file dialog targets
  const uploadInputRef = useRef(null);
  const [expanded, setExpanded] = useState({});          // project name → bool
  const [samples, setSamples] = useState({});            // project name → [sample]
  const [checkedKeys, setCheckedKeys] = useState({});    // key → {project, ...sample}  (batch selection)
  // Filter for the Projects sample list. Its check-all acts on what this
  // filter leaves visible, never on the whole project.
  const [projSampleFilter, setProjSampleFilter] = useState("");
  const [openResults, setOpenResults] = useState({});    // key → bool (inline results expanded)
  const [sampleResults, setSampleResults] = useState({}); // key → {loading, status, present, files}
  const [vsnpResults, setVsnpResults] = useState({});    // key → {loading, step1_present, files, step2} (cross-tool)
  const [activeRun, setActiveRun] = useState(null);      // {project, sample} currently running
  const [queueInfo, setQueueInfo] = useState({ total: 0, done: 0 }); // batch progress
  const [taxon, setTaxon] = useState("");
  const [taxonPresets, setTaxonPresets] = useState(TAXON_PRESETS_FALLBACK); // from /api/taxa
  const [newTaxon, setNewTaxon] = useState("");          // "add search name" input
  const [addingTaxon, setAddingTaxon] = useState(false);
  const [krakenOnly, setKrakenOnly] = useState(false);   // Kraken2 + Krona only, no read parsing
  const [noBlast, setNoBlast] = useState(false);         // Kraken2 + read parsing only, no assembly/BLAST
  const [krakenDb, setKrakenDb] = useState("");
  const [blastDb, setBlastDb] = useState("nt");
  // Known Kraken2 DBs — remembered (saved_kraken_dbs) + discovered on disk —
  // from /api/kraken-dbs. Powers the quick-switch dropdowns in Settings and
  // the run form.
  const [knownDbs, setKnownDbs] = useState([]);   // [{name, path, size_bytes, missing?}]
  const [newDbPath, setNewDbPath] = useState("");  // "Add database location" input
  const [addDbError, setAddDbError] = useState("");
  const [running, setRunning] = useState(false);
  const [jobId, setJobId] = useState(null);
  const [jobStatus, setJobStatus] = useState("idle"); // idle | running | succeeded | failed
  const [logLines, setLogLines] = useState([]);
  const [settingsDraft, setSettingsDraft] = useState({});
  // Server-side folder browser for picking the projects root.
  const [folderBrowser, setFolderBrowser] = useState({ open: false, path: "", parent: null, entries: [], loading: false, error: "" });
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
    loadKnownDbs();
    fetch("./api/taxa")
      .then((r) => r.json())
      .then((d) => {
        if (Array.isArray(d.taxa) && d.taxa.length) setTaxonPresets(d.taxa);
      })
      .catch(() => {});
    loadProjects();
    // The Results table is where finished work is now read, so it
    // has to reflect the run that just ended.
    results.reload();
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

  // Persist a new taxon search name to the shared config/taxa.yaml and select it.
  async function addTaxon() {
    const name = newTaxon.trim();
    if (!name || addingTaxon) return;
    setAddingTaxon(true);
    try {
      const res = await fetch("./api/taxa", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name }),
      });
      if (res.ok) {
        const d = await res.json();
        if (Array.isArray(d.taxa)) setTaxonPresets(d.taxa);
        setTaxon(name);
        setNewTaxon("");
      }
    } catch (_) {
      /* leave the input as-is so the user can retry */
    } finally {
      setAddingTaxon(false);
    }
  }

  // Create a new project (personal root by default). Projects share the same
  // on-disk layout vSNP uses, so a project made here also shows up in vSNP GUI.
  async function createProject() {
    const name = newProjectName.trim();
    if (!name || creatingProject) return;
    setCreatingProject(true);
    try {
      const res = await fetch("./api/projects", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name }),
      });
      if (!res.ok) {
        const detail = await res.json().catch(() => ({}));
        window.alert(`Could not create project: ${detail.detail || res.status}`);
        return;
      }
      const created = await res.json().catch(() => ({}));
      setNewProjectName("");
      loadProjects();
      // Auto-expand the new project and open the Add-samples panel so the user
      // can import / drop / download reads right away. Fetch its (empty)
      // sample + input lists so the row doesn't hang on "Loading samples…".
      if (created.name) {
        const n = created.name;
        setExpanded((e) => ({ ...e, [n]: true }));
        setActiveProject(n);
        await Promise.all([fetchSamples(n), loadInputs(n)]);
      }
    } finally {
      setCreatingProject(false);
    }
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

  // Keep the Inputs pane pointed at a real project: default to the first one
  // and recover if the current target disappears (e.g. after a refresh). This
  // means the import/upload/download controls are always available as long as
  // at least one project exists — no need to hunt for "the active project".
  useEffect(() => {
    if (!projects.length) {
      if (activeProject) setActiveProject("");
      return;
    }
    if (!activeProject || !projects.find((p) => p.name === activeProject)) {
      const first = projects[0].name;
      setActiveProject(first);
      if (inputsByProj[first] === undefined) loadInputs(first);
    }
  }, [projects]);

  function fetchSamples(name) {
    return fetch(`./api/projects/${encodeURIComponent(name)}/samples`)
      .then((r) => r.json())
      .then((data) => setSamples((s) => ({ ...s, [name]: data })))
      .catch(() => setSamples((s) => ({ ...s, [name]: [] })));
  }

  function toggleProject(name) {
    const isExpanded = expanded[name];
    setExpanded((e) => ({ ...e, [name]: !isExpanded }));
    // Target this project for the Inputs pane whenever the user opens it.
    setActiveProject(name);
    if (!isExpanded) {
      if (!samples[name]) fetchSamples(name);
      loadInputs(name);
    }
  }

  // Make a project the Inputs-pane target without toggling its expansion.
  function selectProject(name) {
    setActiveProject(name);
    if (inputsByProj[name] === undefined) loadInputs(name);
  }

  // ---- Sample loading (import / upload / SRA) -------------------------------
  function loadInputs(name) {
    return fetch(`./api/projects/${encodeURIComponent(name)}/inputs`)
      .then((r) => r.json())
      .then((data) => setInputsByProj((m) => ({ ...m, [name]: data })))
      .catch(() => setInputsByProj((m) => ({ ...m, [name]: { files: [], count: 0, total_bytes: 0 } })));
  }

  const setStat = (name, msg) => setAddStatus((m) => ({ ...m, [name]: msg }));

  async function refreshAfterLoad(name) {
    await Promise.all([fetchSamples(name), loadInputs(name)]);
    loadProjects();
  }

  async function linkLocal(name) {
    const path = (addPath[name] || "").trim();
    if (!path) return;
    setStat(name, "Linking…");
    try {
      const res = await fetch(`./api/projects/${encodeURIComponent(name)}/link-local`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ path }),
      });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) { setStat(name, `Import failed: ${data.detail || res.status}`); return; }
      setStat(name, `Linked ${data.linked} file${data.linked === 1 ? "" : "s"}.`);
      setAddPath((m) => ({ ...m, [name]: "" }));
      await refreshAfterLoad(name);
    } catch (e) {
      setStat(name, `Import failed: ${e.message}`);
    }
  }

  function pickFiles(name) {
    uploadProjRef.current = name;
    uploadInputRef.current?.click();
  }

  async function uploadFiles(name, fileList) {
    const files = Array.from(fileList || []).filter((f) => f.name.endsWith(".fastq.gz"));
    if (!name || !files.length) return;
    const fd = new FormData();
    files.forEach((f) => fd.append("files", f));
    setStat(name, `Uploading ${files.length} file${files.length === 1 ? "" : "s"}…`);
    try {
      const res = await fetch(`./api/projects/${encodeURIComponent(name)}/upload`, { method: "POST", body: fd });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) { setStat(name, `Upload failed: ${data.detail || res.status}`); return; }
      setStat(name, `Uploaded ${data.uploaded} file${data.uploaded === 1 ? "" : "s"}.`);
      await refreshAfterLoad(name);
    } catch (e) {
      setStat(name, `Upload failed: ${e.message}`);
    }
  }

  function parseAccessions(text) {
    return (text || "").split(/[\s,]+/).map((s) => s.trim()).filter(Boolean);
  }

  async function sraDownload(name) {
    const accessions = parseAccessions(sraText[name]);
    if (!accessions.length) return;
    setStat(name, `Resolving ${accessions.length} accession${accessions.length === 1 ? "" : "s"}…`);
    setShowLogs(true);
    try {
      const res = await fetch(`./api/projects/${encodeURIComponent(name)}/sra/download`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ accessions }),
      });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) { setStat(name, `Download failed: ${data.detail || res.status}`); return; }
      setStat(name, "Downloading… progress shows in the Pipeline Log below.");
      setSraText((m) => ({ ...m, [name]: "" }));
      // Stream the download log into the shared log panel; refresh when done.
      setJobId(data.job_id);
      setJobStatus("running");
      setLogLines([]);
      streamLogUntilDone(data.job_id, null, () => {
        setStat(name, "Download finished — see samples below.");
        refreshAfterLoad(name);
      });
    } catch (e) {
      setStat(name, `Download failed: ${e.message}`);
    }
  }

  async function deleteInput(name, filename) {
    if (!window.confirm(`Remove ${filename} from this project's download/ folder?`)) return;
    try {
      await fetch(`./api/projects/${encodeURIComponent(name)}/inputs/${encodeURIComponent(filename)}`, { method: "DELETE" });
      await refreshAfterLoad(name);
    } catch (e) {
      setStat(name, `Delete failed: ${e.message}`);
    }
  }

  // --- Sample selection / results (per-sample, decoupled from a single job) ---
  const sampleKey = (project, s) => `${project}::${s.sample}`;
  const isActive = (project, s) =>
    activeRun && activeRun.project === project && activeRun.sample === s.sample;

  /* Samples currently VISIBLE in a project, i.e. after the filter box. The
     check-all below must use this and nothing else: a "select all" that also
     queues samples the user cannot see is how people accidentally run hundreds
     of samples instead of the handful they filtered to. */
  function visibleSamples(project) {
    const q = projSampleFilter.trim().toLowerCase();
    const list = samples[project] || [];
    return q ? list.filter((s) => String(s.sample || "").toLowerCase().includes(q)) : list;
  }

  function checkAllState(project) {
    const vis = visibleSamples(project);
    const on = vis.filter((s) => checkedKeys[sampleKey(project, s)]).length;
    return { total: vis.length, on, checked: vis.length > 0 && on === vis.length,
             indeterminate: on > 0 && on < vis.length };
  }

  function toggleCheckAllVisible(project, checked) {
    const vis = visibleSamples(project);
    setCheckedKeys((m) => {
      const next = { ...m };
      vis.forEach((s) => {
        const k = sampleKey(project, s);
        if (checked) next[k] = { project, ...s };
        else delete next[k];
      });
      return next;
    });
  }

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
          no_blast: noBlast,
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

  // --- Server-side folder browser (projects root + add-database) ----------
  function browseDirs(path) {
    setFolderBrowser((s) => ({ ...s, loading: true, error: "" }));
    fetch(`./api/browse-dirs?path=${encodeURIComponent(path || "")}`)
      .then((r) => (r.ok ? r.json() : r.json().then((e) => { throw new Error(e.detail || "Cannot open folder"); })))
      .then((d) => setFolderBrowser((s) => ({ ...s, path: d.path, parent: d.parent, entries: d.entries, loading: false })))
      .catch((err) => setFolderBrowser((s) => ({ ...s, loading: false, error: err.message })));
  }
  function openFolderBrowser(target = "projects_root") {
    setFolderBrowser({ open: true, target, path: "", parent: null, entries: [], loading: true, error: "" });
    browseDirs(target === "add_db" ? (newDbPath || "") : (settingsDraft.projects_root || ""));
  }
  function chooseFolder() {
    if (folderBrowser.target === "add_db") {
      setNewDbPath(folderBrowser.path);
    } else {
      setSettingsDraft((d) => ({ ...d, projects_root: folderBrowser.path }));
    }
    setFolderBrowser((s) => ({ ...s, open: false }));
  }

  function loadKnownDbs() {
    fetch("./api/kraken-dbs")
      .then((r) => (r.ok ? r.json() : { databases: [] }))
      .then((data) => setKnownDbs(data.databases || []))
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
        saved_project_roots: settingsDraft.saved_project_roots,
        saved_kraken_dbs: settingsDraft.saved_kraken_dbs,
      }),
    })
      .then((r) => r.json())
      .then(() => {
        setKrakenDb(settingsDraft.kraken_db || "");
        setBlastDb(settingsDraft.blast_db || "nt");
        loadKnownDbs();
        loadProjects();
      })
      .catch(() => {});
  }

  // Explicitly add one database location to the saved list. The backend
  // validates it is a real Kraken2 DB (hash/opts/taxo.k2d) before saving, so
  // a typo is an error here instead of a failed run later.
  function addDbLocation() {
    const p = (newDbPath || "").trim();
    if (!p) return;
    setAddDbError("");
    fetch("./api/kraken-dbs/add", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path: p }),
    })
      .then((r) => (r.ok ? r.json() : r.json().then((e) => { throw new Error(e.detail || "Could not add database"); })))
      .then((d) => {
        setNewDbPath("");
        setSettingsDraft((prev) => ({
          ...prev,
          saved_kraken_dbs: d.saved,
          // First database added with none active: make it the active one.
          kraken_db: prev.kraken_db || p,
        }));
        if (!(settingsDraft.kraken_db || "").trim()) setKrakenDb(p);
        loadKnownDbs();
      })
      .catch((err) => setAddDbError(err.message));
  }

  // Remove a database location. If it was the active one, the backend clears
  // the active DB too — a removed database shouldn't stay silently selected.
  function removeSavedDb(p) {
    fetch("./api/kraken-dbs/remove", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path: p }),
    })
      .then((r) => r.json())
      .then((d) => {
        setSettingsDraft((prev) => ({
          ...prev,
          saved_kraken_dbs: d.saved,
          kraken_db: d.cleared_active ? "" : prev.kraken_db,
        }));
        if (d.cleared_active) setKrakenDb("");
        loadKnownDbs();
      })
      .catch(() => {});
  }

  function persistRoots(next) {
    const merged = { ...settingsDraft, ...next };
    setSettingsDraft(merged);
    fetch("./api/config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ projects_root: merged.projects_root, saved_project_roots: merged.saved_project_roots }),
    }).then(() => fetch("./api/config").then((r) => r.json()).then(setSettingsDraft)).catch(() => {});
  }
  function saveCurrentLocation() {
    const cur = (settingsDraft.projects_root || "").trim();
    const list = settingsDraft.saved_project_roots || [];
    if (!cur || list.includes(cur)) return;
    persistRoots({ saved_project_roots: [...list, cur] });
  }
  function removeSavedLocation(p) {
    persistRoots({ saved_project_roots: (settingsDraft.saved_project_roots || []).filter((r) => r !== p) });
  }
  function jumpToLocation(p) {
    if (p) persistRoots({ projects_root: p });
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
      {/* Single hidden file input shared by every project's "Choose files".
          uploadProjRef tracks which project the dialog was opened for. */}
      <input
        ref={uploadInputRef}
        type="file"
        multiple
        accept=".fastq.gz,application/gzip"
        style={{ display: "none" }}
        onChange={(e) => {
          const files = Array.from(e.target.files);
          e.target.value = "";
          if (uploadProjRef.current) uploadFiles(uploadProjRef.current, files);
        }}
      />
      {/* ── Header ─────────────────────────────────────────────── */}
      <header className="app-header">
        <div className="app-brand">
          <img className="app-logo" src="./krona_icon.svg" alt="Krona taxonomy chart icon" />
          <div>
            <h1>
              Kraken ID Parse <span className="version-tag">v{APP_VERSION}</span>
            </h1>
            <p>Classify and isolate reads for species-level identification and contamination screening</p>
          </div>
        </div>
        <div className="header-actions">
          <ThemeToggle />
          <div className="status-pill">
            <span className="dot" data-state={jobStatus} />
            <span>{statusText}</span>
          </div>
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
                <label className="form-label">Kraken2 databases</label>
                {(settingsDraft.saved_kraken_dbs || []).length ? (
                  (settingsDraft.saved_kraken_dbs || []).map((p) => {
                    const info = knownDbs.find((d) => d.path === p);
                    return (
                      <div key={p} style={{ display: "flex", alignItems: "center", gap: 6, fontSize: 12, marginBottom: 3 }}>
                        <span style={{ flex: 1, wordBreak: "break-all", opacity: p === settingsDraft.kraken_db ? 1 : 0.8 }}>
                          {p}
                          {info && info.size_bytes ? ` (${(info.size_bytes / 1073741824).toFixed(1)} GB)` : ""}
                          {info && info.missing ? " ⚠ missing" : ""}
                          {p === settingsDraft.kraken_db ? "  ← active" : ""}
                        </span>
                        <button type="button" className="ghost" title="Use this database"
                          onClick={() => setSettingsDraft((d) => ({ ...d, kraken_db: p }))}
                          disabled={p === settingsDraft.kraken_db}>Use</button>
                        <button type="button" className="ghost" title="Remove from this list (does not delete anything on disk)"
                          onClick={() => removeSavedDb(p)}>✕</button>
                      </div>
                    );
                  })
                ) : (
                  <div className="form-hint">No databases yet — add one below.</div>
                )}
                <div style={{ display: "flex", gap: 6, marginTop: 6 }}>
                  <input
                    style={{ flex: 1 }}
                    placeholder="/path/to/kraken2_db  (directory containing hash.k2d)"
                    value={newDbPath}
                    onChange={(e) => { setNewDbPath(e.target.value); setAddDbError(""); }}
                    onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); addDbLocation(); } }}
                  />
                  <button type="button" className="ghost" onClick={() => openFolderBrowser("add_db")}>Browse…</button>
                  <button type="button" onClick={addDbLocation} disabled={!(newDbPath || "").trim()}>Add</button>
                </div>
                {addDbError ? (
                  <div className="form-hint" style={{ color: "var(--danger, #b00020)" }}>{addDbError}</div>
                ) : null}
                <div className="form-hint">
                  Only databases added here (and in the per-run picker) appear in the dropdowns — nothing is
                  auto-discovered. “Use” makes one active for new runs; remember to click Save.
                </div>
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
                <div style={{ display: "flex", gap: 6 }}>
                  <input
                    style={{ flex: 1 }}
                    value={settingsDraft.projects_root || ""}
                    onChange={(e) => setSettingsDraft((d) => ({ ...d, projects_root: e.target.value }))}
                  />
                  <button type="button" className="ghost" onClick={openFolderBrowser}>Browse…</button>
                </div>
                <div style={{ marginTop: 8 }}>
                  <label className="form-label">Saved locations</label>
                  <span style={{ display: "inline-flex", alignItems: "center", gap: "6px", flexWrap: "wrap" }}>
                    <select value="" onChange={(e) => jumpToLocation(e.target.value)}
                      disabled={!(settingsDraft.saved_project_roots && settingsDraft.saved_project_roots.length)}>
                      <option value="">{settingsDraft.saved_project_roots && settingsDraft.saved_project_roots.length ? "↦ Jump to a saved location…" : "No saved locations yet"}</option>
                      {(settingsDraft.saved_project_roots || []).map((r) => (<option key={r} value={r}>{r}</option>))}
                    </select>
                    <button type="button" className="ghost" onClick={saveCurrentLocation}
                      disabled={!settingsDraft.projects_root || (settingsDraft.saved_project_roots || []).includes(settingsDraft.projects_root)}>★ Save current</button>
                    <button type="button" className="ghost" onClick={() => removeSavedLocation(settingsDraft.projects_root)}
                      disabled={!(settingsDraft.saved_project_roots || []).includes(settingsDraft.projects_root)}>Remove</button>
                  </span>
                </div>
                <div className="form-hint">New projects are created under this root. Projects in the site’s shared projects folder are always visible. Click Save to apply.</div>
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
              <div className="row">
                <input
                  placeholder="New project name (e.g. LSDV_India)"
                  value={newProjectName}
                  onChange={(e) => setNewProjectName(e.target.value.replace(/\s+/g, "_"))}
                  onKeyDown={(e) => { if (e.key === "Enter") createProject(); }}
                  disabled={creatingProject}
                  title="Spaces become underscores. Letters, digits, _ - . are allowed. Created under your personal projects and visible in vSNP too."
                />
                <button onClick={createProject} disabled={creatingProject || !newProjectName.trim()}>
                  {creatingProject ? "Creating…" : "Create"}
                </button>
              </div>
              <div className="form-hint" style={{ marginTop: -4, marginBottom: 8 }}>
                Created under your personal projects root — also visible in vSNP GUI. Add FASTQs to the project’s <code>download/</code> folder.
              </div>
              <div className="list project-list">
                {projectsLoading && <div className="loading-text">Loading projects…</div>}
                {!projectsLoading && projects.length === 0 && (
                  <div className="note">No projects found. Check Settings for the projects path.</div>
                )}
                {projects.map((proj) => (
                  <div
                    key={proj.name}
                    className={`list-item ${activeRun?.project === proj.name || activeProject === proj.name ? "active" : ""}`}
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
                          <div className="empty-msg" style={{ paddingLeft: 4 }}>
                            No FASTQ files yet — add some from the <strong>Inputs</strong> pane on the right.
                          </div>
                        )}
                        {(samples[proj.name] || []).length > 0 && (
                          <div className="sample-item" style={{ display: "flex", alignItems: "center", gap: 8, flexWrap: "wrap" }}>
                            <input
                              type="checkbox"
                              ref={(el) => { if (el) el.indeterminate = checkAllState(proj.name).indeterminate; }}
                              checked={checkAllState(proj.name).checked}
                              disabled={!visibleSamples(proj.name).length}
                              onChange={(e) => toggleCheckAllVisible(proj.name, e.target.checked)}
                              title="Select / deselect every sample shown here (honours the filter)"
                            />
                            <span className="muted" style={{ fontSize: 11 }}>
                              {checkAllState(proj.name).on} of {checkAllState(proj.name).total} selected for run
                            </span>
                            <input
                              type="search"
                              placeholder="Filter samples…"
                              value={projSampleFilter}
                              onChange={(e) => setProjSampleFilter(e.target.value)}
                              style={{ flex: "1 1 140px", minWidth: 120, fontSize: 12 }}
                            />
                          </div>
                        )}
                        {visibleSamples(proj.name).map((s) => {
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

            {/* RIGHT — Inputs (add samples) + batch selection, stacked */}
            <div style={{ display: "flex", flexDirection: "column", gap: 20, minWidth: 0 }}>
              {/* Inputs pane — mirrors vSNP: acts on the active project */}
              <section className="panel">
                <div className="panel-header">
                  <h2>Inputs</h2>
                  {projects.length > 0 && (
                    <select
                      value={activeProject}
                      onChange={(e) => selectProject(e.target.value)}
                      title="Project to add FASTQ files to"
                      style={{ width: "auto", maxWidth: "60%", padding: "6px 10px" }}
                    >
                      {projects.map((p) => (
                        <option key={p.name} value={p.name}>{p.name}</option>
                      ))}
                    </select>
                  )}
                </div>
                {!activeProject ? (
                  <div className="empty-msg">
                    Create a project first (top of the Projects panel), then import, upload, or download FASTQ files into it.
                  </div>
                ) : (
                  <div className="input-columns">
                    {/* ── LEFT: Bring Your Own FASTQ ───────────────────────── */}
                    <div className="input-column">
                      <h3>Bring Your Own FASTQ</h3>

                      {/* Import from a server path */}
                      <div className="row" style={{ margin: 0 }}>
                        <input
                          placeholder="/path/to/folder or .fastq.gz file"
                          value={addPath[activeProject] || ""}
                          onChange={(e) => setAddPath((m) => ({ ...m, [activeProject]: e.target.value }))}
                          onKeyDown={(e) => { if (e.key === "Enter") linkLocal(activeProject); }}
                        />
                        <button className="ghost action" onClick={() => linkLocal(activeProject)} disabled={!(addPath[activeProject] || "").trim()}>Link</button>
                      </div>
                      <div className="form-hint">Symlinks every .fastq.gz found — no copying.</div>

                      {/* Upload / drag & drop */}
                      <div className="block">
                        <h3>Upload / Drag &amp; Drop</h3>
                        <div
                          className="dropzone"
                          onDragOver={(e) => e.preventDefault()}
                          onDrop={(e) => { e.preventDefault(); uploadFiles(activeProject, e.dataTransfer.files); }}
                        >
                          <button type="button" onClick={() => pickFiles(activeProject)}>Choose Files</button>
                          <span className="drop-hint">Or drop FASTQ.GZ files here</span>
                        </div>
                        {addStatus[activeProject] && <div className="note" style={{ marginBottom: 0 }}>{addStatus[activeProject]}</div>}
                      </div>

                      {/* Files already in download/ */}
                      {inputsByProj[activeProject]?.files?.length > 0 && (
                        <div className="block">
                          <h3 style={{ display: "flex", alignItems: "center", gap: 8 }}>
                            <span style={{ flex: 1 }}>
                              Files in download/
                              <span className="muted" style={{ marginLeft: 6, fontWeight: 400, fontSize: 12 }}>
                                ({inputsByProj[activeProject].count}, {fmtSize(inputsByProj[activeProject].total_bytes)})
                              </span>
                            </span>
                            <button className="ghost" style={{ fontSize: 11, padding: "2px 8px" }} onClick={() => loadInputs(activeProject)} title="Refresh">Refresh</button>
                          </h3>
                          <div className="input-files">
                            {inputsByProj[activeProject].files.map((f) => (
                              <div key={f.name} className="input-file-row">
                                <span className="file-name" title={f.name} style={{ flex: 1 }}>{f.name}</span>
                                <span className="file-size">{fmtSize(f.size)}</span>
                                <button className="ghost" style={{ fontSize: 11, padding: "2px 7px" }} title="Remove from download/" onClick={() => deleteInput(activeProject, f.name)}>✕</button>
                              </div>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>

                    {/* ── RIGHT: SRA Download ──────────────────────────────── */}
                    <div className="input-column">
                      <h3>SRA Download</h3>
                      <textarea
                        rows={6}
                        placeholder={"SRR/ERR/DRR or SRX/SRS/PRJNA accessions\n(one per line)"}
                        value={sraText[activeProject] || ""}
                        onChange={(e) => setSraText((m) => ({ ...m, [activeProject]: e.target.value }))}
                        style={{ resize: "vertical", fontFamily: "inherit" }}
                      />
                      <button
                        style={{ width: "100%" }}
                        onClick={() => sraDownload(activeProject)}
                        disabled={!parseAccessions(sraText[activeProject]).length || running}
                      >
                        Download{parseAccessions(sraText[activeProject]).length ? ` (${parseAccessions(sraText[activeProject]).length})` : ""}
                      </button>
                      <div className="form-hint">Runs in the background; progress appears in the Pipeline Log.</div>
                    </div>
                  </div>
                )}
              </section>

              {/* Batch selection pane */}
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
          <div className="row-grid row-grid-single">
            {/* LEFT — configure & run */}
            <section className="panel">
              <h2>Configure &amp; Run</h2>

              <div className="form-section">
                <label className="checkbox-label" style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
                  <input
                    type="checkbox"
                    checked={krakenOnly}
                    onChange={(e) => { setKrakenOnly(e.target.checked); if (e.target.checked) setNoBlast(false); }}
                    disabled={running}
                  />
                  <span>Kraken only (Krona graph, no read parsing)</span>
                </label>
                <div className="note" style={{ marginTop: 4 }}>
                  Runs Kraken2 and produces the Krona graph only — skips read parsing, assembly, and BLAST. No target taxon needed.
                </div>
                <label className="checkbox-label" style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer", marginTop: 10 }}>
                  <input
                    type="checkbox"
                    checked={noBlast}
                    onChange={(e) => { setNoBlast(e.target.checked); if (e.target.checked) setKrakenOnly(false); }}
                    disabled={running}
                  />
                  <span>Parse reads only (skip BLAST)</span>
                </label>
                <div className="note" style={{ marginTop: 4 }}>
                  Runs Kraken2 and extracts the target taxon's reads, then stops — skips assembly, BLAST, and coverage. Leaves the parsed FASTQ.gz reads. Requires a target taxon.
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">Target Taxon</label>
                <select
                  value={taxon}
                  onChange={(e) => setTaxon(e.target.value)}
                  disabled={running || krakenOnly}
                >
                  <option value="">Select a target taxon…</option>
                  {taxonPresets.map((p) => (
                    <option key={p} value={p}>{p}</option>
                  ))}
                </select>
                <div className="add-taxon-row" style={{ display: "flex", gap: 6, marginTop: 8 }}>
                  <input
                    placeholder="Add a new taxon to the list…"
                    value={newTaxon}
                    onChange={(e) => setNewTaxon(e.target.value)}
                    onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); addTaxon(); } }}
                    disabled={running || addingTaxon}
                    style={{ flex: 1 }}
                  />
                  <button
                    type="button"
                    className="ghost"
                    onClick={addTaxon}
                    disabled={running || addingTaxon || !newTaxon.trim()}
                  >
                    {addingTaxon ? "Adding…" : "+ Add"}
                  </button>
                </div>
                <div className="note" style={{ marginTop: 4 }}>
                  New taxa are saved to the shared list (this install’s config/taxa.yaml) and appear in this dropdown and the vSNP GUI.
                </div>
              </div>

              <div className="form-section">
                <label className="form-label">
                  Kraken2 DB path
                  {!krakenDb && <span style={{ color: "var(--danger)", marginLeft: 6, fontSize: 11 }}>⚠ not configured</span>}
                </label>
                {knownDbs.length ? (
                  <select
                    value={knownDbs.some((d) => d.path === krakenDb) ? krakenDb : ""}
                    onChange={(e) => { if (e.target.value) setKrakenDb(e.target.value); }}
                    disabled={running}
                    style={{ marginBottom: 6 }}
                  >
                    <option value="">— switch to a saved database (managed in Settings) —</option>
                    {knownDbs.map((d) => (
                      <option key={d.path} value={d.path}>
                        {d.name}{d.missing ? " ⚠ missing" : ""} — {d.path}
                      </option>
                    ))}
                  </select>
                ) : null}
                <input
                  placeholder="/path/to/kraken2_db  (directory containing hash.k2d)"
                  value={krakenDb}
                  onChange={(e) => setKrakenDb(e.target.value)}
                  disabled={running}
                />
              </div>

              <div className="form-section">
                <label className="form-label">BLAST DB path (or name)</label>
                <input
                  placeholder="nt  or  /path/to/blast_db"
                  value={blastDb}
                  onChange={(e) => setBlastDb(e.target.value)}
                  disabled={running || krakenOnly || noBlast}
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

          </div>
        )}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Pipeline Log                                   */}

        {/* ════════════════════════════════════════════════════════ */}
        {/* SECTION: Results — every completed sample, not just the last  */}
        {/* ════════════════════════════════════════════════════════ */}
        <div className="row-header">
          <h2>Results</h2>
          <button className="ghost" onClick={() => setShowResultsPane(!showResultsPane)}>
            {showResultsPane ? "Hide" : "Show"}
          </button>
        </div>
        {showResultsPane && (
          <div className="row-grid row-grid-split">
            {/* LEFT — Current Run (live status for the batch in flight) */}
<section className="panel">
              <div className="panel-header">
                <h2>Current Run</h2>
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
            {/* RIGHT — every completed sample, searchable (vSNP Step 1 model) */}
            <ResultsPane
              project={activeProject}
              results={results}
              columns={RESULT_COLUMNS}
              labels={{ entity: "sample", sampleHeader: "Sample" }}
            />
          </div>
        )}
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

      <CitationFooter tool="kraken_id_parse_gui" />

      {folderBrowser.open && (
        <div
          onClick={() => setFolderBrowser((s) => ({ ...s, open: false }))}
          style={{ position: "fixed", inset: 0, background: "rgba(0,0,0,0.45)", display: "flex", alignItems: "center", justifyContent: "center", zIndex: 1000 }}
        >
          <div
            onClick={(e) => e.stopPropagation()}
            style={{ background: "var(--panel, #fff)", color: "inherit", borderRadius: 10, width: "min(640px, 92vw)", maxHeight: "80vh", display: "flex", flexDirection: "column", boxShadow: "0 10px 40px rgba(0,0,0,0.3)" }}
          >
            <div style={{ padding: "12px 16px", borderBottom: "1px solid var(--border, #ddd)", fontWeight: 700 }}>
              {folderBrowser.target === "add_db" ? "Select a Kraken2 database folder" : "Select a projects root"}
            </div>
            <div style={{ padding: "10px 16px", display: "flex", gap: 6, alignItems: "center" }}>
              <button type="button" className="ghost" disabled={!folderBrowser.parent || folderBrowser.loading} onClick={() => browseDirs(folderBrowser.parent)}>↑ Up</button>
              <input
                style={{ flex: 1 }}
                value={folderBrowser.path}
                onChange={(e) => setFolderBrowser((s) => ({ ...s, path: e.target.value }))}
                onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); browseDirs(folderBrowser.path); } }}
              />
              <button type="button" className="ghost" onClick={() => browseDirs(folderBrowser.path)}>Go</button>
            </div>
            <div style={{ flex: 1, overflow: "auto", padding: "0 16px", minHeight: 160 }}>
              {folderBrowser.loading ? (
                <div className="note" style={{ padding: 12 }}>Loading…</div>
              ) : folderBrowser.error ? (
                <div className="note" style={{ padding: 12, color: "var(--danger, #c00)" }}>{folderBrowser.error}</div>
              ) : folderBrowser.entries.length === 0 ? (
                <div className="note" style={{ padding: 12 }}>No sub-folders here.</div>
              ) : (
                folderBrowser.entries.map((e) => (
                  <div
                    key={e.path}
                    onClick={() => browseDirs(e.path)}
                    style={{ padding: "7px 8px", cursor: "pointer", borderRadius: 6, display: "flex", gap: 8, alignItems: "center" }}
                    onMouseEnter={(ev) => (ev.currentTarget.style.background = "var(--panel-2, #f0f0f0)")}
                    onMouseLeave={(ev) => (ev.currentTarget.style.background = "transparent")}
                  >
                    <span>📁</span><span>{e.name}</span>
                  </div>
                ))
              )}
            </div>
            <div style={{ padding: "12px 16px", borderTop: "1px solid var(--border, #ddd)", display: "flex", justifyContent: "flex-end", gap: 8 }}>
              <button type="button" className="ghost" onClick={() => setFolderBrowser((s) => ({ ...s, open: false }))}>Cancel</button>
              <button type="button" onClick={chooseFolder} disabled={folderBrowser.loading || !folderBrowser.path}>Select this folder</button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
