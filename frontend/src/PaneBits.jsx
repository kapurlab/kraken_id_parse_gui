import { useEffect, useRef, useState } from "react";
import Elapsed from "./Elapsed";
import { askOnPage } from "./pageAsk.jsx";

// Small pieces shared with the vSNP GUI so the two tools look and behave the
// same: the pane Expand button, a button that shows it is working, the disk
// space line, and byte formatting.

const ExpandIcon = () => (
  <svg viewBox="0 0 16 16" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="1.6"
    strokeLinecap="round" strokeLinejoin="round">
    <path d="M2.5 6V2.5H6M10 2.5h3.5V6M13.5 10v3.5H10M6 13.5H2.5V10" />
  </svg>
);

const CollapseIcon = () => (
  <svg viewBox="0 0 16 16" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="1.6"
    strokeLinecap="round" strokeLinejoin="round">
    <path d="M6 2.5V6H2.5M13.5 6H10V2.5M10 13.5V10h3.5M2.5 10H6v3.5" />
  </svg>
);

// Expand on a pane's header: that pane over most of the window, and back
// (Esc, Collapse, or a click outside). Same control as the vSNP GUI.
export function PaneExpandButton({ expanded, onToggle, what }) {
  return (
    <button
      type="button"
      className="ghost action pane-expand"
      aria-pressed={expanded}
      onClick={onToggle}
      title={expanded
        ? "Put this pane back in its place on the page (Esc)"
        : `Show ${what} over most of the window; Esc puts it back`}
    >
      {expanded ? <CollapseIcon /> : <ExpandIcon />}
      {expanded ? "Collapse" : "Expand"}
    </button>
  );
}

// State + Esc handling for one expanded pane at a time.
export function usePaneExpand() {
  const [expanded, setExpanded] = useState("");
  useEffect(() => {
    if (!expanded) return undefined;
    const onKey = (e) => {
      if (e.key !== "Escape" || e.defaultPrevented) return;
      if (document.querySelector(".page-ask-backdrop, .modal-backdrop")) return;
      setExpanded("");
    };
    window.addEventListener("keydown", onKey);
    const before = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    return () => { window.removeEventListener("keydown", onKey); document.body.style.overflow = before; };
  }, [expanded]);
  const paneClass = (base, id) => (expanded === id ? `${base} is-expanded` : base);
  const expandButton = (id, what) => (
    <PaneExpandButton
      expanded={expanded === id}
      onToggle={() => setExpanded((x) => (x === id ? "" : id))}
      what={what}
    />
  );
  return { expanded, setExpanded, paneClass, expandButton };
}

/**
 * A button that acknowledges the click at once: it disables itself and shows
 * `busyLabel` with a ticking counter until the async onClick settles, so a
 * slow server action never looks like a click that did nothing.
 */
export function BusyButton({ onClick, busyLabel, children, disabled, onError, ...rest }) {
  const [busy, setBusy] = useState(false);
  const mounted = useRef(true);
  useEffect(() => () => { mounted.current = false; }, []);
  return (
    <button
      {...rest}
      disabled={disabled || busy}
      aria-busy={busy ? "true" : undefined}
      onClick={async (e) => {
        if (busy) return;
        setBusy(true);
        try {
          await onClick(e);
        } catch (err) {
          if (onError) onError(err);
          else askOnPage({ notice: true, message: err?.message || String(err) });
        } finally {
          if (mounted.current) setBusy(false);
        }
      }}
    >
      {busy ? (<>{busyLabel || "Working…"} <Elapsed /></>) : children}
    </button>
  );
}

export function formatBytes(bytes) {
  if (!bytes) return "0 B";
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / 1024 / 1024).toFixed(1)} MB`;
  return `${(bytes / 1024 / 1024 / 1024).toFixed(2)} GB`;
}

export function DiskSpace({ disk, style }) {
  if (!disk || !disk.total) return null;
  const pctFree = disk.free / disk.total;
  const low = pctFree < 0.1;
  return (
    <span
      className="muted"
      style={{ fontSize: 12, color: low ? "var(--danger)" : undefined, ...(style || {}) }}
      title={`${formatBytes(disk.free)} free of ${formatBytes(disk.total)} (${Math.round(pctFree * 100)}% free)`}
    >
      {formatBytes(disk.free)} free of {formatBytes(disk.total)}
      {low ? " — low space" : ""}
    </span>
  );
}
