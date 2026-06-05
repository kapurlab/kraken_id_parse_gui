# Handoff — Shared taxa.yaml, "skip BLAST" parse-only mode, and vSNP auto-import

**Date:** 2026-06-04
**Author:** tks5563@psu.edu (with Claude Code)
**Repos:**
- kraken: `/srv/kapurlab/tools/kraken_id_parse_gui` (github.com/kapurlab/kraken_id_parse_gui)
- vSNP:   `/srv/kapurlab/tools/vsnp_gui` (github.com/kapurlab/vsnp_gui)

Picks up after `HANDOFF_inputs_pane_styling.md`. This change spans **both** GUIs
(confirmed with the user: both GUIs get the taxa-from-yaml + add-name + skip-BLAST
controls; auto-import fires for any run that produced parsed reads).

---

## ⚠️ READ FIRST — how to make changes actually appear in the OOD app

Editing files in `/srv/kapurlab/tools/...` is NOT enough. Changes only show up
after you do ALL of this (do it for EVERY GUI change — it cost real time once):

1. **Rebuild the frontend** for any `frontend/src` edit — `dist/` is gitignored
   and is what uvicorn serves:
   ```
   cd /srv/kapurlab/tools/<repo>/frontend && npm run build
   ```
2. **Commit AND push** to `origin/<branch>`:
   ```
   cd /srv/kapurlab/tools/<repo>
   git add <source files> config/taxa.yaml   # NOT env/, __pycache__/, backend/jobs/, scratch files
   git commit -m "..."
   git push origin <branch>
   ```
3. **Relaunch the OOD session** (from the dashboard) and, for a **dev** app,
   **select the matching branch**.

Why both commit *and* push are required:
- **Dev OOD apps** (`*_dev`) run `git worktree add --detach <tmp> origin/<branch>`
  then `npm run build` — they build from the **pushed** GitHub ref, so they never
  see uncommitted/unpushed working-tree edits. (See `*_dev/template/before.sh.erb`
  / `script.sh.erb` under `/var/www/ood/apps/sys/`.)
- **Prod OOD apps** `cd /srv/kapurlab/tools/<repo>/backend` and serve the working
  tree directly, so they reflect edits after a session restart + frontend rebuild
  — but the team standard is to test on a dev branch first, then merge to `main`.

This work shipped on feature branches (kraken `feature/create-projects` @ `6ff97d7`,
vsnp `feature/step1-run-kraken` @ `99852f7`); merge to `main` before pointing prod
users at it.

---

## What was built

### 1. Taxon presets now come from a shared YAML (was hardcoded in JS)

New file **`kraken_id_parse_gui/config/taxa.yaml`** — a plain YAML sequence of
taxon search names. Single source of truth read by **both** GUIs. Editable by
hand or via the GUI "Add a search name…" control.

Parsing is **dependency-free** (a tiny `- name` line parser) so it works whether
or not PyYAML is installed — important because the vSNP env does **not** list
PyYAML in requirements.txt. The file is still valid YAML (verified with PyYAML).

### 2. "Parse reads only (skip BLAST)" mode

A new pipeline mode that runs Kraken2 + taxonomic **read parsing** and then
**stops** — skipping assembly, BLAST, and coverage. It leaves the parsed
`*_R1.fastq.gz` / `*_R2.fastq.gz` reads for the target taxon. This is distinct
from the existing "Kraken only (Krona)" mode, which skips parsing entirely.

### 3. vSNP auto-imports parsed reads into project inputs

When a Kraken run launched **from the vSNP GUI** finishes successfully and
produced parsed reads (full **or** parse-only mode), the backend copies them into
`<project>/download/` (the Step-1 inputs folder), renamed so vSNP treats them as
a **new, distinct sample** — ready to re-run through vSNP with no manual copying.

**The rename is the crux.** vSNP keys a sample on the text left of the first
underscore, so `13-1941-6_Mycobacterium_tuberculosis_complex_R1.fastq.gz` would
collapse back to sample `13-1941-6` and **overwrite the original run**. The
importer replaces every underscore in the sample-identifying stem with `-`,
keeping only the `_R1`/`_R2` read-tag underscore:

```
13-1941-6_Mycobacterium_tuberculosis_complex_R1.fastq.gz
  ->  13-1941-6-Mycobacterium-tuberculosis-complex_R1.fastq.gz
  vSNP sample = 13-1941-6-Mycobacterium-tuberculosis-complex   (unique)
```

---

## Files changed

### kraken repo
| File | Change |
|---|---|
| `config/taxa.yaml` | **NEW** — shared taxon list |
| `bin/kraken_id_parse.py` | `--no-blast` flag + validation + early-exit after parsing (before assembly/BLAST) |
| `backend/app/main.py` | `_TAXA_YAML` + `_read_taxa`/`_write_taxa`; `GET`/`POST /api/taxa`; `no_blast` on `RunPayload`; pass `--no-blast`; job label |
| `frontend/src/App.jsx` | taxa loaded from `./api/taxa`; "Add a search name" control; "Parse reads only (skip BLAST)" checkbox (mutually exclusive with "Kraken only"); sends `no_blast` |
| `frontend/dist/` | rebuilt (`index-8d32ed75.js`) |

### vSNP repo
| File | Change |
|---|---|
| `backend/app/main.py` | `_KRAKEN_TAXA_YAML` + `_read_kraken_taxa`/`_write_kraken_taxa`; `GET`/`POST /api/kraken/taxa`; `_dash_delimited_import_name` + `_import_parsed_reads`; `mode="parse_only"` (→ `--no-blast`); `finalize_callback` on `kraken_run` that auto-imports on success |
| `frontend/src/App.jsx` | taxa loaded from `./api/kraken/taxa`; "Add a search name" control; "Parse reads only (skip BLAST)" mode card; refreshes Inputs pane on run success; updated success note |
| `frontend/dist/` | rebuilt (`index-a9033af6.js`) |

The auto-import hook uses the **existing** `finalize_callback` mechanism in
`vsnp_gui/backend/app/jobs.py` (soft-fail per the T-07 policy — a copy failure
is logged, never masks the run result).

---

## Verified
- Both frontends build clean (`npm run build`).
- Both backends `py_compile` clean.
- `taxa.yaml` add/round-trip preserves order, dedupes case-insensitively, and
  stays valid YAML (incl. names with `:`).
- Rename logic produces unique, non-colliding vSNP sample names.

## NOT yet done
- **Live end-to-end run** on real FASTQs through OOD (build/compile only).
- **Not committed/pushed** to git in either repo.
- vSNP backend was syntax-checked with the kraken env's python (the vSNP env
  `/home/vxk1/miniforge3/envs/vsnp3` is not readable by tks5563). A full import
  under the vSNP env is worth doing in a real session.

## How to deploy / test
1. Kraken GUI: served from `frontend/dist/` (already rebuilt). Restart the
   backend/OOD session to pick up `main.py` changes.
2. vSNP GUI: `frontend/dist/` rebuilt; restart its OOD session (uvicorn serves
   dist on startup) to pick up `main.py`.
3. Smoke test taxa endpoints:
   - `GET ./api/taxa` (kraken) and `GET ./api/kraken/taxa` (vSNP) → same list.
   - Add a name in either GUI → appears in the other.
4. In vSNP: open a Step-1 sample → Kraken modal → pick "Parse reads only" + a
   taxon → Run. On success the parsed reads should appear in the Inputs pane
   (download/) with the `-`-delimited name.
