# Kraken GUI OOD handoff for wgs3

## Current state

Local branch:

```bash
feature/ood-html-reporting
```

Base branch:

```bash
origin/feature/ood-web-gui
```

wgs3 currently has:

```bash
/srv/kapurlab/tools/kraken_id_parse_gui
branch: feature/ood-web-gui
commit: 78b3b39
```

This means wgs3 does **not** yet have the local HTML/PDF/report-results fixes.

## What changed in this branch

- Hardened GUI job execution so pipeline commands run with `shell=False`.
- Kept job logs/session metadata under `backend/jobs/`.
- Added structured reporting under `bin/reporting/`.
- Added `run_manifest.json`, `report.html`, and `report.pdf` generation.
- Replaced brittle LaTeX-only reporting path with an HTML-first report plus browser PDF export.
- Made the PDF layout closely resemble the old LaTeX report:
  - FASTQ quality
  - Bracken pie chart
  - Kraken summary
  - detailed lineage
  - assembly
  - BLAST identification
  - coverage figure
- Fixed report HTML asset links when opened from the GUI.
- Curated the GUI Results pane so it shows only primary outputs by default:
  - Report PDF
  - Report HTML
  - latest stats workbook
  - latest Krona report
  - BLAST summaries
  - de novo/reference-guided FASTA
  - latest coverage graph PDF
  - Bracken workbook
  - pipeline log
- Preserved intermediate files in the run folder for audit/debugging, but hid them from the default GUI results list.
- Added Matplotlib style fallback so older environments do not fail on `seaborn-v0_8-colorblind`.
- Relaxed the backend uvicorn requirement for Python 3.8 compatibility.

## What transfers

Transfer via git:

```bash
backend/
bin/
conda_setup/
frontend/
ood/
README.md
```

Do **not** transfer local user/runtime state:

```bash
backend/jobs/
backend/app/__pycache__/
bin/__pycache__/
bin/reporting/__pycache__/
frontend/node_modules/
env/
~/.config/kraken_id_parse_gui/config.json
```

The local config file is per-machine/per-user and should not be copied to wgs3.

## wgs3 paths verified

Tool source:

```bash
/srv/kapurlab/tools/kraken_id_parse_gui
```

OOD app source in repo:

```bash
ood/apps/kraken_id_parse_gui/
ood/apps/kraken_id_parse_gui_dev/
```

OOD deployed app locations:

```bash
/var/www/ood/apps/sys/kraken_id_parse_gui/
/var/www/ood/apps/sys/kraken_id_parse_gui_dev/
```

Databases currently verified on wgs3:

```bash
/srv/kapurlab/databases/kraken2/k2_standard_08gb
/srv/kapurlab/databases/blast/ref_prok_rep_genomes
```

I did not find `nt_viruses` under `/srv/kapurlab` during this check. If viral work is a primary use case, download/install BLAST `nt_viruses` on wgs3 and set the GUI BLAST DB path to that for viral runs.

## Recommended handoff path

### 1. Commit and push this branch

From the local worktree:

```bash
cd /Users/vivekkapur/kraken/kraken_ood_html_report
git status
git add backend bin conda_setup frontend ood README.md TOD_WGS3_OOD_HANDOFF.md
git commit -m "Add OOD HTML/PDF reporting and curated results"
git push origin feature/ood-html-reporting
```

### 2. Tod tests through the OOD dev app first

The dev OOD app already supports a branch field. Use:

```bash
feature/ood-html-reporting
```

The dev app creates a temporary worktree from `origin/<branch>`, builds the frontend, and starts uvicorn from that branch.

### 3. If dev test passes, update the shared wgs3 repo

On wgs3:

```bash
cd /srv/kapurlab/tools/kraken_id_parse_gui
git fetch origin
git checkout feature/ood-html-reporting
git pull --ff-only origin feature/ood-html-reporting
```

Build frontend:

```bash
cd /srv/kapurlab/tools/kraken_id_parse_gui/frontend
npm install
npm run build
```

Update/install env if needed:

```bash
cd /srv/kapurlab/tools/kraken_id_parse_gui
conda env update -p ./env -f conda_setup/environment.yml --prune
```

If the env update is risky, create a new env first and point the OOD script at it.

### 4. Deploy OOD app files if they changed

From repo root on wgs3:

```bash
sudo rsync -av ood/apps/kraken_id_parse_gui/ \
  /var/www/ood/apps/sys/kraken_id_parse_gui/

sudo rsync -av ood/apps/kraken_id_parse_gui_dev/ \
  /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
```

The app should be available at:

```text
/pun/sys/dashboard/batch_connect/sys/kraken_id_parse_gui/session_contexts/new
```

## User config on wgs3

The GUI writes per-user config here:

```bash
~/.config/kraken_id_parse_gui/config.json
```

For Tod, either set paths in the GUI Settings panel or seed:

```json
{
  "projects_root": "/home/tks5563/projects",
  "shared_projects_root": "/srv/kapurlab/projects",
  "kraken_db": "/srv/kapurlab/databases/kraken2/k2_standard_08gb",
  "blast_db": "/srv/kapurlab/databases/blast/ref_prok_rep_genomes"
}
```

For viral runs, use `nt_viruses` once installed, for example:

```json
{
  "blast_db": "/srv/kapurlab/databases/blast/nt_viruses"
}
```

## Suggested smoke tests

Orbivirus:

```text
Taxon: Orbivirus
Kraken DB: /srv/kapurlab/databases/kraken2/k2_standard_08gb
BLAST DB: /srv/kapurlab/databases/blast/nt_viruses
```

MTB/bacterial:

```text
Taxon: Mycobacterium tuberculosis complex
Kraken DB: /srv/kapurlab/databases/kraken2/k2_standard_08gb
BLAST DB: /srv/kapurlab/databases/blast/ref_prok_rep_genomes
```

Expected primary outputs in the GUI:

```text
Report PDF
Report HTML
Run statistics workbook
Krona taxonomy report
BLAST summary
Consensus BLAST summary
De novo assembly FASTA
Reference-guided FASTA
Coverage graph PDF
Bracken results workbook
Pipeline log
```

