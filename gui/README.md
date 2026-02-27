# Kraken ID Parse GUI

Local web + desktop interface for running the Kraken ID Parse pipeline.

## Requirements

- macOS or Linux
- Python 3.9+
- Node.js 18+
- Kraken ID Parse pipeline (the GUI calls `bin/run_with_config.py` from this repo)

## Setup

### 1. Install backend dependencies

```bash
cd gui/backend
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Install frontend dependencies

```bash
cd gui/frontend
npm install
```

### 3. Install Electron (optional — for desktop app)

```bash
cd gui/electron
npm install
```

## Running

### Browser mode (recommended for getting started)

```bash
cd gui
./start_gui.sh
```

This starts the Flask backend (port 8000) and Vite dev server (port 5173), then opens the GUI in your browser.

### Electron desktop app

```bash
cd gui
./start_electron.sh
```

### Linux

```bash
cd gui
./start_gui_linux.sh
```

### Individual components (for development)

```bash
./start_backend.sh    # Flask API only (port 8000)
./start_frontend.sh   # Vite dev server only (port 5173)
```

## Settings

Configure in the GUI Settings panel:

- **Kraken repo path** (required) — root of this repository (e.g. `~/kraken_id_parse_gui`)
- **Projects root** (required) — directory where projects and run outputs are stored
- **Conda env path** (optional) — adds the environment's `bin/` to PATH when running the pipeline

## Workflow

1. Create a project
2. Link a FASTQ folder (symlinks into the project's `download/` directory)
3. Choose a preset and run
4. Watch logs in real time and review outputs

## Sample / Test Data

Sample data files are not included in the repository to keep it lightweight. To get test data for running the GUI:

### Orbivirus test (Illumina paired-end)

```bash
# Download from SRA
sra_number="SRR9598511"
wget -O "${sra_number}.fastq.gz" \
  "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra_number}/${sra_number}"
fasterq-dump -S ${sra_number}.fastq.gz

# Prepare files
rm ${sra_number}.fastq.gz
mv ${sra_number}.fastq.gz_1.fastq ${sra_number}_R1.fastq
mv ${sra_number}.fastq.gz_2.fastq ${sra_number}_R2.fastq
pigz *fastq
```

Then in the GUI: create a project, link the folder containing these FASTQs, select the **orbivirus** preset, and run.

### Mycobacterium tuberculosis test (Illumina paired-end)

```bash
sra_number="SRR28623786"
wget -O "${sra_number}.fastq.gz" \
  "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra_number}/${sra_number}"
fasterq-dump -S ${sra_number}.fastq.gz

rm ${sra_number}.fastq.gz
mv ${sra_number}.fastq.gz_1.fastq ${sra_number}_R1.fastq
mv ${sra_number}.fastq.gz_2.fastq ${sra_number}_R2.fastq
pigz *fastq
```

Then in the GUI: create a project, link the folder, select the **mtb** preset, and run.

### VCF test data (for Step 2 / vSNP workflows)

A small set of zero-coverage VCFs for MTBC lineage testing is available separately. Contact the maintainers or see the pipeline's test documentation for download instructions.

## Project Structure

```
gui/
├── backend/             # Flask/FastAPI backend (Python)
│   ├── app/main.py      #   API routes
│   ├── app/config.py    #   Settings management
│   ├── app/jobs.py      #   Job tracking
│   ├── app/projects.py  #   Project CRUD
│   ├── app/refs.py      #   Reference management
│   └── app/sra.py       #   SRA download support
├── frontend/            # React/Vite frontend
│   └── src/App.jsx      #   Main application component
├── electron/            # Electron desktop wrapper
│   └── main.js          #   Electron main process
├── docs/                # User and developer documentation
├── assets/icons/        # Application icons
├── start_gui.sh         # Launch browser mode (backend + frontend)
├── start_electron.sh    # Launch Electron desktop app
└── start_gui_linux.sh   # Launch on Linux
```

