# Kraken ID Parse GUI

Local web + desktop interface for running the Kraken ID Parse pipeline.

## Requirements

- macOS or Linux
- Python 3.9+
- Node.js 18+
- Kraken ID Parse repo (this GUI calls its `bin/run_with_config.py`)

## Quick Start

```bash
cd kraken_gui
./start_gui.sh
```

Electron desktop app:

```bash
./start_electron.sh
```

## Settings

Required:

- Kraken repo path (e.g. `~/kraken_id_parse_gui`)
- Projects root (where runs and outputs live)

Optional:

- Conda env path (adds `bin` to PATH when running)

## Workflow

1. Create a project
2. Link a FASTQ folder (symlinks into the project download folder)
3. Choose a preset and run
4. Watch logs and review outputs

