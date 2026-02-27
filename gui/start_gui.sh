#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Auto-detect conda base
CONDA_BASE=""
if [ -n "${CONDA_EXE:-}" ]; then
  CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
fi

if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  export PATH="$CONDA_BASE/bin:$PATH"
fi

pick_free_port() {
  local start=$1
  local end=$2
  for p in $(seq "$start" "$end"); do
    if ! lsof -iTCP:${p} -sTCP:LISTEN >/dev/null 2>&1; then
      echo "$p"
      return 0
    fi
  done
  return 1
}

BACKEND_PORT="$(pick_free_port 8000 8015)" || { echo "No free backend port found (8000-8015)." >&2; exit 1; }
VITE_PORT="$(pick_free_port 5173 5195)" || { echo "No free Vite port found (5173-5195)." >&2; exit 1; }

# Print detected paths for easy copy/paste into Settings
echo "================================================"
echo "  Kraken ID Parse GUI — Detected Paths"
echo "================================================"
echo "  GUI root:        $ROOT_DIR"
if [ -n "$CONDA_BASE" ]; then
  KRAKEN_ENV="$CONDA_BASE/envs/kraken_id_parse"
  if [ -d "$KRAKEN_ENV" ]; then
    echo "  Conda env path:  $KRAKEN_ENV"
  else
    echo "  Conda env path:  (kraken_id_parse env not found — set in Settings)"
  fi
else
  echo "  Conda env path:  (conda not detected — open a conda-enabled shell)"
fi
echo "  Backend port:    $BACKEND_PORT"
echo "  Vite port:       $VITE_PORT"
echo "================================================"
echo ""

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda not found in PATH. Please open a conda-enabled shell." >&2
fi

if ! command -v npm >/dev/null 2>&1; then
  echo "npm not found. Install Node.js (brew install node)." >&2
fi

# Backend
(
  cd "$ROOT_DIR/backend"
  if [ ! -d .venv ]; then
    python3 -m venv .venv
  fi
  source .venv/bin/activate
  pip install -r requirements.txt
  uvicorn app.main:app --reload --port "$BACKEND_PORT"
) &
BACK_PID=$!

# Frontend
(
  cd "$ROOT_DIR/frontend"
  if [ ! -d node_modules ]; then
    npm install
  fi
  VITE_API_URL="http://localhost:${BACKEND_PORT}" npm run dev -- --port "$VITE_PORT" --strictPort
) &
FRONT_PID=$!

sleep 4
open "http://localhost:${VITE_PORT}" || true

echo "Backend PID: $BACK_PID"
echo "Frontend PID: $FRONT_PID"

trap 'kill $BACK_PID $FRONT_PID' INT TERM
wait
