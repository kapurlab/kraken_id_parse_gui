#!/bin/bash
set -euo pipefail

PORT="${PORT:-8000}"

cd "$(dirname "$0")/backend"
if [ ! -d .venv ]; then
  python3 -m venv .venv
fi
source .venv/bin/activate
pip install -r requirements.txt
exec uvicorn app.main:app --reload --port "$PORT"
