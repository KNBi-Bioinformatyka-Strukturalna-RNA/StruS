#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNAPOLIS_URL="https://github.com/tzok/rnapolis-py"
RNAPOLIS_DIR="$SCRIPT_DIR/rnapolis-py"

if command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
else
  echo "ERROR: python/python3 not found in PATH" >&2
  exit 1
fi

if [ ! -d "$SCRIPT_DIR/.venv" ]; then
  "$PYTHON_BIN" -m venv "$SCRIPT_DIR/.venv"
fi

if [ -x "$SCRIPT_DIR/.venv/bin/python" ]; then
  VENV_PY="$SCRIPT_DIR/.venv/bin/python"
elif [ -x "$SCRIPT_DIR/.venv/Scripts/python.exe" ]; then
  VENV_PY="$SCRIPT_DIR/.venv/Scripts/python.exe"
else
  echo "ERROR: Cannot find Python inside .venv" >&2
  exit 1
fi

"$VENV_PY" -m pip install --upgrade pip setuptools wheel
"$VENV_PY" -m pip install -r "$SCRIPT_DIR/requirements.txt"

if [ -d "$RNAPOLIS_DIR/.git" ]; then
  git -C "$RNAPOLIS_DIR" pull --ff-only
else
  git clone "$RNAPOLIS_URL" "$RNAPOLIS_DIR"
fi

"$VENV_PY" -m pip install -e "$RNAPOLIS_DIR"

SCRIPT_DIR="$SCRIPT_DIR" "$VENV_PY" - <<'PY'
import os
import re
from pathlib import Path

strus_path = Path(os.environ["SCRIPT_DIR"]) / "StruS.py"
text = strus_path.read_text(encoding="utf-8")
expected = 'ANNOTATOR = str(SCRIPT_DIR / "rnapolis-py" / "src" / "rnapolis" / "annotator.py")'

if expected not in text:
  text = re.sub(r'ANNOTATOR\s*=\s*".*?"', expected, text)
  strus_path.write_text(text, encoding="utf-8")
PY

chmod +x "$SCRIPT_DIR/StruS" || true

echo
echo "Installation finished."
echo "Use: $SCRIPT_DIR/StruS RTBS target.pdb prediction.pdb"
