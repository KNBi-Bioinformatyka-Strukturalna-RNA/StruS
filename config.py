import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_BIN = sys.executable
ANNOTATOR = str(SCRIPT_DIR / "rnapolis-py" / "src" / "rnapolis" / "annotator.py")
CONVERTER = str(SCRIPT_DIR / "annotation_converter.py")
RTBS = str(SCRIPT_DIR / "RTBS.py")
MBR = str(SCRIPT_DIR / "mbr_matrix.json")