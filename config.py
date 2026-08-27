import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_BIN = sys.executable
ANNOTATOR = str(SCRIPT_DIR / "rnapolis-py" / "src" / "rnapolis" / "annotator.py")
CONVERTER = str(SCRIPT_DIR / "annotation_converter.py")
RTBS = str(SCRIPT_DIR / "RTBS.py")
MBR = str(SCRIPT_DIR / "mbr_matrix.json")
FR3D = str(SCRIPT_DIR / "fr3d-python" / "fr3d" / "classifiers" / "NA_pairwise_interactions.py")
FR3D_TO_DBN = str(SCRIPT_DIR / "fr3d_to_dbn.py")
ANNOTATOR_TO_DBN = str(SCRIPT_DIR / "annotator_to_dbn.py")
CONVERTER_FROM_DBN = str(SCRIPT_DIR / "converter_from_dbn.py")
MOLECULE_FILTER = str(SCRIPT_DIR / "rnapolis-py" / "src" / "rnapolis" / "molecule_filter.py")
RNAPOLIS_SRC = str(SCRIPT_DIR / "rnapolis-py" / "src")