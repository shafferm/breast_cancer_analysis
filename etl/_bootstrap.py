"""
Bootstrap module: ensures the project root is on sys.path.
Import this at the top of any ETL module before importing config.*.
"""

import sys
from pathlib import Path

_project_root = str(Path(__file__).resolve().parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)
