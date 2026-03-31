"""
Shared harmonization utilities.
"""

import logging
from typing import Any, Optional

import etl._bootstrap  # noqa: F401
from config.settings import (
    DAYS_PER_MONTH, MENOPAUSE_AGE_CUTOFF,
    PAM50_MAP, RECEPTOR_STATUS_MAP, STAGE_MAP,
    HISTOLOGY_MAP, VARIANT_CLASS_MAP,
)

logger = logging.getLogger(__name__)

NULL_SENTINELS = frozenset([
    "", "NA", "N/A", "na", "n/a", "[Not Available]", "[Not Evaluated]",
    "[Not Applicable]", "[Discrepancy]", "[Unknown]", "null", "NULL",
    "NaN", "nan", "--", ".", "unknown", "Unknown",
])


def clean_value(val) -> Optional[str]:
    if val is None:
        return None
    s = str(val).strip()
    if s in NULL_SENTINELS:
        return None
    return s


def clean_float(val) -> Optional[float]:
    s = clean_value(val)
    if s is None:
        return None
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def clean_int(val) -> Optional[int]:
    f = clean_float(val)
    if f is None:
        return None
    return int(f)


def days_to_months(days) -> Optional[float]:
    d = clean_float(days)
    if d is None:
        return None
    return round(d / DAYS_PER_MONTH, 2)


def normalize_pam50(raw: str) -> str:
    s = clean_value(raw)
    if s is None:
        return "Unknown"
    return PAM50_MAP.get(s, PAM50_MAP.get(s.strip(), "Unknown"))


def normalize_receptor_status(raw: str) -> str:
    s = clean_value(raw)
    if s is None:
        return "Unknown"
    return RECEPTOR_STATUS_MAP.get(s, RECEPTOR_STATUS_MAP.get(s.title(), "Unknown"))


def derive_tnbc(er: str, pr: str, her2: str) -> Optional[int]:
    statuses = [er, pr, her2]
    if any(s == "Unknown" for s in statuses):
        return None
    if all(s == "Negative" for s in statuses):
        return 1
    return 0


def normalize_stage(raw: str) -> Optional[str]:
    s = clean_value(raw)
    if s is None:
        return None
    return STAGE_MAP.get(s, STAGE_MAP.get(s.upper(), None))


def normalize_histology(raw: str) -> Optional[str]:
    s = clean_value(raw)
    if s is None:
        return None
    return HISTOLOGY_MAP.get(s, s)


def infer_menopausal_status(age: float | None) -> str:
    a = clean_float(age)
    if a is None:
        return "Unknown"
    return "Pre" if a < MENOPAUSE_AGE_CUTOFF else "Post"


def normalize_menopause_status(raw: str) -> str | None:
    """Map verbose cBioPortal MENOPAUSE_STATUS strings to Pre/Post/Peri/Unknown."""
    s = clean_value(raw)
    if s is None:
        return None
    key = s.lower()
    if key.startswith("pre"):
        return "Pre"
    if key.startswith("post"):
        return "Post"
    if key.startswith("peri"):
        return "Peri"
    if key.startswith("indeterminate"):
        return "Unknown"
    return "Unknown"


def normalize_vital_status(raw: str) -> Optional[str]:
    s = clean_value(raw)
    if s is None:
        return None
    mapping = {
        "Alive": "Alive", "alive": "Alive", "LIVING": "Alive", "Living": "Alive", "0": "Alive",
        "Dead": "Deceased", "dead": "Deceased", "Deceased": "Deceased", "DECEASED": "Deceased", "1": "Deceased",
        "Died of Disease": "Deceased", "Died of Other Causes": "Deceased",
    }
    return mapping.get(s, None)


def normalize_variant_classification(raw: str) -> Optional[str]:
    s = clean_value(raw)
    if s is None:
        return None
    return VARIANT_CLASS_MAP.get(s, s)


def compute_npi(tumor_size_cm: float | None, tumor_grade: int | None,
                lymph_nodes_positive: int | None) -> Optional[float]:
    size = clean_float(tumor_size_cm)
    grade = clean_int(tumor_grade)
    nodes = clean_int(lymph_nodes_positive)
    if size is None or grade is None or nodes is None:
        return None
    if nodes == 0:
        ln_score = 1
    elif nodes <= 3:
        ln_score = 2
    else:
        ln_score = 3
    return round(0.2 * size + grade + ln_score, 2)


def parse_tcga_barcode(barcode: str) -> dict[str, Any]:
    parts = barcode.split("-")
    result = {
        "barcode": barcode,
        "patient_id": "-".join(parts[:3]) if len(parts) >= 3 else barcode,
    }
    if len(parts) >= 4:
        sample_code = parts[3][:2]
        sample_type_map = {
            "01": "Primary Solid Tumor", "02": "Recurrent Solid Tumor",
            "06": "Metastatic", "10": "Blood Derived Normal", "11": "Solid Tissue Normal",
        }
        result["sample_type_code"] = sample_code
        result["sample_type"] = sample_type_map.get(sample_code, f"Unknown ({sample_code})")
        result["tumor_id"] = "-".join(parts[:4])
    else:
        result["tumor_id"] = result["patient_id"]
    return result
