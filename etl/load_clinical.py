"""
Load clinical data (patients, tumors, treatments) from TCGA-BRCA and METABRIC.

Usage:  python -m etl.load_clinical
"""

import logging
from pathlib import Path
from typing import Any

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import (
    TCGA_CDR_FILE, TCGA_CLINICAL_PATIENT, TCGA_CLINICAL_DRUG,
    METABRIC_CLINICAL_PATIENT, METABRIC_CLINICAL_SAMPLE,
)
from etl.db import get_connection, bulk_insert, log_harmonization
from etl.harmonize import (
    clean_value, clean_float, clean_int, days_to_months,
    normalize_pam50, normalize_receptor_status, normalize_vital_status,
    normalize_stage, normalize_histology, derive_tnbc,
    infer_menopausal_status, compute_npi,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def read_cbioportal_file(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)


# ---- TCGA ----

def load_tcga_clinical_cdr(conn: Connection) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    logger.info("Loading TCGA-BRCA clinical from CDR file...")
    if not TCGA_CDR_FILE.exists():
        logger.warning(f"CDR not found: {TCGA_CDR_FILE}; trying biotab")
        return load_tcga_clinical_biotab(conn)

    df = pl.read_excel(TCGA_CDR_FILE, sheet_name="TCGA-CDR")
    df = df.filter(pl.col("type") == "BRCA")
    logger.info(f"  [TCGA-BRCA] {len(df):,} cases in CDR")

    patients, tumors = [], []
    for r in df.to_dicts():
        pid = clean_value(r.get("bcr_patient_barcode"))
        if not pid:
            continue
        age = clean_float(r.get("age_at_initial_pathologic_diagnosis"))
        os_ev = clean_int(r.get("OS"))
        vital = {1: "Deceased", 0: "Alive"}.get(os_ev)
        er = normalize_receptor_status(r.get("er_status_by_ihc", ""))
        pr = normalize_receptor_status(r.get("pr_status_by_ihc", ""))
        her2 = normalize_receptor_status(r.get("her2_status_by_ihc", ""))
        patients.append(dict(
            patient_id=pid, study_source="TCGA-BRCA", age_at_diagnosis=age, sex="F",
            race=clean_value(r.get("race")), ethnicity=clean_value(r.get("ethnicity")),
            menopausal_status=infer_menopausal_status(age), vital_status=vital,
            os_months=days_to_months(r.get("OS.time")), os_event=os_ev,
            dfs_months=days_to_months(r.get("PFI.time")), dfs_event=clean_int(r.get("PFI")),
        ))
        tumors.append(dict(
            tumor_id=f"{pid}-01", patient_id=pid,
            histological_type=normalize_histology(r.get("histological_type", "")),
            tumor_grade=clean_int(r.get("neoplasm_histologic_grade")),
            tumor_stage=normalize_stage(r.get("ajcc_pathologic_tumor_stage", "")),
            tumor_size_cm=None, lymph_nodes_positive=None, lymph_nodes_examined=None,
            er_status=er, pr_status=pr, her2_status=her2, tnbc_status=derive_tnbc(er, pr, her2),
            pam50_subtype=normalize_pam50(r.get("paper_BRCA_Subtype_PAM50", "")),
            integrative_cluster=None, nottingham_prognostic_index=None,
        ))
    tx = _load_tcga_treatments({p["patient_id"] for p in patients})
    logger.info(f"  [TCGA-BRCA] {len(patients):,} patients, {len(tumors):,} tumors, {len(tx):,} treatments")
    return patients, tumors, tx


def load_tcga_clinical_biotab(conn: Connection) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    logger.info("Loading TCGA-BRCA from biotab files...")
    if not TCGA_CLINICAL_PATIENT.exists():
        logger.error(f"Not found: {TCGA_CLINICAL_PATIENT}")
        return [], [], []
    df = pl.read_csv(TCGA_CLINICAL_PATIENT, separator="\t", skip_rows_after_header=2,
                     infer_schema_length=0)
    patients, tumors = [], []
    for r in df.to_dicts():
        pid = clean_value(r.get("bcr_patient_barcode"))
        if not pid:
            continue
        age = clean_float(r.get("age_at_initial_pathologic_diagnosis"))
        vital = normalize_vital_status(r.get("vital_status", ""))
        os_ev = 1 if vital == "Deceased" else 0 if vital == "Alive" else None
        dtd, dtf = clean_float(r.get("death_days_to")), clean_float(r.get("last_contact_days_to"))
        os_m = days_to_months(dtd) if os_ev == 1 and dtd else days_to_months(dtf)
        er = normalize_receptor_status(r.get("er_status_by_ihc", ""))
        pr = normalize_receptor_status(r.get("pr_status_by_ihc", ""))
        her2 = normalize_receptor_status(r.get("her2_status_by_ihc", ""))
        sz = clean_float(r.get("tumor_size"))
        gr = clean_int(r.get("neoplasm_histologic_grade"))
        ln = clean_int(r.get("lymph_node_examined_count"))
        patients.append(dict(
            patient_id=pid, study_source="TCGA-BRCA", age_at_diagnosis=age,
            sex=clean_value(r.get("gender", "F")),
            race=clean_value(r.get("race_list") or r.get("race")),
            ethnicity=clean_value(r.get("ethnicity")),
            menopausal_status=infer_menopausal_status(age), vital_status=vital,
            os_months=os_m, os_event=os_ev, dfs_months=None, dfs_event=None,
        ))
        tumors.append(dict(
            tumor_id=f"{pid}-01", patient_id=pid,
            histological_type=normalize_histology(r.get("histological_type", "")),
            tumor_grade=gr, tumor_stage=normalize_stage(r.get("ajcc_pathologic_tumor_stage", "")),
            tumor_size_cm=sz, lymph_nodes_positive=ln,
            lymph_nodes_examined=clean_int(r.get("lymph_nodes_examined_count")),
            er_status=er, pr_status=pr, her2_status=her2, tnbc_status=derive_tnbc(er, pr, her2),
            pam50_subtype="Unknown", integrative_cluster=None,
            nottingham_prognostic_index=compute_npi(sz, gr, ln),
        ))
    tx = _load_tcga_treatments({p["patient_id"] for p in patients})
    logger.info(f"  [TCGA-BRCA] {len(patients):,} patients, {len(tumors):,} tumors, {len(tx):,} treatments")
    return patients, tumors, tx


def _load_tcga_treatments(pids: set[str]) -> list[dict[str, Any]]:
    if not TCGA_CLINICAL_DRUG.exists():
        return []
    df = pl.read_csv(TCGA_CLINICAL_DRUG, separator="\t", skip_rows_after_header=2,
                     infer_schema_length=0)
    tx_map = {"Chemotherapy": "Chemotherapy", "Hormone Therapy": "Hormone Therapy",
              "Targeted Molecular therapy": "Targeted Therapy",
              "Immunotherapy": "Other", "Ancillary": "Other"}
    out, ctr = [], 0
    for r in df.to_dicts():
        pid = clean_value(r.get("bcr_patient_barcode"))
        if not pid or pid not in pids:
            continue
        ctr += 1
        out.append(dict(
            treatment_id=f"TCGA-TX-{ctr:06d}", patient_id=pid,
            treatment_type=tx_map.get(clean_value(r.get("pharmaceutical_therapy_type", "")), "Other"),
            treatment_name=clean_value(r.get("pharmaceutical_therapy_drug_name")), received=1,
        ))
    return out


# ---- METABRIC ----

def load_metabric_clinical(conn: Connection) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    logger.info("Loading METABRIC clinical...")
    if not METABRIC_CLINICAL_PATIENT.exists():
        logger.error(f"Not found: {METABRIC_CLINICAL_PATIENT}")
        return [], [], []
    df_p = read_cbioportal_file(METABRIC_CLINICAL_PATIENT)
    if METABRIC_CLINICAL_SAMPLE.exists():
        df_s = read_cbioportal_file(METABRIC_CLINICAL_SAMPLE)
        df = df_p.join(df_s, on="PATIENT_ID", how="left", suffix="_s")
    else:
        df = df_p
    logger.info(f"  [METABRIC] {len(df):,} rows")

    patients, tumors, treatments = [], [], []
    txc = 0
    for r in df.to_dicts():
        pid = clean_value(r.get("PATIENT_ID"))
        if not pid:
            continue
        age = clean_float(r.get("AGE_AT_DIAGNOSIS"))
        os_m = clean_float(r.get("OS_MONTHS"))
        os_raw = str(r.get("OS_STATUS") or "")
        if "DECEASED" in os_raw.upper() or os_raw.startswith("1"):
            os_ev, vital = 1, "Deceased"
        elif "LIVING" in os_raw.upper() or os_raw.startswith("0"):
            os_ev, vital = 0, "Alive"
        else:
            os_ev, vital = None, None
        rfs_m = clean_float(r.get("RFS_MONTHS") or r.get("DFS_MONTHS"))
        rfs_raw = str(r.get("RFS_STATUS") or r.get("DFS_STATUS") or "")
        rfs_ev = 1 if "Recurred" in rfs_raw or rfs_raw.startswith("1") else \
                 0 if "Not Recurred" in rfs_raw or rfs_raw.startswith("0") else None

        meno_raw = clean_value(r.get("INFERRED_MENOPAUSAL_STATE", ""))
        meno = {"Pre": "Pre", "Post": "Post", "pre": "Pre", "post": "Post"}.get(
            meno_raw, "Unknown") if meno_raw else infer_menopausal_status(age)

        patients.append(dict(
            patient_id=pid, study_source="METABRIC", age_at_diagnosis=age, sex="F",
            race=None, ethnicity=None, menopausal_status=meno, vital_status=vital,
            os_months=os_m, os_event=os_ev, dfs_months=rfs_m, dfs_event=rfs_ev,
        ))

        er = normalize_receptor_status(r.get("ER_STATUS") or r.get("ER_IHC") or "")
        pr = normalize_receptor_status(r.get("PR_STATUS", ""))
        her2 = normalize_receptor_status(r.get("HER2_STATUS") or r.get("HER2_SNP6") or "")
        sid = clean_value(r.get("SAMPLE_ID") or f"{pid}-T")
        sz = clean_float(r.get("TUMOR_SIZE"))
        gr = clean_int(r.get("GRADE") or r.get("TUMOR_GRADE"))
        npi = clean_float(r.get("NOTTINGHAM_PROGNOSTIC_INDEX") or r.get("NPI"))
        ln = clean_int(r.get("LYMPH_NODES_EXAMINED_POSITIVE"))

        tumors.append(dict(
            tumor_id=sid, patient_id=pid,
            histological_type=normalize_histology(r.get("TUMOR_TYPE") or r.get("HISTOLOGICAL_SUBTYPE") or ""),
            tumor_grade=gr, tumor_stage=normalize_stage(r.get("TUMOR_STAGE", "")),
            tumor_size_cm=sz, lymph_nodes_positive=ln, lymph_nodes_examined=None,
            er_status=er, pr_status=pr, her2_status=her2, tnbc_status=derive_tnbc(er, pr, her2),
            pam50_subtype=normalize_pam50(r.get("CLAUDIN_SUBTYPE") or r.get("PAM50_SUBTYPE") or ""),
            integrative_cluster=clean_value(r.get("INTCLUST") or r.get("INTEGRATIVE_CLUSTER") or ""),
            nottingham_prognostic_index=npi if npi else compute_npi(sz, gr, ln),
        ))

        for ttype, field in [("Chemotherapy", "CHEMOTHERAPY"), ("Hormone Therapy", "HORMONE_THERAPY"), ("Radiation", "RADIO_THERAPY")]:
            tv = clean_value(r.get(field, ""))
            if tv is not None:
                txc += 1
                treatments.append(dict(
                    treatment_id=f"MB-TX-{txc:06d}", patient_id=pid, treatment_type=ttype,
                    treatment_name=None, received=1 if tv.upper() in ("YES", "1", "Y") else 0,
                ))

    logger.info(f"  [METABRIC] {len(patients):,} patients, {len(tumors):,} tumors, {len(treatments):,} treatments")
    return patients, tumors, treatments


# ---- MAIN ----

def load_clinical() -> None:
    with get_connection() as conn:
        for loader, _study in [(load_tcga_clinical_cdr, "TCGA-BRCA"), (load_metabric_clinical, "METABRIC")]:
            ps, ts, txs = loader(conn)
            if ps:
                bulk_insert(conn, "patients", ps)
                bulk_insert(conn, "tumors", ts)
                if txs:
                    bulk_insert(conn, "treatments", txs)

        log_harmonization(conn, "patients", "TCGA-BRCA", "os_months",
                          "CDR OS.time (days) / 30.44", "OS.time", "PFI recommended for BRCA.")
        log_harmonization(conn, "patients", "METABRIC", "os_months",
                          "Used OS_MONTHS directly", "OS_MONTHS", None)

        for t in ["patients", "tumors", "treatments"]:
            n = conn.execute(text(f"SELECT count(*) FROM {t}")).scalar()
            logger.info(f"  {t}: {n:,} rows")


if __name__ == "__main__":
    load_clinical()
