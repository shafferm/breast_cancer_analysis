"""
Load copy number (GISTIC2) from TCGA-BRCA and METABRIC.
Both used Affymetrix SNP6.0 - most directly comparable data type.

Usage:  python -m etl.load_copy_number
"""

import logging
from pathlib import Path
from typing import Any

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_CNA, METABRIC_CNA, INSERT_BATCH_SIZE
from etl.db import get_connection, log_harmonization

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _tumor_ids(conn: Connection) -> set[str]:
    return {r[0] for r in conn.execute(text("SELECT tumor_id FROM tumors")).fetchall()}


def _load_matrix(path: Path, study: str) -> pl.DataFrame:
    if not path.exists():
        logger.warning(f"Not found: {path}")
        return pl.DataFrame()
    df = pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(["Entrez_Gene_Id"])
    gene_col = df.columns[0]
    if df[gene_col].is_duplicated().any():
        df = df.unique(subset=[gene_col], keep="first", maintain_order=True)
    logger.info(f"  [{study}] {df.shape[0]} genes x {df.shape[1]-1} samples")
    return df


def _melt_insert(conn: Connection, df: pl.DataFrame, study: str,
                 existing: set[str]) -> None:
    gene_col = df.columns[0]
    sample_cols = df.columns[1:]

    if study == "TCGA-BRCA":
        col_map = {}
        for c in sample_cols:
            tid = "-".join(c.split("-")[:3]) + "-01"
            if tid in existing:
                col_map[c] = tid
    else:
        col_map = {c: c for c in sample_cols if c in existing}

    logger.info(f"  [{study}] Matched {len(col_map)}/{len(sample_cols)} samples")
    if not col_map:
        return

    matched = list(col_map.keys())
    sql = text("INSERT OR IGNORE INTO copy_number (tumor_id, gene_symbol, gistic_value, log2_ratio) VALUES (:tumor_id, :gene_symbol, :gistic_value, :log2_ratio)")
    total, nonzero, batch = 0, 0, []

    for i, row in enumerate(df.iter_rows(named=True)):
        gene = str(row[gene_col]).strip()
        if not gene or gene in ("nan", "?"):
            continue
        for c in matched:
            try:
                val = int(float(row[c]))
            except (ValueError, TypeError):
                continue
            if val not in (-2, -1, 0, 1, 2):
                continue
            if val != 0:
                nonzero += 1
            batch.append({"tumor_id": col_map[c], "gene_symbol": gene, "gistic_value": val, "log2_ratio": None})
            if len(batch) >= INSERT_BATCH_SIZE:
                conn.execute(sql, batch)
                total += len(batch)
                batch = []
        if (i + 1) % 2000 == 0:
            logger.info(f"  [{study}] {i+1}/{len(df)} genes, {total:,} rows")

    if batch:
        conn.execute(sql, batch)
        total += len(batch)
    conn.commit()
    logger.info(f"  [{study}] {total:,} CNA rows ({nonzero:,} non-diploid)")


def load_copy_number() -> None:
    with get_connection() as conn:
        existing = _tumor_ids(conn)
        logger.info(f"{len(existing)} tumors in DB")
        for path, study in [(TCGA_CNA, "TCGA-BRCA"), (METABRIC_CNA, "METABRIC")]:
            if path.exists():
                df = _load_matrix(path, study)
                if not df.is_empty():
                    _melt_insert(conn, df, study, existing)
                    log_harmonization(conn, "copy_number", study, "gistic_value",
                                      "GISTIC2 discrete values from Affymetrix SNP6.0",
                                      "data_cna.txt", "Directly comparable across studies")
            else:
                logger.warning(f"{study} CNA not found: {path}")
        n = conn.execute(text("SELECT count(*) FROM copy_number")).scalar()
        logger.info(f"Copy number load complete: {n:,} rows")


if __name__ == "__main__":
    load_copy_number()
