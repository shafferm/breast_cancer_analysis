"""
Load gene expression from TCGA-BRCA (RNA-seq) and METABRIC (microarray).

Usage:  python -m etl.load_expression
"""

import logging
from pathlib import Path
from typing import Any

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_EXPRESSION, METABRIC_EXPRESSION, INSERT_BATCH_SIZE
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
    logger.info(f"  [{study}] {df.shape[0]} genes x {df.shape[1]-1} samples")
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else df.columns[0]
    if df[gene_col].is_duplicated().any():
        df = df.unique(subset=[gene_col], keep="first", maintain_order=True)
    return df


def _melt_insert(conn: Connection, df: pl.DataFrame, study: str,
                 platform: str, existing: set[str]) -> None:
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else df.columns[0]
    entrez_col = "Entrez_Gene_Id" if "Entrez_Gene_Id" in df.columns else None
    skip_cols = {gene_col, entrez_col, "Entrez_Gene_Id"} - {None}
    sample_cols = [c for c in df.columns if c not in skip_cols]

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

    # Cast matched columns to float and precompute per-gene mean and std
    df = df.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in matched])
    df = df.with_columns([
        pl.mean_horizontal(matched).alias("__mu"),
        pl.std_horizontal(matched).alias("__sigma"),
    ])

    sql = text("INSERT OR IGNORE INTO gene_expression (tumor_id, gene_symbol, entrez_id, expression_value, expression_z_score, platform) VALUES (:tumor_id, :gene_symbol, :entrez_id, :expression_value, :expression_z_score, :platform)")
    total, batch = 0, []

    for i, row in enumerate(df.iter_rows(named=True)):
        gene = str(row[gene_col]).strip()
        if not gene or gene in ("nan", "?"):
            continue
        entrez = None
        if entrez_col:
            try:
                v = row.get(entrez_col)
                if v is not None:
                    entrez = int(float(v))
            except (ValueError, TypeError):
                pass
        mu = row["__mu"]
        sigma = row["__sigma"]
        has_var = sigma is not None and sigma > 0
        for c in matched:
            v = row[c]
            if v is None:
                continue
            z = round((v - mu) / sigma, 4) if has_var else 0.0
            batch.append({"tumor_id": col_map[c], "gene_symbol": gene, "entrez_id": entrez,
                          "expression_value": v, "expression_z_score": z, "platform": platform})
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
    logger.info(f"  [{study}] {total:,} expression rows inserted")


def load_expression() -> None:
    with get_connection() as conn:
        existing = _tumor_ids(conn)
        logger.info(f"{len(existing)} tumors in DB")
        for path, study, plat in [(TCGA_EXPRESSION, "TCGA-BRCA", "RNA-seq"),
                                   (METABRIC_EXPRESSION, "METABRIC", "microarray")]:
            if path.exists():
                df = _load_matrix(path, study)
                if not df.is_empty():
                    _melt_insert(conn, df, study, plat, existing)
                    log_harmonization(conn, "gene_expression", study, "expression_z_score",
                                      f"Z-score per gene across all {study} samples", "expression_value",
                                      "Enables rank-based cross-study comparison")
            else:
                logger.warning(f"{study} expression not found: {path}")
        n = conn.execute(text("SELECT count(*) FROM gene_expression")).scalar()
        logger.info(f"Expression load complete: {n:,} rows")


if __name__ == "__main__":
    load_expression()
