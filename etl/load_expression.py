"""
Load raw gene expression into the gene_expression table (OPTIONAL).

This module is NOT part of the default pipeline (run_all.py). ComBat batch
correction reads expression data directly from the flat files instead.

Run this standalone if you need raw expression values queryable in the DB:

    python -m etl.load_expression
    python -m etl.load_expression --all-genes

Only genes present in BOTH platforms are loaded by default. Pass --all-genes
to load every gene from each study.
"""

import argparse
import logging
from pathlib import Path

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_EXPRESSION, METABRIC_EXPRESSION
from etl.db import get_connection, log_harmonization, drop_indexes, recreate_indexes

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

BULK_BATCH_SIZE = 50_000


def _tumor_ids(conn: Connection) -> set[str]:
    return {r[0] for r in conn.execute(text("SELECT tumor_id FROM tumors")).fetchall()}


def _gene_col(df: pl.DataFrame) -> str:
    return "Hugo_Symbol" if "Hugo_Symbol" in df.columns else df.columns[0]


def _load_matrix(path: Path, study: str) -> pl.DataFrame:
    if not path.exists():
        logger.warning(f"  [{study}] Not found: {path}")
        return pl.DataFrame()
    df = pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)
    logger.info(f"  [{study}] {df.shape[0]:,} genes x {df.shape[1]-1:,} samples")
    gc = _gene_col(df)
    if df[gc].is_duplicated().any():
        df = df.unique(subset=[gc], keep="first", maintain_order=True)
    return df


def _shared_genes(tcga_df: pl.DataFrame, metabric_df: pl.DataFrame) -> set[str]:
    """Return gene symbols present in both platforms. Empty set = no filtering."""
    if tcga_df.is_empty() or metabric_df.is_empty():
        return set()
    tcga_genes = set(tcga_df[_gene_col(tcga_df)].to_list())
    met_genes = set(metabric_df[_gene_col(metabric_df)].to_list())
    return tcga_genes & met_genes


def _melt_insert(conn: Connection, df: pl.DataFrame, study: str,
                 platform: str, existing: set[str], shared_genes: set[str]) -> None:
    gene_col = _gene_col(df)
    entrez_col = "Entrez_Gene_Id" if "Entrez_Gene_Id" in df.columns else None
    skip_cols = {gene_col, entrez_col, "Entrez_Gene_Id"} - {None}
    sample_cols = [c for c in df.columns if c not in skip_cols]

    # Build column-name -> tumor_id mapping
    if study == "TCGA-BRCA":
        col_map = {}
        for c in sample_cols:
            tid = "-".join(c.split("-")[:3]) + "-01"
            if tid in existing:
                col_map[c] = tid
    else:
        col_map = {c: c for c in sample_cols if c in existing}

    logger.info(f"  [{study}] Matched {len(col_map):,} / {len(sample_cols):,} samples")
    if not col_map:
        return

    matched = list(col_map.keys())

    # Filter to shared genes (skip invalid gene symbols)
    df = df.filter(
        pl.col(gene_col).is_not_null()
        & ~pl.col(gene_col).is_in(["", "nan", "?"])
    )
    if shared_genes:
        before = len(df)
        df = df.filter(pl.col(gene_col).is_in(list(shared_genes)))
        logger.info(f"  [{study}] Filtered to {len(df):,} / {before:,} shared genes")

    # Cast matched columns to float
    n = len(matched)
    df = df.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in matched])

    # Compute per-gene mean and std (vectorized)
    df = df.with_columns(pl.mean_horizontal(matched).alias("__mu"))
    if n > 1:
        df = df.with_columns(
            (pl.mean_horizontal([(pl.col(c) - pl.col("__mu")) ** 2 for c in matched]) * n / (n - 1))
            .sqrt()
            .alias("__sigma")
        )
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("__sigma"))

    # Resolve entrez IDs before melt
    if entrez_col:
        df = df.with_columns(
            pl.col(entrez_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).alias("__entrez")
        )
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Int64).alias("__entrez"))

    # Unpivot (melt) from wide to long
    index_cols = [gene_col, "__entrez", "__mu", "__sigma"]
    long = df.unpivot(
        on=matched,
        index=index_cols,
        variable_name="__sample",
        value_name="expression_value",
    )

    # Drop nulls and map sample columns to tumor_ids
    long = long.drop_nulls(subset=["expression_value"])
    mapping_df = pl.DataFrame({"__sample": list(col_map.keys()), "tumor_id": list(col_map.values())})
    long = long.join(mapping_df, on="__sample", how="inner")

    # Compute z-score vectorized
    long = long.with_columns(
        pl.when((pl.col("__sigma").is_not_null()) & (pl.col("__sigma") > 0))
        .then(((pl.col("expression_value") - pl.col("__mu")) / pl.col("__sigma")).round(4))
        .otherwise(0.0)
        .alias("expression_z_score")
    )

    # Select final columns and add platform
    long = long.select([
        "tumor_id",
        pl.col(gene_col).alias("gene_symbol"),
        pl.col("__entrez").alias("entrez_id"),
        "expression_value",
        "expression_z_score",
        pl.lit(platform).alias("platform"),
    ])

    # Bulk insert via SQLAlchemy
    sql = text(
        "INSERT OR IGNORE INTO gene_expression "
        "(tumor_id, gene_symbol, entrez_id, expression_value, expression_z_score, platform) "
        "VALUES (:tumor_id, :gene_symbol, :entrez_id, :expression_value, :expression_z_score, :platform)"
    )
    rows = long.to_dicts()
    saved_indexes = drop_indexes(conn, "gene_expression")
    conn.execute(text("PRAGMA foreign_keys = OFF"))
    total = 0
    for i in range(0, len(rows), BULK_BATCH_SIZE):
        batch = rows[i:i + BULK_BATCH_SIZE]
        conn.execute(sql, batch)
        total += len(batch)
        if total % (BULK_BATCH_SIZE * 5) == 0:
            logger.info(f"  [{study}] {total:,} / {len(rows):,} rows inserted")
    conn.commit()
    conn.execute(text("PRAGMA foreign_keys = ON"))
    recreate_indexes(conn, saved_indexes)
    logger.info(f"  [{study}] {total:,} expression rows inserted")


def load_expression(all_genes: bool = False) -> None:
    with get_connection() as conn:
        existing = _tumor_ids(conn)
        logger.info(f"{len(existing):,} tumors in DB")

        # Load both matrices to find shared genes
        tcga_df = _load_matrix(TCGA_EXPRESSION, "TCGA-BRCA") if TCGA_EXPRESSION.exists() else pl.DataFrame()
        metabric_df = _load_matrix(METABRIC_EXPRESSION, "METABRIC") if METABRIC_EXPRESSION.exists() else pl.DataFrame()

        if all_genes:
            shared = set()
            logger.info("--all-genes: loading all genes from each study")
        else:
            shared = _shared_genes(tcga_df, metabric_df)
            if shared:
                logger.info(f"{len(shared):,} genes shared between platforms")
            else:
                logger.info("Only one platform available — loading all genes")

        for df, study, plat in [(tcga_df, "TCGA-BRCA", "RNA-seq"),
                                 (metabric_df, "METABRIC", "microarray")]:
            if not df.is_empty():
                _melt_insert(conn, df, study, plat, existing, shared)
                log_harmonization(conn, "gene_expression", study, "expression_z_score",
                                  f"Z-score per gene across all {study} samples", "expression_value",
                                  "Enables rank-based cross-study comparison")
            else:
                logger.warning(f"  [{study}] Expression not found")

        n = conn.execute(text("SELECT count(*) FROM gene_expression")).scalar()
        logger.info(f"Expression load complete: {n:,} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load raw gene expression data (optional)")
    parser.add_argument("--all-genes", action="store_true",
                        help="Load all genes, not just those shared between studies")
    args = parser.parse_args()
    load_expression(all_genes=args.all_genes)
