"""
Load copy number (GISTIC2) from TCGA-BRCA and METABRIC.
Both used Affymetrix SNP6.0 - most directly comparable data type.

By default, only genes present in BOTH datasets are loaded. Pass --all-genes
to load every gene from each study.

Usage:
    python -m etl.load_copy_number
    python -m etl.load_copy_number --all-genes
"""

import argparse
import logging
from pathlib import Path

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_CNA, METABRIC_CNA
from etl.db import get_connection, log_harmonization, bulk_executemany, drop_indexes, recreate_indexes

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

VALID_GISTIC = {-2, -1, 0, 1, 2}


def _tumor_ids(conn: Connection) -> set[str]:
    return {r[0] for r in conn.execute(text("SELECT tumor_id FROM tumors")).fetchall()}


def _load_matrix(path: Path, study: str) -> pl.DataFrame:
    if not path.exists():
        logger.warning(f"  [{study}] Not found: {path}")
        return pl.DataFrame()
    df = pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(["Entrez_Gene_Id"])
    gene_col = df.columns[0]
    if df[gene_col].is_duplicated().any():
        df = df.unique(subset=[gene_col], keep="first", maintain_order=True)
    logger.info(f"  [{study}] {df.shape[0]:,} genes x {df.shape[1]-1:,} samples")
    return df


def _shared_genes(tcga_df: pl.DataFrame, metabric_df: pl.DataFrame) -> set[str]:
    """Return gene symbols present in both CNA files. Empty set = no filtering."""
    if tcga_df.is_empty() or metabric_df.is_empty():
        return set()
    tcga_genes = set(tcga_df[tcga_df.columns[0]].to_list())
    met_genes = set(metabric_df[metabric_df.columns[0]].to_list())
    return tcga_genes & met_genes


def _melt_insert(conn: Connection, df: pl.DataFrame, study: str,
                 existing: set[str], shared_genes: set[str]) -> None:
    gene_col = df.columns[0]
    sample_cols = df.columns[1:]

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

    # Filter invalid genes
    df = df.filter(
        pl.col(gene_col).is_not_null()
        & ~pl.col(gene_col).is_in(["", "nan", "?"])
    )

    # Filter to shared genes
    if shared_genes:
        before = len(df)
        df = df.filter(pl.col(gene_col).is_in(list(shared_genes)))
        logger.info(f"  [{study}] Filtered to {len(df):,} / {before:,} shared genes")

    # Cast to Int64 and unpivot
    df = df.with_columns([pl.col(c).cast(pl.Int64, strict=False) for c in matched])
    long = df.unpivot(
        on=matched,
        index=[gene_col],
        variable_name="__sample",
        value_name="gistic_value",
    )

    # Drop nulls and filter to valid GISTIC values
    long = long.drop_nulls(subset=["gistic_value"])
    long = long.filter(pl.col("gistic_value").is_in(list(VALID_GISTIC)))

    # Map sample columns to tumor_ids
    mapping_df = pl.DataFrame({"__sample": list(col_map.keys()), "tumor_id": list(col_map.values())})
    long = long.join(mapping_df, on="__sample", how="inner")

    nonzero = long.filter(pl.col("gistic_value") != 0).height

    # Build tuples: (tumor_id, gene_symbol, gistic_value, log2_ratio=None)
    tuples_df = long.select("tumor_id", pl.col(gene_col), "gistic_value")
    rows = [
        (t, g, v, None)
        for t, g, v in zip(
            tuples_df["tumor_id"].to_list(),
            tuples_df[gene_col].to_list(),
            tuples_df["gistic_value"].to_list(),
        )
    ]

    logger.info(f"  [{study}] Inserting {len(rows):,} rows ({nonzero:,} non-diploid)...")
    saved_indexes = drop_indexes(conn, "copy_number")
    conn.execute(text("PRAGMA foreign_keys = OFF"))
    sql = "INSERT OR IGNORE INTO copy_number (tumor_id, gene_symbol, gistic_value, log2_ratio) VALUES (?, ?, ?, ?)"
    total = bulk_executemany(conn, sql, rows, label=f"[{study}] ")
    conn.commit()
    conn.execute(text("PRAGMA foreign_keys = ON"))
    recreate_indexes(conn, saved_indexes)
    logger.info(f"  [{study}] {total:,} CNA rows inserted")


def load_copy_number(all_genes: bool = False) -> None:
    with get_connection() as conn:
        existing = _tumor_ids(conn)
        logger.info(f"{len(existing):,} tumors in DB")

        # Load both matrices to find shared genes
        tcga_df = _load_matrix(TCGA_CNA, "TCGA-BRCA") if TCGA_CNA.exists() else pl.DataFrame()
        metabric_df = _load_matrix(METABRIC_CNA, "METABRIC") if METABRIC_CNA.exists() else pl.DataFrame()

        if all_genes:
            shared = set()
            logger.info("--all-genes: loading all genes from each study")
        else:
            shared = _shared_genes(tcga_df, metabric_df)
            if shared:
                logger.info(f"{len(shared):,} genes shared between CNA platforms")
            else:
                logger.info("Only one CNA dataset available — loading all genes")

        for df, study in [(tcga_df, "TCGA-BRCA"), (metabric_df, "METABRIC")]:
            if not df.is_empty():
                _melt_insert(conn, df, study, existing, shared)
                log_harmonization(conn, "copy_number", study, "gistic_value",
                                  "GISTIC2 discrete values from Affymetrix SNP6.0",
                                  "data_cna.txt", "Directly comparable across studies")
            else:
                logger.warning(f"  [{study}] CNA not found")

        n = conn.execute(text("SELECT count(*) FROM copy_number")).scalar()
        logger.info(f"Copy number load complete: {n:,} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load copy number data")
    parser.add_argument("--all-genes", action="store_true",
                        help="Load all genes, not just those shared between studies")
    args = parser.parse_args()
    load_copy_number(all_genes=args.all_genes)
