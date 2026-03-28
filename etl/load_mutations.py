"""
Load somatic mutations from TCGA-BRCA (WES) and METABRIC (panel).

Usage:  python -m etl.load_mutations
"""

import logging
import uuid
from pathlib import Path
from typing import Any

import polars as pl
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_MUTATIONS_MAF, METABRIC_MUTATIONS, TARGET_GENOME_BUILD
from config.metabric_panel_genes import METABRIC_PANEL_GENES
from etl.db import get_connection, bulk_insert, log_harmonization
from etl.harmonize import clean_value, clean_float, clean_int, normalize_variant_classification

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _read_maf(path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)
    logger.info(f"  Read {len(df)} rows from {path.name}")
    return df


def _tumor_ids(conn: Connection) -> set[str]:
    return {r[0] for r in conn.execute(text("SELECT tumor_id FROM tumors")).fetchall()}


def _std_cols(df: pl.DataFrame) -> pl.DataFrame:
    m = {"Hugo_Symbol": "gene", "Chromosome": "chr", "Start_Position": "start",
         "End_Position": "end", "Reference_Allele": "ref", "Tumor_Seq_Allele2": "alt",
         "Variant_Classification": "vclass", "Variant_Type": "vtype",
         "HGVSp_Short": "hgvsp", "HGVSp": "hgvsp_long",
         "Tumor_Sample_Barcode": "barcode", "SAMPLE_ID": "barcode",
         "t_alt_count": "t_alt", "t_ref_count": "t_ref", "t_depth": "t_depth"}
    return df.rename({k: v for k, v in m.items() if k in df.columns})


VALID_VTYPES = {"SNP", "INS", "DEL", "DNP", "TNP", "ONP"}


def load_tcga_mutations(conn: Connection) -> list[dict[str, Any]]:
    logger.info("Loading TCGA-BRCA mutations...")
    if not TCGA_MUTATIONS_MAF.exists():
        logger.warning(f"Not found: {TCGA_MUTATIONS_MAF}")
        return []
    df = _std_cols(_read_maf(TCGA_MUTATIONS_MAF))
    if "t_alt" in df.columns and "t_depth" in df.columns:
        df = df.with_columns(
            (pl.col("t_alt").cast(pl.Float64, strict=False) /
             pl.col("t_depth").cast(pl.Float64, strict=False)).alias("vaf")
        )
    elif "t_alt" in df.columns and "t_ref" in df.columns:
        df = df.with_columns(
            (pl.col("t_alt").cast(pl.Float64, strict=False) /
             (pl.col("t_alt").cast(pl.Float64, strict=False) +
              pl.col("t_ref").cast(pl.Float64, strict=False))).alias("vaf")
        )
    else:
        df = df.with_columns(pl.lit(None).alias("vaf"))

    tumors = _tumor_ids(conn)
    rows, skip = [], 0
    for r in df.to_dicts():
        gene = clean_value(r.get("gene"))
        if not gene:
            continue
        pid = "-".join(str(r.get("barcode", "")).split("-")[:3])
        tid = f"{pid}-01"
        if tid not in tumors:
            skip += 1
            continue
        vt = clean_value(r.get("vtype"))
        rows.append(dict(
            mutation_id=f"TCGA-MUT-{uuid.uuid4().hex[:12]}", tumor_id=tid, gene_symbol=gene,
            chromosome=clean_value(r.get("chr")),
            start_position=clean_int(r.get("start")), end_position=clean_int(r.get("end")),
            reference_genome=TARGET_GENOME_BUILD,
            reference_allele=clean_value(r.get("ref")), variant_allele=clean_value(r.get("alt")),
            variant_classification=normalize_variant_classification(r.get("vclass", "")),
            variant_type=vt if vt in VALID_VTYPES else None,
            hgvsp=clean_value(r.get("hgvsp")) or clean_value(r.get("hgvsp_long")),
            vaf=clean_float(r.get("vaf")), sequencing_scope="WES",
            in_metabric_panel=1 if gene in METABRIC_PANEL_GENES else 0,
        ))
    if skip:
        logger.warning(f"  Skipped {skip} (no matching tumor)")
    panel_n = sum(1 for x in rows if x["in_metabric_panel"])
    logger.info(f"  TCGA: {len(rows)} mutations ({panel_n} in METABRIC panel)")
    return rows


def load_metabric_mutations(conn: Connection) -> list[dict[str, Any]]:
    logger.info("Loading METABRIC mutations...")
    if not METABRIC_MUTATIONS.exists():
        logger.warning(f"Not found: {METABRIC_MUTATIONS}")
        return []
    df = _std_cols(_read_maf(METABRIC_MUTATIONS))
    tumors = _tumor_ids(conn)
    rows, skip = [], 0
    for r in df.to_dicts():
        gene = clean_value(r.get("gene"))
        if not gene:
            continue
        tid = clean_value(r.get("barcode"))
        if not tid or tid not in tumors:
            skip += 1
            continue
        vt = clean_value(r.get("vtype"))
        rows.append(dict(
            mutation_id=f"MB-MUT-{uuid.uuid4().hex[:12]}", tumor_id=tid, gene_symbol=gene,
            chromosome=clean_value(r.get("chr")),
            start_position=clean_int(r.get("start")), end_position=clean_int(r.get("end")),
            reference_genome="GRCh37",
            reference_allele=clean_value(r.get("ref")), variant_allele=clean_value(r.get("alt")),
            variant_classification=normalize_variant_classification(r.get("vclass", "")),
            variant_type=vt if vt in VALID_VTYPES else None,
            hgvsp=clean_value(r.get("hgvsp")) or clean_value(r.get("hgvsp_long")),
            vaf=None, sequencing_scope="PANEL", in_metabric_panel=1,
        ))
    if skip:
        logger.warning(f"  Skipped {skip} (no matching tumor)")
    logger.info(f"  METABRIC: {len(rows)} mutations")
    return rows


def load_mutations() -> None:
    with get_connection() as conn:
        for loader in [load_tcga_mutations, load_metabric_mutations]:
            muts = loader(conn)
            if muts:
                bulk_insert(conn, "mutations", muts)
        log_harmonization(conn, "mutations", "BOTH", "in_metabric_panel",
                          f"TCGA mutations flagged against {len(METABRIC_PANEL_GENES)}-gene panel list",
                          "Hugo_Symbol", "Use WHERE in_metabric_panel=1 for cross-study comparisons")
        log_harmonization(conn, "mutations", "METABRIC", "reference_genome",
                          "Stored as GRCh37. Liftover to GRCh38 needed for coordinate comparisons.",
                          "Start_Position", "Gene-level analysis works without liftover")
        n = conn.execute(text("SELECT count(*) FROM mutations")).scalar()
        logger.info(f"Mutation load complete: {n} total")


if __name__ == "__main__":
    load_mutations()
