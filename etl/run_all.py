"""
Run the full BRCA harmonized ETL pipeline.

TODO: Replace argparse with Click for consistent CLI across all ETL modules.
TODO: Add a console_scripts entry point (e.g. `brca-etl`) so this can run as a command.
TODO: Add tqdm progress bars to bulk insert loops and long-running pipeline steps.

Usage:
    python -m etl.run_all
    python -m etl.run_all --force
    python -m etl.run_all --skip-download   # skip step 1 (data already present)
"""

import argparse
import logging
import time

from sqlalchemy import text

import etl._bootstrap  # noqa: F401
from config.settings import (
    DB_PATH, TCGA_CDR_FILE, TCGA_CLINICAL_PATIENT, TCGA_MUTATIONS_MAF,
    TCGA_EXPRESSION, TCGA_CNA, METABRIC_CLINICAL_PATIENT,
    METABRIC_CLINICAL_SAMPLE, METABRIC_MUTATIONS, METABRIC_EXPRESSION, METABRIC_CNA,
)
from etl.init_db import init_db
from etl.load_clinical import load_clinical
from etl.load_mutations import load_mutations
from etl.load_copy_number import load_copy_number
from etl.liftover_mutations import liftover_mutations
from etl.combat_expression import run_combat
from etl.db import get_connection
from etl.download_data import download_all

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def print_summary() -> None:
    with get_connection() as conn:
        logger.info("=" * 60)
        logger.info("DATABASE SUMMARY")
        logger.info("=" * 60)
        for t in ["patients", "tumors", "treatments", "mutations",
                   "gene_expression", "gene_expression_combat",
                   "copy_number", "methylation",
                   "protein_expression", "harmonization_log"]:
            n = conn.execute(text(f"SELECT count(*) FROM {t}")).scalar()
            logger.info(f"  {t:30s} {n:>12,} rows")

        for label, sql in [
            ("Patients by study",
             "SELECT study_source, count(*) FROM patients GROUP BY study_source"),
            ("Tumors by PAM50",
             "SELECT pam50_subtype, count(*) FROM tumors GROUP BY pam50_subtype ORDER BY count(*) DESC"),
            ("Mutations by scope",
             "SELECT sequencing_scope, count(*) FROM mutations GROUP BY sequencing_scope"),
            ("Mutations by genome build",
             "SELECT reference_genome, count(*) FROM mutations GROUP BY reference_genome"),
        ]:
            rows = conn.execute(text(sql)).fetchall()
            if rows:
                logger.info(f"\n{label}:")
                for r in rows:
                    logger.info(f"  {str(r[0] or 'NULL'):30s} {r[1]:>12,}")
        logger.info("=" * 60)


def _preflight_check() -> None:
    """Warn about missing data files before loading steps begin."""
    expected = [
        (TCGA_CDR_FILE, "TCGA CDR clinical"),
        (TCGA_CLINICAL_PATIENT, "TCGA clinical patient"),
        (TCGA_MUTATIONS_MAF, "TCGA mutations"),
        (TCGA_EXPRESSION, "TCGA expression"),
        (TCGA_CNA, "TCGA copy number"),
        (METABRIC_CLINICAL_PATIENT, "METABRIC clinical patient"),
        (METABRIC_CLINICAL_SAMPLE, "METABRIC clinical sample"),
        (METABRIC_MUTATIONS, "METABRIC mutations"),
        (METABRIC_EXPRESSION, "METABRIC expression"),
        (METABRIC_CNA, "METABRIC copy number"),
    ]
    missing = [(p, label) for p, label in expected if not p.exists()]
    if missing:
        logger.warning(f"Pre-flight: {len(missing)} expected data file(s) missing:")
        for p, label in missing:
            logger.warning(f"  MISSING: {label} — {p}")
        logger.warning("Loading steps will skip missing datasets. Re-run download or check data/.")
    else:
        logger.info("Pre-flight: all expected data files present.")


def run_all(force: bool = False, skip_download: bool = False) -> None:
    start = time.time()
    logger.info(f"BRCA ETL pipeline starting (DB: {DB_PATH})")

    if skip_download:
        logger.info("\n--- STEP 1: Download SKIPPED ---")
    else:
        logger.info("\n--- STEP 1: Download source data ---")
        download_all(force=False)

    _preflight_check()

    logger.info("\n--- STEP 2: Init DB ---")
    init_db(force=force)

    logger.info("\n--- STEP 3: Clinical ---")
    load_clinical()

    logger.info("\n--- STEP 4: Mutations ---")
    load_mutations()

    logger.info("\n--- STEP 5: Copy Number ---")
    load_copy_number()

    logger.info("\n--- STEP 6: Liftover Mutations (GRCh37 -> GRCh38) ---")
    try:
        liftover_mutations()
    except ImportError:
        logger.warning("pyliftover not installed. Skipping liftover. Run: pip install pyliftover")

    logger.info("\n--- STEP 7: ComBat Batch Correction (reads flat files) ---")
    try:
        run_combat()
    except ImportError:
        logger.warning("pyComBat not installed. Skipping ComBat. Run: pip install combat")

    print_summary()
    logger.info(f"Pipeline complete in {time.time() - start:.1f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()
    run_all(force=args.force, skip_download=args.skip_download)
