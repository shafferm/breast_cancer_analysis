"""
ComBat batch correction for cross-study gene expression harmonization.

Uses pyComBat (Python implementation of Johnson et al. 2007) to correct
platform-driven batch effects between TCGA-BRCA (RNA-seq) and METABRIC
(Illumina microarray).

Prerequisites:
    pip install combat

Usage:
    python -m etl.combat_expression
    python -m etl.combat_expression --min-samples 10   # gene filter threshold
    python -m etl.combat_expression --dry-run           # preview only

Design decisions:
    1. Only genes measured on BOTH platforms are included.
    2. Expression values must be log2-transformed before ComBat (TCGA FPKM
       should already be log2(FPKM+1); METABRIC is log2 intensity).
    3. The batch variable is study_source (TCGA-BRCA vs METABRIC).
    4. Corrected values are written to a NEW column (expression_combat) via
       a separate table, preserving the original expression_value untouched.
    5. We use parametric ComBat (par_prior=True) as it is robust and fast
       for two-batch designs.

References:
    Johnson WE, et al. (2007) Adjusting batch effects in microarray
    expression data using empirical Bayes methods. Biostatistics 8:118-127.

    Behdenna A, et al. (2023) pyComBat, a Python tool for batch effects
    correction in high-throughput molecular data using empirical Bayes
    methods. BMC Bioinformatics 24:459.
"""

import logging
import argparse
from collections import defaultdict, Counter

import numpy as np
import pandas as pd  # pandas required by pycombat — not used elsewhere in this module
from sqlalchemy import text
from sqlalchemy.engine import Connection

import etl._bootstrap  # noqa: F401
from config.settings import INSERT_BATCH_SIZE
from etl.db import get_connection, log_harmonization

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Schema extension: combat-corrected expression table
# ---------------------------------------------------------------------------

CREATE_COMBAT_TABLE = """
CREATE TABLE IF NOT EXISTS gene_expression_combat (
    tumor_id            TEXT    NOT NULL,
    gene_symbol         TEXT    NOT NULL,
    expression_combat   REAL    NOT NULL,
    PRIMARY KEY (tumor_id, gene_symbol),
    FOREIGN KEY (tumor_id) REFERENCES tumors (tumor_id) ON DELETE CASCADE
)
"""

CREATE_COMBAT_INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_combat_gene ON gene_expression_combat (gene_symbol)",
    "CREATE INDEX IF NOT EXISTS idx_combat_tumor ON gene_expression_combat (tumor_id)",
]


def _ensure_combat_table(conn: Connection) -> None:
    """Create the combat-corrected expression table if it doesn't exist."""
    conn.execute(text(CREATE_COMBAT_TABLE))
    for idx_sql in CREATE_COMBAT_INDEXES:
        conn.execute(text(idx_sql))
    conn.commit()
    logger.info("  gene_expression_combat table ready.")


# ---------------------------------------------------------------------------
# Build the expression matrix from the database
# ---------------------------------------------------------------------------

def _build_expression_matrix(
    conn: Connection, min_samples: int = 10
) -> tuple[np.ndarray | None, list[str] | None, list[str] | None, list[str] | None]:
    """
    Build a genes x samples numpy expression matrix from the database,
    restricted to genes present in BOTH studies.

    Returns:
        matrix:      numpy ndarray (genes x samples), log2 expression values
        batches:     list of batch labels aligned to matrix columns
        gene_list:   list of gene symbols (row order)
        sample_list: list of tumor_ids (column order)
    """
    logger.info("Building expression matrix from database...")

    # Step 1: Find genes present in both platforms
    shared_genes = conn.execute(text("""
        SELECT gene_symbol
        FROM gene_expression
        GROUP BY gene_symbol
        HAVING COUNT(DISTINCT platform) = 2
    """)).fetchall()
    shared_gene_set = {r[0] for r in shared_genes}
    logger.info(f"  {len(shared_gene_set):,} genes on both platforms")

    if len(shared_gene_set) == 0:
        logger.error("No shared genes found. Cannot run ComBat.")
        return None, None, None, None

    # Step 2: Get sample -> study mapping
    sample_study = {}
    rows = conn.execute(text("""
        SELECT t.tumor_id, p.study_source
        FROM tumors t
        JOIN patients p ON t.patient_id = p.patient_id
    """)).fetchall()
    for r in rows:
        sample_study[r[0]] = r[1]

    # Step 3: Pull expression data for shared genes in chunks
    logger.info("  Loading expression values for shared genes...")

    gene_list = sorted(shared_gene_set)
    data = defaultdict(dict)  # gene -> {tumor_id: value}

    chunk_size = 1000
    for i in range(0, len(gene_list), chunk_size):
        chunk = gene_list[i:i + chunk_size]
        placeholders = ",".join(f":p{j}" for j in range(len(chunk)))
        params = {f"p{j}": g for j, g in enumerate(chunk)}
        rows = conn.execute(text(f"""
            SELECT gene_symbol, tumor_id, expression_value
            FROM gene_expression
            WHERE gene_symbol IN ({placeholders})
        """), params).fetchall()
        for gene, tid, val in rows:
            if tid in sample_study and val is not None:
                data[gene][tid] = val

        if (i + chunk_size) % 5000 == 0:
            logger.info(f"    Loaded {min(i + chunk_size, len(gene_list)):,}/{len(gene_list):,} genes")

    # Step 4: Build the matrix
    all_samples = set()
    for gene_data in data.values():
        all_samples.update(gene_data.keys())

    valid_samples = sorted([s for s in all_samples if s in sample_study])
    logger.info(f"  {len(valid_samples):,} samples with expression data")

    # Filter genes by minimum sample coverage
    filtered_genes = [
        gene for gene in gene_list
        if sum(1 for s in valid_samples if s in data.get(gene, {})) >= min_samples
    ]
    logger.info(f"  {len(filtered_genes):,} genes pass min_samples={min_samples} filter")

    if len(filtered_genes) == 0:
        logger.error("No genes pass the sample coverage filter.")
        return None, None, None, None

    # Build numpy matrix: genes (rows) x samples (columns)
    matrix = np.full((len(filtered_genes), len(valid_samples)), np.nan)
    for i, gene in enumerate(filtered_genes):
        gene_data = data.get(gene, {})
        for j, sample in enumerate(valid_samples):
            if sample in gene_data:
                matrix[i, j] = gene_data[sample]

    batches = [sample_study[s] for s in valid_samples]

    # Drop samples with too many NaN genes (>50% missing)
    max_missing_frac = 0.5
    missing_frac = np.mean(np.isnan(matrix), axis=0)
    keep_mask = missing_frac < max_missing_frac
    dropped = int(np.sum(~keep_mask))
    if dropped > 0:
        logger.warning(f"  Dropped {dropped} samples with >{max_missing_frac*100:.0f}% missing genes")
        matrix = matrix[:, keep_mask]
        valid_samples = [s for s, keep in zip(valid_samples, keep_mask) if keep]
        batches = [sample_study[s] for s in valid_samples]

    # Fill remaining NaN with per-gene mean (ComBat requires complete data)
    nan_count = int(np.sum(np.isnan(matrix)))
    if nan_count > 0:
        logger.info(f"  Imputing {nan_count:,} NaN values with per-gene mean")
        gene_means = np.nanmean(matrix, axis=1)
        for i in range(matrix.shape[0]):
            nan_mask = np.isnan(matrix[i])
            if nan_mask.any():
                matrix[i, nan_mask] = gene_means[i]

    batch_counts = Counter(batches)
    logger.info(f"  Matrix shape: {matrix.shape[0]:,} genes x {matrix.shape[1]:,} samples")
    for study, count in batch_counts.items():
        logger.info(f"    {study}: {count:,} samples")

    return matrix, batches, filtered_genes, valid_samples


# ---------------------------------------------------------------------------
# Run ComBat
# ---------------------------------------------------------------------------

def _run_combat(matrix: np.ndarray, batches: list[str],
                gene_list: list[str], sample_list: list[str]) -> np.ndarray:
    """
    Run pyComBat on the expression matrix.

    Args:
        matrix:      numpy ndarray, genes x samples
        batches:     list of batch labels (one per column)
        gene_list:   list of gene symbols (row order)
        sample_list: list of tumor_ids (column order)

    Returns:
        corrected_matrix: numpy ndarray, same shape as input
    """
    try:
        from combat.pycombat import pycombat
    except ImportError:
        logger.error(
            "pyComBat not installed. Run: pip install combat\n"
            "See: https://github.com/epigenelabs/pyComBat"
        )
        raise

    logger.info("Running ComBat batch correction...")
    logger.info("  Method: parametric empirical Bayes")
    logger.info(f"  Batch variable: study_source ({len(set(batches))} levels)")

    # pycombat requires a pandas DataFrame — single conversion at call boundary
    df_pd = pd.DataFrame(matrix, index=gene_list, columns=sample_list)
    df_corrected_pd = pycombat(df_pd, batches)
    corrected_matrix = df_corrected_pd.to_numpy()

    logger.info("  ComBat complete.")

    # Sanity check: corrected values should have reduced batch variance
    batch_arr = np.array(batches)
    unique_batches = list(set(batches))
    if len(unique_batches) == 2:
        b1, b2 = unique_batches
        before_diff = abs(
            matrix[:, batch_arr == b1].mean() -
            matrix[:, batch_arr == b2].mean()
        )
        after_diff = abs(
            corrected_matrix[:, batch_arr == b1].mean() -
            corrected_matrix[:, batch_arr == b2].mean()
        )
        logger.info(f"  Mean gene-level batch difference: {before_diff:.4f} (before) -> {after_diff:.4f} (after)")

    return corrected_matrix


# ---------------------------------------------------------------------------
# Write corrected values back to the database
# ---------------------------------------------------------------------------

def _write_combat_results(conn: Connection, corrected_matrix: np.ndarray,
                          gene_list: list[str], sample_list: list[str]) -> int:
    """
    Write ComBat-corrected expression values to gene_expression_combat table.
    """
    logger.info("Writing ComBat-corrected values to database...")

    conn.execute(text("DELETE FROM gene_expression_combat"))

    sql = text("""INSERT OR REPLACE INTO gene_expression_combat
             (tumor_id, gene_symbol, expression_combat)
             VALUES (:tumor_id, :gene_symbol, :expression_combat)""")

    batch = []
    total = 0
    n_genes = corrected_matrix.shape[0]

    for i, gene in enumerate(gene_list):
        for j, sample in enumerate(sample_list):
            val = corrected_matrix[i, j]
            if not np.isnan(val):
                batch.append({"tumor_id": sample, "gene_symbol": gene,
                               "expression_combat": round(float(val), 6)})

                if len(batch) >= INSERT_BATCH_SIZE:
                    conn.execute(sql, batch)
                    total += len(batch)
                    batch = []

        if (i + 1) % 2000 == 0:
            logger.info(f"  Written {i+1:,}/{n_genes:,} genes ({total:,} rows)")

    if batch:
        conn.execute(sql, batch)
        total += len(batch)

    conn.commit()
    logger.info(f"  Wrote {total:,} ComBat-corrected expression values")
    return total


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_combat(min_samples: int = 10, dry_run: bool = False) -> None:
    """
    Full ComBat batch correction pipeline:
    1. Build expression matrix from database (shared genes only)
    2. Run pyComBat
    3. Write corrected values to gene_expression_combat table
    """
    with get_connection() as conn:
        _ensure_combat_table(conn)

        matrix, batches, gene_list, sample_list = _build_expression_matrix(conn, min_samples)

        if matrix is None:
            logger.error("Cannot proceed without expression data.")
            return

        corrected_matrix = _run_combat(matrix, batches, gene_list, sample_list)

        if dry_run:
            logger.info("DRY RUN: Skipping database write.")
            n_valid = int(np.sum(~np.isnan(corrected_matrix)))
            logger.info(f"  Would write {n_valid:,} corrected values")
            logger.info(f"  Shape: {corrected_matrix.shape[0]:,} genes x {corrected_matrix.shape[1]:,} samples")
            return

        total = _write_combat_results(conn, corrected_matrix, gene_list, sample_list)

        log_harmonization(
            conn, "gene_expression_combat", "BOTH", "expression_combat",
            f"pyComBat parametric batch correction on {corrected_matrix.shape[0]:,} shared genes "
            f"across {corrected_matrix.shape[1]:,} samples. "
            f"Batch variable: study_source (TCGA-BRCA vs METABRIC). "
            f"{total:,} corrected values written.",
            "gene_expression.expression_value",
            "Original expression_value preserved in gene_expression table. "
            "ComBat assumes log-normal input; TCGA should be log2(FPKM+1), "
            "METABRIC is log2 intensity. Corrected values stored separately "
            "so analyses can use either raw within-study or corrected cross-study values."
        )

        logger.info("")
        logger.info("=" * 50)
        logger.info("COMBAT CORRECTION SUMMARY")
        logger.info("=" * 50)
        logger.info(f"  Shared genes:          {corrected_matrix.shape[0]:,}")
        logger.info(f"  Total samples:         {corrected_matrix.shape[1]:,}")
        logger.info(f"  Corrected values:      {total:,}")
        logger.info(f"  Output table:          gene_expression_combat")
        logger.info("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ComBat batch correction on expression data")
    parser.add_argument("--min-samples", type=int, default=10,
                        help="Minimum samples per gene (default: 10)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview without writing to DB")
    args = parser.parse_args()
    run_combat(min_samples=args.min_samples, dry_run=args.dry_run)
