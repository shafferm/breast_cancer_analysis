"""
Lift over METABRIC mutation coordinates from GRCh37 to GRCh38.

Uses pyliftover (pure-Python UCSC liftOver) to convert point coordinates.
The chain file is auto-downloaded on first use and cached locally.

Prerequisites:
    pip install pyliftover

Usage:
    python -m etl.liftover_mutations
    python -m etl.liftover_mutations --dry-run   # Preview without writing

Notes:
    - pyliftover converts point coordinates only (not ranges). For mutations,
      we lift start_position and end_position independently.
    - Variants that fail to lift (unmappable regions, deleted sequence, etc.)
      are flagged with reference_genome = 'GRCh37_UNMAPPED' so they remain
      in the database but are clearly excluded from coordinate-level analyses.
    - Gene-level analyses (e.g., mutation frequency by gene_symbol) are
      unaffected by liftover status.
"""

import logging
import argparse
from typing import TYPE_CHECKING

from sqlalchemy import text

if TYPE_CHECKING:
    from pyliftover import LiftOver

import etl._bootstrap  # noqa: F401
from etl.db import get_connection, log_harmonization

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def _get_converter() -> "LiftOver":
    """
    Initialize the pyliftover converter.
    Auto-downloads the hg19-to-hg38 chain file on first call (~250 KB).
    Cached in ~/.pyliftover/ after first download.
    """
    try:
        from pyliftover import LiftOver
    except ImportError:
        logger.error(
            "pyliftover not installed. Run: pip install pyliftover\n"
            "Alternatively, install the faster C++ version: pip install liftover"
        )
        raise

    logger.info("Initializing LiftOver converter (hg19 -> hg38)...")
    lo = LiftOver('hg19', 'hg38')
    logger.info("  Chain file loaded.")
    return lo


def _lift_position(lo: "LiftOver", chrom: str, pos: int) -> tuple[str | None, int | None]:
    """
    Lift a single genomic position from hg19 -> hg38.

    pyliftover expects 0-based coordinates internally, but UCSC chain files
    and our database both use 1-based. pyliftover handles this correctly
    when you pass 1-based coords (it subtracts 1 internally for the query
    and adds 1 back in the result).

    Actually, pyliftover uses 0-based throughout, so we convert:
      input:  1-based (from database) -> subtract 1 -> query
      output: 0-based result -> add 1 -> store as 1-based

    Returns (new_chrom, new_pos) or (None, None) if unmappable.
    """
    if chrom is None or pos is None:
        return None, None

    # Ensure 'chr' prefix for UCSC chain format
    chrom_str = chrom if chrom.startswith("chr") else f"chr{chrom}"

    # pyliftover uses 0-based coords
    result = lo.convert_coordinate(chrom_str, pos - 1)

    if result is None or len(result) == 0:
        return None, None

    # Take the first (highest-confidence) mapping
    new_chrom, new_pos_0based, strand, score = result[0]

    # Convert back to 1-based and strip 'chr' prefix to match our schema
    new_chrom_clean = new_chrom.replace("chr", "")
    new_pos = new_pos_0based + 1

    return new_chrom_clean, new_pos


def liftover_mutations(dry_run: bool = False) -> None:
    """
    Lift all METABRIC mutations from GRCh37 to GRCh38.

    Updates: chromosome, start_position, end_position, reference_genome
    For unmappable variants: sets reference_genome = 'GRCh37_UNMAPPED'
    """
    lo = _get_converter()

    with get_connection() as conn:
        # Fetch all METABRIC mutations still on GRCh37
        rows = conn.execute(text("""
            SELECT mutation_id, chromosome, start_position, end_position
            FROM mutations
            WHERE reference_genome = 'GRCh37'
        """)).fetchall()

        total = len(rows)
        if total == 0:
            logger.info("No GRCh37 mutations found. Nothing to lift over.")
            return

        logger.info(f"Lifting {total:,} METABRIC mutations from GRCh37 -> GRCh38...")

        lifted = 0
        failed = 0
        multi_chrom = 0

        for i, row in enumerate(rows):
            mut_id = row[0]
            chrom = row[1]
            start = row[2]
            end = row[3]

            # Lift start position
            new_chrom, new_start = _lift_position(lo, chrom, start)

            # Lift end position (if available)
            new_end = None
            if end is not None:
                end_chrom, new_end = _lift_position(lo, chrom, end)
                # Sanity check: start and end should land on the same chromosome
                if end_chrom is not None and end_chrom != new_chrom:
                    multi_chrom += 1
                    new_chrom = None  # Mark as unmappable

            if new_chrom is not None and new_start is not None:
                if not dry_run:
                    conn.execute(text("""
                        UPDATE mutations
                        SET chromosome = :chrom,
                            start_position = :start,
                            end_position = :end,
                            reference_genome = 'GRCh38'
                        WHERE mutation_id = :mut_id
                    """), {"chrom": new_chrom, "start": new_start, "end": new_end, "mut_id": mut_id})
                lifted += 1
            else:
                if not dry_run:
                    conn.execute(text("""
                        UPDATE mutations
                        SET reference_genome = 'GRCh37_UNMAPPED'
                        WHERE mutation_id = :mut_id
                    """), {"mut_id": mut_id})
                failed += 1

            if (i + 1) % 5000 == 0:
                logger.info(f"  Processed {i+1:,} / {total:,} (lifted: {lifted:,}, failed: {failed:,})")

        if not dry_run:
            conn.commit()

            log_harmonization(
                conn, "mutations", "METABRIC", "start_position",
                f"Lifted {lifted:,}/{total:,} mutations from GRCh37 to GRCh38 "
                f"using pyliftover (chain: hg19ToHg38.over.chain.gz). "
                f"{failed:,} variants unmappable.",
                "start_position, end_position, chromosome",
                f"Unmappable variants flagged with reference_genome='GRCh37_UNMAPPED'. "
                f"{multi_chrom} had start/end on different chromosomes after liftover."
            )

        # Summary
        logger.info("")
        logger.info("=" * 50)
        logger.info("LIFTOVER SUMMARY")
        logger.info("=" * 50)
        logger.info(f"  Total METABRIC mutations:  {total:,}")
        logger.info(f"  Successfully lifted:       {lifted:,} ({lifted/total*100:.1f}%)")
        logger.info(f"  Failed / unmappable:       {failed:,} ({failed/total*100:.1f}%)")
        if multi_chrom > 0:
            logger.info(f"  Cross-chromosome splits:   {multi_chrom:,}")
        logger.info(f"  Mode:                      {'DRY RUN' if dry_run else 'COMMITTED'}")
        logger.info("=" * 50)

        if not dry_run:
            # Verify
            counts = conn.execute(text("""
                SELECT reference_genome, count(*)
                FROM mutations
                GROUP BY reference_genome
            """)).fetchall()
            logger.info("\nMutations by reference_genome after liftover:")
            for r in counts:
                logger.info(f"  {r[0]:25s} {r[1]:>10,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Liftover METABRIC mutations GRCh37 -> GRCh38")
    parser.add_argument("--dry-run", action="store_true", help="Preview without writing to DB")
    args = parser.parse_args()
    liftover_mutations(dry_run=args.dry_run)
