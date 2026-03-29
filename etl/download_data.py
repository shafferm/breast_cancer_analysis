"""
Download TCGA-BRCA and METABRIC source data from cBioPortal and GDC.

Both cBioPortal tarballs are public S3 downloads (no authentication required).
The TCGA CDR clinical file is fetched from the GDC API with a GitHub mirror fallback.

Usage:
    python -m etl.download_data            # download everything missing
    python -m etl.download_data --tcga     # TCGA-BRCA + CDR only
    python -m etl.download_data --metabric # METABRIC only
    python -m etl.download_data --force    # re-download even if files exist

Called automatically as step 0 of etl.run_all.
"""

import argparse
import logging
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable

import etl._bootstrap  # noqa: F401
from config.settings import TCGA_DIR, METABRIC_DIR

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Download URLs
# ---------------------------------------------------------------------------

TCGA_URL = "https://datahub.assets.cbioportal.org/brca_tcga.tar.gz"
METABRIC_URL = "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz"

CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_FALLBACK = (
    "https://github.com/umccr/TCGA-data-harmonization/raw/refs/heads/master"
    "/docs/PanCanAtlas/TCGA-CDR-SupplementalTableS1.xlsx"
)
CDR_FILENAME = "TCGA-CDR-SupplementalTableS1.xlsx"

# ---------------------------------------------------------------------------
# Files expected after extraction (used for validation warnings)
# ---------------------------------------------------------------------------

TCGA_EXPECTED = [
    "data_mutations.txt",
    "data_mrna_seq_v2_rsem.txt",
    "data_cna.txt",
]
METABRIC_EXPECTED = [
    "data_clinical_patient.txt",
    "data_clinical_sample.txt",
    "data_mutations.txt",
    "data_mrna_illumina_microarray.txt",
    "data_cna.txt",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_expected_files(directory: Path, expected: list[str]) -> bool:
    """True if directory contains all expected files from tarball extraction."""
    return directory.is_dir() and all((directory / f).exists() for f in expected)


def _make_reporthook(label: str) -> Callable:
    """Return a urlretrieve reporthook that logs progress every ~10%."""
    state: dict = {"last_pct": -1}

    def hook(block_num: int, block_size: int, total_size: int) -> None:
        if total_size <= 0:
            return
        downloaded = block_num * block_size
        pct = min(100, int(downloaded * 100 / total_size))
        milestone = pct // 10 * 10
        if milestone > state["last_pct"]:
            state["last_pct"] = milestone
            mb = downloaded / 1_048_576
            total_mb = total_size / 1_048_576
            logger.info(f"  {label}: {mb:.1f} / {total_mb:.1f} MB ({pct}%)")

    return hook


def _download(url: str, dest: Path, label: str) -> None:
    """Download url to dest, logging progress."""
    logger.info(f"  Downloading {label}...")
    logger.info(f"  URL: {url}")
    urllib.request.urlretrieve(url, dest, reporthook=_make_reporthook(label))
    size_mb = dest.stat().st_size / 1_048_576
    logger.info(f"  {label}: {size_mb:.1f} MB downloaded.")


def _extract_tarball(tarball: Path, dest_dir: Path, label: str) -> None:
    """
    Extract tarball into dest_dir, stripping the top-level directory.
    Equivalent to: tar -xzf tarball --strip-components=1 -C dest_dir
    """
    logger.info(f"  Extracting {label}...")
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tarball) as tf:
        for member in tf.getmembers():
            # Strip the top-level directory component
            parts = member.name.split("/", 1)
            if len(parts) == 1:
                continue  # top-level dir entry — skip
            member.name = parts[1]
            if not member.name:
                continue
            tf.extract(member, dest_dir)
    logger.info(f"  {label}: extraction complete.")


def _validate(directory: Path, expected_files: list[str], label: str) -> None:
    """Warn if any expected file is missing after extraction."""
    missing = [f for f in expected_files if not (directory / f).exists()]
    if missing:
        logger.warning(f"  [{label}] Missing expected files: {', '.join(missing)}")
    else:
        logger.info(f"  [{label}] All expected files present.")


# ---------------------------------------------------------------------------
# Individual downloaders
# ---------------------------------------------------------------------------

def download_metabric(force: bool = False) -> None:
    """Download and extract the METABRIC cBioPortal tarball."""
    logger.info("METABRIC download:")
    if not force and _has_expected_files(METABRIC_DIR, METABRIC_EXPECTED):
        logger.info(f"  {METABRIC_DIR} already has expected files — skipping.")
        logger.info("  Pass --force to re-download.")
        return

    METABRIC_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        _download(METABRIC_URL, tmp_path, "METABRIC")
        _extract_tarball(tmp_path, METABRIC_DIR, "METABRIC")
    finally:
        tmp_path.unlink(missing_ok=True)

    _validate(METABRIC_DIR, METABRIC_EXPECTED, "METABRIC")
    logger.info("METABRIC download complete.")


def download_tcga(force: bool = False) -> None:
    """Download and extract the TCGA-BRCA cBioPortal tarball, then fetch the CDR file."""
    logger.info("TCGA-BRCA download:")
    if not force and _has_expected_files(TCGA_DIR, TCGA_EXPECTED):
        logger.info(f"  {TCGA_DIR} already has expected files — skipping tarball.")
        logger.info("  Pass --force to re-download.")
    else:
        TCGA_DIR.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
            tmp_path = Path(tmp.name)

        try:
            _download(TCGA_URL, tmp_path, "TCGA-BRCA")
            _extract_tarball(tmp_path, TCGA_DIR, "TCGA-BRCA")
        finally:
            tmp_path.unlink(missing_ok=True)

        _validate(TCGA_DIR, TCGA_EXPECTED, "TCGA-BRCA")

    download_cdr(force=force)
    logger.info("TCGA-BRCA download complete.")


def download_cdr(force: bool = False) -> None:
    """Download the TCGA Pan-Cancer CDR Excel file (GDC API, with GitHub mirror fallback)."""
    dest = TCGA_DIR / CDR_FILENAME
    if not force and dest.exists():
        logger.info(f"  CDR file already present at {dest} — skipping.")
        return

    TCGA_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("  Downloading TCGA CDR clinical file (GDC API)...")
    try:
        _download(CDR_URL, dest, "CDR")
        # Sanity check: GDC API returns the file if UUID is valid;
        # a tiny response likely means a JSON error message, not the Excel
        if dest.stat().st_size < 50_000:
            raise ValueError(f"CDR file suspiciously small ({dest.stat().st_size} bytes) — likely an error response")
    except Exception as exc:
        logger.warning(f"  GDC API failed ({exc}); trying GitHub mirror...")
        dest.unlink(missing_ok=True)
        _download(CDR_FALLBACK, dest, "CDR (mirror)")

    logger.info(f"  CDR file saved to {dest}")


# ---------------------------------------------------------------------------
# Public entry point (called by run_all)
# ---------------------------------------------------------------------------

def download_all(force: bool = False) -> None:
    """Download all source data. Skips datasets that already have files."""
    download_tcga(force=force)
    download_metabric(force=force)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download TCGA-BRCA and METABRIC source data from cBioPortal / GDC"
    )
    parser.add_argument("--tcga", action="store_true", help="Download TCGA-BRCA only")
    parser.add_argument("--metabric", action="store_true", help="Download METABRIC only")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if files already exist")
    args = parser.parse_args()

    do_tcga = args.tcga or not (args.tcga or args.metabric)
    do_metabric = args.metabric or not (args.tcga or args.metabric)

    if do_tcga:
        download_tcga(force=args.force)
    if do_metabric:
        download_metabric(force=args.force)
