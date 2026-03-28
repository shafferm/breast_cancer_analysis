"""
Initialize the SQLite database from the DDL schema file.

Usage:
    python -m etl.init_db
    python -m etl.init_db --force
"""

import argparse
import logging
import sys

from sqlalchemy import text

import etl._bootstrap  # noqa: F401
from config.settings import DB_PATH, DDL_PATH
from etl.db import get_connection, execute_script

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def init_db(force: bool = False) -> None:
    if not DDL_PATH.exists():
        logger.error(f"DDL file not found: {DDL_PATH}")
        logger.error("Copy brca_harmonized_schema.sql into the project root.")
        sys.exit(1)

    if DB_PATH.exists():
        if force:
            logger.warning(f"--force: removing existing database at {DB_PATH}")
            DB_PATH.unlink()
        else:
            logger.info(f"Database already exists at {DB_PATH}")
            logger.info("Use --force to drop and recreate.")
            return

    logger.info(f"Creating database at {DB_PATH}")
    ddl = DDL_PATH.read_text()

    with get_connection(DB_PATH) as conn:
        execute_script(conn, ddl)

    with get_connection(DB_PATH) as conn:
        tables = conn.execute(
            text("SELECT name FROM sqlite_master WHERE type='table' AND name != 'sqlite_sequence' ORDER BY name")
        ).fetchall()
        views = conn.execute(
            text("SELECT name FROM sqlite_master WHERE type='view' ORDER BY name")
        ).fetchall()
        indexes = conn.execute(
            text("SELECT count(*) FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'")
        ).scalar()

    logger.info(f"Created {len(tables)} tables: {', '.join(r[0] for r in tables)}")
    logger.info(f"Created {len(views)} views: {', '.join(r[0] for r in views)}")
    logger.info(f"Created {indexes} indexes")
    logger.info("Database initialization complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Initialize the BRCA harmonized database")
    parser.add_argument("--force", action="store_true", help="Drop and recreate the database")
    args = parser.parse_args()
    init_db(force=args.force)
