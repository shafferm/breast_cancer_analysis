"""
Database connection helpers for the BRCA harmonized ETL pipeline.

Uses SQLAlchemy Core (no ORM) for connection management and raw SQL execution.
"""

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator

from sqlalchemy import create_engine, event, text
from sqlalchemy.engine import Connection, Engine

import etl._bootstrap  # noqa: F401
from config.settings import DB_PATH, INSERT_BATCH_SIZE

logger = logging.getLogger(__name__)

_engines = {}


def get_engine(db_path: Path = DB_PATH) -> Engine:
    db_path = str(db_path)
    if db_path not in _engines:
        engine = create_engine(f"sqlite:///{db_path}")

        @event.listens_for(engine, "connect")
        def _set_pragmas(dbapi_conn, connection_record):
            cursor = dbapi_conn.cursor()
            cursor.execute("PRAGMA journal_mode = WAL")
            cursor.execute("PRAGMA foreign_keys = ON")
            cursor.execute("PRAGMA synchronous = NORMAL")
            cursor.close()

        _engines[db_path] = engine
    return _engines[db_path]


@contextmanager
def get_connection(db_path: Path = DB_PATH) -> Generator[Connection, None, None]:
    engine = get_engine(db_path)
    conn = engine.connect()
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def execute_script(conn: Connection, ddl: str) -> None:
    """Execute a multi-statement DDL script."""
    for stmt in ddl.split(";"):
        stmt = stmt.strip()
        if stmt and not stmt.upper().startswith("PRAGMA"):
            conn.execute(text(stmt))


def bulk_insert(conn: Connection, table: str, rows: list[dict[str, Any]],
                batch_size: int = INSERT_BATCH_SIZE) -> int:
    if not rows:
        return 0
    columns = list(rows[0].keys())
    placeholders = ", ".join(f":{c}" for c in columns)
    col_names = ", ".join(columns)
    sql = text(f"INSERT OR IGNORE INTO {table} ({col_names}) VALUES ({placeholders})")
    total = 0
    for i in range(0, len(rows), batch_size):
        batch = rows[i : i + batch_size]
        conn.execute(sql, batch)
        total += len(batch)
        if total % (batch_size * 10) == 0 and total > 0:
            logger.info(f"  {table}: inserted {total:,} / {len(rows):,} rows")
    conn.commit()
    logger.info(f"  {table}: inserted {total:,} rows total")
    return total


def log_harmonization(conn: Connection, table_name: str, study_source: str,
                      field_name: str, transformation: str,
                      source_field: str | None = None,
                      notes: str | None = None) -> None:
    conn.execute(
        text("""INSERT INTO harmonization_log
           (table_name, study_source, field_name, transformation, source_field, notes)
           VALUES (:table_name, :study_source, :field_name, :transformation, :source_field, :notes)"""),
        {"table_name": table_name, "study_source": study_source, "field_name": field_name,
         "transformation": transformation, "source_field": source_field, "notes": notes},
    )


def table_row_count(conn: Connection, table: str) -> int:
    return conn.execute(text(f"SELECT count(*) FROM {table}")).scalar()


def truncate_table(conn: Connection, table: str) -> None:
    conn.execute(text(f"DELETE FROM {table}"))
    conn.commit()
    logger.info(f"  Truncated {table}")
