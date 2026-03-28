"""
Shared pytest fixtures for the BRCA harmonized ETL test suite.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import etl._bootstrap  # noqa: F401
from config.settings import DDL_PATH
from etl.db import get_connection, execute_script


@pytest.fixture
def db_conn(tmp_path):
    """Schema-initialized SQLite connection for a temp database."""
    db_path = tmp_path / "test.db"
    with get_connection(db_path) as conn:
        execute_script(conn, DDL_PATH.read_text())
        yield conn
