"""
Tests for etl/db.py helper functions.

Run with: python -m pytest tests/test_db.py -v
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import etl._bootstrap  # noqa: F401
from sqlalchemy import text

from config.settings import INSERT_BATCH_SIZE
from etl.db import bulk_insert, log_harmonization, table_row_count, truncate_table


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _insert_patient(conn, patient_id: str, study_source: str = "TCGA-BRCA") -> None:
    conn.execute(text(
        "INSERT INTO patients (patient_id, study_source) VALUES (:pid, :src)"
    ), {"pid": patient_id, "src": study_source})


# ---------------------------------------------------------------------------
# bulk_insert
# ---------------------------------------------------------------------------

class TestBulkInsert:
    def test_empty_list_returns_zero(self, db_conn):
        result = bulk_insert(db_conn, "patients", [])
        assert result == 0

    def test_inserts_rows(self, db_conn):
        rows = [
            {"patient_id": f"PT-{i}", "study_source": "TCGA-BRCA"}
            for i in range(3)
        ]
        bulk_insert(db_conn, "patients", rows)
        assert table_row_count(db_conn, "patients") == 3

    def test_returns_count(self, db_conn):
        rows = [
            {"patient_id": f"PT-{i}", "study_source": "METABRIC"}
            for i in range(5)
        ]
        result = bulk_insert(db_conn, "patients", rows)
        assert result == 5

    def test_or_ignore_skips_duplicate(self, db_conn):
        row = {"patient_id": "PT-DUP", "study_source": "TCGA-BRCA"}
        bulk_insert(db_conn, "patients", [row])
        bulk_insert(db_conn, "patients", [row])
        assert table_row_count(db_conn, "patients") == 1

    def test_batches_correctly(self, db_conn):
        n = INSERT_BATCH_SIZE + 1
        rows = [
            {"patient_id": f"PT-{i:05d}", "study_source": "TCGA-BRCA"}
            for i in range(n)
        ]
        result = bulk_insert(db_conn, "patients", rows)
        assert result == n
        assert table_row_count(db_conn, "patients") == n


# ---------------------------------------------------------------------------
# table_row_count
# ---------------------------------------------------------------------------

class TestTableRowCount:
    def test_empty_table(self, db_conn):
        assert table_row_count(db_conn, "patients") == 0

    def test_after_insert(self, db_conn):
        _insert_patient(db_conn, "PT-001")
        assert table_row_count(db_conn, "patients") == 1


# ---------------------------------------------------------------------------
# truncate_table
# ---------------------------------------------------------------------------

class TestTruncateTable:
    def test_clears_rows(self, db_conn):
        _insert_patient(db_conn, "PT-001")
        _insert_patient(db_conn, "PT-002")
        assert table_row_count(db_conn, "patients") == 2
        truncate_table(db_conn, "patients")
        assert table_row_count(db_conn, "patients") == 0


# ---------------------------------------------------------------------------
# log_harmonization
# ---------------------------------------------------------------------------

class TestLogHarmonization:
    def test_writes_entry(self, db_conn):
        before = db_conn.execute(
            text("SELECT count(*) FROM harmonization_log")
        ).scalar()
        log_harmonization(
            db_conn, "patients", "TCGA-BRCA", "vital_status",
            "Mapped LIVING/DECEASED to Alive/Deceased"
        )
        after = db_conn.execute(
            text("SELECT count(*) FROM harmonization_log")
        ).scalar()
        assert after == before + 1

    def test_optional_fields_null(self, db_conn):
        log_harmonization(
            db_conn, "tumors", "METABRIC", "pam50_subtype",
            "Normalized PAM50 labels"
        )
        row = db_conn.execute(
            text("SELECT source_field, notes FROM harmonization_log "
                 "ORDER BY log_id DESC LIMIT 1")
        ).fetchone()
        assert row[0] is None
        assert row[1] is None

    def test_fields_stored_correctly(self, db_conn):
        log_harmonization(
            db_conn, "mutations", "BOTH", "chromosome",
            "GRCh37 -> GRCh38 liftover",
            source_field="start_position",
            notes="Used pyliftover"
        )
        row = db_conn.execute(
            text("SELECT table_name, study_source, field_name, transformation, source_field, notes "
                 "FROM harmonization_log ORDER BY log_id DESC LIMIT 1")
        ).fetchone()
        assert row[0] == "mutations"
        assert row[1] == "BOTH"
        assert row[2] == "chromosome"
        assert row[3] == "GRCh37 -> GRCh38 liftover"
        assert row[4] == "start_position"
        assert row[5] == "Used pyliftover"
