"""
Tests for the BRCA harmonized ETL pipeline.

Run with: python -m pytest tests/ -v
"""

import sys
from pathlib import Path

# Ensure project root is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import etl._bootstrap  # noqa: F401
from sqlalchemy import text
from sqlalchemy.exc import IntegrityError

from etl.harmonize import (
    clean_value, clean_float, clean_int, days_to_months,
    normalize_pam50, normalize_receptor_status, derive_tnbc,
    normalize_stage, normalize_histology, infer_menopausal_status,
    normalize_vital_status, normalize_variant_classification,
    compute_npi, parse_tcga_barcode,
)
from config.metabric_panel_genes import METABRIC_PANEL_GENES


# ---------------------------------------------------------------------------
# clean_value
# ---------------------------------------------------------------------------

class TestCleanValue:
    def test_none(self):
        assert clean_value(None) is None

    def test_na_string(self):
        assert clean_value("NA") is None

    def test_not_available(self):
        assert clean_value("[Not Available]") is None

    def test_empty_string(self):
        assert clean_value("") is None

    def test_valid_string(self):
        assert clean_value("Hello") == "Hello"

    def test_whitespace_stripped(self):
        assert clean_value("  Hello  ") == "Hello"

    def test_nan_string(self):
        assert clean_value("nan") is None


# ---------------------------------------------------------------------------
# clean_float / clean_int
# ---------------------------------------------------------------------------

class TestCleanNumeric:
    def test_float_valid(self):
        assert clean_float("3.14") == 3.14

    def test_float_na(self):
        assert clean_float("NA") is None

    def test_float_none(self):
        assert clean_float(None) is None

    def test_float_integer_string(self):
        assert clean_float("42") == 42.0

    def test_int_valid(self):
        assert clean_int("3") == 3

    def test_int_from_float_string(self):
        assert clean_int("3.7") == 3

    def test_int_na(self):
        assert clean_int("NA") is None


# ---------------------------------------------------------------------------
# days_to_months
# ---------------------------------------------------------------------------

class TestDaysToMonths:
    def test_one_year(self):
        result = days_to_months(365)
        assert abs(result - 11.99) < 0.1

    def test_none(self):
        assert days_to_months(None) is None

    def test_zero(self):
        assert days_to_months(0) == 0.0


# ---------------------------------------------------------------------------
# PAM50
# ---------------------------------------------------------------------------

class TestNormalizePam50:
    def test_luminal_a(self):
        assert normalize_pam50("Luminal A") == "LumA"

    def test_basal_like(self):
        assert normalize_pam50("Basal-like") == "Basal"

    def test_her2_enriched(self):
        assert normalize_pam50("HER2-enriched") == "Her2"

    def test_empty(self):
        assert normalize_pam50("") == "Unknown"

    def test_none(self):
        assert normalize_pam50(None) == "Unknown"

    def test_lowercase(self):
        assert normalize_pam50("lumA") == "LumA"


# ---------------------------------------------------------------------------
# Receptor status
# ---------------------------------------------------------------------------

class TestReceptorStatus:
    def test_positive(self):
        assert normalize_receptor_status("Positive") == "Positive"

    def test_negative(self):
        assert normalize_receptor_status("Negative") == "Negative"

    def test_not_evaluated(self):
        assert normalize_receptor_status("[Not Evaluated]") == "Unknown"

    def test_empty(self):
        assert normalize_receptor_status("") == "Unknown"

    def test_none(self):
        assert normalize_receptor_status(None) == "Unknown"


# ---------------------------------------------------------------------------
# TNBC derivation
# ---------------------------------------------------------------------------

class TestDeriveTnbc:
    def test_triple_negative(self):
        assert derive_tnbc("Negative", "Negative", "Negative") == 1

    def test_er_positive(self):
        assert derive_tnbc("Positive", "Negative", "Negative") == 0

    def test_unknown(self):
        assert derive_tnbc("Negative", "Negative", "Unknown") is None


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------

class TestNormalizeStage:
    def test_stage_iia(self):
        assert normalize_stage("Stage IIA") == "IIA"

    def test_stage_x(self):
        assert normalize_stage("Stage X") is None

    def test_short_form(self):
        assert normalize_stage("IIIC") == "IIIC"

    def test_none(self):
        assert normalize_stage(None) is None


# ---------------------------------------------------------------------------
# Vital status
# ---------------------------------------------------------------------------

class TestVitalStatus:
    def test_living(self):
        assert normalize_vital_status("LIVING") == "Alive"

    def test_numeric_one(self):
        assert normalize_vital_status("1") == "Deceased"

    def test_deceased(self):
        assert normalize_vital_status("Deceased") == "Deceased"

    def test_none(self):
        assert normalize_vital_status(None) is None


# ---------------------------------------------------------------------------
# Variant classification
# ---------------------------------------------------------------------------

class TestVariantClassification:
    def test_missense_variant(self):
        assert normalize_variant_classification("missense_variant") == "Missense_Mutation"

    def test_stop_gained(self):
        assert normalize_variant_classification("stop_gained") == "Nonsense_Mutation"

    def test_passthrough(self):
        assert normalize_variant_classification("Custom_Type") == "Custom_Type"

    def test_none(self):
        assert normalize_variant_classification(None) is None


# ---------------------------------------------------------------------------
# Menopausal status
# ---------------------------------------------------------------------------

class TestMenopausalStatus:
    def test_pre(self):
        assert infer_menopausal_status(45) == "Pre"

    def test_post(self):
        assert infer_menopausal_status(55) == "Post"

    def test_boundary(self):
        assert infer_menopausal_status(50) == "Post"

    def test_none(self):
        assert infer_menopausal_status(None) == "Unknown"


# ---------------------------------------------------------------------------
# NPI
# ---------------------------------------------------------------------------

class TestComputeNpi:
    def test_low_risk(self):
        # 2cm, grade 2, 0 nodes: 0.4 + 2 + 1 = 3.4
        assert compute_npi(2.0, 2, 0) == 3.4

    def test_high_risk(self):
        # 3cm, grade 3, 5 nodes: 0.6 + 3 + 3 = 6.6
        assert compute_npi(3.0, 3, 5) == 6.6

    def test_mid_nodes(self):
        # 1cm, grade 1, 2 nodes: 0.2 + 1 + 2 = 3.2
        assert compute_npi(1.0, 1, 2) == 3.2

    def test_missing(self):
        assert compute_npi(None, 2, 0) is None
        assert compute_npi(2.0, None, 0) is None
        assert compute_npi(2.0, 2, None) is None


# ---------------------------------------------------------------------------
# TCGA barcode parsing
# ---------------------------------------------------------------------------

class TestTcgaBarcode:
    def test_full_barcode(self):
        bc = parse_tcga_barcode("TCGA-A2-A0YF-01A-11R-A034-13")
        assert bc["patient_id"] == "TCGA-A2-A0YF"
        assert bc["tumor_id"] == "TCGA-A2-A0YF-01A"
        assert bc["sample_type"] == "Primary Solid Tumor"
        assert bc["sample_type_code"] == "01"

    def test_normal_sample(self):
        bc = parse_tcga_barcode("TCGA-A2-A0YF-11A-11R-A034-13")
        assert bc["sample_type"] == "Solid Tissue Normal"

    def test_short_barcode(self):
        bc = parse_tcga_barcode("TCGA-A2-A0YF")
        assert bc["patient_id"] == "TCGA-A2-A0YF"
        assert bc["tumor_id"] == "TCGA-A2-A0YF"


# ---------------------------------------------------------------------------
# METABRIC panel
# ---------------------------------------------------------------------------

class TestMetabricPanel:
    def test_key_genes_present(self):
        for gene in ["TP53", "PIK3CA", "GATA3", "CDH1", "MAP3K1", "PTEN", "BRCA1", "BRCA2"]:
            assert gene in METABRIC_PANEL_GENES, f"{gene} missing from panel"

    def test_panel_size(self):
        assert len(METABRIC_PANEL_GENES) > 100

    def test_is_frozenset(self):
        assert isinstance(METABRIC_PANEL_GENES, frozenset)


# ---------------------------------------------------------------------------
# DB initialization
# ---------------------------------------------------------------------------

class TestDbInit:
    def test_init_creates_tables(self, tmp_path):
        from config.settings import DDL_PATH
        from etl.db import get_connection, execute_script

        db_path = tmp_path / "test.db"
        ddl = DDL_PATH.read_text()

        with get_connection(db_path) as conn:
            execute_script(conn, ddl)
            tables = [r[0] for r in conn.execute(text(
                "SELECT name FROM sqlite_master WHERE type='table' AND name != 'sqlite_sequence'"
            )).fetchall()]

        expected = {"patients", "tumors", "treatments", "mutations",
                    "gene_expression", "copy_number", "methylation",
                    "protein_expression", "harmonization_log"}
        assert set(tables) == expected

    def test_foreign_key_enforcement(self, tmp_path):
        from config.settings import DDL_PATH
        from etl.db import get_connection, execute_script

        db_path = tmp_path / "test_fk.db"
        ddl = DDL_PATH.read_text()

        with get_connection(db_path) as conn:
            execute_script(conn, ddl)
            try:
                conn.execute(text(
                    "INSERT INTO tumors (tumor_id, patient_id) VALUES ('FAKE', 'NONEXISTENT')"
                ))
                assert False, "FK constraint should have raised an error"
            except IntegrityError:
                pass  # Expected

    def test_check_constraint(self, tmp_path):
        from config.settings import DDL_PATH
        from etl.db import get_connection, execute_script

        db_path = tmp_path / "test_check.db"
        ddl = DDL_PATH.read_text()

        with get_connection(db_path) as conn:
            execute_script(conn, ddl)
            try:
                conn.execute(text(
                    "INSERT INTO patients (patient_id, study_source) VALUES ('TEST', 'INVALID')"
                ))
                assert False, "CHECK constraint should have raised an error"
            except IntegrityError:
                pass  # Expected


# ---------------------------------------------------------------------------
# normalize_histology
# ---------------------------------------------------------------------------

class TestNormalizeHistology:
    def test_idc(self):
        assert normalize_histology("Infiltrating Ductal Carcinoma") == "IDC"

    def test_ilc(self):
        assert normalize_histology("Infiltrating Lobular Carcinoma") == "ILC"

    def test_metaplastic(self):
        assert normalize_histology("Metaplastic Carcinoma") == "Metaplastic"

    def test_short_form(self):
        assert normalize_histology("IDC") == "IDC"

    def test_passthrough(self):
        assert normalize_histology("Rare Type") == "Rare Type"

    def test_none(self):
        assert normalize_histology(None) is None
