"""
Central configuration for the BRCA harmonized ETL pipeline.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Project layout
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = PROJECT_ROOT / "brca_harmonized.db"
DDL_PATH = PROJECT_ROOT / "schema" / "brca_harmonized_schema.sql"

DATA_DIR = PROJECT_ROOT / "data"
TCGA_DIR = DATA_DIR / "tcga"
METABRIC_DIR = DATA_DIR / "metabric"

# ---------------------------------------------------------------------------
# TCGA-BRCA source files
# ---------------------------------------------------------------------------
TCGA_CDR_FILE = TCGA_DIR / "TCGA-CDR-SupplementalTableS1.xlsx"
TCGA_CLINICAL_PATIENT = TCGA_DIR / "nationwidechildrens.org_clinical_patient_brca.txt"
TCGA_CLINICAL_DRUG = TCGA_DIR / "nationwidechildrens.org_clinical_drug_brca.txt"
TCGA_CLINICAL_CBIO = TCGA_DIR / "data_clinical_patient.txt"
TCGA_PANCAN_SUBTYPES = TCGA_DIR / "TCGASubtype.20170308.tsv"
TCGA_CLINICAL_DRUG_GDC = TCGA_DIR / "nationwidechildrens.org_clinical_drug_public_brca.txt"
TCGA_CLINICAL_RADIATION_GDC = TCGA_DIR / "nationwidechildrens.org_clinical_radiation_public_brca.txt"
TCGA_MUTATIONS_MAF = TCGA_DIR / "data_mutations.txt"
TCGA_EXPRESSION = TCGA_DIR / "data_mrna_seq_v2_rsem.txt"
TCGA_CNA = TCGA_DIR / "data_cna.txt"

# ---------------------------------------------------------------------------
# METABRIC source files (cBioPortal download)
# ---------------------------------------------------------------------------
METABRIC_CLINICAL_PATIENT = METABRIC_DIR / "data_clinical_patient.txt"
METABRIC_CLINICAL_SAMPLE = METABRIC_DIR / "data_clinical_sample.txt"
METABRIC_MUTATIONS = METABRIC_DIR / "data_mutations.txt"
METABRIC_EXPRESSION = METABRIC_DIR / "data_mrna_illumina_microarray.txt"
METABRIC_CNA = METABRIC_DIR / "data_cna.txt"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DAYS_PER_MONTH = 30.44
MENOPAUSE_AGE_CUTOFF = 50
TARGET_GENOME_BUILD = "GRCh38"
INSERT_BATCH_SIZE = 5000

# ---------------------------------------------------------------------------
# Vocabulary maps
# ---------------------------------------------------------------------------
PAM50_MAP = {
    "Luminal A": "LumA", "Luminal B": "LumB", "HER2-enriched": "Her2",
    "Basal-like": "Basal", "Normal-like": "Normal",
    "LumA": "LumA", "LumB": "LumB", "Her2": "Her2", "Basal": "Basal", "Normal": "Normal",
    "lumA": "LumA", "lumB": "LumB", "her2": "Her2", "basal": "Basal", "normal": "Normal",
    "claudin-low": "Basal", "NC": "Unknown", "": "Unknown",
}

RECEPTOR_STATUS_MAP = {
    "Positive": "Positive", "positive": "Positive",
    "Negative": "Negative", "negative": "Negative",
    "Equivocal": "Equivocal", "equivocal": "Equivocal",
    "Indeterminate": "Unknown", "[Not Evaluated]": "Unknown",
    "[Not Available]": "Unknown", "": "Unknown",
}

STAGE_MAP = {
    "Stage I": "I", "Stage IA": "IA", "Stage IB": "IB",
    "Stage II": "II", "Stage IIA": "IIA", "Stage IIB": "IIB",
    "Stage III": "III", "Stage IIIA": "IIIA", "Stage IIIB": "IIIB", "Stage IIIC": "IIIC",
    "Stage IV": "IV", "Stage X": None,
    "I": "I", "IA": "IA", "IB": "IB", "II": "II", "IIA": "IIA", "IIB": "IIB",
    "III": "III", "IIIA": "IIIA", "IIIB": "IIIB", "IIIC": "IIIC", "IV": "IV",
    "0": None, "1": "I", "2": "II", "3": "III", "4": "IV",
}

HISTOLOGY_MAP = {
    "Infiltrating Ductal Carcinoma": "IDC", "Infiltrating Lobular Carcinoma": "ILC",
    "Medullary Carcinoma": "Medullary", "Metaplastic Carcinoma": "Metaplastic",
    "Mixed Histology (please specify)": "Mixed", "Mucinous Carcinoma": "Mucinous",
    "Other, specify": "Other",
    "Breast Invasive Ductal Carcinoma": "IDC", "Breast Invasive Lobular Carcinoma": "ILC",
    "Breast Mixed Ductal and Lobular Carcinoma": "Mixed",
    "Breast Invasive Mixed Mucinous Carcinoma": "Mucinous",
    "Invasive Ductal Carcinoma": "IDC", "Invasive Lobular Carcinoma": "ILC",
    "IDC": "IDC", "ILC": "ILC", "MIXED": "Mixed",
    "BREAST_METAPLASTIC_CARCINOMA": "Metaplastic", "BREAST_IDC_AND_ILC": "Mixed",
}

VARIANT_CLASS_MAP = {
    "Missense_Mutation": "Missense_Mutation", "missense_variant": "Missense_Mutation",
    "Nonsense_Mutation": "Nonsense_Mutation", "stop_gained": "Nonsense_Mutation",
    "Frame_Shift_Del": "Frame_Shift_Del", "frameshift_variant": "Frame_Shift_Del",
    "Frame_Shift_Ins": "Frame_Shift_Ins",
    "Splice_Site": "Splice_Site", "splice_acceptor_variant": "Splice_Site",
    "splice_donor_variant": "Splice_Site",
    "In_Frame_Del": "In_Frame_Del", "inframe_deletion": "In_Frame_Del",
    "In_Frame_Ins": "In_Frame_Ins", "inframe_insertion": "In_Frame_Ins",
    "Translation_Start_Site": "Translation_Start_Site",
    "Nonstop_Mutation": "Nonstop_Mutation",
    "Silent": "Silent", "synonymous_variant": "Silent",
    "5'UTR": "5'UTR", "3'UTR": "3'UTR", "Intron": "Intron",
    "RNA": "RNA", "IGR": "IGR", "5'Flank": "5'Flank", "3'Flank": "3'Flank",
}
