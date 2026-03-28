# BRCA Harmonized Database

A cross-study breast cancer database harmonizing TCGA-BRCA and METABRIC into a
unified SQLite schema with multi-omic data (clinical, mutations, gene expression,
copy number).

## Quick Start

```bash
# 1. Clone and set up
git clone <repo-url>
cd brca_harmonized
conda env create -f environment.yaml
conda activate brca_analysis

# 2. Run the full pipeline (downloads data automatically on first run)
python -m etl.run_all --force
```

## Pipeline Steps

| Step | Module                    | Description                                       |
|------|---------------------------|---------------------------------------------------|
| 1    | `etl.download_data`       | Download source data from cBioPortal / GDC        |
| 2    | `etl.init_db`             | Create SQLite database from DDL                   |
| 3    | `etl.load_clinical`       | Load patients, tumors, treatments                 |
| 4    | `etl.load_mutations`      | Load somatic mutations (WES + panel)              |
| 5    | `etl.load_expression`     | Load gene expression with z-scores                |
| 6    | `etl.load_copy_number`    | Load GISTIC2 copy number alterations              |
| 7    | `etl.liftover_mutations`  | Lift METABRIC coords GRCh37 -> GRCh38             |
| 8    | `etl.combat_expression`   | ComBat batch correction across platforms           |

Each step can be run independently:

```bash
python -m etl.download_data          # Download source data only
python -m etl.load_clinical          # Just clinical data
python -m etl.liftover_mutations     # Just liftover
python -m etl.combat_expression      # Just ComBat
python -m etl.run_all --skip-expr    # Skip expression (fast iteration)
python -m etl.run_all --skip-download  # Skip download (data already present)
```

## Data Acquisition

Source data is downloaded automatically when running `python -m etl.run_all`.
To download without running the pipeline:

```bash
python -m etl.download_data            # all datasets
python -m etl.download_data --tcga     # TCGA-BRCA + CDR only
python -m etl.download_data --metabric # METABRIC only
python -m etl.download_data --force    # re-download even if files exist
```

The tables below document what is fetched and from where.

### TCGA-BRCA

Files are downloaded from cBioPortal (tarball) and the GDC API (CDR). Placed in `data/tcga/`.

| File | Source | Auto-downloaded |
|------|--------|-----------------|
| `TCGA-CDR-SupplementalTableS1.xlsx` | [GDC API](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) | Yes |
| `data_mutations.txt` | [cBioPortal brca_tcga](https://www.cbioportal.org/study?id=brca_tcga) | Yes |
| `data_mrna_seq_v2_rsem.txt` | cBioPortal brca_tcga | Yes |
| `data_cna.txt` | cBioPortal brca_tcga | Yes |

### METABRIC

Files are downloaded from cBioPortal (tarball). Placed in `data/metabric/`.

| File | Contents | Auto-downloaded |
|------|----------|-----------------|
| `data_clinical_patient.txt` | Patient demographics, survival | Yes |
| `data_clinical_sample.txt` | Tumor characteristics, receptor status | Yes |
| `data_mutations.txt` | Somatic mutations (~173-gene panel) | Yes |
| `data_mrna_illumina_microarray.txt` | Gene expression (Illumina HT-12 v3) | Yes |
| `data_cna.txt` | GISTIC2 copy number | Yes |

## Schema Overview

9 tables, 4 views, 27 indexes. See `schema/brca_harmonized_schema.sql` for full DDL.

```
patients  ──1:N──  tumors  ──1:N──  gene_expression
    │                 │               gene_expression_combat
    │                 │               copy_number
    └──1:N── treatments               mutations
                                      methylation (TCGA only)
                                      protein_expression (TCGA only)

harmonization_log  (audit trail for all transformations)
```

## Cross-Study Harmonization Notes

**Gene Expression**: TCGA is RNA-seq, METABRIC is microarray. Raw values are
not directly comparable. Use `gene_expression.expression_z_score` for rank-based
comparison, or `gene_expression_combat.expression_combat` for ComBat-corrected values.

**Mutations**: TCGA has whole-exome coverage (~20K genes), METABRIC has a
targeted ~173-gene panel. Use `WHERE in_metabric_panel = 1` for valid cross-study
mutation frequency comparisons. METABRIC coordinates are lifted from GRCh37 to
GRCh38 in step 6.

**Copy Number**: Both studies used Affymetrix SNP6.0 arrays. GISTIC2 discrete
values are directly comparable without additional normalization.

**Clinical**: Survival times standardized to months. TCGA days converted via
`days / 30.44`. Use PFI (progression-free interval) as the recommended endpoint
for BRCA per Liu et al. (Cell, 2018).

## Project Structure

```
brca_harmonized/
├── README.md
├── environment.yaml
├── requirements.txt
├── setup.py
├── .gitignore
├── schema/
│   └── brca_harmonized_schema.sql
├── config/
│   ├── __init__.py
│   ├── settings.py
│   └── metabric_panel_genes.py
├── etl/
│   ├── __init__.py
│   ├── _bootstrap.py
│   ├── db.py
│   ├── download_data.py
│   ├── harmonize.py
│   ├── init_db.py
│   ├── load_clinical.py
│   ├── load_mutations.py
│   ├── load_expression.py
│   ├── load_copy_number.py
│   ├── liftover_mutations.py
│   ├── combat_expression.py
│   └── run_all.py
├── tests/
│   ├── __init__.py
│   └── test_harmonize.py
└── data/
    ├── tcga/           # Place TCGA downloads here
    └── metabric/       # Place METABRIC downloads here
```

## Running Tests

```bash
conda activate brca_analysis
python -m pytest tests/ -v
```

## References

- Curtis C, et al. (2012) The genomic and transcriptomic architecture of 2,000
  breast tumours reveals novel subgroups. Nature 486:346-52.
- Pereira B, et al. (2016) The somatic mutation profiles of 2,433 breast cancers
  refines their genomic and transcriptomic landscapes. Nat Commun 7:11479.
- Liu J, et al. (2018) An integrated TCGA Pan-Cancer clinical data resource to
  drive high-quality survival outcome analytics. Cell 173:400-416.
- Johnson WE, et al. (2007) Adjusting batch effects in microarray expression data
  using empirical Bayes methods. Biostatistics 8:118-127.
- Behdenna A, et al. (2023) pyComBat, a Python tool for batch effects correction.
  BMC Bioinformatics 24:459.
