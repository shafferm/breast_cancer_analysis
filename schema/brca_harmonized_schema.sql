-- ============================================================================
-- BRCA Harmonized Schema: TCGA-BRCA + METABRIC
-- Target: SQLite 3.x (via SQLAlchemy Core)
-- Generated: 2026-03-25
-- Note: PRAGMAs (journal_mode, foreign_keys, synchronous) are set by the
--       engine connection event in etl/db.py, not in this file.
-- ============================================================================

-- ============================================================================
-- CORE CLINICAL
-- ============================================================================

CREATE TABLE IF NOT EXISTS patients (
    patient_id          TEXT    PRIMARY KEY,                     -- Prefixed by study (e.g. TCGA-A2-A0YF, MB-0001)
    study_source        TEXT    NOT NULL CHECK (study_source IN ('TCGA-BRCA', 'METABRIC')),
    age_at_diagnosis    REAL,                                   -- Age in years at initial diagnosis
    sex                 TEXT    CHECK (sex IN ('M', 'F')),
    race                TEXT,                                   -- TCGA only
    ethnicity           TEXT,                                   -- TCGA only
    menopausal_status   TEXT    CHECK (menopausal_status IN ('Pre', 'Post', 'Peri', 'Unknown')),
    vital_status        TEXT    CHECK (vital_status IN ('Alive', 'Deceased')),
    os_months           REAL,                                   -- Overall survival in months
    os_event            INTEGER CHECK (os_event IN (0, 1)),     -- 1 = death, 0 = censored
    dfs_months          REAL,                                   -- Disease-free / relapse-free survival in months
    dfs_event           INTEGER CHECK (dfs_event IN (0, 1)),    -- 1 = recurrence or death, 0 = censored
    created_at          TEXT    NOT NULL DEFAULT (datetime('now')),
    updated_at          TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_patients_study     ON patients (study_source);
CREATE INDEX idx_patients_vital     ON patients (vital_status);
CREATE INDEX idx_patients_pam50     ON patients (study_source, vital_status);


CREATE TABLE IF NOT EXISTS tumors (
    tumor_id                    TEXT    PRIMARY KEY,
    patient_id                  TEXT    NOT NULL REFERENCES patients (patient_id) ON DELETE CASCADE,
    histological_type           TEXT,                           -- IDC, ILC, mixed, mucinous, etc.
    tumor_grade                 INTEGER CHECK (tumor_grade IN (1, 2, 3)),
    tumor_stage                 TEXT,                           -- AJCC pathologic stage (I, IA, IB, II, IIA, IIB, III, IIIA, IIIB, IIIC, IV)
    tumor_size_cm               REAL,                           -- Largest dimension in cm
    lymph_nodes_positive        INTEGER,
    lymph_nodes_examined        INTEGER,
    er_status                   TEXT    CHECK (er_status IN ('Positive', 'Negative', 'Unknown')),
    pr_status                   TEXT    CHECK (pr_status IN ('Positive', 'Negative', 'Unknown')),
    her2_status                 TEXT    CHECK (her2_status IN ('Positive', 'Negative', 'Equivocal', 'Unknown')),
    tnbc_status                 INTEGER CHECK (tnbc_status IN (0, 1)),   -- Derived: 1 if ER-/PR-/HER2-
    pam50_subtype               TEXT    CHECK (pam50_subtype IN ('LumA', 'LumB', 'Her2', 'Basal', 'Normal', 'Unknown')),
    integrative_cluster         TEXT,                           -- IntClust 1-10 (native in METABRIC, imputable for TCGA)
    nottingham_prognostic_index REAL,                           -- NPI score (native in METABRIC, computable for TCGA)
    created_at                  TEXT    NOT NULL DEFAULT (datetime('now')),
    updated_at                  TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_tumors_patient     ON tumors (patient_id);
CREATE INDEX idx_tumors_er          ON tumors (er_status);
CREATE INDEX idx_tumors_pam50       ON tumors (pam50_subtype);
CREATE INDEX idx_tumors_subtype     ON tumors (er_status, pr_status, her2_status);
CREATE INDEX idx_tumors_intclust    ON tumors (integrative_cluster);


CREATE TABLE IF NOT EXISTS treatments (
    treatment_id    TEXT    PRIMARY KEY,
    patient_id      TEXT    NOT NULL REFERENCES patients (patient_id) ON DELETE CASCADE,
    treatment_type  TEXT    NOT NULL CHECK (treatment_type IN ('Chemotherapy', 'Hormone Therapy', 'Radiation', 'Surgery', 'Targeted Therapy', 'Other')),
    treatment_name  TEXT,                                       -- Drug name or regimen (TCGA has detail, METABRIC typically NULL)
    received        INTEGER NOT NULL CHECK (received IN (0, 1)),
    created_at      TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_treatments_patient ON treatments (patient_id);
CREATE INDEX idx_treatments_type    ON treatments (treatment_type);


-- ============================================================================
-- GENOMIC / MOLECULAR (BOTH STUDIES)
-- ============================================================================

CREATE TABLE IF NOT EXISTS gene_expression (
    tumor_id            TEXT    NOT NULL REFERENCES tumors (tumor_id) ON DELETE CASCADE,
    gene_symbol         TEXT    NOT NULL,                       -- HUGO gene symbol
    entrez_id           INTEGER,                               -- NCBI Entrez gene ID
    expression_value    REAL,                                  -- Log2-transformed normalized expression
    expression_z_score  REAL,                                  -- Z-score relative to study normals/diploid
    platform            TEXT    NOT NULL CHECK (platform IN ('RNA-seq', 'microarray')),
    PRIMARY KEY (tumor_id, gene_symbol)
);

CREATE INDEX idx_expr_gene      ON gene_expression (gene_symbol);
CREATE INDEX idx_expr_tumor     ON gene_expression (tumor_id);
CREATE INDEX idx_expr_platform  ON gene_expression (platform);


CREATE TABLE IF NOT EXISTS copy_number (
    tumor_id        TEXT    NOT NULL REFERENCES tumors (tumor_id) ON DELETE CASCADE,
    gene_symbol     TEXT    NOT NULL,
    gistic_value    INTEGER CHECK (gistic_value IN (-2, -1, 0, 1, 2)),   -- homdel / hetloss / diploid / gain / amp
    log2_ratio      REAL,                                                 -- Continuous log2 copy number ratio
    PRIMARY KEY (tumor_id, gene_symbol)
);

CREATE INDEX idx_cn_gene    ON copy_number (gene_symbol);
CREATE INDEX idx_cn_tumor   ON copy_number (tumor_id);
CREATE INDEX idx_cn_gistic  ON copy_number (gistic_value);


CREATE TABLE IF NOT EXISTS mutations (
    mutation_id             TEXT    PRIMARY KEY,
    tumor_id                TEXT    NOT NULL REFERENCES tumors (tumor_id) ON DELETE CASCADE,
    gene_symbol             TEXT    NOT NULL,
    chromosome              TEXT,
    start_position          INTEGER,                           -- Genomic start position (harmonized to GRCh38)
    end_position            INTEGER,
    reference_genome        TEXT    NOT NULL DEFAULT 'GRCh38',
    reference_allele        TEXT,
    variant_allele          TEXT,
    variant_classification  TEXT,                               -- Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, etc.
    variant_type            TEXT    CHECK (variant_type IN ('SNP', 'INS', 'DEL', 'DNP', 'TNP', 'ONP')),
    hgvsp                   TEXT,                               -- Protein change in HGVSp notation
    vaf                     REAL,                               -- Variant allele frequency
    sequencing_scope        TEXT    NOT NULL CHECK (sequencing_scope IN ('WES', 'PANEL')),
    in_metabric_panel       INTEGER NOT NULL DEFAULT 0 CHECK (in_metabric_panel IN (0, 1))   -- Flag for cross-study comparability
);

CREATE INDEX idx_mut_tumor      ON mutations (tumor_id);
CREATE INDEX idx_mut_gene       ON mutations (gene_symbol);
CREATE INDEX idx_mut_class      ON mutations (variant_classification);
CREATE INDEX idx_mut_panel      ON mutations (in_metabric_panel);
CREATE INDEX idx_mut_chr_pos    ON mutations (chromosome, start_position);


-- ============================================================================
-- TCGA-ONLY DATA TYPES
-- ============================================================================

CREATE TABLE IF NOT EXISTS methylation (
    tumor_id        TEXT    NOT NULL REFERENCES tumors (tumor_id) ON DELETE CASCADE,
    probe_id        TEXT    NOT NULL,                           -- Illumina probe ID (cg...)
    gene_symbol     TEXT,                                      -- Nearest gene
    beta_value      REAL    CHECK (beta_value >= 0.0 AND beta_value <= 1.0),
    chromosome      TEXT,
    position        INTEGER,
    PRIMARY KEY (tumor_id, probe_id)
);

CREATE INDEX idx_meth_gene  ON methylation (gene_symbol);
CREATE INDEX idx_meth_tumor ON methylation (tumor_id);


CREATE TABLE IF NOT EXISTS protein_expression (
    tumor_id            TEXT    NOT NULL REFERENCES tumors (tumor_id) ON DELETE CASCADE,
    antibody_id         TEXT    NOT NULL,
    protein_name        TEXT,
    expression_value    REAL,                                  -- Normalized RPPA expression
    PRIMARY KEY (tumor_id, antibody_id)
);

CREATE INDEX idx_rppa_tumor     ON protein_expression (tumor_id);
CREATE INDEX idx_rppa_protein   ON protein_expression (protein_name);


-- ============================================================================
-- METADATA / AUDIT
-- ============================================================================

CREATE TABLE IF NOT EXISTS harmonization_log (
    log_id          INTEGER PRIMARY KEY AUTOINCREMENT,
    table_name      TEXT    NOT NULL,
    study_source    TEXT    NOT NULL CHECK (study_source IN ('TCGA-BRCA', 'METABRIC', 'BOTH')),
    field_name      TEXT    NOT NULL,
    transformation  TEXT    NOT NULL,                           -- Description of mapping/conversion logic
    source_field    TEXT,                                       -- Original field name in source data
    notes           TEXT,                                       -- Caveats, assumptions, data loss warnings
    created_at      TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_hlog_table ON harmonization_log (table_name);
CREATE INDEX idx_hlog_study ON harmonization_log (study_source);


-- ============================================================================
-- CONVENIENCE VIEWS
-- ============================================================================

-- Flat patient + tumor view for quick exploratory queries
CREATE VIEW IF NOT EXISTS v_patient_tumor AS
SELECT
    p.patient_id,
    p.study_source,
    p.age_at_diagnosis,
    p.sex,
    p.menopausal_status,
    p.vital_status,
    p.os_months,
    p.os_event,
    p.dfs_months,
    p.dfs_event,
    t.tumor_id,
    t.histological_type,
    t.tumor_grade,
    t.tumor_stage,
    t.tumor_size_cm,
    t.lymph_nodes_positive,
    t.er_status,
    t.pr_status,
    t.her2_status,
    t.tnbc_status,
    t.pam50_subtype,
    t.integrative_cluster,
    t.nottingham_prognostic_index
FROM patients p
JOIN tumors t ON p.patient_id = t.patient_id;


-- Cross-study mutation analysis restricted to shared gene panel
CREATE VIEW IF NOT EXISTS v_mutations_cross_study AS
SELECT
    m.mutation_id,
    p.study_source,
    m.tumor_id,
    m.gene_symbol,
    m.variant_classification,
    m.variant_type,
    m.hgvsp,
    m.vaf,
    m.sequencing_scope
FROM mutations m
JOIN tumors t ON m.tumor_id = t.tumor_id
JOIN patients p ON t.patient_id = p.patient_id
WHERE m.in_metabric_panel = 1;


-- Per-gene mutation frequency by study
CREATE VIEW IF NOT EXISTS v_mutation_frequency AS
SELECT
    p.study_source,
    m.gene_symbol,
    COUNT(DISTINCT m.tumor_id) AS n_mutated_samples,
    COUNT(DISTINCT t.tumor_id) AS n_total_samples_in_study,
    ROUND(
        CAST(COUNT(DISTINCT m.tumor_id) AS REAL) /
        CAST(COUNT(DISTINCT t.tumor_id) AS REAL) * 100, 2
    ) AS mutation_pct
FROM mutations m
JOIN tumors t ON m.tumor_id = t.tumor_id
JOIN patients p ON t.patient_id = p.patient_id
WHERE m.in_metabric_panel = 1
GROUP BY p.study_source, m.gene_symbol;


-- Treatment summary per patient
CREATE VIEW IF NOT EXISTS v_treatment_summary AS
SELECT
    patient_id,
    MAX(CASE WHEN treatment_type = 'Chemotherapy'      AND received = 1 THEN 1 ELSE 0 END) AS had_chemo,
    MAX(CASE WHEN treatment_type = 'Hormone Therapy'    AND received = 1 THEN 1 ELSE 0 END) AS had_hormone,
    MAX(CASE WHEN treatment_type = 'Radiation'          AND received = 1 THEN 1 ELSE 0 END) AS had_radiation,
    MAX(CASE WHEN treatment_type = 'Surgery'            AND received = 1 THEN 1 ELSE 0 END) AS had_surgery,
    MAX(CASE WHEN treatment_type = 'Targeted Therapy'   AND received = 1 THEN 1 ELSE 0 END) AS had_targeted
FROM treatments
GROUP BY patient_id;


-- ============================================================================
-- SEED HARMONIZATION LOG WITH KNOWN TRANSFORMATIONS
-- ============================================================================

INSERT INTO harmonization_log (table_name, study_source, field_name, transformation, source_field, notes) VALUES
    ('patients', 'TCGA-BRCA', 'os_months',
     'days_to_death / 30.44 (or days_to_last_follow_up / 30.44 if censored)',
     'demographic.days_to_death, demographic.days_to_last_follow_up',
     'Use TCGA-CDR (Liu et al., Cell 2018) curated endpoints. PFI recommended over OS for BRCA.'),

    ('patients', 'TCGA-BRCA', 'patient_id',
     'Retained TCGA barcode as-is (e.g. TCGA-A2-A0YF)',
     'submitter_id',
     'First 12 characters of TCGA barcode uniquely identify the patient.'),

    ('patients', 'METABRIC', 'patient_id',
     'Prefixed with MB- to match PATIENT_ID from cBioPortal clinical file',
     'PATIENT_ID',
     'Original METABRIC IDs like MB-0001 are already prefixed.'),

    ('patients', 'TCGA-BRCA', 'menopausal_status',
     'Inferred from age: Pre if age < 50, Post if age >= 50',
     'age_at_diagnosis',
     'This is an approximation. METABRIC records menopausal status directly.'),

    ('patients', 'METABRIC', 'race',
     'Set to NULL - field not available in METABRIC',
     NULL,
     'METABRIC cohort is predominantly UK/Canadian, individual-level race/ethnicity not collected.'),

    ('tumors', 'TCGA-BRCA', 'tumor_id',
     'Used TCGA sample barcode (first 16 characters)',
     'submitter_id (sample-level)',
     'Format: TCGA-XX-XXXX-01A. The -01 suffix indicates primary solid tumor.'),

    ('tumors', 'METABRIC', 'nottingham_prognostic_index',
     'Used NPI value directly from cBioPortal clinical file',
     'NOTTINGHAM_PROGNOSTIC_INDEX',
     'For TCGA, NPI can be computed: (0.2 * tumor_size_cm) + tumor_grade + lymph_node_score.'),

    ('tumors', 'BOTH', 'tnbc_status',
     'Derived: 1 if er_status = Negative AND pr_status = Negative AND her2_status = Negative',
     'er_status, pr_status, her2_status',
     'Applied consistently to both studies after receptor status harmonization.'),

    ('gene_expression', 'TCGA-BRCA', 'expression_value',
     'Log2(FPKM + 1) from GDC harmonized RNA-seq pipeline',
     'htseq_fpkm files',
     'TCGA values are RNA-seq based. NOT directly comparable to METABRIC microarray values.'),

    ('gene_expression', 'METABRIC', 'expression_value',
     'Log2 intensity values from Illumina HT-12 v3 arrays as provided by cBioPortal',
     'data_mRNA_median_Zscores.txt / raw expression file',
     'Microarray platform. Use z-scores or ComBat for cross-study comparison.'),

    ('gene_expression', 'BOTH', 'gene_symbol',
     'Mapped to HUGO gene symbols using GENCODE v36 annotation',
     'various probe/gene mappings',
     'METABRIC Illumina probe-to-gene mappings may need updating. Ambiguous probes excluded.'),

    ('mutations', 'TCGA-BRCA', 'start_position',
     'Coordinates already in GRCh38 from GDC harmonized MAF files',
     'Start_Position (MAF)',
     'GDC re-aligned all TCGA data to GRCh38.'),

    ('mutations', 'METABRIC', 'start_position',
     'Lifted over from GRCh37 to GRCh38 using UCSC liftOver or CrossMap',
     'Start_Position (cBioPortal mutation file)',
     'Original METABRIC coordinates are GRCh37. Unmappable variants flagged and excluded.'),

    ('mutations', 'METABRIC', 'in_metabric_panel',
     'Set to 1 for all METABRIC mutations (all are within the targeted panel by definition)',
     NULL,
     'For TCGA, in_metabric_panel = 1 only if the gene appears in the METABRIC ~173-gene panel list.'),

    ('mutations', 'BOTH', 'variant_classification',
     'Standardized to MAF specification vocabulary: Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, etc.',
     'Variant_Classification',
     'METABRIC cBioPortal files use same MAF vocabulary. Verified for consistency.'),

    ('copy_number', 'BOTH', 'gistic_value',
     'Discrete CNA values (-2 to +2) from GISTIC2 output',
     'data_CNA.txt (cBioPortal) / focal_data_by_genes.txt (GDC)',
     'Both studies used Affymetrix SNP6.0 arrays. Most directly comparable data type across studies.');
