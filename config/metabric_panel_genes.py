"""
METABRIC targeted sequencing panel: ~173 cancer-related genes.

Source: Pereira et al., Nat Commun 2016 (Supplementary Table).
The somatic mutation profiles of 2,433 breast cancers.

NOTE: This hardcoded list is a fallback. The preferred approach is to
derive the panel gene set from the actual METABRIC mutation data file,
which contains 173 unique genes — only 59 of which overlap with this list.
See load_mutations.py for the data-driven derivation.
"""

METABRIC_PANEL_GENES = frozenset([
    "AKT1", "AKT2", "AKT3", "ARAF", "ARID1A", "ARID1B", "ARID2",
    "ATM", "ATR", "ATRX", "AXIN1", "BAP1", "BRAF", "BRCA1", "BRCA2",
    "CASP8", "CBFB", "CCND1", "CCNE1", "CDH1", "CDK12", "CDK4",
    "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CEBPA",
    "CHD4", "CHEK2", "CIC", "CREBBP", "CTCF", "CTNNB1", "CUL3",
    "DAXX", "DDR2", "DNMT3A", "EGFR", "EP300", "ERBB2", "ERBB3",
    "ERBB4", "ESR1", "EZH2", "FAM46C", "FANCA", "FANCC", "FANCD2",
    "FANCI", "FAT1", "FBXW7", "FGFR1", "FGFR2", "FGFR3", "FLT3",
    "FOXA1", "FOXL2", "FOXP1", "GATA3", "GNA11", "GNA13", "GNAQ",
    "GNAS", "GPS2", "HGF", "HIST1H3B", "HNF1A", "HRAS", "IDH1",
    "IDH2", "IGF1R", "JAK1", "JAK2", "JAK3", "KDM5C", "KDM6A",
    "KIT", "KMT2A", "KMT2C", "KMT2D", "KRAS", "MAP2K1", "MAP2K2",
    "MAP2K4", "MAP3K1", "MAP3K13", "MDM2", "MDM4", "MED12", "MEN1",
    "MET", "MLH1", "MRE11", "MSH2", "MSH6", "MTOR", "MYC", "MYCN",
    "MYD88", "NBN", "NCOR1", "NF1", "NF2", "NFE2L2", "NOTCH1",
    "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NRAS", "NSD1", "NTRK1",
    "NTRK3", "PALB2", "PBRM1", "PDGFRA", "PDGFRB", "PHF6", "PIK3CA",
    "PIK3CB", "PIK3R1", "PMS2", "POLD1", "POLE", "PPP2R1A", "PRDM1",
    "PTCH1", "PTEN", "PTPN11", "RAD21", "RAD50", "RAD51C", "RAF1",
    "RB1", "RHOA", "RICTOR", "RNF43", "RPTOR", "RUNX1", "SETD2",
    "SF3B1", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1",
    "SMO", "SOCS1", "SOX9", "SPOP", "SRC", "STAG2", "STK11",
    "TAF1", "TBX3", "TERT", "TET2", "TNFAIP3", "TP53", "TSC1",
    "TSC2", "U2AF1", "VHL", "WT1", "XBP1", "ZFP36L1",
])
