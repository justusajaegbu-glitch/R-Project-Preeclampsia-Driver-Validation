# R-Project
Single-cell transcriptomics (Seurat) and machine learning integration of placental and PBMC gene expression.

# Do Random Forest-Ranked Drivers of Preeclampsia Replicate in an Independent Decidua Cohort?

## 1. Objective

**Question:** Do the genes ranked highest by a Random Forest on integrated placenta and PBMC single-cell data (RPS24, RPLP1, RPL35, APRT, SP100, DYNLT1) change in the same direction in independent preeclamptic decidua basalis (GSE60438; 42 controls vs 35 preeclampsia), and do they diagnose preeclampsia beyond gestational age?

**Prediction (stated before validation):** RPS24, RPLP1, RPL35 and APRT, which are higher in preeclampsia in the discovery data, will be higher in preeclamptic decidua, and SP100 and DYNLT1, which are lower, will be lower (limma, adjusted for gestational age, maternal age and fetal sex; FDR < 0.05). Together they will separate the groups (cross-validated AUC well above 0.5).

**Outcome:** Not supported. Only APRT differed, in the opposite direction, and the panel added no information beyond gestational age (Section 7).

## 2. Dataset

| Role | Source | Details |
|---|---|---|
| Discovery | Placenta, GSE226601 (Admati) | Single-cell, integrated with Zhou PBMC data |
| Discovery | PBMC, Zhou et al. | Pre-processed GitHub subsets (raw data CNP0003201 is access-restricted) |
| Validation | GSE60438 | Illumina BeadArray, decidua basalis, 77 samples: 42 normotensive controls, 35 pre-eclamptic |

The integrated Seurat object is `combined_study_light.rds`.

## 3. Tools

R 4.4.2, Seurat, Harmony, randomForest (500 trees, fixed seed), SHAP, GEOquery, limma, pROC, pheatmap. Power analysis uses base R `power.t.test`.

## 4. Workflow

1. **Integration:** QC and normalisation (Seurat), batch correction (Harmony), PCA/UMAP.
2. **Pseudobulk:** aggregated per condition; log2FC after library-size normalisation (CPM).
3. **Machine learning:** Random Forest ranked by Mean Decrease Gini; SHAP for interpretation.
4. **Validation (GSE60438):** probes collapsed to one per gene (highest mean); statistical tests and ROC-AUC below.

Flowchart: <img width="816" height="454" alt="image" src="https://github.com/user-attachments/assets/acac5177-1da6-4042-a7bb-f9c4cb7ad73b" />


## 5. Discovery results

Top Random Forest drivers (Mean Decrease Gini): RPS24 (0.012), RPLP1 (0.010), APRT (0.010), SP100, RPL35, RACK1 (about 0.008).

CPM-normalised pseudobulk log2FC (preeclampsia vs control):

| APRT | RPLP1 | RPS24 | RPL35 | ATP6V1G1 | HMGB1 | CCNI | SP100 | DYNLT1 |
|---|---|---|---|---|---|---|---|---|
| +0.48 | +0.33 | +0.30 | +0.11 | 0.00 | -0.02 | -0.15 | -0.92 | -1.18 |

## 6. Validation strategy (GSE60438)

- **Sample size:** 42 vs 35. At Bonferroni alpha = 0.05/9 = 0.0056, power reaches 80% only for |d| of 0.85 or more (0.65 at alpha = 0.05).
- **Per-gene tests:** Welch t-test and Wilcoxon, with BH and Bonferroni correction across the 9 genes.
- **Adjusted model:** limma on all genes, with gestational age, maternal age and fetal sex as covariates.
- **Gene-set test:** `camera` and `fry` ask whether the 9 drivers shift together.
- **Effect size:** Hedges' g with 95% CI, and adjusted log2FC.
- **Diagnostic performance (ROC-AUC):** per-gene AUC with DeLong 95% CI (above 0.5 means higher in preeclampsia); a pre-specified discovery score (each gene weighted by its discovery direction, nothing fitted); and a 9-gene logistic panel with repeated stratified 5-fold cross-validation (20 repeats), compared with gestational age plus maternal age.

## 7. Validation results

Full tables: `results/validation_stats_GSE60438.csv`, `results/roc_per_gene_GSE60438.csv`, `results/roc_panel_cv_GSE60438.csv`.

**Per gene** (AUC above 0.5 = higher in preeclampsia; FDR from Wilcoxon, BH):

| Gene | Function | Discovery log2FC | Validation AUC (95% CI) | FDR | Agreement |
|---|---|---|---|---|---|
| RPS24 | 40S ribosomal protein | +0.30 | 0.48 (0.34-0.61) | 0.81 | Not replicated |
| RPLP1 | Ribosomal stalk protein | +0.33 | 0.45 (0.32-0.58) | 0.80 | Not replicated |
| RPL35 | 60S ribosomal protein | +0.11 | 0.47 (0.33-0.60) | 0.80 | Not replicated |
| APRT | Purine salvage enzyme | +0.48 | **0.23 (0.12-0.34)** | **0.0004** | **Opposite direction** |
| SP100 | Interferon-inducible nuclear body protein | -0.92 | 0.58 (0.45-0.71) | 0.50 | Opposite direction, not significant |
| DYNLT1 | Dynein light chain, trafficking | -1.18 | 0.50 (0.37-0.63) | 1.00 | Not replicated |
| CCNI | Cyclin I | -0.15 | 0.53 (0.40-0.67) | 0.80 | Not replicated |
| HMGB1 | Alarmin (DAMP), TLR4/RAGE signalling | -0.02 | 0.42 (0.29-0.55) | 0.50 | No discovery effect |
| ATP6V1G1 | V-ATPase subunit | 0.00 | 0.59 (0.46-0.72) | 0.50 | No discovery effect |

APRT is the only gene that separates the groups. It is lower in preeclampsia (Hedges' g about -1.0; adjusted log2FC about -0.2), the opposite of its discovery direction. All other effects are small and their CIs include chance. Nine of the ten top drivers were tested (RACK1 was not in the validation set).

**Diagnostic performance:**

| Model | AUC |
|---|---|
| Discovery-direction score (no fitting) | 0.37 (95% CI 0.25-0.50) |
| 9-gene panel, cross-validated | 0.70 (range 0.65-0.72 over 20 repeats) |
| Gestational age + maternal age, cross-validated | 0.88 (0.84-0.91) |
| Gestational age + maternal age + 9-gene panel, cross-validated | 0.85 (0.77-0.92) |
| Gestational age alone (not cross-validated, lower in preeclampsia) | 0.92 |

The discovery signature does not predict preeclampsia in decidua; its point estimate is below chance. The fitted panel reaches a modest 0.70, mostly through APRT, and adds nothing to gestational age plus maternal age (0.85 vs 0.88).

## 8. Biological interpretation

- **Ribosomal genes (RPS24, RPLP1, RPL35):** very abundant transcripts, so small relative shifts can rank high in a classifier without being large biological effects. Their discovery changes were modest (+0.1 to +0.3), came from pooled counts with no replicates, and are absent in decidua (AUC 0.45-0.48).
- **APRT:** a purine salvage enzyme that supports nucleotide supply in metabolically active tissue. It is the one robust difference in decidua, but lower in preeclampsia, against higher in the discovery data. This may reflect tissue differences or gestational age, which is strongly confounded with disease in this cohort (AUC 0.92 alone).
- **SP100 and DYNLT1:** the largest discovery changes (about 2-fold down), with no decidual difference. SP100 sits in the interferon response, alongside IFI35 (a top Random Forest gene not validated here).
- **HMGB1:** an alarmin released from stressed cells, previously linked to preeclampsia. It showed no effect in either dataset.
- **Overall:** the discovery drivers do not carry over to decidua basalis. A systemic-to-placental link, a tissue mismatch (blood and placenta vs decidua) and low power for small effects cannot be separated with these data.

## 9. Limitations

1. Raw PBMC data was restricted, so pre-processed GitHub subsets were used.
2. Discovery effects come from pooled counts without replicates (no CIs). The control group had about 48 times more total counts (304M vs 6.4M), so values must be library-size normalised.
3. The validation cohort is a different tissue from the discovery data.
4. Gestational age is strongly confounded with disease in GSE60438, so expression-based diagnostic claims are uncertain.
5. AUC CIs are wide at 42 vs 35 samples, so the result is "no evidence of diagnostic value here", not "no diagnostic value".
6. Random Forest rankings vary between runs, so a fixed seed is required.
7. Ligand-receptor incompatibility prevented CellPhoneDB, so a centroid-based signalling analysis was used.



## 10. Reproducibility

Run the validation script inside `Validation.R/` using **R 4.4.2**[cite: 4, 9]. `GSE60438` is downloaded directly via `GEOquery::getGEO("GSE60438")`[cite: 4]. A fixed random seed is set prior to the Random Forest and cross-validation steps to ensure consistent rankings, and all output tables and plots are exported to `results/`[cite: 4, 5, 7].

```text
R-Project/Preeclampsia-Driver-Validation
├── README.md
├── Validation.R/        # R analysis and validation script
└── results/             # Exported ROC/CV CSV tables and PNG figures
