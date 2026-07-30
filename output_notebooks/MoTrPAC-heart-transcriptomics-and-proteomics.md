# Heart Transcriptomic & Proteomic Analysis

This workflow adapts the WAT-SC analytical framework to heart transcriptomics and proteomics data across six sequential steps: 

1. Data audit and PCA
2. Differential expression 
3. Pathway enrichment (FGSEA) 
4. Weighted gene co-expression network analysis (WGCNA) 
5. WGCNA module functional enrichment (ORA)
6. Global protein differential analysis. 

In this markdown file, the folder and file uses path in this [Google Drive](https://drive.google.com/drive/folders/17yrQChPWkv6SgnYBCuWPmEEZtFaZSgLT?usp=sharing)

Each step reads its inputs from official MoTrPAC package data and/or the output folder of the preceding step, and writes a numbered `results/0N_.../` folder containing processed tables, figures, and a saved R object for downstream reuse. Official MoTrPAC helper functions (`run_WGCNA.R`, `rank_genes.R`, `fgsea2.R`, `fora2.R`) are downloaded fresh from the published GitHub repositories and sourced unchanged, with SHA-256 checksums logged for provenance; only project-specific adaptations (e.g., tissue swap, top-5,000 variance prefilter) are layered on top.

**Packages used across this workflow**

| Category | Packages |
|---|---|
| Official MoTrPAC data | `MotrpacRatTraining6moData` (installed via `remotes::install_github("MoTrPAC/MotrpacRatTraining6moData")`, data-only, no vignettes) |
| Official MoTrPAC helper code (downloaded, not installed) | `run_WGCNA.R`, `rank_genes.R`, `fgsea2.R`, `fora2.R` from `MoTrPAC/MotrpacRatTraining6moWAT`; `msigdb_pathways.R`, `WGCNA_MODULE_ORA.Rmd` from `MoTrPAC/MotrpacRatTraining6moWATData`; `proteomics_differential_analysis.R` from `MoTrPAC/MotrpacRatTraining6mo` |
| Differential expression / statistics | `edgeR`, `limma` |
| Pathway enrichment | `fgsea`, `msigdbr` |
| Co-expression network | `WGCNA`, `dynamicTreeCut`, `Biobase` |
| Data wrangling | `dplyr`, `tidyr`, `tibble`, `readr`, `purrr`, `stringr`, `forcats`, `data.table` |
| Plotting | `ggplot2`, `ggrepel`, `scales` |
| Provenance / utility | `digest` (SHA-256 checksums of downloaded source files), `assertthat`, `remotes` |

> [!NOTE]
> The complete `MotrpacRatTraining6mo` analysis package was deliberately **not** installed for the proteomics step, because its `proteomics_training_da()` function imports the CRAN package `metap`, which is unrelated to the timewise analysis used here and can fail to install in Colab. Only the one required source file was downloaded directly.

<br>

## Step 1 — Data audit & PCA (`results/01_rat_heart_data_audit/`)
Confirms sample balance (5 rats × 2 sexes × 5 training groups = 50 heart RNA-seq samples), checks for official outliers, and runs a PCA on the top 5,000 variable normalized transcripts.

- **Inputs:** official objects `PHENO`, `OUTLIERS`, `TRNSCRPT_META`, `TRNSCRPT_HEART_NORM_DATA`
- **Key outputs:** `dataset_inventory.csv`, `sample_counts_by_group.csv`, `heart_rnaseq_sample_metadata.csv`, `heart_rnaseq_pca_coordinates.csv`, `rat_heart_transcriptomics_audit_objects.rds`
- **Key plot:** `heart_rnaseq_pca.png` — shows near-complete separation of male and female samples along PC1 (11.8%), confirming sex as the dominant axis of transcriptomic variation before any modeling.

## Step 2 — Differential expression (`results/02_rat_heart_differential_expression/`)
edgeR `filterByExpr` + TMM normalization, then limma `voomWithQualityWeights` across three contrast families: baseline sex (sedentary M vs. F), within-sex training time course (1/2/4/8 weeks vs. sedentary, each sex), and sex-by-training interaction.

- **Inputs:** official objects `PHENO`, `OUTLIERS`, `TRNSCRPT_META`, `TRNSCRPT_HEART_RAW_COUNTS`, `TRNSCRPT_HEART_NORM_DATA`, `FEATURE_TO_GENE`
- **Key outputs:** **`all_differential_expression_results.csv`** (consumed by Steps 3 & 5), `significant_gene_counts.csv`, `top20_genes_each_contrast.csv`, `female_vs_male_8w_response.csv`
- **Key plots:** `voom_mean_variance.png`, `sample_quality_weights.png`, `volcano_M_SED_F_SED.png` (870 significant transcripts at baseline), `volcano_F_8W_F_SED.png`, `volcano_M_8W_M_SED.png`, `significant_transcript_counts.png`, `female_vs_male_8w_response.png` (ρ = 0.221 concordance)

## Step 3 — Pathway enrichment / FGSEA (`results/03_rat_heart_GitHub_FGSEA/`)
Runs the official `rank_genes()` and `fgsea2()` wrappers unchanged against MSigDB v7.5.1 GO gene sets (GO:BP/MF/CC), with the official ≥85% pathway-retention and ≥15-gene filters, globally BH-adjusted within each contrast family.

- **Inputs:** Step 2's `all_differential_expression_results.csv`
- **Key outputs:** **`filtered_MSigDB_GO_pathways_85pct.csv`** and **`all_GitHub_aligned_FGSEA_GO_results.csv`** (both consumed by Step 5), `significant_GO_pathway_counts.csv`, `top20_GO_pathways_each_contrast.csv`, `cardiac_related_GO_pathways.csv`
- **Key plots:** `baseline_sex_top_GO_BP.png` (127 significant GO terms at baseline), `female_8w_top_GO_BP.png`, `male_8w_top_GO_BP.png`, `interaction_8w_top_GO_BP.png`, `significant_GO_pathway_counts.png`, `training_timecourse_top_pathways_heatmap.png`, `cardiac_related_pathways_heatmap.png`

## Step 4 — WGCNA (`results/04_rat_heart_GitHub_WGCNA_top5000/`)
Selects the 5,000 highest-variance normalized heart transcripts (of ~14,000 detected), then calls the official `run_WGCNA()` function unchanged (signed network, bicor, power selected at R² ≥ 0.90, unsigned TOM, dynamic tree cut, module merge height 0.15).

- **Inputs:** official objects `PHENO`, `OUTLIERS`, `TRNSCRPT_META`, `TRNSCRPT_HEART_NORM_DATA`, `FEATURE_TO_GENE`
- **Key outputs:** **`official_WGCNA_modules.csv`** and **`significant_module_trait_correlations.csv`** (both consumed by Step 5), `module_eigengenes_wide.csv`, `top10_hub_genes_per_module.csv`, `module_contrast_results.csv`
- **Key plots:** `official_run_WGCNA_diagnostics.pdf` (+ 4 page PNGs: soft-threshold selection, gene clustering, dendrogram/module colors, eigengene clustering), `module_trait_correlation_heatmap.png`, `module_eigengene_group_heatmap.png`
- **Result:** 11 co-expression modules; T1/turquoise and T2/blue show the strongest opposing sex associations (ρ = ∓0.866, FDR = 2.0×10⁻¹⁴)

## Step 5 — WGCNA module functional enrichment / ORA (`results/05_rat_heart_GitHub_WGCNA_module_ORA/`)
Intersects Step 3's already-filtered GO pathways with the Step 4 top-5,000 WGCNA gene universe (≥15 genes post-intersection) and runs the official `fora2()` scaled over-representation test per module.

- **Inputs:** Step 2's `all_differential_expression_results.csv`; Step 3's `filtered_MSigDB_GO_pathways_85pct.csv` and `all_GitHub_aligned_FGSEA_GO_results.csv`; Step 4's `official_WGCNA_modules.csv`
- **Key outputs:** `all_module_GO_ORA_results.csv`, `significant_module_GO_ORA_results.csv` (107 significant module–GO associations), `module_ORA_significant_counts.csv`, `module_tables/module_T*_entrez_genes.csv`
- **Key plots:** `module_ORA_significant_counts.png`, `module_ORA_top_GO_BP_heatmap.png`, `module_ORA_top_terms_dotplot.png`, `T1_T2_functional_comparison.png`, `selected_modules_ORA_heatmap.png`
- **Result:** T1 → DNA repair/cell-cycle; T2 → ATP/nucleotide metabolism; T4 → immune activation; T6 → ECM organization; T7 → vascular-remodeling signaling; T8 → mitochondrial respiration

## Step 6 — Global proteomics (`results/06_rat_heart_GitHub_global_proteomics/`)
Reproduces the official 8 sex-stratified `proteomics_timewise_da()` contrasts (validated directly against `PROT_HEART_DA`), then extends to the same 13-contrast structure used in Step 2 for direct RNA–protein comparison.

- **Inputs:** official objects `PROT_HEART_NORM_DATA`, `PROT_HEART_DA`, `PHENO`, `OUTLIERS`, `FEATURE_TO_GENE_FILT`
- **Key outputs:** `official_timewise_protein_DA_reproduced.csv` vs. `official_PROT_HEART_DA_reference.csv` → `official_vs_reproduced_DA_comparison.csv` (validation that our reproduction matches the published pipeline), `all_13_contrast_protein_results.csv`, `significant_proteins.csv` (1,249 significant at baseline: 633 male-high / 616 female-high), `protein_significant_counts.csv`, `female_male_8w_protein_logFC_comparison.csv`
- **Key plots:** `protein_PCA.png`, `protein_sample_correlation_heatmap.png`, `official_vs_reproduced_logFC.png`, `baseline_sex_protein_volcano.png`, `female_8w_protein_volcano.png`, `male_8w_protein_volcano.png`, `female_male_8w_protein_logFC_comparison.png` (ρ = 0.24 concordance)

<br>

**Cross-step summary of significant results**

| Step | Layer | Baseline sex-difference (sedentary M vs. F) | Largest training effect |
|---|---|---|---|
| 2 | Transcripts | 870 significant (575 ↑male / 295 ↑female) | 256 transcripts, male 8-week |
| 3 | GO pathways | 127 significant (75 ↑male / 52 ↑female) | Female peak 4wk (114 pathways) / male peak 2wk (105 pathways) |
| 4 | Co-expression modules | 2 of 11 modules significant | 0 modules reach significance for training or interaction |
| 5 | Module–GO terms | 107 significant module–GO associations (non-directional) | — |
| 6 | Proteins | 1,249 significant (633 ↑male / 616 ↑female) | 424 proteins, male 8-week |

- **Scripts / notebooks (each `.R` script has a matching Colab/Drive `.ipynb`):** `01_rat_heart_data_audit`, `02_rat_heart_differential_expression.R`, `03_rat_heart_GitHub_aligned_FGSEA_FINAL.R`, `04_rat_heart_GitHub_WGCNA_top5000_FIXED.R`, `05_rat_heart_GitHub_WGCNA_module_ORA_FIXED.R`, `06_rat_heart_GitHub_global_proteomics_NO_METAP_FIXED.R`
- **processed results and figures:** `results/0N_.../` folders (see per-step tables above); mirrors the repo's `data_processed/` + `plots/` convention used elsewhere — consider relocating or symlinking into those folders for consistency once finalized
- **Drive link (all `results/` folders, figures, and `.rds` objects):** https://drive.google.com/drive/folders/17yrQChPWkv6SgnYBCuWPmEEZtFaZSgLT?usp=drive_link