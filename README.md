# DECODE-2026-Team-MOTRPAC

## Project Snapshot 

This project adapts and extends the core analytical framework used to identify sex-specific molecular differences in the published MoTrPAC rat endurance-training study of subcutaneous white adipose tissue (Many et al, 2024). Using heart tissue from the same parent study, we apply comparable methods to examine baseline male-versus-female differences and sex-specific responses to endurance exercise training.

The planned workflow includes matched-sample review, differential analysis, and pathway enrichment using heart RNA-sequencing and global proteomics data. Additional analyses, including WGCNA, phenotype correlations, RNA–protein integration, and cardiac-disease pathway mapping, may be completed as time and data availability permit.

By applying a published sexual-dimorphism and exercise-training analysis framework to a different tissue, the project will evaluate whether cardiac gene expression, protein abundance, and biological pathways show shared or tissue-specific patterns. It will also produce a transparent and reusable workflow that can support technical research and broader scientific interpretation.

Potential users include exercise and cardiovascular researchers, computational biologists, MoTrPAC data users, trainees, and others interested in reproducible multi-omic analysis and sex as a biological variable.

## Research Question

### Problem Statement 
Endurance exercise training (ExT) benefits cardiovascular health, but sex-specific molecular mechanisms driving cardiac adaptation remain poorly characterized. MoTrPAC generated the largest public multi-omic resource on this question, profiling 18 solid tissues and blood from male and female rats undergoing endurance training across four time points. Applying this pipeline to adipose tissue revealed profound sexual dimorphism at rest and in the dynamic response to training, with sex-specific signatures preserved or amplified with exercise. No independent replication of this sexual-dimorphism pipeline has been applied to heart tissue, despite matching data being available. 

### Questions

Primary Research Question: 
- Does sedentary heart tissue show sex-specific differences in gene expression and enriched pathways, similar to adipose tissue?

Secondary Research Question: 
- How does the cardiac transcriptomic/proteomic response to ExT differ between males and females across the 1/2/4/8-week time course?

### Objectives

- Reproduce the baseline male-versus-female analytical workflow reported for subcutaneous white adipose tissue from Many et al. (2024).
- Apply the same workflow to matched heart RNA-sequencing and global proteomics data.
- Identify baseline sex differences in cardiac gene expression, protein abundance, and pathway enrichment.
- Characterize within-sex molecular responses across the endurance-training time course.

## Workflow 

<p align="center">
  <img src="images/Pipeline.svg" alt="MoTrPAC analysis workflow" width="250">
</p>

## Set Up

### Background Review

- Read the published paper and relevant supplementary methods. https://doi.org/10.1038/s42255-023-00959-9
- Review the original published GitHub repositories, R package vignettes, and analysis scripts. 

### Software and Repository Setup

- Install the required R packages.
- Python is used for selected analysis steps. These scripts may be converted to R to support a consistent analysis environment.
- The published analysis code is available through two related MoTrPAC R packages. Scripts covering data processing and statistical analysis are provided in the `data-raw/` directory and in the package vignettes and articles for `MotrpacRatTraining6moWATData`. Additional helper functions used for analysis and figure generation are available in the companion `MotrpacRatTraining6moWAT` package. The paper’s Supplementary Methods also provides a summary of the main R and Bioconductor packages used in the workflow.

### Data Preparation

- Download the required MoTrPAC phenotype, transcriptomics, and proteomics files. See the raw_files folder and the file descriptions below for the specific inputs used in this project.
- Confirm filenames, formats, and expected folder locations.
- Review the phenotype data dictionary and sample identifiers.
- Match samples across transcriptomic and proteomic datasets.
- Confirm that the input files are normalized and ready for analysis.

### Validation

- Update file paths for the local or cloud environment.
- Confirm that all required packages load successfully.
- Test that each input file can be read and contains the expected columns.
- Run a small validation step using the published workflow.
- Document software versions, path changes, exclusions, and other deviations from the published analysis.
  
### Raw Data Files and Data Dictionaries (data_raw)

- **Phenotypic data file:** `motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt`<br>
  **Description:** This file contains the full set of phenotypic variables collected across study animals and biospecimens. Variables include information related to animal testing, study registration, calculated measures, and specimen processing, as defined in the accompanying data dictionary. The documentation first highlights the key variables needed for basic analyses and then provides additional detail.     
  **Documentation:** [MoTrPAC Phenotype Technical Guide](https://motrpac-data.org/technical-guides/phenotype)

- **Phenotypic data dictionary:** `motrpac_pass1b-06_pheno_dictionary_merged_v4.0.txt`<br>
  **Description:** This data dictionary defines the fields included in the phenotypic data tables, along with their scope and allowable values. It is intended to make the phenotypic data easier to interpret and use in analysis.   
**Documentation:** [MoTrPAC Phenotype Technical Guide](https://motrpac-data.org/technical-guides/phenotype)

- **Processed transcriptomic data file:** `motrpac_pass1b-06_t58-heart_transcript-rna-seq_normalized-log-cpm.txt`<br>
  **Description:** This file contains normalized RNA-sequencing measurements from rat heart tissue in the MoTrPAC PASS1B-06 endurance-training study. Gene expression is reported as log-transformed counts per million, which accounts for differences in sequencing depth across samples and reduces the influence of very large count values. The table provides sample-level gene-expression data that can be used to compare expression patterns across training groups, time points, and sexes.

- **Analysis output:** `pass1b-06_t58-heart_transcript-rna-seq_training-dea-fdr_20211008.txt`<br>
  **Description:** This file contains statistical results from a time-wise differential gene expression analysis of RNA-sequencing data generated from rat heart tissue in the MoTrPAC PASS1B-06 endurance-training study.

### Ontology for Pathway Enrichment and Interpretation Step (data_msigdb)

-  **Ontology:** `msigdb_go_basic_2021-12-15.obo`<br>
   **Description:** This file contains the December 15, 2021 Gene Ontology structure used to define and organize the biological-process, molecular-function, and cellular-component gene sets included in the pathway-enrichment analysis.

## Methods and Provenance

### Team Approach 

Three team members contributed to different components of the analysis workflow. Individual names may appear in code filenames, plot filenames, and output directories to indicate the analysis component each person developed or completed.

#### Differential Analysis - Jade

Jade contributed to the replication and tissue-adapted differential-analysis workflow. She reviewed and followed the published WAT-SC transcriptomics code, including the data-processing, filtering, differential-analysis, and volcano-plot workflow. Her work included completing the baseline male-versus-female transcriptomic and proteomic analyses associated with Figure 2a and preparing the resulting differential-analysis outputs for downstream pathway enrichment.

#### Pathway Enrichment and Cross-Omic Comparison - Viveka

Viveka contributed to the article and methodology review and to the pathway-enrichment and cross-omic comparison components of the workflow. Her work included analyses supporting Figures 2b and 2c, which compare transcriptomic and proteomic pathway-enrichment results and identify pathways shared across the two omics layers. Viveka also created the dashboard. 

#### WGCNA Analysis — Upama

Upama led the heart transcriptomics data audit and WGCNA analysis. This work included sample and feature review, selection of the 5,000 most variable features, co-expression module construction, module–trait analysis, hub-gene identification, and Gene Ontology over-representation analysis.

### Analysis Methods

#### Differential Analysis Approach

The differential-analysis workflow began with the published WAT-SC transcriptomics pipeline and traced the analysis from the original molecular data through sample and feature preparation, statistical testing, and volcano-plot generation. The published workflow was first reproduced using subcutaneous white adipose tissue data and then adapted to the corresponding heart transcriptomics and proteomics datasets.

Baseline male-versus-female comparisons were performed using sedentary samples. Separate transcriptomic and proteomic analyses identified genes and proteins that differed between males and females. The resulting ranked differential-expression and differential-abundance results were used to generate Figure 2a-style volcano plots and served as inputs for downstream pathway-enrichment analysis.

Minor differences were observed between the reproduced WAT-SC plots and the published figures. The analysis is therefore described as a reproduction and tissue adaptation of the documented analytical workflow rather than an exact recreation of every published result.

Notes: 
#### Figure 2a Replication — WAT-SC Transcriptomics

To reproduce the WAT-SC transcriptomics analysis used for Figure 2a, the published GitHub workflow was followed from the original count data through differential analysis and volcano-plot generation.

1. **Load the source transcriptomics data**  
   The analysis begins with `TRNSCRPT_WATSC_RAW_COUNTS` from the `MotrpacRatTraining6moData` R package.

2. **Prepare the expression data**  
   The script `data-raw/TRNSCRPT_EXP.R` is used to prepare the transcriptomics dataset. This step includes:
   - sample metadata and gene annotation;
   - expression filtering;
   - outlier removal; and
   - covariate processing.

3. **Create the processed expression object**  
   The resulting filtered raw counts and associated metadata are saved as `data/TRNSCRPT_EXP.rda`.

4. **Run differential analysis**  
   The vignette `vignettes/WAT_DA.Rmd` performs the differential-analysis workflow. This analysis calls the `limma_full()` function from the companion MoTrPAC analysis package.

5. **Save the differential-analysis results**  
   The differential-analysis output is saved as `data/TRNSCRPT_DA.rda`.

6. **Generate the volcano plots**  
   The differential-analysis results are used by `vignettes/articles/volcano_plots.Rmd` to generate the WAT-SC transcriptomics volcano plots associated with Figure 2a.



**Primary inputs:**

- [Heart normalized transcriptomics data](raw_files/motrpac_pass1b-06_t58-heart_transcript-rna-seq_normalized-log-cpm.txt)
- [Phenotype data](raw_files/motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt)
- [Phenotype data dictionary](raw_files/motrpac_pass1b-06_pheno_dictionary_merged_v4.0.txt)

**Processed outputs:**

- [Heart transcriptomics differential-expression results](data_processed/heart_transcriptomics_da.csv)
- [Heart proteomics differential-abundance results](data_processed/heart_proteomics_da.csv)
- [WAT-SC transcriptomics differential-expression results](data_processed/watsc_transcriptomics_da.csv)
- [WAT-SC proteomics differential-abundance results](data_processed/watsc_proteomics_da.csv)
- [Liver transcriptomics differential-expression results](data_processed/liver_transcriptomics_da.csv)

**Plots:**

The plots labeled `WATSC` reproduce the published findings for subcutaneous white adipose tissue shown in Figures 2a–2c. The liver and heart plots were generated by applying the same analytical methods to those tissues, allowing the team to compare the published adipose-tissue findings with tissue-adapted results.

-  Figure 2a: Differential analysis identifies individual transcripts and proteins that differ between sedentary males and females.
-  Figure 2b: Pathway enrichment compares whether transcriptomic and proteomic results point to the same sex-associated biological processes.
-  Figure 2c: Cross-omic integration highlights pathways that are significantly enriched in both the RNA and protein analyses.

Naming Convention Note: 
-  WATSC: Replication of the published subcutaneous white adipose tissue analysis.
-  Liver: Tissue-adapted analysis using the same methods applied to liver data.
-  Heart: Tissue-adapted analysis using the same methods applied to heart data.


### Pathway Enrichment and Cross-Omic Comparison Approach

The pathway-enrichment workflow began with ranked transcriptomic and proteomic results from the baseline male-versus-female differential analyses. Fast Gene Set Enrichment Analysis was used to determine whether genes or proteins associated with sex differences were concentrated within defined biological pathways.

Transcriptomic and proteomic enrichment results were compared to assess whether the two molecular layers showed similar pathway direction and strength. This supported the Figure 2b-style comparison. Pathways meeting the selected significance criteria in both omics layers were then identified to support the Figure 2c-style shared-pathway analysis.

The resulting enrichment tables and saved R objects were used for pathway comparison, visualization, and biological interpretation.

**Inputs from differential analysis:**

- [Heart transcriptomics differential-expression results](data_processed/heart_transcriptomics_da.csv)
- [Heart proteomics differential-abundance results](data_processed/heart_proteomics_da.csv)
- [WAT-SC transcriptomics differential-expression results](data_processed/watsc_transcriptomics_da.csv)
- [WAT-SC proteomics differential-abundance results](data_processed/watsc_proteomics_da.csv)

**Pathway reference:**

- [Gene Ontology reference file](raw_files/msigdb_go_basic_2021-12-15.obo)

**Processed outputs:**

- [Heart transcriptomics FGSEA results](data_processed/heart_transcriptomics_fgsea.csv)
- [Heart transcriptomics FGSEA R object](data_processed/heart_transcriptomics_fgsea.rds)
- [Heart proteomics FGSEA results](data_processed/heart_prot_fgsea.csv)
- [WAT-SC transcriptomics FGSEA R object](data_processed/watsc_transcriptomics_fgsea.rds)
- [WAT-SC proteomics FGSEA R object](data_processed/watsc_proteomics_fgsea.rds)

**Plots:**

The plots labeled `WATSC` reproduce the published findings for subcutaneous white adipose tissue shown in Figures 2a–2c. The liver and heart plots were generated by applying the same analytical methods to those tissues, allowing the team to compare the published adipose-tissue findings with tissue-adapted results.

-  Figure 2a: Differential analysis identifies individual transcripts and proteins that differ between sedentary males and females.
-  Figure 2b: Pathway enrichment compares whether transcriptomic and proteomic results point to the same sex-associated biological processes.
-  Figure 2c: Cross-omic integration highlights pathways that are significantly enriched in both the RNA and protein analyses.

Naming Convention Note: 
-  WATSC: Replication of the published subcutaneous white adipose tissue analysis.
-  Liver: Tissue-adapted analysis using the same methods applied to liver data.
-  Heart: Tissue-adapted analysis using the same methods applied to heart data.
-  

#### WGCNA Analysis Approach

The WGCNA analysis began with processed heart transcriptomics data and associated phenotype information. The data were first reviewed through a sample and data audit to confirm sample matching, group sizes, feature annotations, and potential outliers. This step produced a reviewed expression matrix and an analysis-ready set of phenotype traits.

To reduce computational complexity and focus the network analysis on genes showing meaningful variation across samples, the 5,000 most variable transcriptomic features were selected. These features were used to construct a weighted gene co-expression network, which grouped genes with similar expression patterns into modules.

For each module, a module eigengene was calculated to summarize the overall expression pattern of the genes assigned to that module. Module eigengenes were then evaluated in relation to sex, training status, training duration, and other available phenotype traits. Significant module–trait relationships were identified, and highly connected hub genes were selected to represent each module.

Gene Ontology over-representation analysis was then performed using the genes assigned to each module. This final step identified enriched biological processes and supported interpretation of the molecular functions represented by the heart co-expression modules.

**Data-audit inputs and outputs:**

- [Heart RNA-seq sample metadata](MoTrPAC_rat_heart_data_audit/heart_rnaseq_sample_metadata.csv)
- [Sample matching checks](MoTrPAC_rat_heart_data_audit/sample_matching_checks.csv)
- [Sample counts by group](MoTrPAC_rat_heart_data_audit/sample_counts_by_group.csv)
- [Outlier record](MoTrPAC_rat_heart_data_audit/official_heart_rnaseq_outliers.csv)
- [PCA plot](MoTrPAC_rat_heart_data_audit/heart_rnaseq_pca.png)
- [Completed audit objects](MoTrPAC_rat_heart_data_audit/rat_heart_transcriptomics_audit_objects.rds)

**WGCNA inputs and outputs:**

- [Selected 5,000 variable features](MoTrPAC_rat_heart_GitHub_WGCNA/selected_top5000_variable_features.csv)
- [Final gene-module assignments](MoTrPAC_rat_heart_GitHub_WGCNA/official_WGCNA_modules.csv)
- [Module eigengenes](MoTrPAC_rat_heart_GitHub_WGCNA/official_WGCNA_module_eigengenes_long.csv)
- [Module–trait correlations](MoTrPAC_rat_heart_GitHub_WGCNA/module_trait_correlations.csv)
- [Significant module–trait correlations](MoTrPAC_rat_heart_GitHub_WGCNA/significant_module_trait_correlations.csv)
- [Top hub genes by module](MoTrPAC_rat_heart_GitHub_WGCNA/top10_hub_genes_per_module.csv)
- [Official WGCNA code snapshot](MoTrPAC_rat_heart_GitHub_WGCNA/official_run_WGCNA_snapshot.R)
- [Saved WGCNA result](MoTrPAC_rat_heart_GitHub_WGCNA/official_run_WGCNA_result.rds)

**Gene Ontology module-enrichment outputs:**

- [Module Entrez gene lists](MoTrPAC_rat_heart_differential_expression/module_entrez_gene_lists.csv)
- [Complete module GO results](MoTrPAC_rat_heart_differential_expression/all_module_GO_ORA_results.csv)
- [Significant module GO results](MoTrPAC_rat_heart_differential_expression/significant_module_GO_ORA_results.csv)
- [Module ORA summary](MoTrPAC_rat_heart_differential_expression/module_ORA_summary.csv)
- [Top GO terms by module](MoTrPAC_rat_heart_differential_expression/top10_GO_terms_per_module.csv)

**WGCNA and Related Plot Folders**

The WGCNA-related plots are organized into five folders that reflect the major stages of the heart transcriptomics workflow.

- **`MoTrPAC_rat_heart_data_audit`** — Contains plots used to review the heart RNA-seq data before analysis, including sample-level checks, group balance, principal component analysis, and potential outlier review.

- **`MoTrPAC_rat_heart_differential_expression`** — Contains plots related to transcriptomic differential-expression analysis and baseline male-versus-female comparisons.

- **`MoTrPAC_rat_heart_transcriptomic_WGCNA`** — Contains the primary WGCNA visualizations, including network construction, module identification, module eigengenes, module–trait relationships, and hub-gene summaries.

- **`MoTrPAC_rat_heart_GitHub_WGCNA_module_ORA`** — Contains Gene Ontology over-representation analysis plots used to interpret the biological functions represented by individual WGCNA modules.

- **`MoTrPAC_rat_heart_GitHub_aligned_FGSEA`** — Contains pathway-level comparison plots that align FGSEA results with WGCNA findings or compare related enrichment results across analysis components.

### Dashboard

### Limitations and Notes
-  The published differential-abundance analysis was reproduced using the documented methodology; however, the replication was not exact. The initial plots showed minor differences (~1%) from the published results. These variations may reflect differences in software or package versions, input-file processing, filtering criteria, statistical settings, or other implementation details.
-  For the WGCNA analysis, the 5,000 genes with the greatest variability across samples were selected from approximately 14,000 available genes. This reduced computational complexity and focused the co-expression network on genes showing the strongest variation across study samples.

### External Code and AI Use 
-  Mermaid.ai was used to create the pipeline graphic in the README documentation.
-  AI-assisted tools were used to convert Python to R.
-  AI-assisted tools were used to help write the .html code for the dashboard.
-  AI-assisted tools were used to edit and cross-check the README documentation, improve readability and accessibility, and help verify connections between technical steps.

### Citation
Many, G.M., Sanford, J.A., Sagendorf, T.J. et al. Sexual dimorphism and the multi-omic response to exercise training in rat subcutaneous white adipose tissue. Nat Metab 6, 963–979 (2024). https://doi.org/10.1038/s42255-023-00959-9 

### Software Licenses
This repository includes original analysis code as well as code adapted from the published MoTrPAC workflow. Third-party R and Python packages remain subject to their own software licenses and should be used according to the terms provided by their developers.

### Future Research Questions 

- Do male and female hearts show distinct co-expression modules (WGCNA) linked to phenotypic measures (VO2max, body composition)?
- Do the sex-specific pathways identified overlap with known cardiac disease pathways (e.g., hypertrophic/dilated cardiomyopathy, mitochondrial dysfunction)?

## Contributors and Roles
- Xiuqi Jade Li,
- Viveka Patil,
- Upama
- Sarah Myer, MPH, Data Steward, RENCI (smyer@renci.org)
    Documentation and Writing; Team Facilitator
-  Ebuka
