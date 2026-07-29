# DECODE-2026-Team-MOTRPAC

## Project Snapshot 

This project adapts and extends the core analytical framework used to identify sex-specific molecular differences in the published MoTrPAC rat endurance-training study of subcutaneous white adipose tissue. Using heart tissue from the same parent study, we apply comparable methods to examine baseline male-versus-female differences and sex-specific responses to endurance exercise training.

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

- Reproduce the baseline male-versus-female analytical workflow reported for subcutaneous white adipose tissue.**
- Apply the same workflow to matched heart RNA-sequencing and global proteomics data.**
- Identify baseline sex differences in cardiac gene expression, protein abundance, and pathway enrichment.**
- Characterize within-sex molecular responses across the endurance-training time course.**

## Raw Data Files and Data Dictionaries

- **Phenotypic data file:** `motrpac_pass1b-06_pheno_viallabel_data_merged_v4.0.txt`  
  **Description:** This file contains the full set of phenotypic variables collected across study animals and biospecimens. Variables include information related to animal testing, study registration, calculated measures, and specimen processing, as defined in the accompanying data dictionary. The documentation first highlights the key variables needed for basic analyses and then provides additional detail.    
  **Documentation:** [MoTrPAC Phenotype Technical Guide](https://motrpac-data.org/technical-guides/phenotype)

 - **Phenotypic data dictionary:** 'motrpac_pass1b-06_pheno_dictionary_merged_v4.0.txt'
  **Description:** This data dictionary defines the fields included in the phenotypic data tables, along with their scope and allowable values. It is intended to make the phenotypic data easier to interpret and use in analysis.  
**Documentation:** [MoTrPAC Phenotype Technical Guide](https://motrpac-data.org/technical-guides/phenotype)

- **Processed transcriptomic data file:** 'motrpac_pass1b-06_t58-heart_transcript-rna-seq_normalized-log-cpm.txt'
  **Description:** This file contains normalized RNA-sequencing measurements from rat heart tissue in the MoTrPAC PASS1B-06 endurance-training study. Gene expression is reported as log-transformed counts per million, which accounts for differences in sequencing depth across samples and reduces the influence of very large count values. The table provides sample-level gene-expression data that can be used to compare expression patterns across training groups, time points, and sexes.

- **Analysis output:** 'pass1b-06_t58-heart_transcript-rna-seq_training-dea-fdr_20211008.txt'
  **Description:** This file contains statistical results from a time-wise differential gene expression analysis of RNA-sequencing data generated from rat heart tissue in the MoTrPAC PASS1B-06 endurance-training study. 

### Future Research Questions 

- Do male and female hearts show distinct co-expression modules (WGCNA) linked to phenotypic measures (VO2max, body composition)?
- Do the sex-specific pathways identified overlap with known cardiac disease pathways (e.g., hypertrophic/dilated cardiomyopathy, mitochondrial dysfunction)?

## Reference Paper

Many, G.M., Sanford, J.A., Sagendorf, T.J. et al. Sexual dimorphism and the multi-omic response to exercise training in rat subcutaneous white adipose tissue. Nat Metab 6, 963–979 (2024). https://doi.org/10.1038/s42255-023-00959-9
