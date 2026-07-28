# DECODE-2026-Team-MOTRPAC

## Problem Statement

Endurance exercise training (ExT) benefits cardiovascular health, but sex-specific molecular mechanisms driving cardiac adaptation remain poorly characterized. MoTrPAC generated the largest public multi-omic resource on this question, profiling 18 solid tissues and blood from male and female rats undergoing endurance training across four time points. Applying this pipeline to adipose tissue revealed profound sexual dimorphism at rest and in the dynamic response to training, with sex-specific signatures preserved or amplified with exercise. No independent replication of this sexual-dimorphism pipeline has been applied to heart tissue, despite matching data being available. 

## Research Questions

- Does sedentary heart tissue show sex-specific differences in gene expression and enriched pathways, similar to adipose tissue?
- How does the cardiac transcriptomic/proteomic response to ExT differ between males and females across the 1/2/4/8-week time course?
- Do male and female hearts show distinct co-expression modules (WGCNA) linked to phenotypic measures (VO2max, body composition)?
- Do the sex-specific pathways identified overlap with known cardiac disease pathways (e.g., hypertrophic/dilated cardiomyopathy, mitochondrial dysfunction)?

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

## Reference Paper

Many, G.M., Sanford, J.A., Sagendorf, T.J. et al. Sexual dimorphism and the multi-omic response to exercise training in rat subcutaneous white adipose tissue. Nat Metab 6, 963–979 (2024). https://doi.org/10.1038/s42255-023-00959-9
