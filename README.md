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

### Processed Data Files (data_processed)
-  **Proteomics pathway-enrichment output:** `heart_prot_fgsea.csv`<br>
   **Description:** This file contains FGSEA pathway-enrichment results from the heart proteomics analysis, including enrichment scores, adjusted significance values, and leading proteins for each pathway.

-  **Proteomics differential-abundance results:** `heart_proteomics_da.csv`<br>
   **Description:** This file contains differential-abundance results from the heart proteomics analysis.

-  **Transcriptomics differential-expression results:** `heart_transcriptomics_da.csv`<br>
   **Description:** This file contains differential gene-expression results from the heart transcriptomics analysis.

-  **Transcriptomics differential-expression results:** `liver_transcriptomics_da.csv`<br>
   **Description:** his file contains differential gene-expression results from the liver transcriptomics analysis.

-  **Transcriptomics differential-expression results:** `liver_transcriptomics_da.csv`<br>
   **Description:** This file contains differential gene-expression results from the white adapose tissue transcriptomics analysis and was generated by reproducing the code and analytical workflow used in the published study.

  ### Plot Outputs (plots)
The plots labeled `WATSC` reproduce the published findings for subcutaneous white adipose tissue shown in Figures 2a–2c. The liver and heart plots were generated by applying the same analytical methods to those tissues, allowing the team to compare the published adipose-tissue findings with tissue-adapted results.

-  Figure 2a: Differential analysis identifies individual transcripts and proteins that differ between sedentary males and females.
-  Figure 2b: Pathway enrichment compares whether transcriptomic and proteomic results point to the same sex-associated biological processes.
-  Figure 2c: Cross-omic integration highlights pathways that are significantly enriched in both the RNA and protein analyses.

Naming Convention Note: 
-  WATSC: Replication of the published subcutaneous white adipose tissue analysis.
-  Liver: Tissue-adapted analysis using the same methods applied to liver data.
-  Heart: Tissue-adapted analysis using the same methods applied to heart data.

## Methods and Provenance

### Team Approach 

Three team members worked on the analysis. You will see references to individual team members in the naming conventions of the code and plots. 

#### Jade

#### Viveka

#### Upama

### Analysis Steps

### Plots

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
