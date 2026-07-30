# DECODE-2026-Team-MOTRPAC

## Project Snapshot 

This project adapts and extends the core analytical framework used to identify sex-specific molecular differences in the published MoTrPAC rat endurance-training study of subcutaneous white adipose tissue (Many et al, 2024). 

Using heart tissue from the same parent study, we apply comparable methods to examine baseline male-versus-female differences and sex-specific responses to endurance exercise training. The planned workflow includes matched-sample review, differential analysis, and pathway enrichment using heart RNA-sequencing and global proteomics data. Additional analyses include weighted gene co-expression network analysis (WGCNA), phenotype correlations, RNA–protein integration, and cardiac-disease pathway mapping, as time and data availability permit.

By applying a published sexual-dimorphism and exercise-training analysis framework to a different tissue, the project will evaluate whether cardiac gene expression, protein abundance, and biological pathways show shared or tissue-specific patterns. It will also produce a transparent and reusable workflow that can support technical research and broader scientific interpretation.

Potential users include exercise and cardiovascular researchers, computational biologists, MoTrPAC data users, trainees, and others interested in reproducible multi-omic analysis and sex as a biological variable.

<br>

## Research Questions

Endurance exercise training (ExT) benefits cardiovascular health, but sex-specific molecular mechanisms driving cardiac adaptation remain poorly characterized. MoTrPAC generated a large public multi-omic resource that could be used to address this question. Analysis of adipose tissues revealed substantial sexual dimorphism at rest and in the dynamic response to training, with sex-specific signatures preserved or amplified with exercise. No independent replication of this sexual-dimorphism pipeline has been applied to heart tissue, despite matching data being available. 

**Primary Research Question:**
- Does sedentary heart tissue show sex-specific differences in gene expression and enriched biological pathways, similar to those observed in subcutaneous white adipose tissue?

**Secondary Research Question:**
- How does the cardiac transcriptomic/proteomic response to endurance exercise training differ between males and females across the 1-, 2-, 4-, and 8-week time course?
- Do male and female hearts show distinct co-expression modules (WGCNA) linked to phenotypic measures (VO2max, body composition)?
- Do the sex-specific pathways identified overlap with known cardiac disease pathways (e.g., hypertrophic/dilated cardiomyopathy, mitochondrial dysfunction)?

<br>

## Objectives

1. Reproduce the baseline male-versus-female analytical workflow reported for WAT-SC from Many et al. (2024).
2. Apply the same workflow to matched heart RNA-sequencing and global proteomics data.
3. Identify baseline sex differences in cardiac gene expression, protein abundance, and pathway enrichment.
4. Characterize within-sex molecular responses across the endurance-training time course.
5. Explore heart co-expression modules, phenotype associations, and hub genes.

<br>

## Analysis Strategy

<p align="center">
  <img src="images/Pipeline.svg" alt="DeCoDe Institute MoTrPAC team analysis workflow" width="250">
</p>

The project is organized into three related workflows:

**1. Reproduce the Many et al. WAT-SC analysis.**
Follow the published data-processing, differential-analysis, pathway-enrichment, and cross-omic comparison workflow.

**2. Apply the validated workflow to heart tissue.**
Adapt the WAT-SC workflow to matched heart transcriptomics and proteomics data to examine baseline sex differences and exercise-training responses.

**3. Conduct heart-specific exploratory analyses.**
Use WGCNA, module–trait analysis, hub-gene identification, and Gene Ontology enrichment to characterize cardiac co-expression patterns.

<br>

## Project Status

| Workflow                              | Purpose                                      | Status    |
| ------------------------------------- | -------------------------------------------- | --------- |
| Pathway reference construction        | Recreate the GO data                         | Complete  |
| WAT-SC workflow reproduction          | Reproduce the published Many et al. analysis | Complete  |
| Heart baseline sex-difference analysis | Compare sedentary males and females         | Complete  |
| Heart exercise-training response      | Compare response across training time course | To verify |
| Heart WGCNA             | Identify co-expression modules and phenotype relationships | Complete  |
| Cardiac-disease pathways mapping       | Relate results to cardiac disease biology   | To verify |
| Dashboard                             | Present results interactively                | In progress |

<br>


## Analysis Workflows and Methods

### 0. Shared preparation: pathway-reference construction
To reproduce the pathway-enrichment analysis, we reconstructed the Gene Ontology pathway collection used by the published workflow. The archived December 15, 2021 Gene Ontology structure was processed to create the MSigDB pathway subset used as the reference for FGSEA.

**Implementation and outputs**

* script: `code/msigdb-pathway-exploration-Jade-20260729.R`
* archived GO structure from 2021-12-15: `data_msigdb/msigdb_go_basic_2021-12-15.obo`
* subset of msigdb pathways used in our FGSEA: `data_msigdb/msigdb_pathways.rds`



### 1. Reproduction of the Many et al. WAT-SC analysis

The published WAT-SC workflow was followed from the original transcriptommics and proteomics data through sample and feature preparation, statistical testing, pathway enrichment, and figure generation. The source analysis is available through two related MoTrPAC R packages (`MotrpacRatTraining6moWATData` and `MotrpacRatTraining6moWAT`). We examined data and code in both packages to recreate the workflow. This involves three main analytical steps:

1. Differential analysis identifies individual transcripts and proteins that differ between sedentary males and females
2. Pathway enrichment evaluates whether transcriptomic and proteomic results identify similar sex-associated biological processes
3. Cross-omic integration identifies pathways significantly enriched in both the RNA and protein analyses

**Implementation and outputs**

- input data: pulled from `MotrpacRatTraining6moData`
- analysis script: `watsc_*.R` files in [`code/`](code/),
- stand-alone rendered notebook: [`output_notebooks/MoTrPAC_watsc-MvsF-SED.html`](output_notebooks/MoTrPAC_watsc-MvsF-SED.html)
- processed results: `watsc_*.csv` or `watsc_*.rds` files in [`data_processed/`](data_processed/), 
- figures: `watsc_*.png` or `watsc_*.html` files [`plots/`](plots/)


### 2. Heart baseline sex-difference analysis

The WAT-SC analytical pipeline was adapted to corresponding heart transcriptomics and proteomics datasets for male and female sedentary samples. Separate transcriptomic and proteomic differential analyses identified genes and proteins that differ between males and females. The resulting ranked differential-expression and differential-abundance results were then used as inputs for pathway-enrichment analysis and cross-omic comparison.

**Implementation and outputs**

- input data: pulled from `MotrpacRatTraining6moData`
- analysis script: `heart_*.R` files in [`code/`](code/),
- stand-alone rendered notebook: [`output_notebooks/MoTrPAC_heart-MvsF-SED.html`](output_notebooks/MoTrPAC_heart-MvsF-SED.html)
- processed results: `heart_*.csv` or `heart_*.rds` files in [`data_processed/`](data_processed/), 
- figures: `heart_*.png` or `heart_*.html` files [`plots/`](plots/)


### 3. Heart WGCNA Analysis

The WGCNA analysis begins with processed heart transcriptomics data and associated phenotype information. The data are reviewed to confirm sample matching, group sizes, feature annotations, and potential outliers.

To reduce computational complexity and focus the network analysis on genes showing meaningful variation across samples, the 5,000 most variable transcriptomic features were selected. These features were used to construct a weighted gene co-expression network, which grouped genes with similar expression patterns into modules.

For each module, a module eigengene was calculated to summarize the overall expression pattern of the genes assigned to that module. Module eigengenes were then evaluated in relation to sex, training status, training duration, and other available phenotype traits. Significant module–trait relationships were identified, and highly connected hub genes were selected to represent each module.

Gene Ontology over-representation analysis was then performed using the genes assigned to each module. This final step identified enriched biological processes and supported interpretation of the molecular functions represented by the heart co-expression modules.


**Implementation and outputs**

TBA for Heart WGCNA Analysis


### 4. Dashboard

An interactive dashboard is in development to support exploration of differential-analysis, pathway-enrichment, and cross-omic results. 

**Implementation and outputs**

- dashboard source files are available in `output_dashboard/`


<br>


## Repository Guide

```text
.
├── code                         # R and Python analysis scripts
├── data_raw                     # source molecular and phenotype data --- did we use this?
├── data_msigdb                  # archived GO and MSigDB pathway reference
├── data_processed               # df containing generated result 
├── output_dashboard             # interactive dashboard files
├── output_notebooks             # rendered standalone Quarto reports
├── plots                        # static and interactive figures organized by workflow
└── images                       # README and project documentation images
```


Have a section on naming convention --- will need to adjust after clean up
-  WATSC: Replication of the published subcutaneous white adipose tissue analysis.
-  Liver: Tissue-adapted analysis using the same methods applied to liver data.
-  Heart: Tissue-adapted analysis using the same methods applied to heart data.



## Reproducing this Analysis

### Background Review

- Read the [published paper](https://doi.org/10.1038/s42255-023-00959-9) and relevant supplementary methods
- Review the original published GitHub repositories, R package vignettes, and analysis scripts


### Software Setup

- Install the required R packages. Versions of packages for Analysis Workflows are included at the end of the Quarto documents (see `outputs_notebooks`) 
- Python is used for selected analysis steps. These scripts may be converted to R to support a consistent analysis environment
- Update file paths for the local and cloud environment as needed


<br>



## Limitations and Reflections
-  The published differential-analysis workflow was reproduced using the documented methodology, but the reproduction was not exact. These variations may reflect differences in software or package versions, input-file processing, filtering criteria, statistical settings, or other implementation details.
-  The WGCNA analysis selected the 5,000 genes with the greatest variability from approximately 14,000 available genes. This reduced computational complexity and focused the network on genes showing the strongest variation across study samples.
- Some analyses remain in progress or are planned. Preliminary outputs should not be interpreted as finalized findings.

<br>


## Provenance and Team Contributions

Team members contributed to different components of this project. Individual names may appear in code filenames and output directories to indicate the analysis component each person developed or completed. _Team members listed alphabetically by last name_

* Xiuqi "Jade" Li: transcriptomics analysis pipeline, visualization.
* Viveka Patil: article and methodology review, proteomics analysis pipeline, and dashboard development.
* Upama: heart transcriptomics data audit, WGCNA co-expression, module functional enrichment, and global protein differential analysis.
* Sarah Myer, MPH, Data Steward, RENCI (smyer@renci.org): documentation and writing support; team facilitator.
* Ebuka: [Role to be added.]

<br>



## Generative AI Disclosure
-  Mermaid.ai was used to create the pipeline graphic in the README documentation.
-  AI-assisted tools were used to convert Python to R.
-  AI-assisted tools were used to help write the .html code for the dashboard.
-  AI-assisted tools were used to edit and cross-check the README documentation, improve readability and accessibility, and help verify connections between technical steps.

<br>


## Citation
Many, G.M., Sanford, J.A., Sagendorf, T.J. et al. Sexual dimorphism and the multi-omic response to exercise training in rat subcutaneous white adipose tissue. Nat Metab 6, 963–979 (2024). https://doi.org/10.1038/s42255-023-00959-9 

<br>

## Software Licenses
This repository includes original analysis code as well as code adapted from the published MoTrPAC workflow. Third-party R and Python packages remain subject to their own software licenses and should be used according to the terms provided by their developers.

<br>

