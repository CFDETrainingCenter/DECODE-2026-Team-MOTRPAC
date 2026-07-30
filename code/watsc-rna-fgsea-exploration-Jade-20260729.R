# Goal ----
# Examine the code used for the sexual-dimorphism analysis in the WAT-SC
# Article: https://www.nature.com/articles/s42255-023-00959-9
# GitHub: https://github.com/MoTrPAC/MotrpacRatTraining6moWATData/
# Try to replicate FGSEA --- output: transcriptomic part of Fig 2B, 2C


# From examining the GitHub repo, it looks like FGSEA involves the following files
# (1) `data-raw/msigdb_pathways.R` which generates `data/MSIGDB_PAYWAYS.rda`
#     it also uses the function `msigdbr2` from MotrpacRatTraining6moWAT
# (2) `vignettes/MAIN_FGSEA.Rmd` which uses `fgsea2` from MotrpacRatTraining6moWAT
# (3) `vignettes/artciles/plot_FGSEA_male_vs_female.Rmd` which creates Fig 2B and 2C
#     along with proteomics data



# Load libraries -----
library(tidyverse)
library(here)

library(fgsea)                   # fgseaMultilevel()
library(BiocParallel)            # setting BPPARAM in fgseaMultilevel()


# Import data ----
## DA results
watsc_rna_da <- read_csv(here("data_processed", "watsc_transcriptomics_da.csv"))
# for this analysis, we'll need feature_ID, gene_symbol, entrez_gene, logFC, P.Value, and contrast
dim(watsc_rna_da)             # 16404    10


## msigdb pathways 
msigdb_pathways <- readRDS(here("data_msigdb", "msigdb_pathways.rds"))
msigdb_pathways



# Expand 1-to-many transcript-to-gene mappings ----
# during our DA analysis, we collapsed entrez_gene and added ";" in between
watsc_rna_da %>%
  filter(str_detect(entrez_gene, ";"))

# now we want to give each transcript-to-gene mapping its own row
watsc_rna_da_entrez_id <- watsc_rna_da %>%
  separate_rows(entrez_gene, gene_symbol, sep = ";") %>%       # 16,665 × 10
  filter(!is.na(entrez_gene), entrez_gene != "NA") %>%         # 15,537 × 10
  mutate(entrez_gene = as.character(entrez_gene))

head(watsc_rna_da_entrez_id)



# Create an Entrez-ID-to-symbol lookup -----
# for later to translate FGSEA leading-edge IDs.
gene_symbol_by_entrez_id <- watsc_rna_da_entrez_id %>%
  select(entrez_gene, gene_symbol) %>%
  distinct() %>%
  distinct(entrez_gene, .keep_all = TRUE) %>%
  deframe()

head(gene_symbol_by_entrez_id)



# Calculate the transcript-level ranking statistics -----
# positive values represent positive differential expression and 
# negative values represent negative differential expression:
transcript_ranking_statistics <- watsc_rna_da_entrez_id %>%
  mutate(ranking_statistic = -log10(P.Value) * sign(logFC))


# check if any non-finite ranks were produced
any(!is.finite(transcript_ranking_statistics$ranking_statistic))



# Collapse transcript statistics to one rank per Entrez gene ----
# this should reproduce rank_gene(), aka remove identical gene/statistic pairs,
# then average the remaining transcript statistics for each gene and contrast
gene_level_ranking_table <- transcript_ranking_statistics %>%
  mutate(contrast, entrez_gene, ranking_statistic,
         original_row_order = row_number(),
         .keep = "none") %>%
  distinct(contrast, entrez_gene, ranking_statistic, .keep_all = TRUE) %>%
  group_by(contrast, entrez_gene) %>%
  summarise(ranking_statistic = mean(ranking_statistic),
            number_of_distinct_transcript_statistics = n(),
            gene_first_seen = min(original_row_order)) %>%
  ungroup() %>%
  arrange(gene_first_seen)


# fgsea requires one named numeric vector per contrast
ranked_genes_by_contrast <- gene_level_ranking_table %>%
  group_nest(contrast) %>%
  mutate(ranked_genes = map(
    data,
    ~ .x %>%
      # preserve first-seen order when ranks are tied. 
      # the original data.table helper behaved this way, and fgsea tie order can matter.
      arrange(desc(ranking_statistic), gene_first_seen) %>%
      select(entrez_gene, ranking_statistic) %>%
      deframe()
  )) %>%
  select(contrast, ranked_genes) %>%
  deframe()



# Filter pathways to genes tested in this tissue ----
watsc_rna_genes_tested <- watsc_rna_da_entrez_id %>%
  distinct(entrez_gene) %>%
  pull(entrez_gene)


# set filtering criteria
min_pathway_coverage <- 0.85
min_pathway_size_after_filtering <- 15L

tissue_filtered_go_pathways <- msigdb_pathways %>%
  mutate(entrez_gene = map(entrez_gene, ~ intersect(as.character(.x), watsc_rna_genes_tested)),
         pathway_size_after_tissue_filtering = lengths(entrez_gene),
         fraction_of_original_pathway_retained = pathway_size_after_tissue_filtering / set_size) %>%
  filter(fraction_of_original_pathway_retained >= min_pathway_coverage,
         pathway_size_after_tissue_filtering >= min_pathway_size_after_filtering)

nrow(tissue_filtered_go_pathways)       # 3515 retained
nrow(msigdb_pathways_selected)         # 5106 original


# convert the filtered pathway table to the named list required by fgsea
entrez_ids_by_go_pathway <- tissue_filtered_go_pathways %>%
  select(gs_exact_source, entrez_gene) %>%
  deframe()



# Run FGSEA for each contrast ----
fgsea_random_seed <- 0
fgsea_simple_permutations <- 10000L
fgsea_execution_mode <- SerialParam()

fgsea_results_by_contrast <- purrr::imap(
  ranked_genes_by_contrast,
  ~ {
    message("Running FGSEA for contrast: ", .y)
    set.seed(fgsea_random_seed)
    
    fgseaMultilevel(
      pathways = entrez_ids_by_go_pathway,
      stats = .x,
      nPermSimple = fgsea_simple_permutations,
      BPPARAM = fgsea_execution_mode
    ) %>%
      as_tibble() %>%
      mutate(contrast = .y)
  }
)


combined_fgsea_results <- fgsea_results_by_contrast %>%
  bind_rows()



# Add pathway annotation and reproduce the paper's FDR calculation ----
go_pathway_annotations <- tissue_filtered_go_pathways %>%
  select(gs_subcat, gs_exact_source, gs_description) %>%
  distinct()

annotated_fgsea_results <- combined_fgsea_results %>%
  left_join(go_pathway_annotations, by = c("pathway" = "gs_exact_source")) %>%
  group_by(gs_subcat) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup()



# Translate leading-edge Entrez IDs to gene symbols ----
fgsea_results <- annotated_fgsea_results %>%
  mutate(leadingEdge_genes = map(leadingEdge,
                                 ~ unname(gene_symbol_by_entrez_id[as.character(.x)]))) %>%
  select(pathway, gs_subcat, gs_description,
         pval, padj, log2err, ES, NES, size, contrast, 
         leadingEdge, leadingEdge_genes)


# Export file ----
saveRDS(fgsea_results,
        file = here("data_processed", "watsc_transcriptomics_fgsea.rds"))

write_csv(fgsea_results,
          file = here("data_processed", "watsc_transcriptomics_fgsea.csv"))


# Sanity check -----
fgsea_results %>%
  filter(NES > 2) %>%
  arrange(gs_subcat, desc(NES)) %>%
  view()

## the output seems to align with that in Figure 2c,
## since the proteomics dataset decided what was graphed in Fig 2c, 
## we should probably combine the datasets before graphing


# Clear environment ----
rm(list = ls())

