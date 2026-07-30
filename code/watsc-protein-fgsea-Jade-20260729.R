# Goal ----
# Using Viveka's pipeline for proteomics FGSEA
# specifically: heart-protein-fgsea-viveka-20260729.R

# modifications: 
# (1) filter for pathways with ≥85% proteome coverage and ≥15 measured genes.
# (2) run `fgseaMultilevel()` with `seed = 0`, `nPermSimple = 10000`, and `nproc = 1`.
# (3) add `leadingEdge_genes`



# Load libraries ----
library(tidyverse)
library(here)

library(data.table)
library(fgsea)


# Read in msigdb_pathway ----
msigdb_pathways <- readRDS(here("data_msigdb", "msigdb_pathways.rds"))
msigdb_pathways



# Define functions (copied from MotrpacRatTraining6moWAT just as in Viveka's code) ----
.rank_genes <- function(x, genes, metric) {
  setDT(x)
  
  # Subset to non-missing genes
  x <- subset(x, subset = !is.na(get(genes)))
  if (nrow(x) == 0) {
    stop(sprintf("All entries in '%s' are missing.", genes))
  }
  
  # Calculate ranking metric for each gene
  x[, rank_metric := eval(parse(text = metric))]
  cols <- c(genes, "rank_metric")
  x <- unique(x[, cols, with = FALSE])
  x <- x[, lapply(.SD, mean, na.rm = TRUE), by = genes] # mean by gene
  
  stats <- sort(deframe(x), decreasing = TRUE) # sorting doesn't actually matter
  names(stats) <- as.character(names(stats))
  return(stats)
}


rank_genes <- function(x, genes, metric = "-log10(P.Value)*sign(logFC)",
                       contrast_column = "contrast") {
  if (!(genes %in% colnames(x))) {
    stop(sprintf("'%s' is not a column in x.", genes))
  }
  if (!(contrast_column %in% colnames(x))) {
    warning(sprintf("'%s' is not a column in x. Output will be a vector.",
                    contrast_column))
    stats <- .rank_genes(x, genes = genes, metric = metric)
  } else {
    x <- split.data.frame(x = x, f = x[[contrast_column]], drop = TRUE)
    stats <- lapply(x, .rank_genes, genes = genes, metric = metric)
  }
  return(stats)
}



fgsea2 <- function(pathways, stats, gene_column = "entrez_gene", adjust.method = "BH",
                   adjust.globally = TRUE, seed = NULL, ...) {
  # Check input
  adjust.method <- match.arg(adjust.method, choices = p.adjust.methods)
  
  if (!is.list(stats)) {
    stop("`stats` must be a named list.")
  }
  
  # List of pathways to test
  setDT(pathways)
  paths <- deframe(pathways[, list(gs_exact_source, get(gene_column))])
  
  # Run FGSEA on each contrast
  res <- lapply(names(stats), function(contrast_i) {
    message(contrast_i)
    set.seed(seed = seed)
    res_i <- fgseaMultilevel(pathways = paths, stats = stats[[contrast_i]], ...)
    res_i[["contrast"]] <- contrast_i
    return(res_i)
  })
  res <- rbindlist(res)
  setDT(res)
  
  pathways <- pathways[, list(gs_subcat, gs_exact_source, gs_description)]
  
  res <- merge(res, pathways, sort = FALSE, by.x = "pathway", by.y = "gs_exact_source")
  
  # p-value adjustment
  by <- "gs_subcat"
  if (!adjust.globally) {
    by <- c(by, "contrast")
  }
  
  res[, padj := p.adjust(pval, method = adjust.method), by = by]
  
  return(res)
}



# Load your proteomics DA ---
prot_dea <- read_csv(here("data_processed", "watsc_proteomics_da.csv"))



# Construct the Entrez-to-symbol map used for leading-edge lables ----
prot_entrez_to_symbol <- prot_dea %>%
  select(entrez_gene, gene_symbol) %>%
  filter(!is.na(entrez_gene)) %>%
  distinct() %>%
  mutate(entrez_gene = as.character(entrez_gene)) %>%
  deframe()

head(prot_entrez_to_symbol)



# Apply the repositories proteome-coverage pathway filter ----
# keep pathways for which at leats 85% of the original MSigDB members were measured
# with at least 15 measured memebers remaining
# see https://github.com/MoTrPAC/MotrpacRatTraining6moWATData/blob/master/vignettes/MAIN_FGSEA.Rmd
# line 104-110

prot_msigdb <- msigdb_pathways %>%
  mutate(entrez_gene = map(entrez_gene, 
                           ~ intersect(as.character(.x), names(prot_entrez_to_symbol))),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85,                  # at least 85% of the original set remains
         set_size_post >= 15)            # at least 15 measured members remaining

table(prot_msigdb$gs_subcat)             # how many gene sets remain?
#GO:BP GO:CC GO:MF 
#230    97    77 

# pretty close to that reported in the Rmd file (linked above)



# Rank proteins and run FGSEA with the reference setting ---
prot_ranks <- rank_genes(prot_dea, genes = "entrez_gene")

prot_fgsea <- fgsea2(pathways = prot_msigdb, stats = prot_ranks,
                     seed = 0, nPermSimple = 10000,
                     adjust.globally = TRUE, nproc = 1)

prot_fgsea_cleaned <- prot_fgsea %>%
  # map Entrez IDs in leading edge subset to gene symbols
  mutate(leadingEdge_genes = map(leadingEdge,
                                 ~ na.omit(prot_entrez_to_symbol[as.character(.x)]))) %>% 
  # rearrange columns
  select(pathway, gs_subcat, gs_description,
         pval, padj, log2err, ES, NES, size,
         contrast, leadingEdge, leadingEdge_genes)



# Export FGSEA data ----
saveRDS(prot_fgsea_cleaned,
        file = here("data_processed", "watsc_proteomics_fgsea.rds"))

# note, csv won't retain the list columns (i.e., leadingEdge and leadingEdge_genes)
write_csv(prot_fgsea_cleaned,
          file = here("data_processed", "watsc_proteomics_fgsea.csv"))



# Clear environment ----
rm(list = ls())
