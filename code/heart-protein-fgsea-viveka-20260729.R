library(data.table)
library(tibble)
library(fgsea)

library(tidyverse)
library(here)


msigdb_pathways <- readRDS(here("data_msigdb", "msigdb_pathways.rds"))

rank_genes <- function(x, genes, metric = "-log10(P.Value)*sign(logFC)", contrast_column = "contrast") {
  if (!(genes %in% colnames(x))) stop(sprintf("'%s' is not a column in x.", genes))
  if (!(contrast_column %in% colnames(x))) {
    warning(sprintf("'%s' is not a column in x. Output will be a vector.", contrast_column))
    stats <- .rank_genes(x, genes = genes, metric = metric)
  } else {
    x <- split.data.frame(x = x, f = x[[contrast_column]], drop = TRUE)
    stats <- lapply(x, .rank_genes, genes = genes, metric = metric)
  }
  return(stats)
}

.rank_genes <- function(x, genes, metric) {
  setDT(x)
  x <- subset(x, subset = !is.na(get(genes)))
  if (nrow(x) == 0) stop(sprintf("All entries in '%s' are missing.", genes))
  x[, rank_metric := eval(parse(text = metric))]
  cols <- c(genes, "rank_metric")
  x <- unique(x[, cols, with = FALSE])
  x <- x[, lapply(.SD, mean, na.rm = TRUE), by = genes]
  stats <- sort(deframe(x), decreasing = TRUE)
  names(stats) <- as.character(names(stats))
  return(stats)
}

fgsea2 <- function(pathways, stats, gene_column = "entrez_gene",
                   adjust.method = "BH", adjust.globally = TRUE, seed = NULL, ...) {
  adjust.method <- match.arg(adjust.method, choices = p.adjust.methods)
  if (!is.list(stats)) stop("`stats` must be a named list.")
  setDT(pathways)
  paths <- deframe(pathways[, list(gs_exact_source, get(gene_column))])
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
  by <- "gs_subcat"
  if (!adjust.globally) by <- c(by, "contrast")
  res[, padj := p.adjust(pval, method = adjust.method), by = by]
  return(res)
}


# Load your proteomics DA ---
prot_dea <- read_csv(here("data_processed", "heart_proteomics_da.csv"))

# Check for semicolon-collapsed multi-gene IDs
sum(grepl(";", prot_dea$entrez_gene), na.rm = TRUE)
# if nonzero:
# prot_dea <- prot_dea %>% separate_rows(entrez_gene, sep = ";")


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
# with at least 15 measured memebers remaining, to be consisten with watsc analysis

prot_msigdb <- msigdb_pathways %>%
  mutate(entrez_gene = map(entrez_gene, 
                           ~ intersect(as.character(.x), names(prot_entrez_to_symbol))),
         set_size_post = lengths(entrez_gene),
         ratio = set_size_post / set_size) %>%
  filter(ratio >= 0.85,                  # at least 85% of the original set remains
         set_size_post >= 15)            # at least 15 measured members remaining

table(prot_msigdb$gs_subcat)             # how many gene sets remain?
#GO:BP GO:CC GO:MF 
#129    55    42



# Rank + run FGSEA ---
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


# Export FGSEA results ---
write_csv(prot_fgsea_cleaned, here("data_processed", "heart_proteomics_fgsea.csv"))

saveRDS(prot_fgsea_cleaned, here("data_processed", "heart_proteomics_fgsea.rds"))


prot_fgsea_cleaned %>% filter(padj < 0.05) %>% arrange(padj) %>% head(20)
