library(data.table)
library(tibble)
library(fgsea)

heart_pathways <- readRDS(here("data_msigdb", "msigdb_pathways.rds"))

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

# --- Load your proteomics DA ---
prot_dea <- read_csv(here("data_processed", "heart_proteomics_da.csv"))

# Check for semicolon-collapsed multi-gene IDs
sum(grepl(";", prot_dea$entrez_gene), na.rm = TRUE)
# if nonzero:
# prot_dea <- prot_dea %>% separate_rows(entrez_gene, sep = ";")

# --- Rank + run FGSEA ---
prot_ranks <- rank_genes(prot_dea, genes = "entrez_gene")
heart_prot_fgsea <- fgsea2(pathways = heart_pathways, stats = prot_ranks)

write_csv(heart_prot_fgsea, here("data_processed", "heart_prot_fgsea.csv"))

heart_prot_fgsea %>% filter(padj < 0.05) %>% arrange(padj) %>% head(20)