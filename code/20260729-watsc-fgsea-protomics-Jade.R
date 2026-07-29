# Goal ----
# Using Viveka's pipeline for proteomics analysis to reproduce the WAT-SC figures
# DA: 20260729-heart-analysis-proteomics-viveka.R
# FGSEA: 20260729-heart-fgsea-proteomics-viveka.R


# Load libraries ----
library(tidyverse)
library(here)
library(ggrepel)

library(MotrpacRatTraining6moData)
library(limma)

library(data.table)
library(tibble)
library(fgsea)



# Import data ----
pheno <- PHENO
feature_to_gene <- FEATURE_TO_GENE
outliers <- OUTLIERS

prot_norm_raw <- PROT_WATSC_NORM_DATA
prot_meta_heart <- PROT_META %>% filter(tissue == "WAT-SC")


# Differential analysis ----

## 1. Construct sample metadata ----
p_data <- pheno %>%
  filter(tissue == "WAT-SC") %>%
  select(pid:viallabel, sex, timepoint = group) %>%
  mutate(
    sex = factor(stringr::str_to_title(sex), levels = c("Female", "Male")),
    timepoint = ifelse(timepoint == "control", "SED", toupper(timepoint)),
    timepoint = factor(timepoint, levels = c("SED", paste0(2^(0:3), "W"))),
    exp_group = interaction(substr(sex, 1, 1), timepoint, sep = "_")
  ) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = factor(exp_group, levels = unique(exp_group))) %>%
  # add TMT plex as a batch covariate
  inner_join(select(prot_meta_heart, viallabel, tmt_plex), by = "viallabel") %>%
  as.data.frame()

rownames(p_data) <- p_data$viallabel


## 2. Extract the protein expression matrix ----
prot_mat <- prot_norm_raw %>%
  select(feature_ID, where(is.numeric)) %>%
  column_to_rownames("feature_ID") %>%
  as.matrix()

# Match metadata to samples present in the matrix
p_data_prot <- p_data[colnames(prot_mat), ]
p_data_prot <- p_data_prot[!is.na(p_data_prot$viallabel), ]

prot_mat <- prot_mat[, rownames(p_data_prot)]
all.equal(rownames(p_data_prot), colnames(prot_mat))
dim(prot_mat)


## 3. Attach gene annotations (same as RNA script) ----
f_data <- feature_to_gene %>%
  filter(feature_ID %in% rownames(prot_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  group_by(feature_ID) %>%
  filter(!(grepl("^LOC|^NEWGENE", gene_symbol) & !all(grepl("^LOC|^NEWGENE", gene_symbol)))) %>%
  filter(!(all(grepl("^LOC|^NEWGENE", gene_symbol)) & is.na(entrez_gene) & !all(is.na(entrez_gene)))) %>%
  summarise(across(c(gene_symbol, entrez_gene),
                   ~ ifelse(all(is.na(.x)), NA_character_, paste(.x, collapse = ";")))) %>%
  ungroup() %>%
  as.data.frame()

rownames(f_data) <- f_data$feature_ID
f_data <- f_data[rownames(prot_mat), ]
all.equal(rownames(f_data), rownames(prot_mat))


## 4. Check for outliers ----
outliers %>% filter(assay == "PROT")
# no need to worry ;D


## 5. Differential expression: M_SED vs F_SED, with TMT plex as covariate ----
p_data_prot$tmt_plex <- factor(p_data_prot$tmt_plex)

design <- model.matrix(~ 0 + exp_group + tmt_plex, data = p_data_prot)
colnames(design) <- sub("^exp_group", "", colnames(design))

contrast_matrix <- makeContrasts(contrasts = "M_SED - F_SED", levels = design)

linear_fit <- lmFit(prot_mat, design)
contrast_fit <- contrasts.fit(linear_fit, contrast_matrix)
contrast_fit <- eBayes(contrast_fit, robust = TRUE)

prot_MvF_SED <- topTable(
  fit = contrast_fit, coef = 1, number = Inf,
  adjust.method = "BH", sort.by = "none"
) %>%
  rownames_to_column("feature_ID") %>%
  left_join(f_data, by = "feature_ID") %>%
  mutate(contrast = "M_SED - F_SED", se = logFC / t) %>%
  arrange(adj.P.Val, P.Value, feature_ID)

write_csv(prot_MvF_SED, here("data_processed", "watsc_proteomics_da.csv"))

prot_MvF_SED <- read_csv(here("data_processed", "watsc_proteomics_da.csv"))


## Volcano plot ----
fdr_cutoff <- 0.05

prot_MvF_SED_volcano_df <- prot_MvF_SED %>%
  mutate(plot_p = pmax(adj.P.Val, .Machine$double.xmin),
         neg_log10_fdr = -log10(plot_p),
         direction = case_when(adj.P.Val < fdr_cutoff & logFC > 0 ~ "Higher in male",
                               adj.P.Val < fdr_cutoff & logFC < 0 ~ "Higher in female",
                               TRUE ~ "Not significant"),
         direction = factor(direction,
                            levels = c("Higher in female", "Higher in male", "Not significant")),
         feature_label = case_when(!is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
                                   !is.na(entrez_gene) & entrez_gene != "" ~ entrez_gene,
                                   TRUE ~ feature_ID),
         hover_text = paste0("gene_symbol: ", gene_symbol,
                             "<br>feature_ID: ", feature_ID,
                             "<br>entrez_gene: ", entrez_gene,
                             "<br>logFC: ", signif(logFC, 4),
                             "<br>adj.P.Val: ", format(adj.P.Val, scientific = TRUE, digits = 4))
  )

head(prot_MvF_SED_volcano_df)


## tally changes in expression
direction_counts <- prot_MvF_SED_volcano_df %>%
  dplyr::count(direction) %>%
  mutate(label = paste0(direction, ": ", scales::comma(n))) %>%
  pull(label) %>%
  paste(collapse = "   ")

direction_counts


## add gene labels
prot_of_interest <- c("Apeh", "Gm2a", "Sod1", "Npc2", "Acot1",
                      "Lmnb1", "Rps6ka5", "Col14a1", "Slc2a4", "Adipoq")

prot_MvF_SED_volcano_labels <- prot_MvF_SED_volcano_df %>%
  filter(gene_symbol %in% prot_of_interest)


## create basic volcano plot
ggplot(data = prot_MvF_SED_volcano_df,
       mapping = aes(x = logFC, y = neg_log10_fdr, color = direction,
                     text = hover_text)) +
  geom_point(size = 1, alpha = 0.65) +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linewidth = 0.5,
             linetype = "dashed",
             color = "grey35") +
  geom_label_repel(data = prot_MvF_SED_volcano_labels,
                   mapping = aes(label = feature_label), color = "black",
                   size = 3, min.segment.length = 0, box.padding = 0.5,
                   label.padding = 0.2, label.size = 0.2,
                   max.overlaps = Inf, show.legend = FALSE) +
  scale_x_continuous(limits = c(-11, 15), expand = 0) +
  scale_y_continuous(limits = c(0, 30), expand = 0.05) +
  scale_color_manual(values = c("Higher in female" = "#ff6eff",
                                "Higher in male" = "#5555ff",
                                "Not significant" = "grey70"),
                     drop = FALSE) +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_line(color = "grey75", linewidth = 0.2)) +
  labs(title = "Sedentary male vs female WAT-SC protein expression",
       subtitle = direction_counts,
       x = expression(log[2] * "(fold change): Male - Female"),
       y = expression(-log[10] * "(BH-adjusted p-value)"),
       color = NULL)


ggsave(path = here("plots"), filename = "watsc_prot_MvF_SED_volcano.png",
       width = 6.5, height = 5, unit = "in", dpi = 300, bg = "white")



# FGSEA ----
load(here("data_msigdb", "msigdb_pathways.rda"))
msigdb_pathways_selected

## Define functions ----
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

## Load your proteomics DA ---
prot_dea <- read_csv(here("data_processed", "watsc_proteomics_da.csv"))

# Check for semicolon-collapsed multi-gene IDs
sum(grepl(";", prot_dea$entrez_gene), na.rm = TRUE)
# if nonzero:
# prot_dea <- prot_dea %>% separate_rows(entrez_gene, sep = ";")

## Rank + run FGSEA ---
prot_ranks <- rank_genes(prot_dea, genes = "entrez_gene")
prot_fgsea <- fgsea2(pathways = msigdb_pathways_selected, stats = prot_ranks)


## Export FGSEA data ----
saveRDS(prot_fgsea,
        file = here("data_processed", "watsc_proteomics_fgsea.rds"))
