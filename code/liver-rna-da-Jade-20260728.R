# Goal -----
# Conduct the analysis in `20260728-watsc-exploration-Jade.R` for LIVER

# Load library ----
library(tidyverse)
library(here)

library(ggrepel)
library(plotly)
library(htmlwidgets)
library(tidyHeatmap)

library(MotrpacRatTraining6moData)         # for accessing data
library(edgeR)



# Import data ----
pheno <- PHENO

liver_rna_raw <- TRNSCRPT_LIVER_RAW_COUNTS
rna_meta <- TRNSCRPT_META

outliers <- OUTLIERS
feature_to_gene <- FEATURE_TO_GENE 


# 1. Construct sample metadata ----
p_data <- pheno %>%
  filter(tissue == "LIVER") %>%
  select(pid:viallabel, sex, timepoint = group) %>%
  # clean up columns with group info
  mutate(
    sex = factor(stringr::str_to_title(sex), levels = c("Female", "Male")),
    timepoint = ifelse(timepoint == "control", "SED", toupper(timepoint)),
    timepoint = factor(timepoint, levels = c("SED", paste0(2^(0:3), "W"))),
    exp_group = interaction(substr(sex, 1, 1), timepoint, sep = "_")
  ) %>%
  arrange(sex, timepoint) %>%
  mutate(exp_group = factor(exp_group, levels = unique(exp_group))) %>%
  # add four RNA quality measures
  inner_join(
    select(rna_meta, viallabel, rin = RIN, pct_globin, pct_umi_dup, median_5_3_bias),
    by = "viallabel"
  ) %>%
  as.data.frame()

# set viallabel as rownames
rownames(p_data) <- p_data$viallabel

head(p_data)



# 2. Extract the raw count matrix ----
count_mat <- liver_rna_raw %>%
  select(feature_ID, where(is.numeric)) %>%
  column_to_rownames("feature_ID") %>%
  as.matrix()

# rearrange column order so it matches p_data 
count_mat <- count_mat[, rownames(p_data)]

all.equal(row.names(p_data), colnames(count_mat))



# 3. Attach gene annotations -----
f_data <- feature_to_gene %>%
  # limit to features in the count_mat
  filter(feature_ID %in% rownames(count_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  group_by(feature_ID) %>%
  # for each transcript, remove genes that start with "LOC" or "NEWGENE" unless there are no other genes
  filter(!(grepl("^LOC|^NEWGENE", gene_symbol) & !all(grepl("^LOC|^NEWGENE", gene_symbol)))) %>%
  # if all genes start with "LOC" or "NEWGENE", only keep those with Entrez IDs unless none of them have Entrez IDs
  filter(!(all(grepl("^LOC|^NEWGENE", gene_symbol)) & is.na(entrez_gene) & !all(is.na(entrez_gene)))) %>%
  # collapse duplicates with semi-colons
  summarise(
    across(
      c(gene_symbol, entrez_gene),
      ~ ifelse(all(is.na(.x)), NA_character_, paste(.x, collapse = ";"))
    )
  ) %>%
  ungroup() %>%
  as.data.frame()


# set feature_ID as row names and reorder f_data
rownames(f_data) <- f_data$feature_ID
f_data <- f_data[rownames(count_mat), ]

all.equal(rownames(f_data), rownames(count_mat))



# 4. Filter lowly expressed features -----
dge_raw <- DGEList(
  counts = count_mat,
  samples = p_data,
  genes = f_data,
  group = p_data$exp_group,
  remove.zeros = TRUE                 
)                          ## Removing 11746 rows with all zero counts                        

summary(dge_raw)
head(dge_raw)
dim(dge_raw$counts)           # 21137    50

keep <- filterByExpr(dge_raw)           # determine which genes have sufficiently large counts
table(keep)                             # 14,215 TRUE and 6,922 FALSE

dge_raw <- dge_raw[keep, , keep.lib.sizes = FALSE]
summary(dge_raw)
dim(dge_raw$counts)                     # 14215    50


## caculate TMM factors here to create the log-CPM MDS plot used for outlier inspection.
dge_for_mds <- calcNormFactors(dge_raw, method = "TMM")

dge_for_mds$samples$bid == str_sub(dge_for_mds$samples$viallabel, start = 1, end = 5) 
# bid is the first 5 digits of viallabel

plotMDS(cpm(dge_for_mds, log = TRUE), top = 1000,
        label = dge_for_mds$samples$bid, 
        dim.plot = c(1, 2))
## 90218 is an extreme outlier

plotMDS(cpm(dge_for_mds, log = TRUE), top = 1000,
        label = dge_for_mds$samples$bid, 
        dim.plot = c(1, 3))
## 90416 is an extreme outlier --- discard



# 5. Remove the published RNA outliers -----
glimpse(outliers)
unique(outliers$tissue)
unique(outliers$assay)

outliers %>%
  filter(assay == "TRNSCRPT", tissue == "LIVER")
## only identifies 1...

outlier_vials <- paste0(c("90218", "90416"), "016803")
outlier_vials

dge_raw <- dge_raw[, !colnames(dge_raw) %in% outlier_vials]



## 6. Mean-impute, center, and scale RNA QC covariates -----
##   rin                 RNA integrity number
##   pct_globin          percentage of reads assigned to globin
##   pct_umi_dup         percentage of duplicated UMIs
##   median_5_3_bias     median 5'-to-3' coverage bias

qc_covariates <- c("rin", "pct_globin", "pct_umi_dup", "median_5_3_bias")

dge_raw$samples <- dge_raw$samples %>%
  #  select(-group, -lib.size, -norm.factors) %>%       # needed for voomWithQualityWeights() to work   
  mutate(across(all_of(qc_covariates),
                function(x) {
                  x[is.na(x)] <- mean(x, na.rm = TRUE)
                  as.numeric(scale(x))
                }))


# 7. Calculate TMM normalization factors ----
# recalculated after outlier removal
dge_raw <- calcNormFactors(dge_raw, method = "TMM")

summary(dge_raw)
dge_raw$samples$exp_group



# 8. Differential expression analysis -------

## 8.1 Specify the feature-wise linear model ----
model_formula <- paste(
  "~ 0 + exp_group +",
  paste(qc_covariates, collapse = " + ")
)

model_formula


## 8.2 Convert the model formula into a numeric design matrix ----
design <- model.matrix(
  object = as.formula(model_formula),
  data = dge_raw$samples
)

head(design)

colnames(design) <- sub("^exp_group", "", colnames(design))     # remove prefix
head(design)


## 8.3 Encode the M vs F comparison as the contrast ----
## +1 to the adjusted M_SED group mean,
## -1 to the adjusted F_SED group mean,
##  0 to every training group and QC coefficient.

## consequently, positive logFC values mean higher expression in sedentary males
## and negative logFC values mean higher expression in sedentary females.

contrast_matrix <- makeContrasts(
  contrasts = "M_SED - F_SED",
  levels = design
)

contrast_matrix



## 8.4: Transform counts and estimate precision weights ----
## voomWithQualityWeights() first converts counts to TMM-normalized log2 counts
## per million using the effective library sizes stored in `dge_raw`.
## Voom talks into account RNA-seq variability ~ abundance (low low-count 
## observations are generally noisier than high-count observations)

voom_fit <- voomWithQualityWeights(
  counts = dge_raw,
  design = design,
  var.group = factor(dge_raw$samples$viallabel),
  plot = FALSE
)

summary(voom_fit)

## `voom_fit$E` contains normalized log2 expression values.
head(voom_fit$E)

## `voom_fit$weights` contains the combined observation/sample precision weights.
head(voom_fit$weights)


## 8.5 Fit one weighted linear model per retained feature ----
linear_fit <- lmFit(voom_fit, design)

summary(linear_fit)
linear_fit$method


## 8.6: Replace individual model coefficients with M_SED - F_SED ----
## calculates the adjusted difference between the two SED group coefficients 
## and propagates its uncertainty.
contrast_fit <- contrasts.fit(linear_fit, contrast_matrix)

summary(contrast_fit)


## 8.7: Stabilize feature-wise variance estimates ----
## shares information across all features and
## shrinks individual variance estimates toward a common prior distribution,
## producing the moderated t statistic used for inference.
contrast_fit <- eBayes(contrast_fit, robust = TRUE)
summary(contrast_fit)


## 8.8: Extract all tested features and control the FDR ----
rna_MvF_SED <- topTable(
  fit = contrast_fit,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
) %>%
  mutate(
    # record the comparison explicitly in every result row
    contrast = "M_SED - F_SED",
    # calculate se reported by the original wrapper
    se = logFC / t
  ) %>%
  # sort the exported table by p values, while retaining feature ID as a deterministic tie-breaker. 
  # sorting does not affect BH adjustment.
  arrange(adj.P.Val, P.Value, feature_ID) %>%
  # omit B not in the repository's saved differential-analysis objects.
  select(-any_of("B"))


## export results
write_csv(rna_MvF_SED, here("data_processed", "liver_transcriptomics_da.csv"))




# 9. Volcano-plot! ----

## Wrangle data for volcano plot ----
## the manuscript uses BH-adjusted P value ≥ 0.05 to determine significance
fdr_cutoff <- 0.05

rna_MvF_SED_volcano_df <- rna_MvF_SED %>%
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

head(rna_MvF_SED_volcano_df)


### tally changes in expression
direction_counts <- rna_MvF_SED_volcano_df %>%
  dplyr::count(direction) %>%
  mutate(label = paste0(direction, ": ", scales::comma(n))) %>%
  pull(label) %>%
  paste(collapse = "   ")

direction_counts
## close enough to those reported in the manuscript



## Create basic volcano plot ----
rna_MvF_SED_volcano <- ggplot(data = rna_MvF_SED_volcano_df,
                              mapping = aes(x = logFC, y = neg_log10_fdr, color = direction,
                                            text = hover_text)) +
  geom_point(size = 1, alpha = 0.65) +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linewidth = 0.5,
             linetype = "dashed",
             color = "grey35") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(from = 0, to = 30, by = 10),
                     expand = 0.01) +
  scale_color_manual(values = c("Higher in female" = "#ff6eff",
                                "Higher in male" = "#5555ff",
                                "Not significant" = "grey70"),
                     drop = FALSE) +
  theme_classic(base_size = 10) +
  theme(#legend.position = "top", 
    #legend.background = element_rect(color = "grey35", linewidth = 0.2),
    panel.grid = element_line(color = "grey75", linewidth = 0.2)
  ) +
  labs(title = "Sedentary male vs female LIVER transcript expression",
       subtitle = direction_counts,
       x = expression(log[2] * "(fold change): Male - Female"),
       y = expression(-log[10] * "(BH-adjusted p-value)"),
       color = NULL)

rna_MvF_SED_volcano


## Add gene labels -----
## label features with the five smallest adjusted p values in each direction
rna_MvF_SED_volcano_labels <- rna_MvF_SED_volcano_df %>%
  filter(direction != "Not significant") %>%
  group_by(direction) %>%
  slice_min(order_by = adj.P.Val, n = 5, with_ties = FALSE) %>%
  ungroup()

rna_MvF_SED_volcano +
  geom_label_repel(data = rna_MvF_SED_volcano_labels,
                   mapping = aes(label = feature_label), color = "black",
                   size = 3, min.segment.length = 0, box.padding = 0.5,
                   label.padding = 0.2, label.size = 0.2,
                   max.overlaps = Inf, show.legend = FALSE)

ggsave(path = here("plots"), filename = "liver_rna_MvF_SED_volcano.png", 
       width = 6.5, height = 5, unit = "in", dpi = 300, bg = "white")



### Create plotly ----
rna_MvF_SED_volcano_plotly <- plotly::ggplotly(
  rna_MvF_SED_volcano,
  tooltip = "text"
)

rna_MvF_SED_volcano_plotly

htmlwidgets::saveWidget(
  widget = rna_MvF_SED_volcano_plotly,
  file = here("plots", "liver_rna_MvF_SED_volcano_plotly.html"),
  selfcontained = TRUE
)



# PCA plot -----
# Codex suggested using voom_fit$E, which are filtered, outlier-removed, TMM-normalized log2-expression matrix
pca_input <- t(voom_fit$E)

# genes are already log2 transformed, so center but no need to variance scale them
pca_fit <- prcomp(pca_input,
                  center = TRUE,
                  scale. = FALSE)

pca_summary <- summary(pca_fit)
pca_summary

variance_explained <- pca_summary$importance["Proportion of Variance", ] * 100      
# % variance accounted for by each principle
# it's the same as calculating pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100

pca_df <- as.data.frame(pca_fit$x[, 1:4]) %>%
  rownames_to_column("viallabel") %>%
  left_join(dge_raw$samples %>%
              select(viallabel, pid, bid, sex, timepoint, exp_group),
            by = "viallabel")

ggplot(data = pca_df,
      mapping = aes(x = PC1, y = PC2, color = timepoint, shape = sex)) +
  geom_point(size = 3, alpha = 0.52) +
#  stat_ellipse(mapping = aes(group = exp_group),
#               linewidth = 0.5,
#               linetype = "dashed",
#               show.legend = FALSE) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 10) +
  labs(title = "Liver transcriptomic PCA",
       subtitle = "TMM-normalized voom log2 expression",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Training", shape = "Sex")

ggsave(path = here("plots"), filename = "liver_rna_MvF_SED_PCA.png", 
       width = 6, height = 5, unit = "in", dpi = 300, bg = "white")



# Heatmap -----
# select top differentially expressed genes
n_per_direction <- 50
logfc_cutff <- 1

heatmap_top_genes <- rna_MvF_SED %>%
  filter(adj.P.Val < 0.05,
         abs(logFC) > logfc_cutff) %>%
  dplyr::mutate(direction = if_else(logFC > 0, "Higher in male", "Higher in female")) %>%
  group_by(direction) %>%
  slice_min(order_by = adj.P.Val,
            n = n_per_direction,
            with_ties = FALSE) %>%
  ungroup()


# create df for heatmap
heatmap_df <- as.data.frame(voom_fit$E) %>% 
  rownames_to_column(var = "feature_ID") %>%
  # filter for top differentially expressed genes
  filter(feature_ID %in% heatmap_top_genes$feature_ID) %>%
  pivot_longer(col = -feature_ID, names_to = "viallabel", values_to = "expression") %>%
  # join with sample info
  left_join(voom_fit$targets %>%
              select(viallabel, sex, timepoint, exp_group),
            by = "viallabel") %>%
  # filter for just the SED animals
  filter(timepoint == "SED") %>%
  # join with feature info
  left_join(voom_fit$genes, by = "feature_ID") %>%
  # calculate z-score
  group_by(feature_ID) %>%
  mutate(expression_z = (expression - mean(expression)) / sd(expression), .after = expression) %>%
  ungroup() %>%
  mutate(bid = str_sub(viallabel, start = 1, end = 5)) %>%
  distinct(bid, gene_symbol, .keep_all = TRUE)

head(heatmap_df)


heatmap_de <- tidyHeatmap::heatmap(.data = heatmap_df,
                                   .row = gene_symbol,
                                   .column = bid,
                                   .value = expression_z,
                                   clustering_method_columns = "complete",
                                   clustering_method_rows = "complete",
                                   show_row_names = TRUE,
                                   palette_value = rev(RColorBrewer::brewer.pal(3, "RdBu")),
                                   column_title = "Liver Transcriptomics | Heatmap\n(top 50 DE genes in either direction)") %>%
  annotation_tile (sex, palette = c("#ff6eff", "#5555ff"))

heatmap_de

save_pdf(heatmap_de, here("plots", "liver_rna_MvF_SED_Heatmap.pdf"),
         width = 7, height = 7)



# Clear environment ----
rm(list = ls())




