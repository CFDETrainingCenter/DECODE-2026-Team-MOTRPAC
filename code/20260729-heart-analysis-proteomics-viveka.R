# Goal -----
# HEART PROTEOMICS: M_SED vs F_SED differential analysis + volcano plot

library(tidyverse)
library(here)
library(ggrepel)
library(plotly)
library(htmlwidgets)
library(MotrpacRatTraining6moData)
library(limma)
library(htmlwidgets)

# Import data ----
pheno <- PHENO
feature_to_gene <- FEATURE_TO_GENE
outliers <- OUTLIERS

prot_norm_raw <- PROT_HEART_NORM_DATA
prot_meta_heart <- PROT_META %>% filter(tissue == "HEART")


# 1. Construct sample metadata (same pattern as your RNA script) ----
p_data <- pheno %>%
  filter(tissue == "HEART") %>%
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


# 2. Extract the protein expression matrix ----
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


# 3. Attach gene annotations (same as RNA script) ----
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


# 4. Check for outliers ----
outliers %>% filter(assay == "PROT", tissue == "HEART")


# 5. Differential expression: M_SED vs F_SED, with TMT plex as covariate ----
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

write_csv(prot_MvF_SED, here("data_processed", "heart_proteomics_da.csv"))

prot_MvF_SED <- read_csv(here("data_processed", "heart_proteomics_da.csv"))
# Volcano plot ----
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


### tally changes in expression
direction_counts <- prot_MvF_SED_volcano_df %>%
  dplyr::count(direction) %>%
  mutate(label = paste0(direction, ": ", scales::comma(n))) %>%
  pull(label) %>%
  paste(collapse = "   ")

direction_counts


## Create basic volcano plot ----
prot_MvF_SED_volcano <- ggplot(data = prot_MvF_SED_volcano_df,
                               mapping = aes(x = logFC, y = neg_log10_fdr, color = direction,
                                             text = hover_text)) +
  geom_point(size = 1, alpha = 0.65) +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linewidth = 0.5,
             linetype = "dashed",
             color = "grey35") +
  scale_color_manual(values = c("Higher in female" = "#ff6eff",
                                "Higher in male" = "#5555ff",
                                "Not significant" = "grey70"),
                     drop = FALSE) +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_line(color = "grey75", linewidth = 0.2)) +
  labs(title = "Sedentary male vs female HEART protein expression",
       subtitle = direction_counts,
       x = expression(log[2] * "(fold change): Male - Female"),
       y = expression(-log[10] * "(BH-adjusted p-value)"),
       color = NULL)

prot_MvF_SED_volcano


## Add gene labels -----
prot_MvF_SED_volcano_labels <- prot_MvF_SED_volcano_df %>%
  filter(neg_log10_fdr > 5, abs(logFC) > 0.5)   # adjust thresholds based on your data's spread

prot_MvF_SED_volcano +
  geom_label_repel(data = prot_MvF_SED_volcano_labels,
                   mapping = aes(label = feature_label), color = "black",
                   size = 3, min.segment.length = 0, box.padding = 0.5,
                   label.padding = 0.2, label.size = 0.2,
                   max.overlaps = Inf, show.legend = FALSE)

ggsave(path = here("plots"), filename = "heart_prot_MvF_SED_volcano.png",
       width = 6.5, height = 5, unit = "in", dpi = 300, bg = "white")


### Create plotly ----
prot_MvF_SED_volcano_plotly <- plotly::ggplotly(
  prot_MvF_SED_volcano,
  tooltip = "text"
)

prot_MvF_SED_volcano_plotly

htmlwidgets::saveWidget(
  widget = prot_MvF_SED_volcano_plotly,
  file = here("plots", "heart_prot_MvF_SED_volcano_plotly.html"),
  selfcontained = TRUE
)
