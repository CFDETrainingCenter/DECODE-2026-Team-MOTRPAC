# Goal ----
# Using Viveka's pipeline for proteomics differential analysis to reproduce the WAT-SC protein volano plot
# specifically DA pipeline: heart-proteomics-da-viveka-20260729.R

# one modification: this script uses a different statistical model
# specifically, `~ 0 + exp_group` instead of `~ 0 + exp_group + tmt_plex`, the latter adjusts for TMT 
# decided to try the former, bc that's what's done in `MotrpacRatTraining6moWATData`


# Load libraries ----
library(tidyverse)
library(here)
library(ggrepel)

library(MotrpacRatTraining6moData)
library(limma)


# Import data ----
pheno <- PHENO
feature_to_gene <- FEATURE_TO_GENE
outliers <- OUTLIERS

prot_norm <- PROT_WATSC_NORM_DATA
prot_meta <- PROT_META %>% filter(tissue == "WAT-SC")




# Construct sample metadata ----
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
  as.data.frame()

rownames(p_data) <- p_data$viallabel


# Extract the protein expression matrix from the normalized input ----
prot_mat <- prot_norm %>%
  select(feature_ID, where(is.numeric)) %>%
  column_to_rownames("feature_ID") %>%
  as.matrix()

# subset p_data to just the protein samples
p_data_prot <- p_data %>%
  filter(viallabel %in% colnames(prot_mat))

# rearrange the protein matrix so it matches the order of metadata
prot_mat <- prot_mat[, rownames(p_data_prot)]

all.equal(rownames(p_data_prot), colnames(prot_mat))
dim(prot_mat)           # 9964   60



# Rebuild feature annotations as in data-raw/PROT_EXP.R ----
f_data <- feature_to_gene %>%
  filter(feature_ID %in% rownames(prot_mat)) %>%
  select(feature_ID, gene_symbol, entrez_gene) %>%
  distinct() %>%
  group_by(feature_ID) %>%
  # remove a missing symbol when the same protein has a nonmissing symbol.
  filter(!(is.na(gene_symbol) & any(!is.na(gene_symbol)))) %>%
  ungroup() %>%
  as.data.frame()

# set feature_ID as row names and reorder rows of f_data to match that of prot_mat
rownames(f_data) <- f_data$feature_ID
f_data <- f_data[rownames(prot_mat), ]

all.equal(rownames(f_data), rownames(prot_mat))


# how many gene symbols or entrez IDs are missing?
table(is.na(f_data$gene_symbol))
# FALSE  TRUE
#  9908    56

table(is.na(f_data$entrez_gene))
# FALSE
#  9964

# both matching output in MotrpacRatTraining6moWATData/data-raw/PROT_EXP.R



# Fit the MvF_SED model directly with limma ----
design <- model.matrix(~ 0 + exp_group, data = p_data_prot)
colnames(design) <- sub("^exp_group", "", colnames(design))

contrast_matrix <- makeContrasts("M_SED - F_SED", levels = design)

sample_weights <- arrayWeights(object = prot_mat,
                               design = design,
                               var.group = p_data_prot$viallabel)

linear_fit <- lmFit(object = prot_mat,
                    design = design,
                    weights = sample_weights)

contrast_fit <- contrasts.fit(linear_fit, contrast_matrix)
contrast_fit <- eBayes(contrast_fit,
                       trend = TRUE,
                       robust = TRUE,
                       # PROT_DA was generated before limma changed its robust empirical-Bayes
                       # algorithm. Make the historical behavior explicit for exact replication.
                       legacy = TRUE)

prot_MvF_SED <- topTable(fit = contrast_fit, coef = 1, number = Inf,
                         adjust.method = "BH", sort.by = "none") %>%
  rownames_to_column("feature_ID") %>%
  left_join(f_data, by = "feature_ID") %>%
  mutate(contrast = "M_SED - F_SED",
         se = logFC / t) %>%
  arrange(adj.P.Val, P.Value, feature_ID)


# Export DA output ----
write_csv(prot_MvF_SED, here("data_processed", "watsc_proteomics_da.csv"))



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
# "Higher in female: 2,155   Higher in male: 2,071   Not significant: 5,738"
# matching that reported in Many et al Nature Metab 2024 Fig 2b


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
                   mapping = aes(label = gene_symbol), color = "black",
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
       color = NULL)


ggsave(path = here("plots"), filename = "watsc_protein_MvF_SED_volcano.png",
       width = 6.5, height = 5, unit = "in", dpi = 300, bg = "white")



# Clean environment ----
rm(list = ls())
