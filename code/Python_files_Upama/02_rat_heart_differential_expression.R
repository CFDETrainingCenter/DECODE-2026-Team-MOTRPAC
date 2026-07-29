# ============================================================
# MoTrPAC Rat Heart Hackathon
# Step 2: Paper-matched transcriptomic differential analysis
# Tissue adaptation of Many et al., Nature Metabolism (2024)
# Author: Upama
# ============================================================

options(timeout = 1e5)
options(repos = c(CRAN = "https://cloud.r-project.org"))

message("R version: ", R.version.string)
message("Working directory: ", getwd())

# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------
cran_packages <- c(
  "remotes", "BiocManager", "dplyr", "tidyr", "tibble",
  "readr", "ggplot2", "ggrepel", "stringr", "statmod"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, Ncpus = 2)
}

bioc_packages <- c("edgeR", "limma")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    Ncpus = 2
  )
}

if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  message("Installing official MotrpacRatTraining6moData package...")
  remotes::install_github(
    "MoTrPAC/MotrpacRatTraining6moData",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE,
    force = TRUE
  )
}

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

message(
  "MotrpacRatTraining6moData version: ",
  as.character(packageVersion("MotrpacRatTraining6moData"))
)

# ------------------------------------------------------------
# 2. Load official MoTrPAC objects
# ------------------------------------------------------------
data("PHENO", package = "MotrpacRatTraining6moData")
data("OUTLIERS", package = "MotrpacRatTraining6moData")
data("TRNSCRPT_META", package = "MotrpacRatTraining6moData")
data("TRNSCRPT_HEART_RAW_COUNTS", package = "MotrpacRatTraining6moData")
data("TRNSCRPT_HEART_NORM_DATA", package = "MotrpacRatTraining6moData")
data("FEATURE_TO_GENE", package = "MotrpacRatTraining6moData")

raw_data <- TRNSCRPT_HEART_RAW_COUNTS
norm_data <- TRNSCRPT_HEART_NORM_DATA

output_dir <- file.path(
  "results",
  "02_rat_heart_differential_expression"
)
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "contrast_tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Identify heart RNA-seq samples and construct metadata
# ------------------------------------------------------------
pheno_ids <- unique(c(
  rownames(PHENO),
  if ("viallabel" %in% colnames(PHENO)) {
    as.character(PHENO$viallabel)
  } else {
    character()
  }
))

sample_ids <- intersect(colnames(raw_data), pheno_ids)

if (length(sample_ids) == 0) {
  stop("No heart RNA-seq sample columns matched PHENO.")
}

if (length(sample_ids) != 50) {
  warning(
    "Expected 50 heart RNA-seq samples, but detected ",
    length(sample_ids), ". The analysis will continue after validation."
  )
}

if ("viallabel" %in% colnames(PHENO)) {
  p_data <- PHENO |>
    filter(as.character(viallabel) %in% sample_ids) |>
    distinct(viallabel, .keep_all = TRUE)

  p_data <- p_data[
    match(sample_ids, as.character(p_data$viallabel)),
    ,
    drop = FALSE
  ]
} else {
  p_data <- PHENO[sample_ids, , drop = FALSE]
  p_data$viallabel <- rownames(p_data)
}

if (anyNA(match(sample_ids, as.character(p_data$viallabel)))) {
  stop("At least one RNA-seq sample could not be matched to PHENO.")
}

required_pheno <- c("viallabel", "sex", "group")
missing_pheno <- setdiff(required_pheno, colnames(p_data))

if (length(missing_pheno) > 0) {
  stop(
    "Required PHENO columns are missing: ",
    paste(missing_pheno, collapse = ", ")
  )
}

p_data <- p_data |>
  mutate(
    viallabel = as.character(viallabel),
    sex = factor(
      str_to_title(as.character(sex)),
      levels = c("Female", "Male")
    ),
    timepoint = ifelse(
      as.character(group) == "control",
      "SED",
      toupper(as.character(group))
    ),
    timepoint = factor(
      timepoint,
      levels = c("SED", "1W", "2W", "4W", "8W")
    ),
    exp_group = interaction(
      substr(as.character(sex), 1, 1),
      timepoint,
      sep = "_"
    )
  ) |>
  arrange(sex, timepoint) |>
  mutate(
    exp_group = factor(
      exp_group,
      levels = c(
        "F_SED", "F_1W", "F_2W", "F_4W", "F_8W",
        "M_SED", "M_1W", "M_2W", "M_4W", "M_8W"
      )
    )
  )

# ------------------------------------------------------------
# 4. Add the four RNA quality covariates used by the paper
# ------------------------------------------------------------
required_qc <- c(
  "viallabel", "RIN", "pct_globin",
  "pct_umi_dup", "median_5_3_bias"
)
missing_qc <- setdiff(required_qc, colnames(TRNSCRPT_META))

if (length(missing_qc) > 0) {
  stop(
    "Required TRNSCRPT_META columns are missing: ",
    paste(missing_qc, collapse = ", ")
  )
}

qc_data <- TRNSCRPT_META |>
  transmute(
    viallabel = as.character(viallabel),
    rin = as.numeric(RIN),
    pct_globin = as.numeric(pct_globin),
    pct_umi_dup = as.numeric(pct_umi_dup),
    median_5_3_bias = as.numeric(median_5_3_bias)
  ) |>
  distinct(viallabel, .keep_all = TRUE)

p_data <- p_data |>
  left_join(qc_data, by = "viallabel")

if (anyNA(p_data$viallabel)) {
  stop("Metadata join produced missing vial labels.")
}

fix_covariate <- function(x, variable_name) {
  if (all(is.na(x))) {
    stop("Covariate ", variable_name, " is entirely missing.")
  }

  x[is.na(x)] <- mean(x, na.rm = TRUE)

  if (is.na(sd(x)) || sd(x) == 0) {
    warning("Covariate ", variable_name, " has zero variance; setting it to zero.")
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
}

covariate_names <- c(
  "rin", "pct_globin", "pct_umi_dup", "median_5_3_bias"
)

for (covariate in covariate_names) {
  p_data[[paste0(covariate, "_unscaled")]] <- p_data[[covariate]]
  p_data[[covariate]] <- fix_covariate(
    p_data[[covariate]],
    covariate
  )
}

# ------------------------------------------------------------
# 5. Remove official outliers, if any
# ------------------------------------------------------------
official_outlier_ids <- intersect(
  as.character(OUTLIERS$viallabel),
  p_data$viallabel
)

p_data <- p_data |>
  mutate(
    official_outlier = viallabel %in% official_outlier_ids
  )

analysis_samples <- p_data |>
  filter(!official_outlier)

if (nrow(analysis_samples) == 0) {
  stop("No samples remain after official outlier removal.")
}

# Re-establish exact sample order and row names.
analysis_samples <- analysis_samples |>
  arrange(sex, timepoint)

rownames(analysis_samples) <- analysis_samples$viallabel

# ------------------------------------------------------------
# 6. Construct count matrix and gene annotations
# ------------------------------------------------------------
if (!"feature_ID" %in% colnames(raw_data)) {
  stop("TRNSCRPT_HEART_RAW_COUNTS does not contain feature_ID.")
}

count_mat <- raw_data |>
  select(feature_ID, all_of(analysis_samples$viallabel)) |>
  column_to_rownames("feature_ID") |>
  as.matrix()

storage.mode(count_mat) <- "numeric"

if (any(count_mat < 0, na.rm = TRUE)) {
  stop("Negative values were detected in the raw count matrix.")
}

if (anyNA(count_mat)) {
  stop("Missing values were detected in the raw count matrix.")
}

# Replicate the annotation-cleaning logic used in the adipose workflow.
annotation <- FEATURE_TO_GENE |>
  filter(feature_ID %in% rownames(count_mat)) |>
  select(feature_ID, gene_symbol, entrez_gene) |>
  distinct() |>
  group_by(feature_ID) |>
  filter(
    !(
      grepl("^LOC|^NEWGENE", gene_symbol) &
        !all(grepl("^LOC|^NEWGENE", gene_symbol))
    )
  ) |>
  filter(
    !(
      all(grepl("^LOC|^NEWGENE", gene_symbol)) &
        is.na(entrez_gene) &
        !all(is.na(entrez_gene))
    )
  ) |>
  summarise(
    gene_symbol = ifelse(
      all(is.na(gene_symbol)),
      NA_character_,
      paste(unique(na.omit(gene_symbol)), collapse = ";")
    ),
    entrez_gene = ifelse(
      all(is.na(entrez_gene)),
      NA_character_,
      paste(unique(na.omit(as.character(entrez_gene))), collapse = ";")
    ),
    .groups = "drop"
  )

feature_data <- tibble(feature_ID = rownames(count_mat)) |>
  left_join(annotation, by = "feature_ID") |>
  as.data.frame()

rownames(feature_data) <- feature_data$feature_ID

# ------------------------------------------------------------
# 7. edgeR filtering and TMM normalization
# ------------------------------------------------------------
dge <- DGEList(
  counts = count_mat,
  samples = as.data.frame(analysis_samples),
  genes = feature_data,
  group = analysis_samples$exp_group,
  remove.zeros = TRUE
)

n_before_filter <- nrow(dge)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

n_after_filter <- nrow(dge)
n_removed_filter <- n_before_filter - n_after_filter

dge <- calcNormFactors(dge, method = "TMM")

message(
  "Low-expression filtering removed ",
  n_removed_filter,
  " features; ",
  n_after_filter,
  " features remain."
)

# ------------------------------------------------------------
# 8. Build the same model used by the adipose paper
# ------------------------------------------------------------
design <- model.matrix(
  ~ 0 + exp_group + rin + pct_globin + pct_umi_dup + median_5_3_bias,
  data = analysis_samples
)

colnames(design) <- sub("^exp_group", "", colnames(design))
rownames(design) <- analysis_samples$viallabel

if (!is.null(nonEstimable(design))) {
  stop(
    "The design matrix is not full rank. Non-estimable terms: ",
    paste(nonEstimable(design), collapse = ", ")
  )
}

# ------------------------------------------------------------
# 9. Voom with sample-specific quality weights
# ------------------------------------------------------------
png(
  filename = file.path(plot_dir, "voom_mean_variance.png"),
  width = 1800,
  height = 1400,
  res = 220
)

voom_object <- voomWithQualityWeights(
  dge,
  design = design,
  var.group = analysis_samples$viallabel,
  plot = TRUE,
  method = "genebygene"
)

dev.off()

sample_weights <- tibble(
  viallabel = analysis_samples$viallabel,
  exp_group = as.character(analysis_samples$exp_group),
  sex = as.character(analysis_samples$sex),
  timepoint = as.character(analysis_samples$timepoint),
  sample_weight = as.numeric(voom_object$targets$sample.weights)
)

sample_weight_plot <- ggplot(
  sample_weights,
  aes(
    x = reorder(viallabel, sample_weight),
    y = sample_weight,
    fill = sex
  )
) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "RNA-seq sample-quality weights",
    subtitle = "Weights estimated by voomWithQualityWeights",
    x = "Heart RNA-seq sample",
    y = "Sample weight"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(plot_dir, "sample_quality_weights.png"),
  sample_weight_plot,
  width = 8,
  height = 10,
  dpi = 300
)

base_fit <- lmFit(voom_object, design)

# ------------------------------------------------------------
# 10. Define the paper's three contrast families
# ------------------------------------------------------------
trained_vs_sed <- c(
  "F_1W - F_SED",
  "F_2W - F_SED",
  "F_4W - F_SED",
  "F_8W - F_SED",
  "M_1W - M_SED",
  "M_2W - M_SED",
  "M_4W - M_SED",
  "M_8W - M_SED"
)

mvf_sed <- c(
  "M_SED - F_SED"
)

mvf_exercise_response <- c(
  "(M_1W - M_SED) - (F_1W - F_SED)",
  "(M_2W - M_SED) - (F_2W - F_SED)",
  "(M_4W - M_SED) - (F_4W - F_SED)",
  "(M_8W - M_SED) - (F_8W - F_SED)"
)

contrast_groups <- list(
  trained_vs_SED = trained_vs_sed,
  MvF_SED = mvf_sed,
  MvF_exercise_response = mvf_exercise_response
)

contrast_definitions <- tribble(
  ~contrast_group, ~contrast, ~plain_language, ~positive_logFC_means,
  "trained_vs_SED", "F_1W - F_SED",
  "Female heart: 1 week versus sedentary",
  "Higher expression after 1 week in females",
  "trained_vs_SED", "F_2W - F_SED",
  "Female heart: 2 weeks versus sedentary",
  "Higher expression after 2 weeks in females",
  "trained_vs_SED", "F_4W - F_SED",
  "Female heart: 4 weeks versus sedentary",
  "Higher expression after 4 weeks in females",
  "trained_vs_SED", "F_8W - F_SED",
  "Female heart: 8 weeks versus sedentary",
  "Higher expression after 8 weeks in females",
  "trained_vs_SED", "M_1W - M_SED",
  "Male heart: 1 week versus sedentary",
  "Higher expression after 1 week in males",
  "trained_vs_SED", "M_2W - M_SED",
  "Male heart: 2 weeks versus sedentary",
  "Higher expression after 2 weeks in males",
  "trained_vs_SED", "M_4W - M_SED",
  "Male heart: 4 weeks versus sedentary",
  "Higher expression after 4 weeks in males",
  "trained_vs_SED", "M_8W - M_SED",
  "Male heart: 8 weeks versus sedentary",
  "Higher expression after 8 weeks in males",
  "MvF_SED", "M_SED - F_SED",
  "Sedentary male heart versus sedentary female heart",
  "Higher expression in sedentary males",
  "MvF_exercise_response",
  "(M_1W - M_SED) - (F_1W - F_SED)",
  "Male versus female training response at 1 week",
  "Larger 1-week training response in males",
  "MvF_exercise_response",
  "(M_2W - M_SED) - (F_2W - F_SED)",
  "Male versus female training response at 2 weeks",
  "Larger 2-week training response in males",
  "MvF_exercise_response",
  "(M_4W - M_SED) - (F_4W - F_SED)",
  "Male versus female training response at 4 weeks",
  "Larger 4-week training response in males",
  "MvF_exercise_response",
  "(M_8W - M_SED) - (F_8W - F_SED)",
  "Male versus female training response at 8 weeks",
  "Larger 8-week training response in males"
)

# ------------------------------------------------------------
# 11. Fit each contrast family and globally adjust within family
# ------------------------------------------------------------
run_contrast_group <- function(group_name, contrasts_vector) {
  contrast_matrix <- makeContrasts(
    contrasts = contrasts_vector,
    levels = design
  )

  fit <- contrasts.fit(base_fit, contrast_matrix)
  fit <- eBayes(fit, trend = FALSE, robust = TRUE)

  result_list <- lapply(seq_along(contrasts_vector), function(i) {
    result <- topTable(
      fit,
      coef = i,
      number = Inf,
      p.value = 1,
      sort.by = "none",
      adjust.method = "BH"
    )

    result$feature <- rownames(result)
    result$contrast <- contrasts_vector[i]
    result$contrast_group <- group_name

    result$se <- (
      sqrt(fit$s2.post) *
        fit$stdev.unscaled[, i]
    )[rownames(result)]

    result$adj.P.Val_per_contrast <- result$adj.P.Val
    result
  })

  result <- bind_rows(result_list)

  # This reproduces limma_full(adjust.globally = TRUE) within each
  # related contrast family used by the published WAT workflow.
  result$adj.P.Val <- p.adjust(
    result$P.Value,
    method = "BH"
  )

  result |>
    left_join(
      contrast_definitions,
      by = c("contrast_group", "contrast")
    ) |>
    mutate(
      direction = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Positive",
        adj.P.Val < 0.05 & logFC < 0 ~ "Negative",
        TRUE ~ "Not significant"
      ),
      significant_FDR05 = adj.P.Val < 0.05,
      significant_FDR05_FC1.5 =
        adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)
    )
}

results_by_group <- lapply(
  names(contrast_groups),
  function(group_name) {
    run_contrast_group(
      group_name,
      contrast_groups[[group_name]]
    )
  }
)
names(results_by_group) <- names(contrast_groups)

all_results <- bind_rows(results_by_group)

# ------------------------------------------------------------
# 12. Summaries
# ------------------------------------------------------------
significant_counts <- all_results |>
  group_by(
    contrast_group,
    contrast,
    plain_language,
    positive_logFC_means
  ) |>
  summarise(
    tested_features = n(),
    FDR05_total = sum(significant_FDR05, na.rm = TRUE),
    FDR05_positive = sum(
      significant_FDR05 & logFC > 0,
      na.rm = TRUE
    ),
    FDR05_negative = sum(
      significant_FDR05 & logFC < 0,
      na.rm = TRUE
    ),
    FDR05_FC1.5_total = sum(
      significant_FDR05_FC1.5,
      na.rm = TRUE
    ),
    minimum_P_value = min(P.Value, na.rm = TRUE),
    minimum_global_FDR = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )

top20_each_contrast <- all_results |>
  group_by(contrast_group, contrast) |>
  arrange(adj.P.Val, P.Value, desc(abs(logFC)), .by_group = TRUE) |>
  slice_head(n = 20) |>
  ungroup()

model_inventory <- tibble(
  item = c(
    "Input raw features",
    "Input RNA-seq samples",
    "Official outliers removed",
    "Features after filterByExpr",
    "Normalized-table feature count for comparison",
    "Design-matrix columns",
    "Contrast families",
    "Total contrasts",
    "Data package version"
  ),
  value = c(
    nrow(raw_data),
    length(sample_ids),
    length(official_outlier_ids),
    n_after_filter,
    nrow(norm_data),
    ncol(design),
    length(contrast_groups),
    sum(lengths(contrast_groups)),
    as.character(packageVersion("MotrpacRatTraining6moData"))
  )
)

# ------------------------------------------------------------
# 13. Volcano plots for the four priority comparisons
# ------------------------------------------------------------
safe_filename <- function(x) {
  x |>
    str_replace_all("[^A-Za-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}

make_volcano <- function(result_table, contrast_name) {
  plot_data <- result_table |>
    filter(contrast == contrast_name) |>
    mutate(
      minus_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
      plot_class = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Significant positive",
        adj.P.Val < 0.05 & logFC < 0 ~ "Significant negative",
        TRUE ~ "Not significant"
      ),
      label = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        feature,
        gene_symbol
      )
    )

  label_data <- plot_data |>
    filter(adj.P.Val < 0.05) |>
    arrange(adj.P.Val, P.Value, desc(abs(logFC))) |>
    slice_head(n = 12)

  definition <- contrast_definitions |>
    filter(contrast == contrast_name) |>
    slice(1)

  ggplot(
    plot_data,
    aes(x = logFC, y = minus_log10_p)
  ) +
    geom_point(
      aes(colour = plot_class),
      alpha = 0.65,
      size = 1.3
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.4
    ) +
    labs(
      title = definition$plain_language,
      subtitle = paste0(
        "Global BH adjustment within contrast family; ",
        "positive logFC: ",
        definition$positive_logFC_means
      ),
      x = "log2 fold change",
      y = "-log10(P value)",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

priority_contrasts <- c(
  "M_SED - F_SED",
  "F_8W - F_SED",
  "M_8W - M_SED",
  "(M_8W - M_SED) - (F_8W - F_SED)"
)

for (contrast_name in priority_contrasts) {
  volcano <- make_volcano(all_results, contrast_name)

  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        "volcano_",
        safe_filename(contrast_name),
        ".png"
      )
    ),
    plot = volcano,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 14. Significant-count overview
# ------------------------------------------------------------
count_plot_data <- significant_counts |>
  select(
    contrast,
    plain_language,
    FDR05_positive,
    FDR05_negative
  ) |>
  pivot_longer(
    cols = c(FDR05_positive, FDR05_negative),
    names_to = "direction",
    values_to = "n"
  ) |>
  mutate(
    signed_n = ifelse(
      direction == "FDR05_negative",
      -n,
      n
    ),
    direction = recode(
      direction,
      FDR05_positive = "Positive logFC",
      FDR05_negative = "Negative logFC"
    ),
    plain_language = factor(
      plain_language,
      levels = rev(contrast_definitions$plain_language)
    )
  )

count_plot <- ggplot(
  count_plot_data,
  aes(
    x = signed_n,
    y = plain_language,
    fill = direction
  )
) +
  geom_col() +
  geom_vline(xintercept = 0) +
  labs(
    title = "Significant heart transcripts across published contrast families",
    subtitle = "FDR < 0.05; negative bars indicate negative log2 fold change",
    x = "Number of significant transcripts",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "significant_transcript_counts.png"),
  count_plot,
  width = 11,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Male-versus-female 8-week response comparison
# ------------------------------------------------------------
female_8w <- all_results |>
  filter(contrast == "F_8W - F_SED") |>
  select(
    feature,
    gene_symbol,
    female_logFC = logFC,
    female_P = P.Value
  )

male_8w <- all_results |>
  filter(contrast == "M_8W - M_SED") |>
  select(
    feature,
    male_logFC = logFC,
    male_P = P.Value
  )

response_8w <- female_8w |>
  inner_join(male_8w, by = "feature") |>
  mutate(
    response_class = case_when(
      female_logFC > 0 & male_logFC > 0 ~ "Positive in both",
      female_logFC < 0 & male_logFC < 0 ~ "Negative in both",
      TRUE ~ "Opposite directions"
    )
  )

response_cor <- cor(
  response_8w$female_logFC,
  response_8w$male_logFC,
  method = "spearman",
  use = "complete.obs"
)

response_plot <- ggplot(
  response_8w,
  aes(
    x = female_logFC,
    y = male_logFC
  )
) +
  geom_point(
    aes(colour = response_class),
    alpha = 0.55,
    size = 1.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Female and male 8-week heart-training responses",
    subtitle = paste0(
      "Transcript-level Spearman correlation = ",
      round(response_cor, 3)
    ),
    x = "Female: 8W − sedentary log2 fold change",
    y = "Male: 8W − sedentary log2 fold change",
    colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "female_vs_male_8w_response.png"),
  response_plot,
  width = 7.5,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 16. Save all outputs
# ------------------------------------------------------------
write_csv(
  model_inventory,
  file.path(output_dir, "model_inventory.csv")
)

write_csv(
  contrast_definitions,
  file.path(output_dir, "contrast_definitions.csv")
)

write_csv(
  as.data.frame(analysis_samples) |>
    rownames_to_column("analysis_row"),
  file.path(output_dir, "analysis_sample_metadata.csv")
)

write_csv(
  as.data.frame(design) |>
    rownames_to_column("viallabel"),
  file.path(output_dir, "design_matrix.csv")
)

write_csv(
  sample_weights,
  file.path(output_dir, "sample_quality_weights.csv")
)

write_csv(
  as.data.frame(dge$genes),
  file.path(output_dir, "filtered_gene_annotations.csv")
)

write_csv(
  significant_counts,
  file.path(output_dir, "significant_gene_counts.csv")
)

write_csv(
  top20_each_contrast,
  file.path(output_dir, "top20_genes_each_contrast.csv")
)

write_csv(
  all_results,
  file.path(output_dir, "all_differential_expression_results.csv")
)

write_csv(
  response_8w,
  file.path(output_dir, "female_vs_male_8w_response.csv")
)

for (group_name in names(results_by_group)) {
  write_csv(
    results_by_group[[group_name]],
    file.path(
      table_dir,
      paste0(group_name, "_all_results.csv")
    )
  )
}

for (contrast_name in unique(all_results$contrast)) {
  contrast_table <- all_results |>
    filter(contrast == contrast_name) |>
    arrange(adj.P.Val, P.Value)

  write_csv(
    contrast_table,
    file.path(
      table_dir,
      paste0(
        safe_filename(contrast_name),
        ".csv"
      )
    )
  )
}

saveRDS(
  list(
    dge = dge,
    voom_object = voom_object,
    design = design,
    base_fit = base_fit,
    results_by_group = results_by_group,
    all_results = all_results,
    significant_counts = significant_counts,
    sample_metadata = analysis_samples,
    contrast_definitions = contrast_definitions
  ),
  file.path(
    output_dir,
    "rat_heart_differential_expression_objects.rds"
  ),
  compress = "xz"
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\n================ MODEL INVENTORY ================\n")
print(model_inventory)

cat("\n================ CONTRAST SUMMARY ===============\n")
print(
  significant_counts |>
    select(
      contrast_group,
      contrast,
      FDR05_total,
      FDR05_positive,
      FDR05_negative,
      FDR05_FC1.5_total,
      minimum_global_FDR
    )
)

cat("\n================ IMPORTANT NOTES ================\n")
cat(
  "1. Positive M_SED - F_SED logFC means higher expression in sedentary males.\n"
)
cat(
  "2. Positive trained-vs-SED logFC means higher expression after training.\n"
)
cat(
  "3. Positive interaction logFC means a larger training response in males.\n"
)
cat(
  "4. adj.P.Val is globally BH-adjusted within each published contrast family.\n"
)

cat("\nSUCCESS: Differential-expression outputs saved to:\n")
cat(normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
