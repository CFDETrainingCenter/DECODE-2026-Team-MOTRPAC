# ============================================================
# MoTrPAC Rat Heart Hackathon
# Step 6: GitHub-aligned heart global-proteomics analysis
#
# Official component:
#   The exact model/contrast logic from the official GitHub
#   proteomics_timewise_da() function is embedded below.
#
# Installation adaptation:
#   The complete MotrpacRatTraining6mo package is not installed,
#   because its unrelated proteomics_training_da() function imports
#   the CRAN package metap. That dependency is not required for the
#   timewise analysis used in this step.
#
# Official GitHub source:
#   MoTrPAC/MotrpacRatTraining6mo
#   R/proteomics_differential_analysis.R
#
# Official data:
#   MotrpacRatTraining6moData::PROT_HEART_NORM_DATA
#   MotrpacRatTraining6moData::PROT_HEART_DA
#   MotrpacRatTraining6moData::PHENO
#   MotrpacRatTraining6moData::OUTLIERS
#   MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT
#
# Project extension:
#   A joint 10-group limma model generates the same 13 contrasts
#   used for the rat-heart RNA analysis:
#       8 sex-specific trained-versus-sedentary contrasts
#       1 sedentary male-versus-female contrast
#       4 sex-by-training interaction contrasts
#
# Multiple-testing correction for the project extension:
#   - local_BH_FDR: BH within each individual contrast
#   - global_family_BH_FDR: BH pooled across all proteins and
#     contrasts within each related contrast family
#
# The normalized protein values are log2 TMT ratios.
# No missing-value imputation is used for differential analysis.
# Median imputation is used only for PCA visualization.
# ============================================================

options(timeout = 1e5)
options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(0)

message("R version: ", R.version.string)
message("Working directory: ", getwd())

# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

cran_packages <- c(
  "dplyr", "tidyr", "tibble", "purrr", "readr",
  "ggplot2", "stringr", "forcats", "scales",
  "data.table", "digest"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = 2
  )
}

bioc_packages <- c("limma")

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    Ncpus = 2
  )
}

# Install only the official data package.
# The analysis package is intentionally not installed because its
# DESCRIPTION imports metap, which is unrelated to the timewise
# protein analysis and may fail to install in Colab.
if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  remotes::install_github(
    "MoTrPAC/MotrpacRatTraining6moData",
    upgrade = "never",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
  library(data.table)
  library(digest)
  library(limma)
  library(MotrpacRatTraining6moData)
})

message(
  "Official analysis-package source: GitHub snapshot, not installed"
)
message(
  "MotrpacRatTraining6moData version: ",
  as.character(packageVersion("MotrpacRatTraining6moData"))
)
message(
  "limma version: ",
  as.character(packageVersion("limma"))
)

# ------------------------------------------------------------
# 2. Output directories
# ------------------------------------------------------------
output_dir <- file.path(
  "results",
  "06_rat_heart_GitHub_global_proteomics"
)

plot_dir <- file.path(output_dir, "plots")
source_dir <- file.path(output_dir, "official_source_snapshots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Record the exact official GitHub source
# ------------------------------------------------------------
official_source_url <- paste0(
  "https://raw.githubusercontent.com/",
  "MoTrPAC/MotrpacRatTraining6mo/main/",
  "R/proteomics_differential_analysis.R"
)

official_source_file <- file.path(
  source_dir,
  "proteomics_differential_analysis.R"
)

download.file(
  official_source_url,
  official_source_file,
  mode = "wb",
  quiet = FALSE
)

official_source_record <- tibble(
  source_name = "proteomics_differential_analysis.R",
  source_url = official_source_url,
  local_file = official_source_file,
  sha256 = digest(
    official_source_file,
    algo = "sha256",
    file = TRUE
  ),
  downloaded_at = as.character(Sys.time())
)

write_csv(
  official_source_record,
  file.path(output_dir, "official_GitHub_source_record.csv")
)

# ------------------------------------------------------------
# 4. Load official normalized heart-protein data
# ------------------------------------------------------------
# PROT_HEART_NORM_DATA is the object returned by the official
# load_sample_data("HEART", "PROT", exclude_outliers = TRUE)
# before sample-outlier columns are removed.

protein_data <- as.data.frame(
  MotrpacRatTraining6moData::PROT_HEART_NORM_DATA,
  check.names = FALSE
)

outlier_table_all <- as.data.frame(
  MotrpacRatTraining6moData::OUTLIERS,
  check.names = FALSE
)

heart_prot_outliers_for_loading <- outlier_table_all |>
  filter(
    tissue == "HEART",
    assay == "PROT"
  )

outlier_sample_columns <- intersect(
  as.character(
    heart_prot_outliers_for_loading$viallabel
  ),
  colnames(protein_data)
)

if (length(outlier_sample_columns) > 0) {
  protein_data <- protein_data[
    ,
    setdiff(
      colnames(protein_data),
      outlier_sample_columns
    ),
    drop = FALSE
  ]
}

message(
  "Official HEART-PROT outlier columns removed: ",
  length(outlier_sample_columns)
)

required_id_columns <- c(
  "feature_ID",
  "feature",
  "tissue",
  "assay"
)

missing_id_columns <- setdiff(
  required_id_columns,
  colnames(protein_data)
)

if (length(missing_id_columns) > 0) {
  stop(
    "Official protein data are missing expected columns: ",
    paste(missing_id_columns, collapse = ", ")
  )
}

pheno <- as.data.frame(
  MotrpacRatTraining6moData::PHENO,
  check.names = FALSE
)

required_pheno_columns <- c(
  "viallabel",
  "sex",
  "group"
)

missing_pheno_columns <- setdiff(
  required_pheno_columns,
  colnames(pheno)
)

if (length(missing_pheno_columns) > 0) {
  stop(
    "PHENO is missing expected columns: ",
    paste(missing_pheno_columns, collapse = ", ")
  )
}

sample_columns <- intersect(
  colnames(protein_data),
  as.character(pheno$viallabel)
)

if (length(sample_columns) == 0) {
  stop(
    "No normalized protein sample columns matched PHENO$viallabel."
  )
}

sample_metadata <- pheno |>
  filter(viallabel %in% sample_columns) |>
  distinct(viallabel, .keep_all = TRUE)

sample_metadata <- sample_metadata[
  match(
    sample_columns,
    sample_metadata$viallabel
  ),
  ,
  drop = FALSE
]

if (
  any(is.na(sample_metadata$viallabel)) ||
  !identical(
    as.character(sample_metadata$viallabel),
    sample_columns
  )
) {
  stop(
    "Sample metadata could not be aligned to protein columns."
  )
}

sample_metadata <- sample_metadata |>
  mutate(
    sex = tolower(as.character(sex)),
    group = tolower(as.character(group)),
    group10 = case_when(
      sex == "female" & group == "control" ~ "F_SED",
      sex == "female" & group == "1w" ~ "F_1W",
      sex == "female" & group == "2w" ~ "F_2W",
      sex == "female" & group == "4w" ~ "F_4W",
      sex == "female" & group == "8w" ~ "F_8W",
      sex == "male" & group == "control" ~ "M_SED",
      sex == "male" & group == "1w" ~ "M_1W",
      sex == "male" & group == "2w" ~ "M_2W",
      sex == "male" & group == "4w" ~ "M_4W",
      sex == "male" & group == "8w" ~ "M_8W",
      TRUE ~ NA_character_
    ),
    group10 = factor(
      group10,
      levels = c(
        "F_SED", "F_1W", "F_2W", "F_4W", "F_8W",
        "M_SED", "M_1W", "M_2W", "M_4W", "M_8W"
      )
    )
  )

if (any(is.na(sample_metadata$group10))) {
  stop(
    "Some protein samples could not be assigned to the 10-level design."
  )
}

protein_matrix <- as.matrix(
  protein_data[, sample_columns, drop = FALSE]
)

storage.mode(protein_matrix) <- "double"
rownames(protein_matrix) <- as.character(
  protein_data$feature_ID
)

if (anyDuplicated(rownames(protein_matrix)) > 0) {
  stop(
    "The normalized PROT heart object contains duplicated feature_ID values."
  )
}

write_csv(
  sample_metadata,
  file.path(output_dir, "protein_sample_metadata.csv")
)

# ------------------------------------------------------------
# 5. Protein feature annotation
# ------------------------------------------------------------
collapse_unique <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- x[
    !is.na(x) &
      x != "" &
      x != "NA"
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(unique(x), collapse = ";")
}

feature_map <- as.data.frame(
  MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT
) |>
  filter(
    feature_ID %in% rownames(protein_matrix)
  ) |>
  transmute(
    feature_ID = as.character(feature_ID),
    gene_symbol = as.character(gene_symbol),
    entrez_gene = as.character(entrez_gene),
    ensembl_gene = as.character(ensembl_gene),
    rgd_gene = as.character(rgd_gene)
  ) |>
  group_by(feature_ID) |>
  summarise(
    gene_symbol = collapse_unique(gene_symbol),
    entrez_gene = collapse_unique(entrez_gene),
    ensembl_gene = collapse_unique(ensembl_gene),
    rgd_gene = collapse_unique(rgd_gene),
    .groups = "drop"
  )

protein_annotation <- tibble(
  feature_ID = rownames(protein_matrix)
) |>
  left_join(
    feature_map,
    by = "feature_ID"
  )

write_csv(
  protein_annotation,
  file.path(output_dir, "protein_feature_annotation.csv")
)

# ------------------------------------------------------------
# 6. Inventory and missingness
# ------------------------------------------------------------
protein_missingness_by_sample <- tibble(
  viallabel = colnames(protein_matrix),
  n_features = nrow(protein_matrix),
  n_missing = colSums(is.na(protein_matrix)),
  fraction_missing = colMeans(is.na(protein_matrix))
) |>
  left_join(
    sample_metadata |>
      select(viallabel, sex, group, group10),
    by = "viallabel"
  )

protein_missingness_by_feature <- tibble(
  feature_ID = rownames(protein_matrix),
  n_samples = ncol(protein_matrix),
  n_missing = rowSums(is.na(protein_matrix)),
  fraction_missing = rowMeans(is.na(protein_matrix)),
  mean_abundance = rowMeans(
    protein_matrix,
    na.rm = TRUE
  ),
  sd_abundance = apply(
    protein_matrix,
    1,
    sd,
    na.rm = TRUE
  )
) |>
  left_join(
    protein_annotation,
    by = "feature_ID"
  )

write_csv(
  protein_missingness_by_sample,
  file.path(output_dir, "protein_missingness_by_sample.csv")
)

write_csv(
  protein_missingness_by_feature,
  file.path(output_dir, "protein_missingness_by_feature.csv")
)

outlier_table <- as.data.frame(
  MotrpacRatTraining6moData::OUTLIERS
)

heart_prot_outliers <- outlier_table |>
  filter(
    tissue == "HEART",
    assay == "PROT"
  )

write_csv(
  heart_prot_outliers,
  file.path(output_dir, "official_HEART_PROT_outliers.csv")
)

protein_data_inventory <- tibble(
  item = c(
    "Analysis",
    "Tissue",
    "Assay",
    "Normalized data object",
    "Official DA object",
    "Normalized scale",
    "Number of protein features",
    "Number of samples after official outlier exclusion",
    "Female samples",
    "Male samples",
    "Official heart-PROT outliers listed",
    "Features with complete data",
    "Features mapped to at least one gene symbol",
    "Official timewise contrasts",
    "Project integration-ready contrasts",
    "Project global FDR method",
    "Official timewise implementation",
    "MotrpacRatTraining6moData version",
    "Official source SHA-256"
  ),
  value = c(
    "Rat-heart global-proteomics differential analysis",
    "HEART",
    "PROT",
    "PROT_HEART_NORM_DATA",
    "PROT_HEART_DA",
    "Median-MAD normalized, batch-corrected log2 TMT ratio",
    nrow(protein_matrix),
    ncol(protein_matrix),
    sum(sample_metadata$sex == "female"),
    sum(sample_metadata$sex == "male"),
    nrow(heart_prot_outliers),
    sum(rowSums(is.na(protein_matrix)) == 0),
    sum(
      !is.na(protein_annotation$gene_symbol) &
        protein_annotation$gene_symbol != ""
    ),
    8,
    13,
    "BH pooled within each related contrast family",
    "Embedded from official GitHub source; full package not installed",
    as.character(packageVersion("MotrpacRatTraining6moData")),
    official_source_record$sha256
  )
)

write_csv(
  protein_data_inventory,
  file.path(output_dir, "protein_data_inventory.csv")
)

# ------------------------------------------------------------
# 7. PCA for quality assessment
# ------------------------------------------------------------
# PCA only: retain proteins observed in at least 80% of samples
# and replace remaining missing values with feature medians.
pca_keep <- rowMeans(
  !is.na(protein_matrix)
) >= 0.80

pca_matrix <- protein_matrix[
  pca_keep,
  ,
  drop = FALSE
]

feature_medians <- apply(
  pca_matrix,
  1,
  median,
  na.rm = TRUE
)

for (i in seq_len(nrow(pca_matrix))) {
  missing_values <- is.na(pca_matrix[i, ])

  if (any(missing_values)) {
    pca_matrix[i, missing_values] <- feature_medians[i]
  }
}

pca_variable <- apply(
  pca_matrix,
  1,
  sd
)

pca_matrix <- pca_matrix[
  is.finite(pca_variable) &
    pca_variable > 0,
  ,
  drop = FALSE
]

protein_pca <- prcomp(
  t(pca_matrix),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  protein_pca$sdev^2 /
    sum(protein_pca$sdev^2)
)

pca_coordinates <- as_tibble(
  protein_pca$x[, seq_len(min(10, ncol(protein_pca$x))), drop = FALSE],
  rownames = "viallabel"
) |>
  left_join(
    sample_metadata,
    by = "viallabel"
  )

write_csv(
  pca_coordinates,
  file.path(output_dir, "protein_PCA_coordinates.csv")
)

pca_plot <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    shape = sex,
    fill = group
  )
) +
  geom_point(
    size = 4,
    colour = "black",
    alpha = 0.9
  ) +
  labs(
    title = "Rat-heart global-proteomics PCA",
    subtitle = paste0(
      "Official outliers excluded; PC1 ",
      percent(variance_explained[1], accuracy = 0.1),
      ", PC2 ",
      percent(variance_explained[2], accuracy = 0.1)
    ),
    x = paste0(
      "PC1 (",
      percent(variance_explained[1], accuracy = 0.1),
      ")"
    ),
    y = paste0(
      "PC2 (",
      percent(variance_explained[2], accuracy = 0.1),
      ")"
    ),
    shape = "Sex",
    fill = "Training group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "protein_PCA.png"),
  pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# Sample-correlation heatmap.
sample_cor <- cor(
  protein_matrix,
  use = "pairwise.complete.obs",
  method = "spearman"
)

sample_order <- sample_metadata |>
  arrange(
    sex,
    group10,
    viallabel
  ) |>
  pull(viallabel)

sample_cor <- sample_cor[
  sample_order,
  sample_order,
  drop = FALSE
]

sample_cor_long <- as.data.frame(sample_cor) |>
  rownames_to_column("sample_y") |>
  pivot_longer(
    -sample_y,
    names_to = "sample_x",
    values_to = "rho"
  ) |>
  mutate(
    sample_x = factor(
      sample_x,
      levels = sample_order
    ),
    sample_y = factor(
      sample_y,
      levels = rev(sample_order)
    )
  )

sample_correlation_plot <- ggplot(
  sample_cor_long,
  aes(
    x = sample_x,
    y = sample_y,
    fill = rho
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    limits = c(
      min(sample_cor, na.rm = TRUE),
      1
    )
  ) +
  labs(
    title = "Rat-heart protein sample correlations",
    subtitle = "Spearman correlations; official outliers excluded",
    x = NULL,
    y = NULL,
    fill = "Spearman\nrho"
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    plot_dir,
    "protein_sample_correlation_heatmap.png"
  ),
  sample_correlation_plot,
  width = 12,
  height = 11,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Reproduce the official eight timewise contrasts
# ------------------------------------------------------------
# This function preserves the design and contrast logic from the
# official GitHub proteomics_timewise_da() implementation.
# It uses the already loaded official protein matrix and metadata.

limma_result_standard_error <- function(
  fit_object,
  coefficient_index
) {
  as.numeric(
    fit_object$stdev.unscaled[
      ,
      coefficient_index
    ] *
      sqrt(fit_object$s2.post)
  )
}

proteomics_timewise_da_official_logic <- function(
  x,
  metadata,
  outliers,
  tissue = "HEART",
  assay = "PROT"
) {
  timewise_results <- list()

  covariates <- metadata |>
    filter(viallabel %in% colnames(x)) |>
    transmute(
      viallabel = as.character(viallabel),
      sex = if_else(
        sex == "male",
        "M",
        "F"
      ),
      group = as.character(group),
      tr = factor(
        case_when(
          group == "control" ~ "8_1",
          group == "1w" ~ "1_0",
          group == "2w" ~ "2_0",
          group == "4w" ~ "4_0",
          group == "8w" ~ "8_0",
          TRUE ~ NA_character_
        ),
        levels = c(
          "8_1",
          "1_0",
          "2_0",
          "4_0",
          "8_0"
        )
      )
    )

  rownames(covariates) <- covariates$viallabel

  covariates <- covariates[
    colnames(x),
    ,
    drop = FALSE
  ]

  for (current_sex in c("M", "F")) {
    current_metadata <- covariates |>
      filter(sex == current_sex)

    current_counts <- x[
      ,
      rownames(current_metadata),
      drop = FALSE
    ]

    design_official <- model.matrix(
      ~ 0 + tr,
      data = current_metadata
    )

    contrast_official <- makeContrasts(
      tr1_0 - tr8_1,
      tr2_0 - tr8_1,
      tr4_0 - tr8_1,
      tr8_0 - tr8_1,
      levels = design_official
    )

    colnames(contrast_official) <- c(
      "1w",
      "2w",
      "4w",
      "8w"
    )

    fit_official <- lmFit(
      current_counts,
      design_official
    )

    fit_official <- contrasts.fit(
      fit_official,
      contrast_official
    )

    fit_official <- eBayes(
      fit_official
    )

    for (
      coefficient_index in seq_len(
        ncol(fit_official$coefficients)
      )
    ) {
      current_time <- colnames(
        fit_official$coefficients
      )[coefficient_index]

      limma_result <- topTable(
        fit_official,
        coef = coefficient_index,
        number = Inf,
        sort.by = "none"
      )

      case_samples <- rownames(covariates)[
        covariates$sex == current_sex &
          covariates$group == current_time
      ]

      control_samples <- rownames(covariates)[
        covariates$sex == current_sex &
          covariates$group == "control"
      ]

      current_outliers <- outliers |>
        filter(
          tissue == !!tissue,
          assay == !!assay,
          grepl(
            paste0(
              "control|",
              current_time
            ),
            group
          ),
          grepl(
            paste0(
              "^",
              ifelse(
                current_sex == "M",
                "male",
                "female"
              )
            ),
            group
          )
        ) |>
        pull(viallabel)

      current_result <- tibble(
        feature_ID = rownames(limma_result),
        sex = ifelse(
          current_sex == "M",
          "male",
          "female"
        ),
        comparison_group = current_time,
        assay = assay,
        assay_code = "proteomics",
        tissue = tissue,
        tissue_code = "heart",
        removed_samples = if (
          length(current_outliers) > 0
        ) {
          paste(
            current_outliers,
            collapse = ","
          )
        } else {
          NA_character_
        },
        logFC = as.numeric(
          limma_result$logFC
        ),
        logFC_se = limma_result_standard_error(
          fit_official,
          coefficient_index
        ),
        tscore = as.numeric(
          limma_result$t
        ),
        p_value = as.numeric(
          limma_result$P.Value
        ),
        comparison_average_intensity = rowMeans(
          x[
            ,
            case_samples,
            drop = FALSE
          ],
          na.rm = TRUE
        ),
        reference_average_intensity = rowMeans(
          x[
            ,
            control_samples,
            drop = FALSE
          ],
          na.rm = TRUE
        ),
        numNAs = rowSums(
          is.na(
            x[
              ,
              c(
                case_samples,
                control_samples
              ),
              drop = FALSE
            ]
          )
        )
      )

      timewise_results[[
        paste(
          current_sex,
          current_time,
          sep = "_"
        )
      ]] <- current_result
    }
  }

  bind_rows(timewise_results)
}

official_timewise_reproduced <- (
  proteomics_timewise_da_official_logic(
    x = protein_matrix,
    metadata = pheno,
    outliers = outlier_table_all,
    tissue = "HEART",
    assay = "PROT"
  )
) |>
  as_tibble() |>
  left_join(
    protein_annotation,
    by = "feature_ID"
  )

write_csv(
  official_timewise_reproduced,
  file.path(
    output_dir,
    "official_timewise_protein_DA_reproduced.csv"
  )
)

official_reference <- as_tibble(
  MotrpacRatTraining6moData::PROT_HEART_DA
)

write_csv(
  official_reference,
  file.path(
    output_dir,
    "official_PROT_HEART_DA_reference.csv"
  )
)

comparison_keys <- c(
  "feature_ID",
  "sex",
  "comparison_group"
)

official_comparison <- official_reference |>
  select(
    all_of(comparison_keys),
    reference_logFC = logFC,
    reference_p_value = p_value,
    reference_adj_p_value = adj_p_value,
    reference_selection_fdr = selection_fdr
  ) |>
  inner_join(
    official_timewise_reproduced |>
      select(
        all_of(comparison_keys),
        reproduced_logFC = logFC,
        reproduced_p_value = p_value,
        reproduced_tscore = tscore,
        reproduced_numNAs = numNAs
      ),
    by = comparison_keys
  ) |>
  mutate(
    logFC_difference = (
      reproduced_logFC -
        reference_logFC
    ),
    neglog10_reference_p = -log10(
      pmax(
        reference_p_value,
        .Machine$double.xmin
      )
    ),
    neglog10_reproduced_p = -log10(
      pmax(
        reproduced_p_value,
        .Machine$double.xmin
      )
    )
  )

write_csv(
  official_comparison,
  file.path(
    output_dir,
    "official_vs_reproduced_DA_comparison.csv"
  )
)

safe_cor <- function(x, y, method) {
  complete <- complete.cases(x, y)

  if (sum(complete) < 3) {
    return(NA_real_)
  }

  cor(
    x[complete],
    y[complete],
    method = method
  )
}

official_comparison_summary <- official_comparison |>
  group_by(
    sex,
    comparison_group
  ) |>
  summarise(
    n_matched = n(),
    logFC_Pearson = safe_cor(
      reference_logFC,
      reproduced_logFC,
      "pearson"
    ),
    logFC_Spearman = safe_cor(
      reference_logFC,
      reproduced_logFC,
      "spearman"
    ),
    neglog10P_Pearson = safe_cor(
      neglog10_reference_p,
      neglog10_reproduced_p,
      "pearson"
    ),
    maximum_absolute_logFC_difference = max(
      abs(logFC_difference),
      na.rm = TRUE
    ),
    median_absolute_logFC_difference = median(
      abs(logFC_difference),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  official_comparison_summary,
  file.path(
    output_dir,
    "official_vs_reproduced_DA_summary.csv"
  )
)

official_scatter <- ggplot(
  official_comparison,
  aes(
    x = reference_logFC,
    y = reproduced_logFC
  )
) +
  geom_point(
    alpha = 0.25,
    size = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  facet_grid(
    sex ~ comparison_group
  ) +
  labs(
    title = "Official versus reproduced heart-protein log fold changes",
    subtitle = "MoTrPAC proteomics_timewise_da() reproduction",
    x = "PROT_HEART_DA logFC",
    y = "Reproduced logFC"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(
    plot_dir,
    "official_vs_reproduced_logFC.png"
  ),
  official_scatter,
  width = 12,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Integration-ready joint 13-contrast limma model
# ------------------------------------------------------------
design <- model.matrix(
  ~ 0 + group10,
  data = sample_metadata
)

colnames(design) <- levels(
  sample_metadata$group10
)

rownames(design) <- sample_metadata$viallabel

if (
  !identical(
    rownames(design),
    colnames(protein_matrix)
  )
) {
  stop(
    "The 13-contrast design matrix is not aligned with protein samples."
  )
}

contrast_matrix <- makeContrasts(
  F1_vs_FSED = F_1W - F_SED,
  F2_vs_FSED = F_2W - F_SED,
  F4_vs_FSED = F_4W - F_SED,
  F8_vs_FSED = F_8W - F_SED,

  M1_vs_MSED = M_1W - M_SED,
  M2_vs_MSED = M_2W - M_SED,
  M4_vs_MSED = M_4W - M_SED,
  M8_vs_MSED = M_8W - M_SED,

  MSED_vs_FSED = M_SED - F_SED,

  Interaction_1W = (
    M_1W - M_SED
  ) - (
    F_1W - F_SED
  ),

  Interaction_2W = (
    M_2W - M_SED
  ) - (
    F_2W - F_SED
  ),

  Interaction_4W = (
    M_4W - M_SED
  ) - (
    F_4W - F_SED
  ),

  Interaction_8W = (
    M_8W - M_SED
  ) - (
    F_8W - F_SED
  ),

  levels = design
)

contrast_metadata <- tibble(
  contrast = colnames(contrast_matrix),
  contrast_family = c(
    rep("trained_vs_SED", 8),
    "baseline_sex",
    rep("sex_by_training", 4)
  ),
  plain_language = c(
    "Female heart: 1 week versus sedentary",
    "Female heart: 2 weeks versus sedentary",
    "Female heart: 4 weeks versus sedentary",
    "Female heart: 8 weeks versus sedentary",
    "Male heart: 1 week versus sedentary",
    "Male heart: 2 weeks versus sedentary",
    "Male heart: 4 weeks versus sedentary",
    "Male heart: 8 weeks versus sedentary",
    "Sedentary male heart versus sedentary female heart",
    "Male versus female training response at 1 week",
    "Male versus female training response at 2 weeks",
    "Male versus female training response at 4 weeks",
    "Male versus female training response at 8 weeks"
  ),
  positive_direction = c(
    rep("Higher after training", 8),
    "Higher in sedentary males",
    rep("Larger training response in males", 4)
  ),
  negative_direction = c(
    rep("Lower after training", 8),
    "Higher in sedentary females",
    rep("Larger training response in females", 4)
  )
)

write_csv(
  contrast_metadata,
  file.path(output_dir, "protein_contrast_definitions.csv")
)

protein_fit <- lmFit(
  protein_matrix,
  design
)

protein_fit_contrasts <- contrasts.fit(
  protein_fit,
  contrast_matrix
)

protein_fit_ebayes <- eBayes(
  protein_fit_contrasts
)

group_means <- map_dfc(
  levels(sample_metadata$group10),
  function(current_group) {
    current_samples <- sample_metadata$viallabel[
      sample_metadata$group10 == current_group
    ]

    tibble(
      !!paste0("mean_", current_group) := rowMeans(
        protein_matrix[
          ,
          current_samples,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    )
  }
) |>
  mutate(
    feature_ID = rownames(protein_matrix),
    .before = 1
  )

write_csv(
  group_means |>
    left_join(
      protein_annotation,
      by = "feature_ID"
    ),
  file.path(output_dir, "protein_group_means.csv")
)

extract_contrast_result <- function(
  contrast_name,
  contrast_index
) {
  result_table <- topTable(
    protein_fit_ebayes,
    coef = contrast_name,
    number = Inf,
    sort.by = "none"
  )

  standard_error <- (
    protein_fit_contrasts$stdev.unscaled[
      ,
      contrast_index
    ] *
      sqrt(protein_fit_ebayes$s2.post)
  )

  tibble(
    feature_ID = rownames(result_table),
    contrast = contrast_name,
    logFC = as.numeric(result_table$logFC),
    logFC_se = as.numeric(standard_error),
    AveExpr = as.numeric(result_table$AveExpr),
    tscore = as.numeric(result_table$t),
    p_value = as.numeric(result_table$P.Value),
    local_BH_FDR = as.numeric(result_table$adj.P.Val),
    B = as.numeric(result_table$B)
  )
}

project_13_results <- map2_dfr(
  colnames(contrast_matrix),
  seq_len(ncol(contrast_matrix)),
  extract_contrast_result
) |>
  left_join(
    contrast_metadata,
    by = "contrast"
  ) |>
  group_by(contrast_family) |>
  mutate(
    global_family_BH_FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  ungroup() |>
  mutate(
    fold_change = 2^logFC,
    absolute_fold_change = 2^abs(logFC),
    significant_global_FDR05 = (
      !is.na(global_family_BH_FDR) &
        global_family_BH_FDR < 0.05
    ),
    passes_1_5_fold = (
      !is.na(logFC) &
        abs(logFC) >= log2(1.5)
    ),
    direction = case_when(
      logFC > 0 ~ positive_direction,
      logFC < 0 ~ negative_direction,
      TRUE ~ "No change"
    ),
    minus_log10_global_FDR = -log10(
      pmax(
        global_family_BH_FDR,
        .Machine$double.xmin
      )
    )
  ) |>
  left_join(
    protein_annotation,
    by = "feature_ID"
  ) |>
  arrange(
    factor(
      contrast,
      levels = contrast_metadata$contrast
    ),
    global_family_BH_FDR,
    p_value
  )

write_csv(
  project_13_results,
  file.path(
    output_dir,
    "all_13_contrast_protein_results.csv"
  )
)

significant_proteins <- project_13_results |>
  filter(significant_global_FDR05)

write_csv(
  significant_proteins,
  file.path(output_dir, "significant_proteins.csv")
)

top_proteins_per_contrast <- project_13_results |>
  group_by(contrast) |>
  arrange(
    global_family_BH_FDR,
    p_value,
    desc(abs(logFC)),
    .by_group = TRUE
  ) |>
  slice_head(n = 30) |>
  ungroup()

write_csv(
  top_proteins_per_contrast,
  file.path(output_dir, "top30_proteins_per_contrast.csv")
)

protein_significant_counts <- project_13_results |>
  group_by(
    contrast_family,
    contrast,
    plain_language
  ) |>
  summarise(
    n_tested = n(),
    n_significant = sum(
      significant_global_FDR05,
      na.rm = TRUE
    ),
    n_positive = sum(
      significant_global_FDR05 &
        logFC > 0,
      na.rm = TRUE
    ),
    n_negative = sum(
      significant_global_FDR05 &
        logFC < 0,
      na.rm = TRUE
    ),
    n_significant_and_1_5_fold = sum(
      significant_global_FDR05 &
        passes_1_5_fold,
      na.rm = TRUE
    ),
    minimum_global_FDR = min(
      global_family_BH_FDR,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    minimum_global_FDR = ifelse(
      is.infinite(minimum_global_FDR),
      NA_real_,
      minimum_global_FDR
    )
  )

write_csv(
  protein_significant_counts,
  file.path(output_dir, "protein_significant_counts.csv")
)

# Integration-ready gene-level result.
# When several protein features map to one gene, retain the feature
# with the largest absolute moderated t-score within each contrast.
protein_gene_level_13_results <- project_13_results |>
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) |>
  separate_rows(
    gene_symbol,
    sep = ";"
  ) |>
  mutate(
    gene_symbol = str_trim(gene_symbol)
  ) |>
  filter(
    gene_symbol != ""
  ) |>
  group_by(
    contrast,
    gene_symbol
  ) |>
  slice_max(
    order_by = abs(tscore),
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

write_csv(
  protein_gene_level_13_results,
  file.path(
    output_dir,
    "protein_gene_level_13_contrast_results.csv"
  )
)

# ------------------------------------------------------------
# 10. Main result figures
# ------------------------------------------------------------
count_plot_data <- protein_significant_counts |>
  transmute(
    contrast = factor(
      plain_language,
      levels = rev(
        contrast_metadata$plain_language
      )
    ),
    positive = n_positive,
    negative = -n_negative
  ) |>
  pivot_longer(
    cols = c(positive, negative),
    names_to = "direction",
    values_to = "signed_count"
  ) |>
  mutate(
    direction = recode(
      direction,
      positive = "Positive logFC",
      negative = "Negative logFC"
    )
  )

significant_count_plot <- ggplot(
  count_plot_data,
  aes(
    x = signed_count,
    y = contrast,
    fill = direction
  )
) +
  geom_col() +
  geom_vline(xintercept = 0) +
  labs(
    title = "Significant proteins across rat-heart contrasts",
    subtitle = paste(
      "Global BH FDR pooled within each contrast family;",
      "positive and negative logFC shown separately"
    ),
    x = "Number of significant proteins",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    plot_dir,
    "protein_significant_counts.png"
  ),
  significant_count_plot,
  width = 12,
  height = 8,
  dpi = 300
)

make_volcano <- function(
  result_data,
  contrast_name,
  filename
) {
  plot_data <- result_data |>
    filter(contrast == contrast_name) |>
    mutate(
      significance_class = case_when(
        significant_global_FDR05 &
          passes_1_5_fold &
          logFC > 0 ~ "FDR < 0.05 and FC >= 1.5",
        significant_global_FDR05 &
          passes_1_5_fold &
          logFC < 0 ~ "FDR < 0.05 and FC <= 0.67",
        significant_global_FDR05 ~ "FDR < 0.05",
        TRUE ~ "Not significant"
      ),
      display_label = coalesce(
        gene_symbol,
        feature_ID
      )
    )

  label_data <- plot_data |>
    filter(significant_global_FDR05) |>
    arrange(
      global_family_BH_FDR,
      desc(abs(logFC))
    ) |>
    slice_head(n = 12)

  title_text <- contrast_metadata |>
    filter(contrast == contrast_name) |>
    pull(plain_language)

  volcano_plot <- ggplot(
    plot_data,
    aes(
      x = logFC,
      y = minus_log10_global_FDR,
      fill = significance_class
    )
  ) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.15,
      alpha = 0.75,
      size = 1.8
    ) +
    geom_vline(
      xintercept = c(
        -log2(1.5),
        log2(1.5)
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_text(
      data = label_data,
      aes(label = display_label),
      size = 2.5,
      check_overlap = TRUE,
      vjust = -0.6,
      show.legend = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = "Rat-heart global proteomics; family-level global BH FDR",
      x = "Protein log2 fold change",
      y = "-log10(global family BH FDR)",
      fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(
    file.path(plot_dir, filename),
    volcano_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

make_volcano(
  project_13_results,
  "MSED_vs_FSED",
  "baseline_sex_protein_volcano.png"
)

make_volcano(
  project_13_results,
  "F8_vs_FSED",
  "female_8w_protein_volcano.png"
)

make_volcano(
  project_13_results,
  "M8_vs_MSED",
  "male_8w_protein_volcano.png"
)

make_volcano(
  project_13_results,
  "Interaction_8W",
  "interaction_8w_protein_volcano.png"
)

eight_week_comparison <- project_13_results |>
  filter(
    contrast %in% c(
      "F8_vs_FSED",
      "M8_vs_MSED"
    )
  ) |>
  select(
    feature_ID,
    gene_symbol,
    contrast,
    logFC,
    global_family_BH_FDR
  ) |>
  pivot_wider(
    names_from = contrast,
    values_from = c(
      logFC,
      global_family_BH_FDR
    )
  ) |>
  mutate(
    significant_in_either_sex = (
      global_family_BH_FDR_F8_vs_FSED < 0.05 |
        global_family_BH_FDR_M8_vs_MSED < 0.05
    )
  )

write_csv(
  eight_week_comparison,
  file.path(
    output_dir,
    "female_male_8w_protein_logFC_comparison.csv"
  )
)

eight_week_rho <- safe_cor(
  eight_week_comparison$logFC_F8_vs_FSED,
  eight_week_comparison$logFC_M8_vs_MSED,
  "spearman"
)

eight_week_scatter <- ggplot(
  eight_week_comparison,
  aes(
    x = logFC_F8_vs_FSED,
    y = logFC_M8_vs_MSED,
    fill = significant_in_either_sex
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.15,
    alpha = 0.65,
    size = 1.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0
  ) +
  labs(
    title = "Female and male 8-week heart-protein responses",
    subtitle = paste0(
      "Spearman rho = ",
      signif(eight_week_rho, 3)
    ),
    x = "Female 8W versus sedentary log2FC",
    y = "Male 8W versus sedentary log2FC",
    fill = "Significant in\neither sex"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(
    plot_dir,
    "female_male_8w_protein_logFC_comparison.png"
  ),
  eight_week_scatter,
  width = 8,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Save complete objects and session information
# ------------------------------------------------------------
saveRDS(
  list(
    normalized_protein_matrix = protein_matrix,
    sample_metadata = sample_metadata,
    protein_annotation = protein_annotation,
    official_timewise_reproduced = official_timewise_reproduced,
    official_reference = official_reference,
    official_comparison = official_comparison,
    design = design,
    contrast_matrix = contrast_matrix,
    project_13_results = project_13_results,
    protein_gene_level_13_results = protein_gene_level_13_results,
    protein_significant_counts = protein_significant_counts,
    inventory = protein_data_inventory
  ),
  file.path(
    output_dir,
    "rat_heart_GitHub_global_proteomics_complete_objects.rds"
  ),
  compress = "xz"
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

message("")
message(
  "SUCCESS: GitHub-aligned rat-heart global-proteomics analysis completed."
)
message(
  "Results saved to: ",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
)
message(
  "Protein features: ",
  nrow(protein_matrix)
)
message(
  "Samples after official outlier exclusion: ",
  ncol(protein_matrix)
)
message(
  "Official timewise result rows: ",
  nrow(official_timewise_reproduced)
)
message(
  "Project 13-contrast result rows: ",
  nrow(project_13_results)
)
message(
  "Significant protein-contrast associations: ",
  nrow(significant_proteins)
)
