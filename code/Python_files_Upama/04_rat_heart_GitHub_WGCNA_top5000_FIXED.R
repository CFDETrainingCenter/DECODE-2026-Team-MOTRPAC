# ============================================================
# MoTrPAC Rat Heart Hackathon
# Step 4: GitHub-matched transcriptomic WGCNA with top-5,000 variance prefilter
#
# Network construction is performed by sourcing and calling:
# MoTrPAC/MotrpacRatTraining6moWAT/R/run_WGCNA.R
#
# The official run_WGCNA() algorithm and parameters are unchanged.
# Project-specific input adaptations:
# 1. WAT transcriptomics -> HEART transcriptomics
# 2. Select the 5,000 highest-variance normalized heart transcripts before the call
#
# Author: Upama
# Patch: resume support and data.table-safe metadata alignment
# ============================================================

options(timeout = 1e5)
options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(42)

message("R version: ", R.version.string)
message("Working directory: ", getwd())

# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "remotes", "WGCNA", "dynamicTreeCut", "fastcluster",
  "data.table", "tibble", "dplyr", "tidyr", "readr",
  "ggplot2", "stringr", "forcats", "scales", "digest"
)

bioc_packages <- c(
  "Biobase", "impute", "preprocessCore",
  "GO.db", "AnnotationDbi", "limma"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  message(
    "Installing Bioconductor packages: ",
    paste(missing_bioc, collapse = ", ")
  )

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    Ncpus = 2
  )
}

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  message(
    "Installing CRAN packages: ",
    paste(missing_cran, collapse = ", ")
  )

  install.packages(
    missing_cran,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = 2
  )
}

if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  message("Installing MotrpacRatTraining6moData...")

  remotes::install_github(
    "MoTrPAC/MotrpacRatTraining6moData",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE,
    force = TRUE
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(MotrpacRatTraining6moData)
  library(WGCNA)
  library(dynamicTreeCut)
  library(data.table)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
  library(limma)
  library(digest)
})

allowWGCNAThreads(nThreads = 2)

# ------------------------------------------------------------
# 2. Output folders
# ------------------------------------------------------------
output_dir <- file.path(
  "results",
  "04_rat_heart_GitHub_WGCNA_top5000"
)

plot_dir <- file.path(output_dir, "plots")
module_dir <- file.path(output_dir, "module_tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Download and freeze the official run_WGCNA implementation
# ------------------------------------------------------------
official_source_url <- paste0(
  "https://raw.githubusercontent.com/",
  "MoTrPAC/MotrpacRatTraining6moWAT/",
  "master/R/run_WGCNA.R"
)

official_source_file <- file.path(
  output_dir,
  "official_run_WGCNA_snapshot.R"
)

if (!file.exists(official_source_file)) {
  message("Downloading the official run_WGCNA.R source...")

  download.file(
    official_source_url,
    destfile = official_source_file,
    mode = "wb",
    quiet = FALSE
  )
} else {
  message("Using the previously saved official run_WGCNA.R snapshot.")
}

official_source_sha256 <- digest(
  official_source_file,
  algo = "sha256",
  file = TRUE
)

writeLines(
  c(
    paste0("Source URL: ", official_source_url),
    paste0("SHA-256: ", official_source_sha256),
    paste0("Downloaded/used on: ", Sys.time())
  ),
  file.path(output_dir, "official_run_WGCNA_source_record.txt")
)

# Source the official function directly.
source(official_source_file, local = .GlobalEnv)

if (!exists("run_WGCNA", mode = "function")) {
  stop("The official run_WGCNA function was not loaded.")
}

# ------------------------------------------------------------
# 4. Load heart transcriptomics and metadata
# ------------------------------------------------------------
data("PHENO", package = "MotrpacRatTraining6moData")
data("OUTLIERS", package = "MotrpacRatTraining6moData")
data("TRNSCRPT_META", package = "MotrpacRatTraining6moData")
data(
  "TRNSCRPT_HEART_NORM_DATA",
  package = "MotrpacRatTraining6moData"
)
data("FEATURE_TO_GENE", package = "MotrpacRatTraining6moData")

norm_data <- TRNSCRPT_HEART_NORM_DATA

pheno_ids <- unique(c(
  rownames(PHENO),
  if ("viallabel" %in% colnames(PHENO)) {
    as.character(PHENO$viallabel)
  } else {
    character()
  }
))

sample_ids <- intersect(colnames(norm_data), pheno_ids)

if (length(sample_ids) == 0) {
  stop("No heart RNA-seq sample columns matched PHENO.")
}

if ("viallabel" %in% colnames(PHENO)) {
  sample_metadata <- PHENO |>
    filter(as.character(viallabel) %in% sample_ids) |>
    distinct(viallabel, .keep_all = TRUE)

  sample_metadata <- sample_metadata[
    match(sample_ids, as.character(sample_metadata$viallabel)),
    ,
    drop = FALSE
  ]
} else {
  sample_metadata <- PHENO[sample_ids, , drop = FALSE]
  sample_metadata$viallabel <- rownames(sample_metadata)
}

if (anyNA(match(sample_ids, as.character(sample_metadata$viallabel)))) {
  stop("At least one RNA sample could not be matched to PHENO.")
}

sample_metadata <- sample_metadata |>
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
    training_week = recode(
      as.character(timepoint),
      "SED" = 0,
      "1W" = 1,
      "2W" = 2,
      "4W" = 4,
      "8W" = 8
    ) |> as.numeric(),
    trained = as.numeric(timepoint != "SED"),
    sex_male = as.numeric(sex == "Male"),
    sex_x_training = sex_male * training_week,
    exp_group = interaction(
      substr(as.character(sex), 1, 1),
      timepoint,
      sep = "_"
    ),
    exp_group = factor(
      exp_group,
      levels = c(
        "F_SED", "F_1W", "F_2W", "F_4W", "F_8W",
        "M_SED", "M_1W", "M_2W", "M_4W", "M_8W"
      )
    )
  )

official_outlier_ids <- intersect(
  as.character(OUTLIERS$viallabel),
  sample_metadata$viallabel
)

sample_metadata <- sample_metadata |>
  filter(!viallabel %in% official_outlier_ids)

sample_ids <- sample_metadata$viallabel
rownames(sample_metadata) <- sample_ids

# ------------------------------------------------------------
# 5. Build the log2 ExpressionSet expected by official run_WGCNA()
# ------------------------------------------------------------
if (!"feature_ID" %in% colnames(norm_data)) {
  stop("TRNSCRPT_HEART_NORM_DATA does not contain feature_ID.")
}

expression_matrix <- norm_data |>
  select(feature_ID, all_of(sample_ids)) |>
  column_to_rownames("feature_ID") |>
  as.matrix()

storage.mode(expression_matrix) <- "numeric"

if (anyDuplicated(rownames(expression_matrix))) {
  stop("Duplicate feature_ID values were detected.")
}

if (anyNA(expression_matrix) || any(!is.finite(expression_matrix))) {
  stop(
    "Missing or non-finite normalized values were detected. ",
    "The official run_WGCNA implementation does not impute them."
  )
}

feature_sd <- apply(expression_matrix, 1, sd)

if (any(is.na(feature_sd) | feature_sd == 0)) {
  bad_features <- names(feature_sd)[is.na(feature_sd) | feature_sd == 0]

  stop(
    "Zero-variance or invalid features were detected: ",
    paste(head(bad_features, 20), collapse = ", "),
    ". Remove or resolve these features before WGCNA."
  )
}

# ------------------------------------------------------------
# Project-specific computational adaptation:
# select the 5,000 highest-variance normalized transcripts.
#
# IMPORTANT: This prefilter is not performed inside the official
# run_WGCNA() helper. After this point, the official function and
# all of its network parameters are used unchanged.
# ------------------------------------------------------------
n_features_before_selection <- nrow(expression_matrix)
feature_variance <- apply(expression_matrix, 1, var)

n_features_requested <- 5000L
n_features_selected <- min(
  n_features_requested,
  n_features_before_selection
)

selected_feature_ids <- names(
  sort(feature_variance, decreasing = TRUE)
)[seq_len(n_features_selected)]

selected_feature_table <- tibble(
  feature_ID = selected_feature_ids,
  expression_variance = unname(
    feature_variance[selected_feature_ids]
  ),
  variance_rank = seq_len(n_features_selected)
)

write_csv(
  selected_feature_table,
  file.path(
    output_dir,
    "selected_top5000_variable_features.csv"
  )
)

expression_matrix <- expression_matrix[
  selected_feature_ids,
  ,
  drop = FALSE
]

message(
  "Selected ",
  nrow(expression_matrix),
  " highest-variance transcripts from ",
  n_features_before_selection,
  " normalized heart transcript features."
)

# Preserve all metadata columns supplied with the normalized object.
feature_columns <- setdiff(colnames(norm_data), sample_ids)

feature_metadata <- norm_data |>
  select(all_of(feature_columns)) |>
  as.data.frame()

rownames(feature_metadata) <- feature_metadata$feature_ID

# Add/standardize gene annotations.
annotation_summary <- FEATURE_TO_GENE |>
  filter(feature_ID %in% rownames(expression_matrix)) |>
  select(feature_ID, gene_symbol, entrez_gene) |>
  distinct() |>
  group_by(feature_ID) |>
  summarise(
    annotation_gene_symbol = ifelse(
      all(is.na(gene_symbol)),
      NA_character_,
      paste(unique(na.omit(gene_symbol)), collapse = ";")
    ),
    annotation_entrez_gene = ifelse(
      all(is.na(entrez_gene)),
      NA_character_,
      paste(
        unique(na.omit(as.character(entrez_gene))),
        collapse = ";"
      )
    ),
    .groups = "drop"
  )

feature_metadata <- feature_metadata |>
  rownames_to_column("feature_rowname") |>
  left_join(annotation_summary, by = "feature_ID") |>
  column_to_rownames("feature_rowname")

# ExpressionSet requires sample metadata and feature metadata row names
# to match the assay matrix exactly.
sample_metadata <- sample_metadata[
  colnames(expression_matrix),
  ,
  drop = FALSE
]

feature_metadata <- feature_metadata[
  rownames(expression_matrix),
  ,
  drop = FALSE
]

heart_eset <- ExpressionSet(
  assayData = expression_matrix,
  phenoData = AnnotatedDataFrame(
    as.data.frame(sample_metadata)
  ),
  featureData = AnnotatedDataFrame(
    as.data.frame(feature_metadata)
  )
)

# ------------------------------------------------------------
# 6. Memory estimate
# ------------------------------------------------------------
n_features <- nrow(expression_matrix)
one_dense_matrix_gib <- (
  n_features * n_features * 8
) / (1024^3)

memory_note <- tibble(
  metric = c(
    "Normalized features before variance selection",
    "Features supplied to official run_WGCNA",
    "Samples supplied to official run_WGCNA",
    "Approximate size of one dense feature-by-feature matrix (GiB)",
    "Official outliers removed"
  ),
  value = c(
    n_features_before_selection,
    n_features,
    ncol(expression_matrix),
    round(one_dense_matrix_gib, 2),
    length(official_outlier_ids)
  )
)

write_csv(
  memory_note,
  file.path(output_dir, "WGCNA_memory_estimate.csv")
)

message(
  "The top-5,000 GitHub-matched WGCNA will create multiple dense ",
  n_features, " x ", n_features,
  " matrices. One such matrix is approximately ",
  round(one_dense_matrix_gib, 2), " GiB."
)

# ------------------------------------------------------------
# 7. Run the official GitHub implementation unchanged
# ------------------------------------------------------------
# Exact official settings:
# power = 1:20
# RsquaredCut = 0.90
# signed adjacency
# bicor, pairwise.complete.obs
# unsigned TOM
# cutHeight = 0.99
# deepSplit = 2
# minClusterSize = 30
# merge cutHeight = 0.15
# module prefix T for transcriptomics

run_official_wgcna_with_pdf <- function() {
  pdf(
    file.path(plot_dir, "official_run_WGCNA_diagnostics.pdf"),
    width = 13,
    height = 8,
    onefile = TRUE
  )

  on.exit(dev.off(), add = TRUE)

  run_WGCNA(
    object = heart_eset,
    power = 1:20,
    RsquaredCut = 0.90,
    module_prefix = "T"
  )
}

official_result_file <- file.path(
  output_dir,
  "official_run_WGCNA_result.rds"
)

if (file.exists(official_result_file)) {
  message(
    "Found a saved official WGCNA result. ",
    "Skipping network reconstruction and resuming downstream analysis."
  )

  official_result <- readRDS(
    official_result_file
  )
} else {
  message(
    "Starting official run_WGCNA() on the top-5,000 ",
    "variance-filtered heart features..."
  )

  official_result <- run_official_wgcna_with_pdf()

  # Save immediately so downstream errors never require
  # repeating adjacency/TOM construction.
  saveRDS(
    official_result,
    official_result_file,
    compress = "xz"
  )
}

# ------------------------------------------------------------
# 8. Save official outputs without modification
# ------------------------------------------------------------
official_modules <- official_result$modules
official_MEs_long <- official_result$MEs

if (!"feature_ID" %in% colnames(official_modules)) {
  official_modules$feature_ID <- rownames(official_modules)
}

write_csv(
  official_modules,
  file.path(output_dir, "official_WGCNA_modules.csv")
)

write_csv(
  official_MEs_long,
  file.path(output_dir, "official_WGCNA_module_eigengenes_long.csv")
)

# The result is saved immediately after construction above.
# Rewriting it here is harmless and keeps this section self-contained.
saveRDS(
  official_result,
  official_result_file,
  compress = "xz"
)

# ------------------------------------------------------------
# 9. Downstream analysis: module sizes and eigengene matrix
# ------------------------------------------------------------
module_sizes <- official_modules |>
  count(moduleID, moduleColor, name = "n_features") |>
  arrange(desc(n_features))

write_csv(
  module_sizes,
  file.path(output_dir, "module_sizes.csv")
)

ME_wide <- official_MEs_long |>
  mutate(
    viallabel = as.character(viallabel),
    moduleID = as.character(moduleID)
  ) |>
  select(viallabel, moduleID, ME) |>
  distinct() |>
  pivot_wider(
    names_from = moduleID,
    values_from = ME
  )

ME_matrix <- ME_wide |>
  column_to_rownames("viallabel") |>
  as.matrix()

storage.mode(ME_matrix) <- "numeric"

write_csv(
  ME_wide,
  file.path(output_dir, "module_eigengenes_wide.csv")
)

# ------------------------------------------------------------
# 10. Downstream analysis: bicor module membership and hub genes
# ------------------------------------------------------------
datExpr <- t(expression_matrix)
datExpr <- datExpr[rownames(ME_matrix), , drop = FALSE]

module_membership_list <- list()
membership_index <- 1

module_ids <- unique(as.character(official_modules$moduleID))

for (current_module in module_ids) {
  current_features <- official_modules$feature_ID[
    as.character(official_modules$moduleID) == current_module
  ]

  current_features <- intersect(
    current_features,
    colnames(datExpr)
  )

  if (
    length(current_features) == 0 ||
    !current_module %in% colnames(ME_matrix)
  ) {
    next
  }

  membership_values <- WGCNA::bicor(
    datExpr[, current_features, drop = FALSE],
    ME_matrix[, current_module],
    use = "pairwise.complete.obs"
  )

  membership_values <- as.numeric(membership_values)

  module_membership_list[[membership_index]] <- tibble(
    feature_ID = current_features,
    moduleID = current_module,
    module_membership_bicor = membership_values
  )

  membership_index <- membership_index + 1
}

module_membership <- bind_rows(module_membership_list) |>
  left_join(
    official_modules |>
      select(
        feature_ID,
        moduleID,
        moduleColor,
        everything()
      ),
    by = c("feature_ID", "moduleID")
  )

write_csv(
  module_membership,
  file.path(output_dir, "module_membership_bicor.csv")
)

hub_genes <- module_membership |>
  filter(moduleColor != "grey") |>
  group_by(moduleID, moduleColor) |>
  arrange(
    desc(abs(module_membership_bicor)),
    .by_group = TRUE
  ) |>
  slice_head(n = 10) |>
  ungroup()

write_csv(
  hub_genes,
  file.path(output_dir, "top10_hub_genes_per_module.csv")
)

for (current_module in unique(module_membership$moduleID)) {
  write_csv(
    module_membership |>
      filter(moduleID == current_module) |>
      arrange(desc(abs(module_membership_bicor))),
    file.path(
      module_dir,
      paste0("module_", current_module, "_features.csv")
    )
  )
}

# ------------------------------------------------------------
# 11. Downstream analysis: module–trait correlations
# ------------------------------------------------------------
# PHENO may retain data.table class after filtering. Convert explicitly
# to a base data.frame and align samples by viallabel rather than using
# data.table's character-i join syntax.
sample_metadata_for_me <- as.data.frame(
  sample_metadata
)

metadata_match <- match(
  rownames(ME_matrix),
  as.character(sample_metadata_for_me$viallabel)
)

if (anyNA(metadata_match)) {
  missing_samples <- rownames(ME_matrix)[
    is.na(metadata_match)
  ]

  stop(
    "Module eigengene samples were not found in metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

sample_metadata_for_me <- sample_metadata_for_me[
  metadata_match,
  ,
  drop = FALSE
]

rownames(sample_metadata_for_me) <- as.character(
  sample_metadata_for_me$viallabel
)

if (!identical(
  rownames(sample_metadata_for_me),
  rownames(ME_matrix)
)) {
  stop(
    "Sample metadata could not be aligned to the module eigengene matrix."
  )
}

first_numeric <- function(x) {
  values <- suppressWarnings(as.numeric(as.character(x)))
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    return(NA_real_)
  }

  values[1]
}

trait_candidates <- list(
  vo2max = c(
    "vo2_max",
    "vo2.max.test.vo2_max",
    "vo2_max_test___vo2_max"
  ),
  vo2max_change = c(
    "calculated.variables.vo2_max_change",
    "calculated_variables___vo2_max_change"
  ),
  pct_body_fat_change = c(
    "calculated.variables.pct_body_fat_change",
    "calculated_variables___pct_body_fat_change"
  ),
  pct_body_lean_change = c(
    "calculated.variables.pct_body_lean_change",
    "calculated_variables___pct_body_lean_change"
  ),
  terminal_body_weight = c(
    "terminal.weight.bw",
    "terminal_weight___bw"
  ),
  nmr_fat = c(
    "nmr.testing.nmr_fat",
    "nmr_testing___nmr_fat"
  )
)

trait_data <- data.frame(
  sex_male = sample_metadata_for_me$sex_male,
  training_week = sample_metadata_for_me$training_week,
  trained = sample_metadata_for_me$trained,
  sex_x_training = sample_metadata_for_me$sex_x_training,
  row.names = rownames(ME_matrix)
)

available_phenotype_sources <- list()

for (trait_name in names(trait_candidates)) {
  matching_columns <- trait_candidates[[trait_name]][
    trait_candidates[[trait_name]] %in% colnames(PHENO)
  ]

  if (length(matching_columns) == 0) {
    next
  }

  selected_column <- matching_columns[1]

  pid_trait <- PHENO |>
    filter(pid %in% sample_metadata_for_me$pid) |>
    group_by(pid) |>
    summarise(
      value = first_numeric(.data[[selected_column]]),
      .groups = "drop"
    )

  values <- pid_trait$value[
    match(sample_metadata_for_me$pid, pid_trait$pid)
  ]

  if (
    sum(is.finite(values)) >= 8 &&
    sd(values, na.rm = TRUE) > 0
  ) {
    trait_data[[trait_name]] <- values
    available_phenotype_sources[[trait_name]] <- selected_column
  }
}

write_csv(
  as.data.frame(trait_data) |>
    rownames_to_column("viallabel"),
  file.path(output_dir, "module_trait_data.csv")
)

if (length(available_phenotype_sources) > 0) {
  write_csv(
    tibble(
      trait = names(available_phenotype_sources),
      PHENO_column = unlist(available_phenotype_sources)
    ),
    file.path(output_dir, "phenotype_trait_source_columns.csv")
  )
}

correlation_results <- list()
correlation_index <- 1

for (current_module in colnames(ME_matrix)) {
  for (current_trait in colnames(trait_data)) {
    x <- ME_matrix[, current_module]
    y <- trait_data[, current_trait]

    complete_cases <- is.finite(x) & is.finite(y)
    n_complete <- sum(complete_cases)

    if (
      n_complete >= 8 &&
      sd(x[complete_cases]) > 0 &&
      sd(y[complete_cases]) > 0
    ) {
      test <- suppressWarnings(
        cor.test(
          x[complete_cases],
          y[complete_cases],
          method = "spearman",
          exact = FALSE
        )
      )

      correlation_results[[correlation_index]] <- tibble(
        moduleID = current_module,
        trait = current_trait,
        rho = unname(test$estimate),
        p_value = test$p.value,
        n = n_complete
      )

      correlation_index <- correlation_index + 1
    }
  }
}

module_trait_correlations <- bind_rows(correlation_results) |>
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significant_FDR05 = FDR < 0.05
  ) |>
  left_join(
    module_sizes |>
      select(moduleID, moduleColor, n_features),
    by = "moduleID"
  )

write_csv(
  module_trait_correlations,
  file.path(output_dir, "module_trait_correlations.csv")
)

write_csv(
  module_trait_correlations |>
    filter(significant_FDR05) |>
    arrange(FDR),
  file.path(output_dir, "significant_module_trait_correlations.csv")
)

# Module–trait heatmap
trait_plot_data <- module_trait_correlations |>
  mutate(
    module_label = paste0(
      moduleID,
      " / ",
      moduleColor,
      " (n=",
      n_features,
      ")"
    ),
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

module_order <- trait_plot_data |>
  group_by(module_label) |>
  summarise(
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(max_abs_rho)) |>
  pull(module_label)

trait_plot_data <- trait_plot_data |>
  mutate(
    module_label = factor(
      module_label,
      levels = rev(module_order)
    )
  )

module_trait_heatmap <- ggplot(
  trait_plot_data,
  aes(
    x = trait,
    y = module_label,
    fill = rho
  )
) +
  geom_tile(colour = "white") +
  geom_text(
    aes(label = significance),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    title = "Official WGCNA heart modules and sample traits",
    subtitle = paste(
      "Spearman correlations;",
      "* FDR<0.05, ** FDR<0.01, *** FDR<0.001"
    ),
    x = NULL,
    y = "Module",
    fill = "Spearman\nrho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(plot_dir, "module_trait_correlation_heatmap.png"),
  module_trait_heatmap,
  width = 11,
  height = max(
    6,
    0.35 * length(unique(trait_plot_data$moduleID))
  ),
  dpi = 300,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 12. Downstream analysis: eigengene profiles by sex/time
# ------------------------------------------------------------
ME_group_data <- ME_wide |>
  left_join(
    sample_metadata_for_me |>
      as.data.frame() |>
      rownames_to_column("metadata_row") |>
      select(
        viallabel,
        sex,
        timepoint,
        exp_group
      ),
    by = "viallabel"
  ) |>
  pivot_longer(
    cols = all_of(colnames(ME_matrix)),
    names_to = "moduleID",
    values_to = "ME"
  )

ME_group_summary <- ME_group_data |>
  group_by(
    moduleID,
    sex,
    timepoint,
    exp_group
  ) |>
  summarise(
    mean_ME = mean(ME),
    SE = sd(ME) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) |>
  group_by(moduleID) |>
  mutate(
    scaled_mean_ME = as.numeric(scale(mean_ME))
  ) |>
  ungroup()

write_csv(
  ME_group_summary,
  file.path(output_dir, "module_eigengene_group_summary.csv")
)

eigengene_heatmap <- ggplot(
  ME_group_summary,
  aes(
    x = exp_group,
    y = fct_rev(moduleID),
    fill = scaled_mean_ME
  )
) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    title = "Module eigengene profiles across sex and training time",
    x = NULL,
    y = "Official WGCNA module",
    fill = "Scaled\nmean ME"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(plot_dir, "module_eigengene_group_heatmap.png"),
  eigengene_heatmap,
  width = 11,
  height = max(6, 0.30 * ncol(ME_matrix)),
  dpi = 300,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 13. Downstream extension: same 13 contrasts at module level
# ------------------------------------------------------------
qc_columns <- c(
  "viallabel",
  "RIN",
  "pct_globin",
  "pct_umi_dup",
  "median_5_3_bias"
)

analysis_metadata <- sample_metadata_for_me

if (all(qc_columns %in% colnames(TRNSCRPT_META))) {
  qc_data <- TRNSCRPT_META |>
    transmute(
      viallabel = as.character(viallabel),
      rin = as.numeric(RIN),
      pct_globin = as.numeric(pct_globin),
      pct_umi_dup = as.numeric(pct_umi_dup),
      median_5_3_bias = as.numeric(median_5_3_bias)
    ) |>
    distinct(viallabel, .keep_all = TRUE)

  analysis_metadata <- analysis_metadata |>
    left_join(qc_data, by = "viallabel")

  standardize_covariate <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)

    if (sd(x) == 0) {
      return(rep(0, length(x)))
    }

    as.numeric(scale(x))
  }

  for (column_name in c(
    "rin",
    "pct_globin",
    "pct_umi_dup",
    "median_5_3_bias"
  )) {
    analysis_metadata[[column_name]] <- standardize_covariate(
      analysis_metadata[[column_name]]
    )
  }

  module_design <- model.matrix(
    ~ 0 + exp_group + rin + pct_globin +
      pct_umi_dup + median_5_3_bias,
    data = analysis_metadata
  )
} else {
  module_design <- model.matrix(
    ~ 0 + exp_group,
    data = analysis_metadata
  )
}

colnames(module_design) <- sub(
  "^exp_group",
  "",
  colnames(module_design)
)

rownames(module_design) <- analysis_metadata$viallabel

ME_expression <- t(
  ME_matrix[
    analysis_metadata$viallabel,
    ,
    drop = FALSE
  ]
)

base_fit <- lmFit(
  ME_expression,
  module_design
)

contrast_groups <- list(
  trained_vs_SED = c(
    "F_1W - F_SED",
    "F_2W - F_SED",
    "F_4W - F_SED",
    "F_8W - F_SED",
    "M_1W - M_SED",
    "M_2W - M_SED",
    "M_4W - M_SED",
    "M_8W - M_SED"
  ),
  MvF_SED = c(
    "M_SED - F_SED"
  ),
  MvF_exercise_response = c(
    "(M_1W - M_SED) - (F_1W - F_SED)",
    "(M_2W - M_SED) - (F_2W - F_SED)",
    "(M_4W - M_SED) - (F_4W - F_SED)",
    "(M_8W - M_SED) - (F_8W - F_SED)"
  )
)

run_module_contrasts <- function(
  contrast_group_name,
  contrast_vector
) {
  contrast_matrix <- makeContrasts(
    contrasts = contrast_vector,
    levels = module_design
  )

  fit <- contrasts.fit(
    base_fit,
    contrast_matrix
  )

  fit <- eBayes(
    fit,
    robust = TRUE
  )

  result_list <- lapply(
    seq_along(contrast_vector),
    function(i) {
      result <- topTable(
        fit,
        coef = i,
        number = Inf,
        sort.by = "none",
        adjust.method = "BH"
      )

      result$moduleID <- rownames(result)
      result$contrast <- contrast_vector[i]
      result$contrast_group <- contrast_group_name
      result
    }
  )

  results <- bind_rows(result_list)

  results$global_FDR <- p.adjust(
    results$P.Value,
    method = "BH"
  )

  results |>
    mutate(
      significant_global_FDR05 = global_FDR < 0.05
    )
}

module_contrast_results <- bind_rows(
  lapply(
    names(contrast_groups),
    function(group_name) {
      run_module_contrasts(
        group_name,
        contrast_groups[[group_name]]
      )
    }
  )
)

module_contrast_counts <- module_contrast_results |>
  group_by(contrast_group, contrast) |>
  summarise(
    tested_modules = n(),
    significant_modules = sum(
      significant_global_FDR05,
      na.rm = TRUE
    ),
    positive_modules = sum(
      significant_global_FDR05 & logFC > 0,
      na.rm = TRUE
    ),
    negative_modules = sum(
      significant_global_FDR05 & logFC < 0,
      na.rm = TRUE
    ),
    minimum_global_FDR = min(
      global_FDR,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  module_contrast_results,
  file.path(output_dir, "module_contrast_results.csv")
)

write_csv(
  module_contrast_counts,
  file.path(output_dir, "module_contrast_counts.csv")
)

# ------------------------------------------------------------
# 14. Inventory and complete reusable object
# ------------------------------------------------------------
analysis_inventory <- tibble(
  item = c(
    "Heart RNA-seq samples",
    "Normalized heart transcript features available",
    "Heart transcript features supplied",
    "Feature selection before WGCNA",
    "Network-construction implementation",
    "Power candidates",
    "Scale-free R2 threshold",
    "Adjacency type",
    "Correlation",
    "TOM type",
    "Dynamic-tree cut height",
    "Deep split",
    "Minimum cluster size",
    "Module merge cut height",
    "Module prefix",
    "Detected modules excluding grey",
    "Official source SHA-256",
    "WGCNA package version",
    "Data package version"
  ),
  value = c(
    ncol(expression_matrix),
    n_features_before_selection,
    nrow(expression_matrix),
    "Top 5,000 features by across-sample variance; project-specific computational adaptation",
    "Official MoTrPAC run_WGCNA source, unchanged after variance prefilter",
    "1:20",
    "0.90",
    "signed",
    "bicor; pairwise.complete.obs",
    "unsigned",
    "0.99",
    "2",
    "30",
    "0.15",
    "T",
    sum(module_sizes$moduleColor != "grey"),
    official_source_sha256,
    as.character(packageVersion("WGCNA")),
    as.character(
      packageVersion("MotrpacRatTraining6moData")
    )
  )
)

write_csv(
  analysis_inventory,
  file.path(output_dir, "official_WGCNA_analysis_inventory.csv")
)

saveRDS(
  list(
    official_source_sha256 = official_source_sha256,
    heart_eset = heart_eset,
    official_result = official_result,
    ME_matrix = ME_matrix,
    module_membership = module_membership,
    hub_genes = hub_genes,
    trait_data = trait_data,
    module_trait_correlations = module_trait_correlations,
    module_contrast_results = module_contrast_results
  ),
  file.path(output_dir, "rat_heart_GitHub_WGCNA_top5000_complete_objects.rds"),
  compress = "xz"
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\n================ OFFICIAL WGCNA INVENTORY ================\n")
print(analysis_inventory)

cat("\n================ MODULE SIZES =============================\n")
print(module_sizes)

cat("\n================ SIGNIFICANT MODULE–TRAIT RESULTS =========\n")
print(
  module_trait_correlations |>
    filter(significant_FDR05) |>
    arrange(FDR)
)

cat("\n================ MODULE CONTRAST COUNTS ===================\n")
print(module_contrast_counts)

cat("\nSUCCESS: GitHub-matched top-5,000 WGCNA outputs saved to:\n")
cat(normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
