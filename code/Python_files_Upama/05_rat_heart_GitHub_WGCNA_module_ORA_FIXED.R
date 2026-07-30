# ============================================================
# MoTrPAC Rat Heart Hackathon
# Step 5: GitHub-aligned WGCNA module functional enrichment
#
# Official workflow adapted from:
#   MotrpacRatTraining6moWATData/vignettes/WGCNA_MODULE_ORA.Rmd
#
# Official helper sourced unchanged:
#   MotrpacRatTraining6moWAT/R/fora2.R
#
# Analysis:
#   - Transcript Entrez genes grouped by WGCNA module
#   - Background = all unique Entrez genes represented in
#     non-grey WGCNA modules
#   - Same MSigDB v7.5.1 GO pathways used in Step 3
#   - Step 3 pathways already passed the >=85% retention filter
#     against the full detected transcriptome
#   - For the project-specific top-5,000 WGCNA universe:
#       intersect pathways with mapped WGCNA Entrez genes
#       require post-intersection size >= 15
#
# Important adaptation:
#   The >=85% filter is NOT applied a second time against the
#   top-5,000 variance subset. Reapplying it would remove every
#   pathway because the WGCNA subset is intentionally smaller
#   than the detected transcriptome.
#   - Official fora2():
#       adjust.method = "scale"
#       adjust.globally = TRUE
#
# Important:
#   ORA is non-directional. It identifies biological functions
#   over-represented within each co-expression module.
# ============================================================

options(timeout = 1e5)
options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(0)

message("R version: ", R.version.string)
message("Working directory: ", getwd())

# ------------------------------------------------------------
# 1. Install and load required packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "dplyr", "tidyr", "tibble", "purrr", "readr",
  "ggplot2", "stringr", "forcats", "scales",
  "data.table", "digest", "assertthat"
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

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install(
    "fgsea",
    ask = FALSE,
    update = FALSE,
    Ncpus = 2
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
  library(assertthat)
  library(fgsea)
})

message("fgsea version: ", as.character(packageVersion("fgsea")))

# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------
wgcna_dir <- file.path(
  "results",
  "04_rat_heart_GitHub_WGCNA_top5000"
)

step3_dir <- file.path(
  "results",
  "03_rat_heart_GitHub_FGSEA"
)

step2_dir <- file.path(
  "results",
  "02_rat_heart_differential_expression"
)

modules_file <- file.path(
  wgcna_dir,
  "official_WGCNA_modules.csv"
)

trait_file <- file.path(
  wgcna_dir,
  "significant_module_trait_correlations.csv"
)

de_file <- file.path(
  step2_dir,
  "all_differential_expression_results.csv"
)

step3_pathway_file <- file.path(
  step3_dir,
  "filtered_MSigDB_GO_pathways_85pct.csv"
)

step3_fgsea_file <- file.path(
  step3_dir,
  "all_GitHub_aligned_FGSEA_GO_results.csv"
)

output_dir <- file.path(
  "results",
  "05_rat_heart_GitHub_WGCNA_module_ORA"
)

plot_dir <- file.path(output_dir, "plots")
module_table_dir <- file.path(output_dir, "module_tables")
source_dir <- file.path(output_dir, "official_source_snapshots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(module_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- c(
  modules_file,
  de_file,
  step3_pathway_file,
  step3_fgsea_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Required input files were not found:\n",
    paste(
      normalizePath(
        missing_input_files,
        winslash = "/",
        mustWork = FALSE
      ),
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 3. Download and source the official helper unchanged
# ------------------------------------------------------------
official_sources <- tibble(
  source_name = c(
    "fora2.R",
    "WGCNA_MODULE_ORA.Rmd"
  ),
  source_url = c(
    paste0(
      "https://raw.githubusercontent.com/",
      "MoTrPAC/MotrpacRatTraining6moWAT/",
      "master/R/fora2.R"
    ),
    paste0(
      "https://raw.githubusercontent.com/",
      "MoTrPAC/MotrpacRatTraining6moWATData/",
      "master/vignettes/WGCNA_MODULE_ORA.Rmd"
    )
  )
)

for (i in seq_len(nrow(official_sources))) {
  destination <- file.path(
    source_dir,
    official_sources$source_name[i]
  )

  if (!file.exists(destination)) {
    download.file(
      official_sources$source_url[i],
      destination,
      mode = "wb",
      quiet = FALSE
    )
  }
}

source(
  file.path(source_dir, "fora2.R"),
  local = .GlobalEnv
)

if (!exists("fora2", mode = "function")) {
  stop("The official fora2() function was not loaded.")
}

official_sources <- official_sources |>
  mutate(
    local_file = file.path(
      source_dir,
      source_name
    ),
    sha256 = map_chr(
      local_file,
      ~ digest(.x, algo = "sha256", file = TRUE)
    )
  )

write_csv(
  official_sources,
  file.path(output_dir, "official_source_records.csv")
)

# ------------------------------------------------------------
# 4. Read module, annotation and Step 3 pathway data
# ------------------------------------------------------------
modules_raw <- read_csv(
  modules_file,
  show_col_types = FALSE,
  progress = FALSE
)

de_results <- read_csv(
  de_file,
  show_col_types = FALSE,
  progress = FALSE
)

step3_pathways_raw <- read_csv(
  step3_pathway_file,
  show_col_types = FALSE,
  progress = FALSE
)

step3_fgsea <- read_csv(
  step3_fgsea_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_module_columns <- c(
  "moduleID",
  "moduleColor"
)

missing_module_columns <- setdiff(
  required_module_columns,
  colnames(modules_raw)
)

if (length(missing_module_columns) > 0) {
  stop(
    "Missing WGCNA module columns: ",
    paste(missing_module_columns, collapse = ", ")
  )
}

feature_column <- if (
  "feature_ID" %in% colnames(modules_raw)
) {
  "feature_ID"
} else if (
  "feature" %in% colnames(modules_raw)
) {
  "feature"
} else {
  stop(
    "The WGCNA module table must contain feature_ID or feature."
  )
}

required_de_columns <- c(
  "feature",
  "entrez_gene",
  "gene_symbol"
)

missing_de_columns <- setdiff(
  required_de_columns,
  colnames(de_results)
)

if (length(missing_de_columns) > 0) {
  stop(
    "Missing Step 2 annotation columns: ",
    paste(missing_de_columns, collapse = ", ")
  )
}

required_pathway_columns <- c(
  "gs_subcat",
  "gs_exact_source",
  "gs_description",
  "entrez_gene",
  "set_size"
)

missing_pathway_columns <- setdiff(
  required_pathway_columns,
  colnames(step3_pathways_raw)
)

if (length(missing_pathway_columns) > 0) {
  stop(
    "Missing Step 3 pathway columns: ",
    paste(missing_pathway_columns, collapse = ", ")
  )
}

required_fgsea_columns <- c(
  "pathway",
  "gs_name"
)

missing_fgsea_columns <- setdiff(
  required_fgsea_columns,
  colnames(step3_fgsea)
)

if (length(missing_fgsea_columns) > 0) {
  stop(
    "Missing Step 3 pathway-name columns: ",
    paste(missing_fgsea_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. Prepare transcript-to-Entrez mappings
# ------------------------------------------------------------
collapse_nonempty_unique <- function(x) {
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

de_feature_annotation <- de_results |>
  transmute(
    feature_ID = as.character(feature),
    entrez_gene = as.character(entrez_gene)
  ) |>
  group_by(feature_ID) |>
  summarise(
    de_entrez_gene = collapse_nonempty_unique(
      entrez_gene
    ),
    .groups = "drop"
  )

# Build an Entrez-to-symbol lookup using the same one-to-many
# separation approach used in Step 3.
entrez_to_symbol_table <- de_results |>
  transmute(
    entrez_gene = as.character(entrez_gene),
    gene_symbol = as.character(gene_symbol)
  ) |>
  separate_rows(
    entrez_gene,
    gene_symbol,
    sep = ";"
  ) |>
  mutate(
    entrez_gene = str_trim(entrez_gene),
    gene_symbol = str_trim(gene_symbol)
  ) |>
  filter(
    !is.na(entrez_gene),
    entrez_gene != "",
    entrez_gene != "NA"
  ) |>
  group_by(entrez_gene) |>
  summarise(
    gene_symbol = {
      valid_symbols <- gene_symbol[
        !is.na(gene_symbol) &
          gene_symbol != "" &
          gene_symbol != "NA"
      ]

      if (length(valid_symbols) == 0) {
        NA_character_
      } else {
        first(valid_symbols)
      }
    },
    .groups = "drop"
  )

entrez_to_symbol <- setNames(
  entrez_to_symbol_table$gene_symbol,
  entrez_to_symbol_table$entrez_gene
)

module_entrez_column <- if (
  "entrez_gene" %in% colnames(modules_raw)
) {
  "entrez_gene"
} else if (
  "annotation_entrez_gene" %in% colnames(modules_raw)
) {
  "annotation_entrez_gene"
} else {
  NA_character_
}

modules_prepared <- modules_raw |>
  transmute(
    feature_ID = as.character(
      .data[[feature_column]]
    ),
    moduleID = as.character(moduleID),
    moduleColor = as.character(moduleColor),
    module_entrez_gene = if (
      !is.na(module_entrez_column)
    ) {
      as.character(
        .data[[module_entrez_column]]
      )
    } else {
      NA_character_
    }
  ) |>
  left_join(
    de_feature_annotation,
    by = "feature_ID"
  ) |>
  mutate(
    module_entrez_gene = na_if(
      str_trim(module_entrez_gene),
      ""
    ),
    de_entrez_gene = na_if(
      str_trim(de_entrez_gene),
      ""
    ),
    entrez_gene = coalesce(
      module_entrez_gene,
      de_entrez_gene
    )
  )

n_unmapped_features <- sum(
  is.na(modules_prepared$entrez_gene) |
    modules_prepared$entrez_gene == "" |
    modules_prepared$entrez_gene == "NA"
)

message(
  "WGCNA features without an Entrez mapping: ",
  n_unmapped_features
)

TRNSCRPT_MODULES <- modules_prepared |>
  separate_rows(
    entrez_gene,
    sep = ";"
  ) |>
  mutate(
    entrez_gene = str_trim(
      as.character(entrez_gene)
    )
  ) |>
  filter(
    !is.na(entrez_gene),
    entrez_gene != "",
    entrez_gene != "NA",
    tolower(moduleColor) != "grey"
  ) |>
  distinct(
    feature_ID,
    moduleID,
    moduleColor,
    entrez_gene
  ) |>
  left_join(
    entrez_to_symbol_table,
    by = "entrez_gene"
  )

if (nrow(TRNSCRPT_MODULES) == 0) {
  stop(
    "No mapped Entrez genes remained after module preparation."
  )
}

module_order <- TRNSCRPT_MODULES |>
  distinct(moduleID) |>
  mutate(
    module_number = suppressWarnings(
      as.integer(
        str_extract(moduleID, "\\d+")
      )
    )
  ) |>
  arrange(
    is.na(module_number),
    module_number,
    moduleID
  ) |>
  pull(moduleID)

TRNSCRPT_MODULES <- TRNSCRPT_MODULES |>
  mutate(
    moduleID = factor(
      moduleID,
      levels = module_order
    )
  )

# Official vignette structure: one unique Entrez-gene vector per module.
TRNSCRPT_MOD_LIST <- TRNSCRPT_MODULES |>
  group_by(moduleID) |>
  summarise(
    feature = list(
      unique(entrez_gene)
    ),
    .groups = "drop"
  ) |>
  mutate(
    moduleID = as.character(moduleID)
  ) |>
  tibble::deframe()

TRNSCRPT_UNIVERSE <- unique(
  unlist(
    TRNSCRPT_MOD_LIST,
    use.names = FALSE
  )
)

if (length(TRNSCRPT_UNIVERSE) < 2) {
  stop(
    "The WGCNA Entrez-gene universe contains fewer than two genes."
  )
}

module_gene_counts <- TRNSCRPT_MODULES |>
  mutate(
    moduleID = as.character(moduleID)
  ) |>
  group_by(
    moduleID,
    moduleColor
  ) |>
  summarise(
    n_wgcna_features = n_distinct(feature_ID),
    n_unique_entrez_genes = n_distinct(entrez_gene),
    .groups = "drop"
  ) |>
  mutate(
    moduleID = factor(
      moduleID,
      levels = module_order
    )
  ) |>
  arrange(moduleID) |>
  mutate(
    moduleID = as.character(moduleID)
  )

write_csv(
  module_gene_counts,
  file.path(output_dir, "module_gene_counts.csv")
)

write_csv(
  TRNSCRPT_MODULES |>
    mutate(
      moduleID = as.character(moduleID)
    ) |>
    arrange(
      factor(moduleID, levels = module_order),
      entrez_gene
    ),
  file.path(output_dir, "module_entrez_gene_lists.csv")
)

for (current_module in module_order) {
  write_csv(
    TRNSCRPT_MODULES |>
      mutate(
        moduleID = as.character(moduleID)
      ) |>
      filter(
        moduleID == current_module
      ) |>
      arrange(
        entrez_gene,
        feature_ID
      ),
    file.path(
      module_table_dir,
      paste0(
        "module_",
        current_module,
        "_entrez_genes.csv"
      )
    )
  )
}

# ------------------------------------------------------------
# 6. Intersect Step 3 GO pathways with the WGCNA universe
# ------------------------------------------------------------
# The Step 3 pathway file already contains pathways that passed:
#   - original MSigDB size of 15–300 genes
#   - >=85% retention in the full detected heart transcriptome
#   - >=15 detected genes
#
# Step 4 used a project-specific top-5,000 variance prefilter for
# computational feasibility. This is not a biological detection
# universe. Therefore, the official >=85% retention criterion must
# not be applied again to this smaller subset.
#
# For module ORA, the correct adaptation is:
#   1. retain the Step 3-approved pathway collection;
#   2. intersect every pathway with mapped WGCNA Entrez genes;
#   3. require at least 15 genes after intersection;
#   4. use all mapped WGCNA Entrez genes as the ORA universe.

MSIGDB_PATHWAYS_FROM_STEP3 <- step3_pathways_raw |>
  transmute(
    gs_subcat = as.character(gs_subcat),
    gs_exact_source = as.character(gs_exact_source),
    gs_description = as.character(gs_description),
    entrez_gene = str_split(
      as.character(entrez_gene),
      pattern = ";"
    ),
    original_set_size = as.integer(set_size)
  ) |>
  mutate(
    entrez_gene = map(
      entrez_gene,
      ~ unique(
        str_trim(
          as.character(.x)
        )
      )
    ),
    step3_detected_set_size = lengths(entrez_gene)
  )

WGCNA_PATHWAY_DIAGNOSTICS <- MSIGDB_PATHWAYS_FROM_STEP3 |>
  mutate(
    entrez_gene_in_WGCNA = map(
      entrez_gene,
      intersect,
      y = TRNSCRPT_UNIVERSE
    ),
    WGCNA_set_size = lengths(entrez_gene_in_WGCNA),
    WGCNA_fraction_of_step3 = ifelse(
      step3_detected_set_size > 0,
      WGCNA_set_size / step3_detected_set_size,
      NA_real_
    ),
    passes_minimum_15 = WGCNA_set_size >= 15
  )

write_csv(
  WGCNA_PATHWAY_DIAGNOSTICS |>
    transmute(
      gs_subcat,
      gs_exact_source,
      gs_description,
      original_set_size,
      step3_detected_set_size,
      WGCNA_set_size,
      WGCNA_fraction_of_step3,
      passes_minimum_15
    ),
  file.path(
    output_dir,
    "WGCNA_pathway_intersection_diagnostics.csv"
  )
)

message(
  "WGCNA pathway intersection sizes: median=",
  median(
    WGCNA_PATHWAY_DIAGNOSTICS$WGCNA_set_size,
    na.rm = TRUE
  ),
  ", maximum=",
  max(
    WGCNA_PATHWAY_DIAGNOSTICS$WGCNA_set_size,
    na.rm = TRUE
  ),
  ", pathways with >=15 genes=",
  sum(
    WGCNA_PATHWAY_DIAGNOSTICS$passes_minimum_15,
    na.rm = TRUE
  )
)

TRNSCRPT_MSIGDB <- WGCNA_PATHWAY_DIAGNOSTICS |>
  filter(passes_minimum_15) |>
  transmute(
    gs_subcat,
    gs_exact_source,
    gs_description,
    entrez_gene = entrez_gene_in_WGCNA,
    set_size = original_set_size,
    set_size_step3 = step3_detected_set_size,
    set_size_post = WGCNA_set_size,
    ratio = WGCNA_fraction_of_step3
  )

if (nrow(TRNSCRPT_MSIGDB) == 0) {
  stop(
    paste0(
      "No GO pathway contained at least 15 mapped genes in the ",
      "top-5,000 WGCNA universe. Inspect ",
      "WGCNA_pathway_intersection_diagnostics.csv."
    )
  )
}

pathway_filter_summary <- full_join(
  MSIGDB_PATHWAYS_FROM_STEP3 |>
    count(
      gs_subcat,
      name = "Step3_approved_pathways"
    ),
  TRNSCRPT_MSIGDB |>
    count(
      gs_subcat,
      name = "pathways_with_at_least_15_WGCNA_genes"
    ),
  by = "gs_subcat"
) |>
  mutate(
    across(
      where(is.numeric),
      ~ replace_na(.x, 0)
    ),
    pathways_removed_for_small_WGCNA_overlap = (
      Step3_approved_pathways -
        pathways_with_at_least_15_WGCNA_genes
    )
  )

write_csv(
  pathway_filter_summary,
  file.path(output_dir, "pathway_filter_summary.csv")
)

write_csv(
  TRNSCRPT_MSIGDB |>
    mutate(
      entrez_gene = map_chr(
        entrez_gene,
        ~ paste(.x, collapse = ";")
      )
    ),
  file.path(
    output_dir,
    "WGCNA_intersected_MSigDB_GO_pathways_min15.csv"
  )
)

# Short pathway-name lookup from the completed Step 3 results.
pathway_name_lookup <- step3_fgsea |>
  transmute(
    pathway = as.character(pathway),
    gs_name = as.character(gs_name)
  ) |>
  filter(
    !is.na(pathway),
    pathway != "",
    !is.na(gs_name),
    gs_name != ""
  ) |>
  distinct(
    pathway,
    .keep_all = TRUE
  )

make_short_go_label <- function(
  gs_name,
  pathway,
  wrap_width = 42
) {
  label <- as.character(gs_name)
  pathway <- as.character(pathway)

  missing_name <- (
    is.na(label) |
      label == ""
  )

  label[missing_name] <- pathway[missing_name]

  label |>
    str_remove("^(GOBP_|GOMF_|GOCC_)") |>
    str_remove("^GO[:_ ]?(BP|MF|CC)[:_ ]?") |>
    str_replace_all("_", " ") |>
    str_squish() |>
    str_to_sentence() |>
    str_wrap(width = wrap_width)
}

# ------------------------------------------------------------
# 7. Official module ORA
# ------------------------------------------------------------
message(
  "Running official fora2() for ",
  length(TRNSCRPT_MOD_LIST),
  " modules, ",
  length(TRNSCRPT_UNIVERSE),
  " background genes and ",
  nrow(TRNSCRPT_MSIGDB),
  " GO pathways."
)

TRNSCRPT_MODULE_ORA_raw <- fora2(
  pathways = TRNSCRPT_MSIGDB,
  genes = TRNSCRPT_MOD_LIST,
  universe = TRNSCRPT_UNIVERSE,
  adjust.method = "scale",
  adjust.globally = TRUE
)

TRNSCRPT_MODULE_ORA <- as_tibble(
  TRNSCRPT_MODULE_ORA_raw
)

required_ora_columns <- c(
  "module",
  "pathway",
  "gs_subcat",
  "gs_description",
  "pval",
  "padj",
  "overlap",
  "size",
  "overlapGenes"
)

missing_ora_columns <- setdiff(
  required_ora_columns,
  colnames(TRNSCRPT_MODULE_ORA)
)

if (length(missing_ora_columns) > 0) {
  stop(
    "The official fora2() output is missing expected columns: ",
    paste(missing_ora_columns, collapse = ", ")
  )
}

# The official vignette reports maxOverlap and overlapRatio.
# Check explicitly so a future source-code change cannot pass silently.
required_scaled_columns <- c(
  "maxOverlap",
  "overlapRatio"
)

missing_scaled_columns <- setdiff(
  required_scaled_columns,
  colnames(TRNSCRPT_MODULE_ORA)
)

if (length(missing_scaled_columns) > 0) {
  stop(
    "The official fora2() output is missing scaled-P columns: ",
    paste(missing_scaled_columns, collapse = ", ")
  )
}

module_metadata <- module_gene_counts |>
  rename(
    module = moduleID
  )

TRNSCRPT_MODULE_ORA <- TRNSCRPT_MODULE_ORA |>
  mutate(
    module = as.character(module),
    pathway = as.character(pathway),
    gs_subcat = as.character(gs_subcat),
    scaled_pvalue = as.numeric(padj),
    significant_scaled_p05 = (
      !is.na(scaled_pvalue) &
        scaled_pvalue < 0.05
    ),
    minus_log10_scaled_p = -log10(
      pmax(
        scaled_pvalue,
        .Machine$double.xmin
      )
    ),
    overlapGenes_symbol = map(
      overlapGenes,
      ~ {
        symbols <- entrez_to_symbol[
          as.character(.x)
        ]
        unique(
          na.omit(
            as.character(symbols)
          )
        )
      }
    )
  ) |>
  left_join(
    pathway_name_lookup,
    by = "pathway"
  ) |>
  left_join(
    module_metadata,
    by = "module"
  ) |>
  mutate(
    pathway_label = make_short_go_label(
      gs_name,
      pathway,
      wrap_width = 42
    ),
    module = factor(
      module,
      levels = module_order
    )
  ) |>
  arrange(
    module,
    scaled_pvalue,
    desc(overlapRatio)
  ) |>
  mutate(
    module = as.character(module)
  )

# Preserve the full list-column object.
saveRDS(
  TRNSCRPT_MODULE_ORA,
  file.path(
    output_dir,
    "rat_heart_WGCNA_module_ORA_complete.rds"
  ),
  compress = "xz"
)

ora_csv <- TRNSCRPT_MODULE_ORA |>
  mutate(
    overlapGenes_entrez = map_chr(
      overlapGenes,
      ~ paste(
        as.character(.x),
        collapse = ";"
      )
    ),
    overlapGenes_symbol = map_chr(
      overlapGenes_symbol,
      ~ paste(
        as.character(.x),
        collapse = ";"
      )
    )
  ) |>
  select(-overlapGenes) |>
  relocate(
    module,
    moduleColor,
    gs_subcat,
    pathway,
    gs_name,
    pathway_label,
    gs_description
  ) |>
  relocate(
    overlap,
    maxOverlap,
    overlapRatio,
    .after = size
  ) |>
  relocate(
    scaled_pvalue,
    significant_scaled_p05,
    minus_log10_scaled_p,
    .after = padj
  )

write_csv(
  ora_csv,
  file.path(
    output_dir,
    "all_module_GO_ORA_results.csv"
  )
)

significant_ora <- ora_csv |>
  filter(significant_scaled_p05) |>
  arrange(
    factor(module, levels = module_order),
    scaled_pvalue,
    desc(overlapRatio)
  )

write_csv(
  significant_ora,
  file.path(
    output_dir,
    "significant_module_GO_ORA_results.csv"
  )
)

top10_terms_per_module <- significant_ora |>
  group_by(
    module,
    moduleColor
  ) |>
  arrange(
    scaled_pvalue,
    desc(overlapRatio),
    .by_group = TRUE
  ) |>
  slice_head(n = 10) |>
  ungroup()

write_csv(
  top10_terms_per_module,
  file.path(
    output_dir,
    "top10_GO_terms_per_module.csv"
  )
)

significant_counts <- ora_csv |>
  group_by(
    module,
    moduleColor,
    gs_subcat
  ) |>
  summarise(
    n_tested_pathways = n(),
    n_significant_pathways = sum(
      significant_scaled_p05,
      na.rm = TRUE
    ),
    minimum_scaled_pvalue = min(
      scaled_pvalue,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    minimum_scaled_pvalue = ifelse(
      is.infinite(minimum_scaled_pvalue),
      NA_real_,
      minimum_scaled_pvalue
    ),
    module = factor(
      module,
      levels = module_order
    )
  ) |>
  arrange(module, gs_subcat) |>
  mutate(
    module = as.character(module)
  )

write_csv(
  significant_counts,
  file.path(
    output_dir,
    "module_ORA_significant_counts.csv"
  )
)

module_ORA_summary <- significant_counts |>
  group_by(
    module,
    moduleColor
  ) |>
  summarise(
    n_significant_GO_BP = sum(
      n_significant_pathways[
        gs_subcat == "GO:BP"
      ],
      na.rm = TRUE
    ),
    n_significant_GO_CC = sum(
      n_significant_pathways[
        gs_subcat == "GO:CC"
      ],
      na.rm = TRUE
    ),
    n_significant_GO_MF = sum(
      n_significant_pathways[
        gs_subcat == "GO:MF"
      ],
      na.rm = TRUE
    ),
    n_significant_total = sum(
      n_significant_pathways,
      na.rm = TRUE
    ),
    minimum_scaled_pvalue = min(
      minimum_scaled_pvalue,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  left_join(
    module_gene_counts,
    by = c(
      "module" = "moduleID",
      "moduleColor"
    )
  ) |>
  mutate(
    minimum_scaled_pvalue = ifelse(
      is.infinite(minimum_scaled_pvalue),
      NA_real_,
      minimum_scaled_pvalue
    ),
    module = factor(
      module,
      levels = module_order
    )
  ) |>
  arrange(module) |>
  mutate(
    module = as.character(module)
  )

write_csv(
  module_ORA_summary,
  file.path(
    output_dir,
    "module_ORA_summary.csv"
  )
)

# ------------------------------------------------------------
# 8. Publication-ready figures
# ------------------------------------------------------------
save_empty_plot <- function(
  filename,
  title_text,
  subtitle_text
) {
  plot_object <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = subtitle_text,
      size = 5
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title_text) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  ggsave(
    file.path(plot_dir, filename),
    plot_object,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# Figure 1: significant counts by module and GO category.
count_plot_data <- significant_counts |>
  mutate(
    module = factor(
      module,
      levels = module_order
    )
  )

count_plot <- ggplot(
  count_plot_data,
  aes(
    x = module,
    y = n_significant_pathways,
    fill = gs_subcat
  )
) +
  geom_col() +
  labs(
    title = "Significant GO enrichment across rat-heart WGCNA modules",
    subtitle = paste(
      "Official fora2 scaled P < 0.05;",
      "counts are non-directional"
    ),
    x = "WGCNA module",
    y = "Number of significant pathways",
    fill = "GO category"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    plot_dir,
    "module_ORA_significant_counts.png"
  ),
  count_plot,
  width = 11,
  height = 7,
  dpi = 300
)

# Figure 2: top three significant GO:BP terms per module.
dotplot_data <- significant_ora |>
  filter(gs_subcat == "GO:BP") |>
  group_by(module) |>
  arrange(
    scaled_pvalue,
    desc(overlapRatio),
    .by_group = TRUE
  ) |>
  slice_head(n = 3) |>
  ungroup()

if (nrow(dotplot_data) > 0) {
  pathway_order_dot <- dotplot_data |>
    group_by(pathway_label) |>
    summarise(
      best_scaled_p = min(
        scaled_pvalue,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(desc(best_scaled_p)) |>
    pull(pathway_label)

  dotplot_data <- dotplot_data |>
    mutate(
      module = factor(
        module,
        levels = module_order
      ),
      pathway_label = factor(
        pathway_label,
        levels = unique(pathway_order_dot)
      )
    )

  dotplot_height <- max(
    8,
    0.32 * length(
      unique(dotplot_data$pathway_label)
    )
  )

  top_terms_dotplot <- ggplot(
    dotplot_data,
    aes(
      x = module,
      y = pathway_label,
      size = overlap,
      fill = minus_log10_scaled_p
    )
  ) +
    geom_point(
      shape = 21,
      colour = "black",
      alpha = 0.9
    ) +
    scale_size_continuous(
      range = c(2.5, 9)
    ) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick"
    ) +
    labs(
      title = "Top GO biological processes in each WGCNA module",
      subtitle = "Top three significant terms per module",
      x = "WGCNA module",
      y = NULL,
      size = "Overlap genes",
      fill = "-log10(scaled P)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(
      plot_dir,
      "module_ORA_top_terms_dotplot.png"
    ),
    top_terms_dotplot,
    width = 13,
    height = dotplot_height,
    dpi = 300,
    limitsize = FALSE
  )
} else {
  save_empty_plot(
    "module_ORA_top_terms_dotplot.png",
    "Top GO biological processes in each WGCNA module",
    "No GO:BP term reached scaled P < 0.05."
  )
}

make_module_pathway_heatmap <- function(
  input_data,
  selected_modules,
  title_text,
  filename,
  top_n_pathways = 25
) {
  input_selected <- input_data |>
    filter(
      module %in% selected_modules,
      gs_subcat == "GO:BP"
    )

  significant_candidates <- input_selected |>
    filter(significant_scaled_p05) |>
    group_by(
      pathway,
      pathway_label
    ) |>
    summarise(
      minimum_scaled_pvalue = min(
        scaled_pvalue,
        na.rm = TRUE
      ),
      maximum_overlap_ratio = max(
        overlapRatio,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(
      minimum_scaled_pvalue,
      desc(maximum_overlap_ratio)
    ) |>
    slice_head(n = top_n_pathways)

  if (nrow(significant_candidates) == 0) {
    save_empty_plot(
      filename,
      title_text,
      "No GO:BP term reached scaled P < 0.05."
    )
    return(invisible(NULL))
  }

  heatmap_data <- input_selected |>
    filter(
      pathway %in% significant_candidates$pathway
    ) |>
    mutate(
      module = factor(
        module,
        levels = selected_modules
      ),
      significance_mark = ifelse(
        significant_scaled_p05,
        "*",
        ""
      ),
      display_score = pmin(
        minus_log10_scaled_p,
        20
      )
    )

  pathway_order <- heatmap_data |>
    group_by(pathway_label) |>
    summarise(
      best_scaled_p = min(
        scaled_pvalue,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(best_scaled_p) |>
    pull(pathway_label)

  heatmap_data <- heatmap_data |>
    mutate(
      pathway_label = factor(
        pathway_label,
        levels = rev(pathway_order)
      )
    )

  heatmap_height <- max(
    7,
    0.34 * length(
      unique(heatmap_data$pathway_label)
    )
  )

  heatmap_plot <- ggplot(
    heatmap_data,
    aes(
      x = module,
      y = pathway_label,
      fill = display_score
    )
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = significance_mark),
      size = 4
    ) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick"
    ) +
    labs(
      title = title_text,
      subtitle = paste(
        "* indicates official fora2 scaled P < 0.05;",
        "fill capped at -log10(P) = 20"
      ),
      x = "WGCNA module",
      y = NULL,
      fill = "-log10(scaled P)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )

  ggsave(
    file.path(plot_dir, filename),
    heatmap_plot,
    width = 12,
    height = heatmap_height,
    dpi = 300,
    limitsize = FALSE
  )
}

# Figure 3: top GO:BP terms across all modules.
make_module_pathway_heatmap(
  input_data = ora_csv,
  selected_modules = module_order,
  title_text = "Top GO biological-process enrichment across WGCNA modules",
  filename = "module_ORA_top_GO_BP_heatmap.png",
  top_n_pathways = 30
)

# Figure 4: direct T1/T2 comparison.
t1_t2_modules <- intersect(
  c("T1", "T2"),
  module_order
)

if (length(t1_t2_modules) == 2) {
  make_module_pathway_heatmap(
    input_data = ora_csv,
    selected_modules = t1_t2_modules,
    title_text = "Functional comparison of the principal sex-associated modules",
    filename = "T1_T2_functional_comparison.png",
    top_n_pathways = 25
  )
} else {
  save_empty_plot(
    "T1_T2_functional_comparison.png",
    "Functional comparison of T1 and T2",
    "Both T1 and T2 were not available."
  )
}

# Figure 5: modules with strongest Step 4 trait associations.
selected_modules <- intersect(
  c("T1", "T2", "T5", "T10"),
  module_order
)

if (file.exists(trait_file)) {
  trait_results <- read_csv(
    trait_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  trait_module_column <- intersect(
    c("moduleID", "module", "Module"),
    colnames(trait_results)
  )

  trait_fdr_column <- intersect(
    c("FDR", "fdr", "padj", "adjusted_pvalue"),
    colnames(trait_results)
  )

  if (
    length(trait_module_column) > 0 &&
    length(trait_fdr_column) > 0
  ) {
    trait_selected <- trait_results |>
      transmute(
        module = as.character(
          .data[[trait_module_column[1]]]
        ),
        trait_fdr = as.numeric(
          .data[[trait_fdr_column[1]]]
        )
      ) |>
      filter(
        !is.na(module),
        module %in% module_order
      ) |>
      arrange(trait_fdr) |>
      distinct(module, .keep_all = TRUE) |>
      slice_head(n = 6) |>
      pull(module)

    selected_modules <- unique(
      c(
        intersect(c("T1", "T2"), module_order),
        trait_selected
      )
    )
  }
}

selected_modules <- selected_modules[
  selected_modules %in% module_order
]

if (length(selected_modules) > 0) {
  make_module_pathway_heatmap(
    input_data = ora_csv,
    selected_modules = selected_modules,
    title_text = "Functional enrichment of trait-associated WGCNA modules",
    filename = "selected_modules_ORA_heatmap.png",
    top_n_pathways = 30
  )
} else {
  save_empty_plot(
    "selected_modules_ORA_heatmap.png",
    "Functional enrichment of trait-associated WGCNA modules",
    "No selected module was available."
  )
}

# ------------------------------------------------------------
# 9. Analysis inventory and session information
# ------------------------------------------------------------
analysis_inventory <- tibble(
  item = c(
    "Analysis",
    "Tissue",
    "WGCNA input",
    "WGCNA feature-selection adaptation",
    "Number of non-grey modules",
    "Mapped WGCNA features",
    "Unique Entrez genes in WGCNA universe",
    "Pathway source",
    "GO collections",
    "Original Step 3 pathway-size range",
    "Step 3 transcriptome pathway-retention threshold",
    "Top-5,000 WGCNA pathway adaptation",
    "Minimum post-intersection pathway size",
    "ORA function",
    "P-value method",
    "Global adjustment",
    "Significance threshold",
    "ORA directionality",
    "Official fora2 source SHA-256",
    "fgsea version"
  ),
  value = c(
    "WGCNA module over-representation analysis",
    "Rat heart",
    normalizePath(
      modules_file,
      winslash = "/",
      mustWork = FALSE
    ),
    "Top 5,000 across-sample-variable transcripts from Step 4",
    length(TRNSCRPT_MOD_LIST),
    n_distinct(TRNSCRPT_MODULES$feature_ID),
    length(TRNSCRPT_UNIVERSE),
    "Step 3 MSigDB v7.5.1 filtered GO pathways",
    "GO:BP, GO:CC, GO:MF",
    "15–300 genes",
    "At least 85% of original pathway retained in Step 3",
    "No second 85% filter; intersect with mapped WGCNA universe",
    "15 genes",
    "Official MoTrPAC fora2()",
    "Scaled P value",
    "TRUE across modules within GO subcategory",
    "Scaled P < 0.05",
    "Non-directional",
    official_sources$sha256[
      official_sources$source_name == "fora2.R"
    ],
    as.character(packageVersion("fgsea"))
  )
)

write_csv(
  analysis_inventory,
  file.path(
    output_dir,
    "WGCNA_ORA_analysis_inventory.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

saveRDS(
  list(
    module_gene_rows = TRNSCRPT_MODULES,
    module_gene_lists = TRNSCRPT_MOD_LIST,
    universe = TRNSCRPT_UNIVERSE,
    filtered_pathways = TRNSCRPT_MSIGDB,
    module_ORA = TRNSCRPT_MODULE_ORA,
    significant_module_ORA = significant_ora,
    module_summary = module_ORA_summary,
    analysis_inventory = analysis_inventory
  ),
  file.path(
    output_dir,
    "rat_heart_GitHub_WGCNA_module_ORA_complete_objects.rds"
  ),
  compress = "xz"
)

message("")
message(
  "SUCCESS: GitHub-aligned WGCNA module ORA completed."
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
  "Modules analyzed: ",
  length(TRNSCRPT_MOD_LIST)
)
message(
  "WGCNA Entrez-gene universe: ",
  length(TRNSCRPT_UNIVERSE)
)
message(
  "GO pathways tested after WGCNA intersection: ",
  nrow(TRNSCRPT_MSIGDB)
)
message(
  "Significant module-pathway associations: ",
  nrow(significant_ora)
)
