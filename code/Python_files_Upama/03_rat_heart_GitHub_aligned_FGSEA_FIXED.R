# ============================================================
# MoTrPAC Rat Heart Hackathon
# Step 3: GitHub-aligned transcriptomic FGSEA
#
# Directly sources the official MoTrPAC helper implementations:
#   - rank_genes.R
#   - fgsea2.R
#
# Reproduces the official transcriptomic FGSEA choices:
#   - MSigDB v7.5.1 GO:BP / GO:MF / GO:CC
#   - original pathway size 15–300
#   - at least 85% pathway membership retained after intersecting
#     with detected transcriptomic Entrez genes
#   - at least 15 detected genes after intersection
#   - signed -log10(P.Value) ranking, averaged by Entrez gene
#   - seed = 0
#   - nPermSimple = 10,000
#   - BH adjustment globally by GO subcategory within each
#     related contrast family
#
# Tissue adaptation:
#   WAT transcriptomics -> HEART transcriptomics
# Patch: corrected nested-list indexing in ranking inventory
# ============================================================

options(timeout = 1e5)
options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(0)

message("R version: ", R.version.string)
message("Working directory: ", getwd())

# ------------------------------------------------------------
# 1. Install packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "remotes", "dplyr", "tidyr", "tibble", "purrr", "readr",
  "ggplot2", "stringr", "forcats", "scales", "data.table", "digest"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
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

# The published workflow used msigdbr 7.5.1.
required_msigdbr_version <- "7.5.1"

install_msigdbr_751 <- function() {
  remotes::install_version(
    "msigdbr",
    version = required_msigdbr_version,
    repos = "https://cloud.r-project.org",
    upgrade = "never",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install_msigdbr_751()
} else if (
  as.character(packageVersion("msigdbr")) != required_msigdbr_version
) {
  message(
    "Replacing msigdbr ",
    as.character(packageVersion("msigdbr")),
    " with version 7.5.1."
  )
  install_msigdbr_751()
}

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
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
})

message("fgsea version: ", as.character(packageVersion("fgsea")))
message("msigdbr version: ", as.character(packageVersion("msigdbr")))

# ------------------------------------------------------------
# 2. Input/output paths
# ------------------------------------------------------------
de_file <- file.path(
  "results",
  "02_rat_heart_differential_expression",
  "all_differential_expression_results.csv"
)

if (!file.exists(de_file)) {
  stop(
    "Step 2 differential-expression result was not found:\n",
    normalizePath(de_file, winslash = "/", mustWork = FALSE)
  )
}

de_results <- read_csv(
  de_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "feature",
  "gene_symbol",
  "entrez_gene",
  "logFC",
  "P.Value",
  "contrast",
  "contrast_group",
  "plain_language",
  "positive_logFC_means"
)

missing_columns <- setdiff(
  required_columns,
  colnames(de_results)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing Step 2 columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

output_dir <- file.path(
  "results",
  "03_rat_heart_GitHub_FGSEA"
)

plot_dir <- file.path(output_dir, "plots")
rank_dir <- file.path(output_dir, "ranked_gene_lists")
source_dir <- file.path(output_dir, "official_source_snapshots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Download and source official helper functions unchanged
# ------------------------------------------------------------
official_sources <- tibble(
  source_name = c(
    "rank_genes.R",
    "fgsea2.R",
    "msigdb_pathways.R"
  ),
  source_url = c(
    paste0(
      "https://raw.githubusercontent.com/",
      "MoTrPAC/MotrpacRatTraining6moWAT/",
      "master/R/rank_genes.R"
    ),
    paste0(
      "https://raw.githubusercontent.com/",
      "MoTrPAC/MotrpacRatTraining6moWAT/",
      "master/R/fgsea2.R"
    ),
    paste0(
      "https://raw.githubusercontent.com/",
      "MoTrPAC/MotrpacRatTraining6moWATData/",
      "master/data-raw/msigdb_pathways.R"
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

# Source the exact official ranking and FGSEA wrappers.
source(
  file.path(source_dir, "rank_genes.R"),
  local = .GlobalEnv
)

source(
  file.path(source_dir, "fgsea2.R"),
  local = .GlobalEnv
)

if (!exists("rank_genes", mode = "function")) {
  stop("Official rank_genes() was not loaded.")
}

if (!exists("fgsea2", mode = "function")) {
  stop("Official fgsea2() was not loaded.")
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
# 4. Reformat one-to-many transcript-to-gene mappings
# ------------------------------------------------------------
# Mirrors the manuscript vignette:
# separate transcript rows containing multiple Entrez IDs / symbols.

de_separated <- de_results |>
  mutate(
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
  )

contrast_info <- de_results |>
  distinct(
    contrast_group,
    contrast,
    plain_language,
    positive_logFC_means
  )

family_order <- c(
  "trained_vs_SED",
  "MvF_SED",
  "MvF_exercise_response"
)

missing_families <- setdiff(
  family_order,
  unique(de_separated$contrast_group)
)

if (length(missing_families) > 0) {
  stop(
    "Missing expected contrast families: ",
    paste(missing_families, collapse = ", ")
  )
}

de_families <- lapply(
  family_order,
  function(family_name) {
    de_separated |>
      filter(contrast_group == family_name) |>
      as.data.frame()
  }
)

names(de_families) <- family_order

# Entrez-to-symbol lookup for leading-edge reporting.
# The official vignette constructs this from the sedentary sex contrast.
entrez_to_symbol_table <- de_families[["MvF_SED"]] |>
  select(entrez_gene, gene_symbol) |>
  filter(
    !is.na(entrez_gene),
    entrez_gene != "",
    entrez_gene != "NA"
  ) |>
  group_by(entrez_gene) |>
  summarise(
    gene_symbol = first(
      gene_symbol[
        !is.na(gene_symbol) &
          gene_symbol != "" &
          gene_symbol != "NA"
      ],
      default = NA_character_
    ),
    .groups = "drop"
  )

entrez_to_symbol <- setNames(
  entrez_to_symbol_table$gene_symbol,
  entrez_to_symbol_table$entrez_gene
)

detected_entrez <- names(entrez_to_symbol)

message(
  "Detected unique transcriptomic Entrez genes: ",
  length(detected_entrez)
)

# ------------------------------------------------------------
# 5. Build MSigDB v7.5.1 GO pathways as in GitHub
# ------------------------------------------------------------
load_go_collection <- function(subcategory_code, subcategory_label) {
  msigdbr(
    species = "Rattus norvegicus",
    category = "C5",
    subcategory = subcategory_code
  ) |>
    transmute(
      gs_subcat = subcategory_label,
      gs_exact_source = as.character(gs_exact_source),
      gs_description = as.character(gs_description),
      entrez_gene = as.character(entrez_gene)
    ) |>
    filter(
      !is.na(gs_exact_source),
      !is.na(entrez_gene),
      entrez_gene != ""
    )
}

msigdb_long <- bind_rows(
  load_go_collection("BP", "GO:BP"),
  load_go_collection("MF", "GO:MF"),
  load_go_collection("CC", "GO:CC")
) |>
  distinct()

# Equivalent structure to official msigdbr2(..., genes="entrez_gene").
MSIGDB_PATHWAYS <- msigdb_long |>
  group_by(
    gs_subcat,
    gs_exact_source,
    gs_description
  ) |>
  summarise(
    entrez_gene = list(unique(entrez_gene)),
    .groups = "drop"
  ) |>
  mutate(
    set_size = lengths(entrez_gene)
  ) |>
  # Exact prefilter from data-raw/msigdb_pathways.R.
  filter(
    set_size >= 15,
    set_size <= 300
  )

# Official transcriptomic coverage filtering:
# retain >=85% of original pathway genes and >=15 detected genes.
TRNSCRPT_MSIGDB <- MSIGDB_PATHWAYS |>
  mutate(
    entrez_gene = map(
      entrez_gene,
      intersect,
      y = detected_entrez
    ),
    set_size_post = lengths(entrez_gene),
    ratio = set_size_post / set_size
  ) |>
  filter(
    ratio >= 0.85,
    set_size_post >= 15
  )

if (nrow(TRNSCRPT_MSIGDB) == 0) {
  stop(
    "No pathways remained after the official 85% coverage and ",
    "15-gene filters."
  )
}

pathway_filter_summary <- bind_rows(
  MSIGDB_PATHWAYS |>
    count(gs_subcat, name = "pathways_before_coverage") |>
    mutate(stage = "Original size 15–300"),
  TRNSCRPT_MSIGDB |>
    count(gs_subcat, name = "pathways_after_coverage") |>
    mutate(stage = "Coverage >=85% and post-size >=15")
)

pathway_filter_wide <- full_join(
  MSIGDB_PATHWAYS |>
    count(gs_subcat, name = "pathways_before_coverage"),
  TRNSCRPT_MSIGDB |>
    count(gs_subcat, name = "pathways_after_coverage"),
  by = "gs_subcat"
) |>
  mutate(
    pathways_removed = (
      pathways_before_coverage -
        pathways_after_coverage
    )
  )

write_csv(
  pathway_filter_wide,
  file.path(output_dir, "pathway_filter_summary.csv")
)

pathway_membership_table <- TRNSCRPT_MSIGDB |>
  mutate(
    entrez_gene = map_chr(
      entrez_gene,
      ~ paste(.x, collapse = ";")
    )
  )

write_csv(
  pathway_membership_table,
  file.path(
    output_dir,
    "filtered_MSigDB_GO_pathways_85pct.csv"
  )
)

# ------------------------------------------------------------
# 6. Rank genes using the official rank_genes()
# ------------------------------------------------------------
ranked_stats_by_family <- lapply(
  de_families,
  function(family_data) {
    rank_genes(
      family_data,
      genes = "entrez_gene"
    )
  }
)

for (family_name in names(ranked_stats_by_family)) {
  family_stats <- ranked_stats_by_family[[family_name]]

  for (contrast_name in names(family_stats)) {
    rank_table <- tibble(
      entrez_gene = names(family_stats[[contrast_name]]),
      rank_metric = as.numeric(
        family_stats[[contrast_name]]
      )
    ) |>
      left_join(
        entrez_to_symbol_table,
        by = "entrez_gene"
      )

    safe_name <- contrast_name |>
      str_replace_all("[^A-Za-z0-9]+", "_") |>
      str_replace_all("^_|_$", "")

    write_csv(
      rank_table,
      file.path(
        rank_dir,
        paste0(safe_name, "_ranked_genes.csv")
      )
    )
  }
}

ranking_inventory <- bind_rows(
  lapply(
    names(ranked_stats_by_family),
    function(family_name) {
      bind_rows(
        lapply(
          names(ranked_stats_by_family[[family_name]]),
          function(contrast_name) {
            tibble(
              contrast_group = family_name,
              contrast = contrast_name,
              ranked_entrez_genes = length(
                ranked_stats_by_family[[family_name]][[contrast_name]]
              )
            )
          }
        )
      )
    }
  )
) |>
  left_join(
    contrast_info,
    by = c("contrast_group", "contrast")
  )

write_csv(
  ranking_inventory,
  file.path(output_dir, "ranking_gene_counts.csv")
)

message(
  "Gene-ranking inventory completed for ",
  nrow(ranking_inventory),
  " contrasts."
)

# ------------------------------------------------------------
# 7. Run official fgsea2() separately for each contrast family
# ------------------------------------------------------------
# This is important: adjust.globally=TRUE pools p-values across
# related contrasts within each family, separately by GO subcategory.

fgsea_by_family <- lapply(
  family_order,
  function(family_name) {
    message("Running FGSEA family: ", family_name)

    result <- fgsea2(
      pathways = as.data.frame(TRNSCRPT_MSIGDB),
      stats = ranked_stats_by_family[[family_name]],
      seed = 0,
      nPermSimple = 10000,
      adjust.globally = TRUE,
      nproc = 1
    )

    result[, contrast_group := family_name]

    result[, leadingEdge_genes := lapply(
      leadingEdge,
      function(x) {
        symbols <- unname(
          entrez_to_symbol[
            as.character(x)
          ]
        )

        unique(
          na.omit(symbols)
        )
      }
    )]

    result
  }
)

names(fgsea_by_family) <- family_order

fgsea_native <- rbindlist(
  fgsea_by_family,
  use.names = TRUE,
  fill = TRUE
)

fgsea_native <- merge(
  fgsea_native,
  as.data.table(contrast_info),
  by = c("contrast_group", "contrast"),
  all.x = TRUE,
  sort = FALSE
)

fgsea_native[, enrichment_direction := fifelse(
  NES > 0,
  "Positive",
  fifelse(NES < 0, "Negative", "Zero")
)]

fgsea_native[, significant_FDR05 := (
  !is.na(padj) &
    padj < 0.05
)]

# Flatten list columns for CSV.
fgsea_csv <- copy(fgsea_native)

fgsea_csv[, leadingEdge_entrez := vapply(
  leadingEdge,
  function(x) paste(x, collapse = ";"),
  character(1)
)]

fgsea_csv[, leadingEdge_symbols := vapply(
  leadingEdge_genes,
  function(x) paste(x, collapse = ";"),
  character(1)
)]

fgsea_csv[, c("leadingEdge", "leadingEdge_genes") := NULL]

# Add coverage information to each result.
pathway_coverage_metadata <- TRNSCRPT_MSIGDB |>
  select(
    pathway = gs_exact_source,
    set_size,
    set_size_post,
    ratio
  ) |>
  distinct()

fgsea_csv <- as_tibble(fgsea_csv) |>
  left_join(
    pathway_coverage_metadata,
    by = "pathway"
  )

# ------------------------------------------------------------
# 8. Summary tables
# ------------------------------------------------------------
pathway_counts <- fgsea_csv |>
  group_by(
    contrast_group,
    contrast,
    plain_language,
    gs_subcat
  ) |>
  summarise(
    tested_pathways = n(),
    FDR05_total = sum(
      significant_FDR05,
      na.rm = TRUE
    ),
    FDR05_positive_NES = sum(
      significant_FDR05 &
        NES > 0,
      na.rm = TRUE
    ),
    FDR05_negative_NES = sum(
      significant_FDR05 &
        NES < 0,
      na.rm = TRUE
    ),
    minimum_p_value = min(
      pval,
      na.rm = TRUE
    ),
    minimum_global_FDR = min(
      padj,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

significant_pathways <- fgsea_csv |>
  filter(significant_FDR05) |>
  arrange(
    contrast_group,
    contrast,
    padj,
    desc(abs(NES))
  )

top20_pathways <- fgsea_csv |>
  group_by(
    contrast_group,
    contrast,
    gs_subcat
  ) |>
  arrange(
    padj,
    pval,
    desc(abs(NES)),
    .by_group = TRUE
  ) |>
  slice_head(n = 20) |>
  ungroup()

# ------------------------------------------------------------
# 9. Plot helpers
# ------------------------------------------------------------
pretty_description <- function(description, pathway) {
  output <- description

  missing_description <- (
    is.na(output) |
      output == ""
  )

  output[missing_description] <- pathway[
    missing_description
  ]

  str_to_sentence(output)
}

make_dotplot <- function(
  contrast_name,
  filename
) {
  plot_data <- fgsea_csv |>
    filter(
      contrast == contrast_name,
      gs_subcat == "GO:BP"
    )

  significant_data <- plot_data |>
    filter(padj < 0.05)

  if (nrow(significant_data) > 0) {
    shown <- bind_rows(
      significant_data |>
        filter(NES > 0) |>
        arrange(padj, desc(NES)) |>
        slice_head(n = 10),
      significant_data |>
        filter(NES < 0) |>
        arrange(padj, NES) |>
        slice_head(n = 10)
    )

    subtitle_text <- paste(
      "GitHub-aligned FGSEA;",
      "global family-level GO:BP FDR < 0.05"
    )
  } else {
    shown <- plot_data |>
      arrange(padj, pval) |>
      slice_head(n = 15)

    subtitle_text <- paste(
      "No GO:BP pathway reached global family-level FDR < 0.05;",
      "lowest-FDR pathways shown"
    )
  }

  title_text <- contrast_info |>
    filter(contrast == contrast_name) |>
    pull(plain_language) |>
    first()

  shown <- shown |>
    mutate(
      pathway_label = pretty_description(
        gs_description,
        pathway
      ),
      pathway_label = fct_reorder(
        pathway_label,
        NES
      ),
      FDR_score = -log10(
        pmax(
          padj,
          .Machine$double.xmin
        )
      )
    )

  plot_object <- ggplot(
    shown,
    aes(
      x = NES,
      y = pathway_label,
      size = FDR_score,
      fill = NES
    )
  ) +
    geom_point(
      shape = 21,
      colour = "black",
      alpha = 0.9
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0
    ) +
    scale_size_continuous(
      range = c(2.5, 8)
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Normalized enrichment score (NES)",
      y = NULL,
      fill = "NES",
      size = "-log10(FDR)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(plot_dir, filename),
    plot_object,
    width = 10,
    height = 7,
    dpi = 300
  )
}

priority_contrasts <- c(
  "M_SED - F_SED",
  "F_8W - F_SED",
  "M_8W - M_SED",
  "(M_8W - M_SED) - (F_8W - F_SED)"
)

priority_files <- c(
  "baseline_sex_top_GO_BP.png",
  "female_8w_top_GO_BP.png",
  "male_8w_top_GO_BP.png",
  "interaction_8w_top_GO_BP.png"
)

for (i in seq_along(priority_contrasts)) {
  make_dotplot(
    priority_contrasts[i],
    priority_files[i]
  )
}

# Significant pathway count plot.
count_plot_data <- pathway_counts |>
  group_by(
    contrast_group,
    contrast,
    plain_language
  ) |>
  summarise(
    positive = sum(FDR05_positive_NES),
    negative = sum(FDR05_negative_NES),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(positive, negative),
    names_to = "direction",
    values_to = "n"
  ) |>
  mutate(
    signed_n = ifelse(
      direction == "negative",
      -n,
      n
    ),
    direction = recode(
      direction,
      positive = "Positive NES",
      negative = "Negative NES"
    ),
    plain_language = factor(
      plain_language,
      levels = rev(
        contrast_info$plain_language
      )
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
    title = "Significant GO pathways across rat-heart contrasts",
    subtitle = paste(
      "GitHub-aligned FGSEA;",
      "BH FDR pooled within each contrast family and GO category"
    ),
    x = "Number of significant pathways",
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
    "significant_GO_pathway_counts.png"
  ),
  count_plot,
  width = 12,
  height = 8,
  dpi = 300
)

make_heatmap <- function(
  data,
  title_text,
  filename,
  x_levels,
  top_n
) {
  candidates <- data |>
    group_by(pathway) |>
    summarise(
      minimum_FDR = min(
        padj,
        na.rm = TRUE
      ),
      maximum_absolute_NES = max(
        abs(NES),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(
      minimum_FDR,
      desc(maximum_absolute_NES)
    ) |>
    slice_head(n = top_n)

  heatmap_data <- data |>
    filter(
      pathway %in% candidates$pathway
    ) |>
    mutate(
      pathway_label = pretty_description(
        gs_description,
        pathway
      ),
      contrast_label = factor(
        contrast,
        levels = x_levels
      ),
      significance_mark = ifelse(
        padj < 0.05,
        "*",
        ""
      )
    )

  pathway_order <- heatmap_data |>
    group_by(pathway_label) |>
    summarise(
      minimum_FDR = min(
        padj,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(minimum_FDR) |>
    pull(pathway_label)

  heatmap_data <- heatmap_data |>
    mutate(
      pathway_label = factor(
        pathway_label,
        levels = rev(pathway_order)
      )
    )

  heatmap_plot <- ggplot(
    heatmap_data,
    aes(
      x = contrast_label,
      y = pathway_label,
      fill = NES
    )
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = significance_mark),
      size = 3.5
    ) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0
    ) +
    labs(
      title = title_text,
      subtitle = paste(
        "* indicates global family-level",
        "FGSEA FDR < 0.05"
      ),
      x = NULL,
      y = NULL,
      fill = "NES"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 40,
        hjust = 1
      ),
      panel.grid = element_blank()
    )

  ggsave(
    file.path(plot_dir, filename),
    heatmap_plot,
    width = 12,
    height = 9,
    dpi = 300
  )
}

priority_data <- fgsea_csv |>
  filter(
    contrast %in% priority_contrasts,
    gs_subcat == "GO:BP"
  )

make_heatmap(
  priority_data,
  "Top GO biological processes in priority contrasts",
  "top_pathways_priority_contrasts_heatmap.png",
  priority_contrasts,
  25
)

training_contrasts <- contrast_info |>
  filter(
    contrast_group == "trained_vs_SED"
  ) |>
  pull(contrast)

training_data <- fgsea_csv |>
  filter(
    contrast %in% training_contrasts,
    gs_subcat == "GO:BP"
  )

make_heatmap(
  training_data,
  "Top training-responsive GO biological processes",
  "training_timecourse_top_pathways_heatmap.png",
  training_contrasts,
  20
)

cardiac_pattern <- paste(
  c(
    "CARDIAC",
    "HEART",
    "MUSCLE CONTRACTION",
    "SARCOMERE",
    "MYOFIBRIL",
    "CALCIUM",
    "MITOCHON",
    "OXIDATIVE PHOSPHORYLATION",
    "RESPIRATORY",
    "ATP",
    "FATTY ACID",
    "ANGIOGEN",
    "EXTRACELLULAR MATRIX",
    "COLLAGEN",
    "AUTOPHAG",
    "INFLAM"
  ),
  collapse = "|"
)

cardiac_results <- fgsea_csv |>
  mutate(
    cardiac_search_text = str_to_upper(
      paste(
        pathway,
        gs_description
      )
    )
  ) |>
  filter(
    gs_subcat == "GO:BP",
    str_detect(
      cardiac_search_text,
      cardiac_pattern
    )
  ) |>
  select(-cardiac_search_text) |>
  arrange(
    contrast_group,
    contrast,
    padj,
    desc(abs(NES))
  )

cardiac_priority <- cardiac_results |>
  filter(
    contrast %in% priority_contrasts
  )

if (nrow(cardiac_priority) > 0) {
  make_heatmap(
    cardiac_priority,
    "Cardiac and metabolic GO pathway signatures",
    "cardiac_related_pathways_heatmap.png",
    priority_contrasts,
    25
  )
}

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------
analysis_inventory <- tibble(
  item = c(
    "Differential-expression contrasts",
    "Related contrast families",
    "MSigDB release",
    "Gene-set collections",
    "Original pathway-size requirement",
    "Detected-membership requirement",
    "Post-intersection minimum size",
    "Ranking metric",
    "Transcript-to-gene aggregation",
    "FGSEA helper",
    "Random seed",
    "nPermSimple",
    "Parallel processes",
    "FDR adjustment",
    "fgsea package version",
    "msigdbr package version"
  ),
  value = c(
    length(unique(de_results$contrast)),
    length(family_order),
    "7.5.1",
    "GO:BP / GO:MF / GO:CC",
    "15–300 genes",
    "At least 85% of original pathway genes retained",
    "At least 15 detected genes",
    "-log10(P.Value) * sign(logFC)",
    "Arithmetic mean by Entrez gene using official rank_genes()",
    "Official MoTrPAC fgsea2() sourced unchanged",
    0,
    10000,
    1,
    paste(
      "BH globally by GO subcategory",
      "within each related contrast family"
    ),
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("msigdbr"))
  )
)

write_csv(
  analysis_inventory,
  file.path(
    output_dir,
    "GitHub_FGSEA_analysis_inventory.csv"
  )
)

write_csv(
  fgsea_csv,
  file.path(
    output_dir,
    "all_GitHub_aligned_FGSEA_GO_results.csv"
  )
)

write_csv(
  significant_pathways,
  file.path(
    output_dir,
    "significant_GitHub_aligned_FGSEA_GO_pathways.csv"
  )
)

write_csv(
  pathway_counts,
  file.path(
    output_dir,
    "significant_GO_pathway_counts.csv"
  )
)

write_csv(
  top20_pathways,
  file.path(
    output_dir,
    "top20_GO_pathways_each_contrast.csv"
  )
)

write_csv(
  cardiac_results,
  file.path(
    output_dir,
    "cardiac_related_GO_pathways.csv"
  )
)

saveRDS(
  list(
    official_sources = official_sources,
    MSIGDB_PATHWAYS = MSIGDB_PATHWAYS,
    TRNSCRPT_MSIGDB = TRNSCRPT_MSIGDB,
    de_families = de_families,
    ranked_stats_by_family = ranked_stats_by_family,
    fgsea_by_family = fgsea_by_family,
    fgsea_native = fgsea_native,
    contrast_info = contrast_info
  ),
  file.path(
    output_dir,
    "rat_heart_GitHub_aligned_FGSEA_objects.rds"
  ),
  compress = "xz"
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\n================ GITHUB-ALIGNED FGSEA INVENTORY ================\n")
print(analysis_inventory)

cat("\n================ PATHWAY FILTER SUMMARY =========================\n")
print(pathway_filter_wide)

cat("\n================ SIGNIFICANT PATHWAY COUNTS =====================\n")
print(
  pathway_counts |>
    select(
      contrast_group,
      contrast,
      gs_subcat,
      FDR05_total,
      FDR05_positive_NES,
      FDR05_negative_NES,
      minimum_global_FDR
    )
)

cat("\nInterpretation:\n")
cat(
  "Positive NES follows positive logFC; negative NES follows negative logFC.\n"
)
cat(
  "M_SED - F_SED: positive = male-enriched; negative = female-enriched.\n"
)
cat(
  "Trained - SED: positive = training-upregulated; negative = training-downregulated.\n"
)
cat(
  "Interaction: positive = larger male response; negative = larger female response.\n"
)

cat("\nSUCCESS: GitHub-aligned FGSEA outputs saved to:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)
