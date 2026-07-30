# Goal ----
# Examine the code used for the sexual-dimorphism analysis in the WAT-SC
# Article: https://www.nature.com/articles/s42255-023-00959-9
# GitHub: https://github.com/MoTrPAC/MotrpacRatTraining6moWATData/
# Specifically, recreate the MSIGDB_PATHWAYS object used by FGSEA
# This object will be used in FGSEA


# Load libraries -----
library(tidyverse)
library(here)

# the paper uses msigdbr 7.5.1 (see Suppl Methods). 
# uncomment and use the following line to install if needed
# remotes::install_version("msigdbr", version = "7.5.1")
library(msigdbr)

# for pulling GO, `ontologyIndex` is in CRAN
library(ontologyIndex)


# Identify the parent collections for each GO subcollections ----
# the paper uses the following 3 subcollections
go_subcollections <- c("GO:BP", "GO:MF", "GO:CC")

go_collection_info <- msigdbr_collections() %>%
  filter(gs_subcat %in% go_subcollections) %>%
  select(gs_cat, gs_subcat) %>%
  distinct() %>%
  # retain the order in `go_subcollections`
  slice(match(go_subcollections, gs_subcat))

go_collection_info
# in msigdbr 7.5.1, GO:BP, GO:MF, and GO:CC are subcollections of C5



# Retrieve the rat pathway membership from MSigDB ----
# define empty vector
go_memberships_long <- vector(
  mode = "list",
  length = nrow(go_collection_info)
)

# define species and msigdbr version
rat_species_name <- "Rattus norvegicus"
msigdb_version <- packageVersion("msigdbr")          # 7.5.1 for replication


# run loop to fpll out genes within each collection & subcollection
for (subcollection_index in 1:nrow(go_collection_info)) {
  parent_collection <- go_collection_info$gs_cat[subcollection_index]
  go_subcollection <- go_collection_info$gs_subcat[subcollection_index]
  
  message("Reading ", go_subcollection, " from MSigDB ",
          msigdb_version, " for ", rat_species_name, ".")
  
  go_memberships_long[[subcollection_index]] <- msigdbr(
    species = rat_species_name,
    category = parent_collection,
    subcategory = go_subcollection
  )
}

class(go_memberships_long)                         # list
head(go_memberships_long)


# turn list into df
go_memberships_long_df <- go_memberships_long %>%
  bind_rows()

class(go_memberships_long_df)                      # data.frame
glimpse(go_memberships_long_df)



# Clean up the gene list for FGSEA ----
# main task: convert pathway-gene rows into one row per pathway
go_pathways <- go_memberships_long_df %>%
  # keep fields needed for FGSEA
  select(gs_subcat, gs_exact_source, gs_description, entrez_gene) %>%
  distinct() %>%
  # construct the Entrez-ID list-column.
  group_by(gs_subcat, gs_exact_source, gs_description) %>%
  summarise(entrez_gene = list(entrez_gene),
            .groups = "drop")

glimpse(go_pathways)



# Read in GO term names ----
# we'll be using archived GO release
go_obo_release_date <- "2021-12-15"
go_obo_url <- paste0("https://release.geneontology.org/",
                     go_obo_release_date,
                     "/ontology/go-basic.obo")

# download go info
go_obo_file_path <- here("data_msigdb", "msigdb_go_basic_2021-12-15.obo")

download.file(
  url = go_obo_url,
  destfile = go_obo_file_path,
  mode = "wb"
)

# load go_ontology
go_ontology <- ontologyIndex::get_OBO(
  file = go_obo_file_path,
  propagate_relationships = "is_a",
  extract_tags = "minimal"
)

glimpse(as.data.frame(go_ontology))


# clean up go_ontology
go_term_names <- as.data.frame(go_ontology) %>%
  filter(!obsolete, str_starts(id, "GO")) %>%
  mutate(gs_exact_source = id, go_description = name, .keep = "none")



# Update and capitalize GO pathway descriptions ----
msigdb_pathways <- go_pathways %>%
  left_join(go_term_names, by = "gs_exact_source") %>%
  mutate(
    gs_description = coalesce(go_description, gs_description),
    first_word = word(gs_description, 1),
    capitalize_description = first_word == str_to_lower(first_word),
    gs_description = if_else(capitalize_description,
                             str_replace(gs_description, "^.", str_to_upper),
                             gs_description)
  ) %>%
  select(-go_description, -first_word, -capitalize_description)



# Retain pathways with a certain number of genes (i.e, 15-300) ----
min_pathway_size <- 15L
max_pathway_size <- 300L

msigdb_pathways_selected <- msigdb_pathways %>%
  mutate(set_size = lengths(entrez_gene)) %>%
  filter(set_size >= min_pathway_size, set_size <= max_pathway_size) %>%
  select(gs_subcat, gs_exact_source, gs_description, entrez_gene, set_size)


# Check pathway counts ----
msigdb_pathways_selected %>%
  count(gs_subcat, name = "observed_pathways")

# matching the number in `MotrpacRatTraining6moWATData/data-raw/msigdb_pathways.R`
# GO:BP GO:CC GO:MF
#  3913   437   756



# Export msigdb_pathways_selected ----
# IMPORTANT! Don't save as a csv file bc of the list-column
saveRDS(msigdb_pathways_selected,
        file = here("data_msigdb", "msigdb_pathways.rds"))


# Clear environment ----
rm(list = ls())

