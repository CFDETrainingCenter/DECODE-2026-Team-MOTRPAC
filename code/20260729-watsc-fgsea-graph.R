# Goal -----



# Load data ----
rna_fgsea <- readRDS(here("data_processed", "watsc_transcriptomics_fgsea.rds"))
protein_fgsea <- readRDS(here("data_processed", "watsc_proteomics_fgsea.rds"))


# Join pathways tested in both omics
# use inner_join -- it would require NES in both the transcriptome and proteome
shared_pathway_results <- inner_join(
  rna_fgsea %>% select(pathway, gs_subcat, gs_description, rna_padj = padj, rna_NES = NES),
  protein_fgsea %>% select(pathway, gs_subcat, gs_description, protein_padj = padj, protein_NES = NES),
  by = c("pathway", "gs_subcat", "gs_description")
) %>%
#  filter(rna_padj < 0.05 | protein_padj < 0.05) %>%
  mutate(significant_in_both = (rna_padj < 0.05 & protein_padj < 0.05),
         quadrant = case_when(
           protein_NES > 0 & rna_NES > 0 ~ "Q1",
           protein_NES < 0 & rna_NES > 0 ~ "Q2",
           protein_NES < 0 & rna_NES < 0 ~ "Q3",
           protein_NES > 0 & rna_NES < 0 ~ "Q4",
           TRUE ~ NA_character_
         ))

# count the number of significant genes in each subcategory
shared_pathway_results %>%
  filter(significant_in_both == TRUE) %>%
  count(gs_subcat, significant_in_both)



# Figure 2B: transcriptomics vs proteomic NES ----
ggplot(data = shared_pathway_results,
       mapping = aes(x = protein_NES, y = rna_NES, color = significant_in_both)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(size = 1, shape = 16, alpha = 0.8) +
  facet_grid(rows = vars(gs_subcat), space = "free_y") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"), breaks = TRUE,
                     name = NULL, labels = "Significantly enriched in both -omes") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        strip.background = element_blank()) +
  labs(title = "Scatter-plots of normalized enrichment scores", 
       subtitle = "WAT-SC of SED male vs female",
       x = "NES: PROT", y = "NES: TRNSCRPT")

ggsave(here("plots", "watsc_rna_protein_NES_MvF_SED_scatterplot.png"),
       height = 7, width = 4, dpi = 300, units = "in")





