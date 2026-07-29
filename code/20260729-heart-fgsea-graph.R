# Goal -----



# Load data ----
rna_fgsea <- readRDS(here("data_processed", "heart_transcriptomics_fgsea.rds"))
protein_fgsea <- read_csv(here("data_processed", "heart_prot_fgsea.csv"))


# Join pathways tested in both omics ----
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
       subtitle = "HEART of SED male vs female",
       x = "NES: PROT", y = "NES: TRNSCRPT")

ggsave(here("plots", "heart_rna_protein_NES_MvF_SED_scatterplot.png"),
       height = 7, width = 4, dpi = 300, units = "in")



# Select the terms displayed in Figure 2C ----
# first rank signficant pathways by the mean absolute NES across the two omics
# within each GO category, retain every non-Q1 pathway and then fill the remaining positions, 
# up to approximately 10, with Q1 terms.

# Jade's quick note: I think there is some issue with this code, will come back to inspect

ranked_significant_pathways <- shared_pathway_results %>%
  filter(significant_in_both) %>%
  mutate(mean_absolute_NES = (abs(protein_NES) + abs(rna_NES)) / 2) %>%
  group_by(gs_subcat) %>%
  arrange(desc(mean_absolute_NES), .by_group = TRUE) %>%
  mutate(number_outside_Q1 = sum(quadrant != "Q1"),
         number_of_Q1_terms_to_keep = 10L - number_outside_Q1) %>%
  group_by(gs_subcat, quadrant) %>%
  mutate(rank_within_quadrant = row_number()) %>%
  ungroup()

selected_figure_2c_pathways <- ranked_significant_pathways %>%
  filter((quadrant != "Q1" & rank_within_quadrant <= number_outside_Q1) |
           (quadrant == "Q1" & rank_within_quadrant <= number_of_Q1_terms_to_keep)) %>%
  pull(pathway)



# Figure 2C: NES values for selected pathways -----
shared_pathway_results_selected <- shared_pathway_results %>%
  filter(pathway %in% selected_figure_2c_pathways) %>%
  select(pathway, gs_subcat, gs_description, rna_NES, protein_NES) %>%
  pivot_longer(cols = c(rna_NES, protein_NES), names_to = "omic", values_to = "NES") %>%
  mutate(omic = recode(omic, "rna_NES" = "TRNSCRPT", "protein_NES" = "PROT")) %>%
  arrange(NES) %>%
  mutate(gs_description = factor(gs_description, levels = unique(gs_description)))


ggplot(data = shared_pathway_results_selected,
       mapping = aes(x = NES, y = gs_description, shape = omic, color = NES)) +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 3.6) +
  facet_grid(rows = vars(gs_subcat), scales = "free_y",
             space = "free_y", switch = "y") +
  scale_y_discrete(name = NULL, position = "right") +
  scale_color_gradientn(colors = c("#ff6eff", "white", "white", "#5555ff"),
                        values = scales::rescale(c(-2.5, -1, 1, 2.9))) +
  scale_shape_manual(name = NULL,
                     values = c("PROT" = 16, "TRNSCRPT" = 17)) +
  guides(color = guide_colorbar(barwidth = unit(0.6, "in"), barheight = unit(0.07, "in"), 
                                frame.colour = "black", ticks.colour = "black", order = 1),
         shape = guide_legend(title = NULL, order = 2, keyheight = unit(7, "pt"),
                              keywidth = unit(7, "pt"))) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, color = "black"),
        panel.spacing.y = unit(3, "pt"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.3, color = "grey92"),
        plot.margin = margin(t = 20, unit = "pt"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  labs(title = "Selected pathways that are significant enriched in both -omes",
       subtitle = "WAT-SC of SED male vs female, transcriptomics and proteomics")


ggsave(here("plots", "heart_rna_protein_pathway_NES_MvF_SED_scatterplot.png"),
       height = 7, width = 7, dpi = 300, units = "in")
