########## This code builds Figure 2 (volcano plot, heatmap and paired violin plots of selected proteins) and exports it as a 177 x 230 mm, 600 dpi TIFF ##########
# Input : differential abundance table of 2.1 (R_objects/Neutrophil_proteome_limma.rds)
# Panel A: volcano plot of log2 FC (CD177pos / CD177neg) against the adjusted P value; exploratory candidates (|log2 FC| > 0.5,
#          adjusted P < 0.2) in red (more abundant in CD177pos) and blue (less abundant); dashed lines mark the exploratory
#          thresholds; the selected proteins are labelled
# Panel B: clustered heatmap of the selected significant and exploratory proteins (log2 abundance, row z-scores; Euclidean
#          distance, Ward's method); columns labelled by donor (D1-D8) and annotated by CD177 subset
# Panel C: paired violin plots of CD177, EPX and MPO (donors connected by lines; brackets give the adjusted P value)
# Selected proteins, figure size and colours are set in 1.2.
# Run with the working directory set to this repository (paths in 0.config.R).

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.2.Functions_proteomics.R")


proteome <- readRDS(file.path(objects_dir, "Neutrophil_proteome_limma.rds"))
protein_info <- proteome$protein_info
base_data <- proteome$base_data
FC_calculation <- proteome$FC_calculation

##### Panel A - volcano plot
volcano_df <- FC_calculation %>%
  mutate(sig = case_when(`log2 FC` > log2fc_cutoff & p_adj < padj_cutoff ~ "More abundant", `log2 FC` < -log2fc_cutoff & p_adj < padj_cutoff ~ "Less abundant", TRUE ~ "Other"),
         label = ifelse(Accession %in% names(selected_protein_labels) & p_adj < padj_cutoff & abs(`log2 FC`) > log2fc_cutoff, selected_protein_labels[Accession], NA_character_))

hemoglobin_labels <- volcano_df %>% filter(!is.na(label) & sig == "Less abundant") %>%
  mutate(nudge_x = ifelse(label == "Hemoglobin subunit beta", 0, -0.15), nudge_y = ifelse(label == "Hemoglobin subunit beta", 0.8, 0), hjust = ifelse(nudge_x < 0, 1, 0.5))

volcano_panel <- ggplot(volcano_df, aes(x = `log2 FC`, y = log10_p_adj)) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(name = NULL, values = c("More abundant" = "#D62728", "Less abundant" = "#1F77B4", "Other" = "grey80"),
                     breaks = c("More abundant", "Less abundant", "Other"),
                     labels = c("More abundant in CD177pos", "Less abundant in CD177pos", "Other proteins")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.3) +
  ### Labels at 7 pt, kept above the P value threshold and outside the fold-change thresholds (more abundant right, less abundant left)
  geom_text_repel(data = volcano_df %>% filter(!is.na(label) & sig == "More abundant" & label != "CD177"), aes(label = label), size = 7 / .pt, max.overlaps = Inf,
                  xlim = c(log2fc_cutoff + 0.2, NA), ylim = c(-log10(padj_cutoff) + 0.3, NA), box.padding = 0.4, point.padding = 0.1, force = 2,
                  segment.color = "grey55", segment.size = 0.3, min.segment.length = 0, seed = 1) +
  ### The hemoglobin points lie close together, so their labels are placed (beta above its point, alpha and delta left of theirs)
  geom_text_repel(data = hemoglobin_labels, aes(label = label, hjust = hjust), size = 7 / .pt, nudge_x = hemoglobin_labels$nudge_x, nudge_y = hemoglobin_labels$nudge_y,
                  direction = "y", force = 0, max.overlaps = Inf, point.padding = 0.1, segment.color = "grey55", segment.size = 0.3, min.segment.length = 0, seed = 1) +
  geom_text_repel(data = subset(volcano_df, label == "CD177"), aes(label = label), size = 7 / .pt, max.overlaps = Inf, nudge_x = -0.6, box.padding = 0.35,
                  point.padding = 0.1, segment.color = "grey55", segment.size = 0.3, min.segment.length = 0, seed = 1) +
  theme_classic(base_size = 8) +
  theme(legend.position = "right", legend.text = element_text(size = 8)) +
  labs(x = "log2 fold change (CD177pos / CD177neg)", y = "-log10 (BH-adjusted P value)")

volcano_panel

##### Panel B - heatmap of the selected proteins (rows and columns clustered)
heatmap_obj <- plot_protein_heatmap(heatmap_accessions, protein_info, base_data, labels = selected_protein_labels)

##### Panel C - paired violin plots of CD177, EPX and MPO
violin_data <- prepare_violin_data(violin_accessions, selected_protein_labels, protein_info, base_data, FC_calculation)
plot_df <- violin_data$plot_df
stat_df <- violin_data$stat_df

violin_panel <- ggplot(plot_df, aes(x = group, y = log2_abundance)) +
  geom_violin(aes(fill = group), alpha = 0.55, color = NA, trim = TRUE) +
  geom_line(aes(group = donor), color = "grey0", alpha = 0.6, linewidth = 0.20) +
  geom_point(aes(color = group), size = 1.8, alpha = 0.9) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colours) +
  scale_color_manual(values = group_colours) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8), legend.position = "none") +
  stat_pvalue_manual(stat_df, label = "p.label", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.01, bracket.size = 0.3, size = 3) +
  labs(x = NULL, y = "log2 abundance")

violin_panel

##### Multipanel TIFF
### ragg writes the resolution into the file (the default quartz tiff device on macOS records 72 dpi)
agg_tiff(file.path(figures_dir, "Figure_CD177_multipanel.tiff"), width = cd177_figure_size_mm[["width"]], height = cd177_figure_size_mm[["height"]], units = "mm",
         res = cd177_figure_dpi, compression = "lzw")
grid.newpage()
layout <- grid.layout(nrow = 5, ncol = 1, heights = unit(c(0.28, 0.02, 0.32, 0.02, 0.32), "npc"))
pushViewport(viewport(layout = layout))
print(volcano_panel, vp = viewport(layout.pos.row = 1))
pushViewport(viewport(layout.pos.row = 3))
grid.draw(heatmap_obj$gtable)
upViewport(1)
print(violin_panel, vp = viewport(layout.pos.row = 5))
upViewport(1)
dev.off()
