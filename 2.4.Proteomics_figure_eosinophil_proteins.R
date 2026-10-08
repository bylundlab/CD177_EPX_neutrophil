########## This code builds Figure 3G-I (clustered heatmap, stacked paired violin plots and donor-wise log2 FC heatmap) of the HPA eosinophil-enriched proteins detected in the neutrophil proteome and exports the panels as PDFs ##########
# Input : differential abundance table of 2.1 (R_objects/Neutrophil_proteome_limma.rds) and the HPA eosinophil-enriched
#         proteins of 2.2 (R_objects/HPA_eosinophil_enriched_proteins.rds)
# Panels: clustered heatmap (G) 83.6 x 78 mm, stacked violin plots (H) 83.6 x 78 mm and donor-wise log2 FC heatmap (I) 45.6 x 100 mm. They
#         are exported at their final size in the multipanel figure (sizes and the 8 pt font in 1.2), so they are placed at 100 %.
# Run with the working directory set to this repository (paths in 0.config.R).

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.2.Functions_proteomics.R")

##### Load the results of 2.1 and 2.2
proteome <- readRDS(file.path(objects_dir, "Neutrophil_proteome_limma.rds"))
protein_info <- proteome$protein_info
base_data <- proteome$base_data
FC_calculation <- proteome$FC_calculation
eos_strict_table_unique <- readRDS(file.path(objects_dir, "HPA_eosinophil_enriched_proteins.rds"))

##### Proteins of the panels
### Strict HPA eosinophil-enriched proteins, excluding EPX, ALOX15 and catalase; labelled with the HPA gene symbols
eos_plot_table <- eos_strict_table_unique %>% filter(!Accession %in% c("P11678", "P16050", "P04040"))
sel_eos <- eos_plot_table$Accession
eos_labels <- setNames(eos_plot_table$Gene_symbol, eos_plot_table$Accession)
eos_labels["P23610"] <- "F8A1"  # the HPA lists this accession under F8A2 and F8A3

##### Stacked paired violin plots
violin_data <- prepare_violin_data(sel_eos, eos_labels, protein_info, base_data, FC_calculation)
plot_df_eos <- violin_data$plot_df
stat_df_eos <- violin_data$stat_df

violin_eos <- ggplot(plot_df_eos, aes(x = group, y = log2_abundance)) +
  geom_violin(aes(fill = group), alpha = 0.55, color = NA, trim = TRUE) +
  geom_line(aes(group = donor), color = "grey0", alpha = 0.6, linewidth = 0.2) +
  geom_point(aes(color = group), size = 1.2, alpha = 0.9) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colours) +
  scale_color_manual(values = group_colours) +
  scale_x_discrete(labels = c("CD177neg" = "neg", "CD177pos" = "pos")) +
  scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.05, 0.25))) +
  stat_pvalue_manual(stat_df_eos, label = "p.label", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.01, bracket.size = 0.3,
                     size = font_size / .pt, vjust = -0.35) +
  labs(x = NULL, y = "log2 abundance") +
  theme_eos +
  theme(legend.position = "none", panel.spacing = unit(2, "mm"))

##### Clustered heatmap
eos_mat_scaled <- scaled_expression_matrix(sel_eos, protein_info, base_data)
eos_row_labels <- unname(eos_labels[rownames(eos_mat_scaled)])
column_annotation <- heatmap_column_annotation(colnames(eos_mat_scaled))
heatmap_cols <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)

heatmap_eos_panel <- pheatmap(eos_mat_scaled, annotation_col = column_annotation$annotation, annotation_colors = list(Group = group_colours), labels_row = eos_row_labels,
                              labels_col = column_annotation$donor_labels, color = heatmap_cols, clustering_distance_rows = "euclidean",
                              clustering_distance_cols = "euclidean", clustering_method = "ward.D2", show_colnames = TRUE, angle_col = 90, annotation_names_col = FALSE,
                              annotation_legend = FALSE, legend = FALSE, fontsize = font_size, fontsize_row = font_size, fontsize_col = font_size, treeheight_row = 10,
                              treeheight_col = 10, border_color = NA, silent = TRUE)

hm_range <- range(eos_mat_scaled, na.rm = TRUE)
heatmap_scale_source <- ggplot(data.frame(x = seq(hm_range[1], hm_range[2], length.out = 100), y = 1), aes(x = x, y = y, fill = x)) +
  geom_tile() +
  scale_fill_gradientn(colours = heatmap_cols, limits = hm_range, name = "Z-score",
                       guide = guide_colorbar(theme = theme(legend.key.width = unit(14, "mm"), legend.key.height = unit(1.5, "mm")))) +
  theme_eos_legend
group_legend_source <- ggplot(data.frame(Group = factor(c("CD177pos", "CD177neg"), levels = c("CD177pos", "CD177neg")), x = 1:2, y = 1), aes(x = x, y = y, fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = group_colours) +
  guides(fill = guide_legend(title = "Group", nrow = 1)) +
  theme_eos_legend +
  theme(legend.key.size = unit(2.5, "mm"))
heatmap_legends <- cowplot::plot_grid(cowplot::get_legend(heatmap_scale_source), cowplot::get_legend(group_legend_source), nrow = 1, rel_widths = c(0.37, 0.63))

fig_heatmap <- cowplot::plot_grid(heatmap_eos_panel$gtable, heatmap_legends, ncol = 1, rel_heights = c(0.88, 0.12)) +
  theme(plot.margin = margin(1, 1, 1, 1, "mm"), plot.background = element_rect(fill = "white", colour = NA))

##### Donor-wise log2 FC heatmap
fc_accessions <- intersect(names(donor_fc_protein_labels), protein_info$Accession)

fc_labels <- donor_fc_protein_labels
fc_labels[c("P13727", "P68871", "P02042", "P69905")] <- c("MBP", "HBB", "HBD", "HBA")

fc_plot_df <- bind_cols(protein_info, base_data) %>% filter(Accession %in% fc_accessions) %>% select(Accession, all_of(sample_columns)) %>%
  pivot_longer(cols = all_of(sample_columns), names_to = "sample", values_to = "abundance") %>%
  mutate(sample = as.integer(sample), donor = sample_info$donor[match(sample, sample_info$sample)], group = sample_info$group[match(sample, sample_info$sample)]) %>%
  select(Accession, donor, group, abundance) %>% pivot_wider(names_from = group, values_from = abundance) %>%
  mutate(log2FC_donor = log2(CD177pos / CD177neg), Protein = unname(fc_labels[Accession]))

fc_summary <- fc_plot_df %>% group_by(Accession, Protein) %>% summarise(median = median(log2FC_donor, na.rm = TRUE), .groups = "drop")
protein_order <- fc_summary %>% arrange(desc(median)) %>% pull(Protein)
fc_plot_df <- fc_plot_df %>% mutate(Protein = factor(Protein, levels = rev(protein_order)), Donor = factor(paste0("D", donor), levels = paste0("D", unique(sample_info$donor))))

fc_cols_original <- rev(brewer.pal(11, "RdYlBu"))
fc_cols <- c(fc_cols_original[1:4], "#F0F0F0", "#FFFFFF", "#F0F0F0", fc_cols_original[8:11])
fc_values <- scales::rescale(c(-5, -4, -3, -2, -0.5, 0, 0.5, 2, 3, 4, 5))

fc_donor_strip <- ggplot(fc_plot_df, aes(x = Donor, y = Protein, fill = log2FC_donor)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradientn(colours = fc_cols, values = fc_values, limits = c(-5, 5), oob = scales::squish, breaks = c(-5, 0, 5), name = expression(log[2]~FC),
                       guide = guide_colorbar(theme = theme(legend.key.width = unit(20, "mm"), legend.key.height = unit(1.5, "mm")))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), position = "right") +
  theme_eos +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank(),
        axis.text.y.right = element_text(color = "black", margin = margin(l = 2)), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.box.spacing = unit(1.5, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(1, 1, 1, 1, "mm"))

save_panel_pdf(fig_heatmap, file.path(figures_dir, "Panel_HPA_heatmap.pdf"), eos_heatmap_size_mm)
save_panel_pdf(violin_eos, file.path(figures_dir, "Panel_violin_eosinophil.pdf"), eos_violin_size_mm)
save_panel_pdf(fc_donor_strip, file.path(figures_dir, "Panel_donor_log2FC_heatmap.pdf"), eos_donor_fc_size_mm)
