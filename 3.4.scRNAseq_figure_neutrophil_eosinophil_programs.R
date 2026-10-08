########## This code computes the EoP and NeP program scores and exports the panels of Figure 5 (UMAPs, score violins, FeaturePlots and DotPlot) for GSE180298, GSE175879 and E-MTAB-11188 ##########
# Cells are grouped by the de novo clusters of each dataset, shown by cluster number only
# GSE180298 (3.1, resolution 0.45), GSE175879 (3.2, resolution 0.15) and the steady-state BM/PB neutrophils of
# E-MTAB-11188 (3.3.5, resolution 0.15).
# Run with the working directory set to this repository (paths in 0.config.R).

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.4.Functions_scRNAseq_gene_sets_and_figures.R")

##### Load final objects, keeping only the normalized data and the UMAP used for plotting
G180 <- readRDS(file.path(objects_dir, "GSE180298_young_final.rds"))
G180 <- DietSeurat(G180, layers = "data", dimreducs = "umap")
BM <- readRDS(file.path(objects_dir, "GSE175879_NCP_Harmony.rds"))
BM <- DietSeurat(BM, layers = "data", dimreducs = "umap.harmony")
EMTAB <- readRDS(file.path(objects_dir, "EMTAB11188_steady_state.rds"))
EMTAB <- DietSeurat(EMTAB, layers = "data", dimreducs = "umap")

##### Program scores
G180 <- AddModuleScore(G180, features = list(intersect(Eos_core_genes, rownames(G180))), name = "Eos_core", assay = "RNA")
G180 <- AddModuleScore(G180, features = list(intersect(NeP_genes, rownames(G180))), name = "NeP", assay = "RNA")
BM <- AddModuleScore(BM, features = list(intersect(Eos_core_genes, rownames(BM))), name = "Eos_core", assay = "RNA")
BM <- AddModuleScore(BM, features = list(intersect(NeP_genes, rownames(BM))), name = "NeP", assay = "RNA")
EMTAB <- AddModuleScore(EMTAB, features = list(intersect(Eos_core_genes, rownames(EMTAB))), name = "Eos_core", assay = "RNA")
EMTAB <- AddModuleScore(EMTAB, features = list(intersect(NeP_genes, rownames(EMTAB))), name = "NeP", assay = "RNA")

cluster_label_angle <- function(clusters) if (length(clusters) > 10) 90 else 0

##### Main figure panels
layout_objects <- list(GSE180298 = G180, GSE175879_NCP = BM, EMTAB11188 = EMTAB)
layout_reductions <- c(GSE180298 = "umap", GSE175879_NCP = "umap.harmony", EMTAB11188 = "umap")
save_layout_panel <- function(plot, file, width_mm, height_mm = layout_row_height_mm) {
  ggsave(file.path(figures_dir, file), plot = plot, width = width_mm * export_scale, height = height_mm * export_scale, units = "mm", device = "pdf", bg = "white")
}

for (dataset in names(layout_objects)) {
  object <- layout_objects[[dataset]]
  UMAP <- plot_umap_cluster_labels(object, reduction = layout_reductions[[dataset]], cols = categorical_colours(levels(object$seurat_clusters)), size_scale = export_scale)
  print(UMAP)
  save_layout_panel(UMAP, paste0(dataset, "_UMAP_clusters.pdf"), layout_widths_mm[["umap"]])
}

for (dataset in names(layout_objects)) {
  object <- layout_objects[[dataset]]
  Vln_scores <- plot_score_violins_stacked(object, features = c("NeP1", "Eos_core1"), labels = c("NeP score", "EoP score"), colours = program_colours,
                                           group.by = "seurat_clusters", size_scale = export_scale, x_label_angle = cluster_label_angle(levels(object$seurat_clusters)))
  print(Vln_scores)
  save_layout_panel(Vln_scores, paste0(dataset, "_Violin_scores.pdf"), layout_widths_mm[["violins"]])
}

for (dataset in names(layout_objects)) {
  FP <- plot_featureplot_grid(layout_objects[[dataset]], genes = unlist(featureplot_genes, use.names = FALSE), reduction = layout_reductions[[dataset]], ncol = 3,
                              show_legend = FALSE, size_scale = export_scale)
  print(FP)
  save_layout_panel(FP, paste0(dataset, "_FeaturePlots.pdf"), layout_widths_mm[["featureplots"]])
}
FP_legend <- plot_featureplot_scale_legend(size_scale = export_scale)
FP_legend
save_layout_panel(FP_legend, "FeaturePlots_legend.pdf", featureplot_legend_size_mm[["width"]], featureplot_legend_size_mm[["height"]])

dotplot_table <- rbind(data.frame(dataset = "GSE180298", dotplot_summary_table(G180, dotplot_gene_groups)), data.frame(dataset = "GSE175879_NCP", dotplot_summary_table(BM, dotplot_gene_groups)),
                       data.frame(dataset = "EMTAB11188", dotplot_summary_table(EMTAB, dotplot_gene_groups)))
write.csv(transform(dotplot_table, avg_expression = round(avg_expression, 3), relative_expression = round(relative_expression, 3), pct_expressing = round(pct_expressing, 2)),
          file.path(results_dir, "DotPlot_percent_expressing.csv"), row.names = FALSE)

dotplot_headers <- list(GSE180298 = c("GSE180298", "CD34+"), GSE175879_NCP = c("GSE175879", "CD66b- CD64dim CD115-\nCD34+ and CD34dim/-"),
                        EMTAB11188 = c("E-MTAB-11188", "CD33+ CD15+ CD193-\nCD3- CD19- CD56- CD34- CD14-"))
DotPlot <- plot_transposed_dotplot(dotplot_table, dotplot_gene_groups, dotplot_group_colours, dotplot_headers, size_scale = export_scale)
DotPlot
save_layout_panel(DotPlot, "DotPlot_all_datasets.pdf", dotplot_size_mm[["width"]], dotplot_size_mm[["height"]])

save_layout_panel(plot_dotplot_legend(size_scale = export_scale), "DotPlot_legend.pdf", dotplot_legend_size_mm[["width"]], dotplot_legend_size_mm[["height"]])
