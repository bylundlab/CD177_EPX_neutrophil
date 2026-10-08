########## Gene sets, figure settings and figure functions of the scRNA-seq analysis ##########

##### Gene sets
### DotPlot gene groups (top to bottom; names as labelled in the figure) and gene name colours
dotplot_gene_groups <- list(
  "Progenitor/Cycling" = c("CD34", "MKI67"),
  "Azurophil granules (NeP)" = c("MPO", "ELANE", "PRTN3"),
  "Specific granules" = c("LTF", "CAMP", "LCN2"),
  "Late maturation" = c("CD177", "FCGR3B", "MME"),
  "Eosinophil (EoP)" = c("EPX", "PRG2", "IL5RA")
)
dotplot_group_colours <- c("Progenitor/Cycling" = "#B279A2", "Azurophil granules (NeP)" = "#F58518", "Specific granules" = "#54A24B",
                           "Late maturation" = "#4C78A8", "Eosinophil (EoP)" = "#E45756")

### Program scores (EoP = Eos_core)
Eos_core_genes <- c("EPX", "PRG2", "PRG3", "ALOX15", "IL5RA")
NeP_genes <- c("MPO", "ELANE", "PRTN3", "AZU1", "MS4A3", "CEBPE") # no CD177/FCGR3B/SELL: those mark later maturation

### FeaturePlot genes
featureplot_genes <- list(neutrophil = c("MPO", "ELANE", "PRTN3"), eosinophil = c("EPX", "IL5RA", "PRG2"))

### Program colours (violin fills): the DotPlot colours of the gene groups the scores are built from
program_colours <- c(NeP = dotplot_group_colours[["Azurophil granules (NeP)"]], EoP = dotplot_group_colours[["Eosinophil (EoP)"]])

##### Figure layout and panel sizes
layout_row_height_mm <- 45
layout_widths_mm <- c(umap = 50, violins = 48, featureplots = 73)
featureplot_legend_size_mm <- c(width = 35, height = 6)
dotplot_size_mm <- c(width = 177, height = 55)
dotplot_legend_size_mm <- c(width = 47, height = 9)
export_scale <- 3

##### Colour palettes
expression_colours <- rev(brewer.pal(11, "RdYlBu"))

### Cluster colours: Dark2 colours in order for up to 8 clusters, interpolated beyond that
categorical_colours <- function(group_levels) {
  if (length(group_levels) <= 8) {
    cols <- brewer.pal(8, "Dark2")[seq_along(group_levels)]
  } else {
    cols <- colorRampPalette(brewer.pal(8, "Dark2"))(length(group_levels))
  }
  names(cols) <- group_levels
  return(cols)
}

##### Figure functions
plot_umap_cluster_labels <- function(
    object,
    reduction = "umap",
    group.by = "seurat_clusters",
    cols = NULL,
    size_scale = 1
) {
  emb <- as.data.frame(Embeddings(object, reduction)[, 1:2])
  colnames(emb) <- c("x", "y")
  emb$group <- object[[group.by, drop = TRUE]]
  bandwidth <- 0.1 * max(diff(range(emb$x)), diff(range(emb$y)))
  label_positions <- do.call(rbind, lapply(split(emb, emb$group, drop = TRUE), function(cells) {
    density <- MASS::kde2d(cells$x, cells$y, h = bandwidth, n = 100)
    peak <- which(density$z == max(density$z), arr.ind = TRUE)[1, ]
    cells[which.min((cells$x - density$x[peak[1]])^2 + (cells$y - density$y[peak[2]])^2), ]
  }))
  x0 <- min(emb$x)
  y0 <- min(emb$y)
  axis_length <- 0.15 * diff(range(emb$x))
  p <- DimPlot(object, reduction = reduction, group.by = group.by, cols = cols, pt.size = 0.1 * size_scale, raster = FALSE) + NoLegend()
  p$layers[[1]]$aes_params$stroke <- 0.17 * size_scale  # point outline scaled with the point size
  p <- p +
    ggrepel::geom_text_repel(data = label_positions, mapping = aes(x = x, y = y, label = group), inherit.aes = FALSE, size = 8 * size_scale / ggplot2::.pt,
                             fontface = "bold", bg.color = "white", bg.r = 0.15, box.padding = 0.1, point.padding = 0, min.segment.length = Inf, seed = 1) +
    annotate("segment", x = x0, xend = x0 + axis_length, y = y0, yend = y0, linewidth = 0.3 * size_scale) +
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + axis_length, linewidth = 0.3 * size_scale) +
    annotate("text", x = x0, y = y0, label = "UMAP 1", hjust = 0, vjust = 1.4, size = 6 * size_scale / ggplot2::.pt) +
    annotate("text", x = x0, y = y0, label = "UMAP 2", hjust = 0, vjust = -0.4, angle = 90, size = 6 * size_scale / ggplot2::.pt) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(plot.title = element_blank(), legend.position = "none", plot.margin = margin(1, 1, 3.5, 3.5, "mm") * size_scale)  # room for the axis labels
  return(p)
}

### Legend of a plot as a separate grob, centred on its page, e.g. a legend exported apart from the plot it belongs to
get_plot_legend <- function(plot) {
  grob <- ggplotGrob(plot + theme(legend.justification = "center"))  # Seurat plots left-justify legends
  return(grob$grobs[[which(grob$layout$name == "guide-box-right")]])
}

### DotPlot values per group
dotplot_summary_table <- function(
    object,
    gene_groups,
    group.by = "seurat_clusters"
) {
  genes <- unlist(gene_groups, use.names = FALSE)
  groups <- droplevels(as.factor(object[[group.by, drop = TRUE]]))
  present <- intersect(genes, rownames(object))
  expression <- as.matrix(FetchData(object, vars = present, layer = "data"))
  n_cells <- as.vector(table(groups))
  avg <- log1p(rowsum(expm1(expression), groups) / n_cells)  # rows in group-level order
  relative <- sweep(avg, 2, pmax(apply(avg, 2, max), .Machine$double.eps), "/")  # 0 for a gene with no expression
  pct <- 100 * rowsum((expression > 0) * 1, groups) / n_cells
  summary_table <- expand.grid(gene = genes, cluster = levels(groups), stringsAsFactors = FALSE)[, c("cluster", "gene")]
  summary_table$n_cells <- as.vector(table(groups)[summary_table$cluster])
  summary_table$gene_group <- rep(names(gene_groups), lengths(gene_groups))[match(summary_table$gene, genes)]
  summary_table$avg_expression <- NA_real_
  summary_table$relative_expression <- NA_real_
  summary_table$pct_expressing <- NA_real_
  in_object <- summary_table$gene %in% present
  index <- cbind(summary_table$cluster[in_object], summary_table$gene[in_object])
  summary_table$avg_expression[in_object] <- avg[index]
  summary_table$relative_expression[in_object] <- relative[index]
  summary_table$pct_expressing[in_object] <- pct[index]
  return(summary_table)
}

### Transposed DotPlot of several datasets side by side
plot_transposed_dotplot <- function(
    summary_table,
    gene_groups,
    group_colours,
    headers,
    max_dot_size = 3,
    group_gap = 0.5,
    size_scale = 1
) {
  genes <- unlist(gene_groups, use.names = FALSE)
  gene_group_index <- rep(seq_along(gene_groups), lengths(gene_groups))
  gene_position <- -(seq_along(genes) + (gene_group_index - 1) * group_gap)  # first gene at the top; groups group_gap rows apart
  bottom_to_top <- order(gene_position)
  label_colours <- group_colours[names(gene_groups)[gene_group_index]]
  summary_table$y <- gene_position[match(summary_table$gene, genes)]
  datasets <- names(headers)
  panels <- lapply(seq_along(datasets), function(i) {
    dataset_table <- summary_table[summary_table$dataset == datasets[i], ]
    dataset_table$cluster <- factor(dataset_table$cluster, levels = unique(dataset_table$cluster))
    p <- ggplot(dataset_table, aes(x = cluster, y = y)) +
      geom_point(aes(size = pct_expressing, colour = relative_expression), stroke = 0.16 * size_scale, na.rm = TRUE) +
      scale_radius(range = c(0, max_dot_size * size_scale), limits = c(0, 100), guide = "none") +
      scale_colour_gradientn(colours = expression_colours, limits = c(0, 1), guide = "none") +
      scale_y_continuous(breaks = gene_position[bottom_to_top], labels = genes[bottom_to_top], expand = expansion(add = 0.7)) +
      labs(title = headers[[i]][1], subtitle = headers[[i]][2], x = NULL, y = NULL) +
      theme_classic(base_size = 7 * size_scale) +
      theme(plot.title = element_text(size = 7 * size_scale, face = "bold", hjust = 0, margin = margin(b = 0.3 * size_scale, unit = "mm")),
            plot.subtitle = element_text(size = 6 * size_scale, hjust = 0, lineheight = 0.95, margin = margin(b = 1 * size_scale, unit = "mm")),
            plot.title.position = "panel", axis.text.x = element_text(size = 7 * size_scale, colour = "black", margin = margin(t = 0.5 * size_scale, unit = "mm")),
            axis.text.y = if (i == 1) element_text(size = 7 * size_scale, face = "italic", colour = label_colours[bottom_to_top], margin = margin(r = 0.5 * size_scale, unit = "mm")) else element_blank(),
            axis.ticks.y = if (i == 1) element_line(linewidth = 0.3 * size_scale) else element_blank(), axis.ticks.x = element_line(linewidth = 0.3 * size_scale),
            axis.ticks.length = unit(0.7 * size_scale, "mm"), axis.line.x = element_line(linewidth = 0.3 * size_scale),
            axis.line.y = if (i == 1) element_line(linewidth = 0.3 * size_scale) else element_blank(),  # y axis only left of the first dataset
            plot.margin = margin(0.5, 1.5, 0.5, 1.5, "mm") * size_scale)  # 3 mm between datasets
    return(p)
  })
  n_clusters <- sapply(datasets, function(dataset) length(unique(summary_table$cluster[summary_table$dataset == dataset])))
  p <- wrap_plots(panels, nrow = 1, widths = n_clusters)  # equal column width in every dataset
  return(p)
}

### Legends of plot_transposed_dotplot() (same max_dot_size and size_scale)
plot_dotplot_legend <- function(
    max_dot_size = 3,
    size_scale = 1
) {
  legend_values <- data.frame(x = 1:5, y = 0, pct = c(0, 25, 50, 75, 100), relative = seq(0, 1, length.out = 5))
  p <- ggplot(legend_values, aes(x = x, y = y)) +
    geom_point(aes(size = pct, colour = relative), stroke = 0.16 * size_scale) +
    scale_radius(range = c(0, max_dot_size * size_scale), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = "Percent Expressed") +
    scale_colour_gradientn(colours = expression_colours, limits = c(0, 1), breaks = c(0, 1), labels = c("0", "max"), name = "Relative Expression") +
    guides(size = guide_legend(order = 1, override.aes = list(colour = "black")),  # colour bar 15 mm long
           colour = guide_colourbar(order = 2, theme = theme(legend.key.width = unit(15 * size_scale, "mm"), legend.key.height = unit(1.5 * size_scale, "mm")))) +
    theme(legend.position = "right", legend.direction = "horizontal", legend.box = "horizontal", legend.title.position = "top", legend.text.position = "bottom",
          legend.title = element_text(size = 7 * size_scale, margin = margin(b = 0.5 * size_scale, unit = "mm")), legend.text = element_text(size = 7 * size_scale),
          legend.key = element_blank(), legend.key.spacing.x = unit(1 * size_scale, "mm"), legend.spacing.x = unit(3 * size_scale, "mm"), legend.margin = margin(0, 0, 0, 0))
  return(get_plot_legend(p))
}

### Expression FeaturePlot 
plot_expression_featureplot <- function(
    object,
    gene,
    reduction = "umap",
    show_legend = TRUE,
    size_scale = 1
) {
  if (gene %in% rownames(object)) {
    p <- FeaturePlot(object, features = gene, reduction = reduction, pt.size = 0.05 * size_scale, order = TRUE, raster = FALSE) +
      scale_colour_gradientn(colours = c("grey90", expression_colours), breaks = scales::breaks_pretty(n = 2)) +
      coord_fixed()  # same shape as the UMAP
    p$layers[[1]]$aes_params$stroke <- 0.085 * size_scale  # point outline scaled with the point size (default 0.5 dominates small points)
  } else {
    p <- ggplot() + annotate("text", x = 0, y = 0, label = "not in object", size = 8 * size_scale / ggplot2::.pt, colour = "grey50") + ggtitle(gene)
  }
  p <- p + theme_void() +
    theme(plot.title = element_text(size = 8 * size_scale, face = "italic", hjust = 0.5), legend.position = if (show_legend) "right" else "none",
          legend.title = element_blank(), legend.text = element_text(size = 7 * size_scale), legend.key.width = unit(1 * size_scale, "mm"),
          legend.key.height = unit(2 * size_scale, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5 * size_scale, "mm"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm") * size_scale)
  return(p)
}

### Grid of expression FeaturePlots (one panel per gene, filled row by row)
plot_featureplot_grid <- function(
    object,
    genes,
    reduction = "umap",
    ncol = 3,
    show_legend = TRUE,
    size_scale = 1
) {
  panels <- lapply(genes, function(gene) plot_expression_featureplot(object, gene = gene, reduction = reduction, show_legend = show_legend, size_scale = size_scale))
  p <- wrap_plots(panels, ncol = ncol)
  return(p)
}

### One colour legend for FeaturePlots
plot_featureplot_scale_legend <- function(
    size_scale = 1
) {
  p <- ggplot(data.frame(value = seq(0, 1, length.out = 256)), aes(x = value, y = 0, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colours = c("grey90", expression_colours), guide = "none") +
    scale_x_continuous(breaks = c(0, 1), labels = c("0", "max"), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Expression") +
    theme_void() +
    theme(axis.title.y = element_text(size = 8 * size_scale, angle = 0, vjust = 0.5, margin = margin(r = 1.5 * size_scale, unit = "mm")),
          axis.text.x = element_text(size = 7 * size_scale, margin = margin(t = 0.5 * size_scale, unit = "mm")),
          plot.margin = margin(0.5, 3, 0.5, 0.5, "mm") * size_scale)  # right margin for the "max" label
  return(p)
}

### Scores per group as stacked violins 
plot_score_violins_stacked <- function(
    object,
    features,
    labels,
    colours,
    group.by,
    size_scale = 1,
    x_label_angle = 0
) {
  p <- VlnPlot(object, features = features, group.by = group.by, stack = TRUE, flip = TRUE, fill.by = "feature", pt.size = 0)
  p$data$feature <- factor(p$data$feature, levels = features, labels = labels)
  p$layers[[1]]$aes_params$linewidth <- 0.25 * size_scale  # thin violin outlines
  x_text <- if (x_label_angle == 90) {
    element_text(size = 8 * size_scale, angle = 90, hjust = 1, vjust = 0.5, colour = "black")
  } else {
    element_text(size = 8 * size_scale, angle = 0, hjust = 0.5, vjust = 1, colour = "black")
  }
  p <- suppressMessages(p +
    scale_fill_manual(values = setNames(unname(colours), labels)) +
    scale_y_continuous(expand = expansion(mult = 0.05), breaks = scales::breaks_pretty(n = 3)) +
    facet_grid(feature ~ ., scales = "free_y", switch = "y")) +  # score names on the left, as axis titles
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(size = 8 * size_scale, face = "plain", angle = 90, margin = margin(r = 1 * size_scale, unit = "mm")),
          axis.text.x = x_text, axis.text.y = element_text(size = 7 * size_scale, colour = "black"), axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.3 * size_scale), axis.ticks.length = unit(0.7 * size_scale, "mm"),
          panel.background = element_rect(fill = NA, colour = "black", linewidth = 0.3 * size_scale), panel.spacing = unit(1 * size_scale, "mm"),
          plot.margin = margin(1, 1, 1, 1, "mm") * size_scale)
  return(p)
}
