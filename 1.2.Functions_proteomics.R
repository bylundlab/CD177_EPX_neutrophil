########## Settings, protein sets and functions of the proteomic analysis ##########
# Sample layout: the 16 abundance columns are 8 donors x 2 sorted fractions; odd columns are CD177pos and even columns
# CD177neg (columns 1/2 = donor 1, 3/4 = donor 2, ...), so the two fractions of a donor are paired.
# Differential abundance: log2 FC is the log2 ratio of the mean CD177pos and the mean CD177neg abundance; P values are from a
# paired limma model on log2 abundances (donor as blocking factor). Exploratory candidates have |log2 FC| > 0.5 and an adjusted
# P value < 0.2; candidates with an adjusted P value < 0.05 are significantly different.

##### Analysis settings
sample_columns <- as.character(1:16)  # abundance columns of the input table
n_annotation_columns <- 8             # protein annotation columns in front of the abundances
log2fc_cutoff <- 0.5
padj_cutoff <- 0.2
eps <- 1e-9                           # added to the abundances before log2

### Donor and group of each abundance column
sample_info <- data.frame(sample = as.integer(sample_columns), donor = (as.integer(sample_columns) + 1) %/% 2,
                          group = ifelse(as.integer(sample_columns) %% 2 == 1, "CD177pos", "CD177neg"))

##### Protein sets (UniProt accession = label)
### Selected proteins: volcano plot, heatmap and violin plots of Figure 2 and the label_text column of the result table
selected_protein_labels <- c("P11678" = "EPX", "Q8N6Q3" = "CD177", "P13727" = "MBP", "Q9Y2Y8" = "MBP2", "Q05315" = "Galectin-10", "P16050" = "ALOX15",
                             "P05164" = "MPO", "P68871" = "Hemoglobin subunit beta", "P02042" = "Hemoglobin subunit delta", "P69905" = "Hemoglobin subunit alpha")

### Donor-wise log2 FC heatmap of Figure 3I: CD177, MBP, EPX, ALOX15, galectin-10, MPO, the hemoglobin subunits and the HPA
### eosinophil-enriched proteins
donor_fc_protein_labels <- c("P11678" = "EPX", "Q8N6Q3" = "CD177", "P13727" = "Major Basic Protein", "Q05315" = "Galectin-10", "P16050" = "ALOX15",
                             "P05164" = "MPO", "P68871" = "Hemoglobin subunit beta", "P02042" = "Hemoglobin subunit delta", "P69905" = "Hemoglobin subunit alpha",
                             "P12724" = "RNASE3", "P49006" = "MARCKSL1", "Q15744" = "CEBPE", "Q96DT0" = "LGALS12", "P19801" = "AOC1",
                             "P23610" = "F8A1", "Q15349" = "RPS6KA2", "Q13643" = "FHL3", "Q6UX27" = "VSTM1")

### Proteins of the violin plots (Figure 2C) and of the heatmap (Figure 2B)
violin_accessions <- c("Q8N6Q3", "P11678", "P05164")
heatmap_accessions <- c("P11678", "Q8N6Q3", "P13727", "Q9Y2Y8", "Q05315", "P16050", "P68871", "P02042", "P69905")

##### Figure settings
cd177_figure_size_mm <- c(width = 177, height = 230)
cd177_figure_dpi <- 600


eos_heatmap_size_mm <- c(width = 83.6, height = 78)
eos_violin_size_mm <- c(width = 83.6, height = 78)
eos_donor_fc_size_mm <- c(width = 45.6, height = 100)
font_size <- 8  # pt, all text of the Figure 3G-I panels


group_colours <- c("CD177neg" = "#0000AC", "CD177pos" = "#48A4FF")

### Theme of the Figure 3G-I panels
theme_eos <- theme_classic(base_size = font_size, base_line_size = 0.25) +
  theme(axis.text = element_text(size = font_size), strip.text = element_text(size = font_size, face = "bold"), strip.background = element_blank(),
        legend.text = element_text(size = font_size))

### Same font size for the separately built heatmap legends
theme_eos_legend <- theme_void(base_size = font_size) +
  theme(legend.position = "bottom", legend.text = element_text(size = font_size), legend.title.position = "left")

##### Analysis functions
### Proteome table with the protein annotation columns and the abundance columns "1"-"16", read from either
### - the Proteome Discoverer table (sheet "Neutrophil proteome"), or
### - Supplementary Data 1 (sheets "Quantified proteins" and "Missing values"; header in row 4, abundance columns named
###   "1 (127N)" ... "16 (135N)"). Both sheets list their proteins in the order of the Proteome Discoverer table.
read_proteome_table <- function(
    file
) {
  sheets <- excel_sheets(file)
  if ("Neutrophil proteome" %in% sheets) {
    return(read_excel(file, sheet = "Neutrophil proteome"))
  }
  if (!all(c("Quantified proteins", "Missing values") %in% sheets)) stop(file, " has neither the sheet \"Neutrophil proteome\" nor the sheets of Supplementary Data 1")
  annotation_columns <- c("Protein FDR Confidence: Combined", "Accession", "Description", "Coverage [%]", "# Peptides", "# PSMs", "# Unique Peptides", "MW [kDa]")
  read_sheet <- function(sheet) {
    tbl <- read_excel(file, sheet = sheet, skip = 3, .name_repair = "minimal")
    abundance_columns <- grep("^[0-9]+ \\(", names(tbl), value = TRUE)  # "1 (127N)" -> "1"
    tbl %>% select(all_of(annotation_columns), all_of(abundance_columns)) %>% rename_with(~ sub(" \\(.*$", "", .x), all_of(abundance_columns))
  }
  proteome_table <- bind_rows(read_sheet("Quantified proteins"), read_sheet("Missing values"))
  return(proteome_table)
}

classify_regulation <- function(
    log2fc,
    padj
) {
  regulation <- case_when(log2fc > log2fc_cutoff & padj < padj_cutoff ~ "Up",
                          log2fc < -log2fc_cutoff & padj < padj_cutoff ~ "Down",
                          TRUE ~ "NS")
  return(regulation)
}

### All tables in one Excel file, one sheet per element of the named list (the names are the sheet names)
export_to_excel <- function(
    file_name,
    sheets
) {
  wb <- createWorkbook()
  for (sheet_name in names(sheets)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, sheets[[sheet_name]])
  }
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

##### Figure functions
save_panel_pdf <- function(
    plot,
    file,
    size_mm
) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    ggsave(file, plot, width = size_mm[["width"]], height = size_mm[["height"]], units = "mm", device = quartz, type = "pdf")
  } else {
    ggsave(file, plot, width = size_mm[["width"]], height = size_mm[["height"]], units = "mm", device = cairo_pdf)
  }
}

### Log2 abundances of the given proteins (one row per accession, in the given order; the samples are the columns).
### scale_rows = TRUE converts every row to z-scores (rows without variance are set to 0)
scaled_expression_matrix <- function(
    accessions,
    protein_info,
    base_data,
    scale_rows = TRUE
) {
  accessions <- unique(accessions)
  accessions <- accessions[!is.na(accessions) & nzchar(accessions)]
  expr_tbl <- bind_cols(protein_info, base_data) %>% filter(Accession %in% accessions) %>% select(Accession, all_of(sample_columns)) %>%
    group_by(Accession) %>% summarise(across(all_of(sample_columns), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(Accession = factor(Accession, levels = accessions)) %>% arrange(Accession)
  mat <- expr_tbl %>% column_to_rownames("Accession") %>% as.matrix()
  mat_log2 <- log2(mat + eps)
  if (scale_rows) {
    mat_scaled <- t(scale(t(mat_log2)))
    mat_scaled[is.na(mat_scaled)] <- 0
    return(mat_scaled)
  }
  return(mat_log2)
}


heatmap_column_annotation <- function(
    sample_names
) {
  idx <- match(as.integer(sample_names), sample_info$sample)
  annotation <- data.frame(Group = factor(sample_info$group[idx], levels = c("CD177pos", "CD177neg")), row.names = sample_names)
  return(list(annotation = annotation, donor_labels = paste0("D", sample_info$donor[idx])))
}

### Clustered heatmap of selected proteins 
plot_protein_heatmap <- function(
    accessions,
    protein_info,
    base_data,
    labels = NULL,
    scale_rows = TRUE
) {
  mat <- scaled_expression_matrix(accessions, protein_info, base_data, scale_rows)
  column_annotation <- heatmap_column_annotation(colnames(mat))
  row_labels <- rownames(mat)
  if (!is.null(labels)) {
    hit <- intersect(names(labels), row_labels)
    row_labels[match(hit, row_labels)] <- labels[hit]
  }
  hm <- pheatmap(mat, annotation_col = column_annotation$annotation, annotation_colors = list(Group = group_colours[c("CD177pos", "CD177neg")]),
                 labels_row = row_labels, labels_col = column_annotation$donor_labels, clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", clustering_method = "ward.D2", show_colnames = TRUE, angle_col = 0, fontsize_row = 8,
                 fontsize_col = 8, cellwidth = 15, treeheight_row = 16, treeheight_col = 15, border_color = NA, silent = TRUE)
  return(hm)
}

### Data of the paired violin plots
prepare_violin_data <- function(
    accessions,
    labels,
    protein_info,
    base_data,
    FC_calculation
) {
  plot_df <- bind_cols(protein_info, base_data) %>% filter(Accession %in% accessions) %>% select(Accession, all_of(sample_columns)) %>%
    pivot_longer(cols = all_of(sample_columns), names_to = "sample_col", values_to = "abundance") %>%
    mutate(sample_col = as.integer(sample_col), group = sample_info$group[match(sample_col, sample_info$sample)],
           donor = sample_info$donor[match(sample_col, sample_info$sample)], log2_abundance = log2(abundance + eps),
           Protein = unname(labels[Accession]), group = factor(group, levels = c("CD177neg", "CD177pos")))
  stat_df <- FC_calculation %>% filter(Accession %in% accessions) %>%
    transmute(Accession, Protein = unname(labels[Accession]), p_adj = p_adj, group1 = "CD177neg", group2 = "CD177pos",
              p.label = case_when(p_adj < 0.001 ~ "P < 0.001", p_adj < 0.05 ~ paste0("P = ", formatC(p_adj, format = "f", digits = 3)), TRUE ~ "ns")) %>%
    left_join(plot_df %>% group_by(Protein) %>% summarise(y.position = max(log2_abundance, na.rm = TRUE) + 0.35, .groups = "drop"), by = "Protein")
  return(list(plot_df = plot_df, stat_df = stat_df))
}
