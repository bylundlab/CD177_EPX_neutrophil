########## This code does the paired differential abundance analysis (limma) of CD177pos vs CD177neg neutrophils in the neutrophil proteome and exports the result tables ##########
# Input : Supplementary Data 1 (sheets "Quantified proteins" and "Missing values") or the Proteome Discoverer table it was made
#         from (Proteomic_base.xlsx, sheet "Neutrophil proteome"): 8 protein annotation columns and the abundances of 16
#         samples (8 donors x CD177pos/CD177neg; sample layout in 1.2). Both give the same results.
# Proteins with a missing abundance in any sample are left out of the analysis and listed in the sheet na_proteins.
# The abundances are used exactly as they are in the input table (no further normalization).
# Fold changes: log2 of the ratio between the mean abundance of the CD177pos and of the CD177neg samples (CD177pos / CD177neg).
# P values: limma model on log2 abundances with donor as blocking factor (paired design), adjusted by Benjamini-Hochberg.
# Exploratory candidates: |log2 FC| > 0.5 and adjusted P < 0.2 (cutoffs in 1.2); significant: adjusted P < 0.05.
# The abundances and limma results are those of Supplementary Data 1.
# Run with the working directory set to this repository (paths in 0.config.R).

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.2.Functions_proteomics.R")

##### Proteome table
### Proteins with a missing abundance in any sample are kept apart (na_proteins)
data_source <- read_proteome_table(proteomics_input_file)
data <- data_source %>% drop_na
na_proteins <- anti_join(data_source, data, by = names(data_source))

### Protein annotation (bound to every calculation) and abundances of the samples
protein_info <- data %>% select(seq_len(n_annotation_columns))
base_data <- data %>% select(all_of(sample_columns))

##### Abundances of the CD177pos and CD177neg samples
CD177pos_abundances <- base_data %>% select(all_of(sample_columns[sample_info$group == "CD177pos"]))
CD177neg_abundances <- base_data %>% select(all_of(sample_columns[sample_info$group == "CD177neg"]))

##### Differential abundance (limma, paired design: donor as blocking factor, CD177neg as reference)
expr_mat_log2 <- log2(as.matrix(base_data) + eps)
group <- factor(sample_info$group, levels = c("CD177neg", "CD177pos"))
donor <- factor(sample_info$donor)
design <- model.matrix(~ donor + group)
fit <- lmFit(expr_mat_log2, design)
fit <- eBayes(fit)
limma_results <- topTable(fit, coef = "groupCD177pos", number = Inf, sort.by = "none")

### Result table: protein annotation, mean abundances of the two groups, log2 FC and the limma statistics
CD177pos_mean <- rowMeans(CD177pos_abundances, na.rm = TRUE)
CD177neg_mean <- rowMeans(CD177neg_abundances, na.rm = TRUE)
FC_calculation <- bind_cols(protein_info, tibble(`CD177pos mean` = CD177pos_mean, `CD177neg mean` = CD177neg_mean, `log2 FC` = log2(CD177pos_mean / CD177neg_mean),
                                                  `average log2 abundance` = limma_results$AveExpr, `moderated t` = limma_results$t,
                                                  p_value = limma_results$P.Value, p_adj = limma_results$adj.P.Val,
                                                  log10_pvalue = -log10(limma_results$P.Value), log10_p_adj = -log10(limma_results$adj.P.Val)))

##### Differential proteins
Upregulated_proteins <- FC_calculation %>% filter(`log2 FC` > log2fc_cutoff & p_adj < padj_cutoff)
Downregulated_proteins <- FC_calculation %>% filter(`log2 FC` < -log2fc_cutoff & p_adj < padj_cutoff)
empty_row <- as_tibble(matrix(NA, nrow = 1, ncol = ncol(FC_calculation), dimnames = list(NULL, colnames(FC_calculation))))
Differential_proteins <- bind_rows(Upregulated_proteins, empty_row, Downregulated_proteins)

### Regulation class of every protein and the labels of the selected proteins, added to the result table
FC_calculation <- FC_calculation %>% mutate(group = classify_regulation(`log2 FC`, p_adj),
                                            label_text = ifelse(Accession %in% names(selected_protein_labels), unname(selected_protein_labels[Accession]), NA_character_))

##### Export the tables
export_to_excel(file.path(results_dir, "Proteomic analysis.xlsx"), list(FC_calculation = FC_calculation, Differential_proteins = Differential_proteins, na_proteins = na_proteins))
saveRDS(list(protein_info = protein_info, base_data = base_data, FC_calculation = FC_calculation), file.path(objects_dir, "Neutrophil_proteome_limma.rds"))
