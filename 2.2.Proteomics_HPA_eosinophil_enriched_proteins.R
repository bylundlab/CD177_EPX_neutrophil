########## This code selects the eosinophil-enriched proteins of the Human Protein Atlas (HPA) among the proteins detected in the neutrophil proteome ##########
# Input : HPA annotation table (proteinatlas.tsv) and HPA immune-cell RNA data (rna_immune_cell.tsv), downloaded when missing;
#         the differential abundance table of 2.1 (R_objects/Neutrophil_proteome_limma.rds)
# Strict eosinophil-enriched genes: eosinophil RNA expression (nTPM) > 0 and at least 4 times the highest nTPM of any other
# immune cell type, classified by the HPA as "Immune cell enriched"; matched to the detected proteins by UniProt accession.
# The HPA files change with each HPA release; the results here were made with HPA version 25.1 (files downloaded on 2026-09-22).
# Run with the working directory set to this repository (paths in 0.config.R).

##### link to configuration and packages
source("0.config.R")
source("1.1.Packages.R")

##### Load the results of 2.1
proteome <- readRDS(file.path(objects_dir, "Neutrophil_proteome_limma.rds"))
protein_info <- proteome$protein_info
FC_calculation <- proteome$FC_calculation

##### Download the HPA tables (skipped when already present)
dir.create(hpa_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(hpa_dir, "proteinatlas.tsv"))) {
  download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", destfile = file.path(hpa_dir, "proteinatlas.tsv.zip"), mode = "wb")
  unzip(file.path(hpa_dir, "proteinatlas.tsv.zip"), exdir = hpa_dir, overwrite = TRUE)
}
if (!file.exists(file.path(hpa_dir, "rna_immune_cell.tsv"))) {
  download.file("https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip", destfile = file.path(hpa_dir, "rna_immune_cell.tsv.zip"), mode = "wb")
  unzip(file.path(hpa_dir, "rna_immune_cell.tsv.zip"), exdir = hpa_dir, overwrite = TRUE)
}

##### HPA annotation table and immune-cell RNA data
hpa <- read_tsv(file.path(hpa_dir, "proteinatlas.tsv"), show_col_types = FALSE)
immune <- read_tsv(file.path(hpa_dir, "rna_immune_cell.tsv"), show_col_types = FALSE)

##### Strict eosinophil-enriched genes
### Immune-cell expression in wide format (one column per cell type)
immune_wide <- immune %>% filter(`Immune cell` != "total PBMC") %>% select(Gene, `Immune cell`, nTPM) %>% pivot_wider(names_from = `Immune cell`, values_from = nTPM)

### Highest expression of the other cell types and the eosinophil ratio to it
other_cells <- setdiff(names(immune_wide), c("Gene", "eosinophil"))
immune_eos <- immune_wide %>% rowwise() %>% mutate(max_other = max(c_across(all_of(other_cells)), na.rm = TRUE), eos_vs_next = eosinophil / max_other) %>% ungroup()

### Strict criterion: eosinophil nTPM > 0 and at least 4 times the highest nTPM of the other cell types; HPA class
### "Immune cell enriched"
eos_strict <- immune_eos %>% filter(eosinophil > 0, eosinophil >= 4 * max_other) %>%
  left_join(hpa %>% select(Gene_symbol = Gene, Ensembl, Uniprot, `RNA blood cell specificity`), by = c("Gene" = "Ensembl")) %>%
  filter(`RNA blood cell specificity` == "Immune cell enriched")

##### Match the strict eosinophil genes to the detected proteins
### One row per UniProt accession (isoform suffixes removed)
eos_strict_uniprot <- eos_strict %>% separate_rows(Uniprot, sep = "[,;]") %>% mutate(Uniprot = trimws(Uniprot), Uniprot = sub("-[0-9]+$", "", Uniprot)) %>% filter(!is.na(Uniprot), Uniprot != "")
eos_strict_detected <- FC_calculation %>% mutate(Accession_match = sub("-[0-9]+$", "", Accession)) %>% inner_join(eos_strict_uniprot, by = c("Accession_match" = "Uniprot"))

### Protein names from the Proteome Discoverer annotation
protein_names <- protein_info %>% transmute(Accession_match = sub("-[0-9]+$", "", Accession), Protein_name = Description) %>% distinct(Accession_match, .keep_all = TRUE)

### Strict eosinophil-enriched proteins, sorted by log2 FC, one row per UniProt accession
eos_strict_table_unique <- eos_strict_detected %>% left_join(protein_names, by = "Accession_match") %>% mutate(Protein_name = sub(" OS=.*$", "", Protein_name)) %>%
  select(Accession, Gene_symbol, Protein_name, eosinophil, max_other, eos_vs_next, `log2 FC`, p_value, p_adj) %>% arrange(desc(`log2 FC`)) %>% distinct(Accession, .keep_all = TRUE)
eos_strict_table_unique

##### Save the table for the eosinophil panels (2.4)
saveRDS(eos_strict_table_unique, file.path(objects_dir, "HPA_eosinophil_enriched_proteins.rds"))
