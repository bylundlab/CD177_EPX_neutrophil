########## Configuration file ##########
# Source this file at the beginning of each script to define all input and output paths. Run every script with the
# working directory set to this repository. Input data are read from data/ and outputs are written to output/; every
# path can be overridden with an environment variable, e.g. Sys.setenv(OUTPUT_DIR = "/some/other/place") before sourcing.

##### Base directories
base_dir <- getwd()
data_dir <- Sys.getenv("DATA_DIR", file.path(base_dir, "data"))
output_dir <- Sys.getenv("OUTPUT_DIR", file.path(base_dir, "output"))

##### Input data
### Proteomics - protein abundances of CD177pos and CD177neg neutrophils from 8 donors: Supplementary Data 1 of the article, or
### the Proteome Discoverer protein table it was made from (Proteomic_base.xlsx, sheet "Neutrophil proteome")
proteomics_input_file <- Sys.getenv("PROTEOMICS_FILE", file.path(data_dir, "Supplementary_Data_1.xlsx"))

### Human Protein Atlas - annotation table and immune-cell RNA data (downloaded by 2.2 when missing)
hpa_dir <- Sys.getenv("HPA_DIR", file.path(data_dir, "HPA"))

### scRNA-seq
# GSE180298 - CD34+ HSPCs, young healthy donors (Cell Ranger filtered matrices; downloaded from GEO by 3.1 when missing)
raw_data_GSE180298_dir <- Sys.getenv("GSE180298_DIR", file.path(data_dir, "GSE180298"))
# GSE175879 - neutrophil-committed progenitors (BD Rhapsody RSEC tables; downloaded from GEO by 3.2 when missing)
raw_data_GSE175879_dir <- Sys.getenv("GSE175879_DIR", file.path(data_dir, "GSE175879"))
# E-MTAB-11188 - bone marrow and blood neutrophils (ENA FASTQs, counted with kallisto | bustools by 3.3.1 and 3.3.2)
raw_data_EMTAB11188_dir <- Sys.getenv("EMTAB11188_DIR", file.path(data_dir, "E-MTAB-11188"))

##### Output directories
objects_dir <- file.path(output_dir, "R_objects")
results_dir <- file.path(output_dir, "Results")
figures_dir <- file.path(output_dir, "Figures")

### Create output directories if missing
for (d in c(objects_dir, results_dir, figures_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
