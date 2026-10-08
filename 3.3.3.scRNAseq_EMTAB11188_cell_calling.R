########## This code calls cells with emptyDrops in the unfiltered count matrices of the 8 E-MTAB-11188 runs and merges them into one Seurat object ##########
# Input : per-run kallisto | bustools exon-only count matrices (E-MTAB-11188/kb_counts, from 3.3.2)

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.3.Functions_scRNAseq_object_generation_and_QC.R")

##### Samples: ENA run, ArrayExpress sample, donor, tissue and condition (8 different individuals)
EMTAB_samples <- data.frame(
  run_accession = c("ERR7425211", "ERR7425212", "ERR7425213", "ERR7425214", "ERR7425188", "ERR7425190", "ERR7425191", "ERR7425189"),
  arrayexpress_sample = paste("Sample", 1:8),
  donor = c("UPN27", "UPN28", "UPN10", "UPN15", "UPN11", "UPN25", "UPN26", "UPN16"),
  tissue = c("PB", "PB", "BM", "BM", "PB", "PB", "PB", "PB"),
  condition = rep(c("steady_state", "GCSF_treated_donor"), each = 4),
  group = c("PB_steady", "PB_steady", "BM_steady", "BM_steady", "PB_GCSF", "PB_GCSF", "PB_GCSF", "PB_GCSF")
)

##### Seurat object generation
EMTAB_list <- lapply(seq_len(nrow(EMTAB_samples)), function(i) {
  seurat_object <- create_seurat_kb_emptyDrops(file.path(raw_data_EMTAB11188_dir, "kb_counts", EMTAB_samples$run_accession[i], "cellranger"),
                                               project = "E-MTAB-11188", lower = 100, fdr = 0.001, seed = 123)
  for (col in c("run_accession", "arrayexpress_sample", "donor", "tissue", "condition", "group")) seurat_object[[col]] <- EMTAB_samples[[col]][i]
  seurat_object$dataset <- "E-MTAB-11188"
  return(seurat_object)
})
names(EMTAB_list) <- EMTAB_samples$run_accession

EMTAB <- merge(EMTAB_list[[1]], y = EMTAB_list[-1], add.cell.ids = names(EMTAB_list))
rm(EMTAB_list)
table(EMTAB$run_accession, EMTAB$group)

##### Save object
saveRDS(EMTAB, file.path(objects_dir, "EMTAB11188_cells_raw.rds"))
