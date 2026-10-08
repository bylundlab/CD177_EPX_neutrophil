########## Functions for Seurat object generation and quality control of the scRNA-seq datasets ##########

### Seurat object generation from a Cell Ranger filtered feature-barcode matrix (.h5)
create_seurat_10X_h5_data <- function(
    path_to_h5_file,
    project,
    donor_oi,
    dataset_oi,
    cohort_oi
) {
  counts <- Read10X_h5(path_to_h5_file)
  # multi-modal h5 files return a list; keep gene expression only
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }
  seurat_object <- CreateSeuratObject(counts = counts, project = project)
  seurat_object$donor <- donor_oi
  seurat_object$dataset <- dataset_oi
  seurat_object$cohort <- cohort_oi
  return(seurat_object)
}

### Seurat object generation from a BD Rhapsody RSEC molecules-per-cell table (.csv.gz: Cell_Index, barcode_sequences,
### then one column per gene); all genes and cells are kept
create_seurat_BD_rhapsody_data <- function(
    path_to_csv_file,
    project,
    dataset_oi
) {
  mols <- data.table::fread(path_to_csv_file, data.table = FALSE)
  cell_ids <- as.character(mols$Cell_Index)
  stopifnot(!anyDuplicated(cell_ids))
  counts <- t(as.matrix(mols[, 3:ncol(mols)]))
  rownames(counts) <- names(mols)[3:ncol(mols)]
  colnames(counts) <- cell_ids
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  cell_metadata <- data.frame(barcode_sequences = mols$barcode_sequences, Cell_Index = mols$Cell_Index, row.names = cell_ids)
  seurat_object <- CreateSeuratObject(counts = counts, meta.data = cell_metadata, project = project, min.cells = 0, min.features = 0)
  seurat_object$dataset <- dataset_oi
  return(seurat_object)
}

### Seurat object generation from an unfiltered kallisto | bustools count matrix (Cell Ranger layout), keeping the
### barcodes called as cells by emptyDrops (ambient cutoff and FDR are set in the analysis script)
create_seurat_kb_emptyDrops <- function(
    path_to_kb_matrix,
    project,
    lower,
    fdr,
    seed = 123
) {
  counts <- Read10X(path_to_kb_matrix)
  set.seed(seed)
  empty_drops <- emptyDrops(counts, lower = lower)
  is_cell <- !is.na(empty_drops$FDR) & empty_drops$FDR < fdr
  seurat_object <- CreateSeuratObject(counts = counts[, is_cell], project = project)
  return(seurat_object)
}

### Cell-cycle scores (Seurat cc.genes) computed separately for each donor; adds S.Score, G2M.Score, Phase and
### CC.Difference (S.Score - G2M.Score)
add_cell_cycle_scores_per_donor <- function(seurat_object) {
  scores <- lapply(unique(seurat_object$donor), function(donor_oi) {
    donor_object <- subset(seurat_object, cells = colnames(seurat_object)[seurat_object$donor == donor_oi])
    donor_object <- NormalizeData(donor_object, verbose = FALSE)
    donor_object <- CellCycleScoring(donor_object, s.features = intersect(cc.genes$s.genes, rownames(donor_object)),
                                     g2m.features = intersect(cc.genes$g2m.genes, rownames(donor_object)), set.ident = FALSE)
    donor_object@meta.data[, c("S.Score", "G2M.Score", "Phase"), drop = FALSE]
  })
  scores <- do.call(rbind, scores)[colnames(seurat_object), ]
  seurat_object$S.Score <- scores$S.Score
  seurat_object$G2M.Score <- scores$G2M.Score
  seurat_object$Phase <- scores$Phase
  seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
  return(seurat_object)
}

### Add published author cell type labels as external annotation (never used for PCA/clustering)
# author metadata rownames are "<barcode>_<donor>"; Seurat cell names are "<donor>_<barcode>-1"
add_author_annotation <- function(
    seurat_object,
    author_metadata,
    label_column,
    new_column = "author_CellType"
) {
  author_metadata$meta_cell <- rownames(author_metadata)
  author_metadata$donor <- sub("^.*_(young[1-5])$", "\\1", author_metadata$meta_cell)
  author_metadata$barcode <- sub("_young[1-5]$", "", author_metadata$meta_cell)
  author_metadata$Seurat_cell <- paste0(author_metadata$donor, "_", author_metadata$barcode, "-1")
  idx <- match(colnames(seurat_object), author_metadata$Seurat_cell)
  seurat_object[[new_column]] <- author_metadata[[label_column]][idx]
  return(seurat_object)
}

### Doublet detection with scDblFinder, run independently per donor
# With several donors, scDblFinder draws from BiocParallel's session-level random stream, which set.seed() does not
# reset: results are reproducible only when no other BiocParallel call (e.g. emptyDrops) ran earlier in the R session
annotate_doublets_scDblFinder_per_donor <- function(
    seurat_object,
    seed = 123
) {
  seurat_object <- JoinLayers(seurat_object)
  sce <- as.SingleCellExperiment(seurat_object)
  set.seed(seed)
  sce <- scDblFinder(sce, samples = "donor", clusters = TRUE)
  seurat_object$scDblFinder.score <- sce$scDblFinder.score
  seurat_object$scDblFinder.class <- sce$scDblFinder.class
  return(seurat_object)
}
