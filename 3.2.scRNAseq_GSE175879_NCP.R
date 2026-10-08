########## This code does QC, cell-cycle scoring, Harmony integration and clustering of neutrophil-committed progenitors (NCP1-4) from the bone marrow of three healthy donors (GSE175879) ##########

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.3.Functions_scRNAseq_object_generation_and_QC.R")

##### Download GEO supplementary files
dir.create(raw_data_GSE175879_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(raw_data_GSE175879_dir, "GSE175879_RAW.tar"))) {
  getGEOSuppFiles("GSE175879", makeDirectory = FALSE, baseDir = raw_data_GSE175879_dir)
}
if (length(list.files(file.path(raw_data_GSE175879_dir, "raw"), pattern = "RSEC_MolsPerCell_with_barcode\\.csv\\.gz$")) < 15) {
  untar(file.path(raw_data_GSE175879_dir, "GSE175879_RAW.tar"), exdir = file.path(raw_data_GSE175879_dir, "raw"))
}

##### Seurat object generation
### GEO sample sheet: FACS-sorted population and author cell type of each donor x Sample Tag
GSE175879_samples <- read_excel(file.path(raw_data_GSE175879_dir, "GSE175879_metadata_singleCells_15sample-rows.xlsx"))
GSE175879_samples$donor <- sub("_.*$", "", GSE175879_samples$`Sample name`)
GSE175879_samples$sample_tag <- sub("^donor[0-9]+_(SampleTag[0-9]+)_hs$", "\\1", GSE175879_samples$`Sample name`)

count_files <- list.files(file.path(raw_data_GSE175879_dir, "raw"), pattern = "RSEC_MolsPerCell_with_barcode\\.csv\\.gz$", full.names = TRUE)
BM_list <- lapply(count_files, function(count_file) {
  sample_info <- GSE175879_samples[GSE175879_samples$`processed data file` == sub("^GSM[0-9]+_", "", basename(count_file)), ]
  stopifnot(nrow(sample_info) == 1)
  seurat_object <- create_seurat_BD_rhapsody_data(count_file, project = "GSE175879", dataset_oi = "GSE175879")
  seurat_object$sample <- paste(sample_info$donor, sample_info$title, sep = "_")
  seurat_object$donor <- sample_info$donor
  seurat_object$sample_tag <- sample_info$sample_tag
  seurat_object$sorting_population <- sample_info$title
  seurat_object$author_cell_type <- sample_info$`characteristics: cell type`
  seurat_object$compartment <- "Bone marrow"
  return(seurat_object)
})
names(BM_list) <- vapply(BM_list, function(x) unique(x$sample), character(1))
BM <- merge(BM_list[[1]], y = BM_list[-1], add.cell.ids = names(BM_list), project = "GSE175879")
BM$sorting_population <- factor(BM$sorting_population, levels = c("NCP1", "NCP2", "NCP3", "NCP4", "cMOP"))
rm(BM_list)
table(BM$donor, BM$sorting_population)

##### Quality control

BM <- JoinLayers(BM)
BM[["percent.mt"]] <- PercentageFeatureSet(BM, pattern = "^MT-")
BM_list <- lapply(c("donor1", "donor2", "donor3"), function(donor_oi) {
  cells <- colnames(BM)[BM$donor == donor_oi & BM$sorting_population != "cMOP"]
  seurat_object <- CreateSeuratObject(counts = LayerData(BM, assay = "RNA", layer = "counts")[, cells, drop = FALSE],
                                      meta.data = BM@meta.data[cells, , drop = FALSE], project = donor_oi,
                                      min.cells = 10, min.features = 200)
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object[["mitoRatio"]] <- seurat_object$percent.mt / 100
  seurat_object[["log10GenesPerUMI"]] <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  return(seurat_object)
})
BM <- merge(BM_list[[1]], y = BM_list[-1])
rm(BM_list)

VlnPlot(BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "donor", layer = "counts", pt.size = 0)

BM <- subset(BM, subset = nFeature_RNA > 805 & percent.mt < 25 & nFeature_RNA < 4762)
BM <- JoinLayers(BM)
table(BM$donor, BM$sorting_population)

##### Clustering
BM <- add_cell_cycle_scores_per_donor(BM)

BM <- NormalizeData(BM, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
BM <- FindVariableFeatures(BM, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
BM <- ScaleData(BM, features = VariableFeatures(BM), vars.to.regress = c("percent.mt", "CC.Difference"), block.size = 500, verbose = FALSE)
BM <- RunPCA(BM, features = VariableFeatures(BM), npcs = 50, verbose = FALSE)
ElbowPlot(BM, ndims = 50) # 50 PCs used

BM <- RunHarmony(BM, group.by.vars = "donor", reduction.use = "pca", dims.use = 1:50,
                 reduction.save = "harmony", project.dim = FALSE, verbose = TRUE)

BM <- FindNeighbors(BM, reduction = "harmony", dims = 1:50, verbose = FALSE)
BM <- FindClusters(BM, resolution = 0.15, algorithm = 1, verbose = FALSE)
BM <- RunUMAP(BM, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony", verbose = FALSE)
DimPlot(BM, reduction = "umap.harmony", group.by = "donor")
DimPlot(BM, reduction = "umap.harmony", group.by = "sorting_population")
DimPlot(BM, reduction = "umap.harmony", group.by = "seurat_clusters", label = TRUE)

### De novo clusters vs FACS-sorted populations (kept separate; never merged)
table(BM$seurat_clusters, BM$sorting_population)

### EPX-high cluster (cluster 3)
Idents(BM) <- "seurat_clusters"
cluster3_markers <- FindMarkers(BM, ident.1 = "3", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(cluster3_markers, 30)

BM$sorting_population <- factor(BM$sorting_population, levels = c("NCP1", "NCP2", "NCP3", "NCP4"))

##### Save object
saveRDS(BM, file.path(objects_dir, "GSE175879_NCP_Harmony.rds"))
