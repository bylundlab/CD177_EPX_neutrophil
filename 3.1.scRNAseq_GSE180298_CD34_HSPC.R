########## This code does donor-specific QC, doublet removal, Harmony integration and clustering of CD34+ bone marrow HSPCs from five young healthy donors (GSE180298) ##########

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.3.Functions_scRNAseq_object_generation_and_QC.R")

##### Download GEO supplementary files
dir.create(raw_data_GSE180298_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(raw_data_GSE180298_dir, "GSE180298_RAW.tar"))) {
  getGEOSuppFiles("GSE180298", makeDirectory = FALSE, baseDir = raw_data_GSE180298_dir)
}
if (length(list.files(raw_data_GSE180298_dir, pattern = "young[1-5]_filtered_feature_bc_matrix\\.h5$")) < 5) {
  untar(file.path(raw_data_GSE180298_dir, "GSE180298_RAW.tar"), exdir = raw_data_GSE180298_dir)
}

##### Seurat object generation
young1 <- create_seurat_10X_h5_data(file.path(raw_data_GSE180298_dir, "GSM5460406_young1_filtered_feature_bc_matrix.h5"), "GSE180298", "young1", "GSE180298", "Young healthy")
young2 <- create_seurat_10X_h5_data(file.path(raw_data_GSE180298_dir, "GSM5460407_young2_filtered_feature_bc_matrix.h5"), "GSE180298", "young2", "GSE180298", "Young healthy")
young3 <- create_seurat_10X_h5_data(file.path(raw_data_GSE180298_dir, "GSM5460408_young3_filtered_feature_bc_matrix.h5"), "GSE180298", "young3", "GSE180298", "Young healthy")
young4 <- create_seurat_10X_h5_data(file.path(raw_data_GSE180298_dir, "GSM5460409_young4_filtered_feature_bc_matrix.h5"), "GSE180298", "young4", "GSE180298", "Young healthy")
young5 <- create_seurat_10X_h5_data(file.path(raw_data_GSE180298_dir, "GSM5460410_young5_filtered_feature_bc_matrix.h5"), "GSE180298", "young5", "GSE180298", "Young healthy")

### Merge donors
G180 <- merge(young1, y = c(young2, young3, young4, young5),
              add.cell.ids = c("young1", "young2", "young3", "young4", "young5"))
rm(young1, young2, young3, young4, young5)
table(G180$donor)

##### Author annotation (external labels only)
young_meta <- read.delim(gzfile(file.path(raw_data_GSE180298_dir, "GSE180298_young_metadata.txt.gz")), header = TRUE, sep = "\t", check.names = FALSE)
G180 <- add_author_annotation(G180, young_meta, label_column = "CellType")

##### Quality control
DefaultAssay(G180) <- "RNA"
G180[["percent.mt"]] <- PercentageFeatureSet(G180, pattern = "^MT-")

VlnPlot(G180, features = c("nFeature_RNA", "percent.mt"), group.by = "donor", pt.size = 0)

G180 <- subset(G180, subset = nFeature_RNA > 200 &
                 ((donor == "young1" & nFeature_RNA < 4000 & percent.mt < 10) |
                  (donor == "young2" & nFeature_RNA < 2700 & percent.mt < 10) |
                  (donor == "young3" & nFeature_RNA < 4000 & percent.mt < 5) |
                  (donor == "young4" & nFeature_RNA < 4000 & percent.mt < 5) |
                  (donor == "young5" & nFeature_RNA < 5000 & percent.mt < 10)))

### Doublet detection and removal
G180 <- annotate_doublets_scDblFinder_per_donor(G180, seed = 123)
table(G180$donor, G180$scDblFinder.class)
prop.table(table(G180$donor, G180$scDblFinder.class), margin = 1)
G180 <- subset(G180, subset = scDblFinder.class == "singlet")

##### Clustering

G180 <- NormalizeData(G180, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
G180 <- FindVariableFeatures(G180, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
G180 <- ScaleData(G180, features = VariableFeatures(G180), vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), verbose = FALSE)
G180 <- RunPCA(G180, features = VariableFeatures(G180), npcs = 30, verbose = FALSE)
ElbowPlot(G180, ndims = 30) # 21 PCs retained

G180 <- RunHarmony(G180, group.by.vars = "donor", reduction.use = "pca", dims.use = 1:21, theta = 1,
                   reduction.save = "harmony", project.dim = FALSE, verbose = TRUE)

G180 <- FindNeighbors(G180, reduction = "harmony", dims = 1:21, verbose = FALSE)
G180 <- FindClusters(G180, resolution = 0.45, algorithm = 1, verbose = FALSE)
G180 <- RunUMAP(G180, reduction = "harmony", dims = 1:21, n.neighbors = 30, min.dist = 0.3, reduction.name = "umap", verbose = FALSE)
DimPlot(G180, reduction = "umap", group.by = "donor")
DimPlot(G180, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

### De novo clusters vs author labels (kept separate; never merged)
table(G180$seurat_clusters, G180$author_CellType)

##### Save object
saveRDS(G180, file.path(objects_dir, "GSE180298_young_final.rds"))
