########## This code does per-sample QC, doublet removal, clustering and removal of low-quality and contaminant clusters for the E-MTAB-11188 neutrophils ##########

##### link to configuration, packages and functions
source("0.config.R")
source("1.1.Packages.R")
source("1.3.Functions_scRNAseq_object_generation_and_QC.R")

##### Load cell-called object
EMTAB <- readRDS(file.path(objects_dir, "EMTAB11188_cells_raw.rds"))
EMTAB <- JoinLayers(EMTAB)
table(EMTAB$run_accession, EMTAB$group)

##### Quality control
run_order <- c("ERR7425213", "ERR7425214", "ERR7425211", "ERR7425212", "ERR7425188", "ERR7425190", "ERR7425191", "ERR7425189")
EMTAB$sample <- factor(EMTAB$run_accession, levels = run_order)

EMTAB[["percent.mt"]] <- PercentageFeatureSet(EMTAB, pattern = "^MT-")
EMTAB[["percent.hb"]] <- PercentageFeatureSet(EMTAB, features = intersect(c("HBA1", "HBA2", "HBB"), rownames(EMTAB)))

VlnPlot(EMTAB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), group.by = "sample", layer = "counts", pt.size = 0, ncol = 2)
FeatureScatter(EMTAB, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "sample", split.by = "sample", ncol = 4, pt.size = 0.1, plot.cor = FALSE)

cells_called <- table(EMTAB$sample)
EMTAB <- subset(EMTAB, subset = nFeature_RNA >= 100 & percent.mt < 25 & percent.hb <= 5)
cells_after_qc <- table(EMTAB$sample)

### Doublet detection (scDblFinder per sample; one donor per run) and removal
EMTAB <- annotate_doublets_scDblFinder_per_donor(EMTAB, seed = 123)
table(EMTAB$sample, EMTAB$scDblFinder.class)
EMTAB <- subset(EMTAB, subset = scDblFinder.class == "singlet")

cbind(called = cells_called, after_qc = cells_after_qc, singlets = table(EMTAB$sample))

##### Clustering

EMTAB <- NormalizeData(EMTAB, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
EMTAB <- FindVariableFeatures(EMTAB, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

s.genes <- intersect(cc.genes$s.genes, rownames(EMTAB))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(EMTAB))
EMTAB <- CellCycleScoring(EMTAB, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
EMTAB$CC.Difference <- EMTAB$S.Score - EMTAB$G2M.Score

EMTAB <- ScaleData(EMTAB, features = VariableFeatures(EMTAB), vars.to.regress = c("percent.mt", "CC.Difference"), verbose = FALSE)
EMTAB <- RunPCA(EMTAB, features = VariableFeatures(EMTAB), npcs = 50, verbose = FALSE)
ElbowPlot(EMTAB, ndims = 50) # 30 PCs retained

EMTAB <- FindNeighbors(EMTAB, reduction = "pca", dims = 1:30, verbose = FALSE)
EMTAB <- FindClusters(EMTAB, resolution = 0.5, algorithm = 1, verbose = FALSE)
EMTAB <- RunUMAP(EMTAB, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
DimPlot(EMTAB, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(EMTAB, reduction = "umap", group.by = "sample", split.by = "group")

##### Low-quality and contaminant clusters
### DEGs per cluster
Idents(EMTAB) <- "seurat_clusters"
markers <- FindAllMarkers(EMTAB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)

VlnPlot(EMTAB, features = c("nFeature_RNA", "percent.mt"), pt.size = 0)

table(EMTAB$seurat_clusters, EMTAB$sample)

### Marker genes of neutrophil stages (as defined for this dataset by Montaldo et al., Nat Immunol 2022) and of the other
### blood and bone marrow lineages
stage_and_lineage_markers <- unique(c(
  "MPO", "ELANE", "AZU1", "DEFA4", "MKI67", "TOP2A", # Precursors
  "LTF", "CAMP", "LCN2", # Early immature (specific granules)
  "MMP9", "CTSB", # Immature (gelatinase granules)
  "SELL", "MME", "CXCR4", # Mature BM
  "NAMPT", "CXCR2", "SOD2", # Mature PB
  "SERPINA1", "CR1", "CD177", "CD14", # G-CSF-induced
  "IFIT1", "ISG15", "GBP1", # IFN response
  "FCGR3B", "S100A8", # Neutrophils
  "CD34", "CEBPA", # GMPs
  "CD300E", "EREG", "VCAN", # Monocytes
  "CLEC9A", "FLT3", # DCs
  "CD3E", "TRAC", "CD8A", # T cells
  "CD19", "MS4A1", # B cells
  "IGKC", "IGHG2", # PCs
  "CLC", "ADGRE1", "EPX", # Eosinophils
  "CPA3", "MS4A2", "FCER1A", "HDC", # Basophils and mast cells
  "HBB" # Erythrocytes
))
DotPlot(EMTAB, features = stage_and_lineage_markers) + RotatedAxis()

### Remove cluster 4 (low quality: median mitochondrial % 16.5, few genes; mostly ERR7425188 and ERR7425191) and
### cluster 14 (contaminants: monocyte genes VCAN, CD14, HLA-DR and some NK/T-cell genes, few neutrophil genes)
EMTAB <- subset(EMTAB, subset = !(seurat_clusters %in% c("4", "14")))
table(EMTAB$seurat_clusters, EMTAB$group)

##### Save object
saveRDS(EMTAB, file.path(objects_dir, "EMTAB11188_final.rds"))
