########## This code re-clusters the steady-state neutrophils of E-MTAB-11188 (healthy BM and PB; G-CSF-treated donors removed) ##########

##### link to configuration and packages
source("0.config.R")
source("1.1.Packages.R")

##### Steady-state subset (BM and PB; G-CSF-treated donors removed)
EMTAB <- readRDS(file.path(objects_dir, "EMTAB11188_final.rds"))
EMTAB <- subset(EMTAB, subset = group %in% c("BM_steady", "PB_steady"))
EMTAB <- DietSeurat(EMTAB, layers = c("counts", "data"))
EMTAB$sample <- droplevels(EMTAB$sample)
table(EMTAB$sample, EMTAB$group)

##### Clustering 
EMTAB <- NormalizeData(EMTAB, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
EMTAB <- FindVariableFeatures(EMTAB, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

### Cell-cycle scores
s.genes <- intersect(cc.genes$s.genes, rownames(EMTAB))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(EMTAB))
EMTAB <- CellCycleScoring(EMTAB, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
EMTAB$CC.Difference <- EMTAB$S.Score - EMTAB$G2M.Score

EMTAB <- ScaleData(EMTAB, features = VariableFeatures(EMTAB), vars.to.regress = c("percent.mt", "CC.Difference"), verbose = FALSE)
EMTAB <- RunPCA(EMTAB, features = VariableFeatures(EMTAB), npcs = 50, verbose = FALSE)
ElbowPlot(EMTAB, ndims = 50) # 30 PCs retained, as for the full dataset

EMTAB <- FindNeighbors(EMTAB, reduction = "pca", dims = 1:30, verbose = FALSE)
EMTAB <- FindClusters(EMTAB, resolution = 0.15, algorithm = 1, verbose = FALSE)
EMTAB <- RunUMAP(EMTAB, reduction = "pca", dims = 1:30, reduction.name = "umap", verbose = FALSE)
DimPlot(EMTAB, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(EMTAB, reduction = "umap", group.by = "sample", split.by = "group")
table(EMTAB$seurat_clusters, EMTAB$group)

##### Save object
saveRDS(EMTAB, file.path(objects_dir, "EMTAB11188_steady_state.rds"))
