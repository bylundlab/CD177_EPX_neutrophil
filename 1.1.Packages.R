##### Packages #####
### scRNA-seq
library(Seurat)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(DropletUtils)
library(GEOquery)
library(hdf5r)
library(patchwork)

### Proteomics
library(limma)
library(readxl)
library(openxlsx)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(grid)
library(ragg)
library(RColorBrewer)

### Data handling and plotting (attached last, so that the dplyr verbs are not masked by the Bioconductor packages)
library(tidyverse)
