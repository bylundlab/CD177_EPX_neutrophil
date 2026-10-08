# Enzymatically active Eosinophil Peroxidase associates with CD177 and alters oxidative activity in a subset of neutrophils

Jordan Popović, Agnes Dahlstrand Rudin, Alan Bäckerholm, Vignesh Venkatakrishnan, Karin Christenson, Johan Bylund

Code for the proteomic and single-cell RNA sequencing analyses of the study:
1. paired proteomic comparison of CD177pos and CD177neg neutrophils sorted from 8 healthy donors (Figure 2, Figure 3G-I and Supplementary Data 1);
2. the neutrophil (NeP) and eosinophil (EoP) programs across the de novo clusters of three public human scRNA-seq datasets of bone marrow progenitors and neutrophils (Figure 5).

The flow cytometry, microscopy, peroxidase activity and immunoprecipitation data (Figures 1, 3A-F, 4 and 6) are not part of this repository.

#### Publication: Frontiers in Immunology (under revision; reference to be added)

#### Data
- Neutrophil proteome: Supplementary Data 1 of the article, saved as `data/Supplementary_Data_1.xlsx`. Its sheets `Quantified proteins` (3,436 protein groups quantified in all 16 samples) and `Missing values` (323 protein groups) contain the annotation and normalized abundances of all 3,759 protein groups from Proteome Discoverer 2.4; `Sample annotation` gives the TMTpro channel, donor and subset of samples 1-16. The Proteome Discoverer protein table itself (`Proteomic_base.xlsx`, sheet `Neutrophil proteome`) can be used instead (`PROTEOMICS_FILE`) and gives the same results. Raw data, search and quantification results: ProteomeXchange/PRIDE PXD085205.
- Human Protein Atlas (HPA) version 25.1: `proteinatlas.tsv` and `rna_immune_cell.tsv` (https://www.proteinatlas.org), downloaded on 2026-09-22. The scripts download the current release when the files are missing, and the tables change with each release.
- scRNA-seq: GEO GSE180298 (CD34+ bone marrow HSPCs; Ainciburu et al., eLife 2023), GEO GSE175879 (neutrophil-committed progenitors; Calzetti et al., Nat Immunol 2022) and ArrayExpress E-MTAB-11188 / ENA PRJEB48995 (bone marrow and blood neutrophils; Montaldo et al., Nat Immunol 2022).

#### Figures and tables made by the code
| Manuscript | Script | Output in `output/` |
| --- | --- | --- |
| Figure 2A-C | 2.3 | `Figures/Figure_CD177_multipanel.tiff` |
| Figure 3G, 3H, 3I | 2.4 | `Figures/Panel_HPA_heatmap.pdf`, `Panel_violin_eosinophil.pdf`, `Panel_donor_log2FC_heatmap.pdf` |
| Figure 5A-I (per dataset: UMAP, NeP/EoP score violins, FeaturePlots) | 3.4 | `Figures/{GSE180298,GSE175879_NCP,EMTAB11188}_UMAP_clusters.pdf`, `_Violin_scores.pdf`, `_FeaturePlots.pdf`, `FeaturePlots_legend.pdf` |
| Figure 5J | 3.4 | `Figures/DotPlot_all_datasets.pdf`, `DotPlot_legend.pdf`, `Results/DotPlot_percent_expressing.csv` |
| Supplementary Data 1 (limma results) | 2.1 | `Results/Proteomic analysis.xlsx` |

### How to run
- Put the input data in `data/` (layout at the end) or point to them with environment variables: `DATA_DIR`, `PROTEOMICS_FILE`, `HPA_DIR`, `GSE180298_DIR`, `GSE175879_DIR`, `EMTAB11188_DIR` (`0.config.R`). Outputs are written to `output/` (`OUTPUT_DIR`): `R_objects/`, `Results/` and `Figures/`.
- Run the scripts in the order below with the working directory set to this repository, each in a fresh R session (e.g. `Rscript 2.1.Proteomics_CD177pos_vs_CD177neg_limma.R`). Every script reads what the earlier ones saved in `output/R_objects/`.
- A fresh session matters for 3.1 and 3.3.4: with several donors, scDblFinder draws its random numbers from BiocParallel's session-level stream, which `set.seed()` does not reset, so any earlier BiocParallel call in the same session (e.g. emptyDrops in 3.3.3) changes the doublet calls.
- Running a figure script with `Rscript` also writes the plots printed for interactive use to `Rplots.pdf` in the working directory (ignored by git); `Rscript -e "source('2.3.Proteomics_figure_CD177.R')"` avoids it.
- The GEO files and the HPA tables are downloaded by the scripts when missing. The E-MTAB-11188 FASTQs (30-46 GB per run) are downloaded and counted with the bash scripts 3.3.1 and 3.3.2 (a few hours per run).

### Software
R 4.6.1 with Seurat 5.5.1 (SeuratObject 5.4.0), harmony 2.0.5, scDblFinder 1.26.7, SingleCellExperiment 1.34.0, DropletUtils 1.32.0, GEOquery 2.80.0, hdf5r 1.3.16, patchwork 1.3.2, limma 3.68.5, readxl 1.5.0.1, openxlsx 4.2.9, pheatmap 1.0.13, ggpubr 1.0.0, ggrepel 0.9.8, cowplot 1.2.0, ragg 1.5.2, RColorBrewer 1.1.3 and tidyverse 2.0.0 (ggplot2 4.0.3, dplyr 1.2.1):
```r
install.packages(c("Seurat", "harmony", "hdf5r", "patchwork", "readxl", "openxlsx", "pheatmap", "ggpubr", "ggrepel", "cowplot", "ragg", "RColorBrewer", "tidyverse", "BiocManager"))
BiocManager::install(c("limma", "scDblFinder", "SingleCellExperiment", "DropletUtils", "GEOquery"))
```
E-MTAB-11188 counting: kb-python 0.30.2 (Python 3.10; `env/kb-python_requirements.txt`), kallisto 0.52.0, bustools 0.45.1 and the human "standard" index of pachterlab/kallisto-transcriptome-indices v1 (`kb ref -d human -i index.idx -g t2g.txt`). On macOS the kallisto bundled with kb-python 0.30.2 needs an HDF5 library that Homebrew no longer ships; pass the official kallisto 0.52.0 release binary with `KALLISTO=/path/to/kallisto`.

### Code description
- `0.config.R`: input and output paths; every path can be overridden with an environment variable

**1. Packages and functions**
- `1.1.Packages.R`: packages loaded by all scripts
- `1.2.Functions_proteomics.R`: sample layout (16 abundance columns = 8 donors x CD177pos/CD177neg; odd columns CD177pos, even columns CD177neg), cutoffs of the exploratory candidates (|log2 FC| > 0.5, adjusted P < 0.2), selected proteins, figure sizes and colours, and the functions for reading the proteome table (Supplementary Data 1 or the Proteome Discoverer table), the result tables, heatmaps, paired violin plots and PDF export
- `1.3.Functions_scRNAseq_object_generation_and_QC.R`: Seurat objects from Cell Ranger h5 files, BD Rhapsody tables and kallisto | bustools matrices (emptyDrops cell calling), cell-cycle scores per donor, author cell type labels, scDblFinder per donor
- `1.4.Functions_scRNAseq_gene_sets_and_figures.R`: gene sets (EoP and NeP scores, FeaturePlot and DotPlot genes), panel sizes of the 177 x 230 mm figure and the figure functions (UMAP with the cluster numbers on the clusters, stacked score violins, FeaturePlots with one shared legend, transposed DotPlot of several datasets)

**2. Proteomics of CD177pos and CD177neg neutrophils**
- `2.1.Proteomics_CD177pos_vs_CD177neg_limma.R`: reads Supplementary Data 1 (or the Proteome Discoverer table); 3,759 proteins, of which the 3,436 with all 16 abundances are analysed (the 323 others are listed in `na_proteins`); log2 FC as the log2 ratio of the mean CD177pos and mean CD177neg abundances; paired limma model on log2 abundances (donor as blocking factor; Benjamini-Hochberg adjusted P values); 10 exploratory candidates (7 more and 3 less abundant in CD177pos neutrophils), 6 of them with adjusted P < 0.05; `Results/Proteomic analysis.xlsx` (sheets `FC_calculation`, `Differential_proteins`, `na_proteins`), whose mean abundances, log2 FC, average log2 abundance, moderated t, P and adjusted P values are those of Supplementary Data 1
- `2.2.Proteomics_HPA_eosinophil_enriched_proteins.R`: strict eosinophil-enriched HPA genes (eosinophil nTPM > 0 and >= 4x the highest nTPM of any other immune cell type; HPA class "Immune cell enriched") matched to the detected proteins by UniProt accession (12 proteins)
- `2.3.Proteomics_figure_CD177.R`: Figure 2; volcano plot (A), clustered heatmap of the selected significant and exploratory proteins (B; Euclidean distance, Ward's method) and paired violin plots of CD177, EPX and MPO (C); `Figures/Figure_CD177_multipanel.tiff` (177 x 230 mm, 600 dpi)
- `2.4.Proteomics_figure_eosinophil_proteins.R`: Figure 3G-I; the HPA eosinophil-enriched proteins (EPX, ALOX15 and catalase excluded) as a clustered heatmap, paired violin plots and a donor-wise log2 FC heatmap; `Figures/Panel_HPA_heatmap.pdf`, `Panel_violin_eosinophil.pdf` and `Panel_donor_log2FC_heatmap.pdf`, exported at their final size with 8 pt text

**3. Neutrophil and eosinophil programs in scRNA-seq data**
- `3.1.scRNAseq_GSE180298_CD34_HSPC.R`: CD34+ bone marrow HSPCs of five young healthy donors; donor-specific QC, scDblFinder per donor, 2,000 variable genes, regression of genes, UMIs and mitochondrial %, Harmony (theta = 1, 21 PCs) and clustering at resolution 0.45 (33,247 cells, 15 clusters); the authors' cell type labels are kept as external annotation only
- `3.2.scRNAseq_GSE175879_NCP.R`: FACS-sorted neutrophil-committed progenitors (NCP1-4; cMOP excluded) of three donors; QC thresholds of the original study (805 < genes < 4762, mitochondrial % < 25), cell-cycle scores per donor, 3,000 variable genes, regression of mitochondrial % and cell-cycle difference, Harmony (50 PCs) and clustering at resolution 0.15 (15,071 cells, 7 clusters)
- `3.3.1.scRNAseq_EMTAB11188_download_FASTQ.sh`: FASTQs of the 8 runs from healthy donors (bone marrow and blood at steady state, blood of G-CSF-treated donors) from ENA, with MD5 check
- `3.3.2.scRNAseq_EMTAB11188_kb_count.sh`: kallisto | bustools counting per run (exon-only "standard" workflow, as the Cell Ranger v2 matrices of GSE180298); unfiltered count matrices
- `3.3.3.scRNAseq_EMTAB11188_cell_calling.R`: cells called per run with emptyDrops (ambient profile < 100 UMIs, FDR < 0.001) and merged (68,024 cells)
- `3.3.4.scRNAseq_EMTAB11188_QC.R`: QC (>= 100 genes, mitochondrial % < 25, hemoglobin % <= 5), scDblFinder per sample, clustering (30 PCs, resolution 0.5) and removal of a low-quality and a monocyte/lymphocyte contaminant cluster (52,105 cells)
- `3.3.5.scRNAseq_EMTAB11188_steady_state.R`: bone marrow and blood cells at steady state (G-CSF-treated donors removed), re-clustered without integration (30 PCs, resolution 0.15; 24,679 cells, 7 clusters)
- `3.4.scRNAseq_figure_neutrophil_eosinophil_programs.R`: Figure 5; EoP (EPX, PRG2, PRG3, ALOX15, IL5RA) and NeP (MPO, ELANE, PRTN3, AZU1, MS4A3, CEBPE) module scores and the panels of the 177 x 230 mm scRNA-seq figure, exported at 3x their final size: for each dataset the UMAP with cluster numbers, the NeP and EoP score violins and the FeaturePlots of MPO, ELANE, PRTN3, EPX, IL5RA and PRG2; one FeaturePlot legend; the DotPlot of all clusters of the three datasets with its legends; `Figures/*.pdf` and `Results/DotPlot_percent_expressing.csv`

### Input data in data/
| Data | Files | Source | Scripts |
| --- | --- | --- | --- |
| Neutrophil proteome | `Supplementary_Data_1.xlsx` (or `Proteomic_base.xlsx` via `PROTEOMICS_FILE`) | Supplementary Data 1 of the article | 2.1 |
| Human Protein Atlas | `HPA/proteinatlas.tsv`, `HPA/rna_immune_cell.tsv` | proteinatlas.org, downloaded by 2.2 | 2.2 |
| GSE180298 | `GSE180298/GSM5460406_young1_filtered_feature_bc_matrix.h5` ... `GSM5460410_young5_...h5`, `GSE180298/GSE180298_young_metadata.txt.gz` | GEO, downloaded by 3.1 | 3.1 |
| GSE175879 | `GSE175879/raw/*_RSEC_MolsPerCell_with_barcode.csv.gz` (15 files), `GSE175879/GSE175879_metadata_singleCells_15sample-rows.xlsx` | GEO, downloaded by 3.2 | 3.2 |
| E-MTAB-11188 | `E-MTAB-11188/fastq/ERR74252*_1.fastq.gz`, `_2.fastq.gz` (8 runs) and the count matrices `E-MTAB-11188/kb_counts/ERR74252*/cellranger/` | ENA PRJEB48995, downloaded by 3.3.1 and counted by 3.3.2 | 3.3.1-3.3.3 |
| kallisto index | `kb_human_standard/index.idx`, `kb_human_standard/t2g.txt` | `kb ref -d human` | 3.3.2 |
