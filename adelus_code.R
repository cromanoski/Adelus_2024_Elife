# Below is the code used to generate the final RDS for this study
# directory locations need to be adjusted to where files are located. 
# Load necessary libraries 
library(hdf5r) 
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86) 
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ggplot2)
library(data.table)
library(dittoSeq) 
library(DoubletFinder)
library(dplyr)
library(GenomeInfoDb)
library(devtools)
library(glmGamPoi)
set.seed(1234)

# Extract gene annotations from EnsDb, and change to UCSC style
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotation) <- "UCSC"

# Load the RNA multiome data 
counts1 <- Read10X_h5("/directory/snt/outs/filtered_feature_bc_matrix.h5") 
fragpath1 <- "/directory/snt_atac_fragments.tsv.gz"
nt = CreateSeuratObject(counts = counts1$`Gene Expression`, assay = "RNA")
counts2 <- Read10X_h5("/directory/sil1b/outs/filtered_feature_bc_matrix.h5")
fragpath2 <- "/directory/sil1b_atac_fragments.tsv.gz"
il1b <- CreateSeuratObject(counts = counts2$`Gene Expression`, assay = "RNA")
counts3 <- Read10X_h5("/directory/CR74-75_1/outs/filtered_feature_bc_matrix.h5") 
fragpath3 <- "/directory_atac_fragments.tsv.gz"
lane1 <- CreateSeuratObject(counts = counts3$`Gene Expression`, assay = "RNA")
counts4 <- Read10X_h5("/directory/CR74-75_2/outs/filtered_feature_bc_matrix.h5")
fragpath4 <- "/directory/CR74-75_2_atac_fragments.tsv.gz"
lane2 <- CreateSeuratObject(counts = counts4$`Gene Expression`, assay = "RNA")
counts5 <- Read10X_h5("/directory/CR74-75_3/outs/filtered_feature_bc_matrix.h5")
fragpath5 <- "/directory/CR74-75_3_atac_fragments.tsv.gz"
lane3 <- CreateSeuratObject(counts = counts5$`Gene Expression`, assay = "RNA")
counts6 <- Read10X_h5("/directory/CR74-75_4/outs/filtered_feature_bc_matrix.h5")
fragpath6 <- "/directory/CR74-75_4_atac_fragments.tsv.gz"
lane4 <- CreateSeuratObject(counts = counts6$`Gene Expression`, assay = "RNA")
counts7 <- Read10X_h5("/directory/CR74-75_5/outs/filtered_feature_bc_matrix.h5")
fragpath7 <- "/directory/CR74-75_5_atac_fragments.tsv.gz"
lane5 <- CreateSeuratObject(counts = counts7$`Gene Expression`, assay = "RNA")

# Filter out low quality cells in RNA objects
nt[["percent.mt"]] <- PercentageFeatureSet(nt, pattern = "^MT-")
nt <- subset(x = nt, subset = nCount_RNA > 500 & percent.mt < 20) 
il1b[["percent.mt"]] <- PercentageFeatureSet(il1b, pattern = "^MT-")
il1b <- subset(x = il1b, subset = nCount_RNA > 500 & percent.mt < 20) 
lane1[["percent.mt"]] <- PercentageFeatureSet(lane1, pattern = "^MT-")
lane1 <- subset(x = lane1, subset = nCount_RNA > 500 & percent.mt < 20) 
lane2[["percent.mt"]] <- PercentageFeatureSet(lane2, pattern = "^MT-")
lane2 <- subset(x = lane2, subset = nCount_RNA > 500 & percent.mt < 20) 
lane3[["percent.mt"]] <- PercentageFeatureSet(lane3, pattern = "^MT-")
lane3 <- subset(x = lane3, subset = nCount_RNA > 500 & percent.mt < 20) 
lane4[["percent.mt"]] <- PercentageFeatureSet(lane4, pattern = "^MT-")
lane4 <- subset(x = lane4, subset = nCount_RNA > 500 & percent.mt < 20) 
lane5[["percent.mt"]] <- PercentageFeatureSet(lane5, pattern = "^MT-")
lane5 <- subset(x = lane5, subset = nCount_RNA > 500 & percent.mt < 20) 

# Import demuxlet files and subset RNA objects for singlets 
nt = importDemux(nt,demuxlet.best = "/directory/snt.demuxlet.concat.best")
nt <- subset(nt, subset = demux.doublet.call == "SNG") 
il1b = importDemux(il1b,demuxlet.best = "/directory/sil1b.demuxlet.concat.best")
il1b <- subset(il1b, subset = demux.doublet.call == "SNG") 
lane1 = importDemux(lane1,demuxlet.best = "/directory/CR74-75_1.demuxlet.concat.best")
lane1 <- subset(lane1, subset = demux.doublet.call == "SNG") 
lane2 = importDemux(lane2,demuxlet.best = "/directory/CR74-75_2.demuxlet.concat.best")
lane2 <- subset(lane2, subset = demux.doublet.call == "SNG") 
lane3 = importDemux(lane3,demuxlet.best = "/directory/CR74-75_3.demuxlet.concat.best")
lane3 = subset(lane3, subset = demux.doublet.call == "SNG") 
lane4 = importDemux(lane4,demuxlet.best = "/directory/CR74-75_4.demuxlet.concat.best")
lane4 <- subset(lane4, subset = demux.doublet.call == "SNG") 
lane5 = importDemux(lane5,demuxlet.best = "/directory/CR74-75_5.demuxlet.concat.best")
lane5 <- subset(lane5, subset = demux.doublet.call == "SNG") 

# Save lane IDs before integrating RNA objects
nt$condition = "nt" 
il1b$condition = "il1b" 
lane1$condition = "lane1"
lane2$condition = "lane2"
lane3$condition = "lane3"
lane4$condition = "lane4"
lane5$condition = "lane5"

# Add lane IDs to barcodes to avoid duplicates during merging and integration 
x = gsub("^", paste0(nt$condition, "_"), colnames(nt))
nt = RenameCells(nt, new.names = x)
x = gsub("^", paste0(il1b$condition, "_"), colnames(il1b))
il1b = RenameCells(il1b, new.names = x)
x = gsub("^", paste0(lane1$condition, "_"), colnames(lane1))
lane1 = RenameCells(lane1, new.names = x)
x = gsub("^", paste0(lane2$condition, "_"), colnames(lane2))
lane2 = RenameCells(lane2, new.names = x)
x = gsub("^", paste0(lane3$condition, "_"), colnames(lane3))
lane3 = RenameCells(lane3, new.names = x)
x = gsub("^", paste0(lane4$condition, "_"), colnames(lane4))
lane4  = RenameCells(lane4, new.names = x)
x = gsub("^", paste0(lane5$condition, "_"), colnames(lane5))
lane5 = RenameCells(lane5, new.names = x)

# Clean up sample IDs
x = gsub("\\..*", "", nt$Sample)
nt = AddMetaData(object = nt, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", il1b$Sample)
il1b = AddMetaData(object = il1b, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane1$Sample)
lane1 = AddMetaData(object = lane1, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane2$Sample)
lane2 = AddMetaData(object = lane2, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane3$Sample)
lane3 = AddMetaData(object = lane3, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane4$Sample)
lane4 = AddMetaData(object = lane4, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane5$Sample)
lane5 = AddMetaData(object = lane5, metadata = x, col.name = "Donor")

# Add treatment based on lane and donor information for each cell
nt$Treatment = "ctrl_6hr"
il1b$Treatment = "il1b_6hr"
lane1_V105_J98_U114 = subset(lane1, subset = Donor == "V105" 
                             | Donor == "J98" 
                             | Donor == "U114")
cells_lane1_V105_J98_U114 = colnames(lane1_V105_J98_U114)
lane1$Treatment <- ifelse(colnames(lane1) %in% cells_lane1_V105_J98_U114, "siSCR_7d", "ctrl_tx_7d")
lane1$Donor = ifelse(lane1$Donor == "TELOS", "WT_Telo", lane1$Donor)
lane2_V105_J98_U114 = subset(lane2, subset = Donor == "V105" 
                             | Donor == "J98" 
                             | Donor == "U114")
lane2_J76_WT_Telos = subset(lane2, subset = Donor == "J76" 
                            | Donor == "TELOS")
cells_lane2_V105_J98_U114 = colnames(lane2_V105_J98_U114)
cells_lane2_J76_WT_Telos = colnames(lane2_J76_WT_Telos)
lane2$Treatment <- ifelse(colnames(lane2) %in% cells_lane2_V105_J98_U114, "siRNA_7d", 
                          ifelse(colnames(lane2) %in% cells_lane2_J76_WT_Telos, "TGFb_7d", "TGFb_4d"))
lane2$Donor = ifelse(lane2$Donor == "TELOS", "WT_Telo", lane2$Donor)
lane3_J98 = subset(lane3, subset = Donor == "J98")
lane3_U114_V105 = subset(lane3, subset = Donor == "U114" 
                         | Donor == "V105")
lane3_M103 = subset(lane3, subset = Donor == "M103")
lane3_J76_WT_Telos = subset(lane3, subset = Donor == "J76" 
                            | Donor == "TELOS")
cells_lane3_J98 = colnames(lane3_J98)
cells_lane3_U114_V105 = colnames(lane3_U114_V105)
cells_lane3_M103 = colnames(lane3_M103)
cells_lane3_J76_WT_Telos = colnames(lane3_J76_WT_Telos)
lane3$Treatment <- ifelse(colnames(lane3) %in% cells_lane3_J98, "siSCR_4d", 
                          ifelse(colnames(lane3) %in% cells_lane3_U114_V105, "ctrl_tx_7d",
                                 ifelse(colnames(lane3) %in% cells_lane3_M103, "TGFb_7d", "il1b_7d")))
lane3$Donor = ifelse(lane3$Donor == "TELOS", "WT_Telo", lane3$Donor)
lane4_J98 = subset(lane4, subset = Donor == "J98")
lane4_D120_NH3 = subset(lane4, subset = Donor == "D120" 
                        | Donor == "TELOS")
lane4_U114_V105 = subset(lane4, subset = Donor == "U114" 
                         | Donor == "V105")
cells_lane4_J98 = colnames(lane4_J98)
cells_lane4_D120_NH3 = colnames(lane4_D120_NH3)
cells_lane4_U114_V105 = colnames(lane4_U114_V105)
lane4$Treatment <- ifelse(colnames(lane4) %in% cells_lane4_J98, "siRNA_4d", 
                          ifelse(colnames(lane4) %in% cells_lane4_D120_NH3, "ctrl_tx_7d",
                                 ifelse(colnames(lane4) %in% cells_lane4_U114_V105, "TGFb_7d", "il1b_4d")))
lane4$Donor = ifelse(lane4$Donor == "TELOS", "NH3", lane4$Donor)
lane5_D120_NH3 = subset(lane5, subset = Donor == "D120" 
                        | Donor == "TELOS")
cells_lane5_D120_NH3 = colnames(lane5_D120_NH3)
lane5$Treatment <- ifelse(colnames(lane5) %in% cells_lane5_D120_NH3, "TGFb_7d", "il1b_7d")
lane5$Donor = ifelse(lane5$Donor == "TELOS", "NH3", lane5$Donor)

# Merge data 
merged = merge(nt, y = c(il1b, lane1, lane2, lane3, lane4, lane5))

# Load genes for regression 
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

# Split RNA object by treatment
list = SplitObject(merged, split.by = "Treatment") 

# Free up memory prior to integration 
rm(nt)
rm(il1b)
rm(lane1)
rm(lane2)
rm(lane3)
rm(lane4)
rm(lane5)

# Generate cell cycle scores for each cell
list = lapply(X = list, FUN = NormalizeData) 
list = lapply(X = list, FUN = CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Regress out cell cycle scores
list <- lapply(X = list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))

# Integrate data 
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 10000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
list <- lapply(X = list, FUN = RunPCA, features = features)
integrated_rna.anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
vitro_by_treatment_more_features_just_rna_cell_cycle_regressed.sct <- IntegrateData(anchorset = integrated_rna.anchors, normalization.method = "SCT", dims = 1:30) 

# Run PCA and find clusters 
integrated_rna <- vitro_by_treatment_more_features_just_rna_cell_cycle_regressed.sct
integrated_rna <- RunPCA(integrated_rna, verbose = FALSE)
integrated_rna <- RunUMAP(integrated_rna, reduction = "pca", dims = 1:30)
integrated_rna <- FindNeighbors(integrated_rna, reduction = "pca", dims = 1:30)
integrated_rna <- FindClusters(integrated_rna, resolution = 0.5)

# Visualize 
DimPlot(integrated_rna, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5)
DimPlot(integrated_rna, reduction = "umap", group.by = "Treatment", label = T, label.size = 3.5, repel = T)
DimPlot(integrated_rna, reduction = "umap", group.by = "Treatment", split.by = "Donor")
DimPlot(integrated_rna, reduction = "umap", group.by = "Donor", label = T, label.size = 3.5, repel = T)
FeaturePlot(integrated_rna, features = c("ERG", "COL1A1"), reduction = "umap")
sct.markers = FindAllMarkers(integrated_rna, assay = "integrated", only.pos = TRUE)

# Explore cell markers 
write.table(sct.markers, "/Users/directory/sct.markers.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Load libraries for ATAC integration 
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(data.table)
library(dittoSeq)
library(DoubletFinder)
library(dplyr)
library(GenomeInfoDb)
library(DoubletFinder)
library(dittoSeq)
library(devtools)
library(glmGamPoi)
set.seed(1234)

# Extract gene annotations from EnsDb, and change to UCSC style
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) = "UCSC"

# Load count data and create seurat objects
counts1 = Read10X_h5("/directory/snt/filtered_feature_bc_matrix.h5")
fragpath1 = "/directory/snt/atac_fragments.tsv.gz"
nt = CreateChromatinAssay(counts = counts1$`Peaks`, sep = c(":", "-"), genome = 'hg38', fragments = fragpath1)
nt = CreateSeuratObject(counts = nt, assay = "ATAC") 
Annotation(object = nt) <- annotations
counts2 = Read10X_h5("/directory/sil1b/filtered_feature_bc_matrix.h5")
fragpath2 = "/directory/sil1b/atac_fragments.tsv.gz"
il1b = CreateChromatinAssay(counts = counts2$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath2)
il1b = CreateSeuratObject(counts = il1b, assay = "ATAC")
Annotation(object = il1b) <- annotations
counts3 = Read10X_h5("/directory/counts/CR74-75_1/filtered_feature_bc_matrix.h5") 
fragpath3 = "/directory/counts/CR74-75_1/atac_fragments.tsv.gz"
lane1 = CreateChromatinAssay(counts = counts3$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath3)
lane1 = CreateSeuratObject(counts = lane1, assay = "ATAC")
Annotation(object = lane1) <- annotations
counts4 = Read10X_h5("/directory/counts/CR74-75_2/filtered_feature_bc_matrix.h5")
fragpath4 = "/directory/counts/CR74-75_2/atac_fragments.tsv.gz"
lane2 = CreateChromatinAssay(counts = counts4$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath4)
lane2 = CreateSeuratObject(counts = lane2, assay = "ATAC")
Annotation(object = lane2) <- annotations
counts5 = Read10X_h5("/directory/counts/CR74-75_3/filtered_feature_bc_matrix.h5")
fragpath5 = "/directory/counts/CR74-75_3/atac_fragments.tsv.gz"
lane3 = CreateChromatinAssay(counts = counts5$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath5)
lane3 = CreateSeuratObject(counts = lane3, assay = "ATAC")
Annotation(object = lane3) <- annotations
counts6 = Read10X_h5("/directory/counts/CR74-75_4/filtered_feature_bc_matrix.h5")
fragpath6 = "/directory/counts/CR74-75_4/atac_fragments.tsv.gz"
lane4 = CreateChromatinAssay(counts = counts6$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath6)
lane4 = CreateSeuratObject(counts = lane4, assay = "ATAC")
Annotation(object = lane4) <- annotations
counts7 = Read10X_h5("/directory/counts/CR74-75_5/filtered_feature_bc_matrix.h5")
fragpath7 = "/directory/counts/CR74-75_5/atac_fragments.tsv.gz"
lane5 = CreateChromatinAssay(counts = counts7$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragpath7)
lane5 = CreateSeuratObject(counts = lane5, assay = "ATAC")
Annotation(object = lane5) <- annotations

# Filter out low quality cells 
nt <- NucleosomeSignal(nt)
nt <- TSSEnrichment(nt, fast = FALSE) 
VlnPlot(object = nt,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
nt <- subset(x = nt, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
il1b <- NucleosomeSignal(il1b)
il1b <- TSSEnrichment(il1b, fast = FALSE) 
VlnPlot(object = il1b, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
il1b <- subset(x = il1b, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
lane1 <- NucleosomeSignal(lane1)
lane1 <- TSSEnrichment(lane1, fast = FALSE) 
VlnPlot(object = lane1,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
lane1 <- subset(x = lane1,subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
lane2 <- NucleosomeSignal(lane2)
lane2 <- TSSEnrichment(lane2, fast = FALSE) 
VlnPlot(object = lane2, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
lane2 <- subset(x = lane2, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
lane3 <- NucleosomeSignal(lane3)
lane3 <- TSSEnrichment(lane3, fast = FALSE) 
VlnPlot(object = lane3,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
lane3 <- subset(x = lane3, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
lane4 <- NucleosomeSignal(lane4)
lane4 <- TSSEnrichment(lane4, fast = FALSE) 
VlnPlot(object = lane4, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
lane4 <- subset(x = lane4, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)
lane5 <- NucleosomeSignal(lane5)
lane5 <- TSSEnrichment(lane5, fast = FALSE) 
VlnPlot(object = lane5, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size = 0)
lane5 <- subset(x = lane5, subset = nucleosome_signal < 2 & TSS.enrichment > 1 & nCount_ATAC > 500)

# Load demuxlet files
nt = importDemux(nt,demuxlet.best = "/directory/snt/snt.demuxlet.concat.best")
nt <- subset(nt, subset = demux.doublet.call == "SNG") 
il1b = importDemux(il1b,demuxlet.best = "/directory/sil1b/sil1b.demuxlet.concat.best")
il1b <- subset(il1b, subset = demux.doublet.call == "SNG") 
lane1 = importDemux(lane1,demuxlet.best = "/directory/CR74-75_1/CR74-75_1.demuxlet.concat.best")
lane1 <- subset(lane1, subset = demux.doublet.call == "SNG") 
lane2 = importDemux(lane2,demuxlet.best = "/directory/CR74-75_2/CR74-75_2.demuxlet.concat.best")
lane2 <- subset(lane2, subset = demux.doublet.call == "SNG") 
lane3 = importDemux(lane3,demuxlet.best = "/directory/CR74-75_3/CR74-75_3.demuxlet.concat.best")
lane3 <- subset(lane3, subset = demux.doublet.call == "SNG") 
lane4 = importDemux(lane4,demuxlet.best = "/directory/CR74-75_4/CR74-75_4.demuxlet.concat.best")
lane4 <- subset(lane4, subset = demux.doublet.call == "SNG") 
lane5 = importDemux(lane5,demuxlet.best = "/directory/CR74-75_5/CR74-75_5.demuxlet.concat.best")
lane5 <- subset(lane5, subset = demux.doublet.call == "SNG") 

# Save lane information 
nt$condition = "nt"
il1b$condition = "il1b"
lane1$condition = "lane1"
lane2$condition = "lane2"
lane3$condition = "lane3"
lane4$condition = "lane4"
lane5$condition = "lane5"

# Make donor IDs meaningful 
x = gsub("\\..*", "", nt$Sample)
nt = AddMetaData(object = nt, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", il1b$Sample)
il1b = AddMetaData(object = il1b, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane1$Sample)
lane1 = AddMetaData(object = lane1, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane2$Sample)
lane2 = AddMetaData(object = lane2, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane3$Sample)
lane3 = AddMetaData(object = lane3, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane4$Sample)
lane4 = AddMetaData(object = lane4, metadata = x, col.name = "Donor")
x = gsub("\\..*", "", lane5$Sample)
lane5 = AddMetaData(object = lane5, metadata = x, col.name = "Donor")

# Add lane IDs to cell barcodes 
x = gsub("^", paste0(nt$condition, "_"), colnames(nt))
nt = RenameCells(nt, new.names = x)
x = gsub("^", paste0(il1b$condition, "_"), colnames(il1b))
il1b = RenameCells(il1b, new.names = x)
x = gsub("^", paste0(lane1$condition, "_"), colnames(lane1))
lane1 = RenameCells(lane1, new.names = x)
x = gsub("^", paste0(lane2$condition, "_"), colnames(lane2))
lane2 = RenameCells(lane2, new.names = x)
x = gsub("^", paste0(lane3$condition, "_"), colnames(lane3))
lane3 = RenameCells(lane3, new.names = x)
x = gsub("^", paste0(lane4$condition, "_"), colnames(lane4))
lane4  = RenameCells(lane4, new.names = x)
x = gsub("^", paste0(lane5$condition, "_"), colnames(lane5))
lane5 = RenameCells(lane5, new.names = x)

# Utilize donor and lane information to assign treatments to each cell barcode in the metadata
nt$Treatment = "ctrl_6hr"
il1b$Treatment = "il1b_6hr"
lane1_V105_J98_U114 = subset(lane1, subset = Donor == "V105" 
                             | Donor == "J98" 
                             | Donor == "U114")
cells_lane1_V105_J98_U114 = colnames(lane1_V105_J98_U114)
lane1$Treatment <- ifelse(colnames(lane1) %in% cells_lane1_V105_J98_U114, "siSCR_7d", "ctrl_tx_7d")
lane1$Donor = ifelse(lane1$Donor == "TELOS", "WT_Telo", lane1$Donor)
lane2_V105_J98_U114 = subset(lane2, subset = Donor == "V105" 
                             | Donor == "J98" 
                             | Donor == "U114")
lane2_J76_WT_Telos = subset(lane2, subset = Donor == "J76" 
                            | Donor == "TELOS")
cells_lane2_V105_J98_U114 = colnames(lane2_V105_J98_U114)
cells_lane2_J76_WT_Telos = colnames(lane2_J76_WT_Telos)
lane2$Treatment <- ifelse(colnames(lane2) %in% cells_lane2_V105_J98_U114, "siRNA_7d", 
                          ifelse(colnames(lane2) %in% cells_lane2_J76_WT_Telos, "TGFb_7d", "TGFb_4d"))
lane2$Donor = ifelse(lane2$Donor == "TELOS", "WT_Telo", lane2$Donor)
lane3_J98 = subset(lane3, subset = Donor == "J98")
lane3_U114_V105 = subset(lane3, subset = Donor == "U114" 
                         | Donor == "V105")
lane3_M103 = subset(lane3, subset = Donor == "M103")
lane3_J76_WT_Telos = subset(lane3, subset = Donor == "J76" 
                            | Donor == "TELOS")
cells_lane3_J98 = colnames(lane3_J98)
cells_lane3_U114_V105 = colnames(lane3_U114_V105)
cells_lane3_M103 = colnames(lane3_M103)
cells_lane3_J76_WT_Telos = colnames(lane3_J76_WT_Telos)
lane3$Treatment <- ifelse(colnames(lane3) %in% cells_lane3_J98, "siSCR_4d", 
                          ifelse(colnames(lane3) %in% cells_lane3_U114_V105, "ctrl_tx_7d",
                                 ifelse(colnames(lane3) %in% cells_lane3_M103, "TGFb_7d", "il1b_7d")))
lane3$Donor = ifelse(lane3$Donor == "TELOS", "WT_Telo", lane3$Donor)
lane4_J98 = subset(lane4, subset = Donor == "J98")
lane4_D120_NH3 = subset(lane4, subset = Donor == "D120" 
                        | Donor == "TELOS")
lane4_U114_V105 = subset(lane4, subset = Donor == "U114" 
                         | Donor == "V105")
cells_lane4_J98 = colnames(lane4_J98)
cells_lane4_D120_NH3 = colnames(lane4_D120_NH3)
cells_lane4_U114_V105 = colnames(lane4_U114_V105)
lane4$Treatment <- ifelse(colnames(lane4) %in% cells_lane4_J98, "siRNA_4d", 
                          ifelse(colnames(lane4) %in% cells_lane4_D120_NH3, "ctrl_tx_7d",
                                 ifelse(colnames(lane4) %in% cells_lane4_U114_V105, "TGFb_7d", "il1b_4d")))
lane4$Donor = ifelse(lane4$Donor == "TELOS", "NH3", lane4$Donor)
lane5_D120_NH3 = subset(lane5, subset = Donor == "D120" 
                        | Donor == "TELOS")
cells_lane5_D120_NH3 = colnames(lane5_D120_NH3)
lane5$Treatment <- ifelse(colnames(lane5) %in% cells_lane5_D120_NH3, "TGFb_7d", "il1b_7d")
lane5$Donor = ifelse(lane5$Donor == "TELOS", "NH3", lane5$Donor)
combined.peaks <- UnifyPeaks(object.list = list(nt, il1b, lane1, lane2, lane3, lane4, lane5), mode = "reduce")

# Quantify common peaks
nt.counts <- FeatureMatrix(Fragments(nt), features = combined.peaks, sep = c(":", "-"), cells = colnames(nt))
il1b.counts <- FeatureMatrix(Fragments(il1b), features = combined.peaks, sep = c(":", "-"), cells = colnames(il1b))
lane1.counts <- FeatureMatrix(Fragments(lane1), features = combined.peaks, sep = c(":", "-"), cells = colnames(lane1))
lane2.counts <- FeatureMatrix(Fragments(lane2), features = combined.peaks, sep = c(":", "-"), cells = colnames(lane2))
lane3.counts <- FeatureMatrix(Fragments(lane3), features = combined.peaks, sep = c(":", "-"), cells = colnames(lane3))
lane4.counts <- FeatureMatrix(Fragments(lane4), features = combined.peaks, sep = c(":", "-"), cells = colnames(lane4))
lane5.counts <- FeatureMatrix(Fragments(lane5), features = combined.peaks, sep = c(":", "-"), cells = colnames(lane5))

# Save metadata 
metadata <- data.frame(c(nt$Treatment, il1b$Treatment, lane1$Treatment, lane2$Treatment, lane3$Treatment, lane4$Treatment, lane5$Treatment), 
                       c(nt$Donor, il1b$Donor, lane1$Donor, lane2$Donor, lane3$Donor, lane4$Donor, lane5$Donor), 
                       c(nt$condition, il1b$condition, lane1$condition, lane2$condition, lane3$condition, lane4$condition, lane5$condition))
colnames(metadata) = c("Treatment", "Donor", "condition")
write.csv(metadata, file = "/directory/metadataTest.csv")

# Overwrite individual peaks files with the common peak set 
nt_assay = CreateChromatinAssay(counts = nt.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(nt))
nt = CreateSeuratObject(counts = nt_assay, assay = "ATAC")
Annotation(object = nt) = annotations
il1b_assay = CreateChromatinAssay(counts = il1b.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(il1b))
il1b = CreateSeuratObject(counts = il1b_assay, assay = "ATAC")
Annotation(object = il1b) = annotations
lane1_assay = CreateChromatinAssay(counts = lane1.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(lane1))
lane1 = CreateSeuratObject(counts = lane1_assay, assay = "ATAC")
Annotation(object = lane1) = annotations
lane2_assay = CreateChromatinAssay(counts = lane2.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(lane2))
lane2 = CreateSeuratObject(counts = lane2_assay, assay = "ATAC")
Annotation(object = lane2) = annotations
lane3_assay = CreateChromatinAssay(counts = lane3.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(lane3))
lane3 = CreateSeuratObject(counts = lane3_assay, assay = "ATAC")
Annotation(object = lane3) = annotations
lane4_assay = CreateChromatinAssay(counts = lane4.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(lane4))
lane4 = CreateSeuratObject(counts = lane4_assay, assay = "ATAC")
Annotation(object = lane4) = annotations
lane5_assay = CreateChromatinAssay(counts = lane5.counts, sep = c(":", "-"), genome = "hg38", fragments = Fragments(lane5))
lane5 = CreateSeuratObject(counts = lane5_assay, assay = "ATAC")
Annotation(object = lane5) = annotations

# Load and add metadata file back to each seurat object 
metadata = read.csv(file = "/directory/metadataTest.csv", row.names = 1)
nt.barcodes = intersect(colnames(nt), rownames(metadata))
nt = AddMetaData(nt, metadata = subset(metadata, selection = rownames(nt) == nt.barcodes))
head(nt@meta.data)
il1b.barcodes = intersect(colnames(il1b), rownames(metadata))
il1b = AddMetaData(il1b, metadata = subset(metadata, selection = rownames(il1b) == il1b.barcodes))
head(il1b@meta.data)
lane1.barcodes = intersect(colnames(lane1), rownames(metadata))
lane1 = AddMetaData(lane1, metadata = subset(metadata, selection = rownames(lane1) == lane1.barcodes))
head(lane1@meta.data)
lane2.barcodes = intersect(colnames(lane2), rownames(metadata))
lane2 = AddMetaData(lane2, metadata = subset(metadata, selection = rownames(lane2) == lane2.barcodes))
head(lane2@meta.data)
lane3.barcodes = intersect(colnames(lane3), rownames(metadata))
lane3 = AddMetaData(lane3, metadata = subset(metadata, selection = rownames(lane3) == lane3.barcodes))
head(lane3@meta.data)
lane4.barcodes = intersect(colnames(lane4), rownames(metadata))
lane4 = AddMetaData(lane4, metadata = subset(metadata, selection = rownames(lane4) == lane4.barcodes))
head(lane4@meta.data)
lane5.barcodes = intersect(colnames(lane5), rownames(metadata))
lane5 = AddMetaData(lane5, metadata = subset(metadata, selection = rownames(lane5) == lane5.barcodes))
head(lane5@meta.data)

# Merge ATAC objects by lane 
merged = merge(nt, y = c(il1b, lane1, lane2, lane3, lane4, lane5))

# Call peaks on every treatment  
peaks.merged <- CallPeaks(merged, macs2.path = "/opt/conda/bin/macs2", group.by = "Treatment", fragment.tempdir = "/directory/", outdir = "/directory/") 
peaks.merged <- keepStandardChromosomes(peaks.merged, pruning.mode = "coarse") # remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks.merged <- subsetByOverlaps(x = peaks.merged, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks = reduce(peaks.merged) # create a unified set of peaks to quantify in each dataset 
macs2_counts.merged <- FeatureMatrix(fragments = Fragments(merged), features = combined.peaks, cells = colnames(merged)) # quantify counts in each peak
merged_assay = CreateChromatinAssay(counts = macs2_counts.merged, sep = c(":", "-"), genome = "hg38", fragments = Fragments(merged))

# Create a new seurat object with new peak set and add annotations and metadata back to the new object
merged = CreateSeuratObject(counts = merged_assay, assay = "ATAC")
Annotation(object = merged) = annotations 
merged.barcodes = intersect(colnames(merged), rownames(metadata))
merged = AddMetaData(merged, metadata = subset(metadata, selection = rownames(merged) == merged.barcodes))
head(merged@meta.data)

# Compute LSI based on new peak coordinates
merged <- RunTFIDF(merged)
merged<- FindTopFeatures(merged)
merged <- RunSVD(merged)
merged <- RunUMAP(merged, reduction = 'lsi', dims = 2:30)

# Integrate by treatment 
object.list <- SplitObject(merged, split.by = "Treatment") 
integration.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = UnifyPeaks(object.list, mode = "reduce"), reduction = "rlsi", dims = 2:30) 
integrated_atac <- IntegrateEmbeddings(anchorset = integration.anchors, reductions = merged[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 1:30, k.weight = 50)
Fragments(integrated_atac) 

# Visualize new object 
integrated_atac <- RunUMAP(integrated_atac, reduction = "integrated_lsi", dims = 2:30)
integrated_atac <- FindNeighbors(integrated_atac, reduction = 'integrated_lsi', dims = 2:30)
integrated_atac <- FindClusters(integrated_atac, algorithm = 3, resolution = 0.38, verbose = FALSE)
DimPlot(integrated_atac, reduction = "umap")
DimPlot(integrated_atac, reduction = "umap", group.by = "Treatment")
DimPlot(integrated_atac, reduction = "umap", group.by = "Treatment", split.by = "Donor")

# Match cell barcodes between RNA and ATAC objects 

# Load integrated RNA object and run PCA
integrated_rna = readRDS("/directory/integrated_rna_sct_with_multiome.RDS")
integrated_rna <- RunPCA(integrated_rna, verbose = FALSE)
integrated_rna <- RunUMAP(integrated_rna, reduction = "pca", dims = 1:30)
integrated_rna <- FindNeighbors(integrated_rna, reduction = "pca", dims = 1:30)
integrated_rna <- FindClusters(integrated_rna, resolution = 0.5)
DimPlot(integrated_rna, reduction = "umap", label = T)

# Filter out cell barcodes that do not match
integrated_atac_filtered = subset(integrated_atac, cells = WhichCells(integrated_rna))
integrated_rna_filtered = subset(integrated_rna, cells = WhichCells(integrated_atac_filtered))

# Combine datasets
integrated_rna_filtered[["ATAC"]] = integrated_atac_filtered[["ATAC"]]
integrated_rna_filtered[["integrated_lsi"]] = integrated_atac_filtered[["integrated_lsi"]] 
integrated_multi = integrated_rna_filtered

# Perform WNN analysis 
integrated_multi <- FindMultiModalNeighbors(
  object = integrated_multi,
  reduction.list = list("pca", "integrated_lsi"), 
  dims.list = list(1:30, 2:30), 
  verbose = TRUE
)

# Build joint UMAP for visualization and find clusters 
integrated_multi <- RunUMAP(
  object = integrated_multi,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)
integrated_multi <- FindClusters(integrated_multi, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE) 

# Call peaks based on the new clusters 
atac = integrated_multi 
DefaultAssay(atac) = "ATAC"
atac.slim = DietSeurat(
  object = atac,
  counts = TRUE,
  data = TRUE,
  scale.data = TRUE,
  features = NULL,
  assays = c("ATAC"),
  dimreducs = c("wnn.umap"), 
  graphs = NULL
)
peaks.atac <- CallPeaks(atac.slim, macs2.path = "/opt/conda/bin/macs2", group.by = "seurat_clusters", fragment.tempdir = "/directory/", outdir = "/directory/") 
peaks.atac <- keepStandardChromosomes(peaks.atac, pruning.mode = "coarse") # remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks.atac <- subsetByOverlaps(x = peaks.atac, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks = reduce(peaks.atac) # Create a unified set of peaks to quantify in each dataset 
macs2_counts.atac <- FeatureMatrix(fragments = Fragments(atac.slim), features = combined.peaks, cells = colnames(atac.slim)) # quantify counts in each peak

# Create a new seurat object with new peaks and add annotations and metadata back 
atac.slim_assay = CreateChromatinAssay(counts = macs2_counts.atac, sep = c(":", "-"), genome = "hg38", fragments = Fragments(atac.slim))
atac.slim = CreateSeuratObject(counts = atac.slim_assay, assay = "ATAC")
Annotation(object = atac.slim) = annotations # Add annotations back 
atac.slim.barcodes = intersect(colnames(atac.slim), rownames(metadata)) # Add metadata back 
atac.slim = AddMetaData(atac.slim, metadata = subset(metadata, selection = rownames(atac.slim) == atac.slim.barcodes))
head(atac.slim@meta.data)

# Compute LSI and reintegrate ATAC data based on new LSI
atac.slim <- RunTFIDF(atac.slim)
atac.slim <- FindTopFeatures(atac.slim)
atac.slim <- RunSVD(atac.slim)
atac.slim <- RunUMAP(atac.slim, reduction = 'lsi', dims = 2:30)

# Re-integrate ATAC dataset by treatment
object.list <- SplitObject(atac.slim, split.by = "Treatment") 
atac.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = UnifyPeaks(object.list, mode = "reduce"), reduction = "rlsi", dims = 2:30) # Find integration anchors
atac_new <- IntegrateEmbeddings(anchorset = atac.anchors, reductions = atac.slim[["lsi"]], new.reduction.name = "new_integrated_lsi", dims.to.integrate = 1:30, k.weight = 50) 

# Overwrite ATAC multiome data with new peak set and new integrated lsi
integrated_multi[["ATAC"]] = atac_new[["ATAC"]] 
integrated_multi[["new_integrated_lsi"]] = atac_new[["new_integrated_lsi"]] 

# Remove fragment files with no cells to speed up downstream analyses
DefaultAssay(integrated_multi) = "ATAC"
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
}
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
}
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
} 
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
}  
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
} 
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
} 
for (i in 1:length(Fragments(integrated_multi))){ 
  if(length(integrated_multi@assays[["ATAC"]]@fragments[[i]]@cells)==0){integrated_multi@assays[["ATAC"]]@fragments[[i]]=NULL}
} 

# Perform final round of WNN analysis 
integrated_multi <- FindMultiModalNeighbors(
  object = integrated_multi,
  reduction.list = list("pca", "new_integrated_lsi"), 
  dims.list = list(1:30, 2:30), 
  verbose = TRUE
)

# Build joint UMAP visualization
integrated_multi <- RunUMAP(
  object = integrated_multi,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)
integrated_multi <- FindClusters(integrated_multi, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE) 
DimPlot(integrated_multi, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)

# Compute the GC content for each peak using RegionStats
integrated_multi = RegionStats(integrated_multi, genome = BSgenome.Hsapiens.UCSC.hg38)

# Link peaks to genes
final_multi = LinkPeaks(object = integrated_multi, peak.assay = "ATAC", expression.assay = "integrated") 


# Integrate public ex vivo data 

# Load necessary libraries
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(hdf5r) 
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86) 
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ggplot2)
library(data.table)
library(dittoSeq) 
library(DoubletFinder)
library(dplyr)
library(GenomeInfoDb)
library(devtools)
library(glmGamPoi)
set.seed(1234)

# Load individual ex vivo samples processed from GEO 
rca11 = Read10X(data.dir = "/directory/GSE131778/rca11/filtered_feature_bc_matrix/")
rca11 = CreateSeuratObject(rca11)
rca11[["percent.mt"]] = PercentageFeatureSet(rca11, pattern = "^MT-")
rca11 <- subset(rca11, subset = nFeature_RNA > 500 & percent.mt < 20)
rca11 = NormalizeData(rca11)
rca11 = FindVariableFeatures(rca11, selection.method = "vst", nfeatures = 3000)
rca11.all.genes = rownames(rca11)
rca11 = ScaleData(rca11, features = rca11.all.genes)
rca11 = RunPCA(rca11, features = VariableFeatures(object = rca11))
rca11 = RunUMAP(rca11, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca11, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca11 <- FindNeighbors(rca11, reduction = 'pca')
rca11 <- FindClusters(rca11)
nExp_poi <- round(0.20*length(colnames(rca11))) # remove 20% (of all) doublets
rca11 <- doubletFinder_v3(rca11, PCs = 1:20, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca11@meta.data)[grepl("DF.classification", colnames(rca11@meta.data))]
rca11 = rca11[, rca11@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca11, file = "/directory/rca11.RDS")
rca12 = Read10X(data.dir = "/directory/GSE131778/rca12/filtered_feature_bc_matrix/")
rca12 = CreateSeuratObject(rca12)
rca12[["percent.mt"]] = PercentageFeatureSet(rca12, pattern = "^MT-")
rca12 <- subset(rca12, subset = nFeature_RNA > 500 & percent.mt < 20)
rca12 = NormalizeData(rca12)
rca12 = FindVariableFeatures(rca12, selection.method = "vst", nfeatures = 3000)
rca12.all.genes = rownames(rca12)
rca12 = ScaleData(rca12, features = rca12.all.genes)
rca12 = RunPCA(rca12, features = VariableFeatures(object = rca12))
rca12 = RunUMAP(rca12, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca12, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats)
rca12 <- FindNeighbors(rca12, reduction = 'pca')
rca12 <- FindClusters(rca12)
nExp_poi <- round(0.20*length(colnames(rca12))) 
rca12 <- doubletFinder_v3(rca12, PCs = 1:20, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca12@meta.data)[grepl("DF.classification", colnames(rca12@meta.data))]
rca12 = rca12[, rca12@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca12, file = "/directory/rca12.RDS")
rca21 = Read10X(data.dir = "/directory/GSE131778/rca21/filtered_feature_bc_matrix/")
rca21 = CreateSeuratObject(rca21)
rca21[["percent.mt"]] = PercentageFeatureSet(rca21, pattern = "^MT-")
rca21 <- subset(rca21, subset = nFeature_RNA > 500 & percent.mt < 20)
rca21 = NormalizeData(rca21)
rca21 = FindVariableFeatures(rca21, selection.method = "vst", nfeatures = 3000)
rca21.all.genes = rownames(rca21)
rca21 = ScaleData(rca21, features = rca21.all.genes)
rca21 = RunPCA(rca21, features = VariableFeatures(object = rca21))
rca21 = RunUMAP(rca21, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca21, PCs = 1:20, sct = F) # pK identification
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) # choose smallest pk
rca21 <- FindNeighbors(rca21, reduction = 'pca')
rca21 <- FindClusters(rca21)
nExp_poi <- round(0.20*length(colnames(rca21))) # remove 20% (of all) doublets
rca21 <- doubletFinder_v3(rca21, PCs = 1:20, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca21@meta.data)[grepl("DF.classification", colnames(rca21@meta.data))]
rca21 = rca21[, rca21@meta.data[, DF.name] == "Singlet"] # subset singlets
saveRDS(rca21, file = "/directory/rca21.RDS")
rca22 = Read10X(data.dir = "/xdisk/cromanoski/madelus/paper/GSE131778/rca22/filtered_feature_bc_matrix/")
rca22 = CreateSeuratObject(rca22)
rca22[["percent.mt"]] = PercentageFeatureSet(rca22, pattern = "^MT-")
rca22 <- subset(rca22, subset = nFeature_RNA > 500 & percent.mt < 20)
rca22 = NormalizeData(rca22)
rca22 = FindVariableFeatures(rca22, selection.method = "vst", nfeatures = 3000)
rca22.all.genes = rownames(rca22)
rca22 = ScaleData(rca22, features = rca22.all.genes)
rca22 = RunPCA(rca22, features = VariableFeatures(object = rca22))
rca22 = RunUMAP(rca22, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca22, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca22 <- FindNeighbors(rca22, reduction = 'pca')
rca22 <- FindClusters(rca22)
nExp_poi <- round(0.20*length(colnames(rca22))) 
rca22 <- doubletFinder_v3(rca22, PCs = 1:20, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca22@meta.data)[grepl("DF.classification", colnames(rca22@meta.data))]
rca22 = rca22[, rca22@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca22, file = "/directory/rca22.RDS")
rca31 = Read10X(data.dir = "/directory/GSE131778/rca31/filtered_feature_bc_matrix/")
rca31 = CreateSeuratObject(rca31)
rca31[["percent.mt"]] = PercentageFeatureSet(rca31, pattern = "^MT-")
rca31 <- subset(rca31, subset = nFeature_RNA > 500 & percent.mt < 20)
rca31 = NormalizeData(rca31)
rca31 = FindVariableFeatures(rca31, selection.method = "vst", nfeatures = 3000)
rca31.all.genes = rownames(rca31)
rca31 = ScaleData(rca31, features = rca31.all.genes)
rca31 = RunPCA(rca31, features = VariableFeatures(object = rca31))
rca31 = RunUMAP(rca31, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca31, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca31 <- FindNeighbors(rca31, reduction = 'pca')
rca31 <- FindClusters(rca31)
nExp_poi <- round(0.20*length(colnames(rca31))) 
rca31 <- doubletFinder_v3(rca31, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca31@meta.data)[grepl("DF.classification", colnames(rca31@meta.data))]
rca31 = rca31[, rca31@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca31, file = "/directory/rca31.RDS")
rca32 = Read10X(data.dir = "/directory/GSE131778/rca32/filtered_feature_bc_matrix/")
rca32 = CreateSeuratObject(rca32)
rca32[["percent.mt"]] = PercentageFeatureSet(rca32, pattern = "^MT-")
rca32 <- subset(rca32, subset = nFeature_RNA > 500 & percent.mt < 20)
rca32 = NormalizeData(rca32)
rca32 = FindVariableFeatures(rca32, selection.method = "vst", nfeatures = 3000)
rca32.all.genes = rownames(rca32)
rca32 = ScaleData(rca32, features = rca32.all.genes)
rca32 = RunPCA(rca32, features = VariableFeatures(object = rca32))
rca32 = RunUMAP(rca32, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca32, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca32 <- FindNeighbors(rca32, reduction = 'pca')
rca32 <- FindClusters(rca32)
nExp_poi <- round(0.20*length(colnames(rca32))) 
rca32 <- doubletFinder_v3(rca32, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca32@meta.data)[grepl("DF.classification", colnames(rca32@meta.data))]
rca32 = rca32[, rca32@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca32, file = "/directory/rca32.RDS")
rca33 = Read10X(data.dir = "/directory/GSE131778/rca33/filtered_feature_bc_matrix/")
rca33 = CreateSeuratObject(rca33)
rca33[["percent.mt"]] = PercentageFeatureSet(rca33, pattern = "^MT-")
rca33 <- subset(rca33, subset = nFeature_RNA > 500 & percent.mt < 20)
rca33 = NormalizeData(rca33)
rca33 = FindVariableFeatures(rca33, selection.method = "vst", nfeatures = 3000)
rca33.all.genes = rownames(rca33)
rca33 = ScaleData(rca33, features = rca33.all.genes)
rca33 = RunPCA(rca33, features = VariableFeatures(object = rca33))
rca33 = RunUMAP(rca33, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca33, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca33 <- FindNeighbors(rca33, reduction = 'pca')
rca33 <- FindClusters(rca33)
nExp_poi <- round(0.20*length(colnames(rca33))) 
rca33 <- doubletFinder_v3(rca33, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca33@meta.data)[grepl("DF.classification", colnames(rca33@meta.data))]
rca33 = rca33[, rca33@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca33, file = "/directory/rca33.RDS")
rca41 = Read10X(data.dir = "/directory/GSE131778/rca41/filtered_feature_bc_matrix/")
rca41 = CreateSeuratObject(rca41)
rca41[["percent.mt"]] = PercentageFeatureSet(rca41, pattern = "^MT-")
rca41 <- subset(rca41, subset = nFeature_RNA > 500 & percent.mt < 20)
rca41 = NormalizeData(rca41)
rca41 = FindVariableFeatures(rca41, selection.method = "vst", nfeatures = 3000)
rca41.all.genes = rownames(rca41)
rca41 = ScaleData(rca41, features = rca41.all.genes)
rca41 = RunPCA(rca41, features = VariableFeatures(object = rca41))
rca41 = RunUMAP(rca41, dims = 1:20)
sweep.res.list <- paramSweep_v3(rca41, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
rca41 <- FindNeighbors(rca41, reduction = 'pca')
rca41 <- FindClusters(rca41)
nExp_poi <- round(0.20*length(colnames(rca41))) 
rca41 <- doubletFinder_v3(rca41, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(rca41@meta.data)[grepl("DF.classification", colnames(rca41@meta.data))]
rca41 = rca41[, rca41@meta.data[, DF.name] == "Singlet"] 
saveRDS(rca41, file = "/directory/rca41.RDS")
patient1ac = Read10X(data.dir = "/directory/GSE159677/patient1ac/filtered_feature_bc_matrix/")
patient1ac = CreateSeuratObject(patient1ac)
patient1ac[["percent.mt"]] = PercentageFeatureSet(patient1ac, pattern = "^MT-")
patient1ac <- subset(patient1ac, subset = nFeature_RNA > 500 & percent.mt < 20)
patient1ac = NormalizeData(patient1ac)
patient1ac = FindVariableFeatures(patient1ac, selection.method = "vst", nfeatures = 3000)
patient1ac.all.genes = rownames(patient1ac)
patient1ac = ScaleData(patient1ac, features = patient1ac.all.genes)
patient1ac = RunPCA(patient1ac, features = VariableFeatures(object = patient1ac))
patient1ac = RunUMAP(patient1ac, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient1ac, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient1ac <- FindNeighbors(patient1ac, reduction = 'pca')
patient1ac <- FindClusters(patient1ac)
nExp_poi <- round(0.20*length(colnames(patient1ac))) 
patient1ac <- doubletFinder_v3(patient1ac, PCs = 1:20, pN = 0.25, pK = 0.001, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient1ac@meta.data)[grepl("DF.classification", colnames(patient1ac@meta.data))]
patient1ac = patient1ac[, patient1ac@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient1ac, file = "/directory/patient1ac.RDS")
patient1pa = Read10X(data.dir = "/directory/GSE159677/patient1pa/filtered_feature_bc_matrix/")
patient1pa = CreateSeuratObject(patient1pa)
patient1pa[["percent.mt"]] = PercentageFeatureSet(patient1pa, pattern = "^MT-")
patient1pa <- subset(patient1pa, subset = nFeature_RNA > 500 & percent.mt < 20)
patient1pa = NormalizeData(patient1pa)
patient1pa = FindVariableFeatures(patient1pa, selection.method = "vst", nfeatures = 3000)
patient1pa.all.genes = rownames(patient1pa)
patient1pa = ScaleData(patient1pa, features = patient1pa.all.genes)
patient1pa = RunPCA(patient1pa, features = VariableFeatures(object = patient1pa))
patient1pa = RunUMAP(patient1pa, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient1pa, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient1pa <- FindNeighbors(patient1pa, reduction = 'pca')
patient1pa <- FindClusters(patient1pa)
nExp_poi <- round(0.20*length(colnames(patient1pa))) 
patient1pa <- doubletFinder_v3(patient1pa, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient1pa@meta.data)[grepl("DF.classification", colnames(patient1pa@meta.data))]
patient1pa = patient1pa[, patient1pa@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient1pa, file = "/directory/patient1pa.RDS")
patient2ac = Read10X(data.dir = "/directory/GSE159677/patient2ac/filtered_feature_bc_matrix/")
patient2ac = CreateSeuratObject(patient2ac)
patient2ac[["percent.mt"]] = PercentageFeatureSet(patient2ac, pattern = "^MT-")
patient2ac <- subset(patient2ac, subset = nFeature_RNA > 500 & percent.mt < 20)
patient2ac = NormalizeData(patient2ac)
patient2ac = FindVariableFeatures(patient2ac, selection.method = "vst", nfeatures = 3000)
patient2ac.all.genes = rownames(patient2ac)
patient2ac = ScaleData(patient2ac, features = patient2ac.all.genes)
patient2ac = RunPCA(patient2ac, features = VariableFeatures(object = patient2ac))
patient2ac = RunUMAP(patient2ac, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient2ac, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient2ac <- FindNeighbors(patient2ac, reduction = 'pca')
patient2ac <- FindClusters(patient2ac)
nExp_poi <- round(0.20*length(colnames(patient2ac)))
patient2ac <- doubletFinder_v3(patient2ac, PCs = 1:20, pN = 0.25, pK = 0.001, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient2ac@meta.data)[grepl("DF.classification", colnames(patient2ac@meta.data))]
patient2ac = patient2ac[, patient2ac@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient2ac, file = "/directory/patient2ac.RDS")
patient2pa = Read10X(data.dir = "/directory/GSE159677/patient2pa/filtered_feature_bc_matrix/")
patient2pa = CreateSeuratObject(patient2pa)
patient2pa[["percent.mt"]] = PercentageFeatureSet(patient2pa, pattern = "^MT-")
patient2pa <- subset(patient2pa, subset = nFeature_RNA > 500 & percent.mt < 20)
patient2pa = NormalizeData(patient2pa)
patient2pa = FindVariableFeatures(patient2pa, selection.method = "vst", nfeatures = 3000)
patient2pa.all.genes = rownames(patient2pa)
patient2pa = ScaleData(patient2pa, features = patient2pa.all.genes)
patient2pa = RunPCA(patient2pa, features = VariableFeatures(object = patient2pa))
patient2pa = RunUMAP(patient2pa, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient2pa, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient2pa <- FindNeighbors(patient2pa, reduction = 'pca')
patient2pa <- FindClusters(patient2pa)
nExp_poi <- round(0.20*length(colnames(patient2pa))) 
patient2pa <- doubletFinder_v3(patient2pa, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient2pa@meta.data)[grepl("DF.classification", colnames(patient2pa@meta.data))]
patient2pa = patient2pa[, patient2pa@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient2pa, file = "/directory/patient2pa.RDS")
patient3ac = Read10X(data.dir = "/directory/GSE159677/patient3ac/filtered_feature_bc_matrix/")
patient3ac = CreateSeuratObject(patient3ac)
patient3ac[["percent.mt"]] = PercentageFeatureSet(patient3ac, pattern = "^MT-")
patient3ac <- subset(patient3ac, subset = nFeature_RNA > 500 & percent.mt < 20)
patient3ac = NormalizeData(patient3ac)
patient3ac = FindVariableFeatures(patient3ac, selection.method = "vst", nfeatures = 3000)
patient3ac.all.genes = rownames(patient3ac)
patient3ac = ScaleData(patient3ac, features = patient3ac.all.genes)
patient3ac = RunPCA(patient3ac, features = VariableFeatures(object = patient3ac))
patient3ac = RunUMAP(patient3ac, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient3ac, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient3ac <- FindNeighbors(patient3ac, reduction = 'pca')
patient3ac <- FindClusters(patient3ac)
nExp_poi <- round(0.20*length(colnames(patient3ac))) 
patient3ac <- doubletFinder_v3(patient3ac, PCs = 1:20, pN = 0.25, pK = 0.001, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient3ac@meta.data)[grepl("DF.classification", colnames(patient3ac@meta.data))]
patient3ac = patient3ac[, patient3ac@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient3ac, file = "/directory/patient3ac.RDS")
patient3pa = Read10X(data.dir = "/directory/GSE159677/patient3pa/filtered_feature_bc_matrix/")
patient3pa = CreateSeuratObject(patient3pa)
patient3pa[["percent.mt"]] = PercentageFeatureSet(patient3pa, pattern = "^MT-")
patient3pa <- subset(patient3pa, subset = nFeature_RNA > 500 & percent.mt < 20)
patient3pa = NormalizeData(patient3pa)
patient3pa = FindVariableFeatures(patient3pa, selection.method = "vst", nfeatures = 3000)
patient3pa.all.genes = rownames(patient3pa)
patient3pa = ScaleData(patient3pa, features = patient3pa.all.genes)
patient3pa = RunPCA(patient3pa, features = VariableFeatures(object = patient3pa))
patient3pa = RunUMAP(patient3pa, dims = 1:20)
sweep.res.list <- paramSweep_v3(patient3pa, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
patient3pa <- FindNeighbors(patient3pa, reduction = 'pca')
patient3pa <- FindClusters(patient3pa)
nExp_poi <- round(0.20*length(colnames(patient3pa))) 
patient3pa <- doubletFinder_v3(patient3pa, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(patient3pa@meta.data)[grepl("DF.classification", colnames(patient3pa@meta.data))]
patient3pa = patient3pa[, patient3pa@meta.data[, DF.name] == "Singlet"] 
saveRDS(patient3pa, file = "/directory/patient3pa.RDS")
carotid1 = Read10X(data.dir = "/directory/GSE155512/carotid1/filtered_feature_bc_matrix/")
carotid1 = CreateSeuratObject(carotid1)
carotid1[["percent.mt"]] = PercentageFeatureSet(carotid1, pattern = "^MT-")
carotid1 <- subset(carotid1, subset = nFeature_RNA > 500 & percent.mt < 20)
carotid1 = NormalizeData(carotid1)
carotid1 = FindVariableFeatures(carotid1, selection.method = "vst", nfeatures = 3000)
carotid1.all.genes = rownames(carotid1)
carotid1 = ScaleData(carotid1, features = carotid1.all.genes)
carotid1 = RunPCA(carotid1, features = VariableFeatures(object = carotid1))
carotid1 = RunUMAP(carotid1, dims = 1:20)
sweep.res.list <- paramSweep_v3(carotid1, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
carotid1 <- FindNeighbors(carotid1, reduction = 'pca')
carotid1 <- FindClusters(carotid1)
nExp_poi <- round(0.20*length(colnames(carotid1))) 
carotid1 <- doubletFinder_v3(carotid1, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(carotid1@meta.data)[grepl("DF.classification", colnames(carotid1@meta.data))]
carotid1 = carotid1[, carotid1@meta.data[, DF.name] == "Singlet"] 
saveRDS(carotid1, file = "/directory/carotid1.RDS")
carotid2 = Read10X(data.dir = "/directory/GSE155512/carotid2/filtered_feature_bc_matrix/")
carotid2 = CreateSeuratObject(carotid2)
carotid2[["percent.mt"]] = PercentageFeatureSet(carotid2, pattern = "^MT-")
carotid2 <- subset(carotid2, subset = nFeature_RNA > 500 & percent.mt < 20)
carotid2 = NormalizeData(carotid2)
carotid2 = FindVariableFeatures(carotid2, selection.method = "vst", nfeatures = 3000)
carotid2.all.genes = rownames(carotid2)
carotid2 = ScaleData(carotid2, features = carotid2.all.genes)
carotid2 = RunPCA(carotid2, features = VariableFeatures(object = carotid2))
carotid2 = RunUMAP(carotid2, dims = 1:20)
sweep.res.list <- paramSweep_v3(carotid2, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
carotid2 <- FindNeighbors(carotid2, reduction = 'pca')
carotid2 <- FindClusters(carotid2)
nExp_poi <- round(0.20*length(colnames(carotid2))) 
carotid2 <- doubletFinder_v3(carotid2, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(carotid2@meta.data)[grepl("DF.classification", colnames(carotid2@meta.data))]
carotid2 = carotid2[, carotid2@meta.data[, DF.name] == "Singlet"] 
saveRDS(carotid2, file = "/directory/carotid2.RDS")
carotid3 = Read10X(data.dir = "/directory/GSE155512/carotid3/filtered_feature_bc_matrix/")
carotid3 = CreateSeuratObject(carotid3)
carotid3[["percent.mt"]] = PercentageFeatureSet(carotid3, pattern = "^MT-")
carotid3 <- subset(carotid3, subset = nFeature_RNA > 500 & percent.mt < 20)
carotid3 = NormalizeData(carotid3)
carotid3 = FindVariableFeatures(carotid3, selection.method = "vst", nfeatures = 3000)
carotid3.all.genes = rownames(carotid3)
carotid3 = ScaleData(carotid3, features = carotid3.all.genes)
carotid3 = RunPCA(carotid3, features = VariableFeatures(object = carotid3))
carotid3 = RunUMAP(carotid3, dims = 1:20)
sweep.res.list <- paramSweep_v3(carotid3, PCs = 1:20, sct = F) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats) 
carotid3 <- FindNeighbors(carotid3, reduction = 'pca')
carotid3 <- FindClusters(carotid3)
nExp_poi <- round(0.20*length(colnames(carotid3)))
carotid3 <- doubletFinder_v3(carotid3, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(carotid3@meta.data)[grepl("DF.classification", colnames(carotid3@meta.data))]
carotid3 = carotid3[, carotid3@meta.data[, DF.name] == "Singlet"] 
saveRDS(carotid3, file = "/directory/carotid3.RDS")

# save samples before integrating
rca11$sample = "rca11"
rca12$sample = "rca12"
rca21$sample = "rca21"
rca22$sample = "rca22"
rca31$sample = "rca31"
rca32$sample = "rca32"
rca33$sample = "rca33"
rca41$sample = "rca41"
patient1ac$sample = "patient1ac"
patient1pa$sample = "patient1pa"
patient2ac$sample = "patient2ac"
patient2pa$sample = "patient2pa"
patient3ac$sample = "patient3ac"
patient3pa$sample = "patient3pa"
carotid1$sample = "carotid1"
carotid2$sample = "carotid2"
carotid3$sample = "carotid3"

# add sample IDs to barcodes 
x = gsub("^", paste0(rca11$sample, "_"), colnames(rca11))
rca11 = RenameCells(rca11, new.names = x)
x = gsub("^", paste0(rca12$sample, "_"), colnames(rca12))
rca12 = RenameCells(rca12, new.names = x)
x = gsub("^", paste0(rca21$sample, "_"), colnames(rca21))
rca21 = RenameCells(rca21, new.names = x)
x = gsub("^", paste0(rca22$sample, "_"), colnames(rca22))
rca22 = RenameCells(rca22, new.names = x)
x = gsub("^", paste0(rca31$sample, "_"), colnames(rca31))
rca31 = RenameCells(rca31, new.names = x)
x = gsub("^", paste0(rca32$sample, "_"), colnames(rca32))
rca32 = RenameCells(rca32, new.names = x)
x = gsub("^", paste0(rca33$sample, "_"), colnames(rca33))
rca33 = RenameCells(rca33, new.names = x)
x = gsub("^", paste0(rca41$sample, "_"), colnames(rca41))
rca41 = RenameCells(rca41, new.names = x)
x = gsub("^", paste0(patient1ac$sample, "_"), colnames(patient1ac))
patient1ac = RenameCells(patient1ac, new.names = x)
x = gsub("^", paste0(patient1pa$sample, "_"), colnames(patient1pa))
patient1pa = RenameCells(patient1pa, new.names = x)
x = gsub("^", paste0(patient2ac$sample, "_"), colnames(patient2ac))
patient2ac = RenameCells(patient2ac, new.names = x)
x = gsub("^", paste0(patient2pa$sample, "_"), colnames(patient2pa))
patient2pa = RenameCells(patient2pa, new.names = x)
x = gsub("^", paste0(patient3ac$sample, "_"), colnames(patient3ac))
patient3ac = RenameCells(patient3ac, new.names = x)
x = gsub("^", paste0(patient3pa$sample, "_"), colnames(patient3pa))
patient3pa = RenameCells(patient3pa, new.names = x)
x = gsub("^", paste0(carotid1$sample, "_"), colnames(carotid1))
carotid1 = RenameCells(carotid1, new.names = x)
x = gsub("^", paste0(carotid2$sample, "_"), colnames(carotid2))
carotid2 = RenameCells(carotid2, new.names = x)
x = gsub("^", paste0(carotid3$sample, "_"), colnames(carotid3))
carotid3 = RenameCells(carotid3, new.names = x)

# merge datasets
merged = merge(rca11, y = c(rca12, rca21, rca22, rca31, rca32, rca33, rca41,
                            patient1ac, patient1pa, patient2ac, patient2pa, patient3ac, patient3pa,
                            carotid1, carotid2, carotid3))

# cell cycle regression and integration 
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

# integrate samples without doublets 
list = SplitObject(merged, split.by = "sample") 
list = lapply(X = list, FUN = CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
list = lapply(X = list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000) 
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
list <- lapply(X = list, FUN = RunPCA, features = features)
integrated_rna.anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
vivo <- IntegrateData(anchorset = integrated_rna.anchors, normalization.method = "SCT", dims = 1:30) 
vivo <- RunPCA(vivo, verbose = FALSE)
vivo <- RunUMAP(vivo, reduction = "pca", dims = 1:30)
vivo <- FindNeighbors(vivo, reduction = "pca", dims = 1:30)
vivo <- FindClusters(vivo, resolution = 0.6)
saveRDS(vivo, file = "/directory/vivo.RDS") 
