# 1 CellRanger Cloud
Cell Ranger Count v6.0.0
Settings:    chemistry: auto, recovered_cells: 10000, include_introns: false, no_bam: true, no_secondary_analysis: false
Reference:    Mouse (mm10) 2020-A

# 2 Input 10X data into lists
library(Seurat)
Samples = list.files(path='./outs/') # outputs from CellRanger
Sample_Meta <- read.csv("../metadata_samples.csv")
Cell_Meta <- read.csv("../metadata_cell.csv")

Samples_RawDat = vector(mode='list',length=length(Samples))
names(Samples_RawDat) = Samples

for(i in Samples){
    Samples_RawDat[[i]] <- Read10X(data.dir = paste("/home-3/lhayes14@jhu.edu/work/lindsay/10X/outs",i,"mats",sep='/'))
    cat(paste(i,', ',sep=''))
}

# 3 Make list of Seurat Objects
Samples_AllDat = vector(mode='list', length=length(Samples))
names(Samples_AllDat) = Samples

for(i in Samples){
    Samples_AllDat[[i]] <- CreateSeuratObject(counts = Samples_RawDat[[i]])
    CellTracking[i,'Raw'] = ncol(Samples_AllDat[[i]])
    Samples_AllDat[[i]]@meta.data$sample = Sample_Meta$id[which(Samples == i)]
    Samples_AllDat[[i]] = RenameCells(Samples_AllDat[[i]], new.names=paste(Sample_Meta$id[which(Samples ==i)], colnames(Samples_AllDat[[i]]), sep='_'))
    cat(paste(i,', ',sep=''))
}

# 4 Add metadata to cells
for(i in Samples){
    MetaInd = match(colnames(Samples_AllDat[[i]]), as.character(Cell_Meta$cellID))
    Samples_AllDat[[i]]@meta.data$group = Cell_Meta$group[MetaInd]
    Samples_AllDat[[i]]@meta.data$core = Cell_Meta$core[MetaInd]
    Samples_AllDat[[i]]@meta.data$IL6 = Cell_Meta$IL6[MetaInd]
    Samples_AllDat[[i]]@meta.data$date = Cell_Meta$date[MetaInd]
}

# 5 Estimate Doublets
library(scDblFinder)
library(Seurat)
set.seed(100)
for(i in Samples){
TMPsce = Seurat::as.SingleCellExperiment(Samples_AllDat[[i]])
TMPsce <- scDblFinder(TMPsce)
Samples_AllDat[[i]]@meta.data$scDblFinder = TMPsce$scDblFinder.class
cat(paste(i,', ',sep=''))
}

# 6 Filter Cells
Samples <- names(Samples_AllDat)
# estimate filters
means <- data.frame(sample = Samples, nCount_RNA = 0, nFeature_RNA=0)
rownames(means) = Samples
for(i in c(1:length(Samples))){
    means$nCount_RNA[i] = mean(Samples_AllDat[[i]]@meta.data$nCount_RNA)
    means$nFeature_RNA[i] = mean(Samples_AllDat[[i]]@meta.data$nFeature_RNA)
    means$nFeature_RNA[i] = mean(Samples_AllDat[[i]]@meta.data$percent)
}

# Filter
CellTracking = data.frame(sample = Samples, Filt_v1 = 0)
FiltDat_v1 = vector(mode='list',length=length(Samples))
names(FiltDat_v1) = Samples

for(i in Samples){
    FiltDat_v1[[i]] = subset(Samples_AllDat[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 600 & nCount_RNA < 150000 & percent.mt < 10 & scDblFinder =='singlet')
    CellTracking[i,'Filt_v1'] = ncol(FiltDat_v1[[i]])
    cat(paste(i,', ',sep=''))
}

# 7 integrate
# nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 600 & nCount_RNA < 150000 & percent.mt < 10 & scDblFinder =='singlet'
# 3000 HVG, resolution 0.2, regress mito

Samples = list.files(path='../outs/')
NormDat = vector(mode='list',length=length(Samples))
names(NormDat) = Samples

NormDat <- lapply(X = FiltDat_v1, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(object.list = NormDat, nfeatures = 3000)
NormDat <- PrepSCTIntegration(object.list = NormDat, anchor.features = features)
NormDat <- lapply(X = NormDat, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = NormDat, normalization.method = "SCT", anchor.features = features, reduction = "rpca")

FinalDat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(FinalDat) <- "integrated"
FinalDat <- FinalDat %>% RunPCA() %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.2)


###################################################################
####  PLOTS

DimPlot(FinalDat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave("umap_clusters.jpg", width = 5, height = 4)

DimPlot(FinalDat, reduction = "umap", group.by = "group", label = FALSE, raster = TRUE, split.by = "group", pt.size = .1) + scale_color_manual(values = c("black", "turquoise2", "red", "purple2"))
ggsave("umap_group.jpg", width = 12, height = 4)

DimPlot(FinalDat, reduction = "umap", group.by = "core", label = FALSE)
ggsave("umap_core.jpg", width = 5, height = 4)

DimPlot(FinalDat, reduction = "umap", group.by = "sample", label = FALSE)
ggsave("umap_samples.jpg", width = 5, height = 4)

DimPlot(FinalDat, reduction = "umap", group.by = "date", label = FALSE)
ggsave("umap_date.jpg", width = 5, height = 4)

FeaturePlot(object = FinalDat, features = "IL6") + scale_color_gradient(low = "yellow", high = "purple",, trans = "log1p")
ggsave("umap_IL6_man.jpg", width = 5, height = 4)

FeaturePlot(FinalDat, feature = "nFeature_RNA", label = FALSE)
ggsave("umap_genes.jpg", width = 5, height = 4)

FeaturePlot(FinalDat, feature = "nCount_RNA", label = FALSE) + scale_color_gradient(trans = "log")
ggsave("umap_counts.jpg", width = 5, height = 4)

FeaturePlot(FinalDat, feature = "percent.mt", label = FALSE)
ggsave("umap_mito.jpg", width = 5, height = 4)

DefaultAssay(FinalDat) <- "RNA"
FeaturePlot(FinalDat, feature = "nFeature_RNA", label = FALSE)
ggsave("umap_genes_RNA.jpg", width = 5, height = 4)

FeaturePlot(FinalDat, feature = "nCount_RNA", label = FALSE) + scale_color_gradient(trans = "log", low = "lightgrey", high = "blue")
ggsave("umap_counts_RNA_log.jpg", width = 5, height = 4)

###################################################################
#### CLUSTERS

GOI <- c("Ccl3", "Ccl4", "Apoe", "Ccl5", "Ccl12", "Ccl2", "Saa3", "Lcn2", "Spp1", "Postn", "Arg1", "Ifit2", "Cxcl10", "Ccr1", "Ifit1", "Oasl1")
FeaturePlot(FinalDat, features = GOI, cols = c("grey", "red"), max.cutoff = 'q75', min.cutoff = 'q10')
ggsave("CLUSTERS.jpg", width = 20, height = 20)


#### Gene expression
FeaturePlot(FinalDat, features = c("Spp1", "Ccr1", "Ccl5"), cols = c("grey", "red"), max.cutoff = 'q75', min.cutoff = 'q10', split.by = "group")
ggsave("umap_DGE_v1.jpg", width = 20, height = 14)


##################################################################
# Sample specific cells

Idents(object=FinalDat) <- "sample"
head(Idents(FinalDat))
for(i in Samples){
    g <- WhichCells(FinalDat, idents = c(i))
    DimPlot(FinalDat, label=F, cells.highlight= list(g), cols.highlight = c("blue"), cols= "grey")
    ggsave(paste(i,'.jpg',sep=''), width = 5, height = 4)
    }

###################################################################
#### Abundance Assay

library(dplyr)
library(ggplot2)
library(cowplot)

# Filter out low IL-6 samples
load("../raw_counts_meta.rda")
meta3 <- meta %>% filter(IL6 > 50)

dat3 <- meta3 %>% group_by(group, sample, seurat_clusters, IL6) %>% tally()
dat3 <- dat3 %>% group_by(group, sample) %>% mutate(total = sum(n))
dat3 <- dat3 %>% mutate(percent = (n/total)*100)
meta3 %>% group_by(group, sample, IL6) %>% tally()
dat3 %>% filter(seurat_clusters=="1") %>% select(IL6, total)
colnames(dat3) <- c("Group", "Sample", "Cluster","IL6", "n", "total", "percent")

Cluster = unique(dat3$Cluster)
for (i in Cluster){
  tmp <- subset(dat3, Cluster ==i)
  res <- aov(percent ~ Group, data = tmp)
  cat(paste(i))
  print(summary(res))
  print(TukeyHSD(res))
}

dat3 %>% filter(Cluster == "0") %>% ggplot(., aes(x = Cluster, y = percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("grey", "turquoise3", "red", "purple2")) +
  geom_point(size = 1, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
  theme_cowplot() + scale_x_discrete(labels = NULL, breaks = NULL) + theme(legend.position = "none") + ylim(20,40)

 
dat3 %>% filter(Cluster %in% c(1:9)) %>% ggplot(., aes(x = Cluster, y = percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("grey", "turquoise3", "red", "purple2")) +
  geom_point(size = 1, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
  theme_cowplot() + scale_x_discrete(labels = NULL, breaks = NULL) + theme(legend.position = "none") + ylim(0,20)

dat3 %>% filter(Cluster %in% c(10:16)) %>% ggplot(., aes(x = Cluster, y = percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("grey", "turquoise3", "red", "purple2")) +
  geom_point(size = 1, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) +
  theme_cowplot() + scale_x_discrete(labels = NULL, breaks = NULL) + theme(legend.position = "none") + ylim(0,5)


#######################
# Stats on abundance across groups for each cluster
Cluster = unique(dat3$Cluster)
for (i in Cluster){
  tmp <- subset(dat3, Cluster ==i)
  res <- aov(percent ~ Group, data = tmp)
  cat(paste(i))
  print(summary(res))
  print(TukeyHSD(res))
}
#######################
