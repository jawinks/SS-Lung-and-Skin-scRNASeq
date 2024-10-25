library(Seurat)

ssc1.data <- Read10X_h5("H:/My Documents/GSE212109_SScLung/GSM6509488_SC281_ssc1_raw_feature_bc_matrix.h5")
ssc1 <- CreateSeuratObject(counts=ssc1.data$'Gene Expression', project="ssc1", min.cells=3, min.features = 200) 

ssc2.data <- Read10X_h5("H:/My Documents/GSE212109_SScLung/GSM6509489_SC284_ssc2_raw_feature_bc_matrix.h5")
ssc2 <- CreateSeuratObject(counts=ssc2.data$'Gene Expression', project="ssc2", min.cells=3, min.features = 200) 

ssc3.data <- Read10X_h5("H:/My Documents/GSE212109_SScLung/GSM6509490_SC286_ssc3_raw_feature_bc_matrix.h5")
ssc3 <- CreateSeuratObject(counts=ssc3.data$'Gene Expression', project="ssc3", min.cells=3, min.features = 200) 

ssc4.data <- Read10X_h5("H:/My Documents/GSE212109_SScLung/GSM6509491_SC293_ssc4_raw_feature_bc_matrix.h5")
ssc4 <- CreateSeuratObject(counts=ssc4.data, project="ssc4", min.cells=3, min.features = 200) 

ssc5.data <- Read10X_h5("H:/My Documents/GSE212109_SScLung/GSM6509497_SC335_ssc5_raw_feature_bc_matrix.h5")
ssc5 <- CreateSeuratObject(counts=ssc5.data, project="ssc5", min.cells=3, min.features = 200) 

ssc1$condition="SSc"
ssc2$condition="SSc"
ssc3$condition="SSc"
ssc4$condition="SSc"
ssc5$condition="SSc"

ssc1[["percent.mt"]] <- PercentageFeatureSet(object = ssc1, pattern = "^MT-")
ssc2[["percent.mt"]] <- PercentageFeatureSet(object = ssc2, pattern = "^MT-")
ssc3[["percent.mt"]] <- PercentageFeatureSet(object = ssc3, pattern = "^MT-")
ssc4[["percent.mt"]] <- PercentageFeatureSet(object = ssc4, pattern = "^MT-")
ssc5[["percent.mt"]] <- PercentageFeatureSet(object = ssc5, pattern = "^MT-")

rm(ssc1.data, ssc2.data, ssc3.data, ssc4.data, ssc5.data)
gc()

ssc <- merge(ssc1, c(ssc2, ssc3, ssc4, ssc5), merge.data=TRUE)
rm(ssc1, ssc2, ssc3, ssc4, ssc5)
gc()

ssc <- subset(ssc, subset = 
                                  nCount_RNA < 20000 & 
                                  nCount_RNA > 1000 & 
                                  nFeature_RNA > 1000 & 
                                  percent.mt < 20)
ssc <- NormalizeData(ssc)
ssc <- FindVariableFeatures(ssc)
ssc <- ScaleData(ssc)
ssc <- RunPCA(ssc)

options(future.globals.maxSize = 8000 * 1024^2)

ssc <- IntegrateLayers(
  object = ssc,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

ssc <- FindNeighbors(ssc, reduction = "integrated.rpca", dims = 1:30)
ssc <- FindClusters(ssc, resolution = 0.5)
ssc <- RunUMAP(ssc, reduction = "integrated.rpca", dims = 1:30)
DimPlot(ssc, reduction = "umap", label = TRUE)
ssc[["RNA"]] <- JoinLayers(ssc[["RNA"]])

FeaturePlot(ssc, features = c("COL1A1", "PDGFRA", "FN1"))
fibroblasts <- subset(ssc, idents = c(1, 10))

Idents(fibroblasts) <- fibroblasts$orig.ident
cpm <- AggregateExpression(fibroblasts, return.seurat = TRUE)
cpm <- NormalizeData(cpm, normalization.method = "RC", scale.factor = 1e6, verbose= TRUE)
ssclung_cpm <- LayerData(cpm, assay = "RNA", layer = "data")
write.csv(ssclung_cpm, "ssc_lung_cpm_updated.csv")
