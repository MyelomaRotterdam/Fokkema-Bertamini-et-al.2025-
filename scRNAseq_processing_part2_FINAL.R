# --------------------------------------------------- #
# Preprocessing of scRNA sequencing data as described in Fokkema, Bertamini et al. 2025
# Script by Cathelijne Fokkema and Luca Bertamini, optimized for speed by Mathijs Sanders and Gregory van Beek 
# Dept. of Hematology, Erasmus MC Cancer Institute, Rotterdam, the Netherlands
# --------------------------------------------------- #

# --------------------------------------------------- #
# READ-ME (IMPORTANT)
# This script is meant to generate the myeloma tumor cell object by in silico selection of tumor cells 
# from the bone marrow of 25 patients . 
# The data from individual samples is first loaded, tumor cells are then subsetted from the total object by selection 
# on "SDC1", "TNFRSF17", "SLAMF7", "JCHAIN","XBP1", "HLA-DR" expression  as well as known highly expressed markers related
# to primary genetic alterations ("CCND1", "CCND2", "CCND3", "MYEOV","MAF","MAFB","MAFA", "FGFR3", "WHSC1") 
# in patients with whole genome sequencing data,  somatic mutations were called in the scRNAseq to exclude normal plasma cells 
# Tumor cells are merged into a final object where subsequent analysis are performed
# --------------------------------------------------- #


# Loading libraries  --------------------------------------------------------
print("Loading packages ...")

packages <- c('Seurat', 'future', 'parallel', 'cowplot', 'ggplot2', 'dplyr', 'readxl')

invisible(suppressMessages(suppressWarnings(lapply(packages, library, character.only = TRUE))))

options(future.rng.onMisuse="ignore")
options(future.globals.maxSize = 2e4 * 1024^2)
options(future.fork.enable = TRUE)

plan(strategy = 'multicore', workers = 30)

input_dir="..."

## MM_01_BM ----------------------------------------------------------------

# Loading and pre-processing MM_01_BM 
MM_01_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_01_BM", "filtered_feature_bc_matrix"))
MM_01_BM <- CreateSeuratObject(counts = MM_01_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_01_BM <- grep(pattern = "^MT-", x=rownames(x=MM_01_BM), value=T)
percent.mito_MM_01_BM <- Matrix::colSums(x = GetAssayData(object = MM_01_BM, slot="counts")[mito.features_MM_01_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_01_BM, slot = "counts"))
MM_01_BM[["percent.mito"]] <- percent.mito_MM_01_BM
VlnPlot(object = MM_01_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_01_BM <- subset(x = MM_01_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3500 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.10)
# remove non tumor cells
select.cells_MM_01_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_01_BM.csv"), header=T, row.names = 1)
select.cells_MM_01_BM <- as.character(select.cells_MM_01_BM$x)
select.cells_MM_01_BM <- gsub("_.*", "", select.cells_MM_01_BM)
MM_01_BM <- subset(x = MM_01_BM, cells = select.cells_MM_01_BM)
# remove IG genes
counts_MM_01_BM <- GetAssayData(MM_01_BM, assay = "RNA") 
counts_IGexclude_MM_01_BM <- counts_MM_01_BM[-which(rownames(counts_MM_01_BM)%in% rownames(counts_MM_01_BM)[str_detect(rownames(counts_MM_01_BM),"^IG")]),]
MM_01_BM <- subset(MM_01_BM, features = rownames(counts_IGexclude_MM_01_BM))
# normalize and find variables
MM_01_BM <- NormalizeData(object = MM_01_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_01_BM <- FindVariableFeatures(object = MM_01_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_01_BM[["Source"]] <- "BM"
MM_01_BM[["Sample_ID"]] <- "MM_01_BM"
MM_01_BM[["ID"]] <- "MM_01"

## MM_02_BM ----------------------------------------------------------------

# Loading and pre-processing MM_02_BM 
MM_02_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_02_BM", "filtered_feature_bc_matrix"))
MM_02_BM <- CreateSeuratObject(counts = MM_02_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_02_BM <- grep(pattern = "^MT-", x=rownames(x=MM_02_BM), value=T)
percent.mito_MM_02_BM <- Matrix::colSums(x = GetAssayData(object = MM_02_BM, slot="counts")[mito.features_MM_02_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_02_BM, slot = "counts"))
MM_02_BM[["percent.mito"]] <- percent.mito_MM_02_BM
VlnPlot(object = MM_02_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_02_BM <- subset(x = MM_02_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_02_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_02_BM.csv"), header=T, row.names = 1)
select.cells_MM_02_BM <- as.character(select.cells_MM_02_BM$x)
select.cells_MM_02_BM <- gsub("_.*", "", select.cells_MM_02_BM)
MM_02_BM <- subset(x = MM_02_BM, cells = select.cells_MM_02_BM)
# remove IG genes
counts_MM_02_BM <- GetAssayData(MM_02_BM, assay = "RNA") 
counts_IGexclude_MM_02_BM<-counts_MM_02_BM[-which(rownames(counts_MM_02_BM)%in% rownames(counts_MM_02_BM)[str_detect(rownames(counts_MM_02_BM),"^IG")]),]
MM_02_BM <- subset(MM_02_BM, features = rownames(counts_IGexclude_MM_02_BM))
# normalize and find variables
MM_02_BM <- NormalizeData(object = MM_02_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_02_BM <- FindVariableFeatures(object = MM_02_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_02_BM[["Source"]] <- "BM"
MM_02_BM[["Sample_ID"]] <- "MM_02_BM"
MM_02_BM[["ID"]] <- "MM_02"

## MM_03_BM ----------------------------------------------------------------

# Loading and pre-processing MM_03_BM 
MM_03_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_03_BM", "filtered_feature_bc_matrix"))
MM_03_BM <- CreateSeuratObject(counts = MM_03_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_03_BM <- grep(pattern = "^MT-", x=rownames(x=MM_03_BM), value=T)
percent.mito_MM_03_BM <- Matrix::colSums(x = GetAssayData(object = MM_03_BM, slot="counts")[mito.features_MM_03_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_03_BM, slot = "counts"))
MM_03_BM[["percent.mito"]] <- percent.mito_MM_03_BM
VlnPlot(object = MM_03_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_03_BM <- subset(x = MM_03_BM, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_03_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_03_BM.csv"), header=T, row.names = 1)
select.cells_MM_03_BM <- as.character(select.cells_MM_03_BM$x)
select.cells_MM_03_BM <- gsub("_.*", "", select.cells_MM_03_BM)
MM_03_BM <- subset(x = MM_03_BM, cells = select.cells_MM_03_BM)
# remove IG genes
counts_MM_03_BM <- GetAssayData(MM_03_BM, assay = "RNA") 
counts_IGexclude_MM_03_BM <- counts_MM_03_BM[-which(rownames(counts_MM_03_BM)%in% rownames(counts_MM_03_BM)[str_detect(rownames(counts_MM_03_BM),"^IG")]),]
MM_03_BM <- subset(MM_03_BM, features = rownames(counts_IGexclude_MM_03_BM))
# normalize and find variables
MM_03_BM <- NormalizeData(object = MM_03_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_03_BM <- FindVariableFeatures(object = MM_03_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_03_BM[["Source"]] <- "BM"
MM_03_BM[["Sample_ID"]] <- "MM_03_BM"
MM_03_BM[["ID"]] <- "MM_03"

## MM_04_BM ----------------------------------------------------------------

# Loading and pre-processing MM_04_BM 
MM_04_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_04_BM", "filtered_feature_bc_matrix"))
MM_04_BM <- CreateSeuratObject(counts = MM_04_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_object18 <- grep(pattern = "^MT-", x=rownames(x=MM_04_BM), value=T)
percent.mito_object18 <- Matrix::colSums(x = GetAssayData(object = MM_04_BM, slot="counts")[mito.features_object18,]) / Matrix::colSums(x = GetAssayData(object = MM_04_BM, slot = "counts"))
MM_04_BM[["percent.mito"]] <- percent.mito_object18
VlnPlot(object = MM_04_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_04_BM <- subset(x = MM_04_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4500 & nCount_RNA >200 & nCount_RNA <40000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_04_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_04_BM.csv"), header=T, row.names = 1)
select.cells_MM_04_BM <- as.character(select.cells_MM_04_BM$x)
select.cells_MM_04_BM <- gsub("_.*", "", select.cells_MM_04_BM)
MM_04_BM <- subset(x = MM_04_BM, cells = select.cells_MM_04_BM)
# remove IG genes
counts_MM_04_BM <- GetAssayData(MM_04_BM, assay = "RNA") 
counts_IGexclude_MM_04_BM <- counts_MM_04_BM[-which(rownames(counts_MM_04_BM)%in% rownames(counts_MM_04_BM)[str_detect(rownames(counts_MM_04_BM),"^IG")]),]
MM_04_BM <- subset(MM_04_BM, features = rownames(counts_IGexclude_MM_04_BM))
# normalize and find variables
MM_04_BM <- NormalizeData(object = MM_04_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_04_BM <- FindVariableFeatures(object = MM_04_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_04_BM[["Source"]] <- "BM"
MM_04_BM[["Sample_ID"]] <- "MM_04_BM"
MM_04_BM[["ID"]] <- "MM_04"

## MM_05_BM ----------------------------------------------------------------

# Loading and pre-processing MM_05_BM 
MM_05_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_05_BM", "filtered_feature_bc_matrix"))
MM_05_BM <- CreateSeuratObject(counts = MM_05_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_05_BM <- grep(pattern = "^MT-", x=rownames(x=MM_05_BM), value=T)
percent.mito_MM_05_BM <- Matrix::colSums(x = GetAssayData(object = MM_05_BM, slot="counts")[mito.features_MM_05_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_05_BM, slot = "counts"))
MM_05_BM[["percent.mito"]] <- percent.mito_MM_05_BM
VlnPlot(object = MM_05_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_05_BM <- subset(x = MM_05_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)

# remove IG genes
counts_MM_05_BM <- GetAssayData(MM_05_BM, assay = "RNA") 
counts_IGexclude_MM_05_BM <- counts_MM_05_BM[-which(rownames(counts_MM_05_BM)%in% rownames(counts_MM_05_BM)[str_detect(rownames(counts_MM_05_BM),"^IG")]),]
MM_05_BM <- subset(MM_05_BM, features = rownames(counts_IGexclude_MM_05_BM))
# normalize and find variables
MM_05_BM <- NormalizeData(object = MM_05_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_05_BM <- FindVariableFeatures(object = MM_05_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_05_BM[["Source"]] <- "BM"
MM_05_BM[["Sample_ID"]] <- "MM_05_BM"
MM_05_BM[["ID"]] <- "MM_05"

## MM_06_BM ----------------------------------------------------------------

# Loading and pre-processing MM_06_BM
MM_06_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_06_BM", "filtered_feature_bc_matrix"))
MM_06_BM <- CreateSeuratObject(counts = MM_06_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_06_BM <- grep(pattern = "^MT-", x=rownames(x=MM_06_BM), value=T)
percent.mito_MM_06_BM <- Matrix::colSums(x = GetAssayData(object = MM_06_BM, slot="counts")[mito.features_MM_06_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_06_BM, slot = "counts"))
MM_06_BM[["percent.mito"]] <- percent.mito_MM_06_BM
VlnPlot(object = MM_06_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_06_BM <- subset(x = MM_06_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3500 & nCount_RNA >200 & nCount_RNA <40000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_06_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_06_BM.csv"), header=T, row.names = 1)
select.cells_MM_06_BM <- as.character(select.cells_MM_06_BM$x)
select.cells_MM_06_BM <- gsub("_.*", "", select.cells_MM_06_BM)
MM_06_BM <- subset(x = MM_06_BM, cells = select.cells_MM_06_BM)
# remove IG genes
counts_MM_06_BM <- GetAssayData(MM_06_BM, assay = "RNA") 
counts_IGexclude_MM_06_BM <- counts_MM_06_BM[-which(rownames(counts_MM_06_BM)%in% rownames(counts_MM_06_BM)[str_detect(rownames(counts_MM_06_BM),"^IG")]),]
MM_06_BM <- subset(MM_06_BM, features = rownames(counts_IGexclude_MM_06_BM))
# normalize and find variables
MM_06_BM <- NormalizeData(object = MM_06_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_06_BM <- FindVariableFeatures(object = MM_06_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_06_BM[["Source"]] <- "BM"
MM_06_BM[["Sample_ID"]] <- "MM_06_BM"
MM_06_BM[["ID"]] <- "MM_06"

## MM_07_BM ----------------------------------------------------------------

# Loading and pre-processing MM_07_BM
MM_07_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_07_BM", "filtered_feature_bc_matrix"))
MM_07_BM <- CreateSeuratObject(counts = MM_07_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_07_BM <- grep(pattern = "^MT-", x=rownames(x=MM_07_BM), value=T)
percent.mito_MM_07_BM <- Matrix::colSums(x = GetAssayData(object = MM_07_BM, slot="counts")[mito.features_MM_07_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_07_BM, slot = "counts"))
MM_07_BM[["percent.mito"]] <- percent.mito_MM_07_BM
VlnPlot(object = MM_07_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_07_BM <- subset(x = MM_07_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_07_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_07_BM.csv"), header=T, row.names = 1)
select.cells_MM_07_BM <- as.character(select.cells_MM_07_BM$x)
select.cells_MM_07_BM <- gsub("_.*", "", select.cells_MM_07_BM)
MM_07_BM <- subset(x = MM_07_BM, cells = select.cells_MM_07_BM)
# remove IG genes
counts_MM_07_BM <- GetAssayData(MM_07_BM, assay = "RNA") 
counts_IGexclude_MM_07_BM <- counts_MM_07_BM[-which(rownames(counts_MM_07_BM)%in% rownames(counts_MM_07_BM)[str_detect(rownames(counts_MM_07_BM),"^IG")]),]
MM_07_BM <- subset(MM_07_BM, features = rownames(counts_IGexclude_MM_07_BM))
# normalize and find variables
MM_07_BM <- NormalizeData(object = MM_07_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_07_BM <- FindVariableFeatures(object = MM_07_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_07_BM[["Source"]] <- "BM"
MM_07_BM[["Sample_ID"]] <- "MM_07_BM"
MM_07_BM[["ID"]] <- "MM_07"

## MM_08_BM ----------------------------------------------------------------

# Loading and pre-processing MM_08_BM (object 1)
MM_08_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_08_BM", "filtered_feature_bc_matrix"))
MM_08_BM <- CreateSeuratObject(counts = MM_08_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_08_BM <- grep(pattern = "^MT-", x=rownames(x=MM_08_BM), value=T)
percent.mito_MM_08_BM <- Matrix::colSums(x = GetAssayData(object = MM_08_BM, slot="counts")[mito.features_MM_08_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_08_BM, slot = "counts"))
MM_08_BM[["percent.mito"]] <- percent.mito_MM_08_BM
VlnPlot(object = MM_08_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_08_BM <- subset(x = MM_08_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_08_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_08_BM.csv"), header=T, row.names = 1)
select.cells_MM_08_BM=as.character(select.cells_MM_08_BM$x)
select.cells_MM_08_BM <- gsub("_.*", "", select.cells_MM_08_BM)
MM_08_BM <- subset(x = MM_08_BM, cells = select.cells_MM_08_BM)
# remove IG genes
counts_MM_08_BM <- GetAssayData(MM_08_BM, assay = "RNA") 
counts_IGexclude_MM_08_BM <-counts_MM_08_BM[-which(rownames(counts_MM_08_BM)%in% rownames(counts_MM_08_BM)[str_detect(rownames(counts_MM_08_BM),"^IG")]),]
MM_08_BM <- subset(MM_08_BM, features = rownames(counts_IGexclude_MM_08_BM)) # exclusion of Ig genes
# normalize and find variables
MM_08_BM <- NormalizeData(object = MM_08_BM, normalization.method = "LogNormalize", scale.factor = 1e4) # normalization
MM_08_BM <- FindVariableFeatures(object = MM_08_BM, selection.method = "vst", nfeatures = 2000) # find varibales
# minimal annotation
MM_08_BM[["Source"]] <- "BM"
MM_08_BM[["Sample_ID"]] <- "MM_08_BM"
MM_08_BM[["ID"]] <- "MM_08"

## MM_09_BM ----------------------------------------------------------------

# Loading and pre-processing MM_9_BM
MM_09_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_09_BM", "filtered_feature_bc_matrix"))
MM_09_BM <- CreateSeuratObject(counts = MM_09_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_09_BM <- grep(pattern = "^MT-", x=rownames(x=MM_09_BM), value=T)
percent.mito_MM_09_BM <- Matrix::colSums(x = GetAssayData(object = MM_09_BM, slot="counts")[mito.features_MM_09_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_09_BM, slot = "counts"))
MM_09_BM[["percent.mito"]] <- percent.mito_MM_09_BM
VlnPlot(object = MM_09_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_09_BM <- subset(x = MM_09_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_09_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_09_BM.csv"), header=T, row.names = 1)
select.cells_MM_09_BM <- as.character(select.cells_MM_09_BM$x)
select.cells_MM_09_BM <- gsub("_.*", "", select.cells_MM_09_BM)
MM_09_BM <- subset(x = MM_09_BM, cells = select.cells_MM_09_BM)
# remove IG genes
counts_MM_09_BM <- GetAssayData(MM_09_BM, assay = "RNA") 
counts_IGexclude_MM_09_BM <- counts_MM_09_BM[-which(rownames(counts_MM_09_BM)%in% rownames(counts_MM_09_BM)[str_detect(rownames(counts_MM_09_BM),"^IG")]),]
MM_09_BM <- subset(MM_09_BM, features = rownames(counts_IGexclude_MM_09_BM))
# normalize and find variables
MM_09_BM <- NormalizeData(object = MM_09_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_09_BM <- FindVariableFeatures(object = MM_09_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_09_BM[["Source"]] <- "BM"
MM_09_BM[["Sample_ID"]] <- "MM_09_BM"
MM_09_BM[["ID"]] <- "MM_09"

## MM_10_BM ----------------------------------------------------------------

# Loading and pre-processing MM_10_BM 
MM_10_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_10_BM", "filtered_feature_bc_matrix"))
MM_10_BM <- CreateSeuratObject(counts = MM_10_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_10_BM <- grep(pattern = "^MT-", x=rownames(x=MM_10_BM), value=T)
percent.mito_MM_10_BM <- Matrix::colSums(x = GetAssayData(object = MM_10_BM, slot="counts")[mito.features_MM_10_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_10_BM, slot = "counts"))
MM_10_BM[["percent.mito"]] <- percent.mito_MM_10_BM
VlnPlot(object = MM_10_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_10_BM <- subset(x = MM_10_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_10_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_10_BM.csv"), header=T, row.names = 1)
select.cells_MM_10_BM <- as.character(select.cells_MM_10_BM$x)
select.cells_MM_10_BM <- gsub("_.*", "", select.cells_MM_10_BM)
MM_10_BM <- subset(x = MM_10_BM, cells = select.cells_MM_10_BM)
# remove IG genes
counts_MM_10_BM <- GetAssayData(MM_10_BM, assay = "RNA") 
counts_IGexclude_MM_10_BM <- counts_MM_10_BM[-which(rownames(counts_MM_10_BM)%in% rownames(counts_MM_10_BM)[str_detect(rownames(counts_MM_10_BM),"^IG")]),]
MM_10_BM <- subset(MM_10_BM, features = rownames(counts_IGexclude_MM_10_BM))
# normalize and find variables
MM_10_BM <- NormalizeData(object = MM_10_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_10_BM <- FindVariableFeatures(object = MM_10_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_10_BM[["Source"]] <- "BM"
MM_10_BM[["Sample_ID"]] <- "MM_10_BM"
MM_10_BM[["ID"]] <- "MM_10"

## MM_11_BM ----------------------------------------------------------------

# Loading and pre-processing MM_11_BM 
MM_11_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_11_BM", "filtered_feature_bc_matrix"))
MM_11_BM <- CreateSeuratObject(counts = MM_11_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_11_BM <- grep(pattern = "^MT-", x=rownames(x=MM_11_BM), value=T)
percent.mito_MM_11_BM <- Matrix::colSums(x = GetAssayData(object = MM_11_BM, slot="counts")[mito.features_MM_11_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_11_BM, slot = "counts"))
MM_11_BM[["percent.mito"]] <- percent.mito_MM_11_BM
VlnPlot(object = MM_11_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_11_BM <- subset(x = MM_11_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_11_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_11_BM.csv"), header=T, row.names = 1)
select.cells_MM_11_BM <- as.character(select.cells_MM_11_BM$x)
select.cells_MM_11_BM <- gsub("_.*", "", select.cells_MM_11_BM)
MM_11_BM <- subset(x = MM_11_BM, cells = select.cells_MM_11_BM)
# remove IG genes
counts_MM_11_BM <- GetAssayData(MM_11_BM, assay = "RNA") 
counts_IGexclude_counts_MM_11_BM <- counts_MM_11_BM[-which(rownames(counts_MM_11_BM)%in% rownames(counts_MM_11_BM)[str_detect(rownames(counts_MM_11_BM),"^IG")]),]
MM_11_BM <- subset(MM_11_BM, features = rownames(counts_IGexclude_counts_MM_11_BM))
# normalize and find variables
MM_11_BM <- NormalizeData(object = MM_11_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_11_BM <- FindVariableFeatures(object = MM_11_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_11_BM[["Source"]] <- "BM"
MM_11_BM[["Sample_ID"]] <- "MM_11_BM"
MM_11_BM[["ID"]] <- "MM_11"

## MM_12_BM ----------------------------------------------------------------

# Loading and pre-processing MM_12_BM
MM_12_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_12_BM", "filtered_feature_bc_matrix"))
MM_12_BM <- CreateSeuratObject(counts = MM_12_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_12_BM <- grep(pattern = "^MT-", x=rownames(x=MM_12_BM), value=T)
percent.mito_MM_12_BM <- Matrix::colSums(x = GetAssayData(object = MM_12_BM, slot="counts")[mito.features_MM_12_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_12_BM, slot = "counts"))
MM_12_BM[["percent.mito"]] <- percent.mito_MM_12_BM
VlnPlot(object = MM_12_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_12_BM <- subset(x = MM_12_BM, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_12_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_12_BM.csv"), header=T, row.names = 1)
select.cells_MM_12_BM <- as.character(select.cells_MM_12_BM$x)
select.cells_MM_12_BM <- gsub("_.*", "", select.cells_MM_12_BM)
MM_12_BM <- subset(x = MM_12_BM, cells = select.cells_MM_12_BM)
# remove IG genes
counts_MM_12_BM <- GetAssayData(MM_12_BM, assay = "RNA") 
counts_IGexclude_MM_12_BM <- counts_MM_12_BM[-which(rownames(counts_MM_12_BM)%in% rownames(counts_MM_12_BM)[str_detect(rownames(counts_MM_12_BM),"^IG")]),]
MM_12_BM <- subset(MM_12_BM, features = rownames(counts_IGexclude_MM_12_BM))
# normalize and find variables
MM_12_BM <- NormalizeData(object = MM_12_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_12_BM <- FindVariableFeatures(object = MM_12_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_12_BM[["Source"]] <- "BM"
MM_12_BM[["Sample_ID"]] <- "MM_12_BM"
MM_12_BM[["ID"]] <- "MM_12"

## MM_13_BM ----------------------------------------------------------------

# Loading and pre-processing MM_13_BM 
MM_13_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_13_BM", "filtered_feature_bc_matrix"))
MM_13_BM <- CreateSeuratObject(counts = MM_13_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_13_BM <- grep(pattern = "^MT-", x=rownames(x=MM_13_BM), value=T)
percent.mito_MM_13_BM <- Matrix::colSums(x = GetAssayData(object = MM_13_BM, slot="counts")[mito.features_MM_13_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_13_BM, slot = "counts"))
MM_13_BM[["percent.mito"]] <- percent.mito_MM_13_BM
VlnPlot(object = MM_13_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_13_BM <- subset(x = MM_13_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_13_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_13_BM.csv"), header=T, row.names = 1)
select.cells_MM_13_BM <- as.character(select.cells_MM_13_BM$x)
select.cells_MM_13_BM <- gsub("_.*", "", select.cells_MM_13_BM)
MM_13_BM <- subset(x = MM_13_BM, cells = select.cells_MM_13_BM)
# remove IG genes
counts_MM_13_BM <- GetAssayData(MM_13_BM, assay = "RNA") 
counts_IGexclude_MM_13_BM <- counts_MM_13_BM[-which(rownames(counts_MM_13_BM)%in% rownames(counts_MM_13_BM)[str_detect(rownames(counts_MM_13_BM),"^IG")]),]
MM_13_BM <- subset(MM_13_BM, features = rownames(counts_IGexclude_MM_13_BM))
# normalize and find variables
MM_13_BM <- NormalizeData(object = MM_13_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_13_BM <- FindVariableFeatures(object = MM_13_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_13_BM[["Source"]] <- "BM"
MM_13_BM[["Sample_ID"]] <- "MM_13_BM"
MM_13_BM[["ID"]] <- "MM_13"

## MM_14_BM ----------------------------------------------------------------

# Loading and pre-processing MM_14_BM 
MM_14_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_14_BM", "filtered_feature_bc_matrix"))
MM_14_BM <- CreateSeuratObject(counts = MM_14_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_14_BM <- grep(pattern = "^MT-", x=rownames(x=MM_14_BM), value=T)
percent.mito_MM_14_BM <- Matrix::colSums(x = GetAssayData(object = MM_14_BM, slot="counts")[mito.features_MM_14_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_14_BM, slot = "counts"))
MM_14_BM[["percent.mito"]] <- percent.mito_MM_14_BM
VlnPlot(object = MM_14_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_14_BM <- subset(x = MM_14_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_14_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_14_BM.csv"), header=T, row.names = 1)
select.cells_MM_14_BM <- as.character(select.cells_MM_14_BM$x)
select.cells_MM_14_BM <- gsub("_.*", "", select.cells_MM_14_BM)
# remove IG genes
MM_14_BM <- subset(x = MM_14_BM, cells = select.cells_MM_14_BM)
counts_MM_14_BM <- GetAssayData(MM_14_BM, assay = "RNA") 
counts_IGexclude_MM_14_BM <- counts_MM_14_BM[-which(rownames(counts_MM_14_BM)%in% rownames(counts_MM_14_BM)[str_detect(rownames(counts_MM_14_BM),"^IG")]),]
MM_14_BM <- subset(MM_14_BM, features = rownames(counts_IGexclude_MM_14_BM))
# normalize and find variables
MM_14_BM <- NormalizeData(object = MM_14_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_14_BM <- FindVariableFeatures(object = MM_14_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_14_BM[["Source"]] <- "BM"
MM_14_BM[["Sample_ID"]] <- "MM_14_BM"
MM_14_BM[["ID"]] <- "MM_14"

## MM_15_BM ----------------------------------------------------------------

# Loading and pre-processing MM_15_BM
MM_15_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_15_BM", "filtered_feature_bc_matrix"))
MM_15_BM <- CreateSeuratObject(counts = MM_15_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_15_BM <- grep(pattern = "^MT-", x=rownames(x=MM_15_BM), value=T)
percent.mito_MM_15_BM <- Matrix::colSums(x = GetAssayData(object = MM_15_BM, slot="counts")[mito.features_MM_15_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_15_BM, slot = "counts"))
MM_15_BM[["percent.mito"]] <- percent.mito_MM_15_BM
VlnPlot(object = MM_15_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_15_BM <- subset(x = MM_15_BM, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_15_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_15_BM.csv"), header=T, row.names = 1)
select.cells_MM_15_BM <- as.character(select.cells_MM_15_BM$x)
select.cells_MM_15_BM <- gsub("_.*", "", select.cells_MM_15_BM)
MM_15_BM <- subset(x = MM_15_BM, cells = select.cells_MM_15_BM)
# remove IG genes
counts_MM_15_BM <- GetAssayData(MM_15_BM, assay = "RNA") 
counts_IGexclude_MM_15_BM <- counts_MM_15_BM[-which(rownames(counts_MM_15_BM)%in% rownames(counts_MM_15_BM)[str_detect(rownames(counts_MM_15_BM),"^IG")]),]
MM_15_BM <- subset(MM_15_BM, features = rownames(counts_IGexclude_MM_15_BM))
# normalize and find variables
MM_15_BM <- NormalizeData(object = MM_15_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_15_BM <- FindVariableFeatures(object = MM_15_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_15_BM[["Source"]] <- "BM"
MM_15_BM[["Sample_ID"]] <- "MM_15_BM"
MM_15_BM[["ID"]] <- "MM_15"

## MM_16_BM ----------------------------------------------------------------

# Loading and pre-processing MM_16_BM 
MM_16_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_16_BM", "filtered_feature_bc_matrix"))
MM_16_BM <- CreateSeuratObject(counts = MM_16_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_16_BM <- grep(pattern = "^MT-", x=rownames(x=MM_16_BM), value=T)
percent.mito_MM_16_BM <- Matrix::colSums(x = GetAssayData(object = MM_16_BM, slot="counts")[mito.features_MM_16_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_16_BM, slot = "counts"))
MM_16_BM[["percent.mito"]] <- percent.mito_MM_16_BM
VlnPlot(object = MM_16_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_16_BM <- subset(x = MM_16_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <12000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_16_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_16_BM.csv"), header=T, row.names = 1)
select.cells_MM_16_BM <- as.character(select.cells_MM_16_BM$x)
select.cells_MM_16_BM <- gsub("_.*", "", select.cells_MM_16_BM)
MM_16_BM <- subset(x = MM_16_BM, cells = select.cells_MM_16_BM)
# remove IG genes
counts_MM_16_BM <- GetAssayData(MM_16_BM, assay = "RNA") 
counts_IGexclude_MM_16_BM <- counts_MM_16_BM[-which(rownames(counts_MM_16_BM)%in% rownames(counts_MM_16_BM)[str_detect(rownames(counts_MM_16_BM),"^IG")]),]
MM_16_BM <- subset(MM_16_BM, features = rownames(counts_IGexclude_MM_16_BM))
# normalize and find variables
MM_16_BM <- NormalizeData(object = MM_16_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_16_BM <- FindVariableFeatures(object = MM_16_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_16_BM[["Source"]] <- "BM"
MM_16_BM[["Sample_ID"]] <- "MM_16_BM"
MM_16_BM[["ID"]] <- "MM_16"

## MM_17_BM ----------------------------------------------------------------

# Loading and pre-processing MM_17_BM 
MM_17_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_17_BM", "filtered_feature_bc_matrix"))
MM_17_BM <- CreateSeuratObject(counts = MM_17_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_17_BM <- grep(pattern = "^MT-", x=rownames(x=MM_17_BM), value=T)
percent.mito_MM_17_BM <- Matrix::colSums(x = GetAssayData(object = MM_17_BM, slot="counts")[mito.features_MM_17_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_17_BM, slot = "counts"))
MM_17_BM[["percent.mito"]] <- percent.mito_MM_17_BM
VlnPlot(object = MM_17_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_17_BM <- subset(x = MM_17_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <25000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_17_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_17_BM.csv"), header=T, row.names = 1)
select.cells_MM_17_BM <- as.character(select.cells_MM_17_BM$x)
select.cells_MM_17_BM <- gsub("_.*", "", select.cells_MM_17_BM)
MM_17_BM <- subset(x = MM_17_BM, cells = select.cells_MM_17_BM)
# remove IG genes
counts_MM_17_BM <- GetAssayData(MM_17_BM, assay = "RNA") 
counts_IGexclude_MM_17_BM <- counts_MM_17_BM[-which(rownames(counts_MM_17_BM)%in% rownames(counts_MM_17_BM)[str_detect(rownames(counts_MM_17_BM),"^IG")]),]
MM_17_BM <- subset(MM_17_BM, features = rownames(counts_IGexclude_MM_17_BM))
# normalize and find variables
MM_17_BM <- NormalizeData(object = MM_17_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_17_BM <- FindVariableFeatures(object = MM_17_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_17_BM[["Source"]] <- "BM"
MM_17_BM[["Sample_ID"]] <- "MM_17_BM"
MM_17_BM[["ID"]] <- "MM_17"

## MM_18_BM ----------------------------------------------------------------

# Loading and pre-processing MM_18_BM 
MM_18_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_18_BM", "filtered_feature_bc_matrix"))
MM_18_BM <- CreateSeuratObject(counts = MM_18_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_18_BM <- grep(pattern = "^MT-", x=rownames(x=MM_18_BM), value=T)
percent.mito_MM_18_BM <- Matrix::colSums(x = GetAssayData(object = MM_18_BM, slot="counts")[mito.features_MM_18_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_18_BM, slot = "counts"))
MM_18_BM[["percent.mito"]] <- percent.mito_MM_18_BM
VlnPlot(object = MM_18_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_18_BM <- subset(x = MM_18_BM, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <15000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_18_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_18_BM.csv"), header=T, row.names = 1)
select.cells_MM_18_BM <- as.character(select.cells_MM_18_BM$x)
select.cells_MM_18_BM <- gsub("_.*", "", select.cells_MM_18_BM)
MM_18_BM <- subset(x = MM_18_BM, cells = select.cells_MM_18_BM)
# remove IG genes
counts_MM_18_BM <- GetAssayData(MM_18_BM, assay = "RNA") 
counts_IGexclude_MM_18_BM <- counts_MM_18_BM[-which(rownames(counts_MM_18_BM)%in% rownames(counts_MM_18_BM)[str_detect(rownames(counts_MM_18_BM),"^IG")]),]
MM_18_BM <- subset(MM_18_BM, features = rownames(counts_IGexclude_MM_18_BM))
# normalize and find variables
MM_18_BM <- NormalizeData(object = MM_18_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_18_BM <- FindVariableFeatures(object = MM_18_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_18_BM[["Source"]] <- "BM"
MM_18_BM[["Sample_ID"]] <- "MM_18_BM"
MM_18_BM[["ID"]] <- "MM_18"

## MM_19_BM ----------------------------------------------------------------

# Loading and pre-processing MM_19_BM
MM_19_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_19_BM", "filtered_feature_bc_matrix"))
MM_19_BM <- CreateSeuratObject(counts = MM_19_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_19_BM <- grep(pattern = "^MT-", x=rownames(x=MM_19_BM), value=T)
percent.mito_MM_19_BM <- Matrix::colSums(x = GetAssayData(object = MM_19_BM, slot="counts")[mito.features_MM_19_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_19_BM, slot = "counts"))
MM_19_BM[["percent.mito"]] <- percent.mito_MM_19_BM
VlnPlot(object = MM_19_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_19_BM <- subset(x = MM_19_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.7)
# remove non tumor cells
select.cells_MM_19_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_19_BM_v2.csv"))
select.cells_MM_19_BM <- as.character(select.cells_MM_19_BM$x)
select.cells_MM_19_BM <- gsub("_.*", "", select.cells_MM_19_BM)
MM_19_BM <- subset(x = MM_19_BM, cells = select.cells_MM_19_BM)
# remove IG genes
counts_MM_19_BM <- GetAssayData(MM_19_BM, assay = "RNA") 
counts_IGexclude_MM_19_BM <- counts_MM_19_BM[-which(rownames(counts_MM_19_BM)%in% rownames(counts_MM_19_BM)[str_detect(rownames(counts_MM_19_BM),"^IG")]),]
MM_19_BM <- subset(MM_19_BM, features = rownames(counts_IGexclude_MM_19_BM))
# normalize and find variables
MM_19_BM <- NormalizeData(object = MM_19_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_19_BM <- FindVariableFeatures(object = MM_19_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_19_BM[["Source"]] <- "BM"
MM_19_BM[["Sample_ID"]] <- "MM_19_BM"
MM_19_BM[["ID"]] <- "MM_19"

## MM_20_BM ----------------------------------------------------------------

# Loading and pre-processing MM_20_BM 
MM_20_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_20_BM", "filtered_feature_bc_matrix"))
MM_20_BM <- CreateSeuratObject(counts = MM_20_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_20_BM <- grep(pattern = "^MT-", x=rownames(x=MM_20_BM), value=T)
percent.mito_MM_20_BM <- Matrix::colSums(x = GetAssayData(object = MM_20_BM, slot="counts")[mito.features_MM_20_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_20_BM, slot = "counts"))
MM_20_BM[["percent.mito"]] <- percent.mito_MM_20_BM
VlnPlot(object = MM_20_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_20_BM <- subset(x = MM_20_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_20_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_20_BM_v2.csv"))
select.cells_MM_20_BM <- as.character(select.cells_MM_20_BM$x)
select.cells_MM_20_BM <- gsub("_.*", "", select.cells_MM_20_BM)
MM_20_BM <- subset(x = MM_20_BM, cells = select.cells_MM_20_BM)
# remove IG genes
counts_MM_20_BM <- GetAssayData(MM_20_BM, assay = "RNA") 
counts_IGexclude_MM_20_BM<- counts_MM_20_BM[-which(rownames(counts_MM_20_BM)%in% rownames(counts_MM_20_BM)[str_detect(rownames(counts_MM_20_BM),"^IG")]),]
MM_20_BM <- subset(MM_20_BM, features = rownames(counts_IGexclude_MM_20_BM))
# normalize and find variables
MM_20_BM <- NormalizeData(object = MM_20_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_20_BM <- FindVariableFeatures(object = MM_20_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_20_BM[["Source"]] <- "BM"
MM_20_BM[["Sample_ID"]] <- "MM_20_BM"
MM_20_BM[["ID"]] <- "MM_20"

## MM_21_BM ----------------------------------------------------------------

# Loading and pre-processing MM_21_BM
MM_21_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_21_BM", "filtered_feature_bc_matrix"))
MM_21_BM <- CreateSeuratObject(counts = MM_21_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_21_BM <- grep(pattern = "^MT-", x=rownames(x=MM_21_BM), value=T)
percent.mito_MM_21_BM <- Matrix::colSums(x = GetAssayData(object = MM_21_BM, slot="counts")[mito.features_MM_21_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_21_BM, slot = "counts"))
MM_21_BM[["percent.mito"]] <- percent.mito_MM_21_BM
VlnPlot(object = MM_21_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_21_BM <- subset(x = MM_21_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_21_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_21_BM_v2.csv"))
select.cells_MM_21_BM <- as.character(select.cells_MM_21_BM$x)
select.cells_MM_21_BM <- gsub("_.*", "", select.cells_MM_21_BM)
MM_21_BM <- subset(x = MM_21_BM, cells = select.cells_MM_21_BM)
# remove IG genes
counts_MM_21_BM <- GetAssayData(MM_21_BM, assay = "RNA") 
counts_IGexclude_MM_21_BM <- counts_MM_21_BM[-which(rownames(counts_MM_21_BM)%in% rownames(counts_MM_21_BM)[str_detect(rownames(counts_MM_21_BM),"^IG")]),]
MM_21_BM <- subset(MM_21_BM, features = rownames(counts_IGexclude_MM_21_BM))
# normalize and find variables
MM_21_BM <- NormalizeData(object = MM_21_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_21_BM <- FindVariableFeatures(object = MM_21_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_21_BM[["Source"]] <- "BM"
MM_21_BM[["Sample_ID"]] <- "MM_21_BM"
MM_21_BM[["ID"]] <- "MM_21"

## MM_22_BM ----------------------------------------------------------------

# Loading and pre-processing MM_22_BM 
MM_22_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_22_BM", "filtered_feature_bc_matrix"))
MM_22_BM <- CreateSeuratObject(counts = MM_22_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_22_BM <- grep(pattern = "^MT-", x=rownames(x=MM_22_BM), value=T)
percent.mito_MM_22_BM <- Matrix::colSums(x = GetAssayData(object = MM_22_BM, slot="counts")[mito.features_MM_22_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_22_BM, slot = "counts"))
MM_22_BM[["percent.mito"]] <- percent.mito_MM_22_BM
VlnPlot(object = MM_22_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_22_BM <- subset(x = MM_22_BM, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_22_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_22_BM.csv"), header=T, row.names = 1)
select.cells_MM_22_BM <- as.character(select.cells_MM_22_BM$x)
select.cells_MM_22_BM <- gsub("_.*", "", select.cells_MM_22_BM)
MM_22_BM <- subset(x = MM_22_BM, cells = select.cells_MM_22_BM)
# remove IG genes
counts_MM_22_BM <- GetAssayData(MM_22_BM, assay = "RNA") 
counts_IGexclude_MM_22_BM <- counts_MM_22_BM[-which(rownames(counts_MM_22_BM)%in% rownames(counts_MM_22_BM)[str_detect(rownames(counts_MM_22_BM),"^IG")]),]
MM_22_BM <- subset(MM_22_BM, features = rownames(counts_IGexclude_MM_22_BM))
# normalize and find variables
MM_22_BM <- NormalizeData(object = MM_22_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_22_BM <- FindVariableFeatures(object = MM_22_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_22_BM[["Source"]] <- "BM"
MM_22_BM[["Sample_ID"]] <- "MM_22_BM"
MM_22_BM[["ID"]] <- "MM_22"

## MM_23_BM ----------------------------------------------------------------

# Loading and pre-processing MM_23_BM 
MM_23_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_23_BM", "filtered_feature_bc_matrix"))
MM_23_BM <- CreateSeuratObject(counts = MM_23_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_23_BM <- grep(pattern = "^MT-", x=rownames(x=MM_23_BM), value=T)
percent.mito_MM_23_BM <- Matrix::colSums(x = GetAssayData(object = MM_23_BM, slot="counts")[mito.features_MM_23_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_23_BM, slot = "counts"))
MM_23_BM[["percent.mito"]] <- percent.mito_MM_23_BM
VlnPlot(object = MM_23_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_23_BM <- subset(x = MM_23_BM, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.05)
# remove non tumor cells
select.cells_MM_23_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_23_BM.csv"), header=T, row.names = 1)
select.cells_MM_23_BM <- as.character(select.cells_MM_23_BM$x)
select.cells_MM_23_BM <- gsub("_.*", "", select.cells_MM_23_BM)
MM_23_BM <- subset(x = MM_23_BM, cells = select.cells_MM_23_BM)
# remove IG genes
counts_MM_23_BM <- GetAssayData(MM_23_BM, assay = "RNA") 
counts_IGexclude_MM_23_BM <- counts_MM_23_BM[-which(rownames(counts_MM_23_BM)%in% rownames(counts_MM_23_BM)[str_detect(rownames(counts_MM_23_BM),"^IG")]),]
MM_23_BM <- subset(MM_23_BM, features = rownames(counts_IGexclude_MM_23_BM))
# normalize and find variables
MM_23_BM <- NormalizeData(object = MM_23_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_23_BM <- FindVariableFeatures(object = MM_23_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_23_BM[["Source"]] <- "BM"
MM_23_BM[["Sample_ID"]] <- "MM_23_BM"
MM_23_BM[["ID"]] <- "MM_23"

## MM_24_BM ----------------------------------------------------------------

# Loading and pre-processing MM_24_BM 
MM_24_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_24_BM", "filtered_feature_bc_matrix"))
MM_24_BM <- CreateSeuratObject(counts = MM_24_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_24_BM <- grep(pattern = "^MT-", x=rownames(x=MM_24_BM), value=T)
percent.mito_MM_24_BM <- Matrix::colSums(x = GetAssayData(object = MM_24_BM, slot="counts")[mito.features_MM_24_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_24_BM, slot = "counts"))
MM_24_BM[["percent.mito"]] <- percent.mito_MM_24_BM
VlnPlot(object = MM_24_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_24_BM <- subset(x = MM_24_BM, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.15)
# remove non tumor cells
select.cells_MM_24_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_24_BM.csv"), header=T, row.names = 1)
select.cells_MM_24_BM <- as.character(select.cells_MM_24_BM$x)
select.cells_MM_24_BM <- gsub("_.*", "", select.cells_MM_24_BM)
MM_24_BM <- subset(x = MM_24_BM, cells = select.cells_MM_24_BM)
# remove IG genes
counts_MM_24_BM <- GetAssayData(MM_24_BM, assay = "RNA") 
counts_IGexclude_MM_24_BM <-counts_MM_24_BM[-which(rownames(counts_MM_24_BM)%in% rownames(counts_MM_24_BM)[str_detect(rownames(counts_MM_24_BM),"^IG")]),]
MM_24_BM <- subset(MM_24_BM, features = rownames(counts_IGexclude_MM_24_BM))
# normalize and find variables
MM_24_BM <- NormalizeData(object = MM_24_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_24_BM <- FindVariableFeatures(object = MM_24_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_24_BM[["Source"]] <- "BM"
MM_24_BM[["Sample_ID"]] <- "MM_24_BM"
MM_24_BM[["ID"]] <- "MM_24"

## MM_25_BM ----------------------------------------------------------------

# Loading and pre-processing MM_25_BM 
MM_25_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_25_BM", "filtered_feature_bc_matrix"))
MM_25_BM <- CreateSeuratObject(counts = MM_25_BM, min.cells = 3, min.features = 200, project = "Plasmacells")
# remove low quality cells
mito.features_MM_25_BM <- grep(pattern = "^MT-", x=rownames(x=MM_25_BM), value=T)
percent.mito_MM_25_BM <- Matrix::colSums(x = GetAssayData(object = MM_25_BM, slot="counts")[mito.features_MM_25_BM,]) / Matrix::colSums(x = GetAssayData(object = MM_25_BM, slot = "counts"))
MM_25_BM[["percent.mito"]] <- percent.mito_MM_25_BM
VlnPlot(object = MM_25_BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
MM_25_BM <- subset(x = MM_25_BM, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
# remove non tumor cells
select.cells_MM_25_BM <- read.csv(file=file.path(input_dir, "barcodes", "barcodes_MM_25_BM.csv"), header=T, row.names = 1)
select.cells_MM_25_BM <- as.character(select.cells_MM_25_BM$x)
select.cells_MM_25_BM <- gsub("_.*", "", select.cells_MM_25_BM)
MM_25_BM <- subset(x = MM_25_BM, cells = select.cells_MM_25_BM)
# remove IG genes
counts_MM_25_BM <- GetAssayData(MM_25_BM, assay = "RNA") 
counts_IGexclude_MM_25_BM <- counts_MM_25_BM[-which(rownames(counts_MM_25_BM)%in% rownames(counts_MM_25_BM)[str_detect(rownames(counts_MM_25_BM),"^IG")]),]
MM_25_BM <- subset(MM_25_BM, features = rownames(counts_IGexclude_MM_25_BM))
# normalize and find variables
MM_25_BM <- NormalizeData(object = MM_25_BM, normalization.method = "LogNormalize", scale.factor = 1e4)
MM_25_BM <- FindVariableFeatures(object = MM_25_BM, selection.method = "vst", nfeatures = 2000)
# minimal annotation
MM_25_BM[["Source"]] <- "BM"
MM_25_BM[["Sample_ID"]] <- "MM_25_BM"
MM_25_BM[["ID"]] <- "MM_25"

# create sample list ------------------------------------------------------

#to modify
data_list=list(MM_01_BM, MM_02_BM, MM_03_BM, MM_04_BM, MM_05_BM, MM_06_BM, MM_07_BM, 
                    MM_08_BM, MM_09_BM, MM_10_BM, MM_11_BM, MM_12_BM, MM_13_BM, MM_14_BM, MM_15_BM,
                    MM_16_BM, MM_17_BM, MM_18_BM, MM_19_BM, MM_20_BM, MM_21_BM, MM_22_BM, MM_23_BM,
                    MM_24_BM, MM_25_BM )

names(data_list) <- c("MM_01_BM",  "MM_02_BM",  "MM_03_BM",   "MM_04_BM",  "MM_05_BM",  "MM_06_BM", "MM_07_BM",
                           "MM_08_BM", "MM_09_BM", "MM_10_BM", "MM_11_BM", "MM_12_BM", "MM_13_BM", "MM_14_BM", "MM_15_BM",
                           "MM_16_BM", "MM_17_BM", "MM_18_BM", "MM_19_BM", "MM_20_BM", "MM_21_BM", "MM_22_BM", "MM_23_BM",
                           "MM_24_BM", "MM_25_BM" )


# load and add metadata-----------------------------------------------------

meta_data=read_xlsx(file.path(input_dir, "/meta_data_scrnaseq.xlsx"))

colnames(meta_data)

for(i in seq_along(data_list)) {
  meta_data_extracted <- data_list[[i]]@meta.data
  meta_data_extracted_rows <- rownames(meta_data_extracted)
  
  meta_data_extracted = meta_data_extracted%>%
    mutate(barcode = rownames(meta_data_extracted)) # metadata need to have the same sample name
  
  meta_data_extracted <- meta_data_extracted%>% 
    dplyr::select(orig.ident,nCount_RNA, nFeature_RNA, percent.mito, ID, Sample_ID, Source,
                  barcode)%>%
    left_join(meta_data%>% 
                dplyr::select(ID,Age, Gender, ISS, CTC_perc, CTC_level_007, BMMC_perc, Isotype, ISS,
                              Cytogenetics,`t(11_14)_CCND1`, `t(4_14)_MMSET`, `t(14_16)_MAF`, 
                              `t(14_20)_MAFB`, `t(8_14)_MAFA` ,`del17p`,`amp1q`, `del1p`, `del13q`, `MYC_t_amp`),
              by=c("ID"))
  
  rownames(meta_data_extracted) <- meta_data_extracted_rows
  
  data_list[[i]]@meta.data <- meta_data_extracted
}


# merge all samples-------------------------------------------------------------------------

BM_HvsL <- merge(x = data_list[[1]], 
                               y = c(data_list[2:length(data_list)]), 
                               merge.data = TRUE)
# perform SCT, run PCA and UMAP -------------------------------------------

BM_HvsL <- SCTransform(BM_HvsL, method = "glmGamPoi", 
                                        verbose = TRUE, variable.features.n = 3000)
BM_HvsL <- RunPCA(BM_HvsL, npcs = 30, verbose = TRUE)
BM_HvsL <- RunUMAP(BM_HvsL, dims = 1:30, reduction = 'pca', verbose = TRUE)

BM_HvsL <- FindNeighbors(BM_HvsL, dims = 1:30, verbose = TRUE)
BM_HvsL <- FindClusters(BM_HvsL, resolution = c(0.1, 0.2, 0.4, 0.5), verbose = TRUE)





