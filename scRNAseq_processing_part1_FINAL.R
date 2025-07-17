
# --------------------------------------------------- #
# Processing of scRNA sequencing data as described in Fokkema, Bertamini et al. 2025
# Script by Cathelijne Fokkema and Luca Bertamini, optimized for speed by Mathijs Sanders and Gregory van Beek 
# Dept. of Hematology, Erasmus MC Cancer Institute, Rotterdam, the Netherlands
# --------------------------------------------------- #

# --------------------------------------------------- #
# READ-ME (IMPORTANT)
# This script is meant to generate the myeloma tumor cells object by in silico selection of tumor cells from the 7 patients paired
# bone marrow and peripheral blood cells.
# The data from individual samples is first loaded, tumor cells are then subsetted from the total object by selection 
# on "SDC1", "TNFRSF17", "SLAMF7", "JCHAIN","XBP1", "HLA-DR" expression  as well as known highly expressed markers related
# to primary genetic alterations ("CCND1", "CCND2", "CCND3", "MYEOV","MAF","MAFB","MAFA", "FGFR3", "WHSC1") 
# in patients with whole genome sequencing data,  somatic mutations were called in the scRNAseq to exclude normal plasma cells 
# Tumor cells are merged into a final object where subsequent analysis are performed
# --------------------------------------------------- #


### Loading libraries -------------------------------------------------------------------------
print("Loading packages ...")

packages <- c('Seurat', 'future', 'parallel', 'cowplot', 'ggplot2', 'dplyr')

invisible(suppressMessages(suppressWarnings(lapply(packages, library, character.only = TRUE))))

options(future.rng.onMisuse="ignore")
options(future.globals.maxSize = 2e4 * 1024^2)
options(future.fork.enable = TRUE)

plan(strategy = 'multicore', workers = 30)

input_dir="..." 

# MM_01_BM  ---------------------------------------------------------------
MM_01_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_01_BM", "filtered_feature_bc_matrix"))
MM_01_BM <- MM_01_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_01_BM), perl = TRUE),]
MM_01_BM <- CreateSeuratObject(MM_01_BM, project = 'Plasmacells')
select.cells_MM_01_BM <- fread(file.path(input_dir, "barcodes", "barcodes_MM_01_BM.csv"), header = TRUE, sep = ',')
MM_01_BM <- subset(MM_01_BM, cells = select.cells_MM_01_BM$x)
MM_01_BM[["Source"]] <- "BM"
MM_01_BM[["Sample_ID"]] <- "MM_01_BM"
MM_01_BM[["ID"]] <- "MM_01"

# MM_01_CTC  ---------------------------------------------------------------
MM_01_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_01_CTC", "filtered_feature_bc_matrix"))
MM_01_CTC <- MM_01_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_01_CTC), perl = TRUE),]
MM_01_CTC <- CreateSeuratObject(MM_01_CTC, project = 'Plasmacells')
select.cells_MM_01_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_01_CTC.csv"), header = TRUE, sep = ',')
MM_01_CTC <- subset(MM_01_CTC, cells = select.cells_MM_01_CTC$x)
MM_01_CTC[["Source"]] <- "CTC"
MM_01_CTC[["Sample_ID"]] <- "MM_01_CTC"
MM_01_CTC[["ID"]] <- "MM_01"

# MM_02_BM ---------------------------------------------------------------
MM_02_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_02_BM", "filtered_feature_bc_matrix"))
MM_02_BM <- MM_02_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_02_BM), perl = TRUE),]
MM_02_BM <- CreateSeuratObject(MM_02_BM, project = 'Plasmacells')
select.cells_MM_02_BM <- fread(file.path(input_dir, "barcodes", "barcodes_MM_02_BM.csv"), header = TRUE, sep = ',')
MM_02_BM <- subset(MM_02_BM, cells = select.cells_MM_02_BM$x)
MM_02_BM[["Source"]] <- "BM"
MM_02_BM[["Sample_ID"]] <- "MM_02_BM"
MM_02_BM[["ID"]] <- "MM_02"

# MM_02_CTC  ---------------------------------------------------------------
MM_02_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_02_CTC", "filtered_feature_bc_matrix"))
MM_02_CTC <- MM_02_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_02_CTC), perl = TRUE),]
MM_02_CTC <- CreateSeuratObject(MM_02_CTC, project = 'Plasmacells')
select.cells_MM_02_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_02_CTC.csv"), header = TRUE, sep = ',')
MM_02_CTC <- subset(MM_02_CTC, cells = select.cells_MM_02_CTC$x)
MM_02_CTC[["Source"]] <- "CTC"
MM_02_CTC[["Sample_ID"]] <- "MM_02_CTC"
MM_02_CTC[["ID"]] <- "MM_02"

# MM_03_BM ---------------------------------------------------------------
MM_03_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_03_BM", "filtered_feature_bc_matrix"))
MM_03_BM <- MM_03_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_03_BM), perl = TRUE),]
MM_03_BM <- CreateSeuratObject(MM_03_BM, project = 'Plasmacells')
select.cells_MM_03_BM<- fread(file.path(input_dir, "barcodes", "barcodes_MM_03_BM.csv"), header = TRUE, sep = ',')
MM_03_BM <- subset(MM_03_BM, cells = select.cells_MM_03_BM$x)
MM_03_BM[["Source"]] <- "BM"
MM_03_BM[["Sample_ID"]] <- "MM_03_BM"
MM_02_BM[["ID"]] <- "MM_03"


# MM_03_CTC  ---------------------------------------------------------------
MM_03_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_03_CTC", "filtered_feature_bc_matrix"))
MM_03_CTC <- MM_03_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_03_CTC), perl = TRUE),]
MM_03_CTC <- CreateSeuratObject(MM_03_CTC, project = 'Plasmacells')
select.cells_MM_03_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_03_CTC.csv"), header = TRUE, sep = ',')
MM_03_CTC <- subset(MM_03_CTC, cells = select.cells_MM_03_CTC$x)
MM_03_CTC[["Source"]] <- "CTC"
MM_03_CTC[["Sample_ID"]] <- "MM_03_CTC"
MM_03_CTC[["ID"]] <- "MM_03_CTC"

# MM_04_BM  ---------------------------------------------------------------
MM_04_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_04_BM", "filtered_feature_bc_matrix"))
MM_04_BM <- MM_04_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_04_BM), perl = TRUE),]
MM_04_BM <- CreateSeuratObject(MM_04_BM, project = 'Plasmacells')
select.cells_MM_04_BM<- fread(file.path(input_dir, "barcodes", "barcodes_MM_04_BM.csv"), header = TRUE, sep = ',')
MM_04_BM <- subset(MM_04_BM, cells = select.cells_MM_04_BM$x)
MM_04_BM[["Source"]] <- "BM"
MM_04_BM[["Sample_ID"]] <- "MM_04_BM"
MM_04_BM[["ID"]] <- "MM_04"

# MM_04_CTC  ---------------------------------------------------------------
MM_04_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_04_CTC", "filtered_feature_bc_matrix"))
MM_04_CTC <- MM_04_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_04_CTC), perl = TRUE),]
MM_04_CTC <- CreateSeuratObject(MM_04_CTC, project = 'Plasmacells')
select.cells_MM_04_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_04_CTC.csv"), header = TRUE, sep = ',')
MM_04_CTC <- subset(MM_04_CTC, cells = select.cells_MM_04_CTC$x)
MM_04_CTC[["Source"]] <- "CTC"
MM_04_CTC[["Sample_ID"]] <- "MM_04_CTC"
MM_04_CTC[["ID"]] <- "MM_04"

# MM_05_BM  ---------------------------------------------------------------
MM_05_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_05_BM", "filtered_feature_bc_matrix"))
MM_05_BM <- MM_05_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_05_BM), perl = TRUE),]
MM_05_BM <- CreateSeuratObject(MM_05_BM, project = 'Plasmacells')
select.cells_MM_05_BM<- fread(file.path(input_dir, "barcodes", "barcodes_MM_05_BM.csv"), header = TRUE, sep = ',')
MM_05_BM <- subset(MM_05_BM, cells = select.cells_MM_05_BM$x)
MM_05_BM[["Source"]] <- "BM"
MM_05_BM[["Sample_ID"]] <- "MM_05_BM"
MM_05_BM[["ID"]] <- "MM_05"

# MM_05_CTC  ---------------------------------------------------------------
MM_05_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_05_CTC", "filtered_feature_bc_matrix"))
MM_05_CTC <- MM_05_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_05_CTC), perl = TRUE),]
MM_05_CTC <- CreateSeuratObject(MM_05_CTC, project = 'Plasmacells')
select.cells_MM_05_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_05_CTC.csv"), header = TRUE, sep = ',')
MM_05_CTC <- subset(MM_05_CTC, cells = select.cells_MM_05_CTC$x)
MM_05_CTC[["Source"]] <- "CTC"
MM_05_CTC[["Sample_ID"]] <- "MM_05_CTC"
MM_05_CTC[["ID"]] <- "MM_05"

# MM_06_BM  ---------------------------------------------------------------
MM_06_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_06_BM", "filtered_feature_bc_matrix"))
MM_06_BM <- MM_06_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_06_BM), perl = TRUE),]
MM_06_BM <- CreateSeuratObject(MM_06_BM, project = 'Plasmacells')
select.cells_MM_06_BM<- fread(file.path(input_dir, "barcodes", "barcodes_MM_06_BM.csv"), header = TRUE, sep = ',')
MM_06_BM <- subset(MM_06_BM, cells = select.cells_MM_06_BM$x)
MM_06_BM[["Source"]] <- "BM"
MM_06_BM[["Sample_ID"]] <- "MM_06_BM"
MM_06_BM[["ID"]] <- "MM_06"

# MM_06_CTC  ---------------------------------------------------------------
MM_06_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_06_CTC", "filtered_feature_bc_matrix"))
MM_06_CTC <- MM_06_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_06_CTC), perl = TRUE),]
MM_06_CTC <- CreateSeuratObject(MM_06_CTC, project = 'Plasmacells')
select.cells_MM_06_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_06_CTC.csv"), header = TRUE, sep = ',')
MM_06_CTC <- subset(MM_06_CTC, cells = select.cells_MM_06_CTC$x)
MM_06_CTC[["Source"]] <- "CTC"
MM_06_CTC[["Sample_ID"]] <- "MM_06_CTC"
MM_06_CTC[["ID"]] <- "MM_06"

# MM_07_BM  ---------------------------------------------------------------
MM_07_BM <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_07_BM", "filtered_feature_bc_matrix"))
MM_07_BM <- MM_07_BM[!grepl('^IGH|^IGK|^IGL', rownames(MM_07_BM), perl = TRUE),]
MM_07_BM <- CreateSeuratObject(MM_07_BM, project = 'Plasmacells')
select.cells_MM_07_BM<- fread(file.path(input_dir, "barcodes", "barcodes_MM_07_BM.csv"), header = TRUE, sep = ',')
MM_07_BM <- subset(MM_07_BM, cells = select.cells_MM_07_BM$x)
MM_07_BM[["Source"]] <- "BM"
MM_07_BM[["Sample_ID"]] <- "MM_07_BM"
MM_07_BM[["ID"]] <- "MM_07"

# MM_07_CTC  ---------------------------------------------------------------
MM_07_CTC <- Read10X(data.dir = file.path(input_dir,  "count files", "MM_07_CTC", "filtered_feature_bc_matrix"))
MM_07_CTC <- MM_07_CTC[!grepl('^IGH|^IGK|^IGL', rownames(MM_07_CTC), perl = TRUE),]
MM_07_CTC <- CreateSeuratObject(MM_07_CTC, project = 'Plasmacells')
select.cells_MM_07_CTC <- fread(file.path(input_dir, "barcodes", "barcodes_MM_07_CTC.csv"), header = TRUE, sep = ',')
MM_07_CTC <- subset(MM_07_CTC, cells = select.cells_MM_07_CTC$x)
MM_07_CTC[["Source"]] <- "CTC"
MM_07_CTC[["Sample_ID"]] <- "MM_07_CTC"
MM_07_CTC[["ID"]] <- "MM_07"



# merge -------------------------------------------------------------------

BMvsCTC <- merge(x = MM_01_BM, 
             y = list(MM_01_CTC, MM_02_BM, MM_02_CTC, MM_03_BM, MM_03_CTC,
                      MM_04_BM, MM_04_CTC, MM_05_BM, MM_05_CTC, MM_06_BM, MM_06_CTC, 
                      MM_07_BM, MM_07_CTC),
             merge.data = TRUE)


#  SCT transform, PCA< UMAP, clustering-------------------------------------------------------------------------

BMvsCTC <- SCTransform(BMvsCTC, method = "glmGamPoi", verbose = FALSE, variable.features.n = 3000)

BMvsCTC <- RunPCA(BMvsCTC, npcs = 30, verbose = FALSE)

BMvsCTC <- RunUMAP(BMvsCTC, dims = 1:30, reduction = 'pca', verbose = FALSE)

BMvsCTC <- FindNeighbors(BMvsCTC, dims = 1:30, verbose = FALSE)

BMvsCTC <- FindClusters(BMvsCTC, resolution = 0.5, verbose = FALSE)




