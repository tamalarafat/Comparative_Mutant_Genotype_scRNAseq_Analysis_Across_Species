# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

###
# WT C. hirsuta
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
OX_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_2nd_ALL_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
OX_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the data table to ortho gene ids
OX_DF_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_2 = prepare_ortho_data(input_data = OX_data_2E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_3 = prepare_ortho_data(input_data = OX_data_3E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_7 = prepare_ortho_data(input_data = OX_data_7E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

###
# A. thaliana - BLS::STM
###

# WT AThaliana BLS::STM Experiment 1 - leaf 5 and 6
BLS_1 <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_STM_RNA_1ST_ALL/filtered_feature_bc_matrix/")

# WT AThaliana BLS::STM Experiment 2 - leaf 5 and 6
BLS_2 <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_STM_RNA_2nd_ALL_2/filtered_feature_bc_matrix/")


# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(BLS_1)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# SAM data
BLS_1 <- BLS_1[thaliana_ortho_genes, ]

BLS_2 <- BLS_2[thaliana_ortho_genes, ]

# remove the missing genes from the data
OX_DF_1 <- OX_DF_1[thaliana_ortho_genes, ]
OX_DF_2 <- OX_DF_2[thaliana_ortho_genes, ]
OX_DF_3 <- OX_DF_3[thaliana_ortho_genes, ]
OX_DF_7 <- OX_DF_7[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(BLS_1)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)


##### Remove the protoplasting induced genes
OX_DF_1 <- OX_DF_1[genes_to_keep, ]
OX_DF_2 <- OX_DF_2[genes_to_keep, ]
OX_DF_3 <- OX_DF_3[genes_to_keep, ]
OX_DF_7 <- OX_DF_7[genes_to_keep, ]

# SAM data
BLS_1 <- BLS_1[genes_to_keep, ]

BLS_2 <- BLS_2[genes_to_keep, ]

###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.cells = 13, min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 2000)


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_DF_2, project = "OX_2E", min.cells = 21, min.features = 200)

# Add metadata information to the seurat object
OX_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_2E <- subset(OX_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_2E[["percent.mt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_2E[["percent.pt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_2E <- subset(OX_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_2E <- NormalizeData(OX_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_2E <- FindVariableFeatures(OX_2E, selection.method = "vst", nfeatures = 2000)


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_DF_3, project = "OX_3E", min.cells = 8, min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_3E <- NormalizeData(OX_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_3E <- FindVariableFeatures(OX_3E, selection.method = "vst", nfeatures = 2000)


###
# OX - 7 E
###

# First replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_DF_7, project = "OX_7E", min.cells = 18, min.features = 200)

# Add metadata information to the seurat object
OX_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-7", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_7E <- subset(OX_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_7E[["percent.mt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_7E[["percent.pt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_7E <- subset(OX_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_7E <- NormalizeData(OX_7E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_7E <- FindVariableFeatures(OX_7E, selection.method = "vst", nfeatures = 2000)


###
# BLS::STM - 1 E
###

# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
STM_1E <- CreateSeuratObject(counts = BLS_1, project = "STM_1E", min.cells = 11, min.features = 200)

# Add metadata information to the seurat object
STM_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "BLS-STM-1", "BLS", "Leaf")

# Remove cells with a total count more than 110000
STM_1E <- subset(STM_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
STM_1E[["percent.mt"]] <- PercentageFeatureSet(STM_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
STM_1E[["percent.pt"]] <- PercentageFeatureSet(STM_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
STM_1E <- subset(STM_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
STM_1E <- NormalizeData(STM_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
STM_1E <- FindVariableFeatures(STM_1E, selection.method = "vst", nfeatures = 2000)


###
# BLS::STM - 2 E
###

# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
STM_2E <- CreateSeuratObject(counts = BLS_2, project = "STM_2E", min.cells = 13, min.features = 200)

# Add metadata information to the seurat object
STM_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "BLS-STM-2", "BLS", "Leaf")

# Remove cells with a total count more than 110000
STM_2E <- subset(STM_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
STM_2E[["percent.mt"]] <- PercentageFeatureSet(STM_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
STM_2E[["percent.pt"]] <- PercentageFeatureSet(STM_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
STM_2E <- subset(STM_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
STM_2E <- NormalizeData(STM_2E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
STM_2E <- FindVariableFeatures(STM_2E, selection.method = "vst", nfeatures = 2000)


# Let's create a list to store all the seurat objects
seurat.objects <- list(O1 = OX_1E, O2 = OX_2E, O3 = OX_3E, O7 = OX_7E, S1 = STM_1E, S2 = STM_2E)

# Integration of the replicates - find features from the datasets to anchor cells from different sources
anchFeatures <- SelectIntegrationFeatures(object.list = seurat.objects)

fileGenerator(anchFeatures, "anchor_features_seurat.txt")

ingAnchors <- FindIntegrationAnchors(object.list = seurat.objects, dims = 1:50, anchor.features = anchFeatures)

# To keep records of all the genes in the integrated assay, create a feature set with all of the genes from different replicates.
features_integrated <- unique(c(rownames(OX_1E), rownames(OX_2E), rownames(OX_3E), rownames(OX_7E), rownames(STM_1E), rownames(STM_2E)))

# Integrate the replicates
integrated.data <- IntegrateData(anchorset = ingAnchors, dims = 1:50, verbose = T, features.to.integrate = features_integrated)

# Setting the default assay to "integrated"
DefaultAssay(integrated.data) <- "integrated"

# Gene level scaling - standardization
integrated.data <- ScaleData(integrated.data, verbose = FALSE)

# Run PCA
integrated.data <- RunPCA(integrated.data, npcs = 50, verbose = FALSE)

# Run UMAP and tSNE
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:50, n.components = 2)

integrated.data <- RunTSNE(integrated.data, reduction = "pca", dims = 1:50, dim.embed = 2)

# Find neighbours and clusters
integrated.data <- FindNeighbors(integrated.data, reduction = "pca", dims = 1:50)

for (i in seq(0.1, 1.2, 0.1)) {
  integrated.data <- FindClusters(integrated.data, resolution = i, n.start = 50, n.iter = 50)
}

# Save the integrated file as output
save(integrated.data, file = "integrated_wt_hirsta_stm_seurat.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_wt_hirsuta_stm_seurat.txt")
