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
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.features = 200)

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

# Extract the count matrix
OX_W1 <- GetAssayData(OX_1E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W1) <- paste("O1", colnames(OX_W1), sep = "_")


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_DF_2, project = "OX_2E", min.features = 200)

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

# Extract the count matrix
OX_W2 <- GetAssayData(OX_2E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W2) <- paste("O2", colnames(OX_W2), sep = "_")


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_DF_3, project = "OX_3E", min.features = 200)

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

# Extract the count matrix
OX_W3 <- GetAssayData(OX_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W3) <- paste("O3", colnames(OX_W3), sep = "_")

###
# OX - 7 E
###

# First replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_DF_7, project = "OX_7E", min.features = 200)

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

# Extract the count matrix
OX_W7 <- GetAssayData(OX_7E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W7) <- paste("O7", colnames(OX_W7), sep = "_")


###
# BLS::STM - 1 E
###

# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
STM_1E <- CreateSeuratObject(counts = BLS_1, project = "STM_1E", min.features = 200)

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

# Extract the count table from the seurat object
BLS_S1 <- GetAssayData(STM_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(BLS_S1) <- paste("S1", colnames(BLS_S1), sep = "_")



###
# BLS::STM - 2 E
###

# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
STM_2E <- CreateSeuratObject(counts = BLS_2, project = "STM_2E", min.features = 200)

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

# Extract the count table from the seurat object
BLS_S2 <- GetAssayData(STM_2E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(BLS_S2) <- paste("S2", colnames(BLS_S2), sep = "_")

###
# Create a liger object for OX with all the replicates and select HVGs
###

# Creating liger object
WT_OX <- createLiger(list(WO1 = OX_W1, WO2 = OX_W2, WO3 = OX_W3, WO7 = OX_W7), remove.missing = F)

# Normalization of the data
WT_OX <- rliger::normalize(WT_OX)

# Selecting a set of highly variable genes
WT_OX <- selectGenes(WT_OX, num.genes = 3000, combine = "intersection", do.plot = FALSE) # 1782 genes

fileGenerator(WT_OX@var.genes, fileName = "Shared_highly_variable_genes_between_replicates_CH.txt")

###
# Create a liger object for COL with all the replicates and select HVGs
###

# Creating liger object
WT_COL <- createLiger(list(WC1 = BLS_S1, WC2 = BLS_S2), remove.missing = F)

# Normalization of the data
WT_COL <- rliger::normalize(WT_COL)

# Selecting a set of highly variable genes
WT_COL <- selectGenes(WT_COL, num.genes = 3000, do.plot = FALSE, combine = "intersection")

fileGenerator(WT_COL@var.genes, fileName = "Shared_highly_variable_genes_between_replicates_AT_BLS_STM.txt")

###
# Create a liger object for COL with all the replicates and select HVGs - single dataset of SAM
# We have to pick HVGs manually
###

# Combine the HVGs
HVGs_combined = union(WT_OX@var.genes, WT_COL@var.genes)

fileGenerator(HVGs_combined, "HVG_intersect_reps_union_species_without_mincells.txt")

###
# Creating a Liger object, pre-processing, and performing integration
###

# Lets create the liger object
WT_Species <- createLiger(list(WO1 = OX_W1, WO2 = OX_W2, WO3 = OX_W3, WO7 = OX_W7, BS1 = BLS_S1, BS2 = BLS_S2), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
WT_Species@cell.data$Replicates <- WT_Species@cell.data$dataset
WT_Species@cell.data$Replicates <- factor(WT_Species@cell.data$Replicates, levels = c("WO1", "WO2", "WO3", "WO7", "BS1", "BS2"), labels = c("WT-OX-1", "WT-OX-2", "WT-OX-3", "WT-OX-7", "BLS-STM-1", "BLS-STM-2"))

# Add species information
WT_Species@cell.data$Species <- str_sub(WT_Species@cell.data$dataset, 1, nchar(WT_Species@cell.data$dataset) - 1)
WT_Species@cell.data$Species <- factor(WT_Species@cell.data$Species, levels = c("WO", "BS"), labels = c("Hirsuta", "Thaliana"))

# Add genotype information
WT_Species@cell.data$Genotype <- str_sub(WT_Species@cell.data$dataset, 1, nchar(WT_Species@cell.data$dataset) - 1)
WT_Species@cell.data$Genotype <- factor(WT_Species@cell.data$Genotype, levels = c("WO", "BS"), labels = c("WT", "BLS-STM"))

# Normalization of the data
WT_Species <- rliger::normalize(WT_Species)

# Selecting a set of highly variable genes
WT_Species@var.genes <- HVGs_combined

WT_Species <- scaleNotCenter(WT_Species)

# Check which datasets are we integrating
table(WT_Species@cell.data$dataset)

# Run liger integration - factorization of the matrices
WT_Species <- optimizeALS(WT_Species, k = 50, nrep = 10, lambda = 5)

# Quantile normalization of the data - integration in the shared space
WT_Species <- quantile_norm(WT_Species)

# Run liger implemented UMAP
WT_Species <- runUMAP(WT_Species)

Liger_object_K_50 <- WT_Species

#
save(Liger_object_K_50, file = "integrated_wt_hirsta_stm_liger.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_wt_hirsuta_apex_liger.txt")
