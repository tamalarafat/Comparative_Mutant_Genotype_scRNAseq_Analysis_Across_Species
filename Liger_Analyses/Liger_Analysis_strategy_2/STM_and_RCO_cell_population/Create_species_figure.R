# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###############################

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_BLS_STM/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

# Lets get the metadata file
md = integrated.data@meta.data

############################ Set the resolution parameter ############################
Idents(integrated.data) <- integrated.data$RNA_snn_res.0.3
Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))
paste("What is the active ident?", paste(levels(integrated.data@active.ident), collapse = ", "))
###########################################################################

# Get the cell ids from STM cell cluster : All cell ids, stm expressing cell ids; AT1G62360 = STM

DefaultAssay(integrated.data) <- "RNA"

# Pick the colors

temp_sps = levels(integrated.data@meta.data[["Species"]])
temp_reps = levels(integrated.data@meta.data[["Replicates"]])

# set the color base
base_col = "#f2edee" # The base color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A") # Manually picked colors for each species class
sps_col = grp_col[c(1:length(temp_sps))]

temp_col_list = list()

# Lets create color gradient for the replicates for each group of species
for (items in c(1:length(temp_sps))){
  temp_col = colorRampPalette(c(base_col, sps_col[items]))(length(grep(temp_sps[items], unique(str_c(integrated.data@meta.data[["Replicates"]], integrated.data@meta.data[["Species"]])))) + 2)[-c(1:2)]
  temp_col_list[[items]] = temp_col
}

rep_col = rapply(temp_col_list, c) # we can also use c(list, recursive = TRUE) or unlist(list = use.names = FALSE)

# For any gene ID, this script will find total number of cells in the cell atlas expressing that gene
# STM = "AT1G62360"

interested_gene = "AT1G62360"

# Let's get the cell count plot
cell_subset = subset(integrated.data, subset =  AT1G62360 > 0)
total_cell_prop = as.data.frame(table(cell_subset@meta.data$Species))
total_cell_prop$cluster = "All cells \nexpressing STM"
total_cell_prop = total_cell_prop[, c(1, 3, 2)]
colnames(total_cell_prop) <- c("species", "cluster", "total_cell_count")

# Lets get the cell proportion visualization
p1 <- ggplot(data = total_cell_prop, aes(x = cluster, y = total_cell_count, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_text(aes(x = cluster, y = total_cell_count, label = total_cell_count),  position = position_dodge(width = 0.8), size = 8, vjust = -0.5, fontface = "bold") +
  xlab("Cell clusters") + ylab("Number of cells") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "#71D0F5FF"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 26, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 26, face = "bold", colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22, face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) + 
  scale_fill_manual(values = sps_col) + guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))

p1

#####
#####
#####

# This will subset the data based on the provided cluster ident and gene ID
Cluster_subset = subset(integrated.data, idents = 13, subset = AT1G62360 > 0)

# This part checks the number of stm expressing cells in cluster 15 where majority of the STM expressing cells are grouped together - visualized per replicate
cell_prop = as.data.frame(table(Cluster_subset@meta.data$Species, Cluster_subset@meta.data$RNA_snn_res.0.3))
colnames(cell_prop) <- c("species", "cluster", "cell_count")
cell_prop = cell_prop[cell_prop$cluster == 13, ]

# Lets get the cell proportion visualization
p2 <- ggplot(data = cell_prop, aes(x = cluster, y = cell_count, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_text(aes(x = cluster, y = cell_count, label = cell_count),  position = position_dodge(width = 0.8), size = 8, vjust = -0.5, fontface = "bold") +
  xlab("Cell clusters") + ylab("Number of cells") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "#71D0F5FF"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 26, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 26, face = "bold", colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22, face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) + 
  scale_fill_manual(values = sps_col) + guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))



#####
#####
#####



# This part checks the number of stm expressing cells in cluster 15 where majority of the STM expressing cells are grouped together - visualized per replicate
count_n_cluster = merge(total_cell_prop[ , c("species", "total_cell_count")], cell_prop, by = "species")

count_n_cluster = count_n_cluster %>% dplyr::select(species, cluster, cell_count, total_cell_count) %>% 
  dplyr::mutate(frac_detected = round((cell_count / total_cell_count) * 100, 1)) %>%
  mutate(frac_label = str_c(frac_detected, "%")) %>%
  mutate(count_label = str_c(cell_count, frac_label, sep = ",\n")) %>%
  mutate_at(vars(species), as.factor)

# Lets get the cell proportion visualization
p3 <- ggplot(data = count_n_cluster, aes(x = cluster, y = frac_detected, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_text(aes(x = cluster, y = frac_detected, label = frac_label),  position = position_dodge(width = 0.8), size = 8, vjust = -0.5, fontface = "bold") + 
  xlab("Cell clusters") + ylab("Number of cells") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "#71D0F5FF"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 26, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 26, face = "bold", colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22, face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) + 
  scale_fill_manual(values = sps_col) + guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))



arranged_fig <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(filename = "Species_STM_expressing_cells_BLS_STM_CH.png", plot = arranged_fig, width = 32, height = 22, dpi = 300, bg = "white")

