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

# Let's get the cell count plot
cell_subset = subset(integrated.data, subset =  Species == "A.thaliana")

clusterLabels = str_sort(levels(Idents(cell_subset)), numeric = TRUE)

# Lets create a empty vector to store the count information
temp = c()

for (i in c(1:length(clusterLabels))){
  cellNames = WhichCells(cell_subset, idents = clusterLabels[i])
  cellCount = sum(GetAssayData(cell_subset, assay = "RNA", slot = "counts")["AT1G62360", cellNames] != 0)
  temp = c(temp, cellCount)
}

cell_prop = as.data.frame(table(Idents(cell_subset)))
colnames(cell_prop) <- c("cluster", "total_cell")

cell_prop = cell_prop %>% mutate(detected = temp) %>% 
  mutate(non_detected = total_cell - detected) %>% 
  mutate(frac_detected = round((detected / sum(detected)) * 100, 1)) %>% 
  mutate(frac_label = str_c(frac_detected, "%")) %>% 
  mutate(count_label = str_c(detected, frac_label, sep = ",\n"))


long_cell_prop = cell_prop %>% dplyr::select(-c(frac_detected, frac_label, count_label)) %>% 
  gather("detection", "cell_count", detected, non_detected) %>% 
  group_by(cluster) %>% 
  mutate(cluster_prop = round((cell_count/total_cell) * 100, 2)) %>% 
  mutate_at(vars(detection), as.factor)

# frac_detected = In each cluster, the proportion of cells in which expression of the given gene is detected (cells with the genes expression in a cluster / total cells with that genes expression)
# cluster_prop = ration of detected and non detected cells in each cluster

p1 <- ggplot(data = long_cell_prop, aes(x = fct_rev(cluster), y = cell_count)) + 
  geom_bar(stat = "identity", position = "stack", aes(fill = detection), width = 0.6) + 
  geom_text(data = cell_prop, aes(x = cluster, y = (total_cell + 60), label = count_label),  size = 6, fontface = "bold") + 
  coord_flip() + 
  xlab("Cell clusters") + ylab("Number of cells (cells in which expression was detected, proportion of detected cells)") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "#71D0F5FF"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.text = element_text(size = 22, face = "bold", colour = "black"),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 22, face = "bold")) + 
  scale_fill_manual(values = c("#00A08A", "#F98400"), breaks = c("non_detected", "detected"), labels = c("not detected", "detected")) + 
  guides(fill = guide_legend(title = "Expression"))

p1

p2 <- ggplot(data = long_cell_prop, aes(x = fct_rev(cluster), y = cluster_prop)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = detection), width = 0.6) + 
  coord_flip() + 
  xlab("Cell clusters") + ylab("Fraction of cells with detected expression") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "#71D0F5FF"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.text = element_text(size = 22, face = "bold", colour = "black"),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 22, face = "bold")) + 
  scale_fill_manual(values = c("#00A08A", "#F98400"), breaks = c("non_detected", "detected"), labels = c("not detected", "detected")) + 
  guides(fill = guide_legend(title = "Expression"))

arranged_fig <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(filename = "STM_count_n_proportion_BLS_STM.png", plot = arranged_fig, width = 32, height = 22, dpi = 300, bg = "white")


