
# ============================================================
# 03. Cluster marker identification and cell composition
# ============================================================
#
# Purpose:
# - Identify marker genes for each cluster
# - Visualize cluster marker expression
# - Quantify the number and proportion of cells per sample
#
# Clustering resolution: 1.2
# ============================================================


# Run this script from the project root directory.
# Input and output files are referenced using relative paths.

#-----------------------------------
#Load libraries
#-------------------------------
suppressMessages(require(Seurat))
suppressMessages(require(dplyr))
library(Matrix)

#remotes::install_github("rnabioco/clustifyrdata")
library(clustifyrdata)
options(Seurat.object.assay.version = "v3")

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

#seurat_integrated 
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

# Identify positively expressed marker genes for each cluster.
# Markers are identified using the Wilcoxon rank-sum test.
# Only positive markers with log2 fold-change >= 0.25 are retained.


markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25, test.use = "wilcox")                   
write.csv(markers, "results/All_maker_identification_per_cluster_resolution_1.2.csv", row.names = F)

#----------------------------------------
#Visualize with heatmap
#----------------------------------
markers <- read.csv("results/All_maker_identification_per_cluster_resolution_1.2.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() 

# Visualize selected cluster markers using normalized RNA expression.

DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated<- ScaleData(object = seurat_integrated)
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

png("figures/Heatmap_Res1.2.png", res =300, height = 15, width=17, units = "in")

DoHeatmap(seurat_integrated, features = top5$gene) + NoLegend()
#remove \+ NoLegend() to add the legend
dev.off()
#-----------------------------------------------
#Identify the number of cells per cluster
#----------------------------------------------------

# Calculate the number and proportion of cells in each cluster
# for each sample.

n_cells2 <- FetchData(seurat_integrated, 
                          vars = c("ident", "sample")) %>% 
dplyr::count(ident,  sample)
    
n_cell.p <- n_cells2 %>% group_by(ident) %>% summarize(total_count = sum(n))
    
df <- left_join(n_cells2 , n_cell.p,  by = "ident") %>%
      mutate(proportion = n / total_count)

write.csv(df, "results/number_cells_sampleOnly_cluster_1.2_resolution.csv", row.names = FALSE)

sessionInfo()

      
