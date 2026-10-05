
####Based on the cluster tree I need to identify cluster marker genes. This can be done by FindAllMarkers
#Best resolution is 1.2 
#https://satijalab.org/seurat/articles/seurat5_essential_commands

setwd("~/Documents/Analysis/") # set it to your preferred directory

#-----------------------------------
#Load libraries
#-------------------------------
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
suppressMessages(require(dplyr))
library('irlba')
library(Matrix)

#remotes::install_github("rnabioco/clustifyrdata")
library(clustifyrdata)
options(Seurat.object.assay.version = "v3")

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

#seurat_integrated 
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25, test.use = "wilcox")                   
write.csv(markers, "results/All_maker_identification_per_cluster_resolution_1.2.csv", row.names = F)

#----------------------------------------
#Visualize with heatmap
#----------------------------------
markers <- read.csv("results/All_maker_identification_per_cluster_resolution_1.2_8Jan.csv")
top_markers <- markers %>%  group_by(cluster) %>%  top_n(5)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

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

####To identify the number of cells per cluster per genotype 
n_cells2 <- FetchData(seurat_integrated, 
                          vars = c("ident", "sample")) %>% 
dplyr::count(ident,  sample)
    
n_cell.p <- n_cells2 %>% group_by(ident) %>% summarize(total_count = sum(n))
    
df <- left_join(n_cells2 , n_cell.p,  by = "ident") %>%
      mutate(proportion = n / total_count)

write.csv(df, "results/number_cells_sampleOnly_cluster_1.2_resolution.csv", row.names = FALSE)

sessionInfo()

      