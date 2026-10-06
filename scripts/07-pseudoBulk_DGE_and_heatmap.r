#-----------------------------------
# Single-cell DGE using MAST at resolution 0
# Resolution-0 single-cell differential expression using MAST, comparing WT versus Floxed cells.
#-----------------------------------

#-----------------------------------
#Load libraries
#-----------------------------------
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
suppressMessages(require(Matrix))
suppressMessages(require(scCustomize))
suppressMessages(require(viridis))
suppressMessages(require("MAST"))
#-------------------------------

# Set the working directory to the project root before running,
# or use relative paths throughout the script.

options(Seurat.object.assay.version = "v3")
annotations <- read.csv("~/Documents/Analysis/data/annotations.csv")

# Load Seurat object
seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")


#-----------------------------------
# Set cluster identity to resolution 0
#-----------------------------------


Idents(object = seurat_integrated) <- "integrated_snn_res.0" 

#-----------------------------------
# Create cluster + sample identity
#-----------------------------------


seurat_integrated@meta.data$clusters_condition <- paste(as.vector(seurat_integrated@meta.data$integrated_snn_res.0), as.vector(seurat_integrated@meta.data$sample), sep="_") #creat a new column

#table (seurat_integrated@meta.data$clusters_condition)
#Idents(seurat_integrated) <- seurat_integrated$clusters_condition #set this as the cell identity 

## #####this will add the data on the results folder but for each file alone so better to add the data to different folder

output <- "results/MAST_DGE/DGE_resolution_0_cluster_Alldata_Nogrouping_"


Idents(object = seurat_integrated) <- "clusters_condition"

#-----------------------------------
# DGE: WT versus Floxed
#-----------------------------------


for (i in 0){ #or however many clusters you have
  ident1 <- paste0(i,"_WT")  
  ident2 <- paste0(i,"_Floxed") 
  condition.diffgenes <- FindMarkers(seurat_integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,assay = 'RNA',
                                     test.use ="MAST" ) 
  write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
}

#-----------------------------------
# Read DGE results
#-----------------------------------


degenes <- read.csv("results/MAST_DGE/DGE_resolution_0_cluster_Alldata_Nogrouping_0.csv")

#-----------------------------------
# Prepare fold-change columns
#-----------------------------------

degenes$log2FoldChange <- degenes$avg_log2FC

names(degenes)[names(degenes) == "p_val"] <- "pvalue"
names(degenes)[names(degenes) == "p_val_adj"] <- "padj"
names(degenes)[names(degenes) == "X"] <- "gene"

names(degenes)[names(degenes)=="p_val"] <- "pvalue"
names(degenes)[names(degenes)=="p_val_adj"] <- "padj"
names(degenes)[names(degenes)=="X"] <- "gene"

#-----------------------------------
# Select top differentially expressed genes
#-----------------------------------

#-----------------------------------
# Select top differentially expressed genes
#-----------------------------------

res_top25 <- degenes %>%
  slice_max(
    order_by = log2FoldChange,
    n = 25,
    with_ties = FALSE
)

res_lower25 <- degenes %>%
  slice_min(
    order_by = log2FoldChange,
    n = 25,
    with_ties = FALSE
)

top50DEG <- rbind(
  res_top25,
  res_lower25
)


res_top50 <- degenes %>%
  slice_max(
    order_by = log2FoldChange,
    n = 50,
    with_ties = FALSE
)

res_lower50 <- degenes %>%
  slice_min(
    order_by = log2FoldChange,
    n = 50,
    with_ties = FALSE
)

top100DEG <- rbind(
  res_top50,
  res_lower50
)

#-----------------------------------
# Heatmap
#-----------------------------------


DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.0"

seurat_integrated <- ScaleData(object = seurat_integrated)


png("figures/Heatmap_top50_DEG_ALl_RESOLUTION_0.png", res=300, width = 4, height =8  , units = "in") 

DoHeatmap(seurat_integrated, features = as.character(top50DEG$gene), group.by = "sample") 

dev.off()

#-----------------------------------
# Violin plot of genes from DGE
#-----------------------------------

png("figures/dittoViloinPlot_Pak3_Pseudobulk.png", res=300, width = 3.2, units = "in",  height = 1.8)
dittoPlot(seurat_integrated, "Pak3", group.by = "group", jitter.size =0.00001, 
jitter.width = 0.075, color.panel = c("black", "gray85"), 
vlnplot.lineweight = 0.3, vlnplot.width = 1)+
  theme_bw(base_size = 8) + ###to change the font on x and y axis, will add the boundaries with a grid
   theme(panel.grid = element_blank()) ###Remove the grid generated with theme_bw()

dev.off()
