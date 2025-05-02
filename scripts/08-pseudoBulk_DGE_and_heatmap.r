
#This script to perform pseudobulk analysis
#To to this first set set the resolution or cluster identity =0

setwd("~/Documents/Analysis/") # set it to your preferred directory

#-----------------------------------
#Load libraries
# Load libraries
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
setwd("~/Documents/Analysis/") 
options(Seurat.object.assay.version = "v3")
annotations <- read.csv("~/Documents/Analysis/data/annotations.csv")
#To to this first set set the resolution or cluster identity =0

Idents(object = seurat_integrated) <- "integrated_snn_res.0" 

seurat_integrated@meta.data$clusters_condition <- paste(as.vector(seurat_integrated@meta.data$integrated_snn_res.0), as.vector(seurat_integrated@meta.data$sample), sep="_") #creat a new column
#table (seurat_integrated@meta.data$clusters_condition)
#Idents(seurat_integrated) <- seurat_integrated$clusters_condition #set this as the cell identity 

###Working function #####this will put the data on the results folder but for each file alone so better to add the data to different folder
output <- "results/MAST_DGE/DGE_resolution_0_cluster_Alldata_Nogrouping_"


Idents(object = seurat_integrated) <- "clusters_condition"

for (i in 0){ #or however many clusters you have
  ident1 <- paste0(i,"_WT")  
  ident2 <- paste0(i,"_Floxed") 
  condition.diffgenes <- FindMarkers(seurat_integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,assay = 'RNA',
                                     test.use ="MAST" ) 
  write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
}

degenes <- read.csv("results/MAST_DGE/DGE_resolution_0_cluster_Alldata_Nogrouping_0.csv")

degenes$FC <- exp(degenes$avg_log2FC) ###first get the value of the fc
degenes$log2FoldChange <- log2(degenes$FC) #change it to log2
names(degenes)[names(degenes)=="p_val"] <- "pvalue"

names(degenes)[names(degenes)=="p_val_adj"] <- "padj"
names(degenes)[names(degenes)=="X"] <- "gene"

res_top25 <- degenes %>% top_n(n = 25, wt = log2FoldChange) ###3select the top50 DEG
res_lower25 <- degenes %>% top_n(n = 25, wt = - (log2FoldChange))
top50DEG <- rbind(res_top25 ,res_lower25 )

res_top50 <- degenes %>% top_n(n = 50, wt = log2FoldChange) ###3select the top50 DEG
res_lower50 <- degenes %>% top_n(n = 50, wt = -(log2FoldChange))
top100DEG <- rbind(res_top50 ,res_lower50 )


Idents(object = seurat_integrated) <- "integrated_snn_res.0"
DefaultAssay(seurat_integrated) <- "RNA"
DEG <- ScaleData(object = seurat_integrated)

png("figures/Heatmap_top50_DEG_ALl_RESOLUTION_0.png", res=300, res=300, width = 4, height =8  , units = "in") 

DoHeatmap(DEG, features = as.character(top50DEG$gene), group.by = "sample") 

dev.off()

