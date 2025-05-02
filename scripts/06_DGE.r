
git status#title: "DGE"
#author: "Mirvat"
#date: "2023-12-22"


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

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

DefaultAssay(seurat_integrated) <- "RNA"

#set the resolution or cluster identity you want to do differential gene expression

Idents(object = seurat_integrated) <- "celltype" #set the resolution or cluster identity

seurat_integrated@meta.data$clusters_condition <- paste(as.vector(seurat_integrated@meta.data$celltype), as.vector(seurat_integrated@meta.data$sample), sep="_") #create a new column

Idents(seurat_integrated) <- "seurat_integrated$clusters_condition" #set this as the cell identity 
table(seurat_integrated@meta.data$clusters_condition)
dir.create("results/MAST_DGE/")
output <- "results/MAST_DGE/DGE_resolution_celltype_" 
Idents(seurat_integrated) <- seurat_integrated@meta.data$clusters_condition
#e,g write the cell types you need to do DGE for them 
A = c(  "CAF","Tacstd2 hi", "Fmnl2 hi",  "TA",  "Secretory Pr"  ,   "Cachd1 hi"  ,     "Slc14a1 hi" ,"Slc12a2 hi",      "Lgr5 ISC","Pak3 hi"  ,      "Clu hi"  ,
"EC-1", "TA","T cells","Macrophages","Proliferative", "Paneth", "Neutrophils")
for (i in A){ 
  #for (i in 0:32){ #or however many clusters you have if you want to do it on resolution
  ident1 <- paste0(i,"_WT")  #paste0(i,"_ctrl")
  ident2 <- paste0(i,"_Floxed")  #paste0(i,"_BCG") 
  condition.diffgenes <- FindMarkers(seurat_integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,assay = 'RNA',
                                     test.use ="MAST" )  %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c(gene = "gene_name"))
  write.csv(condition.diffgenes, file = paste0(output, i, ".csv"))
  write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
  
  
  upFl <- condition.diffgenes %>%
    dplyr::filter(avg_log2FC < 0)
  downFl <- condition.diffgenes %>%
    dplyr::filter(avg_log2FC > 0)
  
  write.csv(upFl, file = paste0(output, i, "_upFl.csv"))
  write.csv(downFl, file = paste0(output, i, "_downFl.csv"))
}



sessionInfo()

