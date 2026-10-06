
#-----------------------------------
#06_DGE
#-----------------------------------

#author: "Mirvat"
#date: "2023-12-22"

####### use this after AssignCellTypeAndScproportiopn.r

# Set the working directory to the project root before running,
# or use relative paths throughout the script.

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

# -----------------------------------
# Load data
# -----------------------------------

options(Seurat.object.assay.version = "v3")
annotations <- read.csv("~/Documents/Analysis/data/annotations.csv")

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

DefaultAssay(seurat_integrated) <- "RNA"

# -----------------------------------
# Set cell type identity
# -----------------------------------

#set the resolution or cluster identity you want to do differential gene expression

Idents(object = seurat_integrated) <- "celltype" #set the resolution or cluster identity

# -----------------------------------
# Create celltype + sample identity
# -----------------------------------


seurat_integrated@meta.data$clusters_condition <- paste(as.vector(seurat_integrated@meta.data$celltype), as.vector(seurat_integrated@meta.data$sample), sep="_") #create a new column

Idents(seurat_integrated) <- "seurat_integrated$clusters_condition" #set this as the cell identity 
table(seurat_integrated@meta.data$clusters_condition)

# -----------------------------------
# Create output directory
# -----------------------------------


dir.create(
  "results/MAST_DGE/",
  showWarnings = FALSE,
  recursive = TRUE
)  #safer if the directory doesn't already exist.

output <- "results/MAST_DGE/DGE_resolution_celltype_" 

# -----------------------------------
# Set identity for DGE
# -----------------------------------


Idents(seurat_integrated) <- seurat_integrated@meta.data$clusters_condition
# -----------------------------------
# Cell types to analyse 
# MAST differential expression
# WT versus Floxed within each cell type
# -----------------------------------

A = c(  "CAF","Tacstd2 hi", "Fmnl2 hi",  "TA",  "Secretory Pr"  ,   "Cachd1 hi"  ,     "Slc14a1 hi" ,"Slc12a2 hi",      "Lgr5 ISC","Pak3 hi"  ,      "Clu hi"  ,
"EC-1","T cells","Macrophages","Proliferative", "Paneth", "Neutrophils")
for (i in A){ 
  #for (i in 0:32){ #or however many clusters you have if you want to do it on resolution
  ident1 <- paste0(i,"_WT")  #paste0(i,"_ctrl")
  ident2 <- paste0(i,"_Floxed")  #paste0(i,"_BCG") 
  condition.diffgenes <- FindMarkers(seurat_integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,assay = 'RNA',
                                     test.use ="MAST" )  %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c(gene = "gene_name"))
  
  write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
  
 # -----------------------------------
  # Genes higher in Floxed
  # -----------------------------------
  
  # Because ident.1 = WT and ident.2 = Floxed,
  # avg_log2FC > 0 → higher in WT

#avg_log2FC < 0 → higher in Floxed
   #to identify genes  up in Smad4 Floxed 
  upFl <- condition.diffgenes %>%
    dplyr::filter(avg_log2FC < 0)
  #to identify genes  up in Smad4 WT
upWT <- condition.diffgenes %>%
    dplyr::filter(avg_log2FC > 0)
  
  write.csv(upFl, file = paste0(output, i, "_upFl.csv"))
  write.csv(upWT , file = paste0(output, i, "_downFl.csv"))
}

# -----------------------------------
# Plot top 10 marker genes per cluster in a heatmap
# -----------------------------------


Markers_v1_1.2 <- read.table("Analysis/All_maker_identification_per_cluster_resolution_1.2.csv", sep=",", header=T)
Markers_v1_1.2<-  Markers_v1_1.2 %>%  mutate(Diff=pct.1-pct.2) 

top10_1.2 <- Markers_v1_1.2 %>%
  group_by(cluster) %>%  slice_max(order_by = Diff,
    n = 10, with_ties = FALSE ) %>%  ungroup()




DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "celltype"
seurat_integrated <- ScaleData(seurat_integrated)

jpg("figures/Heatmap_SC_labelled_celltype_top10.jpg", res=300,width = 20, height = 25, units = "in")

DoHeatmap(seurat_integrated , features = top10_1.2$gene, slot = "scale.data")
dev.off ()


sessionInfo()

