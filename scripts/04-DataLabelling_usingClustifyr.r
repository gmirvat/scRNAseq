
# https://github.com/rnabioco/clustifyr
#remotes::install_github("rnabioco/clustifyrdata")
#https://rnabioco.github.io/clustifyrdata/articles/download_refs.html
# to use this label transfer, please check the available datasets to use based on your own experiment
library(clustifyrdata)
library(clustifyr)

#After picking the best resolution
#use it to predict cell type 
# use of ref_MCA Mouse Cell Atlas
#some databases are comming with the library(clustifyr)

seurat_integrated <- clustify(
  input = seurat_integrated,
  cluster_col = "integrated_snn_res.1.2",
  ref_mat = ref_MCA,
  seurat_out = TRUE
)

metadata <- seurat_integrated@meta.data
#the type is the column with the predicated cell type

pred.CB <- as.data.frame (table(seurat_integrated@meta.data$type, seurat_integrated@meta.data$integrated_snn_res.1.2))
#This only gives on prediction per cluster and if not predicting any thing it will be 0

pred.CB1 <- pred.CB  %>% group_by (Var1) %>% filter(Freq >0)
names(pred.CB1)[names(pred.CB1)=="Var1"] <- "cell_type_CB"
names(pred.CB1)[names(pred.CB1)=="Var2"] <- "Cluster_id"

pred.CB1 <- pred.CB1 %>% group_by(Cluster_id) %>%  arrange(desc (Freq))

#if the expected cell types as the predicted, then proceed below, otherwise check other dataset or use the topmarker percluster to identify the cell type 
#Name the clusters based on clustify R 
#or you can keep type if it is the only one used, is for clustify
seurat_integrated$seuratClu.Res.1.2 <- seurat_integrated@meta.data$type

###plotting based on classifyR 
#to plot 
Idents(object = seurat_integrated) <- "seuratClu.Res.1.2"
  DimPlot(object = seurat_integrated, 
         reduction = "umap", 
         label = TRUE,
         label.size = 5,
         repel = TRUE)


 




