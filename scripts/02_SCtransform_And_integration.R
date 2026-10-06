# ============================================================
# 02. SCTransform normalization and dataset integration
# ============================================================
#
# Purpose:
# - Assess cell-cycle effects
# - Assess potential batch effects
# - Perform SCTransform normalization
# - Regress selected technical/cell-cycle covariates
# - Integrate experimental groups
# - Perform PCA, UMAP and clustering
# - Evaluate clustering resolution using clustree
# - Assess batch effects after integration
#
# Raw data and intermediate Seurat objects are not included
# in this repository. See README for data availability.
# Expected project structure:
# data/
# results/
# figures/

# ============================================================


# Set the working directory to the project root before running,
# or use relative paths throughout the script.

#---------------------------------------------
# Load libraries
#---------------------------------------------

suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
suppressMessages(require(sctransform))
library(clustree)
library('irlba')
library(Matrix)

options(Seurat.object.assay.version = "v3")
load(file = "results/seurat_filtered_fom_Raw_0.15Mitocondria_nGeneLimit_doubletRemoved_20N.RData")
#------------cell cycle Analysis ------------------------#
#To find the contribution of the cell cycle genes to PCA
###To visualize cells based on cell cycle gene 
#####For more information on how to get cell cycle marker see Cell cycle script.r
cell_cycle_markers <- read.csv("data/cell_cycle_makers.csv")

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes, 
                                 set.ident = TRUE)
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
# Scale the counts , this step take sometime
seurat_phase <- ScaleData(seurat_phase)

#To find the contribution of the cell cycle genes to PCA
seurat_phase  <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), ndims.print = 1:10, #PCs to print genes for
                        nfeatures.print = 10)

seurat_phase <- RunPCA(seurat_phase, npcs = 30, features = c(s_genes, g2m_genes))

# Plot the PCA colored by cell cycle phase

#type="cairo" when running on the cluster
par(mfrow = c(1, 2))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")


DimPlot(seurat_phase)

# view cell cycle scores and phase assignments
head(seurat_phase[[]])

seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score

seurat_phase2 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"))
# cell cycle effects strongly mitigated in PCA
#seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase)) #This command took a long time to run

seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("CC.Difference"))

DimPlot(seurat_phase2)

#plot after cell cycle genes
DimPlot(seurat_phase2,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")


DimPlot(seurat_phase2,
        reduction = "pca",
        group.by= "seq_folder")

#plotting after regressin cell cycle difference

DimPlot(seurat_phase3,
        reduction = "pca",
        group.by= "seq_folder")


DimPlot(seurat_phase3)

DimPlot(seurat_phase3,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

par(mfrow =c(1, 2))
a <- DimPlot(seurat_phase,
             reduction = "pca",
             group.by= "seq_folder", 
             split.by= "seq_folder")
b <- DimPlot(seurat_phase, 
             reduction = "pca",
             group.by= "seq_folder")
plot_grid(a, b)

rm(seurat_phase3)
rm(seurat_phase2)

#----------------See if there is a batch effet ------------
###New ,Use this 
seurat_batch <- NormalizeData(filtered_seurat, verbose = FALSE)
seurat_batch <- FindVariableFeatures(seurat_batch, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Run the standard workflow for visualization and clustering
seurat_batch <- ScaleData(seurat_batch , verbose = FALSE)
seurat_batch <- RunPCA(seurat_batch, npcs = 30, verbose = FALSE)
seurat_batch <- RunUMAP(seurat_batch, reduction = "pca", dims = 1:30)
seurat_batch <- RunTSNE(seurat_batch, reduction = "pca", dims = 1:30)


#BactEffect_before

DimPlot(seurat_batch, reduction = "umap", group.by = "seq_folder", split.by = "seq_folder" )

DimPlot(seurat_batch, reduction = "umap", group.by = "seq_folder")

DimPlot(seurat_batch, reduction = "pca", group.by = "seq_folder")

DimPlot(seurat_batch, reduction = "tsne", group.by = "seq_folder")
rm(seurat_batch)

# ----------------SCTransform and Integration ---------------------------#
#SCTransform
##Based on the BATCH WE NEED TO DO integration
#Integration in Seurat is  designed to mitigate batch effects in single-cell RNA-seq (scRNAseq) datasets.
#Remove batch effect with integration 
#Note: for the checking the data and batch effect  use nfeature 2000, but for the analysis use nfeature 3000. 

#Note  with SCTransform, the more PCs we choose the more variation is accounted for when performing the clustering, but it take a lot longer to perform the clustering. Therefore for this analysis, we will use the first 40 PCs to generate the cluster

# Forty PCs were used for downstream dimensionality reduction and clustering
# based on evaluation of the PCA variance structure.

#------------------- SCTransform -----------------------
#increases the memory limit
options(future.globals.maxSize = 4000 * 1024^2)
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("Floxed", "WT")]

#######
#Normalization (NormalizeData) – Optional when using SCTransform, but fine to do here before scoring the cell cycle.

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE) #display progress bar for normalization procedure
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]]$CC.Difference <- split_seurat[[i]]$S.Score - split_seurat[[i]]$G2M.Score
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("CC.Difference", "mitoRatio","nUMI"), variable.features.n=3000)
}
#this command is required to calculated the CC difference
#seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score
#variable.features.n the default is 3000. 
# Check which assays are stored in objects
#----------------Integration --------------------#
# Select the most variable features to use for integration
#Integrate samples using shared highly variable genes
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
#Get Negative Binomial regression parameters per gene

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

#Now, perform CCA (Perform canonical correlation analysis (CCA),
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
#Finally, integrate across conditions.
# Integrate across conditions

seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT") #this step requires around 1-1.5 hr on this data set
#on the cluster it needs 1/2 hour 
#This would often be a good place to save the R object.                                  
# Save integrated seurat object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")

rm(split_seurat)


######The commands below is after running the PCA in the UMAP_Rscript

#seurat_integrated <- readRDS(file ="results/Analysis_V1/integrated_seurat_Smad4.rds") # this object doesn't have the PCA

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

DimHeatmap(seurat_integrated, 
           dims = 1:12, 
           cells = 500, 
           balanced = TRUE)


# Plot PCA
#PCAPlot(seurat_integrated,
#split.by = "sample")

PCAPlot(seurat_integrated,
        split.by = "sample",
        group.by = "sample")

PCAPlot(seurat_integrated) 


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


DimPlot(seurat_integrated,reduction = "umap" ) 



DimPlot(seurat_integrated, split.by = "sample", group.by = "sample",reduction = "umap" ) 

######TSNE if required########
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


DimPlot(seurat_integrated, reduction = "tsne") 


DimPlot(seurat_integrated, split.by = "sample", reduction = "tsne", group.by = "sample") 

# Explore heatmap of PCs

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

ElbowPlot(object = seurat_integrated, 
          ndims = 40)
# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
write.csv(co1, "results/PCVariation_cumulative90.csv")

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
write.csv(co2, "results/PCVariation0_1.csv")
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  #theme_bw() #will make grid back ground
  theme_classic()

#Cluster the cells

#Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’ [Seurat - Guided Clustering Tutorial].

#We will use the FindClusters() function to perform the graph-based clustering. The resolution is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets. For my data set, I am having between arount 8000 cells/ sample. I will use slightly higher (0.6- 3)

#The FindClusters() function allows us to enter a series of resolutions and will calculate the "granularity" of the clustering. This is very helpful for testing which resolution works for moving forward without having to run the function for each resolution.

# Determine the K-nearest neighbor graph

seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Use the first 40 PCs to construct the cell-cell nearest-neighbor graph.

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2, 2.2))

#---- See how many cluster per resuloution ########
# Assign identity of clusters
#e.g integrated_snn_res.0.8 use the command below
#you can type integrated_snn_res.2 to view clusters based on res 2
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP that is calculated early
#Seurat_integrated_resolution_0.8 

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)


# Assign identity of clusters


Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

#Seurat_integrated_resolution0.4.

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)


saveRDS(seurat_integrated, "results/integrated_seurat_Resolution.rds")


#Single-cell RNA-seq clustering analysis lesson 8

DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample", reduction = "umap")  + NoLegend()

# Explore whether clusters segregate by cell cycle phase


DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "Phase", reduction = "umap")  + NoLegend()


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order  = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


#The sort.cell argument will plot the positive cells above the negative cells, while the min.cutoff argument will determine the threshold for shading. A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene will not exhibit any purple shading (completely gray).
#Batch_effect_after_integration_UMAP_split

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder",
        split.by = "seq_folder")

#Batch_effect_after_integration_UMAP.
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder")


#Batch_effect_after_integration_PCA_split
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder",
        split.by = "seq_folder")

#Batch_effect_after_integration_PCA
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder")


#Batch_effect_after_integration_tsne_split
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder",
        split.by = "seq_folder")
#Batch_effect_after_integration_tsne
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder")


#--------------------------------------------
#Identify the best resolution using clustree
#---------------------------------------
library(clustree)


clustree(seurat_integrated, prefix = "integrated_snn_res.")

# To visulaize the stem cell clusters Lgr5 

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Lgr5", node_colour_aggr = "median")

png("figures/ClusterTree_CD3D_Recluster_conf.png",res=300,width = 15, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd3d", node_colour_aggr = "median")

dev.off()



 saveRDS(seurat_integrated, "results/integrated_seurat_Resolution.rds")

sessionInfo()

  
  
  
