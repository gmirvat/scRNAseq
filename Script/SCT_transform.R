#title: "SCT_integration_analysis"
#author: "Mirvat"
#date: "21/11/2023"
#Submitted batch job 5311005
###In this script, I did regress mitochondrial genes. I only regress the diffece in cell cycle genes. 
#I split the seurat and then peformed the analysis. 
#the Normalization and the scaling were OK in the time required for command excution, however, the regress 
#I am using the raw Counts from the 10X genomics
#repeat the analysis as in V3  exactly

#module load R/default
setwd("/mnt/beegfs/hassan/mirvat/Smad4-SCRNAseq/SeuratAnalysis")
load(file = "results/seurat_filtered_fom_Raw_0.15Mitocondria_nGeneLimit_20N.RData")



# Load libraries
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

cell_cycle_makers <- read.csv("/mnt/beegfs/hassan/mirvat/Smad4-SCRNAseq/SeuratAnalysis/data/cell_cycle_makers.csv")

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
seurat_phase  <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), ndims.print = 6:10, #PCs to print genes for
                        nfeatures.print = 10)


# Plot the PCA colored by cell cycle phase

#seurat_phase <- RunPCA(seurat_phase) #for all variation

seurat_phase <- RunPCA(seurat_phase, features = c(s_genes, g2m_genes))

# Plot the PCA colored by cell cycle phase


png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_before.png", res = 300, height = 4, width =8, units ="in", type="cairo")       
#type="cairo" when running on the cluster
par(mfrow = c(1, 2))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
DimPlot(seurat_phase)

dev.off()

# view cell cycle scores and phase assignments
head(seurat_phase[[]])

seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score

seurat_phase2 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"))
# cell cycle effects strongly mitigated in PCA
#seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase)) #This command took a long time to run

seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("CC.Difference"))

png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after.png", res = 300, height = 3, width =4, units ="in", type="cairo")  
DimPlot(seurat_phase2)
dev.off()
#DimPlot(seurat_phase2)
png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after_Pca.png", res = 300, height = 3, width =4, units ="in", type="cairo")  
DimPlot(seurat_phase2,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()
png("figures/Analysis_V1/PCA_Cell_Cycle_differenceCC_after_Pca.png", res = 300, height = 3, width =4, units ="in", type="cairo")  
DimPlot(seurat_phase3,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()



png("figures/Analysis_V1/BactEffect.png", res = 300, height = 3, width =8, units ="in", type="cairo")
par(mfrow =c(1, 2))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "seq_folder", 
        split.by= "seq_folder")
DimPlot(seurat_phase, 
        reduction = "pca",
        group.by= "seq_folder")

dev.off()
##################
###New ,Use this 
seurat_batch <- NormalizeData(filtered_seurat, verbose = FALSE)
seurat_batch <- FindVariableFeatures(seurat_batch, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Run the standard workflow for visualization and clustering
seurat_batch <- ScaleData(seurat_batch , verbose = FALSE)
seurat_batch <- RunPCA(seurat_batch, npcs = 30, verbose = FALSE)
seurat_batch <- RunUMAP(seurat_batch, reduction = "pca", dims = 1:30)
seurat_batch <- RunTSNE(seurat_batch, reduction = "pca", dims = 1:30)

png("figures/Analysis_V1/BactEffect_before_1.png", res = 300, height = 3, width =12, units ="in", type="cairo")
DimPlot(seurat_batch, reduction = "umap", group.by = "seq_folder", split.by = "seq_folder" )
dev.off()

png("figures/Analysis_V1/BactEffect_before_UMAP_2.png", res = 300, height = 3, width =4, units ="in", type="cairo")
DimPlot(seurat_batch, reduction = "umap", group.by = "seq_folder")
dev.off()
png("figures/Analysis_V1/BactEffect_before_PCA_2.png", res = 300, height = 3, width =4, units ="in", type="cairo")
DimPlot(seurat_batch, reduction = "pca", group.by = "seq_folder")
dev.off()
png("figures/Analysis_V1/BactEffect_before_tsne_2.png", res = 300, height = 3, width =4, units ="in", type="cairo")
DimPlot(seurat_batch, reduction = "tsne", group.by = "seq_folder")
dev.off()

####SCTransform
options(future.globals.maxSize = 4000 * 1024^2)
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("Floxed", "WT")]
#Note: I will keep the mitocondrial genes as these are BAL samples and the mitocondrial genes could be important
#for (i in 1:length(split_seurat)) {
#split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE) #display progress bar for normalization procedure
#split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
#split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"))
#}
######
#I need to regress https:/satijalab.org/seurat/articles/cell_cycle_vignette.html
#The procedure above removes all signal associated with cell cycle. In some cases, we’ve found that this can negatively impact downstream analysis, particularly in differentiating processes (like murine hematopoiesis), where stem cells are quiescent and differentiated cells are proliferating (or vice versa). In this case, regressing out all cell cycle effects can blur the distinction between stem and progenitor cells as well.

#As an alternative, we suggest regressing out the difference between the G2M and S phase scores. This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase among proliferating cells (which are often uninteresting), will be regressed out of the data

#######
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE) #display progress bar for normalization procedure
 split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]]$CC.Difference <- split_seurat[[i]]$S.Score - split_seurat[[i]]$G2M.Score
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("CC.Difference", "mitoRatio","nUMI"), variable.features.n=3000)
}

#variable.features.n the default is 3000. 
# Check which assays are stored in objects
#split_seurat$'TBT-04301005b'@assays

# Select the most variable features to use for integration
#Integrate samples using shared highly variable genes

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
#Get Negative Binomial regression parameters per gene
#Using 2000 genes, 4658 cells

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

#Now, we are going to perform CCA (Perform canonical correlation analysis (CCA),
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

#Finally, we can integrate across conditions.
# Integrate across conditions

seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT") #this step requires around 1-1.5 hr on this data set
#on the cluster it needs 1/2 hour 
#This would often be a good place to save the R object.                                  
# Save integrated seurat object


saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_Smad4.rds")

print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)


######The commands below is after running the PCA in the UMAO_Rscript

seurat_integrated <- readRDS(file ="results/Analysis_V1/integrated_seurat_count.rds") # this object doesn't have the PCA

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

png("figures/Analysis_V1/PCA_MostVariableGenes_split.png", res=300,width = 9, height = 4, units = "in", type="cairo" )

DimHeatmap(seurat_integrated, 
           dims = 1:12, 
           cells = 500, 
           balanced = TRUE)
dev.off()

# Plot PCA
#PCAPlot(seurat_integrated,
#split.by = "sample")
png("figures/Analysis_V1/PCA_After_Integration_split.png", res=300,width = 20, height = 4, units = "in", type="cairo" )

PCAPlot(seurat_integrated,
        split.by = "sample",
        group.by = "sample")
dev.off()

png("figures/Analysis_V1/PCA_After_Integration.png", res=300,width = 9, height = 4, units = "in", type="cairo" )

PCAPlot(seurat_integrated) 
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

png("figures/Analysis_V1/UMAP_AfterIntegration_resolution1.png", res=300,width = 9, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated,reduction = "umap" ) 
dev.off()

png("figures/Analysis_V1/UMAP_AfterIntegration_split.png", res=300,width = 9, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated, split.by = "sample", group.by = "sample",reduction = "umap" ) 
dev.off()

####
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
#to acess the cell cycle genes

png("figures/Analysis_V1/PCA_MostVariableGenes_split.png", res=300,width = 9, height = 4, units = "in", type="cairo" )

DimHeatmap(seurat_integrated, 
           dims = 1:12, 
           cells = 500, 
           balanced = TRUE)
dev.off()

# Plot PCA
#PCAPlot(seurat_integrated,
#split.by = "sample")
png("figures/Analysis_V1/PCA_After_Integration_split.png", res=300,width = 20, height = 4, units = "in", type="cairo" )

PCAPlot(seurat_integrated,
        split.by = "sample",
        group.by = "sample")
dev.off()

png("figures/Analysis_V1/PCA_After_Integration.png", res=300,width = 9, height = 4, units = "in", type="cairo" )

PCAPlot(seurat_integrated) 
dev.off()

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

png("figures/Analysis_V1/UMAP_AfterIntegration.png", res=300,width = 9, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated,reduction = "umap" ) 
dev.off()

png("figures/Analysis_V1/UMAP_AfterIntegration_split.png", res=300,width = 9, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated, split.by = "sample", group.by = "sample",reduction = "umap" ) 
dev.off()
######TSNE if required########
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

png("figures/Analysis_V1/tsne_AfterIntegration.png", res=300,width = 9, height = 4, units = "in", type="cairo")  
DimPlot(seurat_integrated, reduction = "tsne") 
dev.off()
png("figures/Analysis_V1/tsne_AfterIntegration_split.png", res=300,width = 9, height = 4, units = "in", type="cairo")  
DimPlot(seurat_integrated, split.by = "sample", reduction = "tsne", group.by = "sample") 
dev.off()
#############

# Explore heatmap of PCs

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

###########
png("figures/Analysis_V1/elbow_Plot_for_the_40PCA.png", res=300, height=5, width=6, units="in", type="cairo")
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()


# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
write.csv(co1, "results/Analysis_V1/PCVariation_cumulative90.csv")

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
write.csv(co2, "results/Analysis_V1/PCVariation0_1.csv")
# Minimum of the two calculation


pcs <- min(co1, co2)
pcs

plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

png("figures/Analysis_V1/elbow_Plot_rank_for_the_40PCA.png", res=300, height=6, width=6, units="in", type="cairo")
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  #theme_bw() #will make grid back ground
  theme_classic()
dev.off()
#Cluster the cells

#Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’ [Seurat - Guided Clustering Tutorial].

#We will use the FindClusters() function to perform the graph-based clustering. The resolution is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets. For my data set, I am having between arount 8000 cells/ sample. I will use slightly higher (0.6- 3)

#The FindClusters() function allows us to enter a series of resolutions and will calculate the "granularity" of the clustering. This is very helpful for testing which resolution works for moving forward without having to run the function for each resolution.

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2, 2.2, 2.4, 2.6, 2.8 , 3, 3.2, 3.4 ,3.6,  4))

######
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP that is calculated early
png("figures/Analysis_V1/Seurat_integrated_resolution_0.8.png",res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

# Assign identity of clusters


Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

png("figures/Analysis_V1/Seurat_integrated_resolution_2.4.png",res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_count.rds")

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1.8"

# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_1.8.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.3"

png("figures/Analysis_V1/Seurat_integrated_resolution_3.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.2"

png("figures/Analysis_V1/Seurat_integrated_resolution_2.png",res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1"

png("figures/Analysis_V1/Seurat_integrated_resolution_1.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.3.2"

png("figures/Analysis_V1/Seurat_integrated_resolution_3.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Assign identity of clusters

Idents(object = seurat_integrated) <- "integrated_snn_res.2"

#Single-cell RNA-seq clustering analysis lesson 8 



png("figures/Analysis_V1/UMAP_of_cells_in_each_cluster_by_sample_afterIntergration.png",res=300,width = 12, height = 5, units = "in",type="cairo")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample", reduction = "umap")  + NoLegend()
dev.off()



# Explore whether clusters segregate by cell cycle phase

png("figures/Analysis_V1/UMAP_of_clusters_segregate_by_cell_cycle_phase_afterIntergration.png", res=300,width = 12, height = 5, units = "in", type="cairo")

DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "Phase", reduction = "umap")  + NoLegend()
dev.off()

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.2"

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
png("figures/Analysis_V1/Segregation_of_clusters_by_various_sources_of_uninteresting_variation.png", res=300, width=12, height=10, units = "in", type="cairo")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order  = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()

#The sort.cell argument will plot the positive cells above the negative cells, while the min.cutoff argument will determine the threshold for shading. A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene will not exhibit any purple shading (completely gray).

png("figures/Analysis_V1/Batch_effect_after_integration_UMAP_split.png", res= 300, height =4, width =8, units ="in", type="cairo")

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()
png("figures/Analysis_V1/Batch_effect_after_integration_UMAP.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder")
dev.off()

png("figures/Analysis_V1/Batch_effect_after_integration_PCA_split.png", res= 300, height =4, width =20, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()


png("figures/Analysis_V1/Batch_effect_after_integration_PCA.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()

png("figures/Analysis_V1/Batch_effect_after_integration_tsne_split.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()
png("figures/Analysis_V1/Batch_effect_after_integration_tsne.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder")
dev.off()




##########
library(clustree)
png("figures/Analysis_V1/ClusterTree_4_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo")

clustree(seurat_integrated, prefix = "integrated_snn_res.")
dev.off()

png("figures/Analysis_V1/ClusterTree_CD3D_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd3d", node_colour_aggr = "median")

dev.off()

png("figures/Analysis_V1/ClusterTree_CD14_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd14", node_colour_aggr = "median")

dev.off()

png("figures/Analysis_V1/ClusterTree_CD8A_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd8a", node_colour_aggr = "median")

dev.off()

png("figures/Analysis_V1/ClusterTree_CCR7_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ccr7", node_colour_aggr = "median")

dev.off()


#CD45RA (CD4 em TCell)
png("figures/Analysis_V1/ClusterTree_CD45RA_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ptprc", node_colour_aggr = "median")
dev.off()

#NK Tcell, Ncam1
png("figures/Analysis_V1/ClusterTree_CD56_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ncam1", node_colour_aggr = "median")
dev.off()

#NK Tcell

png("figures/Analysis_V1/ClusterTree_SERPINA1_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = " Serpina1", node_colour_aggr = "median")
dev.off()
#Large Macrophages (CD206)
#MRC1

png("figures/Analysis_V1/ClusterTree_MRC1_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Mrc1", node_colour_aggr = "median")
dev.off()


###############

columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

###Plot UMAP for each PC 
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
png("figures/Analysis_V1/UMAPplot_for_each_of_the_PCs.png", res= 300, height =8, width =12, units ="in", type="cairo")

map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

dev.off()

saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_count.rds")



