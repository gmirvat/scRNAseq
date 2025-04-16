#this to plot few figues after the second one in which I did not regress on cell cycle genes
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
suppressMessages(require(sctransform))
library(clustree)
library('irlba')
library(Matrix)
setwd("/Users/mirvat/Documents/Bioinformatics/March13_2023/Smad4_Single_Cell/Analysis/")
options(Seurat.object.assay.version = "v3")


seurat_integrated <- readRDS("results/Analysis_V2/integrated_seurat_Resolution_pc20_noCCRegress.rds")

library(clustree)
png("figures/Analysis_V2/ClusterTree_4_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo")

clustree(seurat_integrated, prefix = "integrated_snn_res.")
dev.off()

png("figures/Analysis_V2/ClusterTree_CD3D_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd3d", node_colour_aggr = "median")

dev.off()


# Assign identity of clusters
Idents(object = seurat_integrated) <-  "integrated_snn_res.1"

# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_1.png", res=300,width = 4, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

png("figures/Analysis_V2/Seurat_integrated_resolution_1.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4, split.by = "sample")
dev.off()


Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_1.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_1.2_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
dev.off()


Idents(object = seurat_integrated) <- "integrated_snn_res.2.2"

# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_2.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_2.2_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
dev.off()


Idents(object = seurat_integrated) <- "integrated_snn_res.2"

# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_2.0.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Plot the UMAP
png("figures/Analysis_V2/Seurat_integrated_resolution_2.0_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
dev.off()
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.3"

png("figures/Analysis_V2/Seurat_integrated_resolution_3.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.2"

png("figures/Analysis_V2/Seurat_integrated_resolution_2.png",res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1"

png("figures/Analysis_V2/Seurat_integrated_resolution_1.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.3.2"

png("figures/Analysis_V2/Seurat_integrated_resolution_3.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Assign identity of clusters

Idents(object = seurat_integrated) <- "integrated_snn_res.2"

#Single-cell RNA-seq clustering analysis lesson 8 



png("figures/Analysis_V2/UMAP_of_cells_in_each_cluster_by_sample_afterIntergration.png",res=300,width = 12, height = 5, units = "in",type="cairo")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample", reduction = "umap")  + NoLegend()
dev.off()



# Explore whether clusters segregate by cell cycle phase
png("figures/Analysis_V2/UMAP_of_clusters_segregate_by_cell_cycle_phase_afterIntergration.png", res=300,width = 12, height = 5, units = "in", type="cairo")

DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "Phase", reduction = "umap")  + NoLegend()
dev.off()

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.2"

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
png("figures/Analysis_V2/Segregation_of_clusters_by_various_sources_of_uninteresting_variation.png", res=300, width=12, height=10, units = "in", type="cairo")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order  = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()

#The sort.cell argument will plot the positive cells above the negative cells, while the min.cutoff argument will determine the threshold for shading. A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene will not exhibit any purple shading (completely gray).

png("figures/Analysis_V2/Batch_effect_after_integration_UMAP_split.png", res= 300, height =4, width =8, units ="in", type="cairo")

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()
png("figures/Analysis_V2/Batch_effect_after_integration_UMAP.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by= "seq_folder")
dev.off()

png("figures/Analysis_V2/Batch_effect_after_integration_PCA_split.png", res= 300, height =4, width =20, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()


png("figures/Analysis_V2/Batch_effect_after_integration_PCA.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()

png("figures/Analysis_V2/Batch_effect_after_integration_tsne_split.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder",
        split.by = "seq_folder")
dev.off()
png("figures/Analysis_V2/Batch_effect_after_integration_tsne.png", res= 300, height =4, width =8, units ="in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "tsne",
        group.by= "seq_folder")
dev.off()




##########
library(clustree)
png("figures/Analysis_V2/ClusterTree_4_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo")

clustree(seurat_integrated, prefix = "integrated_snn_res.")
dev.off()
######1.2 is the most stable 

Idents(object = seurat_integrated) <- "integrated_snn_res1.2"

png("figures/Analysis_V2/Seurat_integrated_resolution_1.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()
#Single-cell RNA-seq clustering analysis lesson 8 



png("figures/Analysis_V2/UMAP_of_cells_in_each_cluster_by_sample_afterIntergration_resolution1.2.png",res=300,width = 12, height = 5, units = "in",type="cairo")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample", reduction = "umap")  + NoLegend()
dev.off()


png("figures/Analysis_V2/ClusterTree_CD3D_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd3d", node_colour_aggr = "median")

dev.off()

png("figures/Analysis_V2/ClusterTree_CD14_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd14", node_colour_aggr = "median")

dev.off()
################
png("figures/Analysis_V2/ClusterTree_CD8A_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd8a", node_colour_aggr = "median")

dev.off()



png("figures/Analysis_V2/ClusterTree_CCR7_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ccr7", node_colour_aggr = "median")

dev.off()


#CD45RA (CD4 em TCell)
png("figures/Analysis_V2/ClusterTree_CD45RA_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ptprc", node_colour_aggr = "median")
dev.off()

#NK Tcell, Ncam1
png("figures/Analysis_V2/ClusterTree_CD56_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Ncam1", node_colour_aggr = "median")
dev.off()

#NK Tcell

png("figures/Analysis_V2/ClusterTree_SERPINA1_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = " Serpina1", node_colour_aggr = "median")
dev.off()
#Large Macrophages (CD206)
#MRC1

png("figures/Analysis_V2/ClusterTree_MRC1_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Mrc1", node_colour_aggr = "median")
dev.off()

png("figures/Analysis_V2/ClusterTree_lgr5_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Lgr5", node_colour_aggr = "median")
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
png("figures/Analysis_V2/UMAPplot_for_each_of_the_PCs.png", res= 300, height =8, width =12, units ="in", type="cairo")

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
