# ------------------------------------------------------
# 05-AssignCellTypeAndScproportiopn.r
#------------------------------------------------------



#In this script 
#1. Assign cell type based on the top marker genes
#2. UMAP plot
#3. See if the proportion of the cell type differ between the genotype
#4. stack plot for the celll type
#5. e.g dotplots for different cell types, and feature plot for genes 
#Note: the plots are in the figures folder. for publication Use png to save, for this repo I used png or jpg to decrease file size

#author: "Mirvat"
#date: "2024-01-19"
  

# Run this script from the project root directory.
# Input and output files are referenced using relative paths.

# upload integrated_seurat_Resolution.rds
# Load libraries
#---------------------------
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
library('irlba')
library(Matrix)
library(scCustomize)
library(viridis)
library(readxl)
#devtools::install_github("rpolicastro/scProportionTest")
library("scProportionTest")
library(dittoSeq)

options(Seurat.object.assay.version = "v3")
seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

# ------------------------------------------------------------
# Set assay and clustering identity
# ------------------------------------------------------------

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.2.2"

# ------------------------------------------------------------
# Assign cell types based on top marker genes
# ------------------------------------------------------------

# I used Resolution 1.2 for marker and cluster 6 has 2.2 groups
# that I detected from 2.2.
#
# To do the cell label, I will use resolution 2.2.
# I will follow the cluster tree to get the original cluster.
#
# Note: I decided to get 2.2 as this cluster get stable from
# 2.2 onwards.
# Based on 1.2 labelling use the cluster tree
#
# Name the cell types as a celltype, create an empty column
# and fill it

seurat_integrated@meta.data$celltype <-NULL

seurat_integrated <- RenameIdents(
  object = seurat_integrated,

  "8"  = "Tacstd2 hi",
  "7"  = "Tacstd2 hi",
  "14" = "Tacstd2 hi",

  "0"  = "Fmnl2 hi",
  "12" = "Fmnl2 hi",

  "15" = "TA",
  "5"  = "TA",

  "2"  = "Secretory Pr",
  "26" = "Secretory Pr",

  "1"  = "Cachd1 hi",

  "6"  = "Slc14a1 hi",

  "22" = "EClpr",
  "10" = "ECepr",

  "19" = "EC-2",
  "16" = "EC-2",

  "3"  = "Pak3 hi",
  "25" = "Pak3 hi",

  "4"  = "Clu hi",

  "24" = "Slc12a2 hi",
  "17" = "Slc12a2 hi",

  "9"  = "Lgr5 ISC",

  "11" = "EC-1",

  "13" = "TA",

  "21" = "T cells",
  "27" = "T cells",

  "20" = "Macrophages",

  "18" = "Proliferative",

  "23" = "Paneth",

  "28" = "Neutrophils",

  "29" = "CAF"
)


# ISC means Intestinal Stem cell


# ------------------------------------------------------------
# Store cell-type identities in metadata
# ------------------------------------------------------------



seurat_integrated[["celltype"]] <- Idents(object = seurat_integrated)
table(seurat_integrated@meta.data$celltype)

# ------------------------------------------------------------
# Save labelled Seurat object
# ------------------------------------------------------------

saveRDS(seurat_integrated, "results/seurat_integrated_Resolution.rds")

# ------------------------------------------------------------
# Plot cell types by genotype
# ------------------------------------------------------------

Idents(object = seurat_integrated) <- "celltype"
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)

#plot by genotype 

jpg("figures/DimPlotcluster_identity_celltype.jpg",res=300, height = 5, width=3.3, units = "in") #png can also be used, see below

DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE, split.by = "group")

dev.off()

# ------------------------------------------------------------
# Test whether cell-type proportions differ between groups
# ------------------------------------------------------------

# Use scProportionTest


prop_test <- sc_utils(seurat_integrated)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WT", sample_2 = "Floxed", 
  sample_identity = "sample")

Test <- prop_test@results$permutation
write.csv(Test, "StatisticsForTheproprtionofCells_res_celltype.csv", row.names = F)

# ------------------------------------------------------------
# Plot cell-type proportion differences
# ------------------------------------------------------------

png("figures/SC_proportionOfcelltype_log2FD_1.png", res =300, height = 4, width=6, units = "in")
permutation_plot(prop_test,log2FD_threshold = 1)
dev.off()

# ------------------------------------------------------------
# Stacked bar plot of cell-type proportions
# ------------------------------------------------------------

p <- print(dittoBarPlot(seurat_integrated, "celltype", group.by = "group",
                        data.out = TRUE))
png("figures/StackplotPercentageOfCells_Smad4_celltype.png",res=300, height = 5, width=3.3, units = "in")

p
dev.off()

# ------------------------------------------------------------
# Feature plots
# ------------------------------------------------------------


Idents(object = seurat_integrated) <- "celltype"
png("figures/Id1_FeaturePlot.png", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Id1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
dev.off()

png("figures/Pak3_FeaturePlot.png", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Pak3"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
dev.off()

# ------------------------------------------------------------
# Plot markers for the main cell types
# ------------------------------------------------------------

# Plot scaled expression across all cell types



#So we plot the scaled expression across all cell types
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"
#Another marker for marcrophages is , "lyz2"
gene.list  <- c("Epcam","Dcn", "Col1a1", "Cd3e", "Cd3g", "Cd74", "H2-Ab1", "S100a8", "S100a9")

jpg("figures/Main_cell_Markers_b.jpg", res=300, width = 8, units = "in",  height = 8)
FeaturePlot_scCustom(seurat_integrated, features = gene.list, 
        label = FALSE, order = TRUE, reduction = "umap", figure_plot = TRUE, min.cutoff = 'q10')
dev.off()

sessionInfo()
  
  
  
  
