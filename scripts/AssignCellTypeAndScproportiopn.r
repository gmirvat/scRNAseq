
#In this script 
#1. Assign cell type based on the top marker genes
#2. UMAP plot
#3. See if the proportion of the cell type differ between the genotype
#4. stack plot for the celll type
#5. e.g dotplots for different cell types, and feature plot for genes 
#Note: the plots are in the figures folder 

#author: "Mirvat"
#date: "2024-01-19"
  

setwd("~/Documents/Analysis/") # set it to your preferred directory

#Just for the analysis I did this I upload integrated_seurat_Resolution.rds
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

options(Seurat.object.assay.version = "v3")
seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.2.2"
#### To add the cell type based on top genes and literature search 
#I used Resolution 1.2 for marker and cluster 6 has 2.2 groups that I detected from 2.2. So to do the cell label, I will use resolution 2.2. I will follow the cluster tree to get the original cluster see figure Analysis_V1/ClusterTree_4_Recluster.png. 
#Note: I decided to get 2.2 as this cluster get stable from 2.2 on wards

####Based on 1.2 labelling  use the cluster tree 
#Name the cell types as a celltype, creat an empty column and fill it 
seurat_integrated@meta.data$celltype <-NULL
Idents(object = seurat_integrated) <- "integrated_snn_res.2.2" 
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "8"	= "Tacstd2 hi",   #"0"	= "ISC-like Tacstd2 hi", #Tacstd2 hi. ISC-L1
                                  "7"	= "Tacstd2 hi",#"0"	= "ISC-like Tacstd2 hi",
                                  "14"	= "Tacstd2 hi",#"0"	= "ISC-like Tacstd2 hi"
                                  "0"="Fmnl2 hi", #"1"="ISC-like Fmnl2 hi", ISC-L2
                                  "12"="Fmnl2 hi", #"1"="ISC-like Fmnl2 hi", ISC-L3 Nav3 the high is associated with worse prognosis
                                  "15"="TA", #lots of Ribosomal genes, "2"="TA", #lots of Ribosomal genes
                                  "5"="TA", #lots of Ribosomal genes, "2"="TA", Slc12a2, 
                                  "2"="Secretory Pr",  #"3"="Secretory Pr", 
                                  "26"="Secretory Pr",  #"3"="Secretory Pr", 
                                  "1"="Cachd1 hi",#"4"="ISC-like Cachd1 hi",  
                                  "6"="Slc14a1 hi", #"5"="ISC-like Slc14a1", SC-L4
                                  "22"="EClpr",#Split "6"="EClpr ECepr",#Split 
                                  "10"="ECepr",#Split ,"6"="EClpr ECepr",#Split 
                                  "19"="EC-2",# "7"="EC-2",
                                  "16"="EC-2", #"7"="EC-2",
                                  "3"="Pak3 hi",#"8"="ISC-like Pak3 hi",
                                  "25"="Pak3 hi",#"8"="ISC-like Pak3 hi",
                                  "4"="Clu hi",#"9"="ISC-like Clu hi", ISC-L5
                                  "24"="Slc12a2 hi",#"10"="ISC-like Kcnq1 hi", #Slc12a2 also ISC-L6
                                  "17"="Slc12a2 hi",#"10"="ISC-like Kcnq1 hi",
                                  "9"="Lgr5 ISC",#"11"="Lgr5 ISC", 
                                  "11"="EC-1", #"12"="EC-1", 
                                  "13"="TA",# "13"="TA",  #CD177 marker and Iftim2
                                  "21"="T cells",# #"14"="T cells", 
                                  "27"="T cells",# #"14"="T cells", 
                                  "20"="Macrophages",#"15"="Macrophages",
                                  "18"="Proliferative", # "16"="Proliferative",
                                  "23"="Paneth", #"17"="Paneth",
                                  "28"="Neutrophils",#"19"="Neutrophils",
                                  "29"="CAF")#"18"="CAF"
#ISC means Intestinal Stem cell 


seurat_integrated[["celltype"]] <- Idents(object = seurat_integrated)
table(seurat_integrated@meta.data$celltype)
#save 
saveRDS(seurat_integrated, "results/seurat_integrated_Resolution.rds")


#Plot 
Idents(object = seurat_integrated) <- "celltype"
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)

#plot by genotype 
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE, split.by = "group")




DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "celltype"
#Dot plot Enterocyte signature 
# Paper: Single-cell RNA Sequencing Reveals How the Aryl Hydrocarbon Receptor Shapes Cellular Differentiation Potency in the Mouse Colon 
features <- c("Krt19", "Krt8", "Lgals3", "Lypd8", "Selenbp1", "Lgals4", "Ethe1", "Fth1", "Gsdmc4", "Maob", "Ces2a", "Ces2c",
              "Ces2e", "Chp2", "Coro2a", "Cyp2c65", "Cyp2c68", "Dgkq", "Dpep1", "Emp1", "Krt20", "Cdkn1a", "Mxd1", "Sprr2a3", "Fam3b",
              "Clec2h", "Fhl1", "Gde1", "Guca2a", "Il18")
  
Clustered_DotPlot(seurat_object = seurat_integrated, features = features, k=8)

# identify if the proportion of cells are different between the 2 groups
#use "scProportionTest" library

prop_test <- sc_utils(seurat_integrated)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WT", sample_2 = "Floxed", 
  sample_identity = "sample")

Test <- prop_test@results$permutation
write.csv(Test, "StatisticsForTheproprtionofCells_res_celltype_from_1.2.csv", row.names = F)
tiff("figures/Celltype_from1_2/SC_proprtionOfcell_scproprtion_res_celltype_from_1.2_log2FD_1.tiff", res =300, height = 4, width=6, units = "in")
permutation_plot(prop_test,log2FD_threshold = 1)
dev.off()


permutation_plot(log2FD_threshold = 1 ,prop_test)
permutation_plot(prop_test)



library(dittoSeq)
#pdf("figures/PercentageOfCells_With_And_WithoutMACSBeads.png", res =300, height = 5, width=4)
p <- print(dittoBarPlot(seurat_integrated, "celltype", group.by = "group",
                        data.out = TRUE))
tiff("figures/Celltype_from1_2/PercentageOfCells_Smad4_celltype_reso_1.2.tiff",res=300, height = 5, width=3.3, units = "in")

p
dev.off()
p




DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "celltype"
#TA, and stem cell with extra Wnt target genes.
features <- c("Lgr5", "Olfm4", "Ascl2", "Mki67", "Tuba1b", "Hist1h1b","Hmgb2", "Cks2", "Axin2","Cd44","Ccnd1", "Ctnnb1", "Lef1", "Bcl9", "Nkd1", "Myc", "Tacstd2", "Fmnl2", "Kcnq1", 
              "SLC12A2", "Cachd1", "Smoc2", "Lrig1", "Ptpro", "Ephb2" )

Clustered_DotPlot(seurat_object = seurat_integrated, features = features, k=8)
tiff("figures/Celltype_from1_2/StemCellmarkers_celltype_from 1.2 resolution.png", res =300, height = 6, width=5)
Clustered_DotPlot(seurat_object = seurat_integrated, features = features, k=8)
dev.off()



DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "celltype"
seurat_integrated <- ScaleData(seurat_integrated)
tiff("figures/Celltype_from1_2/Heatmap_SC_labelled_celltype_top5_only.tiff", res=300,width = 20, height = 25, units = "in")

DoHeatmap(seurat_integrated , features = top10_1.2$gene, slot = "scale.data")

dev.off()


#The analysis performed in the script Difference between WT and Floxedd _Bass Request

SmadWT <- subset(seurat_integrated,  subset = sample == "WT")
SmadFloxed <- subset(x = seurat_integrated, subset = sample == "Floxed")

Idents(object = SmadWT) <- "integrated_snn_res.2.2"
DefaultAssay(SmadWT) <- "RNA"



split seurat_integrated to StemCell_TA, seurat_Paneth, and seurat_Enterocytes, for plotting 



table(seurat_integrated@meta.data$celltype)
Idents(object = seurat_integrated) <- "celltype"
seurat_StemCell_TA <- subset(seurat_integrated, idents = c("Tacstd2 hi", "Fmnl2 hi",  "TA",  "Secretory Pr"  ,   "Cachd1 hi"  ,     "Slc14a1 hi" ,"Slc12a2 hi",      "Lgr5 ISC","Pak3 hi"  ,      "Clu hi"  ))







Idents(object = seurat_integrated) <- "celltype"
seurat_Enterocytes <- subset(seurat_integrated, idents = c("EC-1", "EC-2", "ECepr", "EClpr" ))
DefaultAssay(seurat_Enterocytes) <- "RNA"

#these markers are enterocyte markes from  single servy in the small intestine
Idents(object = seurat_Enterocytes) <- "celltype"
features <- c("Krt19", "Krt8", "Lgals3", "Lypd8", "Selenbp1", "Lgals4", "Ethe1", "Fth1", "Gsdmc4", "Maob", "Ces2a", "Ces2c", "Ces2e", "Chp2", "Coro2a", "Cyp2c65", "Cyp2c68", "Dgkq", "Dpep1", "Emp1", "Krt20", "Cdkn1a", "Mxd1", "Sprr2a3", "Fam3b","Clec2h", "Fhl1", "Gde1", "Guca2a", "Il18", "Reg3g", "Apoa1", "Reg4", "Dmbt1", "Reg3b", "Mt1","Mt2")


tiff("figures/Celltype_from1_2/Heatmap_Enterocyte_subset_labelled_only.tiff", res=300,width = 5, height = 12, units = "in")

DoHeatmap(seurat_Enterocytes , features = features, slot = "scale.data")

dev.off()

tiff("figures/Celltype_from1_2/ClusteredDotPlot_Heatmap_Enterocyte_subset_labelled_only_celltype_from 1.2_resolution.png", res =300, height = 7, width=5, units = "in")
Clustered_DotPlot(seurat_object = seurat_Enterocytes, features = features, k=8)
dev.off()

Clustered_DotPlot(seurat_object = seurat_Enterocytes, features = features, k=8)



features <- c("Reg3b", "Reg3g")
#sample_colors <- c("dodgerblue", "firebrick1")
sample_colors <- c("#999999","#E69F00")
my_levels =c( "ECepr","EClpr", "EC-1", "EC-2")
Idents(seurat_Enterocytes) <- factor(Idents(seurat_Enterocytes), levels= my_levels)

tiff("figures/Celltype_from1_2/Violin_plot_stack_cluster_Reg3g-bEnterocytes_angle90_split_scCustom.png", res=300,width = 4, height = 2, units = "in")
Stacked_VlnPlot(pt.size = 0.05, seurat_Enterocytes, features = features, x_lab_rotate = TRUE,
                colors_use = sample_colors, split.by = "group",  plot_legend = TRUE)+
  theme(axis.text.y = element_text(size = .2))+ #to change the font of the gene name
  guides(color="none") +#to remove the idenity label 
  theme(legend.title = element_text(size=1))
dev.off()

# cd76 Macrophage marker  cd68 Monocyte marker 

features <- c("Cd74", "Cd68")
#sample_colors <- c("dodgerblue", "firebrick1")
sample_colors <- c("#999999","#E69F00")

#Idents(seurat_Enterocytes) <- factor(Idents(seurat_Enterocytes), levels= my_levels)

#tiff("figures/Celltype_from1_2/Violin_plot_stack_cluster_Reg3g-bEnterocytes_angle90_split_scCustom.png", res=300,width = 4, height = 2, units = "in")
Stacked_VlnPlot( seurat_StemCell_TA, features = features, x_lab_rotate = TRUE,
                 colors_use = sample_colors, split.by = "group",  plot_legend = TRUE)+
  theme(axis.text.y = element_text(size = .2))+ #to change the font of the gene name
  guides(color="none") +#to remove the idenity label 
  theme(legend.title = element_text(size=1))
#dev.off()


Plot paneth cell maker alone 

Idents(object = seurat_integrated) <- "celltype"
seurat_Paneth <- subset(seurat_integrated, idents = c("Paneth" ))
DefaultAssay(seurat_Paneth ) <- "RNA"

#these markers are enterocyte markes from  single servy in the small intestine
Idents(object = seurat_Paneth) <- "celltype"
features <- c("Lyz1", "Muc2", "Muc3", "Defa30")
Stacked_VlnPlot( seurat_Paneth, features = features, x_lab_rotate = TRUE,
                 colors_use = sample_colors, split.by = "group",  plot_legend = TRUE)+
  theme(axis.text.y = element_text(size = .2))+ #to change the font of the gene name
  guides(color="none") +#to remove the idenity label 
  theme(legend.title = element_text(size=1))

plots are required for the paper, save them in tiff


Idents(object = seurat_integrated) <- "celltype"
tiff("figures/Celltype_from1_2/ID1_FeaturePlot_labelled_Cell_type_types_Smad4_labelling.tiff", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Id1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
dev.off()

tiff("figures/Celltype_from1_2/Pak3_FeaturePlot_labelled_Cell_type_types_Smad4_labelling.tiff", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Pak3"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
dev.off()
tiff("figures/Celltype_from1_2/Spp1_FeaturePlot_labelled_Cell_type_types_Smad4_labelling.tiff", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Spp1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
dev.off()
tiff("figures/Celltype_from1_2/Lgr5_FeaturePlot_labelled_Cell_type_types_Smad4_labelling.tiff", res =300, height = 4, width=6, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Lgr5"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")

dev.off()


FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Tacstd2"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")
Idents(object = seurat_integrated) <- "celltype"
tiff("figures/Celltype_from1_2/ID1_FeaturePlot_labelled_Cell_type_types_Smad4_labelling_b.tiff", res =300, height = 3, width=5, units = "in")
FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Id1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group", pt.size = 0.2)
dev.off()

FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Id1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group", pt.size = 0.2)



FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Reg3b"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")



FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Reg3g"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")


FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Apoa1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")



FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Sorbs2"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")



FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Ly6a"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")


{#r}
  seurat_integrated <- BuildClusterTree(seurat_integrated, verbose = FALSE, reorder = TRUE, assay="integrated")
  PlotClusterTree(seurat_integrated) 
  Tool(object = seurat_integrated)
  Tool(object = seurat_integrated, slot = "BuildClusterTree")
  
  
  
  FeaturePlot_scCustom(seurat_integrated, 
                       reduction = "umap", 
                       features = c( "Tacstd2"), 
                       order = TRUE,
                       min.cutoff = 'q10', 
                       label = F, split.by = "group")
  FeaturePlot_scCustom(seurat_integrated, 
                       reduction = "umap", 
                       features = c( "Clu"), 
                       order = TRUE,
                       min.cutoff = 'q10', 
                       label = F, split.by = "group")
  
  
  
  
  
  
  sessionInfo()
  
  
  
  