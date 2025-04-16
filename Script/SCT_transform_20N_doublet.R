#title: "SCT_integration_analysis"
#author: "Mirvat"
#date: "21/11/2023"
#Submitted batch job 5311005
###In this script, I did regress mitochondrial genes. I only regress the diffece in cell cycle genes. 
#I split the seurat and then peformed the analysis. 
#the Normalization and the scaling were OK in the time required for command excution, however, the regress 
#I am using the raw Counts from the 10X genomics
#repeat the analysis as in V3  exactly
#####this is running locally #######
#module load R/default
setwd("/Users/mirvat/Documents/Bioinformatics/March13_2023/Smad4_Single_Cell/Analysis/")
load(file = "results/seurat_filtered_fom_Raw_0.15Mitocondria_nGeneLimit_doubletRemoved_20N.RData")



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
seurat_phase  <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), ndims.print = 6:10, #PCs to print genes for
                        nfeatures.print = 10)
#######
#PC_ 6 
#Positive:  Col3a1, Dcn, Col1a2, Col1a1, Sparc, Igfbp7, Bgn, Col4a2, Col5a2, Col6a1 
##Negative:  Nav3, Itpr2, Sema5a, Samd12, Fhod3, Osbpl3, Lgr5, Dkk2, Clu, Chst11 
#PC_ 7 
#Positive:  Top2a, Mki67, Cdk1, Pclaf, Ube2c, Tpx2, Nusap1, Birc5, Cenpf, Kif11 
#Negative:  Tmsb4x, Gstm5, Clu, S100a14, Lcn2, Fth1, Pglyrp1, Gstm1, Psca, Ifitm1 
#PC_ 8 
#Positive:  Anxa3, Anxa1, Ptpn13, Thbs1, Samd4, F3, Myadm, Fam107b, Lamc2, Glis3 
#Negative:  Spock2, Mt3, Pnliprp2, Gulo, Ptn, Tesc, Pnliprp1, Lect2, Ifitm1, Gadd45g 
#PC_ 9 
#Positive:  Ptn, Atp10a, Cadps, Mt3, Dsc3, Inpp4b, Col23a1, 2010109I03Rik, Spns2, Cwh43 
#Negative:  Wfdc18, Pigr, Lrg1, Itpr2, Ltf, Sectm1b, Sox5, Selenop, Thbs1, Nr1h4 
#PC_ 10 
#Positive:  1810065E05Rik, Gsdmc2, Anxa1, Slc26a2, Cyp2c55, Ethe1, Car1, Abca12, Slc8a1, Gsdmc4 
#Negative:  Isg15, Slfn5, Gbp2, Ifit1, Iigp1, Ifi47, Oasl2, Ifi44, Igtp, Ifit2 
######
# Plot the PCA colored by cell cycle phase

#seurat_phase <- RunPCA(seurat_phase) #for all variation

seurat_phase <- RunPCA(seurat_phase, npcs = 30, features = c(s_genes, g2m_genes))

# Plot the PCA colored by cell cycle phase


png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_before.png", res = 300, height = 4, width =8, units ="in", type="cairo")       
#type="cairo" when running on the cluster
par(mfrow = c(1, 2))
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_before_unsplit.png", res = 300, height = 4, width =8, units ="in", type="cairo") 
DimPlot(seurat_phase)
dev.off()


# view cell cycle scores and phase assignments
head(seurat_phase[[]])

seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score

seurat_phase2 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"))
# cell cycle effects strongly mitigated in PCA
#seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase)) #This command took a long time to run

seurat_phase3 <- ScaleData(seurat_phase, vars.to.regress = c("CC.Difference"))

png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after_regress_cell_cycle_genes.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase2)
dev.off()

png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after_regress_cell_cycle_genes_seperated.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase2,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()
#DimPlot(seurat_phase2)
png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after_Pca.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase2,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()
png("figures/Analysis_V1/PCA_Cell_Cycle_differenceCC_after_Pca.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase3,
        reduction = "pca",
        group.by= "seq_folder")
dev.off()
png("figures/Analysis_V1/PCA_Cell_Cycle_differenceCC.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase3)
dev.off()

png("figures/Analysis_V1/PCA_Cell_Cycle_Genes_after_regress_cell_cycle_diffference_seperated.png", res = 300, height = 3, width =8, units ="in", type="cairo")  
DimPlot(seurat_phase3,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()
#close device that might be open
graphics.off()
png("figures/Analysis_V1/BactEffect.png", res = 300, height = 3, width =8, units ="in", type="cairo")
par(mfrow =c(1, 2))
a <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "seq_folder", 
        split.by= "seq_folder")
b <- DimPlot(seurat_phase, 
        reduction = "pca",
        group.by= "seq_folder")
plot_grid(a, b)
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

png("figures/Analysis_V1/BactEffect_before_1.png", res = 300, height = 5, width =12, units ="in", type="cairo")
DimPlot(seurat_batch, reduction = "umap", group.by = "seq_folder", split.by = "seq_folder" )
dev.off()

png("figures/Analysis_V1/BactEffect_before_UMAP_2.png", res = 300, height = 5, width =4, units ="in", type="cairo")
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
#this command is required to calculated the CC difference
#seurat_phase$CC.Difference <- seurat_phase$S.Score - seurat_phase$G2M.Score
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

seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT") #this step was quick

#This would often be a good place to save the R object.                                  
# Save integrated seurat object


saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_Smad4.rds")

#print(x = seurat_integrated[["pca"]], 
#      dims = 1:10, 
 #     nfeatures = 5)


######The commands below is after running the PCA in the UMAO_Rscript

seurat_integrated <- readRDS(file ="results/Analysis_V1/integrated_seurat_count.rds") # this object doesn't have the PCA

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
#PC_ 1 
#Positive:  Cdh13, Mdga2, Malat1, Exoc4, Ptprt, Rbms3, Znrf3, Tenm4, Dkk2, Nav3 
#Cers6, Naaladl2, Tcf4, Ptn, Dach1, Sema3e, Ube2e2, Spock2, Ank3, Limch1 
#Nav2, Pcdh7, Tspan5, Thsd7a, Inpp4b, Gpcpd1, Prox1, Atp10a, Ifitm3, Prex2 
#Negative:  Lgals4, Lypd8, Phgr1, Plac8, Ceacam1, Muc3, 2200002D01Rik, Lgals3, Actb, Krt20 
#AA467197, Krt19, Emp1, Dmbt1, Muc13, Clca4a, Tspan1, Clca4b, Guca2a, Sprr2a3 
#Serpinb1a, Pmp22, Ahnak, Mgat4c, Myo15b, Ly6a, Cdh17, Cdhr5, Car1, Prdx6 
#PC_ 2 
#Positive:  Cdh13, Mdga2, Exoc4, Ptprt, Mecom, Naaladl2, Rbms3, Malat1, Cers6, Nav3 
#Tenm4, Znrf3, Dkk2, Vcl, Fhit, Pard3b, Samd12, Ube2e2, Neat1, Grip1 
#Zbtb20, Sema3e, Thsd7a, Tcf4, Ext1, Ceacam1, Wwox, Magi1, Lypd8, Ank3 
#Negative:  Rps29, Tmsb4x, S100a6, Rps28, Rpl41, Fth1, Rpl39, Rps24, Ifitm3, Rps27a 
#Fau, Rpl32, Rpl38, Rps15a, Rps8, Rpl37a, Rpl13, Rps27, Rpl37, Rps18 
#Rpl35, Rplp0, Rpl36, Rpl27a, Rpl23, Eef1a1, Rps19, Rps20, Mt1, Rpsa 
#PC_ 3 
#Positive:  Srgn, Cd52, Fcer1g, Tyrobp, Lyz2, Ctss, H2-Aa, Laptm5, H2-Ab1, Cd74 
#Rgs1, Pla2g7, H2-Eb1, Zeb2, Tgfbi, C1qa, C1qb, C1qc, Ptprc, Cybb 
#Fth1, Lgals1, Ctsc, Cxcl2, Clec4n, Fxyd5, Msn, Apoe, Lcp1, Bcl2a1b 
#Negative:  S100a6, Ifitm3, Krt8, Krt18, Rpl41, Clu, S100a11, Rps28, Ifitm1, Pglyrp1 
#Gstm5, Hsp90ab1, Lcn2, Rpl7, Mmp7, Wfdc18, Eef1a1, Tmsb10, Rpl36, Krt19 
#Gpx2, Rplp0, Mt3, Rps14, Rpl37, Stmn1, Rps20, Rpl15, Ifitm2, Ly6e 
#PC_ 4 
#Positive:  St6galnac3, Arhgap15, Ikzf3, Cd3g, Itk, Ptpn22, Cd3e, Srgn, Ctla4, Dock10 
#Emb, Dock2, Rgs1, Cd3d, Cd247, Cxcr6, Trbc2, Cd69, Camk4, Sytl3 
#Grap2, Cytip, Ptprc, Skap1, Fli1, AU020206, Hmgb2, Icos, Trdc, Il2rb 
#Negative:  S100a6, Clu, Gstm5, Lcn2, Cd74, S100a11, Pglyrp1, Actb, Ifitm3, Thbs1 
#Fth1, Nupr1, Ftl1, Dkk2, Fn1, Cldn4, Wfdc18, Anxa1, S100a14, Prdx5 
#Ltf, Lyz2, Krt8, Krt18, Tyrobp, Fcer1g, Ptn, Tgfbi, Tmsb4x, Ly6e 
#PC_ 5 
#Positive:  Fth1, Guca2a, Muc13, 1810065E05Rik, Lyz2, Ctss, H2-Aa, Mt3, Abcb1a, Cyp2c55 
#Rps29, Fcer1g, Tyrobp, Car1, C1qa, C1qb, H2-Ab1, C1qc, Mt1, H2-Eb1 
#Dmbt1, Pla2g7, Lgals4, Clec2h, Slc8a1, Tcf7l2, Pigr, Rps8, Vim, Dach1 
#Negative:  Clu, S100a6, Tmsb4x, St6galnac3, S100a11, Cd3g, Ikzf3, Thbs1, Srgn, Arhgap15 
#Pglyrp1, Itk, Cd3e, Anxa1, Nupr1, Ly6a, Gstm5, Lcn2, Dkk2, Ptpn22 
#S100a10, S100a14, Tmsb10, Dock2, Ltf, Dock10, Ctla4, Ahnak, Cldn4, Cxcr6 
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





png("figures/Analysis_V1/PCA_After_Integration.png", res=300,width = 9, height = 4, units = "in", type="cairo" )

PCAPlot(seurat_integrated) 
dev.off()


png("figures/Analysis_V1/UMAP_AfterIntegration.png", res=300,width = 5, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated,reduction = "umap" ) 
dev.off()

png("figures/Analysis_V1/UMAP_AfterIntegration_split.png", res=300,width = 9, height = 4, units = "in", type="cairo")  

DimPlot(seurat_integrated, split.by = "sample", group.by = "sample",reduction = "umap" ) 
dev.off()
######TSNE if required########
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

png("figures/Analysis_V1/tsne_AfterIntegration.png", res=300,width = 5, height = 4, units = "in", type="cairo")  
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
#
#PC_ 1 
#Positive:  Cdh13, Mdga2, Malat1, Exoc4, Ptprt 
#Negative:  Lgals4, Lypd8, Phgr1, Plac8, Ceacam1 
#PC_ 2 
#Positive:  Cdh13, Mdga2, Exoc4, Ptprt, Mecom 
#Negative:  Rps29, Tmsb4x, S100a6, Rps28, Rpl41 
#PC_ 3 
#Positive:  Srgn, Cd52, Fcer1g, Tyrobp, Lyz2 
#Negative:  S100a6, Ifitm3, Krt8, Krt18, Rpl41 
#PC_ 4 
#Positive:  St6galnac3, Arhgap15, Ikzf3, Cd3g, Itk 
#Negative:  S100a6, Clu, Gstm5, Lcn2, Cd74 
#PC_ 5 
#Positive:  Fth1, Guca2a, Muc13, 1810065E05Rik, Lyz2 
#Negative:  Clu, S100a6, Tmsb4x, St6galnac3, S100a11 
#PC_ 6 
#Positive:  Mt3, Ptn, Pnliprp2, Spock2, Pnliprp1 
#Negative:  Gpx2, Thbs1, Slc12a2, Anxa1, Krt19 
#PC_ 7 
#Positive:  Gm42418, Sparc, Col3a1, Col4a1, Col1a1 
#Negative:  Mdga2, Cdh13, Tmsb4x, Rps29, Lgals4 
#PC_ 8 
#Positive:  Ly6a, Col1a2, Col1a1, Col3a1, Dcn 
#Negative:  Gstm5, Lcn2, Clu, Wfdc18, Cyp2c55 
#PC_ 9 
#Positive:  Gm42418, Camk1d, Cmss1, Hexb, 2010109I03Rik 
#Negative:  Dkk2, Rps29, Ptprt, Mdga2, 1810065E05Rik 
#PC_ 10 
#Positive:  Reg3b, Reg3g, Krt19, Lgals4, Guca2a 
#Negative:  Mxd1, Clca4a, Ifrd1, Gprc5a, March3 
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

png("figures/Analysis_V1/Seurat_integrated_resolution.4.png",res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_Smad4_Resolution.rds")

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

Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_1.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()

Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_1.2_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
dev.off()


Idents(object = seurat_integrated) <- "integrated_snn_res.2.2"

# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_2.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_2.2_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
dev.off()


Idents(object = seurat_integrated) <- "integrated_snn_res.2"

# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_2.0.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()


# Plot the UMAP
png("figures/Analysis_V1/Seurat_integrated_resolution_2.0_split.png", res=300,width = 9, height = 5, units = "in", type="cairo")
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        split.by = "sample")
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
######1.2 is the most stable 

Idents(object = seurat_integrated) <- "integrated_snn_res1.2"

png("figures/Analysis_V1/Seurat_integrated_resolution_1.2.png", res=300,width = 9, height = 5, units = "in", type="cairo")
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4)
dev.off()
#Single-cell RNA-seq clustering analysis lesson 8 



png("figures/Analysis_V1/UMAP_of_cells_in_each_cluster_by_sample_afterIntergration_resolution1.2.png",res=300,width = 12, height = 5, units = "in",type="cairo")
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample", reduction = "umap")  + NoLegend()
dev.off()


png("figures/Analysis_V1/ClusterTree_CD3D_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd3d", node_colour_aggr = "median")

dev.off()

png("figures/Analysis_V1/ClusterTree_CD14_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

clustree(seurat_integrated, prefix = "integrated_snn_res.",
         node_colour = "Cd14", node_colour_aggr = "median")

dev.off()
################
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

png("figures/Analysis_V1/ClusterTree_lgr5_Recluster.png",res=300,width = 20, height = 20, units = "in", type="cairo") # type="cairo" for cluster

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

#saveRDS(seurat_integrated, "results/Analysis_V1/integrated_seurat_Smad4_Resolution.rds")



