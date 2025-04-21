
#Smad4_single Cell RNA-seq Creating Seurat object and QC with visualization

#References: 
# 1. https://satijalab.org/seurat/articles/seurat5_essential_commands

# 2. Germain P, Lun A, Garcia Meixide C, Macnair W, Robinson M (2022). “Doublet identification in single-cell sequencing data using scDblFinder.” f1000research. doi:10.12688/f1000research.73600.2.
####the most important part is 
options(Seurat.object.assay.version = "v3") #before start reading the files


setwd ("~/Documents/scRNA-seq/" )# set it to your preferred directory

#-----------------------------------------
#upload the libraries 
# -----------------------------------------
library(tidyverse)
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
library(SingleCellExperiment)
library(Seurat)
library(scDblFinder)

#-----------------------------------------
# create a v3 assay
#-----------------------------------------
options(Seurat.object.assay.version = "v3")

#-----------------------------------------
#Read the seq folders out put and creat Seurat object
#-----------------------------------------
Floxed <- Read10X(data.dir="Analysis/rawCounts/Floxed/raw_feature_bc_matrix/",  strip.suffix = TRUE)
Floxed  <- CreateSeuratObject(counts = Floxed,
                              min.features = 200, #cutoff from Seurat tutorial
                              min.cells = 3,
                              project = "Floxed")
#to view the object type 
Floxed
# read the second folder 
WT <- Read10X(data.dir = "rawCounts/WT/raw_feature_bc_matrix/",  strip.suffix = TRUE)
WT <- CreateSeuratObject(counts = WT,
                         min.features = 200, #cutoff from Seurat tutorial
                         min.cells = 3,
                        project = "WT")
##to view the object type 
WT 
#-----------------------------------------
#to identify the doublets 
#create SingleCellExperiment from the created CreateSeuratObject
#Find the doublets usisng scDblfinder()
#check the results
#return back the SingleCellExperiment to Seurat object
#-----------------------------------------

Floxed.sce <- as.SingleCellExperiment(Floxed)
Floxed.sce <- scDblFinder(Floxed.sce, verbose=FALSE)
# how many singlet and doublet 
table(Floxed.sce$scDblFinder.class)
#Return back to Seurat
Floxed  <- as.Seurat(Floxed.sce)


WT.sce <- as.SingleCellExperiment(WT)
WT.sce <- scDblFinder(WT.sce, verbose=FALSE)
table(WT.sce$scDblFinder.class)#
#-----------------------------------------
#Seurat automatically creates some metadata for each of the cells when you use the Read10X() function to read in data.
#This information is stored in the meta.data slot within the Seurat object.
#metadata coulmns 
#orig.ident: this often contains the sample identity if known, but will default to "SeuratProject"
#nCount_RNA: number of UMIs per cell
#nFeature_RNA: number of genes detected per cell
### Creat one merged object 
# Explore the metadata
#-----------------------------------------
head(Floxed@meta.data)
head(WT@meta.data)

#-----------------------------------------
#merged_seurat
#-----------------------------------------
merged_seurat <- merge(x = WT, 
                           y = c(Floxed))

merged_seurat

merged_seurat@meta.data$ident <- NULL # remove this column to avoid errors after that 
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

#Generating quality metrics
View(merged_seurat@meta.data) #created automatically in Seurat
  

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
#-----------------------------------------
# Compute  Mitochondrial Ratio usisng The PercentageFeatureSet() function 
#The pattern provided ("^MT-") works for human gene names, "^mt-" for mice
# -----------------------------------------
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-") 
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
    
# -----------------------------------------
# Create metadata dataframe
# -----------------------------------------
metadata <- merged_seurat@meta.data
# -----------------------------------------
# Add cell IDs to metadata
# -----------------------------------------
metadata$cells <- rownames(metadata)  
# -----------------------------------------
# Create sample column based on original identity
# -----------------------------------------
metadata$sample <- metadata$orig.ident
# -----------------------------------------
# Rename columns for consistency
# -----------------------------------------
metadata <- metadata %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

# -----------------------------------------
# Check number of unique cells
# -----------------------------------------
length(unique(metadata$cells))  

# -----------------------------------------
# Save updated metadata to Seurat object
# -----------------------------------------
merged_seurat@meta.data <- metadata

# -----------------------------------------
# Save the Seurat object for future use
# -----------------------------------------
save(merged_seurat, file = "results/merged_filtered_seurat_withdoublets.RData")
# -----------------------------------------
#1. Assessing the quality metrics
#2. Cell counts
#3. UMI counts per cell
#4. Genes detected per cell
#5. Complexity (novelty score)
#6. Mitochondrial counts ratio
#7. Doublets 
# -----------------------------------------

#--------- Visualize QC before cleaning and applying the cutoff # ----------

# Visualize the number of cell counts per sample
metadata <- merged_seurat@meta.data 
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

plot(metadata$nUMI, metadata$mitoRatio, pch=20, cex=0.5)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

#----------Based on the quality matrix #------------
# Filter out low quality cells using selected thresholds - these will change with experiment
#---------------------------------
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) &
                            (nGene < 7000) &
                            (scDblFinder.class != "doublet")&
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.15))
filtered_seurat

#---------------------------------
#Gene level filtering 
#many genes with zero counts. 
#These genes can dramatically reduce the average expression for a cell, remove them from the data
#---------------------------------
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data
#inspect
table(metadata_clean$sample)

#---------------------------------
#TO visualize the data following the filer 
#1. use the metadata_clean
## Repeat the above plots
#---------------------------------

#---------------------------------
####Save results 
#---------------------------------
save(filtered_seurat, file="results/seurat_filtered_fom_Raw_0.15Mitocondria_nGeneLimit_doubletRemoved_20N.RData")
#---------------------------------
