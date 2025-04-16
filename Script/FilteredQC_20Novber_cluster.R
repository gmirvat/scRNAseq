suppressMessages(require(tidyverse))
suppressMessages(require(Seurat))
suppressMessages(require(SingleCellExperiment))

setwd("/mnt/beegfs/hassan/mirvat/Smad4-SCRNAseq/SeuratAnalysis")
options(Seurat.object.assay.version = "v3")

Floxed <- Read10X(data.dir="/mnt/beegfs/hassan/mirvat/Smad4-SCRNAseq/40-887618867/01_analysis/cellranger_count/Floxed/raw_feature_bc_matrix/",  strip.suffix = TRUE)
Floxed  <- CreateSeuratObject(counts = Floxed,
                              min.features = 200, #cutoff from Seurat tutorial
                              min.cells = 3,
                              project = "Floxed")
Floxed
#An object of class Seurat 
#23263 features across 71626 samples within 1 assay 
#Active assay: RNA (23263 features, 0 variable features)
# 2 layers present: counts, data
#class(Floxed[["RNA"]])

WT <- Read10X(data.dir = "/mnt/beegfs/hassan/mirvat/Smad4-SCRNAseq/40-887618867/01_analysis/cellranger_count/WT/raw_feature_bc_matrix",  strip.suffix = TRUE)
WT <- CreateSeuratObject(counts = WT,
                         min.features = 200, #cutoff from Seurat tutorial
                         min.cells = 3,
                         project = "WT")
WT 


#merged_seurat <- WT
merged_seurat <- merge(x = WT, 
                       y = c(Floxed))
add.cell.id = c("WT", "Floxed")

merged_seurat
#An object of class Seurat 
#24506 features across 117814 samples within 1 assay 
#Active assay: RNA (24506 features, 0 variable features)
# 2 layers present: counts.WT, counts.Floxed

merged_seurat@meta.data$ident <- NULL # remove this column to avoid errors after that 


# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
table(merged_seurat@meta.data$orig.ident)



# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-") #The pattern provided ("^MT-") works for human gene names. 
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100 


# Create metadata dataframe
metadata <- merged_seurat@meta.data


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="results/merged_filtered_seurat_smad4_20Nov.RData")
metadata$sample <- metadata$orig.ident
# Create sample column if called differently
#metadata$sample <- NA
#metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
#metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
                           
# Create .RData object to load at any time
save(merged_seurat, file="results/merged_filtered_seurat_smad4_20Nov.RData")
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


length(unique(metadata$cells))
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="results/merged_filtered_seurat_smad4_20Nov.RData")


# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (nGene < 7000) &
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.15))


filtered_seurat
# create a v3 assay
#assay.v3 <- CreateAssayObject(counts = pbmc.counts)
# Extract counts
#counts <- GetAssayData(object = filtered_seurat, slot = filtered_seurat[["RNA"]]$counts)
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0



# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]



# Reassign to filtered Seurat object
#Finally, take those filtered counts and create a new Seurat object for downstream analysis.

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

filtered_seurat

save(filtered_seurat, file="results/seurat_filtered_fom_Raw_0.15Mitocondria_nGeneLimit_20N.RData")


sessionInfo()
