

#remotes::install_github("rnabioco/clustifyrdata")
library(clustifyrdata)
seurat_integrated <- clustify(
  input = seurat_integrated,
  cluster_col = "integrated_snn_res.2.2",
  ref_mat = ref_MCA,
  seurat_out = TRUE
)
metadata <- seurat_integrated@meta.data
# calculate correlation



res <- clustify(
  input = seurat_integrated,
  metadata = seurat_integrated$type,
  ref_mat = ref_MCA
)

plot_best_call(
  cor_mat = res,
  metadata = pbmc_meta,
  cluster_col = "classified"
)


library(clustifyrdata)
#https://github.com/rnabioco/clustifyrdatahub
seurat_integrated <- clustify(
  input = seurat_integrated,
  cluster_col = "integrated_snn_res.2.2",
  ref_mat = ref_mouse.rnaseq,
  seurat_out = TRUE
)

metadata <- seurat_integrated@meta.data
table(metadata$type.clustify)
library(clustifyr)
library(usethis)
library(SingleR)
ref_mouse.rnaseq <- average_clusters(
  SingleR::mouse.rnaseq$data,
  SingleR::mouse.rnaseq$types,
  if_log = F
)

d1 <- c(
  paste0("[", "ref_MCA", "](", "https://github.com/rnabioco/clustifyrdata/raw/main/data/ref_MCA.rda", ")"),
  "Mouse Cell Atlas",
  dim(clustifyrdata::ref_MCA)[2],
  dim(clustifyrdata::ref_MCA)[1],
  "mouse",
  paste0("[from](", "https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8", ")")
)

usethis::use_data(ref_mouse.rnaseq, compress = "xz", overwrite = TRUE)
