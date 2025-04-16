cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)


# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

View(cluster0_ann_markers)


# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


conserved_markers <- map_dfr(c(0:5), get_conserved) ###Based on the current 
conserved_markers0_5 <- conserved_markers

conserved_markers7_20 <- map_dfr(c(7:21), get_conserved) 

#clustertree map, I will use resolution 2.2  
#if use group that has 1 and 2 categories 
#names(conserved_markers)[names(conserved_markers)=="1_avg_logFC"] <-"G1_avg_logFC"
#names(conserved_markers)[names(conserved_markers)=="2_avg_logFC"] <-"G2_avg_logFC"

top10_0to5 <- conserved_markers0_5 %>% 
  dplyr::mutate(avg_fc = (WT_avg_log2FC + Floxed_avg_log2FC) /2) %>% 
  dplyr::group_by(cluster_id) %>% 
  top_n(n = 10,  wt = avg_fc)



write.csv(conserved_markers, "results/Analysis_V1/Conserved_maker_identification_per_cluster_resolution_2.2.csv", row.names = F)

#saveRDS(seurat_integrated, "results/NormalControl_G1_G2/integrated_seurat_count_G1_G2.rds")

conserved_markers <- map_dfr(c(6), get_conserved)
#Warning message:
#In FindMarkers.default(object = data.use, slot = data.slot, counts = counts,  :
#                         No features pass logfc.threshold threshold; returning empty data.frame
# remove only.pos = TRUE from  get_conservd
annotations <- read.csv("data/annotations.csv")
get_conserved2 <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample")%>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
Idents(object = seurat_integrated) <- "integrated_snn_res.2.2"
conserved_markers <- map_dfr(c(6), get_conserved2)

conserved_markers
write.csv(conserved_markers, "results/Analysis_V1/Conserved_marker_cluster6_comparingtoAll.csv")
#conserved_markers
#cluster_id    gene     WT_p_val WT_avg_log2FC WT_pct.1 WT_pct.2 WT_p_val_adj
#1          6   Mdga2 2.840176e-07    -0.3997712    0.316    0.452 5.682625e-03
#2          6    Tspo 3.358119e-08    -0.2940823    0.381    0.534 6.718924e-04
#3          6     B2m 6.422903e-10    -0.3767957    0.437    0.621 1.285094e-05
#4          6 S100a10 7.965801e-07    -0.2701116    0.526    0.658 1.593798e-02
#5          6   Plac8 1.621453e-08    -0.4554718    0.121    0.293 3.244203e-04
#6          6  Lgals4 7.295212e-05    -0.3856197    0.237    0.350 1.000000e+00
#Floxed_p_val Floxed_avg_log2FC Floxed_pct.1 Floxed_pct.2 Floxed_p_val_adj     max_pval
#1 5.489966e-28        -0.3729133        0.434        0.617     1.098432e-23 2.840176e-07
#2 9.908842e-21        -0.2533129        0.093        0.239     1.982561e-16 3.358119e-08
#3 1.387757e-17        -0.2640024        0.227        0.392     2.776623e-13 6.422903e-10
#4 1.051744e-16        -0.2578369        0.223        0.378     2.104330e-12 7.965801e-07
#5 8.379076e-16        -0.2921500        0.061        0.169     1.676486e-11 1.621453e-08
#6 4.626946e-11        -0.2897006        0.197        0.302     9.257594e-07 7.295212e-05
#minimump_p_val
#1   1.097993e-27
#2   1.981768e-20
#3   2.775513e-17
#4   2.103489e-16
#5   1.675815e-15
#6   9.253893e-11
#description
#1 MAM domain containing glycosylphosphatidylinositol anchor 2 [Source:MGI Symbol;Acc:MGI:2444706]
#2                                          translocator protein [Source:MGI Symbol;Acc:MGI:88222]
#3                                          beta-2 microglobulin [Source:MGI Symbol;Acc:MGI:88127]
#4                S100 calcium binding protein A10 (calpactin) [Source:MGI Symbol;Acc:MGI:1339468]
#5                                         placenta-specific 8 [Source:MGI Symbol;Acc:MGI:2445289]
#6                         lectin, galactose binding, soluble 4 [Source:MGI Symbol;Acc:MGI:107536]
