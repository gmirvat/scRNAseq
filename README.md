# scRNAseq
This project analysing scRNAseq from mouse colorectal tumours


#software used for scRNAseq 
## R version 4.2.3 (2023-03-15)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] readxl_1.4.3       viridis_0.6.4      viridisLite_0.4.2  scCustomize_2.0.1 
##  [5] irlba_2.3.5.1      Matrix_1.6-3       RCurl_1.98-1.13    cowplot_1.1.1     
##  [9] scales_1.2.1       lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
## [13] dplyr_1.1.4        purrr_1.0.2        readr_2.1.4        tidyr_1.3.0       
## [17] tibble_3.2.1       ggplot2_3.4.4      tidyverse_2.0.0    Seurat_5.0.1      
## [21] SeuratObject_5.0.1 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] circlize_0.4.15        spam_2.10-0            plyr_1.8.9            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.3         
##   [7] RcppHNSW_0.5.0         listenv_0.9.0          scattermore_1.2       
##  [10] digest_0.6.33          htmltools_0.5.7        fansi_1.0.5           
##  [13] magrittr_2.0.3         paletteer_1.5.0        tensor_1.5            
##  [16] cluster_2.1.4          ROCR_1.0-11            tzdb_0.4.0            
##  [19] globals_0.16.2         matrixStats_1.1.0      timechange_0.2.0      
##  [22] spatstat.sparse_3.0-3  colorspace_2.1-0       ggrepel_0.9.4         
##  [25] xfun_0.41              jsonlite_1.8.7         progressr_0.14.0      
##  [28] spatstat.data_3.0-3    survival_3.5-7         zoo_1.8-12            
##  [31] glue_1.6.2             polyclip_1.10-6        gtable_0.3.4          
##  [34] leiden_0.4.3.1         shape_1.4.6            future.apply_1.11.0   
##  [37] abind_1.4-5            spatstat.random_3.2-1  miniUI_0.1.1.1        
##  [40] Rcpp_1.0.11            xtable_1.8-4           reticulate_1.34.0     
##  [43] dotCall64_1.1-0        htmlwidgets_1.6.2      httr_1.4.7            
##  [46] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
##  [49] farver_2.1.1           pkgconfig_2.0.3        sass_0.4.7            
##  [52] uwot_0.1.16            deldir_1.0-9           janitor_2.2.0         
##  [55] utf8_1.2.4             labeling_0.4.3         tidyselect_1.2.0      
##  [58] rlang_1.1.2            reshape2_1.4.4         later_1.3.1           
##  [61] cellranger_1.1.0       munsell_0.5.0          tools_4.2.3           
##  [64] cachem_1.0.8           ggprism_1.0.4          cli_3.6.1             
##  [67] generics_0.1.3         ggridges_0.5.4         evaluate_0.23         
##  [70] fastmap_1.1.1          yaml_2.3.7             goftest_1.2-3         
##  [73] rematch2_2.1.2         knitr_1.45             fitdistrplus_1.1-11   
##  [76] RANN_2.6.1             pbapply_1.7-2          future_1.33.0         
##  [79] nlme_3.1-163           mime_0.12              formatR_1.14          
##  [82] ggrastr_1.0.2          compiler_4.2.3         rstudioapi_0.15.0     
##  [85] beeswarm_0.4.0         plotly_4.10.3          png_0.1-8             
##  [88] spatstat.utils_3.0-4   bslib_0.5.1            stringi_1.8.1         
##  [91] highr_0.10             RSpectra_0.16-1        lattice_0.22-5        
##  [94] vctrs_0.6.4            pillar_1.9.0           lifecycle_1.0.4       
##  [97] GlobalOptions_0.1.2    spatstat.geom_3.2-7    lmtest_0.9-40         
## [100] jquerylib_0.1.4        RcppAnnoy_0.0.21       data.table_1.14.8     
## [103] bitops_1.0-7           httpuv_1.6.12          patchwork_1.1.3       
## [106] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22    
## [109] gridExtra_2.3          vipor_0.4.5            parallelly_1.36.0     
## [112] codetools_0.2-19       fastDummies_1.7.3      MASS_7.3-60           
## [115] withr_2.5.2            sctransform_0.4.1      parallel_4.2.3        
## [118] hms_1.1.3              grid_4.2.3             snakecase_0.11.1      
## [121] rmarkdown_2.25         Rtsne_0.16             spatstat.explore_3.2-5
## [124] shiny_1.8.0            ggbeeswarm_0.7.2