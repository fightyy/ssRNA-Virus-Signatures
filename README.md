### ssRNA-Virus-Signatures

This is the code for **An expanding universe of mutational signatures and its rapid evolution in single stranded RNA viruses (ssRNA-Virus-Signatures)**.

#### Directory structure

The code for ssRNA-Virus-Signatures is organized into different directories and scripts.

The directory structure is as follows:

- bin :  This directory contains ggplot2 theme that is needed for regenerating plots.
- data : This directory contains data that is needed for regenerating plots and descriptions for the data .
- Figure1-5 : This directory contains the analysis scripts for regenerating main figures  .
- supplementary_figure :  This directory contains the scripts that are related to reproducing the supplementary figures.

#### Description of Scripts and How to Use

##### Figure1-5 and supplementary_figure

- Description: These folders contains R scripts to generate plots used in the article. All data can be found in the `data` directory. You can reproduce the results by running these R scripts like `Rscript figure1.R` and shell scripts like `sh Figure5ABC.sh` 

#### *Dependencies*

R session info:

```R
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.6.1

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
other attached packages:
 [1] ggpmisc_0.5.5             ggpp_0.5.6                cowplot_1.1.3             ggpubr_0.6.0              rstatix_0.7.2            
 [6] cols4all_0.7-1            ggnewscale_0.4.10         ggtreeExtra_1.12.0        ggtree_3.10.1             dendextend_1.17.1        
[11] ggdendro_0.2.0            mSigTools_1.0.7           MutationalPatterns_3.12.0 NMF_0.27                  Biobase_2.62.0           
[16] cluster_2.1.6             rngtools_1.5.2            registry_0.5-1            GenomicRanges_1.54.1      GenomeInfoDb_1.38.6      
[21] IRanges_2.36.0            S4Vectors_0.40.2          BiocGenerics_0.48.1       scales_1.3.0              vroom_1.6.5              
[26] car_3.1-2                 carData_3.0-5             factoextra_1.0.7          proxy_0.4-27              lsa_0.73.3               
[31] SnowballC_0.7.1           ggsci_3.0.1               this.path_2.4.0           lubridate_1.9.3           forcats_1.0.0            
[36] stringr_1.5.1             dplyr_1.1.4               purrr_1.0.2               readr_2.1.5               tidyr_1.3.1              
[41] tibble_3.2.1              ggplot2_3.5.0             tidyverse_2.0.0          

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.15.0       jsonlite_1.8.8          magrittr_2.0.3          rmarkdown_2.26         
  [6] farver_2.1.1            fs_1.6.3                zlibbioc_1.48.0         ragg_1.2.7              vctrs_0.6.5            
 [11] memoise_2.0.1           RCurl_1.98-1.14         base64enc_0.1-3         htmltools_0.5.7         polynom_1.4-1          
 [16] broom_1.0.5             Formula_1.2-5           gridGraphics_0.5-1      pracma_2.4.4            htmlwidgets_1.6.4      
 [21] plyr_1.8.9              cachem_1.0.8            lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
 [26] Matrix_1.6-5            R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.11 digest_0.6.34          
 [31] aplot_0.2.2             colorspace_2.1-0        patchwork_1.2.0         confintr_1.0.2          Hmisc_5.1-2            
 [36] textshaping_0.3.7       labeling_0.4.3          spacesXYZ_1.3-0         fansi_1.0.6             timechange_0.3.0       
 [41] abind_1.4-5             mgcv_1.9-1              compiler_4.3.1          bit64_4.0.5             withr_3.0.0            
 [46] doParallel_1.0.17       htmlTable_2.4.2         backports_1.4.1         viridis_0.6.5           ggsignif_0.6.4         
 [51] quantreg_5.97           MASS_7.3-60.0.1         tools_4.3.1             foreign_0.8-86          ape_5.7-1              
 [56] nnet_7.3-19             glue_1.7.0              nlme_3.1-164            grid_4.3.1              checkmate_2.3.1        
 [61] gridBase_0.4-7          reshape2_1.4.4          generics_0.1.3          gtable_0.3.4            tzdb_0.4.0             
 [66] data.table_1.15.2       hms_1.1.3               utf8_1.2.4              XVector_0.42.0          ggrepel_0.9.5          
 [71] foreach_1.5.2           pillar_1.9.0            yulab.utils_0.1.4       splines_4.3.1           treeio_1.26.0          
 [76] lattice_0.22-5          survival_3.5-8          bit_4.0.5               SparseM_1.81            tidyselect_1.2.0       
 [81] knitr_1.45              gridExtra_2.3           xfun_0.42               stringi_1.8.3           lazyeval_0.2.2         
 [86] ggfun_0.1.4             evaluate_0.23           codetools_0.2-19        ggplotify_0.1.2         cli_3.6.2              
 [91] rpart_4.1.23            systemfonts_1.0.6       munsell_0.5.0           Rcpp_1.0.12             png_0.1-8              
 [96] MatrixModels_0.5-3      ggalluvial_0.12.5       bitops_1.0-7            viridisLite_0.4.2       tidytree_0.4.6         
[101] crayon_1.5.2            rlang_1.1.3            
```

Python packages：

```python
## Requirements
To run the shell scripts, you will need the following Python packages:
- [sigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting)
- [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment)
```

#### *Data source*

We downloaded ssRNA viral sequences from multiple databases including the Bacterial and Viral Bioinformatics Resource Center (BVBRC)(Olson, et al. 2023), the Global Initiative on Sharing All Influenza Data (GISAID)(Shu and McCauley 2017), Nextstrain(Hadfield, et al. 2018) and the literature(Kustin T, Stern A. 2021. Biased Mutation and Selection in RNA Viruses. Mol Biol Evol 38:575-588.). You can download data from the accession column of `data/all_meta_corrected_host.csv` .

