graphics.off()

rm_list <- ls()
rm_list <- rm_list[!(rm_list %in% c("expression_data","tpm_loaded"))]  #in case you want to keep these loaded for multiple runs: note having these in memory breaks the save function, because ggplot2 objects save their environments, including these huge objects
rm(list = rm_list)  #removping  all variables
#set.seed(19)
rm(list = ls())
############## Set to directory for home user
setwd("C:\\Users/JAM526/pancan_fpkmuq/met_analysis/pipeline/")
source("control_panel.R")

#resetting things from control panel in this script
base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis/pipeline"   #directory where everything is based
data_dir <- "data/merged"  #an inner-data directory, with the merged data files
data_dir_1 <- "data"      #the outder directory, with other data contained inside. 
output_dir <- "output"   #where data are stored
pathway_lists_dir <- "pathway_lists"   #directory with lists of pathways.  look at what is contained here to add new gene lists
gene_mapping_file <- "ensembl_hgnc_mapping_20200323.rds"   #file with mapping from ensemble ids to hgnc ids
options_list <- as.list(globalenv())      #saving list of options


source("helper_functions.R")  
source("main_functions.R")

  dc <- list() #initializing object that will be piped from function to function
  dc$options <- options_list  #creates a record of all the options used for this run
  
  #the actual pipeline
  dc <- dc %>% load_metadata() %>% initial_mutates() %>% filter_metadata() %>%  load_filter_rename_gen_data() %>% 
    merge_gen_meta() %>% normalize_data() %>% dim_reduce() %>% cluster() %>% add_relation_points()  %>%
    plot() %>% chi_square()  %>%  run_rfc() %>% volcano() %>% make_heatmap() %>% save_output()
  rm(list = "expression_data")
  dc <- dc[names(dc)[!(names(dc) %in% c("expression_data"))]]  #deleding expression data from DC because it is too large to save multiple times. 
  
  

  graphics.off()
  
  


#list of commands to outpu results
 dc$scatter_color2_lines
 dc$cluster_comp_bp_pvals
 dc$rfc_plot
 dc$volcano
 dc$heatmap_rowmean
 dc$rfc_plot

