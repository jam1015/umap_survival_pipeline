#more packages than necessary
library(lubridate)
library(Matrix.utils)
library(class)
library(pracma)
library(dplyr)
library(stringr)
library(umap)
library(scales)
library(ggplot2)
library(ggthemes)
library(dbscan)
library(MASS)
library(RVAideMemoire)
library(made4)
library(plotly)
library(Rtsne)
library(TCGAbiolinks)
library(randomForest)
#library(org.Hs.eg.db)
library(GenomicFeatures)
library(dplyr)
#library(biomaRt)
library(randomcoloR)
library(magrittr)
library(rlang)
library(survival)
library(survminer)
library(survMisc)
library(rsample)     
library(randomForest)
library(ranger)      
library(caret)       
library(h2o)        
library(grid)
library(ggplotify)
library(gridExtra)
library(tsne)
library(Rtsne)
library(irlba) 
library(SummarizedExperiment)
library(tools)
library(ComplexHeatmap)
library(viridis)
library(EnhancedVolcano)
library(tidyverse)
library(ggrepel)


#options to be changed by shiny user
load_local_gen <- TRUE  #should we load gen data from a local array or should we to tcgabiolinks

three_dee_output <- TRUE #should we do 3d output
plot_black_rings <- FALSE #sholud we plot black rings, otherwise mets ar holoow points
output_dims <- ifelse(three_dee_output, 3 , 2)
#how should the data be normalized? 
whiten <- FALSE   #whitening removes correlation from the dataset; only for tSNE
standardize <-  TRUE   #should we normalize the data to 1
standardize_cols <- FALSE  #normalizes data to 1 column-wise
spherize <- TRUE     #should we center and project onto unit hypersphere
filter_noise <- FALSE           #should we get rid of noise from umap after clustering (group zero)
survival_zero <- FALSE         #should we include the group zero on the survival plot         

#should we use umap or tsne
use_umap <- TRUE  #should we use umap? if not, use tsne
use_umap_wrapper <- TRUE # should we use umap learn, or pure R implementation? 
use_tsne_wrapper <- TRUE # should we use the tsne c++ wrapper or native r_versin

# what umap options should we use
custom.config <- umap.defaults    #initializing umap configuration
custom.config$n_epochs <- 2000    #how much should we run umap for
custom.config$min_dist <- 0.001#10^(-6)    #minimum distance for umap
custom.config$n_neighbors <- 5     #number of neighbors for umap
custom.config$n_components <- output_dims     # how many dimensions we are doing umap into
hdbscan_min_pts <- 5   #min points in a cluster for hdbscan
force_hard_cluster <- TRUE  #should we use knn to force hdbscan noise points to belong to the nearst cluster
k_for_knn <- 1 #how many neighbors to look at when forcing hard clustering


#other umap parameters that one could change.  see documentation of umap package for details.
custom.config$metric <- "euclidean"
custom.config$input <- "data"
custom.config$init <- "spectral"
custom.config$set_op_mix_ratio <- 1
custom.config$local_connectivity <- 1
custom.config$bandwidth <- 1
custom.config$alpha <- 1
custom.config$gamma <- 1
custom.config$negative_sample_rate <- 5
custom.config$a <- NA
custom.config$b <- NA
custom.config$spread <- 1
custom.config$random_state <- NA
custom.config$transform_state <- NA
custom.config$knn_repeats <- 1
custom.config$verbose <- TRUE
custom.config$umap_learn_args <- NA





############### tsne_options if we are not using umap
#########################################
#basic parameters to change for tsne: note that this is a mixture of the parameters from both the tsne and Rtsne package
tsne_perplexity <- 10
tsne_max_iter <- 500
tsne_eta <- 200

#slightly more advanced parameters to change for tsne
tsne_stop_lying_iter <- 250L
tsne_mom_switch_iter <- 250L
tsne_momentum <- .5
tsne_final_momentum <- .8
tsne_exaggeration_factor <- 12
tsne_theta <- .5    # set to zero for pure t-sne, positive <1 for bhtsne
tsne_min_cost <- 0   #only tsne package

#things having to do with normalization that we shouldn't set because we have a normalization function
tsne_pca <- FALSE 
tsne_pca_center <- TRUE
tsne_pca_scale <- TRUE
tsne_normalize <- TRUE
tsne_initial_dims <- NA #should we reduce the dimension with PCA first? not done here because we want to do this in the normalization function
tsne_partial_pca <- FALSE
tsne_whiten <- FALSE   #only tsne package.  should be done in normalization function. 
tsne_Y_init <- NULL

#more tsne_options
tsne_check_duplicates <- FALSE
tsne_verbose <- TRUE
tsne_is_distance <- FALSE  #can pre-compute a distance matrix
tsne_index <- NA            # for nearest neighbors function (unused)
tsne_distance <- NA         #for  nearest neighbors function (unused)
tsne_epoch <- 100            #only tsne package; how long an epoch is 
tsne_epoch_callback <- NULL  #only tsne package; a function to excecute every epoch
tsne_num_threads <- 1


#########################################

#the simple load loads the expression data into memory each time, and deletes it while 
#the non-simple load keeps it and memory but causes problems for saving ggplot objecst which try to save their environment
use_simple_load <- FALSE

#pathway and cancer to use
used_pathway_file <- c("rp.csv")  #which pathway should we use

#these options will select a random genes from the transcriptome
use_rand_sample <- FALSE
rand_sample_number <- 25
if(use_rand_sample){
human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") #the biomart package to get the genes from 
}

filter_dataset <- TRUE        #are we filtering by what dataset we are using
disease_filt <- "BRCA"            #which disease are we using
dataset_filt <- "Lee"    #which dataset should we use
filter_disease <- FALSE #should we filter the disease
filter_sample_type <- FALSE  #should we filter the sample type?
#sample_type_filt <- "Primary Blood Derived Cancer - Peripheral Blood"
sample_type_filt <- c("Metastatic")   #what sample type should we filter for? 



tcga_mets <- FALSE    #special case to look at data that has primary metastasis site in tcga data
tcga_met_file <- "tcga_metastasis.csv" #additional metastasis file infor from tcga; the primary metastasis site


## where the data are: need to change this
base_dir <-  "C:/Users/JAM526/pancan_fpkmuq/met_analysis/pipeline"   #directory where everything is based
data_dir <- "data/merged"  #an inner-data directory where separate datasets are merged
data_dir_1 <- "data"      #the outder directory
output_dir <- "output"   #where data are stored
pathway_lists_dir <- "pathway_lists"   #where we store the lists of pathways

tcga_projects <- getGDCprojects()  #df with all of the tcga projects
## Useful for knowing the names/filenames of the pathways

pathway_table <- read_csv(file.path(base_dir,data_dir_1,pathway_lists_dir,"gene_list_master.csv")) #the master pathway table



#plotting options
chi_sq_var <- "target_site" 
chi_sq_var_quo <- quo("target_site")#what variable do we run chi square on? 
color_selection_1 <- quo(disease)  #what do we do color selection with
border_selection <- quo(sample_type) #really a shape selection; what shape do we represent things with 
color_selection_2 <- quo(Cluster)  #how do we color the second plot (for the clusters)
color_selection_3 <- quo(biopsy_tissue)
shape_selection <- NULL  #quo(disease)     #use a quosure.  for plotting shape selections.
print_p_threshold <- 0.01   #after 
use_tpm <- FALSE   #should we use tpm or fpkm
 
draw_lines <- TRUE #draw lines between matched samples
   #are we jsut seleecting a few datasets
use_plotly <- TRUE


symbol_size <- 1.5
symbol_size_factor <- 1.5   #what symbol size on the plot? 
num_trees <- 500   #number of trees for random forest
fig_res <- 300     #figure resolution for output
subset_for_stats <- FALSE   #when we run statistics should we use a subset?
stats_filter_string <- "Tumor" #if we are subsettinf for stats what filter do we use? 
save_image <- TRUE            # are we saving an image
axis_numbers <- FALSE         # are we displaying axis numbers
show_p_vals <- TRUE           # are we displaying p-values
#should we display every run? 
show_one_iteration  <- FALSE  #if we are looping throug multiple runs do we show all of them in output?

tpm_file <- "joined_tpm_mat.rds"       #where is the tpm file?
fpkm_file <- "joined_fpkm_mat.rds"     #where is the tpm file for all the data?
fisher <- TRUE   #do we use post-hoc fisher or chi-square residuals for output?
font_size <- 14 #font size on outputs
save_data <- TRUE

#how should we add noise: note much better to just raise the "min dist" param of umap
jitter_cluster <- FALSE  #should we add noise to the data before clustering?
jitter_plot <- FALSE     #should we add the clustering noise to the output plot?
jitter_factor <- NULL   #parameter for jitter function
jitter_amount <- .05    #parameter for jitter function


vars_for_heatmap <- 35 #how many genes to include in heatmap





filter_for_stats <- TRUE     # 
stats_filter_elim <- "Primary Tumor" #what we get rid of when filtering for stats
stats_filter_var <- quo(sample_type)  #what variable we are filtering on when filtering for stats

run_stats <- TRUE              #should we we do the chi square analysis 
survival_analysis <- TRUE      #should se do a surival analysis
contingency_analysis <- FALSE  #should we check cluster composition

rfc_mtry <- seq(5, 40, by = 2) #mtry parameter of rfc


out_fig_width <- 6.5     #width of output figs     
out_fig_height <- 6.5     #height of output figs
out_fig_units <- "in"    # what units for output figs


exclude_title_legend <- c("patient_id","sample","tissue","sample_tissue","biopsy_tissue","symbol_size")  #never show these in a title of legend

#these colors files are used at various points in the pipeline
colors_file_1 <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/pipeline/data/color_palettes/pa200_set1_chroma23_lightness23.RDS"
colors_file_2 <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/pipeline/data/color_palettes/aeb_colors.RDS"
colors_file_3 <- "C:/Users/JAM526/pancan_fpkmuq/met_analysis/pipeline/data/color_palettes/out_pal_33.RDS"


shape_scale  <- c(`Metastatic` = 21, `Primary Tumor` = 19, `Solid Tissue Normal` = 23,`Recurrent Tumor` = 22)   #ggplot2 codes for shapes to use 
shape_scale_3d <- c(`Metastatic` = "circle-open", `Primary Tumor` = 'circle', `Solid Tissue Normal` = "diamond-open", `Recurrent Tumor` = "square-open")      #plotly codes for shapes in 3d output


save_output <- TRUE  #should we save the output

volcano_text_threshold <- 4.5 # what is the value of -log10(p) above which we shuld plot labels


run_transpose <- FALSE #if run_transpose is true we run umap on the genes in question, outputting genes.  not well tested.  keep false for now. 








