#the main functions
load_metadata <- function(dc_in)  { #dc_in stands for data container
  
  dc_out <- dc_in
  
  
  #loading the metadata about each sample
  all_metadata <-  readRDS(file.path(base_dir,data_dir,"joined_meta.rds")) #loading a single matrix of all metadata, pre-wrangled
  surv_data <- readRDS(file.path(base_dir,data_dir,"gdc_survival.rds"))  #loading a matrix of all the survival data
  surv_data <- surv_data[,c("sample","OS","OS.time")]    #just getting essential survival data
  all_metadata <- left_join(all_metadata,surv_data,by = "sample") %>% as_tibble()     #joining survival data to all metadata
  
  dc_out$tcga_met_data <- read_csv(file.path(base_dir,data_dir,tcga_met_file))    #tcga metastasis data
  dc_out$all_metadata <- all_metadata   
  
  return(dc_out)
}

initial_mutates <- function(dc_in){
  
  
  #filling in missing tissue with biopsy tissue
  all_metadata <- dc_in$all_metadata
  
  na_inds <- is.na( all_metadata$biopsy_tissue)   # biopsy tissue is the tissue that it was actuall taken from
  all_metadata$biopsy_tissue[na_inds] <- all_metadata$tissue[na_inds]   #tissue represents the "origin tissue".  should probably be labeled that
  
  
  
  #making tables for the first descendant and last ancestor of each sample.
  dc_in$descendant_table <- make_descendant_table_meta(all_metadata)    
  dc_in$ancestor_table <-  make_ancestor_table_meta(all_metadata)
  
  all_metadata <- get_target_site(all_metadata)  #appends target site to  the metadata; the final known location.  could get complicated if there are multiple mets
  
  
  
  #these next two lines can be updated to avoind making this temporary dataframe and also to asign target site to everyone with the same patient id
  temp_met <- left_join(all_metadata, dc_in$tcga_met_data, by = "sample")  #appending tcga met data
  all_metadata$target_site[!is.na(temp_met$new_tumor_event_site)] <- temp_met$new_tumor_event_site[!is.na(temp_met$new_tumor_event_site)] #updating target site
  
  dc_in$all_metadata <- all_metadata  #
  
  
  dc_out <- dc_in
  
  
  return(dc_out)
}



##################################################
filter_metadata <- function(dc_in){
  
  #loading ancestor/descendant tables, calling them "filtered" in advance
  ancestor_table <- dc_in$ancestor_table      
  descendant_table<- dc_in$descendant_table
  
  #loading all the other metadata
  metadata <- dc_in$all_metadata %>% unique()
  metadata <- dplyr::filter(dc_in$all_metadata,!is.na("disease"))
  
  if(tcga_mets){
    #filtering to only the data with primary metastasis site in tcga
    metadata <- metadata %>% filter(metadata$sample %in% dc_in$tcga_met_data$sample)
  }
  
  #these next three if statements filter the datasets (metadata, ancestor, descendant) according to values in control panel
  
  if(filter_dataset){
    
    metadata <- dplyr::filter(metadata,dataset %in% dataset_filt)  
    ancestor_table <- dplyr::filter(ancestor_table,dataset %in% dataset_filt) 
    descendant_table <- dplyr::filter(descendant_table,dataset %in% dataset_filt) 
  }
  if(filter_disease){ 
    metadata <- dplyr::filter(metadata,disease %in% disease_filt)
    ancestor_table <- dplyr::filter(ancestor_table,disease %in% disease_filt)
    descendant_table <- dplyr::filter(descendant_table,disease %in% disease_filt)
    
  }
  if(filter_sample_type){ 
    metadata <- dplyr::filter(metadata,sample_type %in% sample_type_filt)
    ancestor_table <- dplyr::filter(ancestor_table,sample_type %in% sample_type_filt)
    descendant_table <- dplyr::filter(descendant_table,sample_type %in% sample_type_filt)
    
  }
  
  metadata_filtered <- metadata
  ancestor_table_filtered <- ancestor_table
  descendant_table_filtered <- descendant_table
  
  dc_out <- dc_in
  
  dc_out$metadata_filtered <- metadata_filtered
  dc_out$descendant_table_filtered <- descendant_table_filtered
  dc_out$ancestor_table_filtered <- ancestor_table_filtered
  
  return(dc_out)
}

##################################################



load_filter_rename_gen_data <- function(dc_in){
  
  
  
  
  
  #deciding what genetic data loading function to use
  if (load_local_gen){
    if(use_simple_load){
      expression_data <- load_gen_data_simple()
    } else{
      expression_data <- load_gen_data()
    }
    
    ensembl_gene_mapping <- readRDS(file.path(base_dir,data_dir_1,gene_mapping_file)) #reads a mapping of gene ids made with biomaRt
    ensembl_gene_mapping <-  ensembl_gene_mapping %>% 
      group_by(hgnc_symbol) %>% mutate(symbol_count = n()) %>% ungroup()  #getting the number of times each hgnc symbol appears (mapped to multiple ensembl ids) 
    ensembl_gene_mapping <-  ensembl_gene_mapping %>% 
      mutate(unique_code = str_c(hgnc_symbol,str_extract(ensembl_gene_id ,"([1-9]{1,}0*)*$"),sep = "_"))  #appending tail of ensembl id  
    ensembl_gene_mapping <-  ensembl_gene_mapping %>%
      mutate(hgnc_unique  = case_when(symbol_count > 1 ~ unique_code, symbol_count <= 1 ~ hgnc_symbol))   #if there are more than one occurences of a gene symbol use the symbol with the appended ensembl code scale 
    ensembl_gene_mapping %<>% dplyr::select(ensembl_gene_id, hgnc_symbol, hgnc_unique) #selecting ensembl id, blank hgnc symbol, and hgnc symbol with appended suffix
    
    
    ensembl_gene_mapping <-  ensembl_gene_mapping %>% 
      dplyr::filter(ensembl_gene_id %in% colnames(expression_data))   #filtering for only the values we have in our data
    
    if (use_rand_sample) { #if we use a random sample we select a random set of genes from the ensembl mapping file
      
      
      mapped_genes <- ensembl_gene_mapping %>% sample_n(size = rand_sample_number)
      tx <- mapped_genes$hgnc_symbol 
      
    } else{
      
      
      
      #reading the file we chose to use that has the list of the genes 
      tx <- read.csv(file.path(base_dir,data_dir_1,pathway_lists_dir,used_pathway_file),header = F)$V1
      
      #filtering the gene mapping to be what we want
      mapped_genes <- ensembl_gene_mapping %>% dplyr::filter(hgnc_symbol %in% levels(tx))
      
      
    }
    
    #creating a vector whose values are ensembl ids and names are hgnc ids (doing this because this is how I originally wrote it/don't want to break it)
    entrez_to_ensembl <- mapped_genes$ensembl_gene_id
    names(entrez_to_ensembl) <- mapped_genes$hgnc_unique
    
    entrez_to_gene_name <- names(entrez_to_ensembl)
    
    
    #the actual retrieval of the genes we want from the expression matrix
    wanted_fpkm <- expression_data[,entrez_to_ensembl]
    
    #
    
    #labeling with hgnc gene code, making sure we have valid names
    colnames(wanted_fpkm) <- entrez_to_gene_name %>% make.names() %>% str_replace_all(pattern = "[.]",replacement = "_")
    
    #filtering to samples we want based on filtered metadata
    wanted_fpkm <- wanted_fpkm[dc_in$metadata_filtered$sample,]
    
    #removing NA values introduced in merging with other datasets
    wanted_fpkm <- wanted_fpkm[, colSums(is.na(wanted_fpkm)) == 0]
    
    
    #converting expression data to tibble because tibbles work well  
    fpkm_df <- dplyr::as_tibble(wanted_fpkm)
    fpkm_df$sample <- rownames(wanted_fpkm)
    
  } else { #loading data from tcga; should not use; not tested/is slow
    
    wanted_project <- tcga_projects %>% dplyr::filter(tumor == disease_filt) %>% dplyr::select(project_id) 
    rna_query <- GDCquery(project = wanted_project,
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - FPKM",
                          barcode = dc_in$metadata_filtered$sample)
    
    GDCdownload(rna_query, method = "api", files.per.chunk = 10)
    
    data <- GDCprepare(rna_query)  #using summarizedexperiment package
    
    data <- assay(data)
    colnames(data) <- colnames(data) %>% substr(1,15)
    
    wanted_fpkm <- data[entrez_to_ensembl,colnames(data) %in% dc_in$metadata_filtered$sample] %>% t()
    fpkm_df <- dplyr::as_tibble(wanted_fpkm)
    fpkm_df$sample <- rownames(wanted_fpkm)
  }
  
  #deleting expression data to make sure it's not taking up memory
  expression_data <-NULL
  dc_out <- dc_in
  dc_out$fpkm_df <- fpkm_df
  
  
  
  return(dc_out)
}

########################################
merge_gen_meta <- function(dc_in){
  #getting dimensions of metadata
  metadata_end <- dim(dc_in$metadata)[2]
  
  #joining metadata and expression data
  meta_fpkm_joined <- left_join(dc_in$metadata,dc_in$fpkm_df,by = "sample")
  
  # where are genetic data stored; this is a bad practice that I used in my matlab scripts
  expression_indices <- (metadata_end+1):NCOL(meta_fpkm_joined)
  
  #deleting samples with no expressiondata
  has_gen <-  apply(meta_fpkm_joined[,expression_indices], 2, function(x) !any(is.na(x)))
  meta_fpkm_joined <- meta_fpkm_joined %>% dplyr::select(-names(has_gen)[!has_gen])
  
  
  
  
  #saving results to dc_out
  dc_out <- dc_in
  dc_out$expression_indices <- expression_indices
  dc_out$meta_fpkm_joined <- meta_fpkm_joined
 
  
  return(dc_out)  
}

##################################################
normalize_data <- function(dc_in){
  
  dr_matrix <- data.matrix(dc_in$meta_fpkm_joined[,dc_in$expression_indices])
  
  
  if(standardize){  #make rows sum to 1 
    dr_matrix <- sweep(dr_matrix,1,rowSums(dr_matrix),"/")
    dr_matrix[is.na(dr_matrix)] <- 0
  }
  
  if(standardize_cols){ #make columns sum to 1
    dr_matrix <- sweep(dr_matrix,2,colSums(dr_matrix),"/")
    dr_matrix[is.na(dr_matrix)] <- 0
  }
  
  
  #centering mean and projecting onto unit hypersphere
  if(spherize){
    centroid <- colMeans(dr_matrix)
    dr_matrix <- sweep(dr_matrix,2,centroid,"-")
    radius <- sqrt(rowSums(dr_matrix^2))
    dr_matrix <- sweep(dr_matrix,1,radius,"/")
  }
  #are we
  
  #these Run transpose options are poorly tested; do not use until further development
  if( run_transpose){
    dr_matrix <- t(dr_matrix)}
  if(run_transpose){
    colnames(dr_matrix) <- dc_in$meta_fpkm_joined$sample
  } else {
    rownames(dr_matrix) <- dc_in$meta_fpkm_joined$sample 
  }
  
  dc_out <- dc_in
  dc_out$dr_matrix <- dr_matrix
  
  
  return(dc_out)
  
}

##############################################
dim_reduce <- function(dc_in){
  
  dr_matrix <- dc_in$dr_matrix
  
  if (use_umap){
    # selecting umap implementation and running umap; custom.config was set in the control panel
    umap_method <-    ifelse(use_umap_wrapper,  "umap-learn",  "naive")
    embedding <- umap(dr_matrix,config = custom.config, method = umap_method,verbose = TRUE)
    
    # extracting the dataframe with the umap output values, 
    layout_df <- data.frame(embedding$layout)
    
    
    
    
  } else{
    
    
    #options for using tsne if we're not using umap 
    if (use_tsne_wrapper){
      
      #could include the Rtsne_neighbors here
      
      embedding <-   Rtsne(dr_matrix, 
                           dims = output_dims,  
                           initial_dims = dim(dr_matrix)[2],
                           perplexity = tsne_perplexity, 
                           theta = tsne_theta, 
                           check_duplicates = tsne_check_duplicates,
                           pca = tsne_pca, 
                           partial_pca = tsne_partial_pca,
                           max_iter = tsne_max_iter,
                           verbose = tsne_verbose, 
                           is_distance = tsne_is_distance,
                           Y_init = tsne_Y_init,
                           pca_center = tsne_pca_center, 
                           pca_scale = tsne_pca_scale,
                           normalize = tsne_normalize,
                           stop_lying_iter = ifelse(is.null(tsne_Y_init), 
                                                    tsne_stop_lying_iter,0L), 
                           mom_switch_iter = ifelse(is.null(tsne_Y_init), tsne_mom_switch_iter, 0L),
                           momentum = tsne_momentum, 
                           final_momentum = tsne_final_momentum, 
                           eta = tsne_eta,
                           exaggeration_factor = tsne_exaggeration_factor, 
                           num_threads = tsne_num_threads)
      rownames(embedding$Y) <- rownames(dr_matrix)
      layout_df <- data.frame(embedding$Y)
      
    } else {
      embedding <- tsne::tsne(dr_matrix, 
                              initial_config = tsne_Y_init, 
                              k = output_dims, 
                              initial_dims = dim(dr_matrix)[2], 
                              perplexity = tsne_perplexity,
                              max_iter = tsne_max_iter, 
                              min_cost = tsne_min_cost, 
                              epoch_callback = tsne_epoch_callback, 
                              whiten = tsne_whiten,
                              epoch = tsne_epoch)
      rownames(embedding) <- rownames(dr_matrix)
      layout_df <- data.frame(embedding)
      
      
    }
    
    
  }
  
  #getting the filtered and combined clinical and expression data
  meta_fpkm_joined <- dc_in$meta_fpkm_joined 
  
  #joining the output of the dimension reduction to the 
  joined_dr <- bind_cols(meta_fpkm_joined,layout_df) 
  
  
  if(jitter_plot){
    #this section adds random noise to the clusters so that the points aren't all on top of each other. far better to rais the min_dist parameter of umap than to do this. in any case not sure if this is a good implementation.
    layout_df_jittered <- layout_df %>% mutate_if(is.numeric, jitter, amount = jitter_amount, factor = jitter_factor)
    joined_dr <- bind_cols(meta_fpkm_joined,layout_df_jittered) 
    joined_dr <- bind_cols(joined_dr, rename_all(layout_df,.funs = list(paste0),"_original"))
    dc_in$layout_df_jittered <- layout_df_jittered
  } else { 
    
    
  }
  
  dc_out <- dc_in
  dc_out$jitter_plot <- jitter_plot
  dc_out$embedding <- embedding
  dc_out$joined_dr <- joined_dr
  dc_out$layout_df <- layout_df
  
  return(dc_out)
}


####### this function runs clustering on the data
cluster <- function(dc_in){ 
  
  if(three_dee_output){ #extracting matrix to clust based on dimensions we put out
    
    to_clust <- cbind(dc_in$joined_dr$X1,dc_in$joined_dr$X2,dc_in$joined_dr$X3)
    colnames(to_clust) <- c("X1","X2","X3")} else{
      
      to_clust <- cbind(dc_in$joined_dr$X1,dc_in$joined_dr$X2)
      colnames(to_clust) <- c("X1","X2")
      
    }
  
  
  
  if(jitter_cluster){  #the idea behind this section is that we can choose between any combination of jittering the output plot and what the clustering algorithm sees.  AGAIN: better to just raise the min-dist parameter of umap
    
    if(jitter_plot){
      #we've already jittered the input vector
    } else{
      #we haven't jittered the input so we need to here
      to_clust %<>% jitter(amount = jitter_amount,factor = jitter_factor)}
  } else{  #we're not jittering the cluster
    
    if(jitter_plot){
      # we need to replace the clustering values by the un-jittered values
      if(three_dee_output){
        to_clust <- cbind(dc_in$joined_dr$X1_original,dc_in$joined_dr$X2_original,dc_in$joined_dr$X3_original)
        colnames(to_clust) <- c("X1","X2","X3")} else{
          
          to_clust <- cbind(dc_in$joined_dr$X1_original,dc_in$joined_dr$X2_original)
          colnames(to_clust) <- c("X1","X2")
        }
      
      
    } else{
      
      
    }
    
  }
  
  
  #actually running hdbscan
  dc_in$joined_dr$Cluster <- factor(paste("Cluster",
                                          hdbscan(to_clust,minPts = hdbscan_min_pts)$cluster)
  ) 
  
  
  if(force_hard_cluster){
    # use nearest_neighbors to force any "cluster 0" to join the nearest cluster 
    if(  any(dc_in$joined_dr$Cluster == "Cluster 0")){
      knn_df <- cbind(as_tibble(to_clust), as_tibble(dc_in$joined_dr[c("sample","Cluster")])) %>% as_tibble()
      
      #separating successful clusters from what was labeled as noise points
      cluster_zero <- knn_df %>% dplyr::filter(Cluster == "Cluster 0")
      successful_cluster <- knn_df %>% dplyr::filter(!(Cluster == "Cluster 0"))
      train <- successful_cluster %>% select_if(is.numeric) %>% as.matrix()
      test <-  cluster_zero %>% select_if(is.numeric)    %>% as.matrix()
      cl <- successful_cluster$Cluster
      
      #actually running k-nearest neighbors
      assigned_cluster <- knn(train,test,cl,k = k_for_knn)
      cluster_zero$Cluster <- assigned_cluster
      
      #assigning the results of knn to the cluster in the joined_dr tibble.  could probabably have done this and above steps with tidyverse
      dc_in$joined_dr$Cluster[ dc_in$joined_dr$sample %in% cluster_zero$sample  ] <- cluster_zero$Cluster
      
    }
  }
  
  #making cluster not a factor
  dc_in$joined_dr$Cluster <- droplevels(dc_in$joined_dr$Cluster)
  
  
  #explicitly recording the force_hard_cluster option
  dc_in$force_hard_cluster <- force_hard_cluster
  
  dc_out <- dc_in
  
  
  return(dc_out)
}


add_relation_points <- function(dc_in){  #makes ancestor and descendant tables for any input points given ordering of normal tissue, primary tumor, met; adds the actual points, whereas we've already computed which samples the ancestors and descendants are. 
  
  dc_out <- dc_in

  #these two functions make ancestor and descendant tables
  dc_out$ancestor_table_filtered   <- add_ancestor_points(dc_in$joined_dr,dc_in$ancestor_table_filtered)     %>% dplyr::filter(!( is.na(X1)    |is.na(X2)| is.na( X1_ancestor)| is.na(X2_ancestor)))
  dc_out$descendant_table_filtered <- add_descendant_points(dc_in$joined_dr,dc_in$descendant_table_filtered) %>% dplyr::filter(!( is.na(X1)    |is.na(X2)| is.na( X1_descendant)| is.na(X2_descendant)))
  
  
  return(dc_out) 
}


plot <- function(dc_in) {
  
  #getting the pathway files and their names from the dable describing them as such 
  pathway_files <- pathway_table$pathway_files
  names(pathway_files) <- pathway_table$pathway_names
  
  
  axis_label <- ifelse(use_umap,"UMAP","t-SNE")  #label axes with DR methods we used
  #getting values from dc in
  joined_dr <- dc_in$joined_dr
  dr_matrix <- dc_in$dr_matrix
  ancestor_table_filtered <- dc_in$ancestor_table_filtered
  
  #getting name of pathway, or if we sued a random sample
  if(use_rand_sample){
    pway_name_string <- paste(as.character(rand_sample_number), "Randomly Selected Genes")} else{
      pway_name_string <- names(pathway_files)[pathway_files %in% used_pathway_file] 
    }
  
   #variables to select for output graph
  merge_selection_vector <- c("sample","X1_ancestor", "X2_ancestor", "X3_ancestor","sample_ancestor","ancestor_type"  )
  
  #doing a join to make df for plot
  joined_dr_for_plot <- left_join(joined_dr, dplyr::select(ancestor_table_filtered, any_of(merge_selection_vector)),by = "sample") 
  
  #adding symbol size to tibble rather than leaving it just in environment
  joined_dr_for_plot$symbol_size <- symbol_size 
  
  #checking if there is only one thing, then we can make a title out of it because there is no need to label it
  is_one <- (apply(joined_dr_for_plot,2,function(x) {length(unique(x))} ) == 1) & !(colnames(joined_dr_for_plot) %in% colnames(dr_matrix)  )
  
  
  
  title <-(joined_dr_for_plot[,is_one]) #getting variables we are using for title
  title %<>% dplyr::select_if(~!all(is.na(.)))  #getting rid of things that are all the same because they are NA
  title <- dplyr::select(title, - !!((exclude_title_legend[exclude_title_legend %in% colnames(title)])))  #removing variables we explicitly said in control panel we don't want to use for tilte or have in legend
  
  
  title_string <- paste(axis_label, 
                        pway_name_string, 
                        apply( title[1,] , 1 , paste , collapse = " " )) #actually making base of title strings
  
  
  #Did we apply these normalizations to data? if so include in title
  if(spherize) { title_string <- paste(title_string, "Spherized")}  
  if(standardize) { title_string <- paste(title_string, "Standardized")}
  if(whiten) { title_string <- paste(title_string, "Whitened")}
  
  dc_out <- dc_in
  
  
  #getting color palettes and shortening them to have exactly the colors we will be using 
  pal1 <- readRDS(colors_file_1) %>% unname()
  pal2 <- readRDS(colors_file_2) %>% unname()
  pal3 <- readRDS(colors_file_3) %>% unname()
  
  pal_length_1 <- factor(joined_dr_for_plot[[quo_name(color_selection_1)]]) %>% levels() %>% length()
  pal_length_2 <- factor(joined_dr_for_plot[[quo_name(color_selection_2)]]) %>% levels() %>% length()
  pal_length_3 <- factor(joined_dr_for_plot[[quo_name(color_selection_3)]]) %>% levels() %>% length()
  
  pal1 <- pal1[1:pal_length_1]
  pal2 <- pal2[1:pal_length_2]
  pal3 <- pal3[1:pal_length_3]
  
  dc_out$three_dee <- three_dee_output
  
  
  env_list  <- as.list.environment(environment())
  if(three_dee_output){
    
    
    
    if(plot_black_rings){  #should not use this option; was hoping to have option to choose between black rings indicating met or hollow circle indicating met; turns out the black rings are are hard to do well in plotly; could try it but the other way is better
      dc_out$scatter_color1 <- plot_black_rings_3d(color_slxn = color_selection_1 , palette =  pal1, data_in = joined_dr_for_plot)
      dc_out$scatter_color2 <- plot_black_rings_3d(color_slxn = color_selection_2 , palette =  pal2, data_in = joined_dr_for_plot)
      dc_out$scatter_color3 <- plot_black_rings_3d(color_slxn = color_selection_3 , palette =  pal3, data_in = joined_dr_for_plot)
    } else {
      dc_out$scatter_color1 <- plot_hollow_rings_3d(color_slxn = color_selection_1 , palette =  pal1, data_in = joined_dr_for_plot)
      dc_out$scatter_color2 <- plot_hollow_rings_3d(color_slxn = color_selection_2 , palette =  pal2, data_in = joined_dr_for_plot)
      dc_out$scatter_color3 <- plot_hollow_rings_3d(color_slxn = color_selection_3 , palette =  pal3, data_in = joined_dr_for_plot)
      
    }
    
    
  } else{   #did 2d dimension reduction 
    #have caption string
    caption_string <- paste0("Parameters: n-Neighbors = ",custom.config$n_neighbors ,", Min. Dist. = ",custom.config$min_dist,", n-Epochs = ",custom.config$n_epochs)
    
    if(plot_black_rings){
      dc_out$scatter_color1 <- plot_black_rings_2d(color_slxn = color_selection_1, palette = pal1,data_in = joined_dr_for_plot)
      dc_out$scatter_color2 <- plot_black_rings_2d(color_slxn = color_selection_2, palette = pal2,data_in = joined_dr_for_plot)
      dc_out$scatter_color3 <- plot_black_rings_2d(color_slxn = color_selection_3, palette = pal3,data_in = joined_dr_for_plot)
    } else{
      
      dc_out$scatter_color1 <- plot_hollow_rings_2d(color_slxn = color_selection_1 , palette =  pal1,data_in = joined_dr_for_plot)
      dc_out$scatter_color2 <- plot_hollow_rings_2d(color_slxn = color_selection_2 , palette =  pal2,data_in = joined_dr_for_plot)
      dc_out$scatter_color3 <- plot_hollow_rings_2d(color_slxn = color_selection_3 , palette =  pal3,data_in = joined_dr_for_plot)
    }
    
    
    
  }
  
  
  
  dc_out$scatter_color1_lines <- dc_out$scatter_color1
  dc_out$scatter_color2_lines <- dc_out$scatter_color2
  dc_out$scatter_color3_lines <- dc_out$scatter_color3
  
  if(draw_lines){
    #drawing the lines connecting points to their ancestors and descendants
    dc_out$scatter_color1_lines <- plot_lines_proper(joined_dr_for_plot,dc_out$scatter_color1_lines)
    dc_out$scatter_color2_lines <- plot_lines_proper(joined_dr_for_plot,dc_out$scatter_color2_lines)
    dc_out$scatter_color3_lines <- plot_lines_proper(joined_dr_for_plot,dc_out$scatter_color3_lines)
    
  }
  
  dc_out$title_string <- title_string
  
  return(dc_out)}


## running the post_hoc frequency analysis 
chi_square <- function(dc_in){
  
  
  joined_dr <- dc_in$joined_dr
  
  
  if(filter_for_stats){  #applying fiter we use before runing chi square
    joined_dr <- joined_dr %>% dplyr::filter(!!stats_filter_var != stats_filter_elim) 
  }
  
  #doing same thing we did to the palette before
  pal1 <- readRDS(colors_file_1) %>% unname()
  pal_length_1 <- factor(joined_dr[[quo_name(chi_sq_var_quo)]]) %>% levels() %>% length()
  pal1 <- pal1[1:pal_length_1]
  
  
  #similar title things to before
  is_one <- (apply(joined_dr,2,function(x) {length(unique(x))} ) == 1) & !(colnames(joined_dr) %in% colnames(dc_in$dr_matrix)  )
  axis_label <- ifelse(use_umap,"UMAP","t-SNE")
  #pway_name_string <- names(pathway_files)[pathway_files %in% used_pathway_file] 
  
  
  #making contingency table of cluster and the chosen variable
  tbl <- table(as.character(joined_dr$Cluster),unlist(lapply(joined_dr[[chi_sq_var]],as.character)))
    title_string <- dc_in$title_string

  
  #converting the above table to matrix
  tbl_mat <- as.matrix(tbl)
  
  #running chi square on the table
  cs <- chisq.test(tbl, correct = TRUE)
  
  #choosing between fisher exact post hoc or chi square resitdual post hoc
  if(fisher) {
    post_hoc <- matrix(NA,dim(tbl_mat)[1],dim(tbl_mat)[2])
    rownames(post_hoc) <- rownames(tbl_mat)
    colnames(post_hoc) <- colnames(tbl_mat)
    
    
    for (e in 1:dim(tbl_mat)[1]){
      for (f in 1:dim(tbl_mat)[2]){
        collapsed <- contingency_collapse(tbl_mat,e,f)
        
        ft_out <- fisher.test(collapsed)
        post_hoc[e,f] <- ft_out$p.value
        
        
      }
    }
    post_hoc_fisher <- post_hoc
    dc_in$post_hoc_fisher <- post_hoc_fisher
  } else {
    post_hoc_res <- cs$stdres
    rownames(post_hoc_res) <- rownames(tbl_mat)
    colnames(post_hoc_res) <- colnames(tbl_mat)
    post_hoc <- pchisq((post_hoc_res)^2,df = 1,lower.tail = FALSE)
    
    dc_in$post_hoc_chisq <- post_hoc
  }
  
  
  #making table for plotting, whereas the other one was for the statistical test
  freq_tab <- dplyr::count(joined_dr,(!!sym(chi_sq_var)),Cluster)
  
  
  
  freq_tab$fisher_exact_p <- post_hoc[cbind( as.character(freq_tab$Cluster),as.character(freq_tab[[chi_sq_var]]))]  #appending the p values, calling it fisher even tbough it could be the residual method
  fish_sig <- freq_tab$fisher_exact_p < print_p_threshold  #do we consider the post-hoc p vals significant
  
  freq_tab %<>% mutate(fisher_exact_p = paste("P =",formatC(.data$fisher_exact_p, format = "e", digits = 2)))  #making the p-values nice to disply
  freq_tab$fisher_exact_p[!fish_sig] <- ""  #deleting strings of non-significant p values
  freq_tab %<>% mutate_if(is.character,tools::toTitleCase)  #makign things title case
  
  title_string <- paste0(" Cluster Composition of ",title_string)
  cs_format <- formatC(cs$p.value, format = "e", digits = 2)  #making p vals print nicely
  caption_string <- paste0("Chi-sq P = ",cs_format,"; P < ",as.character(print_p_threshold) ," displayed above")  #making string for captin
  

  #making actual bar plot
  cluster_comp_bp <- ggplot(freq_tab,aes(x = Cluster,y = n,label = fisher_exact_p ,fill = !! chi_sq_var_quo))+scale_fill_manual(values = toupper(pal1))+
    labs(title = c(title_string),x = "Cluster",y = "Frequency",caption = caption_string, fill =  make_title_case(quo_name(chi_sq_var_quo) )) + 
    geom_bar(stat = "identity") + theme(aspect.ratio=1)
  
  
  if(show_p_vals){
    
    #if we show the p vals, show the p vals
    cluster_comp_bp_pvals <- cluster_comp_bp + geom_text(size = 3, position = position_stack(vjust = 0.5))
    
  }
  
  dc_out <- dc_in
  dc_out$freq_tab <- freq_tab
  dc_out$cluster_comp_bp_pvals <-   cluster_comp_bp_pvals
  dc_out$cluster_comp_bp <-  cluster_comp_bp
  dc_out$chi_square <- cs
  return(dc_out)
  
}


run_rfc <- function(dc_in){   #running random forest classifier; have yet to understand this 
  
  dr_matrix <- dc_in$dr_matrix
  joined_dr <- dc_in$joined_dr
  
  gen_df <- as_tibble(dr_matrix, rownames = "sample")
  info_df <- joined_dr[,c("sample","Cluster")]
  rfc_df <- left_join(gen_df,info_df,by = "sample",keep = FALSE) %>% dplyr::select(-sample) %>% dplyr::filter(Cluster != 0)
  
  hyper_grid <- expand.grid(   #grid of parameters to do with random forest
    mtry       = rfc_mtry,
    node_size  = seq(2, 20, by = 2),
    sampe_size = c(.55, .632, .70, .80),
    OOB_RMSE   = 0
  )
  
  

# rfc_running: legitimately for get ow this works -------------------------------------------------------------

  
  for(i in 1:nrow(hyper_grid)) {
    
    
    # train model
    model <- ranger(
      formula         = Cluster ~ ., 
      data            = rfc_df, 
      num.trees       = 500,
      mtry            = hyper_grid$mtry[i],
      min.node.size   = hyper_grid$node_size[i],
      sample.fraction = hyper_grid$sampe_size[i],
      seed            = 42069
    )
    
    # add OOB error to grid
    hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
  }
  
  RMSE <- hyper_grid$OOB_RMSE
  min_error <- RMSE == min(RMSE)
  params <- hyper_grid[min_error,][1,]
  
  
  OOB_RMSE <- vector(mode = "numeric", length = 100)
  
  for(i in seq_along(OOB_RMSE)) {
    
    rfc_out <- ranger(
      formula         = Cluster ~ ., 
      data            = rfc_df, 
      num.trees       = num_trees,
      mtry            = params$mtry,
      min.node.size   = params$node_size,
      sample.fraction = params$sampe_size,
      importance      = 'impurity'
    )
    
    OOB_RMSE[i] <- sqrt(rfc_out$prediction.error)
    
    
    
    
  }
  
  

# plotting rfc: variable importance ------------------------------------------------------------

  
  rfc_plot_data <-     rfc_out$variable.importance %>% 
    tidy() %>%
    dplyr::arrange(desc(x)) %>%
    dplyr::top_n(15) 
  title_string <- dc_in$title_string
  colnames(rfc_plot_data) <- c("Gene","Importance")
  rfc_plot_data %<>% arrange(-Importance)
  caption_string <- paste0("Parameters: n-trees = ",num_trees ,", Mtry  = ",params$mtry  ,", Min Node Size = ",params$node_size, ", Sample Fraction = ",params$sampe_size )
  rfc_plot <-   ggplot(rfc_plot_data,aes(x = reorder(Gene, Importance), y = Importance)) +
    geom_col() +
    coord_flip()  + labs(title = paste(title_string,"Random Forest Gene Importance"),y = "Predictive Importance", x = "Gene",
                         caption = caption_string)+theme(aspect.ratio=1,text = element_text(size=font_size))+theme(plot.title = element_text(hjust = 0.5))  + theme( text = element_text(size=font_size))
  dc_out <- dc_in
  dc_out$rfc_plot <- rfc_plot
  dc_out$rfc_out <- rfc_out
  
  return(dc_out)
  
  
  
  
}




make_heatmap <- function(dc_in){
  
  joined_dr <- dc_in$joined_dr
  mat <- dc_in$dr_matrix  #matrix for heatmap
  rn <- rownames(mat)    #rownames for heatmap
  
  mat <- mat[order(match(rn,joined_dr$sample)),]  #reorderingmat to match joined dr
  
  rownames(mat) <- joined_dr$patient_id  #replacing sample id with patient id
  
  wanted_variables <-   dc_in$rfc_out$variable.importance %>% 
    tidy() %>%
    dplyr::arrange(desc(x)) %>%
    dplyr::top_n(vars_for_heatmap) %>% pull(names)   #vars for heatmap is a number of variables wer are displaying
  mat <- mat[,colnames(mat) %in% wanted_variables]
  
  #manipulating palettes as we have before; the names of the palettes are important
  pal_length_1 <- factor(joined_dr[[quo_name(color_selection_1)]]) %>% levels() %>% length()
  pal1 <- readRDS(colors_file_1) %>% unname()
  pal1 <- pal1[1:pal_length_1]
  names(pal1) <- factor(joined_dr[[quo_name(color_selection_1)]]) %>% levels()
  
  pal_length_3 <- factor(joined_dr[[quo_name(color_selection_3)]]) %>% levels() %>% length()
  pal3 <- readRDS(colors_file_3) %>% unname() %>% rev()
  pal3 <- pal3[1:pal_length_3]
  names(pal3) <- factor(joined_dr[[quo_name(color_selection_3)]]) %>% levels()
  
  #making the palettee upper case for some reason
  pal3 %<>% toupper()
  pal1 %<>% toupper()
  
  #makign heatmap annotations aout of sample type and target tissue.  could make these variables something you set in control panel to add flexibity
  ha = HeatmapAnnotation(
    `Sample Type` = joined_dr$sample_type,
    `Target Tissue` = joined_dr$target_site,
    col = list( `Sample Type` = pal3,
                `Target Tissue` = pal1
    ),show_legend = TRUE)
  
  
  dc_out <- dc_in  
  
  #setting something to split bty
  split <- as.character(joined_dr$Cluster)
  

  #making two kinds of heamap
  dc_out$heatmap_rowmean <-  Heatmap(t(scale(mat)),name = "SD From Row Mean",column_split = joined_dr$Cluster, top_annotation = ha)
  dc_out$heatmap_colmean <-  Heatmap(scale(t(mat)),name = "SD From Col. Mean",column_split = split, top_annotation = ha)
  
  
  return(dc_out)
  
  
}



volcano <- function(dc_in){ #makes volvano plot
  
  dc_out <- dc_in
  
  mat <- dc_in$dr_matrix 
  mat <- mat - min(mat) #subtracting the minimum value of the matrix for some reason 
  mat <- mat  %>% as_tibble(rownames = "sample")
  
  
  mat_joined <- left_join(dplyr::select(dc_in$joined_dr,sample,Cluster),mat, by = "sample") %>% ungroup()  #adding cluster to the df we're using for volcano
  
  
  full_vec <- mat_joined %>% select_if(is.numeric) %>% as.matrix() %>% c()   #
  increment_value <- min(full_vec[full_vec>0])/100             #finding value to increment by so we can take log of zero vals
  mat_joined <- mat_joined %>% mutate_if(is.numeric, function(x){log2(x + increment_value)})  #log transforming
  mat_joined_nosample <- mat_joined %>% dplyr::select(-sample)   #getting rid of sample
  mat_joined_nosample <- mat_joined_nosample %>% mutate( Cluster = as.integer(str_sub(as.character(Cluster),-1)))  $#chainging cluster value for some reaon
  mat_joined_nosample_nocluster <- mat_joined_nosample %>% dplyr::select(-Cluster)
  
  
  
  #running aov on every gene and outputing 
  aov.models <- mat_joined_nosample_nocluster %>% map(~ aov(mat_joined_nosample$Cluster ~ .x))
  p_vals <- lapply(aov.models,extract_anova_pvalue) %>% as_tibble() %>% pivot_longer(cols = everything(), names_to = "Gene", values_t = "p_val")
  
  
  #adding wor cluter
  mat_joined_nosample <-                               mat_joined_nosample %>%
    mutate(Cluster = paste("Cluster",Cluster, sep = "_"))
  
  
  #making CD so that we could plot horizontal error bars
  gene_SD <-       mat_joined_nosample %>%
    group_by(Cluster) %>% 
    summarize_all(function(x) sd(x)/sqrt(length(x))) %>% 
    pivot_longer(cols = -Cluster, names_to = "Gene", values_to = "SD") %>%
    pivot_wider(names_from = Cluster,values_from = SD) %>% pivot_longer(cols = -Gene,names_to = "Cluster", values_to = "SD")
  
  #pivoting longer for some reason; need to check why i did this
  longer_mat_joined <- mat_joined_nosample %>% pivot_longer(cols = -Cluster, names_to = "Gene", values_to = "value")                  
  #processing to subtract logs by mean (same as taking log of ratio(
  l2fc_mat <-                                             longer_mat_joined %>%  
    group_by(Gene) %>%
    mutate(mean_gene_log2 = mean(value))   %>%
    group_by(Gene,Cluster) %>% 
    summarize(mean_gene_log2 = mean(mean_gene_log2), mean_gene_cluster_log2 = mean(value)) %>% 
    ungroup() %>% 
    mutate(l2fc = mean_gene_cluster_log2 - mean_gene_log2)                             
  
  
  full_volcano_data <- left_join(gene_SD, l2fc_mat) %>% left_join(p_vals)  #joining p vals
  full_volcano_data <- full_volcano_data %>% mutate(upper_val = l2fc + SD, lower_val = l2fc - SD, pval_neglog10 = - log10(p_val))  #adding upper and lower error bar points
  full_volcano_data <- full_volcano_data %>% mutate(   Cluster = str_replace_all(string = Cluster ,pattern = "[_]", replacement = " ")        )   # changing format of cluster string
  
#processing color palette as we did before  
  pal_length_2 <- factor(full_volcano_data[[quo_name(color_selection_2)]]) %>% levels() %>% length()
  pal2 <- readRDS(colors_file_2) %>% unname()
  pal2 <- pal2[1:pal_length_2]
  names(pal2) <- factor(full_volcano_data[[quo_name(color_selection_2)]]) %>% levels() 
  
  
  
  #getting rightmost point for each geneto add label
  full_volcano_data_rightmost <- full_volcano_data %>% group_by(Gene) %>% dplyr::filter(upper_val == max(upper_val))
  
  volcano_text_threshold
  
  #making volcano plot
  volcano <-  full_volcano_data %>% ggplot(mapping = aes(y = pval_neglog10, x = l2fc, color = Cluster, xmin = lower_val, xmax = upper_val, label = Gene)) +
    geom_point() + 
    geom_text_repel(data = dplyr::filter(full_volcano_data_rightmost,pval_neglog10 > 4.5) , color = "black",mapping = aes(y = pval_neglog10, x = l2fc))+ 
    scale_color_manual(values = toupper(pal2)) + 
    labs(x = expression(
      paste(Log[2],"(Ratio to Mean)")), y = expression(
        paste("-",Log[10],"(P Value)")), title = paste("Volcano Plot:", dc_in$title_string))
  
  
  
  dc_out$volcano <- volcano
  
  
  return(dc_out)
  
}






save_output <- function(dc_in){ #saves output and makes timestamp string
  
  dc_out <- dc_in
  if(save_output){
    time_now <- ymd_hms(now()) %>% str_replace_all("[-:\\s]","_")
    file_name <- paste0(c("output",time_now),  collapse = "_")
    file_name <- paste0(file_name, ".RDS")
    saveRDS(object = dc_out,file = file.path(base_dir,output_dir,file_name))
  }
  return(dc_out)
}


