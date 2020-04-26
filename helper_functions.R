# takes a vector and checks for duplicates, INCLUDING the first instance
all_dup <- function (value){
  duplicated(value) | duplicated(value, fromLast = TRUE)
} 




#collapses contingency table to be a 2x2  yes/no contingenccy for one category specified by row and column
contingency_collapse <- function(table,r,c){
  out_mat <- matrix(NA,2,2)
  out_mat[1,1]  <- table[r,c]
  out_mat[1,2]  <- sum(table[r,])-table[r,c]
  out_mat[2,1]  <- sum(table[,c])-table[r,c]
  out_mat[2,2]  <- sum(table) -(out_mat[1,1] + out_mat[1,2] +out_mat[2,1])
  return(out_mat)
}


#exponentiates base two and decrements by 1
exp_dec <- function(x) {
  (2^x)-1
}


#Loads genetic  data we want only if it is not already loaded
#is a function but used like a script
load_gen_data <- function(){
  if(!exists("expression_data")){
    print("data not loaded")
    
    if(use_tpm){ 
      expression_data <<- readRDS(file.path(base_dir,data_dir,tpm_file))
      tpm_loaded <<- TRUE  #asignment operator traces back to first parent environment with 'tpm_loaded' until base environment
    } else {
      expression_data <<- readRDS(file.path(base_dir,data_dir,fpkm_file))
      tpm_loaded <<- FALSE 
    } 
    
    
  } else {   
    print("data are loaded")
    if(use_tpm){
      if(tpm_loaded){
      }else{
        expression_data <<- readRDS(file.path(base_dir,data_dir,tpm_file))
        tpm_loaded <<- TRUE
      }
    } else {
      if(tpm_loaded){
        expression_data <<- readRDS(file.path(base_dir,data_dir,fpkm_file))
        tpm_loaded <<- FALSE
      }else{
      }
      
    }  
  }
  return(expression_data)
}





load_gen_data_simple <- function(){


    
    if(use_tpm){ 
      expression_data <- readRDS(file.path(base_dir,data_dir,tpm_file))
      tpm_loaded <- TRUE  #asignment operator traces back to first parent environment with 'tpm_loaded' until base environment
    } else {
      expression_data <- readRDS(file.path(base_dir,data_dir,fpkm_file))
     
    } 
    
    
     


  
  return(expression_data)
}





#gets the ggplot color hues but excludes the inital red hue, useful for using it for noise points
ggcolor_hue_minus <- function(n) {
  hues = seq(15, 375, length = n + 1)
  full <- hcl(h = hues, l = 65, c = 100)[1:n]
  minus_one <-full[2:numel(full)]
  return(minus_one)
}


make_title_case <- function(str){
  str <- sub(pattern = "_",replacement = " ",fixed = TRUE, x = str)
  str <- toTitleCase(str)
}






#defines descentandt of each tissue type
descendant <- function(str_in){
  
  switch_fun <- function(x) { switch(x, "Solid Tissue Normal" = {"Primary Tumor"},  
                                     "Primary Tumor" =   {"Metastatic"},
                                     "Recurrent Tumor" = {"Metastatic"},
                                     "Metastatic" = {NA},
                                     NA)
  }
  return(sapply(str_in,  switch_fun))
  
}


#ancestor
ancestor <- function(str_in){
  
  switch_fun <- function(x) { switch(x, "Solid Tissue Normal" = {NA},  
                                     "Primary Tumor" =   {"Solid Tissue Normal"},
                                     "Metastatic" = {"Primary Tumor"},
                                     "Recurrent Tumor" = {"Primary Tumor"},
                                     NA)
  }
  return(sapply(str_in,  switch_fun))
  
}



make_tissue_factor <- function(tissue_string_in) {
  out_factor <-  factor(tissue_string_in, levels=c("Solid Tissue Normal", "Primary Tumor", "Recurrent Tumor","Metastatic"), ordered=TRUE)
}


#adding ancestors and descendant tables
add_ancestor_points <- function(data_in,ancestors_in){

  ancestors_in <- ancestors_in %>% dplyr::select(-any_of(c("X1", "X2", "X3","X1_ancestor", "X2_ancestor", "X3_ancestor")))

  sub_data_in <- data_in %>% dplyr::select(sample, any_of(c("X1", "X2", "X3")))
  sub_data_in_ancestor <- sub_data_in %>% rename_all(.funs = ~paste0(.,"_ancestor"))
  
  ancestors_points <- ancestors_in %>% left_join( sub_data_in) %>% left_join(sub_data_in_ancestor) 
  
  
  
  return(ancestors_points)
  
}


add_descendant_points <- function(data_in, descendants_in){
  descendants_in <- descendants_in %>% dplyr::select(-any_of(c("X1", "X2", "X3","X1_ancestor", "X2_ancestor", "X3_ancestor")))
  sub_data_in <- data_in %>% dplyr::select(sample, any_of(c("X1", "X2", "X3")))
  sub_data_in_descendant <- sub_data_in %>% rename_all(.funs = ~paste0(.,"_descendant"))
  
  descendants_points <- descendants_in %>% left_join( sub_data_in) %>% left_join(sub_data_in_descendant) 
  
  
  
  
  return(descendants_points)
}

#these versions of make descendant table do not  deal with the x dimension
make_ancestor_table_meta <- function(data_in){
  
  
  
  
  starting_table <- data_in %>% transmute(sample_type = sample_type , ancestor_type = ancestor(sample_type)   , patient_id = patient_id, sample = sample, disease = disease , dataset = dataset )
  ancestor_table  <- data_in %>% transmute(ancestor_type = sample_type, patient_id = patient_id, sample_ancestor = sample) 
  
  
  ancestor_table = left_join(starting_table, ancestor_table, by = c("ancestor_type", "patient_id"))
  
  return(ancestor_table)
  
}


make_descendant_table_meta <- function(data_in){
  
  starting_table <- data_in %>% transmute(sample_type = sample_type , descendant_type = descendant(sample_type)   , patient_id = patient_id, sample = sample, disease = disease, dataset  = dataset  )
  descendant_table  <- data_in %>% transmute(descendant_type = sample_type, patient_id = patient_id, sample_descendant = sample) 
  
  descendant_table = left_join(starting_table, descendant_table, by = c("descendant_type", "patient_id"))
  
  return(descendant_table)
  
  
  
}




get_ancestor_coordinates <- function(data_in){
  
  starting_table <- data_in %>% transmute(sample_type = sample_type , ancestor_type = ancestor(sample_type)   , patient_id = patient_id, sample = sample  )
  ancestor_table  <- data_in %>% transmute(ancestor_type = sample_type, patient_id = patient_id, sample_ancestor = sample) 
  
  
  ancestor_table = left_join(starting_table, ancestor_table, by = c("ancestor_type", "patient_id"))
  
  return(ancestor_table)
  
}


get_descendant_coordinates <- function(data_in){
  
  starting_table <- data_in %>% transmute(sample_type = sample_type , descendant_type = descendant(sample_type)   , patient_id = patient_id, sample = sample  )
  descendant_table  <- data_in %>% transmute(descendant_type = sample_type, patient_id = patient_id, sample_descendant = sample) 
  
  descendant_table = left_join(starting_table, descendant_table, by = c("descendant_type", "patient_id"))
  
  return(descendant_table)
  
  
  
}





get_target_site <- function(data_in){
  
  temp_target <- data_in %>% dplyr::select(patient_id, sample_type, biopsy_tissue, sample)  
  
  #getting latest known target tissue site for each sample
  
  temp_target <- temp_target %>% 
    mutate(sample_type = make_tissue_factor(sample_type)) %>% 
    group_by(patient_id) %>%
    mutate(target_sample_type = max(sample_type)) %>% 
    group_by(patient_id) %>% 
    dplyr::filter(sample_type  == max(target_sample_type)) %>% 
    dplyr::rename( target_site = biopsy_tissue ) %>% 
    dplyr::select(-c(target_sample_type,sample_type,sample)) %>% distinct()
  
  data_in <- left_join(data_in, temp_target, by = "patient_id")
  return(data_in)
}












plot_mets_proper <- function(joined_dr, scatter_in){
  #  joined_dr <- joined_dr %>% filter(sample_type == "Metastatic")

  
    
    

    
   
    sub_df <- joined_dr %>% dplyr::filter(sample_type == "Metastatic")

    
    
    
        
        scatter_in <-    scatter_in %>% add_trace(data = sub_df, x = ~X1 , y = ~X2, z = ~X3,
                                                  marker = list(color = "rgba(0,0,0,0)" ,symbol = "100",symbols = "100", shape ="100",size = 9,line = list(
                                                    color = 'rgb(0,0,0)',
                                                    width = 3
                                                  )) , 
                                                  name = "metastatic", type = "scatter3d",mode = "marker", showlegend = TRUE)
        
  
   

  return(scatter_in)
}




#### plotting the lines in 3d  #would like to get rid of for loops but it could make the legend situation difficult with plotlyh 



plot_lines_proper <- function(table_in, scatter_in){

  
 
 remove_num <-  function(x) gsub('[0-9]+', '', x)
 

    
    if(three_dee_output){
      
      lines_table  <- table_in %>%  dplyr::select( any_of(c("sample", "sample_ancestor", "X1_ancestor",   "X2_ancestor", "X3_ancestor"  , "X1","X2","X3")) ) %>%
        dplyr::filter(!is.na(sample_ancestor)) %>% mutate(pairing = paste(sample, sample_ancestor, sep = "|")) %>% pivot_longer(cols = c("X1", "X1_ancestor"), names_to = "x1_type", values_to = "X") %>%
        pivot_longer(cols = c("X2", "X2_ancestor"), names_to = "x2_type", values_to = "Y") %>%
        pivot_longer(cols = c("X3", "X3_ancestor"), names_to = "x3_type", values_to = "Z") %>% mutate_at(.vars = c("x1_type","x2_type","x3_type"), .funs =  remove_num) %>% dplyr::filter((x1_type == x2_type) & (x1_type == x3_type) & (x2_type == x3_type)) %>% group2NA(groupNames = "pairing")
    
   
  
    scatter_in <-  scatter_in %>% add_paths(data = lines_table , x = ~X , y = ~Y, z = ~Z,
                                            line = list(color = '#000000', width = 3, opacity = 0) , 
                                            name = "Same Patient", type = "scatter3d",mode = "lines", showlegend = TRUE)
    
    
  } else {
 
    scatter_in <- scatter_in + geom_segment(data = table_in, mapping = aes(x = X1_ancestor,y = X2_ancestor,
                                                                           xend = X1 ,yend = X2, linetype = "")) + 
      labs(linetype = "Same Patient") +
      guides( color = guide_legend(order = 1),shape = guide_legend(order = 2), linetype = guide_legend(order = 3))
  }
  
  
  
  
  
  return(scatter_in)
}



extract_anova_pvalue <- function(anova_in){
  summary(anova_in)[[1]][["Pr(>F)"]][1]
}






plot_black_rings_2d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in) { 
  
  plot <- ggplot(data_in,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color = !!(border_selection) ,fill= !! (color_slxn)),shape = 21, stroke = 2) + 
    labs(title = env$title_string,  y = paste(env$axis_label,'2'), x= paste(env$axis_label,'1') ,caption = env$caption_string, color = make_title_case(quo_name(border_selection)) ,fill = make_title_case(quo_name(color_slxn)))+ 
    theme(aspect.ratio=1,text = element_text(size=font_size)) + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = toupper(palette)) + scale_color_manual(values= c("#000000FF","#FFFFFF00","#FF00FF"))
  return(plot)
}



plot_hollow_rings_2d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in) { 

  plot <- ggplot(data_in,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color = !!(color_slxn) ,shape = !! (border_selection)), stroke = 2) + 
    labs(title = env$title_string,  y = paste(env$axis_label,'2'), x= paste(env$axis_label,'1') ,caption = env$caption_string, shape = make_title_case(quo_name(border_selection)) ,color = make_title_case(quo_name(color_slxn)))+ 
    theme(aspect.ratio=1,text = element_text(size=font_size)) + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = toupper(palette)) + scale_color_manual(values = toupper(palette)) +scale_shape_manual(values = shape_scale)
  return(plot)
}





plot_black_rings_3d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in){

  out <- plot_ly() %>%
    add_markers( data = data_in, x = ~X1,marker = list( sizemode = 'diameter'), y = ~X2, z = ~X3, color = color_slxn, colors = palette,opacity = .8,size = ~symbol_size ,sizes = c(symbol_size*symbol_size_factor,symbol_size*symbol_size_factor)) %>%
    layout(title = paste0(env$title_string,": ",make_title_case(quo_name(color_slxn))),scene = list(xaxis = list(title = paste(env$axis_label,'1')),
                                                                                                    yaxis = list(title = paste(env$axis_label,'2')),
                                                                                                    zaxis = list(title = paste(env$axis_label,'3'))))

  out<- plot_mets_proper(data_in,out)
  return(out)
}




plot_hollow_rings_3d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in){
  
  used_shape_scale <- shape_scale_3d[names(shape_scale_3d) %in% dplyr::pull(data_in,(quo_name(border_selection)))]
  
  out <- plot_ly() %>% 
    add_markers( data = data_in, marker = list( sizemode = 'diameter'),x = ~X1, y = ~X2, z = ~X3, color = color_slxn, colors = palette,opacity = .8, symbol = border_selection ,symbols = used_shape_scale, size = ~symbol_size, sizes = c(symbol_size*symbol_size_factor,symbol_size*symbol_size_factor)) %>%
    layout(title = paste0(env$title_string,": ",make_title_case(quo_name(color_slxn))),scene = list(xaxis = list(title = paste(env$axis_label,'1')),
                                                                                                    yaxis = list(title = paste(env$axis_label,'2')),
                                                                                                    zaxis = list(title = paste(env$axis_label,'3'))))

  
  return(out)
}











