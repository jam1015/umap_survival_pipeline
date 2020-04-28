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
#is a function but used like a script; the <<- asigns to a variable in the global environment
#the idea is that it also appropriately loads tpm or fpkm 
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




#loads the genetic data in a more simple way. probably better to use this. 
load_gen_data_simple <- function(){
    if(use_tpm){ 
      expression_data <- readRDS(file.path(base_dir,data_dir,tpm_file))
      tpm_loaded <<- TRUE  #asignment operator traces back to first parent environment with 'tpm_loaded' until base environment
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

#makies things title cae and does a replacement
make_title_case <- function(str){
  str <- sub(pattern = "_",replacement = " ",fixed = TRUE, x = str)
  str <- toTitleCase(str)
}






#defines descentant of each tissue type
descendant <- function(str_in){
  
  switch_fun <- function(x) { switch(x, "Solid Tissue Normal" = {"Primary Tumor"},  
                                     "Primary Tumor" =   {"Metastatic"},
                                     "Recurrent Tumor" = {"Metastatic"},
                                     "Metastatic" = {NA},
                                     NA)
  }
  return(sapply(str_in,  switch_fun))
  
}


# decines ancestor of each tissue type
ancestor <- function(str_in){
  switch_fun <- function(x) { switch(x, "Solid Tissue Normal" = {NA},  
                                     "Primary Tumor" =   {"Solid Tissue Normal"},
                                     "Metastatic" = {"Primary Tumor"},
                                     "Recurrent Tumor" = {"Primary Tumor"},
                                     NA)
  }
  return(sapply(str_in,  switch_fun))
  
}


#makes an ordered factor out of a string with tissue types, in rough order of ancestor ---> descendant
make_tissue_factor <- function(tissue_string_in) {
  out_factor <-  factor(tissue_string_in, levels=c("Solid Tissue Normal", "Primary Tumor", "Recurrent Tumor","Metastatic"), ordered=TRUE)
}


#adding ancestors and descendant tables: data in is data in is the metadata + output dimiensionally reduced points, ancestors in is an ancestor table
add_ancestor_points <- function(data_in,ancestors_in){

   ancestors_in <- ancestors_in %>% dplyr::select(-any_of(c("X1", "X2", "X3","X1_ancestor", "X2_ancestor", "X3_ancestor"))) #these data are not present, but getting rid of them just in case
  sub_data_in <- data_in %>% dplyr::select(sample, any_of(c("X1", "X2", "X3")))  #making a copy of the data in  points with sample labels
  sub_data_in_ancestor <- sub_data_in %>% rename_all(.funs = ~paste0(.,"_ancestor"))  #adding ancestor label
 
  #the way this works is left joining in sequence: sub_data_in to ancestor_in by shared column sample, and that joined tibble to sub_data_in_ancestor by shared column "sample_ancestor"

   ancestors_points <- ancestors_in %>% left_join( sub_data_in) %>% left_join(sub_data_in_ancestor) 

  return(ancestors_points)
}

#similar to add_ancestor_points but adds descendants instead
add_descendant_points <- function(data_in, descendants_in){
  descendants_in <- descendants_in %>% dplyr::select(-any_of(c("X1", "X2", "X3","X1_ancestor", "X2_ancestor", "X3_ancestor")))
  sub_data_in <- data_in %>% dplyr::select(sample, any_of(c("X1", "X2", "X3")))
  sub_data_in_descendant <- sub_data_in %>% rename_all(.funs = ~paste0(.,"_descendant"))
  
  descendants_points <- descendants_in %>% left_join( sub_data_in) %>% left_join(sub_data_in_descendant) 
  
  
  
  
  return(descendants_points)
}

#making ancestor and descendant tables based on metadata supplied; 
make_ancestor_table_meta <- function(data_in){
  #making a table where the ancestor type of each sample is noted
  starting_table <- data_in %>% transmute(sample_type = sample_type , ancestor_type = ancestor(sample_type)   , patient_id = patient_id, sample = sample, disease = disease , dataset = dataset )
  #calling the ancestor  sample_ancestor 
  ancestor_table  <- data_in %>% transmute(sample_type= sample_type, patient_id = patient_id, sample_ancestor = sample) 
  
  #joint left join on ancestor type and patient id so that each "ancestor_type" and "patient_id" in starting table  is matched with the corresponding value in ancestor table
  ancestor_table = left_join(starting_table, ancestor_table, by = c("ancestor_type" = "sample_type" ,"patient_id" =  "patient_id"))
  
  return(ancestor_table)
  
}

#functions same way as make ancestor table meta
make_descendant_table_meta <- function(data_in){
  #making a table where the descendant type of each sample is noted
  starting_table <- data_in %>% transmute(sample_type = sample_type , descendant_type = descendant(sample_type)   , patient_id = patient_id, sample = sample, disease = disease , dataset = dataset )
  #calling the descendant  sample_descendant 
  descendant_table  <- data_in %>% transmute(sample_type= sample_type, patient_id = patient_id, sample_descendant = sample) 
  
  #joint left join on descendant type and patient id so that each "descendant_type" and "patient_id" in starting table  is matched with the corresponding value in descendant table
  descendant_table = left_join(starting_table, descendant_table, by = c("descendant_type" = "sample_type" ,"patient_id" =  "patient_id"))
  
  return(descendant_table)
  
}



#data in is all_metadata from pipeline
get_target_site <- function(data_in){
 
  temp_target <- data_in %>% dplyr::select(patient_id, sample_type, biopsy_tissue, sample)  
  
  #getting latest known target tissue site for each sample
 
  temp_target <- temp_target %>% 
    mutate(sample_type = make_tissue_factor(sample_type)) %>%   #makes an ordered factor so we can use max and min on it 
    group_by(patient_id) %>% 
    mutate(target_sample_type = max(sample_type)) %>%   # this is adding  the destination; the max of the ordered sample type factor
    group_by(patient_id) %>% 
    dplyr::filter(sample_type  == max(target_sample_type)) %>%   #getting samples that are at their final known locatin 
    dplyr::rename( target_site = biopsy_tissue ) %>% 
    dplyr::select(patient_id,target_site) %>% distinct() #only keeping patient id and target site
  
  data_in <- left_join(data_in, temp_target, by = "patient_id")
  return(data_in)
}





#called proper because there was an old function that is now cone called plot_lines
#table in is something that has information about the metastasis locations, scatter in is the ggplot2 or plotly scatter we are adding the lines connecting samples from same patient with
plot_lines_proper <- function(table_in, scatter_in){

  
 
 remove_num <-  function(x) gsub('[0-9]+', '', x)  #function to replace numbers in a string
 

    
    if(three_dee_output){
      browser()
      #complicated plot the lines in 3d using plotly; there is no good way to do this so I had to do a complex pivoting option and use group2na to separate the things that should be separated by NA values
      lines_table  <- table_in %>%
        dplyr::select( any_of(c("sample", "sample_ancestor", "X1_ancestor",   "X2_ancestor", "X3_ancestor"  , "X1","X2","X3")) ) %>% 
        dplyr::rename(X1_origin = X1 ,   X2_origin = X2, X3_origin = X3) %>%  #adding suffix that will be helpful
        dplyr::filter(!is.na(sample_ancestor)) %>%   #getting rid of the samples that have no ancestors
        pivot_longer(cols = c("X1_origin", "X1_ancestor","X2_origin", "X2_ancestor","X3_origin", "X3_ancestor"), names_to = "coord_type", values_to = "coord") %>% #pivoting longer so that each element of each vector is on it's own row
        tidyr::separate(coord_type,sep = "_", into = c("dim","type")) %>%  #separating the dimension(ie X1,X2,X3) from whether it is an ancestor or the point in question
        pivot_wider(names_from = dim, values_from = coord) %>%  #pivoting wider but now ancestors and origins samples are on separate lines
        group2NA(groupNames = c("sample","sample_ancestor"), ordered = "type")  #separating by NA for plotting
   
  #adding the lines connecting primary to mets using plotly
    scatter_in <-  scatter_in %>% add_paths(data = lines_table , x = ~X1 , y = ~X2, z = ~X3,
                                            line = list(color = '#000000', width = 3, opacity = 0) , 
                                            name = "Same Patient", type = "scatter3d",mode = "lines", showlegend = TRUE)
    
    
  } else {  #adding the ancestor-descendant lines in 2d; much easier in ggplot2
 
    scatter_in <- scatter_in + geom_segment(data = table_in, mapping = aes(x = X1_ancestor,y = X2_ancestor,
                                                                           xend = X1 ,yend = X2, linetype = "")) + 
      labs(linetype = "Same Patient") +
      guides( color = guide_legend(order = 1),shape = guide_legend(order = 2), linetype = guide_legend(order = 3)) #setting order of the legend
  }
  

  return(scatter_in)
}



extract_anova_pvalue <- function(anova_in){
  summary(anova_in)[[1]][["Pr(>F)"]][1]
}





#indicages metastasis by putting a black ring around the point; probably better to used the filled vs hollow option (not this function); color selection indicates what we're coloring by in the data
plot_black_rings_2d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in) { #passing the parent environment so that we don't have to specify many optiond
  
  plot <- ggplot(data_in,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color = !!(border_selection) ,fill= !! (color_slxn)),shape = 21, stroke = 2) + 
    labs(title = env$title_string,  y = paste(env$axis_label,'2'), x= paste(env$axis_label,'1') ,caption = env$caption_string, color = make_title_case(quo_name(border_selection)) ,fill = make_title_case(quo_name(color_slxn)))+ 
    theme(aspect.ratio=1,text = element_text(size=font_size)) + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = toupper(palette)) + scale_color_manual(values= c("#000000FF","#FFFFFF00","#FF00FF"))
  return(plot)
}


#like plot_black_rings_3d but used hollow points to indicate metastasis
plot_hollow_rings_2d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in) { 

  plot <- ggplot(data_in,aes(x=X1,y=X2)) + 
    geom_point(size = symbol_size,aes(color = !!(color_slxn) ,shape = !! (border_selection)), stroke = 2) + 
    labs(title = env$title_string,  y = paste(env$axis_label,'2'), x= paste(env$axis_label,'1') ,caption = env$caption_string, shape = make_title_case(quo_name(border_selection)) ,color = make_title_case(quo_name(color_slxn)))+ 
    theme(aspect.ratio=1,text = element_text(size=font_size)) + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = toupper(palette)) + scale_color_manual(values = toupper(palette)) +scale_shape_manual(values = shape_scale)
  return(plot)
}




#plots the all points in same shape then uses plot_mets_black to add black rings behind the mets; not ideal solution, should use plot_hollow_rings_3d
plot_black_rings_3d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in){

  out <- plot_ly() %>%
    add_markers( data = data_in, x = ~X1,marker = list( sizemode = 'diameter'), y = ~X2, z = ~X3, color = color_slxn, colors = palette,opacity = .8,size = ~symbol_size ,sizes = c(symbol_size*symbol_size_factor,symbol_size*symbol_size_factor)) %>%
    layout(title = paste0(env$title_string,": ",make_title_case(quo_name(color_slxn))),scene = list(xaxis = list(title = paste(env$axis_label,'1')),
                                                                                                    yaxis = list(title = paste(env$axis_label,'2')),
                                                                                                    zaxis = list(title = paste(env$axis_label,'3'))))

  out<- plot_mets_black(data_in,out)
  return(out)
}



#plots hollow rings to indicate metastases, anod other shapes for other sample types in 3d
plot_hollow_rings_3d <- function(color_slxn,env = as.list.environment(parent.frame()),palette,data_in){
  
  used_shape_scale <- shape_scale_3d[names(shape_scale_3d) %in% dplyr::pull(data_in,(quo_name(border_selection)))]
  
  out <- plot_ly() %>% 
    add_markers( data = data_in, marker = list( sizemode = 'diameter'),x = ~X1, y = ~X2, z = ~X3, color = color_slxn, colors = palette,opacity = .8, symbol = border_selection ,symbols = used_shape_scale, size = ~symbol_size, sizes = c(symbol_size*symbol_size_factor,symbol_size*symbol_size_factor)) %>%
    layout(title = paste0(env$title_string,": ",make_title_case(quo_name(color_slxn))),scene = list(xaxis = list(title = paste(env$axis_label,'1')),
                                                                                                    yaxis = list(title = paste(env$axis_label,'2')),
                                                                                                    zaxis = list(title = paste(env$axis_label,'3'))))

  
  return(out)
}

#called "proper" because there was a previous function that did not plot it as nicely; adds black points on top of regular colored points in plotly 3d scatter plot; not ideal
#searched for a way to do it in the main plotly call but had to rely on this hack; uses plotly
plot_mets_black <- function(joined_dr, scatter_in){
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














