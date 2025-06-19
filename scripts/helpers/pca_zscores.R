# Get PCA Z-scores Helper
# Input:
# Data frame in the following wide format:
# |subjectid|signal_names|signal_components|X1|X2|....|X101|
# 
# Xs are the percent of gait cycle or "item" if using v3dR
#
# Result is a 2 item list: a zscores dataframe and another list
# of results from each PCAs for each metric.
#
# "zscores" dataframe:
# |subjectid|signal_names|signal_components|PCnum|Zscore|
#
# Each pca_list element is named by signal and component
#

library(tibble)

not_any_na <- function(x) all(!is.na(x))

pca_zscores <- function(df, scale, rank) {
  
  zscores_out = NULL
  propvar_out = NULL
  cumvar_out = NULL
  pca_list = list()
  
  # create dim names vector length of rank (number of pcs returned)
  PCnum <- NULL
  for (i in 1:rank) {
    PCnum[i] <- paste0("PC",i)
  }
  

  u_sn <- unique(df$signal_names)
  u_sn <- levels(u_sn)[as.numeric(u_sn)]
  
  for (i in 1:length(u_sn))
  {
    u_sc <- unique(df$signal_components)
    u_sc <- levels(u_sc)[as.numeric(u_sc)]
    
    
    for (j in 1:length(u_sc))
    {
      tmp_df <- df %>%
        filter(signal_names == u_sn[i]) %>%
        filter(signal_components == u_sc[j]) %>%
        #select(-c(session, signal_names, signal_components))
        dplyr::select(-c(signal_names, signal_components))
        
        
      
      if (dim(tmp_df)[1] == 0) {
        tmp_df <- NULL
        break()
      }
      
      # prcomp defaults center = TRUE and scale = FALSE
      df_PCA <- prcomp(tmp_df %>% dplyr::select(where(is.numeric)) %>% dplyr::select(where(not_any_na)), center = TRUE, scale. = scale, rank. = rank)

      # Add Variance, Proportion of Variance and Cumulative Variance
      eigs <- "eigs"
      propvar <- "propvar"
      cumvar <- "cumvar"
      
      df_PCA[[eigs]] <- df_PCA$sdev^2
      df_PCA[[propvar]] <- df_PCA$eigs/sum(df_PCA$eigs)
      df_PCA[[cumvar]] <- cumsum(df_PCA$eigs)/sum(df_PCA$eigs)
      

      # Cut Variance, Proportion of Variance and Cumulative Variance
      # to length rank
      df_PCA[[eigs]] <- as.matrix(df_PCA$eigs[1:rank])
      dimnames(df_PCA$eigs)[[1]] <- PCnum
      dimnames(df_PCA$eigs)[[2]] <- 'variance'
      

      df_PCA[[propvar]] <- as.matrix(df_PCA$propvar[1:rank])
      dimnames(df_PCA$propvar)[[1]] <- PCnum
      dimnames(df_PCA$propvar)[[2]] <- 'variance' 
      

      df_PCA[[cumvar]] <- as.matrix(df_PCA$cumvar[1:rank])
      dimnames(df_PCA$cumvar)[[1]] <- PCnum
      dimnames(df_PCA$cumvar)[[2]] <- 'variance' 
      

      # Create seperate df of just Proportion of Variance
      propvar <- as.data.frame(df_PCA$propvar) %>%
        rownames_to_column(var = "PC") %>%
        mutate(signal_components = u_sc[j],
               signal_names =  u_sn[i])
      
      # Create seperate df of just Cumulative Variance
      cumvar <- as.data.frame(df_PCA$cumvar) %>%
        rownames_to_column(var = "PC") %>%
        mutate(signal_components = u_sc[j],
               signal_names =  u_sn[i])
      
      # Create seperate df of just zscores
      zscores <- as.data.frame(df_PCA$x) %>%
        #rownames_to_column(var = "subject") %>%
        mutate(signal_names =  u_sn[i], signal_components = u_sc[j], .before = everything()) #add back in signal_components and signal_names before PCs
      
      # fold back in all data that was provided to the prcomp.
      zscores <- cbind(tmp_df %>% dplyr::select(!where(is.numeric)), zscores)
      
      # bind dfs loop for return
      zscores_out = rbind(zscores_out, zscores)
      propvar_out = rbind(propvar_out, propvar)
      cumvar_out = rbind(cumvar_out, cumvar)
      
      
      # need to dynamically name list.
      list_name = paste(u_sn[i], u_sc[j], sep="_")
      pca_list[[list_name]] <- df_PCA
      
    }
    
  }
  
  out_object = list("zscores" = zscores_out,
                    "propvar" = propvar_out,
                    "cumvar" = cumvar_out,
                    "pca_list" = pca_list)
  
  return(out_object)
}
