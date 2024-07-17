#Library for single cell data resolution improvement
library(Matrix)


heterogenity_select <- function(cells_wide_df, marker_df, heterogenity_factor, p_val, max_genes, select_stat, min_occ = 5) {
  

  if (unique(grepl('[0-9]', colnames(cells_wide_df))) == TRUE) {
    
    cluster_group <-  sort(as.numeric(as.character(unique(colnames(cells_wide_df)))))
    
  } else {
    
    cluster_group <-  as.character(unique(colnames(cells_wide_df)))
    
  }
  
  marker_df <- marker_df[marker_df$p_val < p_val,]
  

  for (cluster in cluster_group) {
    
    cat(paste('\n\n Cluster ', cluster, '- searching marker genes... \n\n' ))

    tmp2 <- as.data.frame(cells_wide_df[, colnames(cells_wide_df) %in% cluster])
    tmp2[tmp2 > 0] <- 1L
    tmp2[tmp2 == 0] <- 0L
    
    
    gen_cor <- unique(marker_df$gene[marker_df$cluster %in% cluster])
    tmp3 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), ]
    rm(tmp2)
    
    
    gen <- c()
    per_obj <- c()
    
    for (i in 1:length(rownames(tmp3))) {
      gen <- c(gen, rownames(tmp3)[i])
      per_obj <- c(per_obj, sum(tmp3[i,])/length(tmp3[i,])*100)
      
    } 
    
    rm(tmp3)
    df <- data.frame(as.character(gen), as.numeric(per_obj))
    colnames(df) <- c('gen', 'per_obj')
    
    if (exists('heterogenity_markers_df') == FALSE) {
      heterogenity_markers_df <- df
      heterogenity_markers_df$cluster <- cluster
    } else {
      tmp_het <- df
      tmp_het$cluster <- cluster
      heterogenity_markers_df <- rbind(heterogenity_markers_df, tmp_het)
     }
    

    genes_CSSG <- as.vector(df$gen[df$per_obj <= heterogenity_factor & df$per_obj >= min_occ])
    
        
    if (length(genes_CSSG) == 0) {
      
      genes_CSSG <- as.vector(df$gen[df$per_obj >= min_occ])
      
    } 
    
  
  
  marker_df_CSSG <- marker_df[marker_df$cluster %in% cluster, ]
  marker_df_CSSG <- marker_df_CSSG[marker_df_CSSG$gene %in% genes_CSSG, ]
  
  
    if (exists('marker_df_tmp') == FALSE) {
      if (select_stat %in% c('avg_logFC', 'avg_log2FC')) {
        marker_df_tmp <- marker_df_CSSG[order(marker_df_CSSG$avg_logFC, decreasing =  TRUE), ]
        marker_df_tmp <- marker_df_tmp[1:as.numeric(max_genes), ]
        marker_df_tmp <- marker_df_tmp[!is.na(marker_df_tmp$cluster),]
      } else {
        marker_df_tmp <- marker_df_CSSG[order(marker_df_CSSG$p_val, decreasing =  FALSE), ]
        marker_df_tmp <- marker_df_tmp[1:as.numeric(max_genes), ]
        marker_df_tmp <- marker_df_tmp[!is.na(marker_df_tmp$cluster),]
      }
      
    } else {
      if (select_stat %in% c('avg_logFC', 'avg_log2FC')) {
        marker_df_tmp1 <- marker_df_CSSG[order(marker_df_CSSG$avg_logFC, decreasing =  TRUE), ]
        marker_df_tmp1 <- marker_df_tmp1[1:as.numeric(max_genes), ]
        marker_df_tmp1 <- marker_df_tmp1[!is.na(marker_df_tmp1$cluster),]
      } else {
        marker_df_tmp1 <- marker_df_CSSG[order(marker_df_CSSG$p_val, decreasing =  FALSE), ]
        marker_df_tmp1 <- marker_df_tmp1[1:as.numeric(max_genes), ]
        marker_df_tmp1 <- marker_df_tmp1[!is.na(marker_df_tmp1$cluster),]
      }
      
      marker_df_tmp <- rbind(marker_df_tmp, marker_df_tmp1)
      
      
      }
    
  
  }
  
  return(list('marker_df' = marker_df_tmp, 'heterogenity_markers_df' = heterogenity_markers_df))
  
  
}

CSSG_markers <- function(cells_wide_df, markers_df, max_combine, loss_val) {
  
      add_avg_val <- function(x, tmp2) {
        
        av <- mean(tmp2[x,])
        return(av)
        
      }
      
  
      cat('\n\n The CSSG start \n it can last several minutes depending on the number of clusters and set parameters \n\n\n')
      
      CPU <- detectCores() - 1
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      
      if (unique(grepl('[0-9]', colnames(cells_wide_df))) == TRUE) {
        
        cluster_group <-  sort(as.numeric(as.character(unique(colnames(cells_wide_df)))))
        
      } else {
        
        cluster_group <-  as.character(unique(colnames(cells_wide_df)))
        
      }
      
      complete_df <- data.frame()

      for (cluster in cluster_group) {
        
        cat(paste('\n\n Cluster ', cluster, '- searching heterogeneity marker genes... \n\n' ))
        
        
        # start data storing 
        
        approved_df <- data.frame()
        
        
        # initisl val 
        tmp_recursive <- NaN
        
        
        gen_cor <- unique(markers_df$gene[markers_df$cluster %in% cluster])
        
        
        tmp2 <- cells_wide_df[, colnames(cells_wide_df) %in% cluster, drop = FALSE]
        tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), , drop = FALSE]

        
        
        tmp3 <- tmp2
        tmp3[tmp3 > 0] <- 1L
        tmp3[tmp3 == 0] <- 0L
        
        tmp2 <- data.frame(row.names = rownames(tmp2), avg = rowMeans(tmp2))
       
        
        # First loop for duble combination
        
        res_df <- tmp3
        
        ncol_res_df <- ncol(res_df)
        perc0 <- rowSums(res_df == 0) / ncol_res_df
        perc1 <- rowSums(res_df == 1) / ncol_res_df
        avg_exp <- unlist(lapply(strsplit(rownames(res_df), split = ' '), function(x) add_avg_val(x, tmp2)))
        multi = 1 - (perc0 + perc1)
        
        # Combine perc0 and perc1 into a sparse matrix (optional)
        last_df <- as.data.frame(cbind(perc0, perc1, avg_exp, multi))
        
        # Filter based on perc1 quantile

        last_df <- last_df[last_df$perc1 > quantile(last_df$perc1, 0.6), ,drop = FALSE]
        
        
        # first results
        
        
        approved_df <- rbind(approved_df, last_df[last_df$perc0 <= loss_val, ,drop = FALSE])
        
        if (nrow(approved_df) > 0) {
          
          approved_df$het <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1)))
          approved_df$het_adj <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1))) - (approved_df$perc0*2)
          approved_df$cluster <- cluster
          
          complete_df <- rbind(approved_df, complete_df)
          
        }
 

        # res_df <- res_df[rownames(res_df) %in% rownames(last_df), ,drop = FALSE]


     
        while (TRUE) {
          
          
          pb <- txtProgressBar(max = nrow(res_df), style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
          
          
          res_df_multi <- foreach(i =  1:nrow(res_df), .options.snow = opts, .combine=rbind) %dopar% {
            
            library(Matrix)
            library(stringr)
            
            
            exclude_rows <- strsplit(rownames(res_df)[i], split = ' ')[[1]]
            add_tmp  <- tmp3[!rownames(tmp3) %in% exclude_rows, ,drop = FALSE]
            
            if (nrow(add_tmp) != 0) {
              
            row_vector =  paste(rownames(add_tmp), rownames(res_df)[rep(i, nrow(add_tmp))])
            
            row_vector <- strsplit(row_vector, split = ' ')
            row_vector <- lapply(row_vector, sort)
            row_vector <- lapply(row_vector, function(x) paste(x, collapse = ' '))
            

            
            # Perform the matrix addition
            results_tmp <- add_tmp + res_df[rep(i, nrow(add_tmp)), , drop = FALSE] 
            rownames(results_tmp) <- row_vector

            } else {
              
              results_tmp <- res_df
            }
          
            # top1
            #
            ncol_results_tmp <- ncol(results_tmp)
            perc0 <- rowSums(results_tmp == 0) / ncol_results_tmp
            perc1 <- rowSums(results_tmp == 1) / ncol_results_tmp
            avg_exp <- unlist(lapply(strsplit(rownames(results_tmp), split = ' '), function(x) add_avg_val(x, tmp2)))
            multi = 1 - (perc0 + perc1)
            
            sorting_df <- data.frame(perc1, avg_exp)
            
            perc0 <- Matrix(perc0, ncol = 1, sparse = TRUE)
            perc1 <- Matrix(perc1, ncol = 1, sparse = TRUE)
            avg_exp <- Matrix(avg_exp, ncol = 1, sparse = TRUE)
            multi <- Matrix(multi, ncol = 1, sparse = TRUE)

           

            results_tmp <- cbind(results_tmp, perc0, perc1, avg_exp, multi)
            colnames(results_tmp)[(ncol(results_tmp)-3):ncol(results_tmp)] <- c("perc0", "perc1", "avg_exp", "multi")
            
            
            rm(perc0, perc1, avg_exp, multi)
            
            
            # Get order of rows based on sorting criteria
            sorted_indices <- order(sorting_df$perc1, sorting_df$avg_exp, decreasing = TRUE)
            
            # Subset the sparse matrix based on the sorted order
            results_tmp <- results_tmp[sorted_indices, , drop = FALSE]
            
            # Further subset to keep only the first row
            
          
            if (nrow(add_tmp) != 0) {
              
              results_tmp[1, , drop = FALSE]
              
            }

            return(results_tmp)
            
         
          }
          
          
          
          res_df_multi <- res_df_multi[!duplicated(rownames(res_df_multi)), , drop = FALSE]
          final_df <- as.data.frame(res_df_multi[,!colnames(res_df_multi) %in% cluster, drop = FALSE])

          
          res_df <- as.matrix(res_df_multi[,colnames(res_df_multi) %in% cluster, drop = FALSE])

          if (max_combine < nrow(res_df)) {
            
            res_df <- res_df[1:max_combine, ,drop = FALSE]
            
          }
          
          res_df <- Matrix(res_df, sparse = TRUE)
          
          
          
          colnames(res_df) <- rep(as.character(cluster), ncol(res_df))
          rm(res_df_multi)
        

            
          if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] == 0) {
            
            approved_df <- rbind(approved_df, final_df[final_df$perc0 == 0 & final_df$multi < 0.20, , drop = FALSE])
            to_exclude <- rownames(final_df)[final_df$perc0 == 0 & final_df$multi < 0.20]
            res_df <- res_df[!rownames(res_df) %in% to_exclude, ,drop = FALSE]
            
            
          } else if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] < as.numeric(loss_val)) {
              
            approved_df <- rbind(approved_df, final_df[final_df$perc0 <= as.numeric(loss_val), ,drop = FALSE])
              
              
          } 
          
          
          
          
          
          if (nrow(approved_df) >= max_combine) {
            
            approved_df$het <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1)))
            approved_df$het_adj <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1))) - (approved_df$perc0*2)
            approved_df$cluster <- cluster
            
            complete_df <- rbind(approved_df, complete_df)
            
            break
            
            
          } else if (nrow(res_df) == 0) {
          
          approved_df <- rbind(approved_df, final_df[final_df$perc0 <= quantile(final_df$perc0, 0.25), ,drop = FALSE])
          approved_df$het <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1)))
          approved_df$het_adj <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1))) - (approved_df$perc0*2)
          approved_df$cluster <- cluster
          
          complete_df <- rbind(approved_df, complete_df)
          
          break
          
          
        }  else if (is.data.frame(tmp_recursive)) {
          
          
          if ((min(tmp_recursive$perc0) == min(final_df$perc0)) & (min(tmp_recursive&multi) <= min(final_df$multi))) {
          
          approved_df <- rbind(approved_df, tmp_recursive[tmp_recursive$perc0 <= quantile(tmp_recursive$perc0, 0.25), ,drop = FALSE])
          approved_df$het <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1)))
          approved_df$het_adj <- (1- ((approved_df$perc1 + ((1-(approved_df$perc1 + approved_df$perc0))*1.25))/(str_count(string = rownames(approved_df), pattern = ' ') + 1))) - (approved_df$perc0*2)
          approved_df$cluster <- cluster
          
          complete_df <- rbind(approved_df, complete_df)
          
          break
          
          } else {
            
            tmp_recursive <- final_df
            
            
          }
          
          
        } else {
            
            tmp_recursive <- final_df
            
            
        }
          
          
        }
        
        
      }  
      
      
      
      
      close(pb)
      stopCluster(cl)  
      
      colnames(complete_df) <- c('loss_val', 'perc_1', 'avg_exp', 'perc_multi', 'hf', 'adj_hf', 'cluster')
      cat('\n\n The CSSG finish \n')
      
      return(complete_df)

      
}

###########################################################################################################################################################







