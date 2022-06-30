#Library for single cell data resolution improvement

heterogenity_select <- function(cells_wide_df, marker_df, heterogenity_factor, p_val, max_genes, select_stat) {
    
  if (unique(grepl('[0-9]', colnames(cells_wide_df))) == TRUE) {
    
    cluster_group <-  sort(as.numeric(as.character(unique(colnames(cells_wide_df)))))
    
  } else {
    
    cluster_group <-  as.character(unique(colnames(cells_wide_df)))
    
  }
  
  marker_df <- marker_df[marker_df$p_val < p_val,]
  

  for (cluster in cluster_group) {
    
    

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
    
    genes_CSSG <- as.vector(df$gen[df$per_obj <= heterogenity_factor])
    
  
  
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

CSSG_markers <- function(cells_wide_df, markers_df, max_combine, loss_pval) {
  
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
      

      for (cluster in cluster_group) {
        
        cat(paste('\n\n Cluster ', cluster, '- searching heterogeneity marker genes... \n\n' ))
        
        tmp2 <- as.data.frame(cells_wide_df[, colnames(cells_wide_df) %in% cluster])
        tmp2[tmp2 > 0] <- 1L
        tmp2[tmp2 == 0] <- 0L
        
        gen_cor <- unique(markers_df$gene[markers_df$cluster %in% cluster])
        tmp3 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), ]
        rm(tmp2)
        
        
        # First loop for duble combination
        
        res_df <- as.matrix(tmp3)
        
        perc0 <- apply(res_df == 0, 1, sum)
        perc0 <- perc0/length(res_df[1,])
        perc1 <- apply(res_df == 1, 1, sum)
        perc1 <- perc1/length(res_df[1,])
        last_df <- data.frame(perc0, perc1)
        rownames(last_df) <- rownames(res_df)
        up_tr_q <- quantile(last_df$perc1, 0.75)
        last_df <- last_df[last_df$perc1 >= up_tr_q,]

      

        res_df <- res_df[rownames(res_df) %in% rownames(last_df),]

        
       
        low = 1
        up = 0
        test = TRUE
        while (test == TRUE) {
          
          iterations <- (length(rownames(res_df)))
          pb <- txtProgressBar(max = iterations, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
          
          
          
          res_df_multi <- foreach(i =  1:nrow(res_df), .options.snow = opts, .combine=rbind) %dopar% {
            genes_df <- tmp3[!rownames(tmp3) %in% strsplit(rownames(res_df)[i], split = ' ')[[1]],]
            rep_df <- res_df[rep(i, nrow(genes_df)),]
            rownames_tmp <- rep(rownames(rep_df)[1], nrow(genes_df))
            rownames(rep_df) <- make.unique(rownames_tmp, sep = '>') 
            results_tmp <- as.matrix(rep_df) + as.matrix(genes_df)
            rownames(results_tmp) <- paste(gsub('>.*', '', rownames(rep_df)), rownames(genes_df))
            return(results_tmp)
            
            
          }
          
          perc0 <- apply(res_df_multi == 0, 1, sum)
          perc0 <- perc0/length(res_df_multi[1,])
          perc1 <- apply(res_df_multi == 1, 1, sum)
          perc1 <- perc1/length(res_df_multi[1,])
            
          final_df <- data.frame(perc0, perc1)
          rownames(final_df) <- rownames(res_df_multi)
            
          #Dedup  
          
          row_vector <- strsplit(rownames(final_df), split = ' ')
          row_vector <- lapply(row_vector, sort)
          row_vector <- lapply(row_vector, toString)
          
          final_df <- as.data.frame(final_df)
          final_df$duplicates <- unlist(row_vector)
        
            
            
          final_df <- final_df[!duplicated(final_df$duplicates),]
          final_df <- final_df[,!colnames(final_df) %in% 'duplicates']
          up_tr_q <- quantile(final_df$perc1, 0.75)
          final_df_up <- final_df[final_df$perc1 >= up_tr_q,]
          final_df_up <- final_df_up[order(final_df_up$perc1,decreasing = TRUE),]
		      final_df_up <- final_df_up[1:as.numeric(max_combine),]
          final_df_up <- drop_na(final_df_up)
		  
		      down_tr_q <- quantile(final_df$perc0, 0.25) 
          final_df_down <- final_df[final_df$perc0 <= down_tr_q,]
          final_df_down <- final_df_down[order(final_df_down$perc0,decreasing = FALSE),]
          final_df_down <- final_df_down[1:as.numeric(max_combine),]
          final_df_down <- drop_na(final_df_down)
         
		      duplicated_rows <- rownames(final_df_up) %in% rownames(final_df_down)
		      final_df <- rbind(final_df_up[!duplicated_rows,], final_df_down)
		  
          res_df <- res_df_multi
          res_df <- res_df[rownames(res_df) %in% rownames(final_df),]
      
        
          
          if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] < as.numeric(loss_pval)) {
            last_df <- rbind(last_df, final_df)
            last_df$het <- (1- ((last_df$perc1 + ((1-(last_df$perc1 + last_df$perc0))*1.25))/(str_count(string = rownames(last_df), pattern = ' ') + 1)))
			      last_df$het_adj <- (1- ((last_df$perc1 + ((1-(last_df$perc1 + last_df$perc0))*1.25))/(str_count(string = rownames(last_df), pattern = ' ') + 1))) - (last_df$perc0*2)
            last_df <- last_df[, !colnames(last_df) %in% 'perc1']
            colnames(last_df) <- c('loss_pval','hf', 'adj_hf')
            test = FALSE
          } else if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] >= low & final_df$perc1[order(final_df$perc1, decreasing = TRUE)][1] <= up) {
            last_df$het <- (1- ((last_df$perc1 + ((1-(last_df$perc1 + last_df$perc0))*1.25))/(str_count(string = rownames(last_df), pattern = ' ') + 1)))
			      last_df$het_adj <- (1- ((last_df$perc1 + ((1-(last_df$perc1 + last_df$perc0))*1.25))/(str_count(string = rownames(last_df), pattern = ' ') + 1))) - (last_df$perc0*2)
            last_df <- last_df[, !colnames(last_df) %in% 'perc1']
            colnames(last_df) <- c('loss_pval','hf', 'adj_hf')
            test = FALSE
          } else {
            last_df <- rbind(last_df, final_df)
            
          }
          
          low <- final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1]
          up <- final_df$perc1[order(final_df$perc1, decreasing = TRUE)][1]
      
          
          
          
          
          
        }
        
          if (exists('complete_df') == FALSE) {
            complete_df <- last_df[last_df$`loss_pval` <= quantile(last_df$`loss_pval`, 0.25),]
            complete_df$cluster <- cluster
          } else {
            last_df <- last_df[last_df$`loss_pval` <= quantile(last_df$`loss_pval`, 0.25),]
            last_df$cluster <- cluster
            complete_df <- rbind(complete_df, last_df)
          }
        
      }  
      
      
      
      
      close(pb)
      stopCluster(cl)  
      
      cat('\n\n The CSSG finish \n')
      
      return(complete_df)

      
}

###########################################################################################################################################################







