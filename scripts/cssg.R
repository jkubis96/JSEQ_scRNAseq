#Library for single cell data resolution improvement

library(Matrix)
library(stringr)

heterogenity_select <- function(cells_wide_df, marker_df, heterogenity_factor, p_val, max_genes, select_stat, min_occ = 5, qmin = .10) {
  
  set.seed(123)

  if (unique(grepl('[0-9]', colnames(cells_wide_df))) == TRUE) {
    
    cluster_group <-  sort(as.numeric(as.character(unique(colnames(cells_wide_df)))))
    
  } else {
    
    cluster_group <-  as.character(unique(colnames(cells_wide_df)))
    
  }
  
  marker_df <- marker_df[marker_df$p_val < p_val,]
  

  for (cluster in cluster_group) {
    
    cat(paste('\n\n Cluster ', cluster, '- searching marker genes... \n\n' ))
    
    gen_cor <- unique(marker_df$gene[marker_df$cluster %in% cluster])
    
    
    tmp2 <- cells_wide_df[, colnames(cells_wide_df) %in% cluster, drop = FALSE]
    
    tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), , drop = FALSE]
    
    
    #expression level 
    
    rmean <- rowMeans(tmp2 > 0)
    rmean <- rmean[rmean > quantile(rmean, qmin)]
    
    if (length(rmean) > 10) {
      
      tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(names(rmean)), , drop = FALSE]
      
    }
    
    
    tmp2[tmp2 > 0] <- 1L
    tmp2[tmp2 == 0] <- 0L

    perc <- (rowSums(tmp2 == 1) / ncol(tmp2)) * 100
    
    perc2 <- perc[perc <= heterogenity_factor & perc >= min_occ]
    

    
    if (exists('heterogenity_markers_df') == FALSE) {
      heterogenity_markers_df <- data.frame(gen = names(perc2), per_obj = perc2)
      heterogenity_markers_df$cluster <- cluster
    } else {
      tmp_het <- data.frame(gen = names(perc2), per_obj = perc2)
      tmp_het$cluster <- cluster
      heterogenity_markers_df <- rbind(heterogenity_markers_df, tmp_het)
    }
    
    marker_df_CSSG <- marker_df[marker_df$cluster %in% cluster, ]
    
    if (length(perc2) > 10) {
      
      marker_df_CSSG <- marker_df_CSSG[marker_df_CSSG$gene %in% names(perc2), ]
      
    } else {
      
      marker_df_CSSG <- marker_df_CSSG[marker_df_CSSG$gene %in% names(perc), ]
      
    }
  
  
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
      
      set.seed(123)
  
      add_avg_val <- function(x, tmp2) {
        
        av <- mean(tmp2[x,])
        return(av)
        
      }
      
      
      check_memory_par <- function(par_object, per_mem = .8) {
        os_type <- Sys.info()["sysname"]
        

        
        if (os_type == "Windows") {

          current_memory <- memory.size()
          max_memory <- memory.limit()
          avaiable <- max_memory - current_memory
          
          object_size_bytes <- 0
          
          for (obj in par_object) {
            
            object_size_bytes <- object_size_bytes + object.size(par_object)
            
            
          }
          
          
          object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2
          
          av_cpu <- (avaiable*per_mem)/object_size_mb
          
          CPU <- detectCores()
          
          if (av_cpu >= CPU) {
            av_cpu <- CPU - 2
          } else if (av_cpu < 1) {
            av_cpu <- 1
            
          }
          
          
        } else if (os_type == "Linux") {

          available_memory_bytes <- as.numeric(system("awk '/MemAvailable:/ {print $2 * 1024}' /proc/meminfo", intern = TRUE))
          avaiable <- available_memory_bytes / (1024^2)
          
          
          object_size_bytes <- 0
          
          for (obj in par_object) {
            
            object_size_bytes <- object_size_bytes + object.size(par_object)
            
            
          }
          
          
          object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2
          
          av_cpu <- (avaiable*per_mem)/object_size_mb
          
          CPU <- detectCores()
          
          if (av_cpu >= CPU) {
            av_cpu <- CPU - 2
          } else if (av_cpu < 1) {
            av_cpu <- 1
            
          }
          
          
        } else if (os_type == "Darwin") {

          free_pages <- as.numeric(system("vm_stat | grep 'Pages free:' | awk '{print $3}' | sed 's/[^0-9]//g'", intern = TRUE))
          
          page_size <- as.numeric(system("sysctl -n hw.pagesize", intern = TRUE))
          
          available_memory_bytes <- free_pages * page_size
          
          avaiable <- available_memory_bytes / (1024^2)
          
          object_size_bytes <- 0
          
          for (obj in par_object) {
            
            object_size_bytes <- object_size_bytes + object.size(par_object)
            
            
          }
          
          
          object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2
          
          av_cpu <- ((avaiable*per_mem)/object_size_mb) - 1
          
          CPU <- detectCores()
          
          
          if (av_cpu >= CPU) {
            av_cpu <- CPU - 2
          } else if (av_cpu <= 1) {
            av_cpu <- 1
          }
          
          
        } else {
          av_cpu <- detectCores() - 2
          print("Unknown or unsupported operating system.")
        }
        
        print(paste("Number of CPU cores allowed for parallel processing:", av_cpu))
        
        return(av_cpu)
      }
      
      
  
      cat('\n\n The CSSG start \n it can last several minutes depending on the number of clusters and set parameters \n\n\n')
      
      
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
 

      
        
        while (TRUE) {
          
          ls <- ls(all.names = TRUE)
          
          CPU <- check_memory_par(par_object = c('res_df', 'tmp3', 'tmp2', 'add_avg_val'))
          cl <- makeCluster(CPU)
          registerDoParallel(cl)
          registerDoSNOW(cl)
          
          
          pb <- txtProgressBar(max = nrow(res_df), style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
          
          
          res_df_multi <- foreach(i =  1:nrow(res_df), .options.snow = opts, .combine=rbind, .inorder=FALSE, .packages = c('Matrix', 'stringr'), .export = c('res_df', 'tmp3', 'tmp2', 'add_avg_val'), .noexport = ls) %dopar% {
         
            
            exclude_rows <- strsplit(rownames(res_df)[i], split = ' ')[[1]]
            add_tmp  <- as.matrix(tmp3[!rownames(tmp3) %in% exclude_rows, ,drop = FALSE])
            
            if (nrow(add_tmp) != 0) {
              
            row_vector =  paste(rownames(add_tmp), rownames(res_df)[rep(i, nrow(add_tmp))])
            
            row_vector <- strsplit(row_vector, split = ' ')
            row_vector <- lapply(row_vector, sort)
            row_vector <- lapply(row_vector, function(x) paste(x, collapse = ' '))
            

            # Perform the matrix addition
            
            results_tmp <- add_tmp + as.matrix(res_df[rep(i, nrow(add_tmp)), , drop = FALSE])
            rownames(results_tmp) <- row_vector

            } else {
              
              results_tmp <- as.matrix(res_df)
            }
          
            # top1
            #
            ncol_results_tmp <- ncol(results_tmp)
            perc0 <- rowSums(results_tmp == 0) / ncol_results_tmp
            perc1 <- rowSums(results_tmp == 1) / ncol_results_tmp
            avg_exp <- unlist(lapply(strsplit(rownames(results_tmp), split = ' '), function(x) add_avg_val(x, tmp2)))
            multi = 1 - (perc0 + perc1)

            sorting_df <- data.frame(perc1, avg_exp)
            
            results_tmp <- cbind(results_tmp, perc0, perc1, avg_exp, multi)
           
            
            rm(perc0, perc1, avg_exp, multi)
            
            
            results_tmp  <- results_tmp[order(sorting_df$perc1, sorting_df$avg_exp, decreasing = TRUE),]
           
          
            if (nrow(add_tmp) != 0) {
              
              results_tmp <- results_tmp[1, , drop = FALSE]
              
            }
            
            gc()
            
            return(results_tmp)
            
          }
          
          
          
          res_df_multi <- res_df_multi[!duplicated(rownames(res_df_multi)), , drop = FALSE]
          final_df <- as.data.frame(res_df_multi[,!colnames(res_df_multi) %in% cluster, drop = FALSE])

          
          res_df <- as.matrix(res_df_multi[,colnames(res_df_multi) %in% cluster, drop = FALSE])

          if (max_combine < nrow(res_df)) {
            
            res_df <- res_df[1:max_combine, ,drop = FALSE]
            
          }
          
          
          
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
          
          close(pb)
          stopCluster(cl) 
          gc()
          
          
        }
        
        
      }  
      
      
      
      colnames(complete_df) <- c('loss_val', 'perc_1', 'avg_exp', 'perc_multi', 'hf', 'adj_hf', 'cluster')
      cat('\n\n The CSSG finish \n')
      
      return(complete_df)

      
}

###########################################################################################################################################################







