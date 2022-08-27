#Cell name change for future
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cluster_nameing <- function(matrix_a, markers) {
  
  
  matrix_a <- as.data.frame(matrix_a)
  colname <- colnames(matrix_a)
  rownames(matrix_a) <- make.unique(toupper(rownames(matrix_a)), sep = '')
  
   if (length(markers_subclass) != 0) {
  
  list.markers <- c()
  for (l.marker in markers) {
    for (marker in l.marker) {
      TF <- (grepl('+', marker, fixed = T))
      if (TF == TRUE)  {
        list.markers <- c(list.markers, textclean::mgsub(marker, c('+'), c('')))
      }
    } 
  } 
  
  
  index = 0
  for (i in colnames(matrix_a)) {
    index = index +1
    rename_df <- as.data.frame(matrix_a[toupper(rownames(matrix_a)) %in% toupper(list.markers),index ,drop = F])
    rename_df <- as.data.frame(rename_df[order(rename_df, decreasing = T), ,drop = F])
    sum <- sum(rename_df)
    if (sum > 0) {
      colnames(matrix_a)[index] <- rownames(rename_df)[1]
    } else colnames(matrix_a)[index] <- 'Unknow'
  }  
  
  
  col <- 0
  for (c.marker in markers) {
    col <- col + 1
    for (marker in c.marker) {
      cell_n <- 0
      for (cell in colnames(matrix_a)) {
        cell_n <- cell_n + 1 
        if (cell %in% textclean::mgsub(marker, c('+'), c(''))) {
          colnames(matrix_a)[cell_n] <- colnames(markers[col])
        }
      }
    }
  }
  
  } else {
  
  colnames(matrix_a) <- paste0('Cluster_', colnames(matrix_a))
  
  }
  
  
  
  average_expression <<- matrix_a
  rm(matrix_a)
  
}

cssg_naming <- function(seurat_project, CSSG_df) {
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  cell_names <- c()
  for (cluster in sort(as.numeric(as.character(unique(seurat_project@active.ident))))) {
    sub_cssg <- CSSG_df[CSSG_df$cluster %in% cluster, ]
    sub_cssg <- sub_cssg[sub_cssg$`loss_pval` == min(sub_cssg$`loss_pval`), ]
    sub_cssg <- sub_cssg[order(sub_cssg$`adj_hf`, decreasing = TRUE),]
    tmp <- subset(seurat_project, features =  gsub(' ', '', strsplit(rownames(sub_cssg)[1], split = ' ')[[1]]))
    tmp <- as.matrix(GetAssayData(tmp, slot = 'data'))
    colnames(tmp) <- seurat_project@active.ident
    for (i in 1:length(colnames(tmp))) {
      if (as.numeric(cluster) %in% colnames(tmp)[i] & colSums(tmp)[i] > 0) {
        tmp <- tmp[order(tmp[,i], decreasing = T), ,drop = F]
        cell_names[i] <- firstup(tolower(rownames(tmp)[1]))
      } else if (as.numeric(cluster) %in% colnames(tmp)[i] & colSums(tmp)[i] == 0) {
        cell_names[i] <- 'Bad!'
      }
    }
  }
  
  return(cell_names)
  
}


subcluster_naming <- function(average_expression, markers_subclass, cell_markers, top_cell_markers) {
  
  if (length(markers_subclass) != 0) {
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    
    cell_names.1 <- c()
    for (i in 1:length(colnames(average_expression))) {
      rename_df <- average_expression[rownames(average_expression) %in% markers_subclass,]
      rename_df <- as.data.frame(rename_df[order(rename_df[,i], decreasing = T), ,drop = F])
      if (sum(rename_df[,i]) > 0) {
        cell_names.1[i] <- toupper(rownames(rename_df[1,]))
      } else {
        cell_names.1[i] <- 'BAD!'
      } 
      
      if (colnames(rename_df)[i] %in% 'Unknow') {
      cell_names.1[i] <- ''
      }
    }
    

    
    cluster <- 0
    cell_names.2 <- c()
    for (col in 1:length(colnames(average_expression))) {
      tmp_names <- cell_markers[cell_markers$cluster %in% cluster,]
      tmp_names <- tmp_names[tmp_names$gen %in% top_cell_markers$gene[top_cell_markers$cluster %in% cluster],]
      tmp_names <- as.data.frame(tmp_names[order(tmp_names$per_obj, decreasing = T), ])
      cluster <- cluster + 1
      if (!toupper(cell_names.1[col]) %in% toupper(tmp_names$gen[1])) {
        cell_names.2[col] <- firstup(tolower(tmp_names$gen[1]))
      } else if (toupper(cell_names.1[col]) %in% toupper(tmp_names$gen[1])) {
        cell_names.2[col] <- firstup(tolower(tmp_names$gen[2]))
      }
    }
    
     sub_names <- paste(cell_names.1, cell_names.2)
    
  } else if (length(markers_subclass) == 0) {
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    cluster <- 0
    cell_names.1 <- c()
    cell_names.2 <- c()
    for (col in 1:length(colnames(average_expression))) {
      tmp_names <- cell_markers[cell_markers$cluster %in% cluster,]
      tmp_names <- tmp_names[tmp_names$gen %in% top_cell_markers$gene[top_cell_markers$cluster %in% cluster],]
      tmp_names <- as.data.frame(tmp_names[order(tmp_names$per_obj, decreasing = T), ])
      cluster <- cluster + 1
      cell_names.1[col] <- toupper(tmp_names$gen[1])
      cell_names.2[col] <- firstup(tolower(tmp_names$gen[2]))
    }
    
    
     sub_names <- paste(cell_names.1, cell_names.2)
    
  }
  
  subclass_marker_list_pheat <- rownames(average_expression)[toupper(rownames(average_expression)) %in% toupper(cell_names.2)]
  subclass_marker_list_pheat <- unique(c(subclass_marker_list_pheat, cell_names[!grepl('Bad!',cell_names)]))
  
  return(list('sub_names' = sub_names , 'used_markers' = subclass_marker_list_pheat))
  
}

name_repairing <- function(seurat_project, average_expression, markers_class, markers_subclass, species) {
  
  if (length(markers_class) != 0) {
  
  #remove empty cells (without markers expression)
  
  class_marker_list <- c()
  for (class in markers_class) {
    class_marker_list <- c(class_marker_list, textclean::mgsub(class, c('+'), c('')))
  }
  
  
  class_marker_list <- rownames(seurat_project)[toupper(rownames(seurat_project)) %in% toupper(class_marker_list)]
  
  
  second_matrix <- average_expression[toupper(rownames(average_expression)) %in% toupper(class_marker_list),]
  
  
  renamed_old.1 <- c()
  renamed_new.1 <- c()
  renamed_old.2 <- c()
  renamed_new.2 <- c()
  
  index_marker <- 0
  for (marker in markers_class) {
    index_marker <- index_marker + 1
    for (col in 1:length(colnames(second_matrix))) {
      second_matrix <- as.data.frame(second_matrix[order(second_matrix[,col], decreasing = T), ,drop = F])
      if (max(second_matrix[,col]) == 0 & !grepl('Unknow', as.character(colnames(second_matrix)[col]))) {
        renamed_old.1 <- c(renamed_old.1, colnames(second_matrix)[col])
		 mark <- c()
		 for (change in colnames(markers_class)) {
          mark <- c(mark, change[grepl(change, colnames(second_matrix)[col])])
          if (length(mark) != 0) {break}
        }
        n_marker = 0
        for (marker_2 in markers_class) {
          n_marker = n_marker + 1
          if (grepl(toupper(rownames(second_matrix)[1]), toupper(list(marker_2)))) {
            colnames(second_matrix)[col] <- gsub(pattern = mark, replacement = 'Unknow', x = colnames(second_matrix)[col])
            renamed_new.1 <- c(renamed_new.1, colnames(second_matrix)[col])
            break}
		}
			
      } else if (grepl(as.character(colnames(markers_class)[index_marker]), as.character(colnames(second_matrix)[col])) & !grepl(as.character(toupper(rownames(second_matrix)[1])), toupper(as.character(list(textclean::mgsub(marker, c('+'), c(''))))))) {
        renamed_old.1 <- c(renamed_old.1, colnames(second_matrix)[col])
        mark <- c()
        for (change in colnames(markers_class)) {
          mark <- c(mark, change[grepl(change, colnames(second_matrix)[col])])
          if (length(mark) != 0) {break}
        }
        n_marker = 0
        for (marker_2 in markers_class) {
          n_marker = n_marker + 1
          if (grepl(toupper(rownames(second_matrix)[1]), toupper(list(marker_2)))) {
            colnames(second_matrix)[col] <- gsub(pattern = mark, replacement = colnames(markers_class)[n_marker], x = colnames(second_matrix)[col])
            renamed_new.1 <- c(renamed_new.1, colnames(second_matrix)[col])
            break
          }
        }
      }
    }
  }
  
  Renamed_idents <- as.character(Idents(seurat_project))
  
  
  if (length(renamed_old.1) != 0) {
    n = 0
    for (name in 1:length(renamed_new.1)) {
      n = n + 1
      Renamed_idents[Renamed_idents %in% renamed_old.1[n]] <- renamed_new.1[n]
    }
  }
  
  
  
  
  if (length(markers_subclass) != 0) {
    subclass_marker_list <- c()
    for (subclass in markers_subclass) {
      subclass_marker_list <- c(subclass_marker_list, toupper(subclass))
    }
    
    subclass_marker_list <- rownames(seurat_project)[toupper(rownames(seurat_project)) %in% toupper(subclass_marker_list)]
    
    
    
    old.names <- as.character(colnames(second_matrix))
    
    
    second_matrix <- average_expression[toupper(rownames(average_expression)) %in% toupper(subclass_marker_list),]
    colnames(second_matrix) <- as.character(old.names)
    
    ###########################################################################################################################################################
    #Renamed function
    

    
    for (col in 1:length(colnames(second_matrix))) {
      second_matrix <- as.data.frame(second_matrix[order(second_matrix[,col], decreasing = T), ,drop = F])
      if ((lengths(gregexpr("\\W+", as.character(colnames(second_matrix)[col]))) + 1 == 4) & grepl('Unknow', as.character(colnames(second_matrix)[col]))) {
       
        mark <- c()
        for (change in markers_subclass) {
          mark <- c(mark, change[grepl(paste(change,' ', sep = ''), colnames(second_matrix)[col])])
          if (length(mark) !=0) {
		  renamed_old.2 <- unique(c(renamed_old.2, colnames(second_matrix)[col]))
		  colnames(second_matrix)[col] <- gsub(pattern = mark, replacement =  '', x = colnames(second_matrix)[col])
		  renamed_new.2 <- c(renamed_new.2, colnames(second_matrix)[col])
		  break}
        }
        
        
      } else if (second_matrix[1,col] == 0) {
        renamed_old.2 <- c(renamed_old.2, colnames(second_matrix)[col])
        renamed_new.2 <- c(renamed_new.2, 'BAD!')
      } else if (!grepl(toupper(rownames(second_matrix)[1]), colnames(second_matrix)[col]) & !grepl('Unknow', colnames(second_matrix)[col])) {
        renamed_old.2 <- unique(c(renamed_old.2, colnames(second_matrix)[col]))
        mark <- c()
        for (change in markers_subclass) {
          mark <- c(mark, change[grepl(paste(change,' ', sep = ''), colnames(second_matrix)[col])])
          if (length(mark) !=0) {break}
        }
        
        colnames(second_matrix)[col] <- gsub(pattern = mark, replacement =  toupper(rownames(second_matrix)[1]), x = colnames(second_matrix)[col])
        renamed_new.2 <- c(renamed_new.2, colnames(second_matrix)[col])
      } 
    }
    
    
    if (length(renamed_old.2) != 0) {
      n = 0
      for (name in 1:length(renamed_new.2)) {
        n = n + 1
        Renamed_idents[Renamed_idents %in% renamed_old.2[n]] <- renamed_new.2[n]
      }
    }
    
  }
  
  } else {
  
	Renamed_idents <- as.character(Idents(seurat_project))
	
  }
  
  part_name_1 <- sub(" .*", "", Renamed_idents)
  part_name_2 <- sub('.*? ', "", Renamed_idents)
  
  
  if (species == 'human') {
    part_name_2 <- toupper(part_name_2)
  } else {
    
    part_name_2 <- str_to_title(part_name_2)
  }
  
  Renamed_idents <- paste(part_name_1, part_name_2)
  
   if (length(markers_subclass) == 0) {subclass_marker_list <- c()}
   if (length(markers_class) == 0) {
   class_marker_list <- c()
   renamed_new.1 <- c()
   renamed_new.2 <- c()}
  
  return(list('Renamed_idents' = Renamed_idents, 'renamed_new.1' = renamed_new.1, 'renamed_new.2' = renamed_new.2, 'subclass_marker_list' = subclass_marker_list, 'class_marker_list' = class_marker_list))
  
}


dim_reuction_pcs <- function(dim_stats) {
  
  
    dim <- 1
    score <- c()
    element <- 0
    for (i in dims$`Elbow$data$stdev`) {
      element <- element + 1
      if (i-i*0.01 > dims$`Elbow$data$stdev`[element+1] & element < 50 | i-i*0.02 > dims$`Elbow$data$stdev`[element+2] & element < 49 | i-i*0.02 > dims$`Elbow$data$stdev`[element+3] & element < 48 | i-i*0.02 > dims$`Elbow$data$stdev`[element+4] & element < 47) {
        dim <- dim + 1
      } else 
        break
    }
    dim <- as.numeric(dim)
    
    return(dim)
  
}



bin_cell_test <- function(renamed_list, p_val) {
  
  subclass_names <- renamed_list$Renamed_idents
  bad <- subclass_names[grepl('BAD!', toupper(as.character(subclass_names)))]
  renamed.subnames <- c(as.character(renamed_list$renamed_new.1), as.character(renamed_list$renamed_new.2))
  renamed.subnames <- subclass_names[as.character(toupper(subclass_names)) %in% as.character(toupper(renamed.subnames))]
  renamed.subnames <- renamed.subnames[!as.character(renamed.subnames) %in% as.character(bad)]
  new.subnames <- subclass_names[!as.character(subclass_names) %in% as.character(bad)]
  new.subnames <- new.subnames[!as.character(new.subnames) %in% as.character(renamed.subnames)]
  bad.subnames <- subclass_names[as.character(subclass_names) %in% as.character(bad)]
  subclass_names[subclass_names %in% bad.subnames] <- 'Undefined'
  
  data <- as.data.frame(summary(as.factor(subclass_names), maxsum = length(unique(subclass_names))))
  colnames(data)[1] <- 'n'
  data$names <- rownames(data)
  data$p_val <- NA
  
  for (n in 1:length(data$n)) {
    bin <- binom.test(data$n[n], sum(data$n), p = 1/sum(data$n),
                      conf.level = 0.9)
    data$p_val[n] <- bin$p.value
  }
  
  
  
  below.names <- data$names[data$p_val > p_val]
  data$test[data$names %in% new.subnames] <- "Good marked types" 
  data$test[data$names %in% renamed.subnames] <- "Renamed"
  data$test[data$names %in% 'Undefined'] <- "Undefined types"
  data$test[data$names %in% below.names] <- "Non-significant"
  
  return(list('data' = data, 'below.names' = below.names, 'bad.subnames' = bad.subnames))
  
}


cell_stat_graph <- function(data) {
  
  
  threshold <- ggplot(data, aes(y = n, x = reorder(names, -n), fill = test, sort = test)) +
    geom_bar(stat = 'identity') +
    ylab("Cells types") +
    xlab("Number of cells")+
    theme_bw() +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
    labs(fill = "Cells threshold") +
    coord_flip() + 
    scale_fill_manual(values = c("Good marked types" = "#00BA38", "Renamed" = "#00BF7D", "Undefined types" = "#F8766D",  "Non-significant" = "gray"))
  
  
  
  return(threshold)
  
} 


aggregation_num <- function(seurat_project) {
  
  
  subset_num <- round(length(colnames(seurat_project))/10000)
  cells_num <- round(length(colnames(seurat_project))/subset_num)
  exp_matrix <- GetAssayData(seurat_project, slot = 'data')
  cols <- as.data.frame(summary(seurat_project@active.ident, maxsum = length(unique(colnames(exp_matrix)))))
  cols <- as.data.frame(cols[order(as.numeric(rownames(cols))), , drop = F])
  colnames(cols)[1] <- 'n'
  
  if (round(length(colnames(seurat_project))) > 14999) {
    
    for (batch in 1:subset_num) {
      if (batch == 1 & round(length(colnames(seurat_project))) > 14999) {
        average_expression <- as.data.frame(exp_matrix[,1:cells_num])
        colnames(average_expression) <- seurat_project@active.ident[1:cells_num]
        average_expression <- t(average_expression)
        average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
        rownames(average_expression) <- average_expression[,1]
        average_expression <- average_expression[,-1]
        average_expression <- t(average_expression)
        
      } else if (batch == subset_num & round(length(colnames(seurat_project))) > 14999) {
        average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))])
        colnames(average_expression_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))]
        print((((batch - 1) * cells_num)+1))
        print(length(colnames(seurat_project)))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
        rownames(average_expression_tmp) <- average_expression_tmp[,1]
        average_expression_tmp <- average_expression_tmp[,-1]
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
        average_expression <- as.data.frame(average_expression)
        average_expression <- t(average_expression)
        average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
        rownames(average_expression) <- average_expression[,1]
        average_expression <- average_expression[,-1]
        average_expression <- t(average_expression)
        average_expression <- average_expression[,order(as.numeric(colnames(average_expression)))]
        average_expression <- as.data.frame(average_expression)
        
        num_col <- 0
        for (col in average_expression) {
          num_col <- num_col + 1
          num_row <- 0
          for (mean in col) {
            num_row <- num_row + 1
            average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
          } 
        }
      } else if (batch > 1 & batch < subset_num & round(length(colnames(seurat_project))) > 14999) {
        average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
        colnames(average_expression_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
        print((((batch - 1) * cells_num)+1))
        print(batch * cells_num)
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
        rownames(average_expression_tmp) <- average_expression_tmp[,1]
        average_expression_tmp <- average_expression_tmp[,-1]
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
      } 
    }
  }
  
  if (round(length(colnames(seurat_project))) < 14999) {
    average_expression <- as.data.frame(exp_matrix)
    colnames(average_expression) <- seurat_project@active.ident
    average_expression <- t(average_expression)
    average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
    rownames(average_expression) <- average_expression[,1]
    average_expression <- average_expression[,-1]
    average_expression <- t(average_expression)
    average_expression <- average_expression[,order(as.numeric(colnames(average_expression)))]
    average_expression <- as.data.frame(average_expression)
    
    num_col <- 0
    for (col in average_expression) {
      num_col <- num_col + 1
      num_row <- 0
      for (mean in col) {
        num_row <- num_row + 1
        average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
        
      } 
      
    }
  }
  
  return(average_expression)
}



aggregation_chr <- function(seurat_project) {
  
  subset_num <- round(length(colnames(seurat_project))/10000)
  cells_num <- round(length(colnames(seurat_project))/subset_num)
  exp_matrix <- GetAssayData(seurat_project, slot = 'data')
  cols <- as.data.frame(summary(seurat_project@active.ident, maxsum = length(unique(colnames(exp_matrix)))))
  cols <- as.data.frame(cols[order(rownames(cols)), , drop = F])
  colnames(cols)[1] <- 'n'
  
  if (round(length(colnames(seurat_project))) > 14999){
    
    for (batch in 1:subset_num) {
      if (batch == 1 & round(length(colnames(seurat_project))) > 14999) {
        average_expression <- as.data.frame(exp_matrix[,1:cells_num])
        colnames(average_expression) <- seurat_project@active.ident[1:cells_num]
        average_expression <- t(average_expression)
        average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
        rownames(average_expression) <- average_expression[,1]
        average_expression <- average_expression[,-1]
        average_expression <- t(average_expression)
      } else if (batch == subset_num & round(length(colnames(seurat_project))) > 14999) {
        average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))])
        colnames(average_expression_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))]
        print((((batch - 1) * cells_num)+1))
        print(length(colnames(seurat_project)))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
        rownames(average_expression_tmp) <- average_expression_tmp[,1]
        average_expression_tmp <- average_expression_tmp[,-1]
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
        average_expression <- as.data.frame(average_expression)
        average_expression <- t(average_expression)
        average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
        rownames(average_expression) <- average_expression[,1]
        average_expression <- average_expression[,-1]
        average_expression <- t(average_expression)
        average_expression <- average_expression[,order(colnames(average_expression))]
        average_expression <- as.data.frame(average_expression)
        
        num_col <- 0
        for (col in average_expression) {
          num_col <- num_col + 1
          num_row <- 0
          for (mean in col) {
            num_row <- num_row + 1
            average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
          } 
        }
        
      } else if (batch > 1 & batch < subset_num & round(length(colnames(seurat_project))) > 14999) {
        average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
        colnames(average_expression_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
        print((((batch - 1) * cells_num)+1))
        print(batch * cells_num)
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
        rownames(average_expression_tmp) <- average_expression_tmp[,1]
        average_expression_tmp <- average_expression_tmp[,-1]
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
      } 
    }
  }
  
  if (round(length(colnames(seurat_project))) < 14999) {
    average_expression <- as.data.frame(exp_matrix)
    colnames(average_expression) <- seurat_project@active.ident
    average_expression <- t(average_expression)
    average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
    rownames(average_expression) <- average_expression[,1]
    average_expression <- average_expression[,-1]
    average_expression <- t(average_expression)
    average_expression <- average_expression[,order(colnames(average_expression))]
    average_expression <- as.data.frame(average_expression)
    
    num_col <- 0
    for (col in average_expression) {
      num_col <- num_col + 1
      num_row <- 0
      for (mean in col) {
        num_row <- num_row + 1
        average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
        
      } 
      
    }
  }
  
  return(average_expression)
}



heterogenity_stats <- function(seurat_project) {
  

  subset_num <- round(length(colnames(seurat_project))/1000)
  cells_num <- round(length(colnames(seurat_project))/subset_num)
  exp_matrix_obl <- GetAssayData(seurat_project, slot = 'data')
  
  if (round(length(colnames(seurat_project))) > 14999){
    
    for (batch in 1:subset_num) {
      if (batch == 1 & round(length(colnames(seurat_project))) > 14999) {
        exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,1:cells_num])
        colnames(exp_matrix_obl_tmp) <- seurat_project@active.ident[1:cells_num]
        mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
        exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
        positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
        sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
        positive_expression_perc <- positive_expression/sum_expression
        names <- colnames(exp_matrix_obl_tmp)
        exp_stat <- data.frame(names,mean_expression, positive_expression_perc)
      } else if (batch == subset_num & round(length(colnames(seurat_project))) > 14999) {
        exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix_obl))])
        colnames(exp_matrix_obl_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix_obl))]
        mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
        exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
        positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
        sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
        positive_expression_perc <- positive_expression/sum_expression
        names <- colnames(exp_matrix_obl_tmp)
        exp_stat_tmp <- data.frame(names,mean_expression, positive_expression_perc)
        exp_stat <- rbind(exp_stat, exp_stat_tmp)
      } else if (batch > 1 & batch < subset_num & round(length(colnames(seurat_project))) > 14999) {
        exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
        colnames(exp_matrix_obl_tmp) <- seurat_project@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
        mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
        exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
        positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
        sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
        positive_expression_perc <- positive_expression/sum_expression
        names <- colnames(exp_matrix_obl_tmp)
        exp_stat_tmp <- data.frame(names,mean_expression, positive_expression_perc)
        exp_stat <- rbind(exp_stat, exp_stat_tmp)
      } 
    }
  } 
  
  if (round(length(colnames(seurat_project))) < 14999){
    exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl)
    colnames(exp_matrix_obl_tmp) <- seurat_project@active.ident
    mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
    exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
    positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
    sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
    positive_expression_perc <- positive_expression/sum_expression
    names <- colnames(exp_matrix_obl_tmp)
    exp_stat <- data.frame(names,mean_expression, positive_expression_perc)
  }
  
  return(exp_stat)
}



hd_cluster_factors <- function(seurat_object, markers_cssg) {
  
  library(umap)
  
  for (cluster in unique(Idents(seurat_object))) {
    
    

    tmp_cssg <- markers_cssg[markers_cssg$cluster %in% cluster,]
    
    tmp_cssg <- tmp_cssg[order(tmp_cssg$adj_hf, decreasing = TRUE),]
    tmp_cssg <- tmp_cssg[tmp_cssg$loss_pval == min(tmp_cssg$loss_pval),]

    gen_cor <- c()
    for (i in 1:nrow(tmp_cssg)) {
      gen_cor <- c(gen_cor, gsub(' ', '', strsplit(rownames(tmp_cssg)[i], split = ' ')[[1]]))
    }
    
    gen_cor <- unique(gen_cor)
    
    seurat_object_sub <- subset(seurat_object, idents = cluster, features = gen_cor)
    
    tmp <- as.data.frame(GetAssayData(seurat_object_sub, slot = 'data'))

    
    tmp_df <- umap(t(as.data.frame(tmp)), method = 'umap-learn')$layout
    colnames(tmp_df) <- c('x', 'y')
    tmp_df <- as.data.frame(tmp_df)
    tmp_df$BARCODES <- rownames(tmp_df)
    
    if (exists('df_factor') == TRUE) {
      df_factor <- rbind(df_factor, tmp_df)
      
      
    }
    else {
      df_factor <- tmp_df
      
      
    }
  
  }

  return(df_factor)
}


hdmap_cordinates <- function(seurat_object, fUMAP) {
  
  
  HDMAP <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
  HDMAP$idents <- Idents(seurat_object)
  HDMAP$BARCODES <- rownames(HDMAP) 
  
  HDMAP$UMAP_1 <- (HDMAP$UMAP_1 - min(HDMAP$UMAP_1)) / (max(HDMAP$UMAP_1) - min(HDMAP$UMAP_1)) *100
  HDMAP$UMAP_2 <- (HDMAP$UMAP_2 - min(HDMAP$UMAP_2)) / (max(HDMAP$UMAP_2) - min(HDMAP$UMAP_2)) *100
  

  
  agg_umap <- aggregate(HDMAP[,1:2], by = list(HDMAP[,3]), FUN = mean)
  colnames(agg_umap) <- c('idents', 'avg_UMAP1', 'avg_UMAP2')
  
  
  
  HDMAP <- merge(HDMAP, agg_umap, by = 'idents', all = T)
  
  HDMAP <- merge(HDMAP, fUMAP, by = 'BARCODES',all.x = TRUE)
  
  HDMAP$x <- (HDMAP$x - min(HDMAP$x)) / (max(HDMAP$x) - min(HDMAP$x)) 
  HDMAP$y <- (HDMAP$y - min(HDMAP$y)) / (max(HDMAP$y) - min(HDMAP$y)) 
  
  HDMAP$HDMAP_1 <- HDMAP$avg_UMAP1 + HDMAP$x
  HDMAP$HDMAP_2 <- HDMAP$avg_UMAP2 + HDMAP$y
  
  
  return(HDMAP)
  
}



DimPlotFactor <- function(HDMAP) {
  
  
  plot <- ggplot() +
    geom_point(mapping = aes(x = HDMAP$HDMAP_1  , y = HDMAP$HDMAP_2 , color = HDMAP$idents)) +
    xlab('HDMAP_1') +
    ylab('HDMAP_2') +
    labs(colour="Cell subtypes") +
    theme_classic()
  
  
  return(plot)
  
}



