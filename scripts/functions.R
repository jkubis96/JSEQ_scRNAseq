#Aggregate differenced markers
all_unique <- function(data, range, value_check_1){
  eval(parse(text= paste0(data,'<<- aggregate(',data,range,' ,list(',
                          data,'$',value_check_1,'),FUN=list)')))
}



#Cell name change
# cell_nameing <- function(matrix_a, markers) {
#   
#   
#   matrix_a <- as.data.frame(matrix_a)
#   colname <- colnames(matrix_a)
#   matrix_a$genes <- toupper(rownames(matrix_a))
#   
#   
#   
#   var_num <- as.numeric(length(colnames(markers)))
#   var_col <- as.numeric(length(colnames(matrix_a)))-1
#   gen <- c()
#   index <- 0
#   tmp <- 0
#   for (marker in 1:var_num) {
#     col <- 0
#     tmp <- c()
#     gen <- as.list(markers[,marker])
#     for (i in gen){
#       tmp <- matrix_a[matrix_a$genes %in% textclean::mgsub(i, c('+','-','*'), c('','','')),]
#       mark_genes_plus <- matrix_a[matrix_a$genes %in% textclean::mgsub(i, c('+'), c('')),]
#       mark_genes_minus <- matrix_a[matrix_a$genes %in% textclean::mgsub(i, c('-'), c('')),]
#       mark_genes_obligatory <- matrix_a[matrix_a$genes %in% textclean::mgsub(i, c('*'), c('')),]
#       for (value in 1:var_col) {
#         if (nrow(mark_genes_plus) == 0) {value_plus <- 1} 
#         else if (nrow(mark_genes_plus) != 0) {value_plus  <- sum(mark_genes_plus[,value], na.rm = T)}
#         if (nrow(mark_genes_minus) == 0) {value_minus <- 0} 
#         else if (nrow(mark_genes_minus) != 0) {value_minus <- sum(mark_genes_minus[,value], na.rm = T)}
#         if (nrow(mark_genes_obligatory) == 0) {value_obligatory <- 1} 
#         else if (nrow(mark_genes_obligatory) != 0) {value_obligatory <- prod(mark_genes_obligatory[,value], na.rm = T)}
#         value_plus <- value_plus*value_obligatory
#         value_minus <- value_minus*value_obligatory
#         if ((value_plus > 0)&(value_minus == 0)){
#           colnames(tmp)[value] <- paste(colnames(tmp)[value], colnames(markers[marker]))
#           colnames(matrix_a) <- colnames(tmp)}
#         
#       } 
#     }  
#   } 
#   
#   
#   
#   #remove genes
#   matrix_a <- matrix_a[, 1:var_col]
#   
#   `%!in%` = Negate(`%in%`)
#   
#   
#   for (rm in colname) {
#     
#     rm <- paste0(rm,' ')
#     colnames(matrix_a) <- gsub(colnames(matrix_a), pattern = rm, replacement = '')
#   }
#   
#   # col <- 0
#   # for (un in colnames(matrix_a)) {
#   #   col <- col + 1
#   #   if (un %!in% (colnames(markers))) {colnames(matrix_a)[col] <-'Unknow'}
#   # }
#   # 
#   
#   exp_matrix <<- matrix_a
#   rm(matrix_a)
#   
# }
# 
# 
# 
# ########################################################################
# 
# 
# 
# 
# cluster_markers <- function(top, markers){
#   
#   
#   var_num <- as.numeric(length(colnames(markers)))
#   
#   {
#     top$cell <- NA
#     index <- 0
#     ind <- 0
#     for (g in top$gene) {
#       index <- index+1
#       for (i in markers[,1:var_num]){
#         ind <- (ind+1)
#         for (c in g) {
#           for (e in i) {
#             if (c %in% e) {
#               index_markes <- ind-((index-1)*var_num)
#               top$cell[index] <- paste(if (!is.na(top$cell[index])) {top$cell[index]}, colnames(markers)[index_markes])
#               print(index_markes)
#               print(index)
#             } 
#           }
#         }
#       }  
#     } 
#   }
#   
#   
#   for (i in 1:length(top$cluster)) {
#     top$cell[i] <- paste(if (!is.na(top$cell[i])) {top$cell[i]}, top$gene[[i]][1])
#     
#   }
#   
#   
#   group_names <<- top
#   
# }





#Cluster nameing
# cluster_nameing<- function(matrix_a, markers) {
#   
#   
#   matrix_a <- as.data.frame(matrix_a)
#   colname <- colnames(matrix_a)
#   rownames(matrix_a) <- make.unique(toupper(rownames(matrix_a)), sep = '')
#   
#   
#   
#   var_num <- as.numeric(length(colnames(markers)))
#   var_col <- as.numeric(length(colnames(matrix_a)))
#   gen <- c()
#   index <- 0
#   tmp <- 0
#   for (marker in 1:var_num) {
#     col <- 0
#     tmp <- c()
#     gen <- as.list(markers[,marker])
#     for (i in gen){
#       tmp <- floor((matrix_a[textclean::mgsub(i, c('+','-','*'), c('','','')),]))
#       mark_genes_plus <- floor((matrix_a[textclean::mgsub(i, c('+'), c('')),]))
#       mark_genes_minus <- floor((matrix_a[matrix_a$genes %in% textclean::mgsub(i, c('-'), c('')),]))
#       mark_genes_obligatory <- floor((matrix_a[textclean::mgsub(i, c('*'), c('')),]))
#       for (value in 1:var_col) {
#         if (nrow(mark_genes_plus) == 0) {value_plus <- 1}
#         else if (nrow(mark_genes_plus) > 1) {value_plus  <- sum(mark_genes_plus[,value], na.rm = T)}
#         if (nrow(mark_genes_minus) == 0) {value_minus <- 0}
#         else if (nrow(mark_genes_minus) > 1) {value_minus <- sum(mark_genes_minus[,value], na.rm = T)}
#         if (nrow(mark_genes_obligatory) == 0) {value_obligatory <- 1}
#         else if (nrow(mark_genes_obligatory) > 1) {value_obligatory <- prod(mark_genes_obligatory[,value], na.rm = T)}
#         value_plus <- value_plus*value_obligatory
#         value_minus <- value_minus*value_obligatory
#         if ((value_plus > 0)&(value_minus == 0)){
#           colnames(tmp)[value] <- paste(colnames(tmp)[value], colnames(markers[marker]))
#           colnames(matrix_a) <- colnames(tmp)}
#         
#       }
#     }
#   }
#   
#   
#   
#   #remove genes
#   matrix_a <- matrix_a[, 1:var_col]
#   
#   `%!in%` = Negate(`%in%`)
#   
#   
#   for (rm in colname) {
#     
#     rm <- paste0(rm,' ')
#     colnames(matrix_a) <- gsub(colnames(matrix_a), pattern = rm, replacement = '')
#   }
#   
#   
#   average_expression <<- matrix_a
#   rm(matrix_a)
#   
# }

###################################################################################



# #Cell name change for future
# cluster_nameing<- function(matrix_a, markers) {
#   
#   
#   matrix_a <- as.data.frame(matrix_a)
#   colname <- colnames(matrix_a)
#   rownames(matrix_a) <- make.unique(toupper(rownames(matrix_a)), sep = '')
#   
#   
#   
#   var_num <- as.numeric(length(colnames(markers)))
#   var_col <- as.numeric(length(colnames(matrix_a)))
#   gen <- c()
#   index <- 0
#   tmp <- 0
#   for (marker in 1:var_num) {
#     col <- 0
#     tmp <- c()
#     gen <- as.list(markers[,marker])
#     for (i in gen){
#       tmp <- floor((matrix_a[textclean::mgsub(i, c('+','-','*'), c('','','')),]))
#       mark_genes_plus <- floor((matrix_a[textclean::mgsub(i, c('+'), c('')),]))
#       mark_genes_minus <- floor((matrix_a[textclean::mgsub(i, c('-'), c('')),]))
#       for (value in 1:var_col) {
#         if (nrow(mark_genes_plus) == 0) {value_plus <- 1}
#         else if (nrow(mark_genes_plus) != 0) {value_plus  <- sum(mark_genes_plus[,value], na.rm = T)}
#         if (nrow(mark_genes_minus) == 0) {value_minus <- 0}
#         else if (nrow(mark_genes_minus) != 0) {value_minus <- sum(mark_genes_minus[,value], na.rm = T)}
#         if ((value_plus > 1 & value_minus < 1)){
#           colnames(tmp)[value] <- paste(colnames(tmp)[value], colnames(markers[marker]))
#           colnames(matrix_a) <- colnames(tmp)}
#         
#       }
#     }
#   }
#   
#   
#   
#   #remove genes
#   matrix_a <- matrix_a[, 1:var_col]
#   
#   `%!in%` = Negate(`%in%`)
#   
#   
#   for (rm in colname) {
#     
#     rm <- paste0(rm,' ')
#     colnames(matrix_a) <- gsub(colnames(matrix_a), pattern = rm, replacement = '')
#   }
#   
#   # col <- 0
#   # for (un in colnames(matrix_a)) {
#   #   col <- col + 1
#   #   if (un %!in% (colnames(markers))) {colnames(matrix_a)[col] <-'Unknow'}
#   # }
#   #
#   
#   average_expression <<- matrix_a
#   rm(matrix_a)
#   
# }





#Cell name change for future
cluster_nameing<- function(matrix_a, markers) {
  
  
  matrix_a <- as.data.frame(matrix_a)
  colname <- colnames(matrix_a)
  rownames(matrix_a) <- make.unique(toupper(rownames(matrix_a)), sep = '')
  
  list.markers <- c()
  for (l.marker in markers) {
    for (marker in l.marker) {
      TF <- (grepl('+', marker, fixed = T))
      if (TF == TRUE)  {
        list.markers <- c(list.markers, textclean::mgsub(marker, c('+'), c('')))
      }
    } 
  } 
  
  
  num <- as.numeric(length(list.markers))
  index = 0
  for (i in colnames(matrix_a)) {
    index = index +1
    rename_df <- c()
    rename_df <- matrix_a[toupper(rownames(matrix_a)) %in% list.markers,]
    rename_df <- as.data.frame(rename_df[order(rename_df[,index], decreasing = T), ,drop = F])
    if (sum(rename_df[,index] > 0)) {
      colnames(matrix_a)[index] <- rownames(rename_df[1,])
    } else colnames(matrix_a)[index] <- 'Unknow'
  }  
  matrix_b <- round(matrix_a, digits = 1) 
  
  col <- 0
  for (c.marker in markers) {
    col <- col + 1
    cell_n <- 0
    for (marker in c.marker) {
      for (cell in colnames(matrix_a)) {
        cell_n <- cell_n +1 
        if (cell %in% textclean::mgsub(marker, c('+'), c(''))) {
          colnames(matrix_a)[cell_n] <- colnames(markers[col])
        }
      }
    }
  }
 
  average_expression <<- matrix_a
  rm(matrix_a)
  # 
}