#Aggregate differenced markers
all_unique <- function(data, range, value_check_1){
  eval(parse(text= paste0(data,'<<- aggregate(',data,range,' ,list(',
                          data,'$',value_check_1,'),FUN=list)')))
}





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
  
  
  index = 0
  for (i in colnames(matrix_a)) {
    index = index +1
    rename_df <- as.data.frame(matrix_a[rownames(matrix_a) %in% list.markers,index , drop = F])
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
        cell_n <- cell_n +1 
        if (cell %in% textclean::mgsub(marker, c('+'), c(''))) {
          new.names <- colnames(markers[col])
          colnames(matrix_a)[cell_n] <- new.names
        }
      }
    }
  }
  
  
  average_expression <<- matrix_a
  rm(matrix_a)
  
}

