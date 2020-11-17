all_unique <- function(data, range, value_check_1){
  eval(parse(text= paste0(data,'<<- aggregate(',data,range,' ,list(',
                          data,'$',value_check_1,'),FUN=list)')))
}



markers_patterning <- function(top, markers){


var_num <- as.numeric(length(colnames(markers)))

{
top$cell <- NA
index <- 0
ind <- 0
for (g in top$gene) {
  index <- index+1
  for (i in markers[,1:var_num]){
    ind <- (ind+1)
    for (c in g) {
      for (e in i) {
        if (c %in% e) {
          index_markes <- ind-((index-1)*var_num)
          top$cell[index] <- paste(if (!is.na(top$cell[index])) {top$cell[index]}, colnames(markers)[index_markes])
          print(index_markes)
          print(index)
        } 
      }
    }
  }  
} 
}


  for (i in 1:length(top$cluster)) {
     top$cell[i] <- paste(if (!is.na(top$cell[i])) {top$cell[i]}, top$gene[[i]][1], top$gene[[i]][2])

   }


  top <<- top
  
}



cell_names <- function(input_seurat, cell_markers){
  
  
  matrix_a <- as.matrix(GetAssayData(input_seurat, slot = 'counts'))
  matrix_a <- as.data.frame(matrix_a)
  colname <- colnames(matrix_a)
  matrix_a$genes <- toupper(rownames(matrix_a))
  
  
  
  
  
  
  var_num <- as.numeric(length(colnames(cell_markers)))
  var_col <- as.numeric(length(colnames(matrix_a)))-1
  gen <- c()
  index <- 0
  tmp <- 0
  for (marker in 1:var_num) {
    col <- 0
    tmp <- c()
    gen <- as.list(cell_markers[,marker])
    for (i in gen){
      tmp <- matrix_a[matrix_a$genes %in% i,]
      for (value in tmp[,1:var_col]) {
        col <-  col + 1
        value  <- sum(value, na.rm = T)
        print(value)
        if (value > 0){
          colnames(tmp)[col] <- paste(colnames(tmp)[col], colnames(cell_markers[marker]))
          colnames(matrix_a) <- colnames(tmp)
          
        } 
      } 
    }  
  } 
  
  #remove genes
  matrix_a <- matrix_a[, 1:var_col]
  
  for (rm in colname) {
    colnames(matrix_a) <- gsub(colnames(matrix_a), pattern = rm, replacement = '')
  }
  
expression_matrix_cells <<- matrix_a
  
  
}