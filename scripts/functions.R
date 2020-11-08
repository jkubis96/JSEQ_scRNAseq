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
