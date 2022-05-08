load_annotation <- function(path) {
  
  options(scipen = 999)
  library(readr)
  library(stringr)
  library(dplyr)
  library(doSNOW)
  library(foreach)
  library(doParallel)

  cat('\n\n Data loading... \n\n')


  GTF <- read_delim(path, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, comment = '#', 
                    col_types = cols(
                      X1 = col_character(),
                      X2 = col_character(),
                      X3 = col_character(),
                      X4 = col_integer(),
                      X5 = col_integer(),
                      X6 = col_character(),
                      X7 = col_character(),
                      X8 = col_character(),
                      X9 = col_character()
                    ))
  
  GTF <- GTF[toupper(GTF$X3) %in% c("GENE", "TRANSCRIPT", "EXON", "CDS", 'MRNA', 'FIVE_PRIM_UTR', 'THREE_PRIM_UTR'), ]

return(GTF)

}

 


create_GTF_df <- function(input) {
    cat('\n\n GTF converting... \n\n')

    df <- input[,1:8]
    
   
    
    if (TRUE %in% unique(grepl('gene_name', input$X9))) {
    #GENECODE & https://www.ensembl.org/index.html
    #############################################################
      
    df$gene_id_check <- grepl('gene_id', input$X9)
    df$gene_id <-  gsub(".*=", "", gsub(";.*", "", gsub(".*gene_id", "", input$X9)))

    df$gene_name_check <- grepl('gene_name', input$X9)
    df$gene_name <-  gsub(".*=", "", gsub(";.*", "", gsub(".*gene_name", "", input$X9)))

    df$transcript_name_check <- grepl('transcript_name', input$X9)
    df$transcript_name <- gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_name", "", input$X9))))
    
    df$transcript_id_check <-grepl('transcript_id', input$X9)
    df$transcript_id <- gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_id", "", input$X9))))
    
  
    
    df$gene_id[df$gene_id == FALSE] <- df$transcript_id[df$gene_id == FALSE]
    
    df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
    
    df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
    
    df$transcript_id[df$transcript_id == FALSE] <- df$gene_id[df$transcript_id == FALSE]
    
    df$transcript_name[df$transcript_name_check == FALSE] <- df$gene_name[df$transcript_name_check == FALSE]
    
    #repaire gene names
    gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id')]
    gene_names <- gene_names[!duplicated(gene_names),]
    gene_names <- gene_names[gene_names$gene_name != gene_names$gene_id, ]
    
    cat('\n\n Unifying the names of genes... \n\n')
    
    iterations <- length(gene_names$gene_name)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    
    CPU <- detectCores() - 2
    
    cl <- makeCluster(CPU)
    
    
    registerDoParallel(cl)
    registerDoSNOW(cl)
    
    
   
  
    df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
      tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
      tmp$gene_name <- gene_names$gene_name[n]
      tmp
      
    }
    
    close(pb)
    stopCluster(cl)  
    
    
    df <- df2 
    
    rm(df2)
    
    
    
    df$lp <- 1:length(df$X1)
    
    #repair duplicated genes names from different loci
    
    gen_list <- c()
    check <- 0
    for (gen in df$gene_name) {
      if (gen != check) {
        gen_list <- c(gen_list, gen)
        check <- gen
      }
      
    }
    
    duplicated_names <- gen_list[duplicated(gen_list)]
    
    if (length(duplicated_names) > 0) {
      dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
      dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
      dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
      dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
      
      genes <- unique(duplicated_names)
      global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
      cat('\n\n Duplicated genes repairing...             \n\n ')
      
      for (strand in c('+','-')) {
        tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
        for (gen in genes){
          tmp <- tmp_strand[tmp_strand$gene_name %in% gen, ]
          group <- 1
          for (i in 1:length(tmp_strand$gene_name)) {
            cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp_strand$gene_name)*100),2), '%           '))
            if (i < length(tmp$gene_name)) {
              tmp1 <- tmp[i,]
              if (i == 1) {
                df_group <- tmp1
                df_group$group <- group}
              tmp2 <- tmp[i+1,]
              if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - 110000 >  tmp1$X5[1]) {
                group <- group + 1
                tmp2$group <- group
                df_group <- rbind(df_group, tmp2)
              } else {
                tmp2$group <- group
                df_group <- rbind(df_group, tmp2)
              }
            } else if (i == length(tmp$gene_name)){
              
              global_df <- rbind(global_df, df_group)
              rm(df_group)
              
            }
          } 
        }
      }
      
      
  
      global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
      
      
      
      
      #rename genes
      
      for (new_name in 1:length(global_df$lp)) {
        cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
        df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
      }
      
      
      
      rm(global_df)
      
    }
    
    

    
    df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
    
    } else if (TRUE %in% unique(grepl('=gene-', input$X9))) {
      
    #NCBI
    ################################################################
    
      
      
      
      df$gene_id_check <-grepl('GeneID:', input$X9)
      df$gene_id <- gsub(".*?\\s", "", gsub(",.*", "", gsub(".*GeneID:", "", input$X9)))
      
      df$gene_name_check <- grepl('=gene-', input$X9)
      df$gene_name <-  gsub(";.*", "", gsub(".*=gene-", "", input$X9))
      
      df$gene_name_check2 <- grepl(';gene=', input$X9)
      df$gene_name2 <-  gsub(";.*", "", gsub(".*;gene=", "", input$X9))
      
      df$transcript_name <- NA
      
      df$transcript_id_check <- grepl('Genbank:', input$X9)
      df$transcript_id <-  gsub(";.*", "", gsub(",.*", "", gsub(".*Genbank:", "", input$X9)))
      
      
      
      df$gene_id[df$gene_id == FALSE] <- df$transcript_id[df$gene_id == FALSE]
      
      df$gene_name[df$gene_name_check == FALSE & df$gene_name_check2 == TRUE] <- df$gene_name2[df$gene_name_check2 == TRUE & df$gene_name_check == FALSE]
      
      df$gene_name_check[df$gene_name_check2 == TRUE] <- TRUE
      
      df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
      
      df$transcript_id[df$transcript_id == FALSE] <- df$gene_id[df$transcript_id == FALSE]
      
      
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id')]
      gene_names <- gene_names[!duplicated(gene_names),]
      gene_names <- gene_names[gene_names$gene_name != gene_names$gene_id, ]
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      #DOBRZE ZROBIONE ALE TRZEBA ZAMIENIÆ I ZBINDOWAC
      
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp

      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- df2 
      
      rm(df2)
      
    
      
      
      df$transcript_name <- df$gene_name
      
      
      
    
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      gen_list <- c()
      check <- 0
      for (gen in df$gene_name) {
        if (gen != check) {
          gen_list <- c(gen_list, gen)
          check <- gen
        }
        
      }
      
      duplicated_names <- gen_list[duplicated(gen_list)]
      
      if (length(duplicated_names) > 0) {
      cat('\n\n Duplicated genes repairing...             \n\n')
        
      dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
      dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
      dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
      dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
      
      genes <- unique(duplicated_names)
      global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
      
      for (strand in c('+','-')) {
        tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
        for (gen in genes){
          tmp <- tmp_strand[tmp_strand$gene_name %in% gen, ]
          group <- 1
          for (i in 1:length(tmp_strand$gene_name)) {
            cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp_strand$gene_name)*100),2), '%           '))
            if (i < length(tmp$gene_name)) {
              tmp1 <- tmp[i,]
              if (i == 1) {
                df_group <- tmp1
                df_group$group <- group}
              tmp2 <- tmp[i+1,]
              if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - 110000 >  tmp1$X5[1]) {
                group <- group + 1
                tmp2$group <- group
                df_group <- rbind(df_group, tmp2)
              } else {
                tmp2$group <- group
                df_group <- rbind(df_group, tmp2)
              }
            } else if (i == length(tmp$gene_name)){
              
              global_df <- rbind(global_df, df_group)
              rm(df_group)
              
            }
          } 
        }
      }
      
      

      
      global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
      
      
      
      #rename genes
      
      for (new_name in 1:length(global_df$lp)) {
        cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
        df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
      }
      
      
      
      rm(global_df)
      
      }
      
    
      

      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check2', 'gene_name2', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
  
  
    } else if (TRUE %in% unique(grepl('gene:', input$X9)) & TRUE %in% unique(grepl('transcript:', input$X9))) {
      
      #CUSTOM
      ################################################################
      
      
    
      df$gene_id <- NA
      
      df$gene_name_check <- grepl('gene:', input$X9)
      df$gene_name <-  gsub(";.*", "", gsub(".*gene:", "", input$X9))
    
      df$transcript_name_check <- grepl('transcript:', input$X9)
      df$transcript_name <-  gsub(";.*", "", gsub(",.*", "", gsub(".*transcript:", "", input$X9)))
      
      df$transcript_id <- NA
      
      
      
      df <- df[!(df$gene_name_check == FALSE & df$transcript_name_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$transcript_name[df$gene_name_check == FALSE ]
      df$transcript_name[df$transcript_name_check == FALSE] <-  df$gene_name[df$transcript_name_check == FALSE ]
      
      
      
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'transcript_name')]
      gene_names <- gene_names[!duplicated(gene_names),]
      gene_names <- gene_names[gene_names$gene_name != gene_names$transcript_name, ]
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
    
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$transcript_name %in% gene_names$transcript_name[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp
        
      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- df2 
      
      rm(df2)
      
      
      
      
      df$transcript_id <- df$gene_name
      df$gene_id <- df$gene_name
      
      
      
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      gen_list <- c()
      check <- 0
      for (gen in df$gene_name) {
        if (gen != check) {
          gen_list <- c(gen_list, gen)
          check <- gen
        }
        
      }
      
      duplicated_names <- gen_list[duplicated(gen_list)]
      
      if (length(duplicated_names) > 0) {
        cat('\n\n Duplicated genes repairing...             \n\n')
        
        dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
        dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
        dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
        dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
        
        genes <- unique(duplicated_names)
        global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
        
        for (strand in c('+','-')) {
          tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
          for (gen in genes){
            tmp <- tmp_strand[tmp_strand$gene_name %in% gen, ]
            group <- 1
            for (i in 1:length(tmp_strand$gene_name)) {
              cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp_strand$gene_name)*100),2), '%           '))
              if (i < length(tmp$gene_name)) {
                tmp1 <- tmp[i,]
                if (i == 1) {
                  df_group <- tmp1
                  df_group$group <- group}
                tmp2 <- tmp[i+1,]
                if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - 110000 >  tmp1$X5[1]) {
                  group <- group + 1
                  tmp2$group <- group
                  df_group <- rbind(df_group, tmp2)
                } else {
                  tmp2$group <- group
                  df_group <- rbind(df_group, tmp2)
                }
              } else if (i == length(tmp$gene_name)){
                
                global_df <- rbind(global_df, df_group)
                rm(df_group)
                
              }
            } 
          }
        }
        
        
        
        
        global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
        
        
        
        #rename genes
        
        for (new_name in 1:length(global_df$lp)) {
          cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
          df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
        }
        
        
        
        rm(global_df)
        
      }
      
      
      
      
      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check2', 'gene_name2', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
      
      
    }
    
    
    
    
    
    colnames(df) <- c('chr','source','annotationType','start','end','score','strand','phase','gene_id','gene_name','transcript_name','transcript_id')

return(df)

}






add_UTR <- function(input, five_prim_utr_length, three_prim_utr_length) {
  
  cat('\n\n UTRs sequence extending...             \n\n')
  
  input$sort_val <- 1:length(input$chr)
  chromosomes <- unique(input$chr)
  
  CDS <- input[toupper(input$annotationType) %in% c("EXON", "CDS"), ]
  
  
  df <- data.frame(matrix(ncol = length(colnames(CDS)), nrow = 0))
  colnames(df) <- colnames(CDS)
  
  for (chr in chromosomes) {
  
  tmp_primary <- CDS[CDS$chr %in% chr,]
  tmp_primary <- tmp_primary[order(tmp_primary$start),]
  tmp_primary$lp <- 1:length(tmp_primary$chr)
  
  
  for (sen in c('+','-')) {
  tmp <- tmp_primary[tmp_primary$strand %in% sen,]
  gen_list <- c()
  check <- 0
  for (gen in tmp$gene_name) {
    if (gen != check) {
      gen_list <- c(gen_list, gen)
      check <- gen
    }
    
  }
  
  if (length(gen_list) > 0) {
    
  
  for (i in 1:length(gen_list)) {
          
          cat('\r',paste('LOC:',chr, '| STRAND:', sen , '| PROGRESS:', round((i/length(gen_list)*100),2), '%           '))
          
                if (i < length(gen_list)) {
                  tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
                  tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
                  min_val <- min(tmp2$start)
                  max_val <- max(tmp1$end)
                } else if (i == length(gen_list)) {
                  tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
                  tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
                  min_val <- 0
                  max_val <- 0
                } 
            
            if (min_val > max_val + (five_prim_utr_length + three_prim_utr_length)/3 & i == 1) {
              tmp1 <- tmp1[tmp1$end < min(tmp2$start),]
              tmp_primary <- tmp_primary[!tmp_primary$lp %in% tmp1$lp,]

              if (sen == '+') {
                
                tmp0 <- tmp1[tmp1$gene_name %in% gen_list[i],]
                tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
                tmp0 <- tmp0[tmp0$lp == min(tmp0$lp),]
                tmp0$end <- as.integer(tmp0$start - 1)
                tmp0$start <- as.integer(tmp0$end - five_prim_utr_length)
                if (tmp0$start < 0) {tmp0$start <- 0}
                tmp0$source <- 'JBIO'
                tmp0$annotationType <- 'five_prim_UTR'
                tmp0$score <- '.'
                tmp0$phase <- '.'
                
                if (min(tmp2$start) > (max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)+2)) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length,0)
                } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/2) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length*0.5,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length*0.5,0)
                } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/3) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length*0.3,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length*0.3,0)
                }
                
                tmp1 <- tmp1[tmp1$end == max(tmp1$end),]
                tmp1 <- tmp1[tmp1$lp == min(tmp1$lp),]
                tmp1$start <- as.integer(tmp1$end + 1)
                tmp1$end <- as.integer(tmp1$end + three_prim_utr_length_cor)
                tmp1$source <- 'JBIO'
                tmp1$annotationType <- 'three_prim_UTR'
                tmp1$score <- '.'
                tmp1$phase <- '.'
                tmp2 <- tmp2[tmp2$start == min(tmp2$start),]
                tmp2 <- tmp2[tmp2$lp == min(tmp2$lp),]
                tmp2$source <- 'JBIO'
                tmp2$annotationType <- 'five_prim_UTR'
                tmp2$score <- '.'
                tmp2$phase <- '.'
                tmp2$end <- as.integer(tmp2$start - 1)
                tmp2$start <- as.integer(tmp2$start - five_prim_utr_length_cor)
                
              } else if (sen == '-') {
                
                tmp0 <- tmp1[tmp1$gene_name %in% gen_list[i],]
                tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
                tmp0 <- tmp0[tmp0$lp == min(tmp0$lp),]
                tmp0$end <- as.integer(tmp0$start - 1)
                tmp0$start <- as.integer(tmp0$end - three_prim_utr_length)
                if (tmp0$start < 0) {tmp0$start <- 0}
                tmp0$source <- 'JBIO'
                tmp0$annotationType <- 'three_prim_UTR'
                tmp0$score <- '.'
                tmp0$phase <- '.'
                
                if (min(tmp2$start) > (max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)+2)) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length,0)
                } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/2) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length*0.5,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length*0.5,0)
                } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/3) {
                  five_prim_utr_length_cor <- round(five_prim_utr_length*0.3,0)
                  three_prim_utr_length_cor <- round(three_prim_utr_length*0.3,0)
                }
                
                tmp1 <- tmp1[tmp1$end == max(tmp1$end),]
                tmp1 <- tmp1[tmp1$lp == min(tmp1$lp),]
                tmp1$start <- as.integer(tmp1$end + 1)
                tmp1$end <- as.integer(tmp1$end + five_prim_utr_length_cor)
                tmp1$source <- 'JBIO'
                tmp1$annotationType <- 'five_prim_UTR'
                tmp1$score <- '.'
                tmp1$phase <- '.'
                tmp2 <- tmp2[tmp2$start == min(tmp2$start),]
                tmp2 <- tmp2[tmp2$lp == min(tmp2$lp),]
                tmp2$source <- 'JBIO'
                tmp2$annotationType <- 'three_prim_UTR'
                tmp2$score <- '.'
                tmp2$phase <- '.'
                tmp2$end <- as.integer(tmp2$start - 1)
                tmp2$start <- as.integer(tmp2$start - three_prim_utr_length_cor)
              }
              
              final <- rbind(tmp0, tmp1, tmp2)
              final <- final[,!colnames(final) %in% 'lp']
              df <- rbind(df, final)
              rm(final)
              
            } else if (min_val > max_val + (five_prim_utr_length + three_prim_utr_length)/3 & i > 1 & i < length(gen_list)) {
              tmp1 <- tmp1[tmp1$end < min(tmp2$start),]
              tmp_primary <- tmp_primary[!tmp_primary$lp %in% tmp1$lp,]

              
              if (min(tmp2$start) > (max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)+2)) {
                five_prim_utr_length_cor <- round(five_prim_utr_length,0)
                three_prim_utr_length_cor <- round(three_prim_utr_length,0)
              } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/2) {
                five_prim_utr_length_cor <- round(five_prim_utr_length*0.5,0)
                three_prim_utr_length_cor <- round(three_prim_utr_length*0.5,0)
              } else if (min(tmp2$start) > max(tmp1$end) + (five_prim_utr_length + three_prim_utr_length)/3) {
                five_prim_utr_length_cor <- round(five_prim_utr_length*0.3,0)
                three_prim_utr_length_cor <- round(three_prim_utr_length*0.3,0)
              }
              
              if (sen == '+') {
                
                tmp1 <- tmp1[tmp1$end == max(tmp1$end),]
                tmp1 <- tmp1[tmp1$lp == min(tmp1$lp),]
                tmp1$start <- as.integer(tmp1$end + 1)
                tmp1$end <- as.integer(tmp1$end + three_prim_utr_length_cor)
                tmp1$source <- 'JBIO'
                tmp1$annotationType <- 'three_prim_UTR'
                tmp1$score <- '.'
                tmp1$phase <- '.'
                tmp2 <- tmp2[tmp2$start == min(tmp2$start),]
                tmp2 <- tmp2[tmp2$lp == min(tmp2$lp),]
                tmp2$source <- 'JBIO'
                tmp2$annotationType <- 'five_prim_UTR'
                tmp2$score <- '.'
                tmp2$phase <- '.'
                tmp2$end <- as.integer(tmp2$start - 1)
                tmp2$start <- as.integer(tmp2$start - five_prim_utr_length_cor)
              
              } else if (sen == '-') {
                
                tmp1 <- tmp1[tmp1$end == max(tmp1$end),]
                tmp1 <- tmp1[tmp1$lp == min(tmp1$lp),]
                tmp1$start <- as.integer(tmp1$end + 1)
                tmp1$end <- as.integer(tmp1$end + five_prim_utr_length_cor)
                tmp1$source <- 'JBIO'
                tmp1$annotationType <- 'five_prim_UTR'
                tmp1$score <- '.'
                tmp1$phase <- '.'
                tmp2 <- tmp2[tmp2$start == min(tmp2$start),]
                tmp2 <- tmp2[tmp2$lp == min(tmp2$lp),]
                tmp2$source <- 'JBIO'
                tmp2$annotationType <- 'three_prim_UTR'
                tmp2$score <- '.'
                tmp2$phase <- '.'
                tmp2$end <- as.integer(tmp2$start - 1)
                tmp2$start <- as.integer(tmp2$start - three_prim_utr_length_cor)
                 
              }
              
              final <- rbind(tmp1, tmp2)
              final <- final[,!colnames(final) %in% 'lp']
              df <- rbind(df, final)
              rm(final)
              
            } else if (i == length(gen_list) & i != 1 & i != 0) {
              
               if (sen == '+') {
                
                tmp0 <- tmp1[tmp1$gene_name %in% gen_list[i],]
                tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
                tmp0 <- tmp0[tmp0$lp == min(tmp0$lp),]
                tmp0$start <- as.integer(tmp0$end + 1)
                tmp0$end <- as.integer(tmp0$start + three_prim_utr_length)
                tmp0$source <- 'JBIO'
                tmp0$annotationType <- 'three_prim_UTR'
                tmp0$score <- '.'
                tmp0$phase <- '.'
              
               } else if (sen == '-') {
                 
                 tmp0 <- tmp1[tmp1$gene_name %in% gen_list[i],]
                 tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
                 tmp0 <- tmp0[tmp0$lp == min(tmp0$lp),]
                 tmp0$start <- as.integer(tmp0$end + 1)
                 tmp0$end <- as.integer(tmp0$start + five_prim_utr_length)
                 tmp0$source <- 'JBIO'
                 tmp0$annotationType <- 'fivee_prim_UTR'
                 tmp0$score <- '.'
                 tmp0$phase <- '.'
                
               }
              
              final <- tmp0
              final <- final[,!colnames(final) %in% 'lp']
              df <- rbind(df, final)
              rm(final)
            
          }
        } 
      }
    }
}

  sum_genes <- as.data.frame(summary(as.factor(df$gene_name), maxsum = length(unique(df$gene_name))))
  colnames(sum_genes) <- 'n'
  sum_genes$gene_name <- rownames(sum_genes)
  sum_genes <- sum_genes$gene_name[sum_genes$n  < 10]
  
  df <- df[df$gene_name %in% sum_genes,]
  
  output <- rbind(input, df)
  output <- output[order(output$sort_val, output$start),]
  output <- output[,!colnames(output) %in% 'sort_val']
  output <- distinct(output)
  
  return(output)
}



prepare_to_reffflat <- function(input) {
  
  cat('\n\n Transcripts repair & extending...             \n\n')
  
  input$sort_val <- 1:length(input$chr)
  chromosomes <- unique(input$chr)
  
  CDS <- input[toupper(input$annotationType) %in% c('GENE', 'EXON', 'CDS', 'TRANSCRIPT', 'MRNA', 'FIVE_PRIM_UTR', 'THREE_PRIM_UTR'), ]
  CDS$annotationType[toupper(CDS$annotationType) %in% c('MRNA', 'GENE', 'TRANSCRIPT')] <- 'transcript'
  CDS$diff <- CDS$end - CDS$start
  
  df <- data.frame(matrix(ncol = length(colnames(CDS)), nrow = 0))
  colnames(df) <- colnames(CDS)
  
  for (chr in chromosomes) {
    
    tmp_primary <- CDS[CDS$chr %in% chr,]
    tmp_primary <- tmp_primary[order(tmp_primary$start),]
    tmp_primary$lp <- 1:length(tmp_primary$chr)
    
    
    for (sen in c('+','-')) {
      tmp <- tmp_primary[tmp_primary$strand %in% sen,]
      gen_list <- c()
      check <- 0
      for (gen in tmp$gene_name) {
        if (gen != check) {
          gen_list <- c(gen_list, gen)
          check <- gen
        }
        
      }
      
      if (length(gen_list) > 0) {
        
      for (i in 1:length(gen_list)) {
        
        cat('\r',paste('LOC:',chr, '| STRAND:', sen , '| PROGRESS:', round((i/length(gen_list)*100),2), '%           '))
       
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            tmp_primary <- tmp_primary[!tmp_primary$lp %in% tmp1$lp,]
            tmp2 <- tmp1[tmp1$annotationType == 'transcript',]
            
            if (!'transcript' %in% tmp1$annotationType) {
              if ('exon' %in% tmp1$annotationType) {
              tmp_t <- tmp1[tmp1$annotationType == 'exon',]
              } else if ('CDS' %in% tmp1$annotationType) {
                
              tmp_t <- tmp1[tmp1$annotationType == 'CDS',]
              tmp_exon_tmp <- tmp1[tmp1$annotationType == 'CDS',]
              tmp_exon_tmp <- distinct(tmp_exon_tmp)
                
                for (cds in 1:length(tmp_exon$chr)) {
                tmp_exon <- tmp_exon_tmp[cds,]
                tmp_exon$source <- 'JBIO'
                tmp_exon$annotationType <- 'exon'
                tmp_exon$score <- '.'
                tmp_exon$phase <- '.'
                tmp_exon <- tmp_exon[,!colnames(tmp_exon) %in% c('lp', 'diff')]
                df <- rbind(df, tmp_exon)
                }
              }
              
            tmp2 <- tmp_t[tmp_t$start == min(tmp_t$start),][1,]
            tmp2$end <- tmp_t$end[tmp_t$end == max(tmp_t$end)][1]
            tmp2$source <- 'JBIO'
            tmp2$annotationType <- 'transcript'
            tmp2$score <- '.'
            tmp2$phase <- '.'
            rm(tmp_t)
              
            }
            
            
            tmp2 <- tmp2[tmp2$diff == max(tmp2$diff),] 
            tmp2 <- tmp2[1,]
            
            if (TRUE %in% grepl('UTR', tmp1$annotationType)) {
              
              if (sen == '+') {
              
              if ('five_prim_UTR' %in% tmp1$annotationType) {
                  tmp_utr <- tmp1[tmp1$annotationType == 'five_prim_UTR',]
               
              
              if (tmp2$start > tmp_utr$start[tmp_utr$start == min(tmp_utr$start)]) {
                  tmp2$start <- as.integer(tmp_utr$start[tmp_utr$start == min(tmp_utr$start)])
                }
                  rm(tmp_utr)
                }
              
              if ('three_prim_UTR' %in% tmp1$annotationType) {
                  tmp_utr <- tmp1[tmp1$annotationType == 'three_prim_UTR',]
              
                
              if (tmp2$end < tmp_utr$end[tmp_utr$end == max(tmp_utr$end)]) {
                  tmp2$end <- as.integer(tmp_utr$end[tmp_utr$end == max(tmp_utr$end)])
                }
                rm(tmp_utr)
                } 
              
            } else  if (sen == '-') {
              
              if ('three_prim_UTR' %in% tmp1$annotationType) {
                tmp_utr <- tmp1[tmp1$annotationType == 'three_prim_UTR',]
                
                
              if (tmp2$start > tmp_utr$start[tmp_utr$start == min(tmp_utr$start)]) {
                  tmp2$start <- as.integer(tmp_utr$start[tmp_utr$start == min(tmp_utr$start)])
                }
                rm(tmp_utr)
                }
              
              if ('five_prim_UTR' %in% tmp1$annotationType) {
                tmp_utr <- tmp1[tmp1$annotationType == 'five_prim_UTR',]
                
              
              if (tmp2$end < tmp_utr$end[tmp_utr$end == max(tmp_utr$end)]) {
                  tmp2$end <- as.integer(tmp_utr$end[tmp_utr$end == max(tmp_utr$end)])
                }
                rm(tmp_utr)
                } 
              
              
            }
            
            
      
            tmp2$source <- 'JBIO'
            tmp2$annotationType <- 'transcript'
            tmp2$score <- '.'
            tmp2$phase <- '.'
            
           
            
            final_transcript <- tmp2[,!colnames(tmp2) %in% c('lp', 'diff')]
 
            
            df <- rbind(df, final_transcript)
            rm(tmp1, tmp2, final_transcript)
            
            }
          }
        }
      } 
    }
  
  
  sum_genes <- as.data.frame(summary(as.factor(df$gene_name), maxsum = length(unique(df$gene_name))))
  colnames(sum_genes) <- 'n'
  sum_genes$gene_name <- rownames(sum_genes)
  sum_genes <- sum_genes$gene_name[sum_genes$n  < 10]
  
  df <- df[df$gene_name %in% sum_genes,]
  
  input <- input[toupper(input$annotationType) %in% c('GENE', "EXON", "CDS", 'TRANSCRIPT', 'MRNA'), ]
  input$annotationType[toupper(input$annotationType) %in% c('MRNA', 'GENE', 'TRANSCRIPT')] <- 'transcript'
  output <- rbind(input, df)
  output <- output[order(output$sort_val, output$start),]
  output <- output[,!colnames(output) %in% 'sort_val']
  output <- distinct(output)
  
  return(output)
}





refflat_create <- function(input) {
  
  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  CPU <- detectCores() - 2
  
  cl <- makeCluster(CPU)

  
  
  registerDoParallel(cl)
  registerDoSNOW(cl)
  
  
  
  
  input <- input[!duplicated(input[,c('chr', 'start', 'end' ,'gene_name', 'strand', 'annotationType')]),]
  input$sort_val <- 1:length(input$chr)
  input$diff <- input$end - input$start
  
  chromosomes <- unique(input$chr)

  
  df <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(df) <- c('geneName', 'name', 'chrom', 'strand', 'txStart','txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds')
  
  cat('\n\n Refflat creating ... \n\n')
  results <- foreach(chr = chromosomes, .packages = c('dplyr'), .options.snow = opts, .combine=rbind) %dopar% {
    
    tmp_primary <- input[input$chr %in% chr,]
    tmp_primary <- tmp_primary[order(tmp_primary$start),]
    tmp_primary$lp <- 1:length(tmp_primary$chr)
    
    
    for (sen in c('+','-')) {
      tmp <- tmp_primary[tmp_primary$strand %in% sen,]
      gen_list <- c()
      check <- 0
      for (gen in tmp$gene_name) {
        if (gen != check) {
          gen_list <- c(gen_list, gen)
          check <- gen
        }
        
      }
      
      
      if (length(gen_list) == 1) {break}
      
      for (i in 1:length(gen_list)) {
        

          check_vector <- tmp$annotationType[tmp$gene_name %in% gen_list[i]]
          if ('exon' %in% check_vector & 'transcript' %in% check_vector) {
          tmp_tx <- tmp[tmp$gene_name %in% gen_list[i],]
          tmp_tx <- tmp_tx[tmp_tx$annotationType == 'transcript',]
          tmp_dif <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff),]
          min_tr <- min(tmp_dif$start)
          max_tr <- max(tmp_dif$end)
          tmp_tx <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff) | tmp_tx$start < min_tr | tmp_tx$end > max_tr ,]
          tmp_ex <- tmp[tmp$gene_name %in% gen_list[i] & tmp$annotationType == 'exon',]

          for (trancript in 1:nrow(tmp_tx)) {
            tmp_tx_tmp <- tmp_tx[trancript,]
            tmp_ex_tmp <- tmp_ex[tmp_ex$start >= tmp_tx_tmp$start & tmp_ex$end <= tmp_tx_tmp$end,]
            
            decission <- TRUE
            trn <- 0
            while (decission) {
              trn <- trn + 1
                dec <- c()
                start <- c()
                end <- c()
                before <- 0
                for(ex in 1:nrow(tmp_ex_tmp)) {
                  if (before <= tmp_ex_tmp$start[ex]) {
                    start <- c(start, as.integer(tmp_ex_tmp$start[ex]))
                    end <- c(end, as.integer(tmp_ex_tmp$end[ex]))
                    before  <-  as.integer(tmp_ex_tmp$end[ex])
                    
                   
      
                  } else if (before > as.integer(tmp_ex_tmp$start[ex])) {
                    
                    dec <- c(dec, ex-1)
                   
                  }
                }
                
                if (length(dec) > 0) {
                tmp_ex_tmp <- tmp_ex_tmp[-dec[1], ]
                decission <- TRUE
                
                } else if (length(dec) == 0) {
                  decission <- FALSE
                  
                }
                
                
                df[nrow(df) + 1,] <- c(as.character(gsub('"', '', tmp_tx_tmp$gene_name[1])), as.character(paste0(gsub('"', '', tmp_ex_tmp$gene_name[1]),'.',trn)), as.character(tmp_ex_tmp$chr[1]), as.character(tmp_ex_tmp$strand[1]), as.integer(tmp_tx_tmp$start[1]), as.integer(tmp_tx_tmp$end[1]), as.integer(min(start)), as.integer(max(end)), as.integer(length(start)), sub('"', '', paste0(start, collapse = ',')), sub('"', '', paste0(end, collapse = ',')))
                df <- distinct(df)
                
                
                } 
                  
                 
                
                }
              
            
            
                }
          
          
              }
            }
                results <- df
        }
  
      
      close(pb)
      stopCluster(cl)  
      return(results)
} 





###############

create_full_GTF <- function(input) {
  
  
  output <- input[,1:8]
  gene_id <- paste0('gene_id ',gsub(' ', '',input$gene_id), ';')
  gene_name <-  paste0('gene_name ', gsub(' ','',input$gene_name), ';')
  transcript_name <- paste0('transcript_name ', gsub(' ', '',input$transcript_name), ';')
  transcript_id <- paste0('transcript_id ', gsub(' ', '',input$transcript_id), ';')
  output$combine <- paste0(gene_id, gene_name, transcript_name, transcript_id)
  
  
  return(output)
  
}




create_reduced_GTF <- function(input) {
  
  output <- input %>%
    dplyr::select(chr,start,end,strand,transcript_name,transcript_id, gene_name, gene_id)
  output$annotationType <- input$annotationType
  output$transcriptType <- NA
  
  
  return(output)  
  
}

