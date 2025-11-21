args <- commandArgs()
set.seed(123)
options(scipen = 999)


#Paths and arguments from env
{
  print(args)
  species  <- args[6]
  annotation <-args[7]
  output <- args[8]

}

#Configuration file 

{
  conf_file <- read.csv(
    file = file.path(getwd(), 'requirements_file/genome.conf'),
    header = FALSE,
    sep = '=',
    row.names = 1,
    comment.char = "#"
  )
  
  extend <- as.logical(conf_file$V2[grep(pattern = 'extend', rownames(conf_file))])
  
  optimize <- as.logical(conf_file$V2[grep(pattern = 'optimize_names', rownames(conf_file))])
  
  sep_factor <- as.numeric(as.character(conf_file$V2[grep(pattern = 'sep_factor', rownames(conf_file))]))
  
  three_prime_utr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'three_prime_utr', rownames(conf_file))]))
  
  five_prime_utr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'five_prim_utr', rownames(conf_file))]))
 
}



GTF <- GTF.tool::load_annotation(annotation)

if (extend) {
  
GTF2 <- GTF.tool::create_GTF_df(GTF, optimize = T, shift = sep_factor)

} else {
  
  GTF2 <- GTF.tool::create_GTF_df(GTF, optimize = optimize, shift = sep_factor)
  
}


GTF3 <- GTF.tool::add_CDS(GTF2)


if (extend) {
  
  if ((sum(toupper(GTF3$annotationType) %in% c('TRANSCRIPT', 'CDS', 'MRNA') & GTF3$gene_type == "protein_coding")) / length(toupper(GTF3$annotationType) %in% c('TRANSCRIPT', 'CDS', 'MRNA')) < 0.20) {
    GTF3 <- GTF.tool::add_UTR(GTF3, 
                              five_prime_utr = five_prime_utr , 
                              three_prime_utr = three_prime_utr, 
                              biotype = NULL, transcript_limit = 300
                              )
    
  } else {
    GTF3 <- GTF.tool::add_UTR(GTF3, 
                              five_prime_utr = five_prime_utr , 
                              three_prime_utr = three_prime_utr,
                              biotype = 'protein_coding')
    
  }
  
  
  GTF4 <- GTF.tool::create_full_GTF(GTF3)
  
} else {
  
  GTF4 <- GTF.tool::create_full_GTF(GTF3)
  
}


write.table(GTF4, file.path(output, 'correct_annotation.gtf'), quote = F, sep = '\t', col.names = F, row.names = F)

GTF5 <- GTF.tool::create_reduced_GTF(GTF3)

write.table(GTF5, file.path(output, 'reduced_annotation.gtf'), quote = F, sep = '\t', col.names = T, row.names = F)

GTF6 <- GTF.tool::refflat_create(GTF3)

write.table(GTF6, file.path(output, 'correct_annotation.refflat'), quote = F, sep = '\t', col.names = F, row.names = F)
