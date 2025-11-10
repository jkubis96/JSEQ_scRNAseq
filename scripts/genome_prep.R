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
  
  sep_factor <- as.numeric(conf_file$V2[grep(pattern = 'sep_factor', rownames(conf_file))])
  
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

GTF4 <- GTF.tool::add_introns(GTF3)

if (extend) {
  
  GTF5 <- GTF.tool::add_UTR(GTF4, five_prime_utr = five_prime_utr , three_prime_utr = three_prime_utr)
  
  GTF6 <- GTF.tool::create_full_GTF(GTF5)
  
} else {
  
  GTF6 <- GTF.tool::create_full_GTF(GTF4)
  
}


write.table(GTF6, file.path(output, 'correct_annotation.gtf'), quote = F, sep = '\t', col.names = F, row.names = F)

GTF7 <- GTF.tool::create_reduced_GTF(GTF5)

write.table(GTF7, file.path(output, 'reduced_annotation.gtf'), quote = F, sep = '\t', col.names = T, row.names = F)

GTF8 <- GTF.tool::refflat_create(GTF5)

write.table(GTF8, file.path(output, 'correct_annotation.refflat'), quote = F, sep = '\t', col.names = F, row.names = F)
