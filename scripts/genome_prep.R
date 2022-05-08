args <- commandArgs()


#Paths and arguments from env
{
  print(args)
  species  <- args[6]
  annotation <-args[7]
  output <- args[8]
  
  functions <- file.path(getwd(), 'scripts/gtf_tool.R')
  source(functions, local = T)
}

#Configuration file 

{
  conf_file <- read.csv(file = file.path(getwd(), 'requirements_file/genome.conf'), header = F, sep = '=', row.names = 1)
  
  human_extend <- as.character(conf_file$V2[grep(pattern = 'human_extend', rownames(conf_file))])
  
  mice_extend <- as.character(conf_file$V2[grep(pattern = 'mice_extend', rownames(conf_file))])
  
  mix_extend <- as.character(conf_file$V2[grep(pattern = 'mix_extend', rownames(conf_file))])
  
  custom_extend <- as.character(conf_file$V2[grep(pattern = 'custom_extend', rownames(conf_file))])
  
  three_prime_utr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'three_prime_utr', rownames(conf_file))]))
  
  five_prime_utr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'five_prim_utr', rownames(conf_file))]))
 
}



GTF <- load_annotation(annotation)

GTF2 <- create_GTF_df(GTF)

if (species == 'human'&  grepl('T', toupper(human_extend))) {
GTF2 <- add_UTR(GTF2, five_prime_utr, three_prime_utr)
} else if (species == 'mice'&  grepl('T', toupper(mice_extend))) {
  GTF2 <- add_UTR(GTF2, five_prime_utr, three_prime_utr)
} else if (species == 'custom'&  grepl('T', toupper(custom_extend))) {
  GTF2 <- add_UTR(GTF2, five_prime_utr, three_prime_utr)
} else if (species == 'mix' &  grepl('T', toupper(human_extend))) {
  GTF2 <- add_UTR(GTF2, five_prime_utr, three_prime_utr)
} 


GTF3 <- prepare_to_reffflat(GTF2)

GTF4 <- refflat_create(GTF3)

write.table(GTF4, file.path(output, 'correct_annotation.refflat'), quote = F, sep = '\t', col.names = F, row.names = F)


GTF5 <- create_full_GTF(GTF3)

write.table(GTF5, file.path(output, 'correct_annotation.gtf'), quote = F, sep = '\t', col.names = F, row.names = F)

GTF6 <- create_reduced_GTF(GTF3)


write.table(GTF6, file.path(output, 'reduced_annotation.gtf'), quote = F, sep = '\t', col.names = T, row.names = F)



