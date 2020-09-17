####script to find overlapped features across cells####

library(RMySQL)
library(dplyr)
library(data.table)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Please specify the first argument as desired feature and second argument as feature size\n 
(if feature size is not applicable (e.g. cpg_islands), just please specify any number \n
 it will be used as counter for output filename)\n", call.=FALSE)
}

####generate the db####
zuse1db = dbConnect(MySQL(), user='i747v', password='igor',
                    dbname='mouse_public', host='zuse1db')

#function for processing input cpg bases data
cpgs_proc <- function(input_cpgs) {
  
  #to get the input data (cpgBases)
  input_cpgs <-  dbSendQuery(zuse1db, "select C.expID, B.chrom, B.pos, B.strand
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        where C.expID between 100310 and 100386") %>% fetch(n=-1) 
  
  input_cpgs <- setDT(input_cpgs) %>% .[, compl_pos := ifelse(strand == 1, pos+1, pos-1)]
  
  #create the complementary df
  cpgs_compl <- data.table(input_cpgs$expID, input_cpgs$chrom, input_cpgs$compl_pos)
  
  #remove extra columns
  input_cpgs <- input_cpgs[,c("strand","compl_pos"):= NULL] 
  
  colnames(cpgs_compl) <- colnames(input_cpgs)
  
  #concatenate the two dataframes and order the rows
  input_cpgs <- rbind(input_cpgs, cpgs_compl) %>% arrange(expID) %>% distinct() 
  
  input_cpgs$pos_d <- input_cpgs$pos
  
  #delete cpgbases complementary data table
  cpgs_compl <- NULL
  
  cat(sprintf("Cpg bases positions for cells 100310:100386 were generated\n"))
  
  
  
}

input_cpgs <- cpgs_proc()
fwrite(input_cpgs, "input_cpgs.tsv")
cat("The input cpg bases file was saved at the same directory - specify path to file as third argument \n 
for the next script running for faster implementation. It will be used the local input data instead of sql query ")

####function for finding the intersection windows from the raw data####
windows_proc <- function(input_cpgs, winsize) {
  
  input_cpgs[, start := pos - pos %% winsize] %>% .[, finish := start + (winsize - 1)] %>%
    .[, c("pos") := NULL] %>% .[, feature_coord := paste(chrom, ":",start, ":", finish)] %>%
    .[, c("chrom", "start", "finish") := NULL] %>% .[, feature_coord := gsub(" ", "", feature_coord)] %>% 
    distinct() %>% setkey(expID)
  
  unique_features_number <- nrow(input_cpgs)
  input_cpgs <- split(input_cpgs, input_cpgs$expID)
  
  cat(sprintf("For each cpg site the consequent feature was found\n"))
  
  working_set1 <- input_cpgs
  working_set2 <- input_cpgs
  res_matrix <- matrix(nrow = length(input_cpgs), ncol = length(input_cpgs))
  for (i in seq_along(working_set1)) {
    for (j in seq_along(working_set1)) {
      res_two_ds <- intersect(working_set1[[i]]$feature_coord, working_set1[[j]]$feature_coord)
      res_matrix[i,j] <- length(res_two_ds)  
    }
  }
  
  res_matrix_pct <- res_matrix / unique_features_number
  colnames(res_matrix_pct) <- names(input_cpgs)
  rownames(res_matrix_pct) <- names(input_cpgs)
  
  return(res_matrix_pct)

}

if (args[1] == "window") {
  matrix_result <- windows_proc(input_cpgs, winsize = as.numeric(args[2]))
  cat(sprintf("Windows coordinates were generated\n"))
}

#process promoters
prom_db_proc <- function(prom_length) {
  #take data for promoters
  output_db <- dbSendQuery(zuse1db, "select chrom, strand, txStart, txEnd
                         from knownGene") %>% fetch(n=-1) %>% setDT() %>% 
    .[, feature_start := ifelse(strand == "+", txStart - prom_length, txEnd+1)] %>%
    .[, feature_end := ifelse(strand == "+", txStart-1, txEnd+1 + prom_length)] %>%
    .[, c("strand", "txStart", "txEnd") := NULL] %>% distinct(.)

  unique_features_number <- nrow(output_db)

  #set keys for data.table class
  setkey(input_cpgs, chrom, pos, pos_d)
  setkey(output_db, chrom, feature_start, feature_end)
  
  #find cpg sites positions in promoters, create concatenated column and delete extra columns
  cells_features <- foverlaps(input_cpgs, output_db, nomatch = NULL) %>%
    .[,c("pos","pos_d"):= NULL] %>% setkey(expID) %>% distinct() %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)] %>%
    .[,c("chrom","feature_start", "feature_end"):= NULL]
  
  cells_features <- split(cells_features, cells_features$expID)
  
  cat(sprintf("For each cpg site the consequent feature was found\n"))
  
  working_set1 <- cells_features
  working_set2 <- cells_features
  res_matrix <- matrix(nrow = length(cells_features), ncol = length(cells_features))
  for (i in seq_along(working_set1)) {
    for (j in seq_along(working_set1)) {
      res_two_ds <- intersect(working_set1[[i]]$feature_coord, working_set1[[j]]$feature_coord)
      res_matrix[i,j] <- length(res_two_ds)  
    }
  }
  
  res_matrix_pct <- res_matrix / unique_features_number
  colnames(res_matrix_pct) <- names(cells_features)
  rownames(res_matrix_pct) <- names(cells_features)
  
  return(res_matrix_pct)
  
  
}

if (args[1] == "promoter") {
  matrix_result <- prom_db_proc(prom_length = as.numeric(args[2]))
  cat(sprintf("Promoters coordinates were generated\n"))
}

#function for processing the cpg islands table
islands_proc <- function() {
  #take data for cpg islands
  output_db <- dbSendQuery(zuse1db, "select chrom, chromStart, chromEnd
                         from cpgIslandExt") %>% fetch(., n=-1) %>% setDT() %>%
    setnames(c("chromStart","chromEnd"), c("feature_start", "feature_end"))
  
  unique_features_number <- nrow(output_db)
  
  #set keys for data.table class
  setkey(input_cpgs, chrom, pos, pos_d)
  setkey(output_db, chrom, feature_start, feature_end)
  
  #find cpg sites positions in promoters, create concatenated column and delete extra columns
  cells_features <- foverlaps(input_cpgs, output_db, nomatch = NULL) %>%
    .[,c("pos","pos_d"):= NULL] %>% setkey(expID) %>% distinct() %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)] %>%
    .[,c("chrom","feature_start", "feature_end"):= NULL]
  
  cells_features <- split(cells_features, cells_features$expID)
  
  cat(sprintf("For each cpg site the consequent feature was found\n"))
  
  working_set1 <- cells_features
  working_set2 <- cells_features
  res_matrix <- matrix(nrow = length(cells_features), ncol = length(cells_features))
  for (i in seq_along(working_set1)) {
    for (j in seq_along(working_set1)) {
      res_two_ds <- intersect(working_set1[[i]]$feature_coord, working_set1[[j]]$feature_coord)
      res_matrix[i,j] <- length(res_two_ds)  
    }
  }
  
  res_matrix_pct <- res_matrix / unique_features_number
  colnames(res_matrix_pct) <- names(cells_features)
  rownames(res_matrix_pct) <- names(cells_features)
  
  return(res_matrix_pct)
  
} 

if (args[1] == "cpg_islands") {
  matrix_result <- islands_proc()
  cat(sprintf("Cpg islands coordinates were generated\n"))
}

#function for processing exons
exons_proc <- function() {
  output_db <- dbSendQuery(zuse1db, "select chrom, exonStarts, exonEnds
                         from knownGene") %>% fetch(n=-1) %>% 
    separate_rows(exonStarts, exonEnds, sep=",", convert = TRUE) %>%
    na.omit() %>% as.data.table() %>%
    setnames(c("exonStarts","exonEnds"), c("feature_start", "feature_end"))
  
  unique_features_number <- nrow(output_db)
  
  #set keys for data.table class
  setkey(input_cpgs, chrom, pos, pos_d)
  setkey(output_db, chrom, feature_start, feature_end)
  
  #find cpg sites positions in promoters, create concatenated column and delete extra columns
  cells_features <- foverlaps(input_cpgs, output_db, nomatch = NULL) %>%
    .[,c("pos","pos_d"):= NULL] %>% setkey(expID) %>% distinct() %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)] %>%
    .[,c("chrom","feature_start", "feature_end"):= NULL]
  
  cells_features <- split(cells_features, cells_features$expID)
  
  cat(sprintf("For each cpg site the consequent feature was found\n"))
  
  working_set1 <- cells_features
  working_set2 <- cells_features
  res_matrix <- matrix(nrow = length(cells_features), ncol = length(cells_features))
  for (i in seq_along(working_set1)) {
    for (j in seq_along(working_set1)) {
      res_two_ds <- intersect(working_set1[[i]]$feature_coord, working_set1[[j]]$feature_coord)
      res_matrix[i,j] <- length(res_two_ds)  
    }
  }
  
  res_matrix_pct <- res_matrix / unique_features_number
  colnames(res_matrix_pct) <- names(cells_features)
  rownames(res_matrix_pct) <- names(cells_features)
  
  return(res_matrix_pct)
 
  
}

if (args[1] == "exons") {
  matrix_result <- exons_proc()
  cat(sprintf("Cpg islands coordinates were generated\n"))
}

introns_proc <- function() {
  
  output_db <- fread("/home/igor/Data/common_matrices/common_exons_introns/introns_psql.tsv", sep='\t')
  
  unique_features_number <- nrow(output_db)
  
  #set keys for data.table class
  setkey(input_cpgs, chrom, pos, pos_d)
  setkey(output_db, chrom, feature_start, feature_end)
  
  #find cpg sites positions in promoters, create concatenated column and delete extra columns
  cells_features <- foverlaps(input_cpgs, output_db, nomatch = NULL) %>%
    .[,c("pos","pos_d"):= NULL] %>% setkey(expID) %>% distinct() %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)] %>%
    .[,c("chrom","feature_start", "feature_end"):= NULL]
  
  cells_features <- split(cells_features, cells_features$expID)
  
  cat(sprintf("For each cpg site the consequent feature was found\n"))
  
  working_set1 <- cells_features
  working_set2 <- cells_features
  res_matrix <- matrix(nrow = length(cells_features), ncol = length(cells_features))
  for (i in seq_along(working_set1)) {
    for (j in seq_along(working_set1)) {
      res_two_ds <- intersect(working_set1[[i]]$feature_coord, working_set1[[j]]$feature_coord)
      res_matrix[i,j] <- length(res_two_ds)  
    }
  }
  
  res_matrix_pct <- res_matrix / unique_features_number
  colnames(res_matrix_pct) <- names(cells_features)
  rownames(res_matrix_pct) <- names(cells_features)
  
  return(res_matrix_pct)
  
}

if (args[1] == "intron") {
  matrix_result <- intron_proc()
  cat(sprintf("Exons coordinates were generated\n"))
}

cat(sprintf("The matrix of overlapped genomic features across all of the cells was contructed\n"))

#mean
mean_value <- mean(matrix_result[upper.tri(matrix_result)]) %>% round(4)

cat(sprintf("Mean value of overlapped features across all of the cells was calculated: %s\n", mean_value))

fwrite(matrix_result, sprintf("overlapping_mtx_%s_%s.tsv", args[1], args[2]))

cat(sprintf("The result was saved as tsv table: overlapping_mtx_prom_%s.rds", args[2]))


