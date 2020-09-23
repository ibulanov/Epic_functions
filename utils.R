library(RMySQL)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

scm_mm10 <- dbConnect(MySQL(), user='i747v', password='igor',
                          dbname="scm_mm10", host='zuse1db')

scm_mm10
#convert string ranges to sql query
input_string <- args[1]
spl_str <- strsplit(input_string, split=',', fixed=FALSE)
colon <- strsplit(spl_str[[1]], split = ":", fixed = TRUE)
result_seq <- c()
for (ee in colon) {
  if (length(ee) > 1){
    tmp_v <- seq(ee[1], ee[2])
    result_seq <- c(result_seq, tmp_v)
  } else {
    tmp_v <- as.integer(ee)
    result_seq <- c(result_seq, tmp_v)
  }
}
result_seq <- toString(result_seq)

#test <- dbSendQuery(mouse_public, sprintf("select expID, chrom, pos
#                    from cpgCoverage
#                    where expID in (%s)", result_seq) %>% fetch(n=-1)

cat(sprintf("select expID, chrom, pos from cpgCoverage where expID in (%s)", result_seq))


#to test taking data based on chr
output_db <- dbSendQuery(mouse_public, "select chrom, chromStart, chromEnd
                         from cpgIslandExt
                         where chrom = 'chr1' ") %>% fetch(., n=-1)

#preparing scNMT data for MySQL db
list_of_scnmt_cells <- list()
for (i in seq_along(sample_metadata$id_met)) {
  cur_cell <- fread(sprintf("/home/igor/Data/scnmt_gastrulation/met/cpg_level/%s.tsv.gz", sample_metadata$id_met[i]))
  cur_cell$expID <- seq(400001,401000)[i]
  cur_cell$coverage <- cur_cell$met_reads + cur_cell$nonmet_reads
  cur_cell$nonmet_reads <- NULL
  cur_cell <- cur_cell[,c(5,1,2,6,3,4)]
  colnames(cur_cell) <- c("expID", "chrom", "pos", "coverage", "mC", "ratio")
  list_of_scnmt_cells[[i]] <- cur_cell
  names(list_of_scnmt_cells)[i] <- sample_metadata$id_met[i]
  print(i)
}

#replace chr from numeric to string chrX
for (i in seq_along(list_of_scnmt_cells)) {
  list_of_scnmt_cells[[i]]$chrom <- paste0("chr", list_of_scnmt_cells[[i]]$chrom)
} 



#save the results
for (i in seq_along(input_cpgs_chr1)) {
  fwrite(input_cpgs_chr1[[i]], 
         sprintf("/home/igor/Data/melissa/inhouse_chr1_plus_scnmt_77/%s.tsv.gz", names(input_cpgs_chr1)[i]),
         sep = '\t', col.names = FALSE)
  print(i)
}

#scnmt_table <- rbindlist(list_of_scnmt_cells)

scnmt <- fread("/home/igor/Data/melissa/scnmt_fordb/scnmt_986.tsv.gz")

protein_coding_genes <- fread("/home/igor/Data/scnmt_gastrulation/met/feature_level/met_dat_anno.tsv")

na_list_cells <- c()
#to check whether cells with full NA regions
for (i in seq_along(input_melissa$met)) {
  na_list_cells[i] <- any(!is.na(input_melissa$met[[i]]))
}

no_na_list_cells <- c()
#to check whether cells with full NA regions
for (i in seq_along(input_melissa$met)) {
  no_na_list_cells[i] <- any(is.na(input_melissa$met[[i]]))
}

length_bpr <- c()
for (i in seq_along(input_melissa$met)) {
  length_bpr[i] <- sum(!is.na(input_melissa$met[[i]]))
}
