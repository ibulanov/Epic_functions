#20/01/21 - process epidermis and dermis cells 

library(RMySQL)
library(dplyr)
library(data.table)

dbconn <- dbConnect(MySQL(), user='igor', password='igor123',
                    dbname="scm_ih_hg19", host='a130-pc-bioinfo-09')

dbListTables(dbconn)

#epidermis
seq = c(15001:15401) %>% toString()
#dermis
seq2 = c(15402:15963) %>% toString()


epidermis <- dbSendQuery(dbconn, sprintf("select C.expID, C.chrom, C.pos, B.strand, C.coverage, C.mC, C.ratio 
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        where C.expID in (%s)", seq)) %>% fetch(n=-1)

dermis <- dbSendQuery(dbconn, sprintf("select C.expID, C.chrom, C.pos, B.strand, C.coverage, C.mC, C.ratio 
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        where C.expID in (%s)", seq2)) %>% fetch(n=-1)
#back up
epi_r <- epidermis
der_r <- dermis

#normalize
epidermis$pos <- ifelse(epidermis$strand == 1, epidermis$pos, epidermis$pos-1)
dermis$pos <- ifelse(dermis$strand == 1, dermis$pos, dermis$pos-1)

#mutate
epidermis$ratio <- ifelse(epidermis$ratio >= 0.5, 1, 0)
epidermis$CpG_end <- epidermis$pos + 1
epidermis <- epidermis[, c(1,2,3,8,7,6,5)]
colnames(epidermis) <- c("expID", "chr",	"CpG_start",	"CpG_end",	"meth_frac",	"count_meth", "total_reads")

dermis$ratio <- ifelse(dermis$ratio >= 0.5, 1, 0)
dermis$CpG_end <- dermis$pos + 1
dermis <- dermis[, c(1,2,3,8,7,6,5)]
colnames(dermis) <- c("expID", "chr",	"CpG_start",	"CpG_end",	"meth_frac",	"count_meth", "total_reads")

#split
epidermis <- split(epidermis, epidermis$expID)
dermis <- split(dermis, dermis$expID)

#save files
setwd("/home/igor/Data/inhouse_epi_dermis_hg19/")
for (f in seq_along(epidermis)) {
  #epidermis[[c]]$expID = NULL
  fwrite(epidermis[[f]], sprintf("%s_cov.tsv.gz", names(epidermis)[f]), sep='\t')
  print(f)
}

for (f in seq_along(dermis)) {
  #epidermis[[c]]$expID = NULL
  fwrite(dermis[[f]], sprintf("%s_cov.tsv.gz", names(dermis)[f]), sep='\t')
  print(f)
}

#creating the true file
true_file <- data.table(cell_id = c(names(epidermis), names(dermis)), epigenotype_id = c(rep(1, 401), rep(2, 562)))
fwrite(true_file, "/home/igor/Data/Epiclomal_results/true_files/963_epi_d_true.txt.gz", sep='\t')


