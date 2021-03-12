### 21/01/2021 - uploading the published datasets to the mysql server

library(RMySQL)

mouse_public = dbConnect(MySQL(), user='i747v', password='igor',
                   dbname="mouse_public", host='zuse1db')

dbListTables(dbconn)

#test to add cpgCoverage for hg19 Zhu dataset.

tmp <- fread("/home/igor/Data/mysql_data/KC_2k/cov_300KC_all_chr_50_Zhu/scBS-2C-6-1_cov.tsv.gz")
tmp$CpG_end = NULL
tmp$expID = 1

tmp = tmp[,c(6,1,2,5,4,3)]
colnames(tmp) <- c("expID","chrom","pos",  "coverage", "mC", "ratio")



#mm9
dbconn = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_mm9", host='a130-pc-bioinfo-09')
scp_ih_hg19
#cpgExp
#cpgCoverage
#cpgBases
#knownGene

# to updaload datasets for mm9 db
#Smallwood 2014
#cpgExp

#create and upload cpgExp file
file_names = list.files("/home/igor/Data/Smallwood2014/")

cell_ids = gsub(".cov.txt.gz", "", file_names)

setwd("/home/igor/Data/mysql_data/scm_pub_mm9/")
mm9_cpgExp <- data.table(expID = seq_along(cell_ids),
                         name = cell_ids,
                         comment = c(rep("MII oocytes",12), rep("mouse ESCs cultured in 2i/LIF or serum/LIF", 39))) 

fwrite(mm9_cpgExp, "/home/igor/Data/mysql_data/scm_pub_mm9/mm9_cpgExp.tsv", sep = '\t')

dbSendStatement(dbconn, "LOAD DATA LOCAL INFILE '/home/igor/Data/mysql_data/scm_pub_mm9/mm9_cpgExp.tsv' INTO TABLE cpgExp;")

#29/01 new cpgExp

pub_mm9_cpgExp <- fread("/home/igor/Data/mysql_data/scm_pub_mm9/mm9_cpgExp.tsv")
cell_ids = strsplit(pub_mm9_cpgExp$V2, "_", fixed = TRUE) %>% unlist() %>% matrix(nrow = 126, byrow = TRUE) %>% as.data.frame()
new_cpgExp = matrix(nrow = 4, ncol = 3) %>% as.data.table()
new_cpgExp$V1 = c(1,13,25,32)
new_cpgExp$V2 = c("MII oocytes","mouse ESCs cultured in 2i/LIF","mouse ESCs cultured in 2i/LIF or serum/LIF","mouse ESCs cultured in serum/LIF")
new_cpgExp$V3 = c("1-12, GSE56879", "13-24, GSE56879", "25-31, GSE56879", "32-51, GSE56879")

fwrite(new_cpgExp, "/home/igor/Data/mysql_data/scm_pub_mm9/mm9_cpgExp.tsv", sep='\t')

dbSendStatement(mm9, "LOAD DATA LOCAL INFILE '/home/igor/Data/mysql_data/scm_pub_mm9/mm9_cpgExp.tsv' INTO TABLE cpgExp;")


#cpgCoverage

setwd("/home/igor/Data/Smallwood2014/")
file_names = list.files(".")[2:52]
cnt = 6001
for (f in seq_along(file_names)) {
  tmp <- fread(file_names[f])
  tmp$V3 <-  NULL
  tmp$expID <-  cnt
  tmp <-  tmp[,c(6,1,2,5,4,3)]
  colnames(tmp) <- c("expID","chrom","pos",  "coverage", "mC", "ratio")
  tmp$ratio <- ifelse(tmp$ratio >= 0.5, 1, 0)
  setkey(tmp, expID, chrom)
  dbWriteTable(scm_pub_mm10, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  cnt = cnt + 1
  print(f)
}

#upload cpgCoverage files
setwd("/home/igor/Data/Smallwood2014/mysql_proc/cpgCoverage")
cells = list.files(".")
for (f in c(2:51)) {
  tmp <- fread(cells[f])
  tmp$V2 <- paste0("chr", tmp$V2)
  colnames(tmp) <- c("expID","chrom","pos",  "coverage", "mC", "ratio") 
  dbWriteTable(dbconn, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  print(f)
}

#check
test = dbSendQuery(dbconn, "select * from cpgCoverage where expID in (30, 51)") %>% fetch(-1)

#17/02 edit real Smallwood mm9 data for Epiclomal
setwd("/home/igor/Data/Smallwood2014/mm9/")
lf = list.files(".")[2:52]

for (f in seq_along(lf)) {
  tmp <- fread(lf[f])
  tmp <- tmp[V3 == "+"]
  tmp$CpG_end <- tmp$V2 + 1
  tmp$total_reads <- tmp$V4 + tmp$V5
  tmp$meth_frac <- tmp$V4 / tmp$total_reads
  tmp <- tmp[complete.cases(tmp),]
  tmp$meth_frac <- ifelse(tmp$meth_frac >= 0.5, 1, 0)
  tmp <- tmp[,c(1,2,8,10,4,9)]
  colnames(tmp) = c('chr', 'CpG_start', 'CpG_end', 'meth_frac', 'count_meth', 'total_reads')
  tmp$chr <- paste0("chr", tmp$chr)
  fwrite(tmp, sprintf("/home/igor/Data/Smallwood2014/mm9/Epi_proc/%s_cov.tsv.gz", f), sep='\t')
  print(f)
}

#for mm9 dataset upload .CpG files
setwd("/home/igor/Data/Smallwood2014/mm9/")
lf = list.files(".")[2:52]
cnt <- 1
for (f in seq_along(lf)) {
  tmp <- fread(lf[f])
  tmp <- tmp[V3 == "+"]
  tmp$expID <- cnt
  tmp$coverage <- tmp$V4 + tmp$V5
  tmp$ratio <- tmp$V4 / tmp$coverage
  tmp <- tmp[complete.cases(tmp),]
  tmp$ratio <- ifelse(tmp$ratio >= 0.5, 1, 0)
  tmp <- tmp[,c(8,1,2,9,4,10)]
  colnames(tmp) <- c("expID","chrom","pos", "coverage", "mC", "ratio") 
  tmp$chrom <- paste0("chr", tmp$chrom)
  dbWriteTable(scm_pub_mm9, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  cnt <- cnt + 1
  print(f)
}

##Argelaguet and Mulqueen 200 cells

scm_pub_mm10 = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_mm10", host='a130-pc-bioinfo-09')

#cpgCoverage
setwd("/home/igor/Data/sci-MET_mouse_cortex_scNMT_new/")
file_names = list.files(".")
cnt = 1
for (f in seq_along(file_names)) {
  tmp <- fread(file_names[f])
  tmp$CpG_end <-  NULL
  tmp$expID <-  cnt
  tmp <-  tmp[,c(6,1,2,5,4,3)]
  colnames(tmp) <- c("expID","chrom","pos",  "coverage", "mC", "ratio")
  tmp$ratio <- ifelse(tmp$ratio >= 0.5, 1, 0)
  setkey(tmp, expID, chrom)
  fwrite(tmp, sprintf("/home/igor/Data/Smallwood2014/mysql_proc/cpgCoverage/%s", file_names[f]), sep='\t', col.names = FALSE)
  cnt = cnt + 1
  print(f)
}

#cpgExp
mm10_cpgExp_scNMT <- data.table(expID = c(1,268, 366, 463, 851), name = c("E4.5-5.5", "E6.5", "E6.75", "E7.5", "PS_VE") , comment = "Argelaguet et al 2019, scNMT-seq, met modality") 
scnmt_cells <- list.files("/home/igor/Data/scnmt_gastrulation/met/cpg_level") %>% gsub(".tsv.gz", "", .)
dbWriteTable(dbconn, "cpgExp", mm10_cpgExp_scNMT, append = TRUE, row.names = FALSE)

#from initial data
setwd("/home/igor/Data/scnmt_gastrulation/met/cpg_level/")

mm_cells <- list.files(".")
mm10_cell_ids = gsub(".tsv.gz", "", mm10_cell_ids)

cnt = 1
for (f in seq_along(mm_cells)) {
  tmp <- fread(mm_cells[f])
  tmp$total_reads <- tmp$met_reads	+ tmp$nonmet_reads
  tmp$expID <- cnt
  tmp <-tmp[,c(7,1,2,6,3,5)]
  colnames(tmp) <- c("expID","chrom","pos", "coverage", "mC", "ratio")
  setkey(tmp, expID, chrom)
  fwrite(tmp, sprintf("/home/igor/Data/mysql_data/scm_pub_mm10/scNMT/%s", mm_cells[f]), sep='\t', col.names = FALSE)
  cnt = cnt + 1
  print(f)
}

#upload scNMT cpgCoverage
setwd("/home/igor/Data/mysql_data/scm_pub_mm10/scNMT/")
scnmt_cells_full = list.files(".")
for (f in seq_along(scnmt_cells_full)) {
  tmp <- fread(scnmt_cells_full[f])
  tmp$V2 <- paste0("chr", tmp$V2)
  colnames(tmp) <- c("expID","chrom","pos", "coverage", "mC", "ratio")
  dbWriteTable(dbconn, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  print(f)
}

#Mulqueen data
setwd("/home/igor/Data/sci-MET_mouse_cortex_scNMT_new/")
mm_cells <- list.files(".")[1:100]
cell_ids <- gsub(".tsv.gz", "", mm_cells)

#from processed epiclomal data - not sure
cnt <- 1
for (f in seq_along(cell_ids)) {
  tmp = fread(mm_cells[f])
  tmp$CpG_end <-  NULL
  tmp$expID <-  cnt
  tmp <-  tmp[,c(6,1,2,5,4,3)]
  colnames(tmp) <- c("expID","chrom","pos",  "coverage", "mC", "ratio")
  tmp$ratio <- ifelse(tmp$ratio >= 0.5, 1, 0)
  setkey(tmp, expID, chrom)
  fwrite(tmp, sprintf("/home/igor/Data/mysql_data/scm_pub_mm10/sci-MET/%s", mm_cells[f]), sep='\t', col.names = FALSE)
  cnt = cnt + 1
  print(f)
}

#prepare mouse cortex data for mysql server
setwd("/home/igor/Data/sci-MET_mouse_cortex_unarch/")
mm_cells <- list.files(".")

cell_ids = strsplit(mm_cells, "_", fixed = TRUE) %>% unlist() %>% matrix(nrow = 100, byrow = TRUE) %>% as.data.frame()
#cell_ids = paste0(cell_ids$V3, "_", cell_ids$V4)

ls = list.dirs("/home/igor/Data/sci-MET_mouse_cortex_unarch/")[2:27]

cnt <- 2000
for (d in seq_along(ls)) {
  setwd(ls[d])
  file.list <- list.files(".")
  tmp <- do.call("rbind", lapply(file.list, FUN = function(file) {
    fread(file, header=TRUE, sep="\t")})) %>% setDT()
  tmp <- subset(tmp, tmp$mc_class %in% c("CGA", "CGC", "CGG", "CGN", "CGT"))
  #tmp$pos <- ifelse(tmp$strand == "-", tmp$pos - 1, tmp$pos)
  tmp[,c(3,4)] <- NULL
  tmp$chr <- paste0("chr", tmp$chr)
  tmp$expID <- cnt
  tmp <- tmp[,c(6,1,2,4,3,5)]
  colnames(tmp) <- c("expID","chrom","pos", "coverage","mC","ratio")
  dbWriteTable(dbconn, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  cnt <- cnt + 1
  print(d)
}

#cpgExp for sci-MET mm10
dbSendQuery(dbconn, "INSERT INTO cpgExp VALUES (2000,'Mouse cortex cells', 'Mulqueen et al, 2018, sci-MET, 2000-8169')")


#hg19 Zhu, Guo from scratch to MySQL
setwd("/home/igor/Data/Zhu_Guo2017_100")
cells = list.files(".")
cnt <- 1
for (f in c(2:length(cells))){
  tmp <- fread(cells[f])
  tmp <- tmp[tmp$Type == "CpG"]
  tmp$expID <- cnt
  tmp <- tmp[,c(11,1,2,5,6,8)]
  colnames(tmp) <- c("expID","chrom","pos", "coverage","mC","ratio")
  setkey(tmp, expID, chrom, pos)
  dbWriteTable(dbconn, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  cnt <- cnt + 1
  print(f)
}

#test mm9 ih 
mouse_public = dbConnect(MySQL(), user='i747v', password='igor',
                         dbname="mouse_public", host='zuse1db')

test = dbSendQuery(mouse_public, "select * from cpgCoverage where expID = 100310") %>% fetch(-1)

test2 = dbSendQuery(mouse_public, "select C.expID, B.chrom, B.pos, B.strand, C.ratio
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                    where expID = 100310") %>% fetch(-1)

#check hg19
dbconn = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_hg19", host='a130-pc-bioinfo-09')

hg19_bases = dbSendQuery(dbconn, "select * from cpgBases where pos = 24334") %>% fetch(-1)

#compare my mm10 db result
mm10 = dbConnect(MySQL(), user='igor', password='igor123',
                 dbname="scm_pub_mm10", host='a130-pc-bioinfo-09')

dbSendQuery(mm10, "select * from cpgBases LIMIT 10 OFFSET 1") %>% fetch(-1)


#new scimet processing  hg19
library(data.table)
library(dplyr)
#scimet_mm <- fread("/home/igor/Data/Mulqueen2018/mouse/GSE112554_MouseCortex_Barcoded_CG.bed.gz")
#scimet_mm <- scimet_mm[scimet_mm$V7 == "Mouse3"]
#test_un <- unique(tmp$BARCODE)

tmp <- fread("/home/igor/Data/Mulqueen2018/hg19/GSE112554_HumanCellLineSplit_Barcoded_CG.bed.gz")

setkey(tmp, BARCODE, CHR, START, Tn5_Tracked_Condition)
hg19_list <- split(tmp, tmp$BARCODE)
mdata <- data.table(barcode = tmp$BARCODE, lib_prep_cond = tmp$Library_Preparation_Condition, tn5_cond = tmp$Tn5_Tracked_Condition)
mdata <- distinct(mdata)

#change the order of Mulqueen cells
setkey(mdata, tn5_cond)
hg19_list <- hg19_list[mdata$barcode]

#26/02 to check if names were assigned correctly
mdata$cell_id <- c(1:641)

#upload cpgExp hg19
cpgExp_hg19 <- fread("/home/igor/Downloads/cpgExp_hg19_final.tsv.gz")
#dbSendStatement(scm_pub_hg19, "LOAD DATA LOCAL INFILE '/home/igor/Data/Mulqueen2018/hg19/cpgExp_hg19_final.tsv' INTO TABLE cpgExp;")
dbWriteTable(scm_pub_hg19, "cpgExp", cpgExp_hg19, append = TRUE, row.names = FALSE)

#for mysql
cnt <- 1001
for (f in seq_along(hg19_list)) {
  tmp <- hg19_list[[f]]
  tmp$PERC_MET <- tmp$PERC_MET / 100
  tmp$expID <- cnt
  tmp$coverage <- tmp$C_MET_COUNT + tmp$C_UNMET_COUNT
  tmp <- tmp[, c(10,1,2,11,5,4)]  
  colnames(tmp) <- c("expID","chrom","pos", "coverage","mC","ratio")
  dbWriteTable(scm_pub_hg19, "cpgCoverage", tmp, append = TRUE, row.names = FALSE)
  cnt <- cnt + 1
  print(f)
}

#01/02 db checks

scm_pub_mm9 = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_mm9", host='a130-pc-bioinfo-09')
scm_pub_mm10 = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_mm10", host='a130-pc-bioinfo-09')
scm_pub_hg19 = dbConnect(MySQL(), user='igor', password='igor123',
                   dbname="scm_pub_hg19", host='a130-pc-bioinfo-09')

dbReadTable(scm_pub_mm9, "cpgExp")
dbReadTable(scm_pub_mm10, "cpgExp")
dbReadTable(scm_pub_hg19, "cpgExp")

#unique expID hg19
test = dbSendQuery(scm_pub_hg19, "SELECT DISTINCT expID from cpgCoverage") %>% fetch(-1)

#DELETE MOUSE CELLS IN ZHOU HG19 DATASET
dbSendQuery(scm_pub_hg19, "DELETE FROM cpgCoverage WHERE expID BETWEEN 400 AND 435") %>% fetch(-1)

#02/02 to add publication names
cpgExp <- dbReadTable(scm_pub_mm9, "cpgExp")
dbSendQuery(scm_pub_mm9, "DELETE FROM cpgExp") %>% fetch(-1)
dbWriteTable(scm_pub_mm9, "cpgExp", cpgExp, append = TRUE, row.names = FALSE)

cpgExp <- dbReadTable(scm_pub_mm10, "cpgExp")
dbSendQuery(scm_pub_mm10, "DELETE FROM cpgExp") %>% fetch(-1)
dbWriteTable(scm_pub_mm10, "cpgExp", cpgExp, append = TRUE, row.names = FALSE)

cpgExp <- dbReadTable(scm_pub_hg19, "cpgExp")
dbSendQuery(scm_pub_hg19, "DELETE FROM cpgExp") %>% fetch(-1)
dbWriteTable(scm_pub_hg19, "cpgExp", cpgExp, append = TRUE, row.names = FALSE)

#create own table
dbSendQuery(scm_pub_mm9, "CREATE TABLE IF NOT EXISTS dataInfo (
    id INT,
    article TEXT,
    link TEXT,
    protocol TEXT,
    n_cells INT
)  ENGINE=INNODB;") %>% fetch(-1)

dbSendQuery(scm_pub_mm10, "CREATE TABLE IF NOT EXISTS dataInfo (
    id INT,
    article TEXT,
    link TEXT,
    protocol TEXT,
    n_cells INT
)  ENGINE=INNODB;") %>% fetch(-1)

dbSendQuery(scm_pub_hg19, "CREATE TABLE IF NOT EXISTS dataInfo (
    id INT,
    article TEXT,
    link TEXT,
    protocol TEXT,
    n_cells INT
)  ENGINE=INNODB;") %>% fetch(-1)

#fill it in
dbSendQuery(scm_pub_mm9, "INSERT INTO dataInfo VALUES (1, '[Smallwood et al, 2014]', 'GSE56879', 'scBS-seq', 51)")

dbSendQuery(scm_pub_mm10, "INSERT INTO dataInfo VALUES (1, '[Argelaguet et al, 2019]', 'GSE121708', 'scNMT-seq', '1140')")
dbSendQuery(scm_pub_mm10, "INSERT INTO dataInfo VALUES (2,  '[Luo et al, 2017]', 'GSE97179', 'MethylC-seq', 3385)")

dbSendQuery(scm_pub_hg19, "INSERT INTO dataInfo VALUES (1,  '[Zhu, Guo, Ren et al, 2018]', 'GSE81233','scBS-seq', 564)")
dbSendQuery(scm_pub_hg19, "INSERT INTO dataInfo VALUES (2,  '[Mulqueen et al, 2018]', 'GSM3072971','sci-MET', 641)")

#add column to the cpgExp table
dbSendQuery(scm_pub_mm9, "ALTER TABLE cpgExp
ADD COLUMN idDataset TEXT AFTER comment;") %>% fetch(-1)

#change cpgExp
cpgExp_mm9 <- dbSendQuery(scm_pub_mm9, "select * from cpgExp") %>% fetch(-1)
cpgExp_mm9$comment = gsub(", GSE56879 [Smallwood et al, 2014]", "", cpgExp_mm9$comment, fixed = TRUE)

dbSendQuery(scm_pub_hg19, "DELETE FROM cpgExp") %>% fetch(-1)
dbWriteTable(scm_pub_hg19, "cpgExp", test, append = TRUE, row.names = FALSE)

#22/02 Upload mm9 Smallwood data to the public dataset.


