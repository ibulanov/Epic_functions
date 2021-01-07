library(data.table)
library(REpiclomal)

ten_cells <- list.files("/Users/igorbulanov/Documents/Data/KC_1700")[1:10]
all_cells <- list.files("/Users/igorbulanov/Documents/Data/KC_chr1/300_cpgCoverage")


for (i in c(1:10)){
  tmp_file <- fread(ten_cells[i])
  tmp_file <- tmp_file[tmp_file$V1=="chr1"]
  tmp_file$V1 <- gsub("chr", "", tmp_file$V1)
  fwrite(tmp_file, sprintf("/Users/igorbulanov/Documents/Data/KC_chr1/%s", ten_cells[i]), sep='\t')
}

ten_cells_chr1 <- data.frame()
for (i in c(1:1700)){
  tmp_file <- fread(ten_cells[i])
  tmp_file <- tmp_file[tmp_file$V1=="chr1"]
  tmp_file$V1 <- gsub("chr", "", tmp_file$V1)
  tmp_file$cell_id <- ten_cells[i]
  ten_cells_chr1 <- rbind(ten_cells_chr1, tmp_file)
  print(i)
}

#all cells
all_cells_chr1 <- data.frame()
for (i in seq_along(all_cells) {
  tmp_file <- fread(all_cells[i])
  tmp_file <- tmp_file[tmp_file$V1=="chr1"]
  tmp_file$V1 <- gsub("chr", "", tmp_file$V1)
  tmp_file$cell_id <- all_cells[i]
  all_cells_chr1 <- rbind(all_cells_chr1, tmp_file)
  print(i)
}

#16 most complete cells
all_cells_chr1 <- data.frame()
for (i in c(1:16)){
  tmp_file <- fread(all_cells[mmm_candidates][i])
  tmp_file <- tmp_file[tmp_file$V1=="chr1"]
  tmp_file$V1 <- gsub("chr", "", tmp_file$V1)
  tmp_file$cell_id <- all_cells[mmm_candidates][i]
  all_cells_chr1 <- rbind(all_cells_chr1, tmp_file)
  print(i)
}

all_cells_chr1$cell_id <- gsub(".tsv.gz", "", all_cells_chr1$cell_id)
all_cells_chr1$pos_d <- all_cells_chr1$V2
all_cells_chr1 <- all_cells_chr1[, c(4,1,2,5,3)]
colnames(all_cells_chr1) <- c("cell_id", "chr", "pos", "pos_d", "ratio")

ten_cells_chr1_test <- all_cells_chr1[c(1:1000),] 

fwrite(all_cells_chr1, "/Users/igorbulanov/Documents/Data/KC_chr1/cpgs_16cells_chr1.tsv", sep='\t')

annos <- fread("/Users/igorbulanov/Documents/Data/Annotations/mm9_annos_200.tsv")

annos_test <- annos[c(1:100),]

fwrite(annos_test, "/Users/igorbulanov/Documents/Data/test/annos_test.tsv", sep='\t')

anno_2k <- fread("/Users/igorbulanov/Documents/Data/Annotations/anno_2k_chr1.tsv")

epidata29 <- load_data(outdir = "/Users/igorbulanov/Documents/Data/",
                 input_CpG_data_file = "/Users/igorbulanov/Documents/Data/KC_chr1/cpgs10cells_chr1_29k.tsv",
                 input_regions_file = "/Users/igorbulanov/Documents/Data/Annotations/annos_mm9_37k_chr.tsv",
                 use_cache = FALSE)

# process sci-MET brain cells
brain_cells <- list.files("/Users/igorbulanov/Documents/Data/ESC77_BRAIN20/")
for (i in c(1:20)) {
  tmp <- fread(brain_cells[i])
  tmp[,c(3:6)] <- NULL
  tmp$chr <- paste0("chr", tmp$chr)
  fwrite(tmp, sprintf("/Users/igorbulanov/Documents/Data/ESC77_BRAIN20/proc/%s", brain_cells[i]), sep='\t')
  print(i)
}

#filter columns by number of NA
real_obs <- data.table(i = c(1:1652), count = NA)
for (i in c(1:1652)) {
  cur_value <- sum(!is.na(mean_meth[,i]))
  real_obs$count[i] <- cur_value
}

setkey(real_obs, count)

mmm_candidates <- real_obs$i[1637:1652]


#save cpg coords by cell separately
list_files = rownames(epiclomal_input$mean_meth_matrix)

for (i in c(1:16)) {
  cell_16_split[[i]]$cell_id <- NULL
  fwrite(cell_16_split[[i]], sprintf("/Users/igorbulanov/Documents/Data/KC_chr1/16_cells/%s.tsv.gz", list_files[i]), sep='\t')
}

cpgCoverage <- split(cpgCoverage, cpgCoverage$expID)
for (i in c(1:300)) {
  cpgCoverage[[i]]$expID <- NULL
  fwrite(cpgCoverage[[i]], sprintf("/Users/igorbulanov/Documents/Data/KC_chr1/300_cpgCoverage/coverage_%s.tsv.gz", names(cpgCoverage)[i]), sep='\t')
}

#all cells remove other chrs and save separately
all_cells_chr1 <- data.frame()
for (i in seq_along(all_cells)) {
  tmp_file <- fread(all_cells[i])
  tmp_file <- tmp_file[tmp_file$V1=="chr1"]
  fwrite(tmp_file, sprintf("/Users/igorbulanov/Documents/Data/KC_chr1/30_cells/%s", all_cells[i]), sep='\t')
}

#compute files size
file_sizes <- fread("/Users/igorbulanov/Documents/Data/KC_chr1/true_clone_memb_all_kc.txt.gz")
file_sizes$size <- file.info(sprintf("%s.tsv.gz", file_sizes$cell_id))$size

for (i in seq_along(command)) {
  tmp <- fread(command[i])
  fwrite(tmp, sprintf("/Users/igorbulanov/Documents/Data/KC_chr1/300_cells/%s", command[i]), sep='\t')
  print(i)
}  

#add 1 to each CpG_end coords
for (f in seq_along(all_cells)) {
  tmp <- fread(all_cells[f])
  tmp$CpG_end <- tmp$CpG_end + 1
  fwrite(tmp, all_cells[f], sep='\t')
  print(f)
}

for (i in seq(1:nrow(input))) {
  input[i,is.na(input[i,])] <- mean(input[i,],na.rm=TRUE)
}

#create table with cell id and cons cluster
cluster_posteriors <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/24_11_epiclomal_results_300KC_chr1/epi_result_done/0_0.95_10000/epi_region/663/cluster_posteriors.tsv.gz")
library(tibble)
cluster_posteriors <- cluster_posteriors %>% column_to_rownames(var = "cell_id")
cells_cl <- data.frame(cell_id = rownames(cluster_posteriors),  cl_n = NA)

for (r in c(1:300)) {
  cells_cl$cl_n[r] <- (match(max(cluster_posteriors[r,]), cluster_posteriors[r,])-1)
}

#rename coverage cells
list_cov_files = list.files("/Users/igorbulanov/Documents/Data/KC_chr1/300_cpgCoverage")


for (f in seq_along(kc300_list)) {
  tmp_f <- fread(kc300_list[f])
  fwrite(tmp_f, sprintf("/Users/igorbulanov/Documents/Data/KC_coverage/coverage_300_all_chr/%s_cov.tsv.gz", cell_ids[f]))
  print(f)
}

#### 04/01/21 Fill in the sparse CpG matrix by imputed CpGs from Epiclomal####

#define cell to cluster assignments table
library(dplyr)
library(tibble)

cluster_map <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/epi_region/944/cluster_MAP.tsv.gz")

input_epiclomal <- column_to_rownames(input_epiclomal, "cell_id") %>% t()
region_map <- column_to_rownames(region_map, "cluster_id") %>% t()

result_comb <- input_epiclomal
for (c in c(1:300)) {
  tmp_table <- data.table(obs = input_epiclomal[,c], imp = region_map[,cluster_map$cluster[c]]) %>% 
    mutate(res = coalesce(obs, imp))
  result_comb[,c] <- tmp_table$res
  print(c)
}

#fwrite(result_comb, "/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/epiclomal_comb.tsv", sep='\t')
fwrite(result_comb, "/home/igor/Data/Epiclomal_results/24_12_KC300_chr_all/5_0.99_20000/visual/epiclomal_comb_5_99_20k_prom1kb.tsv", sep='\t')

#replace 0's by 0.1 in the entire dataset and replace NA by 0 then
test5[test5==0] <- 0.1
test5[test5==0] <- 0.1
test5[is.na(test5)] <- 0

for (i in c(3:8)) {
  test5[,i] <- ifelse(is.na(test5[,i]), 0, test5[,i])
}

#save file for ucsc
for (i in c(1:6)){
  fwrite(list_of_cl[[i]], sprintf("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/%s.tsv",files_id[i]),sep ='\t', col.names = FALSE)
}

## comb of clusters AND cells
#comb_cls_cells <- cbind(region_map, epiclomal_comb)


#find the most variable cells within cluster
cells_id_cl1 <- cluster_map$cell_id[cluster_map$cluster==1] %>% as.character()
cells_id_cl2 <- cluster_map$cell_id[cluster_map$cluster==2] %>% as.character()
cells_id_cl3 <- cluster_map$cell_id[cluster_map$cluster==3] %>% as.character()
cells_id_cl4 <- cluster_map$cell_id[cluster_map$cluster==4] %>% as.character()
cells_id_cl5 <- cluster_map$cell_id[cluster_map$cluster==5] %>% as.character()
cells_id_cl6 <- cluster_map$cell_id[cluster_map$cluster==6] %>% as.character()
cells_id_cl7 <- cluster_map$cell_id[cluster_map$cluster==7] %>% as.character()

cpgs_cells_cl1 <- as.matrix(epiclomal_comb[,cells_id_cl1])
cpgs_cells_cl2 <- as.matrix(epiclomal_comb[,cells_id_cl2])
cpgs_cells_cl3 <- as.matrix(epiclomal_comb[,cells_id_cl3])
cpgs_cells_cl4 <- as.matrix(epiclomal_comb[,cells_id_cl4])
cpgs_cells_cl5 <- as.matrix(epiclomal_comb[,cells_id_cl5])
cpgs_cells_cl6 <- as.matrix(epiclomal_comb[,cells_id_cl6])
cpgs_cells_cl7 <- as.matrix(epiclomal_comb[,cells_id_cl7])


#test-  c(1:nrow-1)
comp_tbl_test <- matrix(nrow = 37, ncol = 2)

for (clm in c(1:37)) {
  tmp_comp <- which(cpgs_cells_cl1[,clm] != cpgs_cells_cl1[,clm+1])
  comp_tbl_test[clm,1] <- paste0(colnames(cpgs_cells_cl1)[clm], ":",colnames(cpgs_cells_cl1)[clm+1], ",")
  comp_tbl_test[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl1 <- matrix(nrow = length(cells_id_cl1), ncol = 2)
for (clm in c(1:length(cells_id_cl1)-1)) {
  tmp_comp <- which(cpgs_cells_cl1[,clm] != cpgs_cells_cl1[,clm+1])
  comp_tbl_cl1[clm,1] <- paste0(colnames(cpgs_cells_cl1)[clm], ",",colnames(cpgs_cells_cl1)[clm+1])
  comp_tbl_cl1[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl2 <- matrix(nrow = length(cells_id_cl2), ncol = 2)
for (clm in c(1:length(cells_id_cl2)-1)) {
  tmp_comp <- which(cpgs_cells_cl2[,clm] != cpgs_cells_cl2[,clm+1])
  comp_tbl_cl2[clm,1] <- paste0(colnames(cpgs_cells_cl2)[clm], ",",colnames(cpgs_cells_cl2)[clm+1])
  comp_tbl_cl2[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl3 <- matrix(nrow = length(cells_id_cl3), ncol = 2)
for (clm in c(1:length(cells_id_cl3)-1)) {
  tmp_comp <- which(cpgs_cells_cl3[,clm] != cpgs_cells_cl3[,clm+1])
  comp_tbl_cl3[clm,1] <- paste0(colnames(cpgs_cells_cl3)[clm], ",",colnames(cpgs_cells_cl3)[clm+1])
  comp_tbl_cl3[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl4 <- matrix(nrow = length(cells_id_cl4), ncol = 2)
for (clm in c(1:length(cells_id_cl4)-1)) {
  tmp_comp <- which(cpgs_cells_cl4[,clm] != cpgs_cells_cl4[,clm+1])
  comp_tbl_cl4[clm,1] <- paste0(colnames(cpgs_cells_cl4)[clm], ",",colnames(cpgs_cells_cl4)[clm+1])
  comp_tbl_cl4[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl5 <- matrix(nrow = length(cells_id_cl5), ncol = 2)
for (clm in c(1:length(cells_id_cl5)-1)) {
  tmp_comp <- which(cpgs_cells_cl5[,clm] != cpgs_cells_cl5[,clm+1])
  comp_tbl_cl5[clm,1] <- paste0(colnames(cpgs_cells_cl5)[clm], ",",colnames(cpgs_cells_cl5)[clm+1])
  comp_tbl_cl5[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl6 <- matrix(nrow = length(cells_id_cl6), ncol = 2)
for (clm in c(1:length(cells_id_cl6)-1)) {
  tmp_comp <- which(cpgs_cells_cl6[,clm] != cpgs_cells_cl6[,clm+1])
  comp_tbl_cl6[clm,1] <- paste0(colnames(cpgs_cells_cl6)[clm], ",",colnames(cpgs_cells_cl6)[clm+1])
  comp_tbl_cl6[clm,2] <- toString(tmp_comp)
}

comp_tbl_cl7 <- matrix(nrow = length(cells_id_cl7), ncol = 2)
for (clm in c(1:length(cells_id_cl7)-1)) {
  tmp_comp <- which(cpgs_cells_cl7[,clm] != cpgs_cells_cl7[,clm+1])
  comp_tbl_cl7[clm,1] <- paste0(colnames(cpgs_cells_cl7)[clm], ",",colnames(cpgs_cells_cl7)[clm+1])
  comp_tbl_cl7[clm,2] <- toString(tmp_comp)
}

#old
#cpg = rownames(epiclomal_comb)[c(22053, 31725,33680,34256,44604,24327,660,39782,28203,5093,6209,25700,41098,20830,4222,10773,36142,15458,40218,45541,27341,30071,1886,4240)]
#new
cpg = rownames(epiclomal_comb)[c(20125, 22053, 31726,16284,44779,660,13827,24328,25389,39782,5093,6209,9800,41098,28203,4257,20831,21461,5485,9926,10773,42557,24186,27342,30069,45863,6501,36145,15456,21941,25682,26081,34143,23479,39396,4240,13072)]

cpg <- rownames(epiclomal_comb)[c(20126,22053,31726,16284,44779,660,13827,24328,25389,4257,45863,36145,15456,39396,4240)]

library(stringr)
cpgs <- data.table(cpg = cpg,
                   chrom = str_split(cpg, ":", simplify = TRUE)[,1],
                   pos = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric(),
                   pos_d = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric()
)

 
####create table of regions####
filtered_regions = data.table(region = filter_regions$V1, 
                              chrom = str_split(filter_regions$V1, ":", simplify = TRUE)[,1],
                              start = str_split(filter_regions$V1, "[[:punct:]]", simplify = TRUE)[,2] %>% as.numeric(),
                              end = str_split(filter_regions$V1, "[[:punct:]]", simplify = TRUE)[,3] %>% as.numeric()
)
filtered_regions$chrom = gsub("chr", "", filtered_regions$chrom)

setkey(cpgs, chrom, pos, pos_d)
setkey(filtered_regions, chrom, start, end)

overlap_cpgs_regions  <- foverlaps(filtered_regions, cpgs, nomatch = NULL)

#create table of best cells and desired regions

epiclomal_best_cells <- result_comb[,c(cl1_cells_id,cl2_cells_id,cl3_cells_id,cl4_cells_id,cl5_cells_id,cl6_cells_id,cl7_cells_id)]

# separate table on columns
#epiclomal_best_cells <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/epiclomal_best_cells.tsv", sep='\t', col.names = TRUE)
rn <- rownames(epiclomal_best_cells)

epiclomal_best_cells <- as.data.frame(epiclomal_best_cells) %>% setDT()
epiclomal_best_cells$rn <- epi_rn

#add additional columns
epiclomal_best_cells$chrom <- str_split(epiclomal_best_cells$rn, ":", simplify = TRUE)[,1]
#epiclomal_best_cells$chrom <- paste0("chr", epiclomal_best_cells$chrom)
epiclomal_best_cells$start <- str_split(epiclomal_best_cells$rn, ":", simplify = TRUE)[,2] %>% as.numeric()
epiclomal_best_cells$end <- epiclomal_best_cells$start + 1 
epiclomal_best_cells <- cbind(epiclomal_best_cells, region_map)

#mutate
epiclomal_best_cells <- epiclomal_best_cells[, c(36:46, 1:35)]

#merge with desired list of regions
setkey(epiclomal_best_cells, chrom, start, end)
setkey(desired_regions, chrom, start,end)
epiclomal_best_cells_regions <- foverlaps(epiclomal_best_cells, desired_regions, nomatch = NULL)
epiclomal_best_cells_regions <- epiclomal_best_cells_regions[,c(1,6,7:49)]
colnames(epiclomal_best_cells_regions)[1:10] <- c("chrom", "start", "end", "cl_1", "cl_2", "cl_3", "cl_4", "cl_5", "cl_6", "cl_7")

fwrite(epiclomal_best_cells_regions, "/home/igor/Data/Epiclomal_results/24_12_KC300_chr_all/5_0.99_20000/visual/epiclomal_best_cells_regions.tsv", sep='\t', col.names = TRUE)

epiclomal_best_cells_regions <- as.matrix(epiclomal_best_cells_regions)

setwd("/home/igor/Data/Epiclomal_results/24_12_KC300_chr_all/5_0.99_20000/visual/cls_cells_ucsc")
count = 1
for (i in c(4:45)) {
  tmp <- epiclomal_best_cells_regions[,c(1,2,3,i)] 
  ident <- sprintf("track type=bedGraph name=%s color=0,0,255 maxHeightPixels=20,12,10 viewLimits=0.0:1.0 windowingFunction=mean",
                   colnames(tmp)[4])
  id <- colnames(tmp)[4]
  ident_dt <- data.table(id = ident, V2 = NA, V3 = NA, V4 = NA)
  tmp <- rbind(ident_dt, tmp, use.names=FALSE)
  fwrite(tmp, sprintf("%s_%s.bed", count,id),sep='\t', col.names = FALSE)
  count = count + 1
  print(i)
}



