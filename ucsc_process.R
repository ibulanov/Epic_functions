library(data.table)

#split list of cl by regions
for (i in seq_along(list_of_cl)) {
  list_of_cl[[i]] <- split(list_of_cl[[i]], list_of_cl[[i]]$region_id)
}

#delete region_id column in each file
for (i in seq_along(list_of_cl)) {
  for (j in seq_along(names(list_of_cl[[1]]))) {
    list_of_cl[[i]][[j]]$region_id <- NULL
  }
}

#save all files
setwd("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc")
for (i in seq_along(list_of_cl)) {
  for (j in seq_along(names(list_of_cl[[1]]))) {
    fwrite(list_of_cl[[i]][[j]], sprintf("%s_%s.tsv",names(list_of_cl)[i], names(list_of_cl[[i]])[j]), sep = '\t', col.names = FALSE)
  }
}

test <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/cl_0_chr1:9647976-9648975.tsv")


# save files per cluster
all_clusters <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/all_clusters_3reg.tsv")

for (f in c(5:10)) {
  tmp <- all_clusters[,c(2,3,4,f)]
  print(head(tmp))
}

#take cells based on their sizes
setwd("/Users/igorbulanov/Documents/Data/KC_coverage/coverage_300_all_chr")
list_cells <- list.files("/Users/igorbulanov/Documents/Data/KC_coverage/coverage_300_all_chr")
cells_sizes <- data.table(cells = list_cells, size = file.info(list_cells)$size)
cells_sizes$cells <- gsub("_cov.tsv.gz", "", cells_sizes$cells)


#new top5 biggest promoters!
library(stringr)
cl_cells_top5r_ucsc <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/top5_diff_chr/epi_cl_cells_top5r.tsv")
cl_cells_top5r_ucsc$chrom <- str_split(cl_cells_top5r_ucsc$rn, ":", simplify = TRUE)[,1]
cl_cells_top5r_ucsc$chrom <- paste0("chr", cl_cells_top5r_ucsc$chrom)
cl_cells_top5r_ucsc$start <- str_split(cl_cells_top5r_ucsc$rn, ":", simplify = TRUE)[,2] %>% as.numeric()
cl_cells_top5r_ucsc$end <- cl_cells_top5r_ucsc$start + 1 

final_top5_ucsc <- as.matrix(final_top5_ucsc)

setwd("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/top5_diff_chr")
count = 1
for (i in c(4:39)) {
  tmp <- final_top5_ucsc[,c(1,2,3,i)] 
  ident <- sprintf("track type=bedGraph name=%s color=0,0,255 maxHeightPixels=20,12,10 viewLimits=0.0:1.0 windowingFunction=mean",
                   colnames(tmp)[4])
  id <- colnames(tmp)[4]
  ident_dt <- data.table(id = ident, V2 = NA, V3 = NA, V4 = NA)
  tmp <- rbind(ident_dt, tmp, use.names=FALSE)
  fwrite(tmp, sprintf("%s_%s.bed", count,id),sep='\t', col.names = FALSE)
  count = count + 1
  print(i)
}

#check how much real cpg obs
input_cells <- input_epiclomal[,c("10134", "10792", "10176", "10293", "10326",
                                  "10094", "10080", "10147", "10044", "10215",
                                  "10210", "11509", "10034", "10143", "11019",
                                  "10964", "10840", "10860", "10249", "10342",
                                  "10425", "10355", "10349", "10339", "10240",
                                  "10300", "11415", "10327", "10413", "10146")]


real_cpgs = data.table(cells = colnames(input_cells), real_obs = NA)
real_cpgs$real_obs = sapply(input_cells, function(x) sum(!is.na(x)))

library(tidyverse)
input_cells %>%
  select(everything()) %>%  # replace to your needs
  summarise_all(funs(sum(!is.na(.))))


####06_01_21 another tactics to find most variable cells and regions for 10k cpgs and gene1kb annos####
library(tibble)
library(dplyr)
library(data.table)

best_dir <- "/home/igor/Data/Epiclomal_results/17_01_200_mm10_sci-MET_scNMT_clustering/5_0.95_10000/epi_region/8/"
input_dir <- "/home/igor/Data/Epiclomal_results/16_01_200_sci-MET_scNMT_mm10_preproc/epiclomal_input/5_0.95_10000/"
output_dir <- "/home/igor/Data/Epiclomal_results/17_01_200_mm10_sci-MET_scNMT_clustering/5_0.95_10000/visual/"

cluster_map <- fread(paste0(best_dir, "cluster_MAP.tsv.gz"), header = TRUE) %>%
  .[,c(1,2)] 
colnames(cluster_map) <- c("cell_id", "cluster")
cluster_map$cluster <- cluster_map$cluster + 1

setwd(input_dir)
input_epiclomal <- fread(paste0(input_dir, "input_Epiclomal_top200_mm10.tsv.gz")) %>%
  column_to_rownames(., "cell_id") %>%
  t()


region_map <- fread(paste0(best_dir, "genotype_MAP.tsv.gz")) %>% 
  column_to_rownames(., "cluster_id") %>%
  t()

#change colnames in region_map
colnames(region_map) <- paste0("cl_", as.numeric(colnames(region_map)) + 1)

result_comb <- input_epiclomal
for (c in c(1:nrow(cluster_map))) {
  tmp_table <- data.table(obs = input_epiclomal[,c], imp = region_map[,cluster_map$cluster[c]]) %>% 
    mutate(res = coalesce(obs, imp))
  result_comb[,c] <- tmp_table$res
  print(c)
}

fwrite(result_comb, paste0(output_dir, "epiclomal_comb_5_95_10k_200_mm10.tsv"),
       sep='\t',
       col.names = TRUE)

#M3C umap
library(M3C)
M3C::umap(result_comb,labels=as.factor(cluster_map$cluster),controlscale=TRUE,scale=3)

#find the most variable cells within cluster
cell_ids <- list()
for(c in c(1:ncol(region_map))) {
  cell_ids[[c]] <- cluster_map$cell_id[cluster_map$cluster==c] %>% as.character()
}
names(cell_ids) <- colnames(region_map)

cpgs_cells <- list()
for(c in c(1:ncol(region_map))) {
  cpgs_cells[[c]] <- as.matrix(result_comb[,cell_ids[[c]]])
}
names(cpgs_cells) <- colnames(region_map)


comp_tbl <- list()
comp_tbl[[1]] <- matrix(nrow = (ncol(cpgs_cells[[1]]) ** 2), ncol = 4)
comp_tbl[[2]] <- matrix(nrow = (ncol(cpgs_cells[[2]]) ** 2), ncol = 4)
comp_tbl[[3]] <- matrix(nrow = (ncol(cpgs_cells[[3]]) ** 2), ncol = 4)
comp_tbl[[4]] <- matrix(nrow = (ncol(cpgs_cells[[4]]) ** 2), ncol = 4)
comp_tbl[[5]] <- matrix(nrow = (ncol(cpgs_cells[[5]]) ** 2), ncol = 4)

for (r in c(1:5)) {
  r_count <-  1
  for (c in c(1:ncol(cpgs_cells[[r]]))) {
    for (l in c(1:ncol(cpgs_cells[[r]]))) {
      tmp_v <- which(cpgs_cells[[r]][,c] != cpgs_cells[[r]][,l])
      comp_tbl[[r]][r_count,1] <- colnames(cpgs_cells[[r]])[c]
      comp_tbl[[r]][r_count,2] <- colnames(cpgs_cells[[r]])[l]
      comp_tbl[[r]][r_count,3] <- toString(tmp_v)
      comp_tbl[[r]][r_count,4] <- length(tmp_v)
      r_count <- r_count + 1
    }
    print(c)
  }
  comp_tbl[[r]] <- as.data.table(comp_tbl[[r]])
  setkey(comp_tbl[[r]], V4)
}

#for 5 cl separately
#to see count of diff cpgs between two cells
table(comp_tbl[[1]]$V4)

most_diff <- subset(comp_tbl[[1]], (V4 %in% c(11:17)))

sort(table(most_diff$V1))

#cl_N_ids = c()

most_diff <- subset(most_diff, (V1 %in% cl_5_ids)) %>% 
  subset(., (V2 %in% cl_5_ids))

unique(most_diff$V3) %>% paste(., collapse = ',') %>% strsplit(., ",") %>% .[[1]] %>% as.numeric(.) %>% table(.)

#create table of most variable cells and regions

cpg <- rownames(result_comb)[c(4872,6129,7271,11687,22054,27441,37767,44783,4242,4715,24906,44581,43996,3245,2722,29128,30446,42424,45982,45540,2695,25188,25723,27517,28493,33678,34219,44711,7722,24326,25149,32695,43673,44670)]

library(stringr)
cpgs <- data.table(cpg = cpg,
                   chrom = str_split(cpg, ":", simplify = TRUE)[,1],
                   pos = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric(),
                   pos_d = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric()
)

filter_regions <- fread("/home/igor/Data/Epiclomal_results/24_12_KC300_chr_all_20k/filter_regions/5_0.99_20000/filtered_regions_KC300_all_cpg.tsv", col.names = "V1")

####create table of regions####
filtered_regions <- data.table(region = filter_regions$V1, 
                              chrom = str_split(filter_regions$V1, ":", simplify = TRUE)[,1],
                              start = str_split(filter_regions$V1, "[[:punct:]]", simplify = TRUE)[,2] %>% as.numeric(),
                              end = str_split(filter_regions$V1, "[[:punct:]]", simplify = TRUE)[,3] %>% as.numeric()
)
filtered_regions$chrom = gsub("chr", "", filtered_regions$chrom)

setkey(cpgs, chrom, pos, pos_d)
setkey(filtered_regions, chrom, start, end)

overlap_cpgs_regions  <- foverlaps(filtered_regions, cpgs, nomatch = NULL)

#create table of best cells and desired regions

epiclomal_best_cells <- result_comb[,c(cl_1_ids, cl_2_ids,cl_3_ids,cl_4_ids,cl_5_ids)]

# separate table on columns
#epiclomal_best_cells <- fread("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/epiclomal_best_cells.tsv", sep='\t', col.names = TRUE)
rn <- rownames(epiclomal_best_cells)

epiclomal_best_cells <- as.data.frame(epiclomal_best_cells) %>% setDT()
epiclomal_best_cells$rn <- rn

#add additional columns
epiclomal_best_cells$chrom <- str_split(epiclomal_best_cells$rn, ":", simplify = TRUE)[,1]
#epiclomal_best_cells$chrom <- paste0("chr", epiclomal_best_cells$chrom)
epiclomal_best_cells$start <- str_split(epiclomal_best_cells$rn, ":", simplify = TRUE)[,2] %>% as.numeric()
epiclomal_best_cells$end <- epiclomal_best_cells$start + 1 
epiclomal_best_cells <- cbind(epiclomal_best_cells, region_map)

#mutate
epiclomal_best_cells <- epiclomal_best_cells[, c(27:35, 1:26)]

desired_regions <- overlap_cpgs_regions[,c(5,1,6,7)] %>% distinct(., .keep_all = TRUE)

#merge with desired list of regions
setkey(epiclomal_best_cells, chrom, start, end)
setkey(desired_regions, chrom, start,end)
epiclomal_best_cells_regions <- foverlaps(epiclomal_best_cells, desired_regions, nomatch = NULL)
epiclomal_best_cells_regions <- epiclomal_best_cells_regions[,c(1,6,7:38)]
colnames(epiclomal_best_cells_regions)[1:8] <- c("chrom", "start", "end", "cl_1", "cl_2", "cl_3", "cl_4", "cl_5")

#colnames(epiclomal_best_cells_regions)[9:13] <- paste0("cl1_", colnames(epiclomal_best_cells_regions)[9:13])

fwrite(epiclomal_best_cells_regions, paste0(output_dir, "epiclomal_best_cells_regions.tsv"),
       sep='\t',
       col.names = TRUE)

epiclomal_best_cells_regions <- as.matrix(epiclomal_best_cells_regions)

setwd("/home/igor/Data/Epiclomal_results/24_12_KC300_chr_all_20k/visual/cls_cells_ucsc")
count = 1
for (i in c(4:34)) {
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


#preparing the bar plot
count_cpgs <- data.table(cell_ids = colnames(input_epiclomal), unmeth = NA, meth = NA) %>% as.matrix()
count_cpgs2 <- data.table(cell_ids = colnames(input_epiclomal), count = NA) %>% as.matrix()

#for count_cpgs
for (c in seq_along(colnames(input_epiclomal))) { 
  count_cpgs[c,2] <- as.numeric(table(input_epiclomal[,c])[["0"]])
  count_cpgs[c,3] <- as.numeric(table(input_epiclomal[,c])[["1"]])
}
count_cpgs <- gather(count_cpgs, meth_status, ratio, unmeth:meth)

#for count_cpgs_2 - real cpgs totally
for (c in seq_along(colnames(input_epiclomal))) { 
  count_cpgs2[c,2] <- as.numeric(sum(!is.na(input_epiclomal[,c])))
}

library(ggpubr)
count_cpgs <- as.data.frame(count_cpgs)
ggbarplot(count_cpgs, "cell_ids", "ratio",
          fill = "meth_status", label = TRUE,
          x.text.angle = 90, width = 0.7)

#hist of real cpgs
count_cpgs2 <- as.data.frame(count_cpgs2)
count_cpgs2$count <- as.numeric(count_cpgs2$count)
gghistogram(count_cpgs2[,2], title = "Count of the real CpG obs")

#count distinct 0s and 1s 
dist_cells_1 <- data.table(cell_ids = colnames(input_epiclomal), dist_0 = NA, dist_1 = NA) %>% as.matrix()
for (c in c(1:300)) {
  tmp <- data.table(real = input_epiclomal[,c], imp = region_map[,cluster_map$cluster[c]])
  tmp <- tmp[complete.cases(tmp)]
  tmp$same <- ifelse(tmp$real != tmp$imp, TRUE, FALSE) 
  dist_var <- tmp$real[tmp$same=="TRUE"]
  dist_cells_1[c, 2] <- ifelse(length(table(dist_var)) > 0 , table(dist_var)[[1]], 0)
  dist_cells_1[c, 3] <- ifelse(length(table(dist_var)) > 1 , table(dist_var)[[2]], 0)
}

dist_cells <- as.data.frame(dist_cells)
dist_cells$dist_0 <- as.numeric(dist_cells$dist_0)
dist_cells$dist_1 <- as.numeric(dist_cells$dist_1)

#barplot of dist values to 0s
dist_cells <- gather(dist_cells, meth_status, count, dist_0:dist_1)
ggbarplot(dist_cells, "cell_ids", "count",
          fill = "meth_status", label = TRUE,
          x.text.angle = 90, width = 0.7)

#hist of distinct values
gghistogram(dist_cells_1$dist_0, title = "Distinguishable CpGs (from 1s to 0s)")
gghistogram(dist_cells_1$dist_1, title = "Distinguishable CpGs (from 0s to 1s)")

#table of proportion across cluster center for 2 type of annos (prom 1kb, gene1kb)
cluster_prop = data.table(annos = c(rep("prom1kb_10k",6), rep("prom1kb_20k",5), rep("gene1kb_10k", 5)),
                          cluster = c(1:6, 1:5, 1:5), prop_0 = NA, prop_1 = NA)

for (c in c(1:6)) {
  cluster_prop$prop_0[c] <- table(prom10k_clusters[,c])[[1]]
  cluster_prop$prop_1[c] <- table(prom10k_clusters[,c])[[2]]
}

#for each annos separately
prop_prom10k = data.table(annos = rep("prom1kb_10k",6), cluster = c(1:6), prop_0 = NA, prop_1 = NA)
for (c in c(1:6)) {
  prop_prom10k$prop_0[c] <- table(prom10k_clusters[,c])[[1]]
  prop_prom10k$prop_1[c] <- table(prom10k_clusters[,c])[[2]]
}   

prop_prom20k = data.table(annos = rep("prom1kb_20k",5), cluster = c(1:5), prop_0 = NA, prop_1 = NA)
for (c in c(1:5)) {
  prop_prom20k$prop_0[c] <- table(region_map[,c])[[1]]
  prop_prom20k$prop_1[c] <- table(region_map[,c])[[2]]
}     

prop_gene1kb = data.table(annos = rep("gene1kb_10k",5), cluster = c(1:5), prop_0 = NA, prop_1 = NA)
for (c in c(1:5)) {
  prop_gene1kb$prop_0[c] <- table(gene10k_clusters[,c])[[1]]
  prop_gene1kb$prop_1[c] <- table(gene10k_clusters[,c])[[2]]
}

####to add cpgs coords for cluster centersk####
library(stringi)

gene10k_clusters <- fread("/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/593/genotype_MAP.tsv.gz") %>%
  column_to_rownames(var = "cluster_id") %>%
  t() %>% as.data.frame()

colnames(gene10k_clusters) <- paste0("cl_", as.numeric(colnames(gene10k_clusters)) + 1)

gene10k_clusters$cpg <- colnames(input_CpG_data)

cpg_split <- stri_split_fixed(gene10k_clusters$cpg, ":")

gene10k_clusters$chr <- NA
gene10k_clusters$pos <- NA
for (i in seq_along(cpg_split)) {
  gene10k_clusters$chr[i] <- cpg_split[[i]][1]
  gene10k_clusters$pos[i] <- cpg_split[[i]][2]
}

gene10k_clusters$cpg <- NULL
gene10k_clusters <- gene10k_clusters[,c(7,8,1:6)]

for (i in c(3:8)) {
  tmp <- gene10k_clusters[,c(1,2,i)] 
  fwrite(tmp, sprintf("%s_%s.txt", "cluster_centers_gene1kb_10k",colnames(gene10k_clusters)[i]),sep='\t', col.names = FALSE)
  print(i)
}

#track the error - input= must be a single character string containing a file name, a system command containing at least one space, a URL starting 'http[s]://', 'ftp[s]://' or 'file://', or, the input data itself containing at least one \n or \r
#Calls: get_cov_data -> fread

#cell_based_methylation_extraction.R
args <- list()
args$output_directory = "/home/igor/Data/Epiclomal_results/14_01_26brain_88scNMT/cell_based_CpGs"
args$data_ID "26brain_88scNMT"
args$path_CpG_coordinates <- "/home/igor/Data/Epiclomal_results/14_01_26brain_88scNMT/CpG_coordinates_in_regions" 
args$path_cell_data <- "/home/igor/Data/sci-MET_mous_cortex_scNMT"
args$cell_ID <- "E6.5_Plate3_B1"
args$data_type <- "novoalign"
args$genome <- "mouse"
args$include_chrY <- 1


true[cell_id == "E4.5-5.5_new_Plate1_A09"] <- NULL

# test = subset(true, !(cell_id %in% c("E4.5-5.5_new_Plate1_A09.tsv.gz",
                                             "E4.5-5.5_new_Plate1_C02.tsv.gz",
                                             "E4.5-5.5_new_Plate1_E02.tsv.gz",
                                             "E4.5-5.5_new_Plate1_E09.tsv.gz",
                                             "E4.5-5.5_new_Plate1_G07.tsv.gz",
                                             "E7.5_Plate2_B3.tsv.gz")))
