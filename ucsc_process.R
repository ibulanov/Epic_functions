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

cluster_map <- fread("/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/593/cluster_MAP.tsv.gz") %>%
  .[,c(1,2)] 
cluster_map$cluster <- cluster_map$cluster + 1
colnames(cluster_map) <- c("cell_id", "cluster")

input_epiclomal <- fread("/home/igor/Data/Epiclomal_results/30_12_KC300_20k_gene1kb_preproc/epiclomal_input/5_0.99_10000/input_Epiclomal_KC300_20k_gene1kb.tsv.gz") %>%
  column_to_rownames(., "cell_id") %>%
  t()


region_map <- fread("/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/593/genotype_MAP.tsv.gz") %>% 
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

fwrite(result_comb, "/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/visual/epiclomal_comb_5_99_10k_gene1kb.tsv",
       sep='\t',
       col.names = TRUE)

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

most_diff <- subset(comp_tbl[[1]], (V4 %in% c(10:17)))

sort(table(most_diff$V1))

most_diff <- subset(most_diff, (V1 %in% c(10934,11467,11494,11399,10300))) %>% 
  subset(., (V2 %in% c(10934,11467,11494,11399,10300)))

cpgs_index <- unique(most_diff$V3) %>% paste(., collapse = ',') %>% strsplit(., ",") %>% .[[1]] %>% as.numeric(.) %>% table(.)

#create table of most variable cells and regions

cpg <- rownames(result_comb)[c(10368,13063,18733,21574,21657,2694,6986,12378,16508,21719,2162,7043,13941,9207,21740,2132,3249,18877,19427,21722,2863,6854,12558,13832,18569,21601)]

library(stringr)
cpgs <- data.table(cpg = cpg,
                   chrom = str_split(cpg, ":", simplify = TRUE)[,1],
                   pos = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric(),
                   pos_d = str_split(cpg, ":", simplify = TRUE)[,2] %>% as.numeric()
)

filter_regions <- fread("/home/igor/Data/Epiclomal_results/30_12_KC300_20k_gene1kb_preproc/filter_regions/5_0.99_10000/filtered_regions_KC300_20k_gene1kb.tsv", col.names = "V1")

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
epiclomal_best_cells <- epiclomal_best_cells[, c(28:36, 1:27)]

desired_regions <- overlap_cpgs_regions[,c(5,1,6,7)] %>% distinct(., .keep_all = TRUE)

#merge with desired list of regions
setkey(epiclomal_best_cells, chrom, start, end)
setkey(desired_regions, chrom, start,end)
epiclomal_best_cells_regions <- foverlaps(epiclomal_best_cells, desired_regions, nomatch = NULL)
epiclomal_best_cells_regions <- epiclomal_best_cells_regions[,c(1,6,7:39)]
colnames(epiclomal_best_cells_regions)[1:8] <- c("chrom", "start", "end", "cl_1", "cl_2", "cl_3", "cl_4", "cl_5")

fwrite(epiclomal_best_cells_regions, "/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/visual/epiclomal_best_cells_regions.tsv",
       sep='\t',
       col.names = TRUE)

epiclomal_best_cells_regions <- as.matrix(epiclomal_best_cells_regions)

setwd("/home/igor/Data/Epiclomal_results/05_01_KC300_5_99_10k_gene1kb_clustering/visual/cls_cells_ucsc")
count = 1
for (i in c(4:35)) {
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
