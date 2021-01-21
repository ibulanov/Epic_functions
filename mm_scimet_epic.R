# 13_01_21
# Use sci_MET in-house and published data for Epiclomal

library(data.table)
library(dplyr)
library(RMySQL)
library(tidyr)
library(stringr)

#prepare mouse cortex data
setwd("/home/igor/Data/Mulqueen2018_mm10_100")
mm_cells <- list.files(".")

cell_ids = strsplit(mm_cells, "_", fixed = TRUE) %>% unlist() %>% matrix(nrow = 100, byrow = TRUE) %>% as.data.frame()
cell_ids = paste0(cell_ids$V3, "_", cell_ids$V4)

dirs = list.dirs("/home/igor/Data/Mulqueen2018_mm10_100")[2:101]

for (d in seq_along(dirs)) {
  setwd(dirs[d])
  file.list = list.files(".")
  dataset <- do.call("rbind", lapply(file.list, FUN = function(file) {
    fread(file, header=TRUE, sep="\t")})) %>% setDT() %>% .[, c("strand", "mc_class") := NULL]
  
  colnames(dataset) = c("chr","CpG_start", "count_meth", "total_reads", "meth_frac")
  dataset$CpG_end = dataset$CpG_start + 1
  dataset$chr = paste0("chr", dataset$chr)
  dataset = dataset[,c(1,2,6,5,3,4)]
  fwrite(dataset, sprintf("/home/igor/Data/Mulqueen2018_mm10_100/%s.tsv.gz", cell_ids[d]), sep='\t')
  print(d)
}

#restructure data mouse cortex data - old conde for wrong code above deleted already.
setwd("/home/igor/Data/sci-MET_mous_cortex_scNMT")
for (f in c(1:26)) {
  tmp <- fread(all_files[f], sep=',')
  tmp <- tmp[,c(1:3,6,5,4)]
  colnames(tmp) <- c("chr", "CpG_start", "CpG_end", "meth_frac", "total_reads", "count_meth")
  fwrite(tmp, all_files[f], sep='\t')
  print(f)
}

for (f in c(1:114)) {
  tmp <- fread(mm_cells[f], sep='\t')
  tmp <- tmp[,c(1,2,3,4,6,5)]
  fwrite(tmp, mm_cells[f], sep='\t')
  print(f)
}
  
#for scNMT data
setwd("/home/igor/Data/scnmt_gastrulation/met/epi_100_high")

mm_cells <- list.files(".")
for (f in c(1:100)) {
  setwd("/home/igor/Data/scnmt_gastrulation/met/epi_100_high")
  tmp <- fread(mm_cells[f])
  tmp$total_reads <- tmp$met_reads	+ tmp$nonmet_reads
  colnames(tmp) <- c("chr", "CpG_start", "count_meth", "nonmet_reads", "meth_frac", "total_reads")
  tmp$CpG_end <- tmp$CpG_start + 1
  tmp <- tmp[,c(1,2,7,5,3,6)]
  tmp$chr = paste0("chr", tmp$chr)
  setwd("/home/igor/Data/sci-MET_mouse_cortex_scNMT/new_high")
  fwrite(tmp, sprintf("%s", mm_cells[f]), sep='\t')
  print(f)
}

# creating the true file
all_files =list.files("/home/igor/Data/sci-MET_mouse_cortex_scNMT")
true_file = data.table(cell_id = all_files, epigenotype_id = NA)
true_file$epigenotype_id[1:100] = 1
true_file$epigenotype_id[101:200] = 2
fwrite(true_file, "/home/igor/Data/sci-MET_mous_cortex_scNMT_true.tsv.gz", sep='\t')

#16/01/21
#Guo2017 20 cells for epiclomal
setwd("/home/igor/Data/Zhu_Guo2017_100/smallest_50/30")
hg_cells = list.files("/home/igor/Data/Zhu_Guo2017_100/smallest_50/30")

cell_ids = strsplit(hg_cells, "_", fixed = TRUE) %>% unlist() %>% matrix(nrow = 30, byrow = TRUE) %>% as.data.frame() 
cell_ids = strsplit(cell_ids$V2, ".", fixed = TRUE) %>% unlist() %>% matrix(nrow = 30, byrow = TRUE) %>% as.data.frame() %>% .[,1]

for (f in seq_along(hg_cells)) {
  tmp <- fread(hg_cells[f])
  tmp <- tmp[tmp$Type == "CpG"]
  tmp$Pos <- ifelse(tmp$Chain == "+", tmp$Pos + 1, tmp$Pos) #right
  tmp$MetRate <- ifelse(tmp$MetRate <= 0.5, 0, 1)
  tmp$CpG_end <- tmp$Pos
  tmp <- tmp[,c(1,2,11,8,6,5)]
  colnames(tmp) <- c("chr","CpG_start","CpG_end","meth_frac","count_meth","total_reads")
  fwrite(tmp, sprintf("/home/igor/Data/Zhu_Guo2017_100/smallest_50/proc_for_epi/%s_cov.tsv.gz", cell_ids[f]), sep='\t')
  print(f)
}

#true for 320 KC and Zhu
cells_350 <- list.files("/home/igor/Data/mysql_data/KC_2k/cov_300KC_all_chr_20_Zhu") %>%
  gsub(".tsv.gz", "", .) %>%
  gsub("_cov", "", .)

true350 <- data.table(cell_id = cells350, epigenotype_id = c(rep(1, 300), rep(2, 50)))
fwrite(true350, "/home/igor/Data/Epiclomal_results/true_files/350KC_Zhu_true.txt.gz", sep="\t")
# 


#UMAP

####21_01_21 UMAP for mm9 and hg19 datasets####
library(tibble)
library(dplyr)
library(data.table)

best_dir <- "/home/igor/Data/Epiclomal_results/17_01_320_KC_Zhu_hg19_clust/5_0.99_10000/epi_region/9/"
input_dir <- "/home/igor/Data/Epiclomal_results/16_01_320_KC_Zhu_hg19_preproc/epiclomal_input/5_0.99_10000/"
output_dir <- "/home/igor/Data/Epiclomal_results/17_01_320_KC_Zhu_hg19_clust/5_0.99_10000/visual/"

cluster_map <- fread(paste0(best_dir, "cluster_MAP.tsv.gz"), header = TRUE) %>%
  .[,c(1,2)] 
colnames(cluster_map) <- c("cell_id", "cluster")
cluster_map$cluster <- cluster_map$cluster + 1

setwd(input_dir)
input_epiclomal <- fread(paste0(input_dir, "input_Epiclomal_320KC_Zhu.tsv.gz")) %>%
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
M3C::umap(result_comb,
          labels=as.factor(cluster_map$cluster),
          dotsize = 2,
          axistextsize = 10,
          legendtextsize = 15
          )







