#05/02 - finding how much sites were transferred to 0s and to 1s
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
rownames(region_map) <- seq(1:nrow(region_map))
region_map <- region_map[,-10]

#change colnames in region_map
colnames(region_map) <- paste0("cl_", as.numeric(colnames(region_map)) + 1)

result_comb <- input_epiclomal
for (c in c(1:nrow(cluster_map))) {
  tmp_table <- data.table(obs = input_epiclomal[,c], imp = region_map[,cluster_map$cluster[c]]) %>% 
    mutate(res = coalesce(obs, imp))
  result_comb[,c] <- tmp_table$res
  print(c)
}

#find the most variable cells within a cluster
cell_ids <- list()
for(c in c(1:ncol(region_map))) {
  cell_ids[[c]] <- cluster_map$cell_id[cluster_map$cluster==c] %>% as.character()
}
names(cell_ids) <- colnames(region_map)

##create cl center and cells for cl 1

#1
cl1_cells <- result_comb[,cell_ids_cl[[1]]]
cl1_cells <- cbind(region_map[,1], cl1_cells)
#2
table_diff_1 <- data.table(cell_id = cell_ids[[1]], diff_0 = NA , diff_1 = NA)
#3
for (f in seq_along(cell_ids[[1]])) {
  tmp <- cl1_cells[,c(1,f+1)] %>% as.data.frame()
  tmp$diff <- ifelse(tmp[,1] == tmp[,2], NA, tmp[,2])
  table_diff_1$diff_0[f] = table(tmp$diff)[["0"]]
  table_diff_1$diff_1[f] = table(tmp$diff)[["1"]]
}

##create cl center and cells for each cluster
#1
cl_cells <- list()
for (c in seq_along(cell_ids_cl)) {
  cl_cells[[c]] <- cbind(region_map[,c],result_comb[,cell_ids_cl[[c]]])
}
#2
table_diff <- list()
for (f in c(1:9)) {
  table_diff[[f]] <- data.table(cell_id = cell_ids[[f]], diff_0 = NA , diff_1 = NA)
}
#3
for (f in c(1:9)) {
  for (r in seq_along(cell_ids[[f]])) {
    tmp <- cl_cells[[f]][,c(1,r+1)] %>% as.data.frame()
    tmp$diff <- ifelse(tmp[,1] == tmp[,2], NA, tmp[,2])
    table_diff[[f]]$diff_0[r] = table(tmp$diff)[["0"]]
    table_diff[[f]]$diff_1[r] = table(tmp$diff)[["1"]]  
  }
}

#the process for finding most variable cpg across cells (?)
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
