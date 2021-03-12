#05/02 - finding how much sites were transferred to 0s and to 1s
library(tibble)
library(dplyr)
library(data.table)

best_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_SW_mm10_CGI_clust/0_0.95_10000/epi_region/462/"
input_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_sw_mm10_CGI_preproc/epiclomal_input/0_0.95_10000/"
out_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_SW_mm10_CGI_clust/0_0.95_10000/custom_visual/"

cluster_map <- fread(paste0(best_dir, "cluster_MAP.tsv.gz"), header = TRUE) %>%
  .[,c(1,2)] 
colnames(cluster_map) <- c("cell_id", "cluster")
cluster_map$cluster <- cluster_map$cluster + 1

setwd(input_dir)
input_epiclomal <- fread(paste0(input_dir, sprintf("input_Epiclomal_32_Sw.tsv.gz"))) %>%
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
for (c in seq_along(cell_ids)) {
  cl_cells[[c]] <- cbind(region_map[,c],result_comb[,cell_ids[[c]]])
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
    table_diff[[f]]$diff_0[r] = table(factor(tmp$diff, levels = 0:1))[[1]]
    table_diff[[f]]$diff_1[r] = table(factor(tmp$diff, levels = 0:1))[[2]]
  }
}

#count total dist 0 and 1
table_diff_total <- rbind(table_diff[[1]], table_diff[[2]], table_diff[[3]], table_diff[[4]], 
                          table_diff[[5]], table_diff[[6]], table_diff[[7]], table_diff[[8]], 
                          table_diff[[9]])

#the process for finding most variable cpg across cells (?)
cpgs_cells <- list()
for(c in c(1:ncol(region_map))) {
  cpgs_cells[[c]] <- as.matrix(result_comb[,cell_ids[[c]]])
}
names(cpgs_cells) <- colnames(region_map)


comp_tbl <- list()
for (t in c(1:ncol(region_map))) {
  comp_tbl[[t]] <- matrix(nrow = (ncol(cpgs_cells[[t]]) ** 2), ncol = 4)
}

for (r in c(1:7,9)) {
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
  comp_tbl[[r]]$V4 <- as.numeric(comp_tbl[[r]]$V4)
  setkey(comp_tbl[[r]], V4)
  
}

#for 5 cl separately
#to see count of diff cpgs between two cells
table(comp_tbl[[5]]$V4)

most_diff <- subset(comp_tbl[[5]], (V4 %in% c(253:273)))

sort(table(most_diff$V1))

cl_5_ids = c("E4.5-5.5_new_Plate2_A10", "E7.5_Plate4_H4", "E6.75_Plate2_G10", "E4.5-5.5_new_Plate2_D10","E4.5-5.5_new_Plate2_A06")

most_diff <- subset(most_diff, (V1 %in% cl_1_ids)) %>% 
  subset(., (V2 %in% cl_1_ids))

unique(most_diff$V3) %>% paste(., collapse = ',') %>% strsplit(., ",") %>% .[[1]] %>% as.numeric(.) %>% table(.)

#to take regions
result_comb_top1 <- result_comb[,cl_1_ids]
result_comb_top1 <- cbind(region_map[,1], result_comb_top1)
rownames(result_comb_top1) <- rownames(result_comb)


