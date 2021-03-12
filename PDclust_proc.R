library(PDclust)
library(tibble)
library(data.table)
library(dplyr)
library(magrittr)

#process Smallwood data

#1
setwd("/home/igor/Data/Mulqueen2018/hg19/Mulqueen_sci-MET_hg19_Epiclomal/")
lf = list.files(".")

cpg_files <- list()
for (f in seq_along(lf)) {
  tmp <- fread(lf[f])
  tmp <- tmp[,c(1:4)]
  tmp$meth_frac <- tmp$meth_frac * 100
  colnames(tmp) <- c("chr","start","end", "meth")
  tmp$end <- tmp$end + 1
  tmp <- as_tibble(tmp)
  tmp$end <- as.integer(tmp$end)
  tmp$meth <- as.integer(tmp$meth)
  cpg_files[[f]] <- tmp
  print(f)
}

names(cpg_files) <- c(1:32)

#2
cpg_files_pairwise <- create_pairwise_master(sw_cells, cores_to_use = 5)

#3
cpg_files_pairwise_matrix <- convert_to_dissimilarity_matrix(cpg_files_pairwise)

#4
cluster_results <- cluster_dissimilarity(cpg_files_pairwise_matrix, num_clusters = 2)

plot(cluster_results$hclust_obj)

test = cluster_results$cluster_assignments
test = setDT(test, keep.rownames = TRUE)
test$rn <- as.numeric(test$rn)
setkey(test, rn)


#Mulqueen dataset

#1
setwd("/home/igor/Data/Mulqueen2018/hg19/Mulqueen_sci-MET_hg19_Epiclomal/")
lf = list.files(".")

cpg_files <- list()
for (f in seq_along(lf)) {
  tmp <- fread(lf[f])
  tmp <- tmp[,c(1:4)]
  tmp$meth_frac <- tmp$meth_frac * 100
  colnames(tmp) <- c("chr","start","end", "meth")
  tmp$end <- tmp$end + 1
  tmp <- as_tibble(tmp)
  tmp$end <- as.integer(tmp$end)
  tmp$meth <- as.integer(tmp$meth)
  cpg_files[[f]] <- tmp
  print(f)
}

names(cpg_files) <- c(1:315)
#2
cpg_files_pairwise <- create_pairwise_master(cpg_files, cores_to_use = 5)
#3
cpg_files_pairwise_matrix <- convert_to_dissimilarity_matrix(cpg_files_pairwise)
#4
cluster_results <- cluster_dissimilarity(cpg_files_pairwise_matrix, num_clusters = 2)

#plot
plot(cluster_results$hclust_obj, cex.lab = 0.7, cex.axis = 0.7, cex.sub = 0.7)

#custom label font size
plot(cluster_results$hclust_obj, xlab="xlab", ylab="ylab", main="main", sub="")

# reduced label size
par(cex=0.3, mar=c(5, 8, 4, 1))
plot(cluster_results$hclust_obj, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(xlab="xlab", ylab="ylab", main="main")
axis(2)

library(pheatmap)
heatmap_pallete <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21)

pheatmap(cpg_files_pairwise_matrix,
         cluster_rows = cluster_results$hclust_obj,
         cluster_cols = cluster_results$hclust_obj,
         treeheight_row = 0,
         border_color = NA,
         color = heatmap_pallete,
         show_colnames = F,
         annotation_col = cluster_results$cluster_assignments,
         fontsize_row = 3)


#create a table as after euclideanClust from Epiclomal
test <- cluster_results$cluster_assignments
test <- setDT(test, keep.rownames = TRUE)
colnames(test) <- c("cell_id", "EuclideanClust_region_num_clusters_1")

for (k in c(2:Max_K)) {
  cluster_results <- cluster_dissimilarity(cpg_files_pairwise_matrix, num_clusters = k)
  tmp_r <- cluster_results$cluster_assignments
  test <- cbind(test,tmp_r)
  colnames(test)[k+1] <- sprintf("EuclideanClust_region_num_clusters_%s", k)
  test$EuclideanClust_best_cluster_2 <- test$EuclideanClust_region_num_clusters_2
}

cluster_results <- cluster_dissimilarity(cpg_files_pairwise_matrix, num_clusters = 2)
test = cluster_results$cluster_assignments
test = setDT(test, keep.rownames = TRUE)
test$EuclideanClust_best_cluster_2 = test$cluster


#test$rn <- as.numeric(test$rn)
#setkey(test, rn)

#another approach to use not default input data, but processed by Epiclomal
library(magrittr)

load("/home/igor/Data/Epiclomal_results/315_Mulq/ERB/24_02_315Mulq_All_ERB245/0_0.99_20000/simple_hclust/CpG_based_imputed.RDa.gz")
input_epiclomal <- t(input_CpG_data)

mulq_cpg_files <- list()

for(c in c(1:ncol(input_epiclomal))) {
  tmp <- data.table(chr = rep(NA, nrow(input_epiclomal)),
                    start = rep(NA, nrow(input_epiclomal)),
                    end = rep(NA, nrow(input_epiclomal)),
                    meth = rep(NA, nrow(input_epiclomal))
                    )
  tmp$chr <-  strsplit(rownames(input_epiclomal), ":")  %>% sapply(extract2, 1) %>% paste0("chr", .)
  tmp$start <- strsplit(rownames(input_epiclomal), ":")  %>% sapply(extract2, 2) %>% as.numeric()
  tmp$end <- tmp$start + 1
  tmp$meth <- input_epiclomal[,c]
  tmp <- tmp[complete.cases(tmp),]
  mulq_cpg_files[[c]] <- tmp
  print(c)
}

names(mulq_cpg_files) <- c(1:315)
