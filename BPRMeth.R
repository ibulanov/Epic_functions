library(Melissa)
library(BPRMeth)
library(data.table)

# Set seed for reproducible results
set.seed(12) 
# Create basis object
basis_obj <- create_rbf_object(M = 4)

# Perform clustering
# Read annotation file, which will create annotation regions as well
anno_dt <- read_anno(file = "/home/igor/Data/zuse1db/melissa/anno_files/anno_prom_1k_37k.tsv", is_centre = TRUE)
cell_names <- names(input_melissa$met)

clustered_profiles_vb <- list()
for (i in c(1:77)) {
  # Read bulk methylation data
  met_cell <- read_met(file = sprintf("/home/igor/Data/zuse1db/melissa/input_cells_3_1/%s", cell_names[i]))
  # Create methylation regions that have CpG coverage
  met_region <- create_region_object(met_dt = met_cell, anno_dt = anno_dt)
  #cluster regions from one cell (VB)
  cluster_cell <- infer_profiles_mle(X = met_region$met, K = 1, model = "bernoulli", 
                                alpha_0 = .5, beta_0 = .1,
                                basis = basis_obj, vb_max_iter = 100, is_verbose = TRUE)
  clustered_profiles_vb[[i]] <- cluster_cell

}

vb_profiles_matrix <- matrix(data = NA, nrow = 5, ncol = 77)
for (i in c(1:77)) {
  vb_profiles_matrix[,i] <- clustered_profiles_vb[[i]]$W
}

colnames(vb_profiles_matrix) <- names(input_melissa$met)

#cluster regions from one cell (MLE)
cluster_310_mle <- cluster_profiles_mle(X = met_region$met, K = 2, model = "bernoulli", 
                                   alpha_0 = .5, beta_0 = .1,
                                   basis = basis_obj, vb_max_iter = 100, is_verbose = TRUE)


#umap of melissa 

w_1 <- melissa_output$W[,,1] %>% as.data.frame()
w_2 <- melissa_output$W[,,2] %>% as.data.frame()
w_1$V6 <- 1
w_2$V6 <- 2
w_main <- rbind(w_1, w_2)
w_data <- w_main[, c(1:5)]
w_label <- w_main[, 6]

#dimreduction
#PCA by factoMineR
library(FactoMineR)
res.pca = PCA(t(infer_prof_sc), scale.unit=TRUE, ncp=5, graph=T)

#UMAP
library(umap)
result_umap <- umap(w_data, config = umap.defaults, method = "umap-learn")
plot_umap(result_umap$layout, w_label)

plot_umap <- function(x, labels,
         main="UMAP visualization of the clustered weights (promoter 1kb, rbf = 4)",
         colors=c("chocolate4", "#17becf", "darkgoldenrod2", "chartreuse4"),
         pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=1) {

  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }

  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

#fill bprmeth result matrix by melissa output matrix
input_pbr_13k_list <- list()


filling_nas_2 <- function(melissa_output, infer_profiles_bases_4_mle) {
  for (q in seq_along(infer_profiles_bases_4_mle)[c(1:40)]) {
    melissa_cl2 <- melissa_output$W[,,2] %>% as.data.frame()
    melissa_cl2$status <- NA
    
    bpr_inf_prof_mle <- infer_profiles_bases_4_mle[[q]]$W %>% as.data.frame()
    bpr_inf_prof_mle$w6 <- "natural"
    
    for (i in seq_along(input_melissa$met[[q]])) {
      if (is.na(input_melissa$met[[q]][[i]])) {
        melissa_cl2[i, 6] <- "imputed"
      } else {
        melissa_cl2[i, 6] <- "natural"
      }
    }
    
    count = 1
    for (j in seq_along(melissa_cl2$status)) {
      if (melissa_cl2$status[j] == "natural") {
        melissa_cl2[j,] <- bpr_inf_prof_mle[count,]
        count = count + 1
      }
    }
    melissa_cl2$status <- NULL
    
    input_pbr_13k_list[[q]] <- melissa_cl2
    print(q)
    
  }
  
  return(input_pbr_13k_list)
}

filling_nas_1 <- function(melissa_output, infer_profiles_bases_4_mle) {
  for (q in seq_along(infer_profiles_bases_4_mle)[c(41:77)]) {
    melissa_cl1 <- melissa_output$W[,,1] %>% as.data.frame()
    melissa_cl1$status <- NA
    
    bpr_inf_prof_mle <- infer_profiles_bases_4_mle[[q]]$W %>% as.data.frame()
    bpr_inf_prof_mle$w6 <- "natural"
    
    for (i in seq_along(input_melissa$met[[q]])) {
      if (is.na(input_melissa$met[[q]][[i]])) {
        melissa_cl1[i, 6] <- "imputed"
      } else {
        melissa_cl1[i, 6] <- "natural"
      }
    }
    
    count = 1
    for (j in seq_along(melissa_cl1$status)) {
      if (melissa_cl1$status[j] == "natural") {
        melissa_cl1[j,] <- bpr_inf_prof_mle[count,]
        count = count + 1
      }
    }
    melissa_cl1$status <- NULL
    
    input_pbr_13k_list[[q]] <- melissa_cl1
    print(q)
    
  }
  
  return(input_pbr_13k_list)
}

#replace mean weights from melissa by profiles weights for the second cluster cells
result_2 <- filling_nas_2(melissa_output, infer_profiles_bases_4_mle = infer_profiles_cells_base4_mle)
#replace mean weights from melissa by profiles weights for the first cluster cells
result_1 <- filling_nas_1(melissa_output, infer_profiles_bases_4_mle = infer_profiles_cells_base4_mle)

#create main result list of df
for (i in seq_along(result_1)[c(1:40)]) {
  result_1[[i]] <- result_2[[i]]
}


#reshape dfs to one columns
for (i in seq_along(result_1)) {
  result_1[[i]] <- tidyr::gather(result_1[[i]])
  result_1[[i]]$key <- NULL
}

#convert values to another class
for (i in seq_along(result_1)) {
  result_1[[i]] <- as.data.frame(result_1[[i]])
}

#create matrix
infer_prof_sc <- matrix(nrow = 69955, ncol = 77) %>% as.data.frame()


for (i in seq_along(result_1)) {
  infer_prof_sc[,i] <- result_1[[i]]
}

#assign the names
names(infer_prof_sc) <- names(input_melissa$met)

####process real data (scNMT-seq)####
input_data <- fread("/home/igor/Data/scnmt_gastrulation/met/cpg_level/E4.5-5.5_new_Plate1_A02.tsv.gz")

#modify annotation file
anno <- fread("/home/igor/Data/zuse1db/melissa/anno_files/anno_prom_1k_37k.tsv")

replace_chrom_names <- function(input) {
  
  colnames(input)[1] <- "chrom"
  input$chrom <- ("chr", "", input$chrom)
  input$chrom[input$chrom == "Un_random"] <- NA
  input$chrom[input$chrom == "1_random"] <- NA
  input$chrom[input$chrom == "3_random"] <- NA
  input$chrom[input$chrom == "4_random"] <- NA
  input$chrom[input$chrom == "5_random"] <- NA
  input$chrom[input$chrom == "7_random"] <- NA
  input$chrom[input$chrom == "8_random"] <- NA
  input$chrom[input$chrom == "9_random"] <- NA
  input$chrom[input$chrom == "X_random"] <- NA
  input$chrom[input$chrom == "Y_random"] <- NA
  input$chrom[input$chrom == "13_random"] <- NA
  input$chrom[input$chrom == "16_random"] <- NA
  input$chrom[input$chrom == "17_random"] <- NA
  input$chrom[input$chrom == "M"] <- "MT"
  
  input$chrom[input$chrom == "1_GL456210_random"] <- NA
  input$chrom[input$chrom == "1_GL456211_random"] <- NA
  input$chrom[input$chrom == "1_GL456212_random"] <- NA
  input$chrom[input$chrom == "1_GL456221_random"] <- NA
  input$chrom[input$chrom == "4_GL456216_random"] <- NA
  input$chrom[input$chrom == "4_GL456350_random"] <- NA
  input$chrom[input$chrom == "4_JH584292_random"] <- NA
  input$chrom[input$chrom == "4_JH584293_random"] <- NA
  input$chrom[input$chrom == "4_JH584294_random"] <- NA
  input$chrom[input$chrom == "4_JH584295_random"] <- NA
  input$chrom[input$chrom == "5_GL456354_random"] <- NA
  input$chrom[input$chrom == "5_JH584296_random"] <- NA
  input$chrom[input$chrom == "5_JH584297_random"] <- NA
  input$chrom[input$chrom == "5_JH584298_random"] <- NA
  input$chrom[input$chrom == "5_JH584299_random"] <- NA
  input$chrom[input$chrom == "7_GL456219_random"] <- NA
  input$chrom[input$chrom == "Un_GL456239"] <- NA
  input$chrom[input$chrom == "Un_GL456372"] <- NA
  input$chrom[input$chrom == "Un_GL456381"] <- NA
  input$chrom[input$chrom == "Un_GL456385"] <- NA
  input$chrom[input$chrom == "Un_JH584304"] <- NA
  input$chrom[input$chrom == "X_GL456233_random"] <- NA
  input$chrom[input$chrom == "Y_JH584303_random"] <- NA
  
  input <- input[complete.cases(input), ]
  
  return(input)
  
}

anno <- replace_chrom_names(anno)

fwrite(anno, "/home/igor/Data/melissa/annos_mm10_prom5kb.tsv", sep='\t')


#process the scnmt cells
scnmt_cells <- list.files(path= "/home/igor/Data/scnmt_gastrulation/met/cpg_level", pattern = "*.tsv.gz", full.names = FALSE)

list_of_scnmt_cells <- list()

for (i in seq_along(filtered_cells)) {
  cur_cell <- fread(sprintf("/home/igor/Data/scnmt_gastrulation/met/cpg_level/%s.tsv.gz", filtered_cells[i]))
  list_of_scnmt_cells[[i]] <- cur_cell[,c(1,2,5)]
  colnames(list_of_scnmt_cells[[i]]) <- c("chrom", "pos", "met_level")
  names(list_of_scnmt_cells)[i] <- filtered_cells[i]
  print(i)
}

list_of_scnmt_cells[sapply(list_of_scnmt_cells, is.null)] <- NULL

for (i in seq_along(list_of_scnmt_cells)) {
  fwrite(list_of_scnmt_cells[[i]], 
         sprintf("/home/igor/Data/melissa/scnmt_same_filter/%s.tsv.gz", names(list_of_scnmt_cells)[i]),
         sep = '\t')
  print(i)
}

input_melissa <- create_melissa_data_obj(
  met_dir = "/home/igor/Data/melissa/scnmt_same_filter",
  anno_file = "/home/igor/Data/melissa/anno_mm10_18k_prom1kb.tsv",
  chrom_size_file = NULL,
  chr_discarded = NULL,
  is_centre = TRUE,
  is_window = FALSE,
  upstream = -5000,
  downstream = 5000,
  cov = 5,
  sd_thresh = -1,
  no_cores = 5
) 

input_melissa <- filter_by_cpg_coverage(input_melissa, min_cpgcov = 3)
input_melissa <- filter_by_variability(input_melissa, min_var = 0.05)
input_melissa <- filter_by_coverage_across_cells(input_melissa, min_cell_cov_prcg = 0.05)

#run melissa
# Set seed for reproducible results
set.seed(15) 
# Create basis object
basis_obj <- create_rbf_object(M = 4)


melissa_obj <- melissa(X = input_melissa$met, K = 3, basis = basis_obj,
                               vb_max_iter = 500, vb_init_nstart = 10, 
                               is_parallel = TRUE, no_cores = 2, is_verbose = TRUE)

#Run bprmeth for each cell 
run_bpr_cells <- function() {
  #create list of inputs for bpr
  list_of_input_objs <- list()
  for (i in seq_along(input_melissa$met)) {
    cur_cell <- list(met = input_melissa$met[[i]],
                     anno = input_melissa$anno_region)
    list_of_input_objs[[i]] <- cur_cell
  }
  
  #delete extra na 
  for (i in seq_along(list_of_input_objs)) {
    list_of_input_objs[[i]]$met <- list_of_input_objs[[i]]$met[!is.na(list_of_input_objs[[i]]$met)]
  }
  
  #run bprmeth for each cell
  result <- list()
  for (i in seq_along(list_of_input_objs)) {
    tmp_cell <- infer_profiles_mle(X = list_of_input_objs[[i]]$met, model = "bernoulli",
                                   basis = basis_obj, 
                                   vb_max_iter = 500,
                                   is_parallel = TRUE,
                                   no_cores = 5,
                                   is_verbose = TRUE
    )
    result[[i]] <- tmp_cell
    print(i)
  }
  return(result)
}

bpr_cells <- run_bpr_cells()


scnmt_filling_nas <- function(input_melissa, output_bpr) {
  
  both_packs_list <- list()
  
  names(output_bpr) <- names(input_melissa$met)
  
  cells_cl <- as.data.frame(melissa_obj$r_nk)
  cells_cl$cluster1[cells_cl$cluster1 >= 0.5] <- "cluster1"
  cells_cl$cluster2[cells_cl$cluster2 >= 0.5] <- "cluster2"
  cells_cl$cluster3[cells_cl$cluster3 >= 0.5] <- "cluster3"
  
  cells_cl$cluster1[cells_cl$cluster1 != "cluster1"] <- NA
  cells_cl$cluster2[cells_cl$cluster2 != "cluster2"] <- NA
  cells_cl$cluster3[cells_cl$cluster3 != "cluster3"] <- NA
  cells_cl <- cells_cl %>% tidyr::unite("result", cluster1:cluster3, na.rm = TRUE, remove = FALSE)
  cells_cl[,c(2:4)] <- NULL
  
  
  for (q in seq_along(output_bpr)) {
    tmp_df <- as.data.frame(output_bpr[[q]]$W)
    tmp_df$w6 <- "natural"
    output_bpr[[q]]$W <- tmp_df
    
    if (cells_cl$result[q] == "cluster1") {
      melissa_cl <- as.data.frame(melissa_obj$W[,,1])
    } else if (cells_cl$result[q] == "cluster2") {
      melissa_cl <- as.data.frame(melissa_obj$W[,,2])
    } else if (cells_cl$result[q] == "cluster3") {
      melissa_cl <- as.data.frame(melissa_obj$W[,,3])
    }
    
    for (i in seq_along(input_melissa$met[[q]])) {
      if (is.na(input_melissa$met[[q]][[i]])) {
        melissa_cl[i, 6] <- "imputed"
      } else {
        melissa_cl[i, 6] <- "natural"
      }
    }
    
    count = 1
    for (j in seq_along(melissa_cl$V6)) {
      if (melissa_cl$V6[j] == "natural") {
        melissa_cl[j,] <- output_bpr[[q]]$W[count,] %>% as.data.frame()
        count = count + 1
      }
    }
    melissa_cl$V6 <- NULL
    
    both_packs_list[[q]] <- melissa_cl
    print(q)
    
  }
  
  return(both_packs_list)
}

result <- scnmt_filling_nas(
  input_melissa = input_melissa, 
  output_bpr = output_bpr)

#reshape dfs to one columns
for (i in seq_along(result)) {
  result[[i]] <- tidyr::gather(result[[i]])
  result[[i]]$key <- NULL
  result[[i]] <- as.data.frame(result[[i]])
}

#create matrix
infer_prof_sc <- matrix(nrow = 35485, ncol = 986) %>% as.data.frame()

for (i in seq_along(result)) {
  infer_prof_sc[,i] <- result[[i]]
}

#assign the names
names(infer_prof_sc) <- NULL
dev_fraction <- names(input_melissa$met) 

sample_metadata <- fread("/home/igor/Data/scnmt_gastrulation/sample_metadata.txt") %>%
  .[, c("id_met", "pass_metQC", "stage")] %>%
  .[complete.cases(.),] %>% .[,subset(., pass_metQC == "TRUE")] %>%
  .[, pass_metQC := NULL] %>% setorder()

infer_prof_sc <- t(infer_prof_sc) %>% as.data.frame()

infer_prof_sc$V35486 <- sample_metadata$stage


#PCA by factoMineR
library(FactoMineR)
res.pca <- PCA(infer_prof_sc, quali.sup=22001, 
              scale.unit=TRUE, ncp=5, graph=FALSE) 
plot.PCA(res.pca, label = "none", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#307332")) 


#UMAP of cluster weights
cluster_1 <- melissa_obj$W[,,1] %>% as.data.frame()
cluster_2 <- melissa_obj$W[,,2] %>% as.data.frame()
cluster_3 <- melissa_obj$W[,,3] %>% as.data.frame()
cluster_4 <- melissa_obj$W[,,4] %>% as.data.frame()

cluster_1$label <- 1
cluster_2$label <- 2
cluster_3$label <- 3
cluster_4$label <- 4

umap_matrix <- rbind(cluster_1, cluster_2, cluster_3, cluster_4)

labels <- umap_matrix[,6]
clusters_data <- umap_matrix[,c(1:5)]

library(umap)
result_umap <- umap(clusters_data, config = umap.defaults, method = "umap-learn")
plot_umap(result_umap$layout, labels)

#barplot of bpr real genomic regions
bpr_regions <- data.frame(cell = rep(NA, 986), bpr = rep(NA, 986))
for (i in seq_along(list_of_input_objs)) {
  bpr_regions$cell[i] <- i   
  bpr_regions$bpr[i] <- length(list_of_input_objs[[i]]$met)
}

barplot(bpr_regions$bpr,  xlab = "cells", ylab = "Number of met profiles by BPRMeth", names.arg = bpr_regions$cell)
