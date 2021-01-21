library(Melissa)
library(BPRMeth)
library(factoextra)
library(RMySQL)
library(dplyr)
library(data.table)

#args
args = commandArgs(trailingOnly=TRUE)
#1 - cell files directory (for storing from db and then this directory will be used for Melissa data object creation)
#2 - annotation file directory - the file with be created from db - the name has to be specified clearly: "/path/to/file/custom_annotations.tsv"
#3 - Number of clusters for Melissa running
#4 - Number of cores for Melissa and BPRMeth
#5 - scNMT sample metadata file directory
#6 - Database name, "mouse_public" or "scm_mm10"
#7 - desired cell IDs (ranges, numbers, e.g. "100310:100384, 100386" or "all")
#8 - To choose desired chr (e.g. either only 1 chr or all, e.g. 'chr1' or 'all') 
#9 - Name of output PCA plot
#10 - promoters length for annotation creating
#11 - Number of RBF
#12 - Directory for the intermediate results


#the custom function for inferring the mean value instead of plot
infer_melissa_mean <- function(melissa_obj, region=1, cluster=1) {
  
  W_Sigma <- list()
  for (cl in seq_along(melissa_obj$pi_k)) {
    W_Sigma[[cl]] <- melissa_obj$W_Sigma[[cl]][[region]]
  }
  cluster_obj <- list(W = melissa_obj$W[region, , ], W_Sigma = W_Sigma, 
                      basis = melissa_obj$basis)
  class(cluster_obj) <- c("cluster_profiles_vb_bernoulli", "cluster_profiles_vb")
  
  aes_xs <- seq(from = -1, to = 1, by = 0.01)
  K <- NCOL(cluster_obj$W)
  ys <- matrix(0, ncol = K, nrow = length(aes_xs))
  aes_ys = x = y = Cluster <- NULL
  ys_low = ys_high <- NULL
  if (methods::is(cluster_obj, "cluster_profiles_mle")) {
    if (methods::is(cluster_obj, "cluster_profiles_mle_gaussian")) {
      for (k in 1:K) {
        ys[, k] <- eval_function(cluster_obj$basis, 
                                 aes_xs, cluster_obj$W[, k])
      }
    } else {
      for (k in 1:K) {
        ys[, k] <- eval_probit_function(cluster_obj$basis, 
                                        aes_xs, cluster_obj$W[, k])
      }
    }
  } else if (methods::is(cluster_obj, "cluster_profiles_vb")) {
    tmp <- BPRMeth:::.predictive_cluster_profile(cluster_obj, aes_xs)
    ys <- tmp$W_pred
    if (methods::is(cluster_obj, "cluster_profiles_vb_binomial") || 
        methods::is(cluster_obj, "cluster_profiles_vb_bernoulli")) {
      ys_low <- ys - ys * (1 - ys)
      ys_high <- ys + ys * (1 - ys)
    } else if (methods::is(cluster_obj, "cluster_profiles_vb_gaussian")) {
      ys_low <- ys - 2 * tmp$W_sd_pred
      ys_high <- ys + 2 * tmp$W_sd_pred
    }
  } else {
    stop("No plotting function for this model!")
  }
  
  dt <- data.table::data.table(aes_xs = numeric(), aes_ys = numeric(), 
                               ys_low = numeric(), ys_high = numeric(), Cluster = numeric())
  if (methods::is(cluster_obj, "cluster_profiles_vb") || methods::is(cluster_obj, 
                                                                     "cluster_profiles_gibbs")) {
    for (k in 1:K) {
      dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs, 
                                             aes_ys = ys[, k],
                                             ys_low = ys_low[, k], 
                                             ys_high = ys_high[, k], 
                                             Cluster = as.factor(k)))
    }
  } else {
    for (k in 1:K) {
      dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs, 
                                             aes_ys = ys[, k],
                                             ys_low = 0,
                                             ys_high = 0,
                                             Cluster = as.factor(k)))
    }
  }
  
  
  clust_met_mean <- mean(dt$aes_ys[dt$Cluster==cluster])
  
  
  
  return(clust_met_mean)
  
}

#connect to the db
dbconn <- dbConnect(MySQL(), user='i747v', password='igor',
                     dbname=args[6], host='zuse1db')


inhouse_cpgs_proc <- function(dbconn) {
  
  if (args[7] == "all") {
    result_seq = ""
  } else {
    #convert string ranges to sql query
    input_string <- args[7]
    spl_str <- strsplit(input_string, split=',', fixed=FALSE)
    colon <- strsplit(spl_str[[1]], split = ":", fixed = TRUE)
    cell_seq <- c()
    for (ee in colon) {
      if (length(ee) > 1){
        tmp_v <- seq(ee[1], ee[2])
        cell_seq <- c(cell_seq, tmp_v)
      } else {
        tmp_v <- as.integer(ee)
        cell_seq <- c(cell_seq, tmp_v)
      }
    }
    
    cell_seq <- toString(cell_seq)
    result_seq <- paste0(sprintf("where C.expID in (%s)", cell_seq))
  }
  
  if (args[8]=="all") {
    chrom <- ""
  } else {
    chrom <- sprintf("where B.chrom in ('%s')", args[8])
  }
  
  if (chrom != "" && result_seq != "") {
    and <- "and"
    result_seq <- gsub("where", "", result_seq)
    result_seq <- gsub(" C.expID", "C.expID", result_seq)
  } else {
    and <- ""
  } 
  
  cat("Performing SQL query...")
  
  #to get the input data (cpgBases)
  input_cpgs <-  dbSendQuery(dbconn, sprintf("select C.expID, B.chrom, B.pos, B.strand, C.ratio
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        %s
                        %s
                        %s", chrom, and, result_seq)) %>% fetch(n=-1) 
  
  input_cpgs <- setDT(input_cpgs) %>% .[, ratio:= ifelse(ratio >= 0.5, 1, 0)] %>% 
    .[, compl_pos := ifelse(strand == 1, pos+1, pos-1)]
  
  #create the complementary df
  cpgs_compl <- data.table(input_cpgs$expID, input_cpgs$chrom, input_cpgs$compl_pos, input_cpgs$ratio)
  
  #remove extra columns
  input_cpgs <- input_cpgs[,c("strand","compl_pos"):= NULL] 
  
  colnames(cpgs_compl) <- colnames(input_cpgs)
  
  #concatenate the two dataframes and order the rows
  input_cpgs <- rbind(input_cpgs, cpgs_compl) %>% arrange(expID) %>% distinct() 
  
  #delete cpgbases complementary data table
  cpgs_compl <- NULL
  
  cat(sprintf("%s CpG sites were found\n", nrow(input_cpgs)))
  
  input_cpgs <- split(input_cpgs, input_cpgs$expID)
  
  for (i in seq_along(input_cpgs)) {
    fwrite(input_cpgs[[i]][,c(2:4)], paste0(args[[1]], names(input_cpgs)[i], ".tsv.gz"), sep = '\t')
  }
  
  cat(sprintf("%s cells were written as separate files to the directory: %s \n", length(input_cpgs), args[1]))
  
}

scnmt_cpgs_proc <- function(dbconn) {
  
  if (args[7] == "all") {
    result_seq = ""
  } else {
    #convert string ranges to sql query
    input_string <- args[7]
    spl_str <- strsplit(input_string, split=',', fixed=FALSE)
    colon <- strsplit(spl_str[[1]], split = ":", fixed = TRUE)
    cell_seq <- c()
    for (ee in colon) {
      if (length(ee) > 1){
        tmp_v <- seq(ee[1], ee[2])
        cell_seq <- c(cell_seq, tmp_v)
      } else {
        tmp_v <- as.integer(ee)
        cell_seq <- c(cell_seq, tmp_v)
      }
    }
    
    cell_seq <- toString(cell_seq)
    result_seq <- paste0(sprintf("where expID in (%s)", cell_seq))
  }
  
  if (args[8]=="all") {
    chrom <- ""
  } else {
    chrom <- sprintf("where chrom in ('%s')", args[8])
  }
  
  if (chrom != "" && result_seq != "") {
    and <- "and"
    result_seq <- gsub("where", "", result_seq)
    result_seq <- gsub(" expID", "expID", result_seq)
  } else {
    and <- ""
  }
  
  cat("Performing SQL query...")
  
  input_cpgs <-  dbSendQuery(dbconn, sprintf("select expID, chrom, pos, ratio
                                             from cpgCoverage
                                             %s
                                             %s
                                             %s", chrom, and, result_seq)) %>% fetch (n=-1) %>% setDT()
  
  
  input_cpgs <- split(input_cpgs, input_cpgs$expID)
  
  for (i in seq_along(input_cpgs)) {
    fwrite(input_cpgs[[i]][,c(2:4)], paste0(args[1], sprintf("/%s.tsv.gz", names(input_cpgs)[i])), sep = '\t')
  }
  
  cat(sprintf("%s cells were written as separate files to the directory: %s \n", length(input_cpgs), args[1]))
  
}

if (args[6] == "mouse_public") {
  inhouse_cpgs_proc(dbconn)

  prom_length <- as.numeric(args[10])
  suppressWarnings(annotations <- dbSendQuery(dbconn, "select chrom, strand, txStart, txEnd
                         from knownGene") %>% fetch(n=-1) %>% setDT() %>%
    .[, feature_start := ifelse(strand == "+", txStart - prom_length, txEnd+1)] %>%
    .[, feature_end := ifelse(strand == "+", txStart-1, txEnd+1 + prom_length)] %>%
    .[, c("strand", "txStart", "txEnd") := NULL] %>% distinct(.) %>%
    .[, strand := "*"] %>% .[, id := rownames(.)] %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)])

  fwrite(annotations, args[2], sep = '\t')

} else if (args[6] == "scm_mm10") {
  scnmt_cpgs_proc(dbconn)

  prom_length <- as.numeric(args[10])
  annotations <- dbSendQuery(dbconn, "select chrom, strand, txStart, txEnd
                         from knownGene") %>% fetch(n=-1) %>% setDT() %>%
    .[, feature_start := ifelse(strand == "+", txStart - prom_length, txEnd+1)] %>%
    .[, feature_end := ifelse(strand == "+", txStart-1, txEnd+1 + prom_length)] %>%
    .[, c("strand", "txStart", "txEnd") := NULL] %>% distinct(.) %>%
    .[, strand := "*"] %>% .[, id := rownames(.)] %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)]

  fwrite(annotations, args[2], sep = '\t')

}

#create melissa data obj
input_melissa <- create_melissa_data_obj(
  met_dir = "/home/igor/Data/" #args[1],
  anno_file = args[2],
  chrom_size_file = NULL,
  chr_discarded = NULL,
  is_centre = TRUE,
  is_window = FALSE,
  upstream = -3000,
  downstream = 1000,
  cov = 1,
  sd_thresh = -1,
  no_cores = as.integer(args[4])
)

cat(sprintf("Melissa object was created: totally %s regions\n", length(input_melissa$met[[1]])))

#filter data
input_melissa <- filter_by_cpg_coverage(input_melissa, min_cpgcov = 1)
cat(sprintf("Melissa object was filtered by cpg coverage: %s regions left\n", length(input_melissa$met[[1]])))
input_melissa <- filter_by_variability(input_melissa, min_var = 0.01)
cat(sprintf("Melissa object was filtered by variability: %s regions left\n", length(input_melissa$met[[1]])))
input_melissa <- filter_by_coverage_across_cells(input_melissa, min_cell_cov_prcg = 0.01)
cat(sprintf("Melissa object was filtered by coverage across cells: %s regions left\n", length(input_melissa$met[[1]])))

####run melissa####

# Set seed for reproducible results
set.seed(15) 
# Create basis object
basis_obj <- create_rbf_object(M = as.integer(args[10]))

cat(sprintf("Running Melissa with supposed %s clusters, %s cores for parallel running\n", args[3], args[4]))

melissa_obj <- melissa(X = input_melissa$met, K = 2, basis = basis_obj,
                       vb_max_iter = 500, vb_init_nstart = 10, 
                       is_parallel = TRUE, no_cores = 4, is_verbose = TRUE)

#saveRDS(melissa_obj, sprintf("melissa_obj_%scl_rbf%s_%s_%s.rds", as.integer(args[3]), args[11], args[6], nrow(melissa_obj$r_nk)))
saveRDS(melissa_obj, "melissa_obj_geo51.rds")

cat(sprintf("Melissa is done. Clustering result:\n"))

print(melissa_obj$pi_k)

####Run bprmeth for each cell#### 

#create list of inputs for bpr
list_of_input_objs <- list()
for (i in seq_along(input_melissa$met)) {
  cur_cell <- list(met = input_melissa$met[[i]], anno = input_melissa$anno_region)
  list_of_input_objs[[i]] <- cur_cell
}

#delete extra na 
for (i in seq_along(list_of_input_objs)) {
  list_of_input_objs[[i]]$met <- list_of_input_objs[[i]]$met[!is.na(list_of_input_objs[[i]]$met)]
}

cat(sprintf("BPRMeth is running\n"))

#run bprmeth for each cell
output_bpr <- list()
for (i in seq_along(list_of_input_objs)) {
  tmp_cell <- infer_profiles_mle(X = list_of_input_objs[[i]]$met, model = "bernoulli",
                                 basis = basis_obj, 
                                 vb_max_iter = 500,
                                 is_parallel = TRUE,
                                 no_cores = 6, #as.integer(args[4]),
                                 is_verbose = TRUE
  )
  output_bpr[[i]] <- tmp_cell
  print(i)
}
saveRDS(output_bpr, "output_bpr_geo51.rds")

cat(sprintf("BPRMeth is done. Time to use output data together from both packages...\n"))

#replacement of real profiles to consensus data
both_packs_list <- list()

#names(output_bpr) <- names(input_melissa$met)

#extract the table to clustered cells
cells_cl <- as.data.frame(melissa_obj$r_nk)

#transform data from numeric chances to just cluster number
for (i in seq_along(cells_cl)) {
  cells_cl[i] <- ifelse(cells_cl[i] > 0.9, i, NA)
}
cells_cl <- cells_cl %>% tidyr::unite("result", cluster1:last(colnames(cells_cl)), na.rm = TRUE, remove = FALSE)
cells_cl[,c(2:ncol(cells_cl))] <- NULL
cells_cl$result <- as.numeric(cells_cl$result)

cat(sprintf("Replacement of consensus met profiles from Melissa by cell profiles: cell processing...\n"))

#the replacement procedure
for (q in seq_along(output_bpr)) {
  #add_col <- as.numeric(args[11]) + 2
  add_col <- 4 + 2
  
  #take the real met profiles
  tmp_df <- as.data.frame(output_bpr[[q]]$W)
  tmp_df[, add_col] <- "natural"
  output_bpr[[q]]$W <- tmp_df
  
  #to define which cluster matrix to use
  melissa_cl <- as.data.frame(melissa_obj$W[,,cells_cl$result[q]])
  
  #to determine which region is real, which is calculated
  for (i in seq_along(input_melissa$met[[q]])) {
    if (is.na(input_melissa$met[[q]][[i]])) {
      melissa_cl[i, add_col] <- "imputed"
    } else {
      melissa_cl[i, add_col] <- "natural"
    }
  }
  
  #to replace only imputed data by profiles from BPRMeth package
  count = 1
  for (j in c(1:nrow(melissa_cl))) {
    if (melissa_cl[,add_col][j] == "natural") {
      melissa_cl[j,] <- output_bpr[[q]]$W[count,]
      count = count + 1
    }
  }
  
  melissa_cl[,add_col] <- NULL
  
  both_packs_list[[q]] <- melissa_cl
  
  print(q)
}


####an old approach to take weights to psa in one column####

# #reshape dfs to one columns (from wide to long)
# for (i in seq_along(both_packs_list)) {
#   both_packs_list[[i]] <- tidyr::gather(both_packs_list[[i]])
#   both_packs_list[[i]]$key <- NULL
#   both_packs_list[[i]] <- as.data.frame(both_packs_list[[i]])
# }
# 
# 
# #create matrix
# infer_prof_sc <- matrix(nrow = nrow(both_packs_list[[1]]), ncol = nrow(cells_cl)) %>% as.data.frame()
# 
# #fill it up 
# for (i in seq_along(both_packs_list)) {
#   infer_prof_sc[,i] <- both_packs_list[[i]]
# }
# 
# #to create long data matrix for initial clusters
# melissa_mtx <- matrix(nrow = nrow(both_packs_list[[1]]), ncol = length(melissa_obj$pi_k))
# 
# #fill it up
# for (i in seq_along(melissa_obj$pi_k)) {
#   melissa_mtx[,i] <- melissa_obj$W[,,i][,1]
# }
# 
# #concatenate tables
# infer_prof_sc_clusters <- cbind(infer_prof_sc, melissa_mtx)
# 
# #for scNMT data - developmental stage defining#
# 
# #to define developmental stage for each cell
# if (args[5] != "no") {
#   sample_metadata <- fread(args[5]) %>%
#     .[, c("id_met", "pass_metQC", "stage")] %>%
#     .[complete.cases(.),] %>% .[,subset(., pass_metQC == "TRUE")] %>%
#     .[, pass_metQC := NULL] %>% setorder()
#   
#   dev_stage <- sample_metadata$stage
#   dev_stage <- c(dev_stage, colnames(melissa_obj$r_nk))
#   
#   #delete stages for that two NAs cells 
#   dev_stage <- dev_stage[-c(41,308)]
#   group <- as.factor(dev_stage)
# }
# 
# colnames(infer_prof_sc_clusters) <- c(rownames(melissa_obj$r_nk), colnames(melissa_obj$r_nk))
# colnames(infer_prof_sc_clusters) <- gsub(".tsv", "", colnames(infer_prof_sc_clusters))
# 
# print(head(infer_prof_sc_clusters[,c(1:10)], 10))
# print(tail(infer_prof_sc_clusters[,c(1:10)], 10))
# 
# #create groups for mix data
# groups <- data.table(cell = colnames(infer_prof_sc_clusters),
#                      set = NA)
# 
# groups$set[1:76] <- "in-house"
# groups$set[77:152] <- "scNMT"
# groups$set[153:155] <- groups$cell[153:155]
# groups <- groups$set %>% as.factor()
# 
# 
# 
# fwrite(infer_prof_sc_clusters, sprintf("replacement_mtx_%s_cl.tsv", as.numeric(args[3])), sep='\t')
# cat("the matrix with replaced values saved as", sprintf("replacement_mtx_%s_cl.tsv", as.numeric(args[3])), "at the same folder where is a script\n")
# 
# cat("Performing PCA\n")
# 
# #Perform PCA
# res_prcomp <- prcomp(t(infer_prof_sc_clusters), scale = TRUE)
# 
# #Plot PCA
# custom_palette <- c("#00AFBB", "#E7B800", "#FC4E07", "#1307fc", "#fc1307")
# 
# fviz_pca_ind(res_prcomp,
#              pointsize = 0.5,
#              labelsize = 4,
#              geom = c("point"),
#              col.ind = groups, # color by groups
#              palette = c(rep("#000000", as.integer(args[3])), custom_palette),
#              addEllipses = FALSE, # Concentration ellipses
#              ellipse.type = "confidence",
#              legend.title = "Groups",
#              repel = FALSE,
#              title = sprintf("PCA of 77 in-house and 77 scNMT cells (RBF = %s)", as.integer(args[11])) #sprintf("PCA of cells based on methylation profiles of promoters 1kb (predefined K = %s)", as.numeric(args[3]))
# )
# 
# pdf(sprintf("%s.pdf", args[9]), width = 9, height = 6)
# plot(pca_plot)
# dev.off()
# 
# #Plotting the melissa region before and after bprmeth package
# 
# #save initial melissa profiles
# for (i in c(1:nrow(melissa_obj$W[,,2]))) {
#   plot_melissa_profiles(melissa_obj = melissa_obj, region = i, title = sprintf("Methylation profiles for region %s by Melissa", i))
#   ggsave(filename = sprintf("melissa_prof_%s.pdf", i), device = "pdf", width = 9, height = 6)                    
#   print(i)
# }
# 
# #for cell 1 plot profiles after replacement by real bprmeth profiles
# melissa_bprmeth_obj_1 <- melissa_obj 
# melissa_bprmeth_obj_1$W[,,3] <- as.matrix(both_packs_list[[1]])
# 
# for (i in c(1:nrow(melissa_obj$W[,,2]))) {
#   plot_melissa_profiles(melissa_obj = melissa_bprmeth_obj_1,
#                         region = i, title = sprintf("Methylation profiles for region %s by Melissa and BPR", i),
#                         )
#   ggsave(filename = sprintf("melissa_bpr_prof_%s.pdf", i), device = "pdf", width = 9, height = 6)                    
#   print(i)
# }
# 
# #plot bernoulli profiles and individual cpg obs
# 
# basis_prof <- create_rbf_object(M = 4)
# basis_mean <- create_rbf_object(M = 0)
# bern_prof <- infer_profiles_vb(X = bernoulli_data, model = "bernoulli",
#                                basis = basis_prof, is_parallel = FALSE)
# bern_mean <- infer_profiles_vb(X = bernoulli_data, model = "bernoulli",
#                                basis = basis_mean, is_parallel = FALSE)
# p <- plot_infer_profiles(region = 3, obj_prof = tmp_cell, 
#                           obs = bernoulli_data, 
#                          title = "Bernoulli profile by BPRMeth")
# print(p)



####a new approach to take mean met value for one region####

#for finding a 1 mean value for each cluster curve
replaced_cells <- list()
cells_table <- matrix(ncol = length(both_packs_list), nrow = nrow(both_packs_list[[1]]))
for (i in seq_along(both_packs_list)) {
  #create empty matrix each region * N clusters per cell
  regions_means <- c() #matrix(ncol = length(melissa_obj$W_Sigma), nrow = length(melissa_cl$V1))
  
  #melissa object for the cell
  melissa_tmp_cell <- melissa_obj
  #replace the weight matrix of which cell it belongs
  melissa_tmp_cell$W[,,cells_cl$result[i]] <- as.matrix(both_packs_list[[i]])
  #put the modified melissa object to the list
  replaced_cells[[i]] <- melissa_tmp_cell
  
  #in this cell to find means of each cluster for each regions 
  for (r in seq_along(replaced_cells[[i]]$W_Sigma[[1]])) {
    cell_region_means <- infer_melissa_mean(melissa_obj = replaced_cells[[i]], region = r, cluster = cells_cl$result[i])
    regions_means <- c(regions_means, cell_region_means)
  }
  cells_table[,i] <- regions_means
  print(i)
}

#for cluster centers
cluster_means <- matrix(ncol = length(melissa_obj$pi_k), nrow = nrow(cells_table))

for (k in seq_along(melissa_obj$pi_k)) {
  for (r in c(1:nrow(cells_table))) {
    mean_1region <- infer_melissa_mean(melissa = melissa_obj, region = r, cluster = k)
    cluster_means[r, k] <- mean_1region 
  }
}

final_table <- cbind(cells_table, cluster_means)
colnames(final_table) <- c(gsub("(.tsv|.txt)", "", names(input_melissa$met), perl = TRUE), seq(length(melissa_obj$pi_k)))

fwrite(final_table, paste0(args[12], sprintf("melissa_bpr_%scl_%srbf_met_values.tsv", args[3], args[11])), sep='\t')



#PCA
#Perform PCA
res_prcomp <- prcomp(t(final_table), scale = TRUE)


#only for smallwood data
groups <- data.table(cell = c(colnames(final_table), "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7"),
                    set = NA)

#groups$set[1:76] <- "in-house"
#groups$set[77:88] <- "MII (ovulated) oocytes"
#groups$set[89:100] <- "ESC cultured in 2i/LIF"
#groups$set[101:119] <- "ESC cultured in serum/LIF"
#groups$set[120:123] <- seq(length(melissa_obj$pi_k))
groups$set[1:629] <- "S1"
groups$set[630:1209] <- "S2"
groups$set[1210:1696] <- "S3"
groups$set[1697:1703] <- seq(1:7)

groups <- groups$set %>% as.factor()

#Plot PCA
custom_palette <- c("#00AFBB", "#E7B800", "#FC4E07", "#1307fc", "#fc1307", "#7d8050", "#1600bb")

fviz_pca_ind(res_prcomp,
             pointsize = 0.5,
             labelsize = 4,
             geom = c("point"),
             col.ind = groups, # color by groups
             palette = c(rep("#000000", 7), custom_palette),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = FALSE,
             title = sprintf("PCA of cells based on methylation profiles of promoters 1kb (predefined K = %s)", 7)
) 

pdf(sprintf("%s.pdf", args[9]), width = 9, height = 6)
plot(pca_plot)
dev.off()
