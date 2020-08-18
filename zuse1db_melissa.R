library(RMySQL)
library(tidyverse)
library(data.table)
library(Melissa)
library(BPRMeth)

dt_obj <- melissa_encode_dt

mouse_public <- dbConnect(MySQL(), user='i747v', password='igor', 
                     dbname='mouse_public', host='zuse1db')

zuse1db_preproc <- function() {
  
  zuse1db <- dbSendQuery(mouse_public, "select C.expID, C.chrom, C.pos, B.strand, C.coverage, C.mC, C.ratio 
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        where C.expID between 100310 and 100386") %>% fetch(n=-1)
  
  #or zuse1db <- fread("/home/igor/Data/common_matrices/input_cpgs_23kk.tsv", sep='\t')
  
  input_cpgs <- setDT(zuse1db) %>% .[, compl_pos := ifelse(strand == 1, pos+1, pos-1)]
  
  #create the complementary df
  cpgs_compl <- data.table(input_cpgs$expID, input_cpgs$chrom, input_cpgs$compl_pos, input_cpgs$ratio)
  
  #remove extra columns
  input_cpgs <- input_cpgs[,c("strand","compl_pos"):= NULL] 
  
  colnames(cpgs_compl) <- colnames(input_cpgs)
  
  #concatenate the two dataframes and order the rows
  input_cpgs <- rbind(input_cpgs, cpgs_compl) %>% arrange(expID) %>% distinct() 
  
  input_cpgs$pos_d <- input_cpgs$pos
  
  #delete cpgbases complementary data table
  cpgs_compl <- NULL
  
  cat(sprintf("Cpg bases positions for cells 100310:100386 were generated\n"))
  
  return(input_cpgs)
  
}

process_zuse1db_full <- function(zuse1db) {
  
  zuse1db <- dbSendQuery(mouse_public, "select C.expID, C.chrom, C.pos, B.strand, C.coverage, C.mC, C.ratio 
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        where C.expID between 100310 and 100386") %>% fetch(n=-1)
  
  zuse1db <-  zuse1db[, compl_pos := ifelse(strand == 1, pos+1, pos-1)] %>% .[, unmet_reads := coverage - mC] %>%
    .[, c("expID", "chrom", "pos", "compl_pos", "ratio", "mC", "unmet_reads")]
  
  #create the complementary df
  cpgs_compl <- data.table(zuse1db$expID, zuse1db$chrom, zuse1db$compl_pos, zuse1db$ratio, zuse1db$mC, zuse1db$unmet_reads)
  
  #remove extra columns
  zuse1db <- zuse1db[, compl_pos := NULL] 
  
  colnames(cpgs_compl) <- colnames(zuse1db)
  
  #concatenate the two dataframes and order the rows
  zuse1db <- rbind(zuse1db, cpgs_compl) %>% arrange(expID) %>% distinct() 
  
  zuse1db$pos_d <- zuse1db$pos
  
  #delete cpgbases complementary data table
  cpgs_compl <- NULL
  
  cat(sprintf("Cpg bases positions for cells 100310:100386 were generated\n"))
  
  return(zuse1db)
  
  
  
}

windows_proc <- function(input_cpgs, winsize) {
  
  input_cpgs[, start := pos - pos %% winsize] %>% .[, finish := start + (winsize - 1)] %>%
    .[, c("pos") := NULL] %>% .[, feature_coord := paste(chrom, ":",start, ":", finish)] %>%
    .[, c("chrom", "start", "finish") := NULL] %>% .[, feature_coord := gsub(" ", "", feature_coord)] %>% 
    distinct() %>% setkey(expID)
  
  unique_features_number <- nrow(input_cpgs)

}

prom_db_proc <- function(prom_length) {
  #take data for promoters
  output_db <- dbSendQuery(mouse_public, "select chrom, strand, txStart, txEnd
                         from knownGene") %>% fetch(n=-1) %>% setDT() %>% 
    .[, feature_start := ifelse(strand == "+", txStart - prom_length, txEnd)] %>%
    .[, feature_end := ifelse(strand == "+", txStart -1, txEnd + prom_length - 1)] %>%
    .[, c("strand", "txStart", "txEnd") := NULL] %>% distinct()
  
  unique_features_number <- nrow(output_db)
  
  #set keys for data.table class
  setkey(input_cpgs, chrom, pos, pos_d)
  setkey(output_db, chrom, feature_start, feature_end)
  
  #find cpg sites positions in promoters, create concatenated column and delete extra columns
  cells_features <- foverlaps(input_cpgs, output_db, nomatch = NULL) %>%
    .[,c("pos_d"):= NULL] %>% setkey(expID) %>% distinct() %>%
    .[, feature_coord := paste(chrom, ":", feature_start, ":", feature_end)] %>%
    .[, feature_coord := gsub(" ", "", feature_coord)] #%>%
    .[,c("chrom","feature_start", "feature_end"):= NULL]
 
  return(cells_features)
   
}
  
cells_features <- prom_db_proc(prom_length)


#### try to mutate the input data####
#to add binary ratio
cells_features <- cells_features %>% 
  mutate(bin_ratio = ifelse(ratio >= 0.5 , 1, 0))

#to normalize pos values to [-1;1]

cells_features$norm_pos <- 2*((cells_features$pos - cells_features$feature_start) /
                                (cells_features$feature_end - cells_features$feature_start))-1


input_melissa <- cells_features[, c(1,4,7,8,10)] %>% .[,c(2,3,1,5,4)] %>% 
  .[, rn_pos := round(norm_pos, digits = 3) * 1000] %>%
  .[, rn := paste(chrom, ":", rn_pos)] %>%
  .[, rn := gsub(" ", "", rn)] %>%
  .[, c("chrom", "rn_pos") := NULL]

#split the table as test input data for melissa
input_melissa_spl <- split(input_melissa, input_melissa$expID)

#test
input_melissa_test <- input_melissa_spl[c(1:5)]


#split every genomic region for each cell
for (i in seq_along(input_melissa_spl)) {
  input_melissa_spl[[i]] <- split(input_melissa_spl[[i]], input_melissa_spl[[i]]$feature_coord)
  print(i)
}

#move rn column to rownames and delete extra columns according to test data
for (i in seq_along(input_melissa_spl)) {
  for (j in seq_along(input_melissa_spl[[i]]))
    input_melissa_spl[[i]][[j]] <- input_melissa_spl[[i]][[j]] %>% remove_rownames() %>% 
      .[, rn := make.unique(rn)] %>% column_to_rownames(var = "rn") %>% select(norm_pos, bin_ratio) %>% as.matrix()
  print(c(i,j))
}

#to extract promoters of each cells to create reduced list of promoters
cells_proms <- input$met
for (i in seq_along(input$met)) {
  cells_proms[[i]] <- names(input$met[[i]]) 
}

# using default function (create_melissa_data_obj) to create an input data
#for met_dir
for (i in seq_along(zuse1db_3_1)) {
  cur_dir <- sprintf("/home/igor/Data/zuse1db/melissa/input_cells_3_1/cell_%s.tsv", names(zuse1db_3_1[i]))
  zuse1db_3_1[[i]]$expID <- NULL
  fwrite(zuse1db_3_1[[i]], cur_dir, sep = '\t')
}

#for anno_file
output_db[, strand := "*"] %>% .[, id:= rownames(output_db)] %>% 
  .[, name:= paste(chrom, ":", feature_start, ":", feature_end)] %>%
  .[, name:= gsub(" ", "", name)]

fwrite(output_db, "/home/igor/Data/zuse1db/melissa/anno_files/anno_prom_1k_39k.tsv", sep='\t')


#### melissa pipeline ####
binarise_files(indir = "/home/igor/Data/zuse1db/melissa/input_cells_unmet/",
               outdir = "/home/igor/Data/zuse1db/melissa/input_cells_unmet/binarised", format = 2, no_cores = NULL)


input_melissa <- create_melissa_data_obj(
  met_dir = "/home/igor/Data/zuse1db/melissa/input_cells_3_1/",
  anno_file = "/home/igor/Data/zuse1db/melissa/anno_files/anno_win_1k.tsv",
  chrom_size_file = NULL,
  chr_discarded = NULL,
  is_centre = TRUE,
  is_window = FALSE,
  upstream = -5000,
  downstream = 5000,
  cov = 5,
  sd_thresh = -1,
  no_cores = 5
) %>% filter_by_cpg_coverage(min_cpgcov = 3) %>% filter_by_variability(min_var = 0.05) %>%
  filter_by_coverage_across_cells(min_cell_cov_prcg = 0.05)
  

library(BPRMeth)
# Create RBF basis object with N RBFs
basis_obj <- create_rbf_object(M = 4)
basis_mean <- create_rbf_object(M = 0)

#partitioning the dataset
set.seed(15)
# Partition to training and test set
input_melissa_prom_1k_37_part <- partition_dataset(dt_obj = input_melissa_prom_1k_37, data_train_prcg = 0.2,
                                region_train_prcg = 1, cpg_train_prcg = 0.4, 
                                is_synth = FALSE)

#RUN MELISSA
set.seed(15)
# Run Melissa with K = 4 clusters
melissa_out_prom1k_3_1 <- melissa(X = input_melissa_prom_1k_3_1$met, K = 2, basis = basis_obj,
                       vb_max_iter = 500, vb_init_nstart = 10, 
                       is_parallel = TRUE, no_cores = 5)

melissa_obj$pi_k

#evaluate clustering performing
melissa_obj <- eval_cluster_performance(melissa_obj, melissa_obj$r_nk)

#plot
plot_melissa_profiles(melissa_obj = melissa_out_prom1k_3_1, region = 5000, 
                      title = "Methylation profiles for region 10")

#run BPRMeth package for each separate cell but for all region

list_of_input_objs <- list()

for (i in seq_along(input_melissa$met)) {
  cur_cell <- list(met = input_melissa$met[[i]],
                   anno = input_melissa$anno_region)
  list_of_input_objs[[i]] <- cur_cell
}


cluster_profiles_cells_f_13k_base0 <- list()
for (i in seq_along(list_of_input_objs)) {
  cur_cluster_13k <- cluster_profiles_vb(X = cell_13k$met, K = 2, model = "bernoulli",
                                alpha_0 = .5, beta_0 = .1,
                                basis = basis_mean, 
                                vb_max_iter = 500, 
                                is_verbose = TRUE
                                )
  cluster_profiles_cells_f_13k_base0[[i]] <- cur_cluster_13k
  print(i)
}

#delete extra na 
for (i in seq_along(list_of_input_objs)) {
  list_of_input_objs[[i]]$met <- list_of_input_objs[[i]]$met[!is.na(list_of_input_objs[[i]]$met)]
}

cl_13k_13k <- cluster_profiles_vb(X = cell_13k$met, K = 2, model = "bernoulli",
                                   alpha_0 = .5, beta_0 = .1,
                                   basis = basis_obj, 
                                   vb_max_iter = 500, 
                                   is_verbose = TRUE
)

#assign 13k annotations region to flexible numbers of met values
for (i in seq_along(list_of_input_objs)) {
  list_of_input_objs[[i]]$anno <- cell_13k$anno
}

#extract values
w_mtx_cells <- data.frame(cell_id = names(input_melissa_obj$met), cluster_1 = NA, cluster_2 = NA)
w_mtx_cells$cluster_1 <- cluster_profiles_cells_f_13k[[i]]$W[1,][[1]]

for (i in seq_along(cluster_profiles_cells_f_13k)) {
  w_mtx_cells$cluster_1[i] <- cluster_profiles_cells_f_13k[[i]]$W[1,][[1]]
  w_mtx_cells$cluster_2[i] <- cluster_profiles_cells_f_13k[[i]]$W[1,][[2]]
}

#try to run infer profiles_vb

fit_profiles <- infer_profiles_vb(X = cell_13k$met, model = "binomial",
                                  basis = basis_obj, is_parallel = TRUE)

#for 77 cells
infer_profiles_bases_0 <- list()
for (i in seq_along(list_of_input_objs)) {
  cur_profile <- infer_profiles_vb(X = list_of_input_objs[[i]]$met, model = "bernoulli",
                                         alpha_0 = .5, beta_0 = .1,
                                         basis = basis_mean, 
                                         vb_max_iter = 500,
                                         is_parallel = TRUE,
                                         no_cores = 5,
                                         is_verbose = TRUE
  )
  infer_profiles_bases_0[[i]] <- cur_profile
  print(i)
}

#detect min length

ll <- c()
for (i in seq_along(infer_profiles_bases_0)) {
  cur_length <- length(infer_profiles_bases_0[[i]]$W)
  ll <- append(cur_length, ll) 
}

#matrix of inferred profiles (cells and genomic regions)
cell_region_weights <- res_matrix <- matrix(nrow = max(ll), ncol = 77)

for (i in seq_along(infer_profiles_bases_0)) {
  cell_region_weights[,i][c(1:length(infer_profiles_bases_0[[i]]$W))] <- infer_profiles_bases_0[[i]]$W
}

cell_region_weights_cut <- cell_region_weights[c(1:min(ll)), c(1:77)]
colnames(cell_region_weights_cut) <- seq(100310, 100386)

#PCA by factoMineR
library(FactoMineR)
res.pca = PCA(t(cell_region_weights_cut), scale.unit=FALSE, ncp=5, graph=T)


