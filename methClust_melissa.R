library(Melissa)

input.data <- fread("/home/igor/input.data/mouse_scimet/cell100386.txt")
example_input <- melissa_encode_dt

head(input.data$met[[2]][[50]])

cat("Number of cells: ", length(input.data$met))

#creating basis object

library(BPRMeth)

# Create RBF basis object with 4 RBFs
basis_obj <- create_rbf_object(M = 4)

#clustering and imputing scBS-seq input.data
set.seed(15)
# Partition to training and test set
input.data <- partition_input.dataset(data = input.data, input.data_train_prcg = 0.2,
                            region_train_prcg = 1, cpg_train_prcg = 0.4, 
                            is_synth = TRUE)

# Run Melissa with K = 4 clusters
melissa_obj <- melissa(X = input.data2, K = 2, basis = basis_obj,
                       vb_max_iter = 20, vb_init_nstart = 1, is_parallel = FALSE)


melissa_obj$pi_k
head(melissa_obj$r_nk)
melissa_obj$W[10, , 3]

#plotiing the methylation profiles

plot_melissa_profiles(melissa_obj = melissa_obj, region = 22, 
                      title = "Methylation profiles for region 22")

plot_melissa_profiles(melissa_obj = melissa_obj, region = 77, 
                      title = "Methylation profiles for region 77")


#evaluating the clustering perfomance

# Run clustering performance
melissa_obj <- eval_cluster_performance(melissa_obj, input.data$opts$C_true)

# ARI metric
cat("ARI: ", melissa_obj$clustering$ari)
# Clustering assignment error metric
cat("Clustering assignment error: ", melissa_obj$clustering$error)


#evaluation the imputation perfomance
imputation_obj <- impute_met_state(obj = melissa_obj, test = input.data$met_test)

melissa_obj <- eval_imputation_performance(obj = melissa_obj, 
                                           imputation_obj = imputation_obj)

# AUC 
cat("AUC: ", melissa_obj$imputation$auc)

# F-measure
cat("F-measure: ", melissa_obj$imputation$f_measure)



