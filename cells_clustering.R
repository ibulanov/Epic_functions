library(Melissa)
library(BPRMeth)
library(factoextra)
library(RMySQL)
library(dplyr)
library(data.table)

#args
args = commandArgs(trailingOnly=TRUE)
#1 - cell files directory (for storing from db and then this directory will be used for Melissa data object creation)
#2 - annotation file
#3 - Number of clusters for Melissa running
#4 - Number of cores for Melissa and BPRMeth
#5 - scNMT sample metadata file - not use now
#6 - Database name
#7 - desired cell IDs (ranges, numbers, e.g. "100310:100384, 100386" or "all")
#8 - To choose desired data (e.g. either only 1 chr or all, e.g. 'chr1' or 'all') 
#9 - Name of output PCA plot


dbconn <- dbConnect(MySQL(), user='i747v', password='igor',
                    dbname=args[6], host='zuse1db')


inhouse_cpgs_proc <- function(dbconn) {
  
  if (args[7] = "all") {
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
  
  #to get the input data (cpgBases)
  input_cpgs <-  dbSendQuery(dbconn, sprintf("select C.expID, B.chrom, B.pos, B.strand, C.ratio
                        from cpgCoverage C
                        inner join cpgBases B on ((C.pos = B.pos) AND (C.chrom = B.chrom))
                        %s
                        %s
                        %s", chrom,and,result_seq)) %>% fetch(n=-1) 
  
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
    fwrite(input_cpgs[[i]][,c(2:4)], paste0(args[1], names(input_cpgs)[i], ".tsv"), sep = '\t')
  }
  
  cat(sprintf("%s cells were written as separate files to the directory: %s \n", length(input_cpgs), args[1]))
  
}

scnmt_cpgs_proc <- function(dbconn, result_seq=result_seq) {
  
  if (args[7] = "all") {
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
  
  input_cpgs <-  dbSendQuery(dbconn, sprintf("select expID, chrom, pos, ratio
                                             from cpgCoverage
                                             %s
                                             %s
                                             %s", chrom, and, result_seq)) %>% fetch (n=-1)
  
}
#
input_cpgs <- cpgs_proc()

#create melissa data obj
input_melissa <- create_melissa_data_obj(
  met_dir = args[1],
  anno_file = args[2],
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

cat(sprintf("Melissa object was created: totally %s regions\n", length(input_melissa$met[[1]])))

#filter data
input_melissa <- filter_by_cpg_coverage(input_melissa, min_cpgcov = 3)
cat(sprintf("Melissa object was filtered by cpg coverage: %s regions left\n", length(input_melissa$met[[1]])))
input_melissa <- filter_by_variability(input_melissa, min_var = 0.05)
cat(sprintf("Melissa object was filtered by variability: %s regions left\n", length(input_melissa$met[[1]])))
input_melissa <- filter_by_coverage_across_cells(input_melissa, min_cell_cov_prcg = 0.05)
cat(sprintf("Melissa object was filtered by coverage across cells: %s regions left\n", length(input_melissa$met[[1]])))

####run melissa####

# Set seed for reproducible results
set.seed(15) 
# Create basis object
basis_obj <- create_rbf_object(M = 4)

cat(sprintf("Running Melissa with supposed %s clusters, %s cores for parallel running\n", args[3], args[4]))

melissa_obj <- melissa(X = input_melissa$met, K = as.integer(args[3]), basis = basis_obj,
                       vb_max_iter = 500, vb_init_nstart = 10, 
                       is_parallel = TRUE, no_cores = as.integer(args[4]), is_verbose = TRUE)

cat(sprintf("Melissa is done. Clustering result:\n"))
print(melissa_obj$pi_k)

####Run bprmeth for each cell#### 

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

cat(sprintf("BPRMeth is running\n"))

#run bprmeth for each cell
output_bpr <- list()
for (i in seq_along(list_of_input_objs)) {
  tmp_cell <- infer_profiles_mle(X = list_of_input_objs[[i]]$met, model = "bernoulli",
                                 basis = basis_obj, 
                                 vb_max_iter = 500,
                                 is_parallel = TRUE,
                                 no_cores = as.integer(args[4]),
                                 is_verbose = TRUE
  )
  output_bpr[[i]] <- tmp_cell
  print(i)
}

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
  #take the real met profiles
  tmp_df <- as.data.frame(output_bpr[[q]]$W)
  tmp_df$w6 <- "natural"
  output_bpr[[q]]$W <- tmp_df
  
  #to define which cluster matrix to use
  melissa_cl <- as.data.frame(melissa_obj$W[,,cells_cl$result[q]])
  
  #to determine which region is real, which is calculated
  for (i in seq_along(input_melissa$met[[q]])) {
    if (is.na(input_melissa$met[[q]][[i]])) {
      melissa_cl[i, 6] <- "imputed"
    } else {
      melissa_cl[i, 6] <- "natural"
    }
  }
  
  #to replace only imputed data by profiles from BPRMeth package
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

#reshape dfs to one columns (from wide to long)
for (i in seq_along(both_packs_list)) {
  both_packs_list[[i]] <- tidyr::gather(both_packs_list[[i]])
  both_packs_list[[i]]$key <- NULL
  both_packs_list[[i]] <- as.data.frame(both_packs_list[[i]])
}

#create matrix
infer_prof_sc <- matrix(nrow = nrow(both_packs_list[[1]]), ncol = nrow(cells_cl)) %>% as.data.frame()

#fill it up 
for (i in seq_along(both_packs_list)) {
  infer_prof_sc[,i] <- both_packs_list[[i]]
}

#to create long data matrix for initial clusters
melissa_mtx <- matrix(nrow = nrow(both_packs_list[[1]]), ncol = length(melissa_obj$pi_k))

#fill it up
for (i in seq_along(melissa_obj$pi_k)) {
  melissa_mtx[,i] <- melissa_obj$W[,,i]
}

#concatenate tables
infer_prof_sc_clusters <- cbind(infer_prof_sc, melissa_mtx)

#for scNMT data - developmental stage defining#

#to define developmental stage for each cell
sample_metadata <- fread("/home/igor/Data/scnmt_gastrulation/sample_metadata.txt") %>%
  .[, c("id_met", "pass_metQC", "stage")] %>%
  .[complete.cases(.),] %>% .[,subset(., pass_metQC == "TRUE")] %>%
  .[, pass_metQC := NULL] %>% setorder()

dev_stage <- sample_metadata$stage
dev_stage <- c(dev_stage, colnames(melissa_obj$r_nk))
group <- as.factor(dev_stage)
colnames(infer_prof_sc_clusters) <- c(rownames(melissa_obj$r_nk), colnames(melissa_obj$r_nk))


cat("Performing PCA\n")
#Perform PCA
res_prcomp <- prcomp(t(infer_prof_sc_clusters), scale = TRUE)

#Plot PCA
pca_plot <- fviz_pca_ind(res_prcomp,
             pointsize = 1.5,
             labelsize = 4,
             geom = c("point"),
             #col.ind = group, # color by groups
             #palette = c("#000000", "#000000", "#000000", "#00AFBB", "#E7B800", "#FC4E07", "#1307fc", "#fc1307"),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = FALSE,
             title = "PCA of cells based on methylation profiles of promoters 1kb"
)

pdf(sprintf("%s.pdf", args[9]), width = 9, height = 6)
plot(pca_plot)
dev.off()

