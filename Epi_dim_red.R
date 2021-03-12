#24/01/21 Dimensional reduction plots of Epiclomal results

library(tibble)
library(dplyr)
library(data.table)
library(M3C)
library(factoextra)
suppressMessages(library(Rtsne))
suppressMessages(library(NMF))
suppressMessages(library(pheatmap))

#1 - 300KC prom1kb
#UMAP


prepare_epiclomal <- function(best_dir, input_dir, out_dir, data_id, true) {
  
  cluster_map <- fread(paste0(best_dir, "cluster_MAP.tsv.gz"), header = TRUE) %>%
    .[,c(1,2)] 
  colnames(cluster_map) <- c("cell_id", "cluster")
  cluster_map$cluster <- cluster_map$cluster + 1
  
  setwd(input_dir)
  input_epiclomal <- fread(paste0(input_dir, sprintf("input_Epiclomal_%s.tsv.gz", data_id))) %>%
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
  
  
  #result_comb <- cbind(result_comb, region_map)
  
  #tSNE
  max_perplexity = floor((ncol(result_comb)-1)/3)
  num_tsne_iter <- 5000
  
  for (perp in unique(floor(seq(8,max_perplexity,length.out=5)))) {
    clusterer <- paste0("tSNE_perp",perp)
    rfile <- paste0(out_dir,"results_", clusterer ,".txt")
    if (!file.exists(rfile)) {    
      print(paste0("Running only tSNE with perplexity ", perp))  
      htsne <- Rtsne(t(result_comb), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)
      png(file=paste0(out_dir,data_id,"_tSNE_perp",perp,".png"))
      plot(htsne$Y, col= as.factor(true$epigenotype_id), pch = as.numeric(factor(cluster_map$cluster)))
      #plot(htsne$Y, col=as.factor(c(cluster_map$cluster, rep(max(cluster_map$cluster)+1, max(cluster_map$cluster)))))
      #title(main=paste0(data_id, ", tSNE (perp", perp, "), V=", format(round(vmeasure, 2), nsmall = 2)))
      title(main=paste0(data_id, ", tSNE (perp", perp, ")"))
      dev.off()
    } else {
      print(paste0("Results already exist in file ", rfile))
    }        
  }

  #NMF
  min_nmf_rank <- 2
  max_nmf_rank <- 1       #15
  num_nmf_runs <- 10
  
  #replace 0s to 0.0001 for NMF
  result_comb_nmf <- result_comb
  result_comb_nmf[result_comb_nmf==0] <- 0.0001
  
  # Run NMF on a range of ranks from 2 to 15
  for (nmf_rank in c(2:7)) {
    nmf_file <- paste0(out_dir,"nmf_rank", nmf_rank, ".rds")
    if (!file.exists(nmf_file)) {    
      print (paste0("Running NMF with rank ", nmf_rank))
      res <- nmf(result_comb_nmf, nmf_rank, nrun=num_nmf_runs, .opt='vP')
      saveRDS(res, nmf_file)
    } else {
      res <- readRDS(file=nmf_file)
    }   
    
    # get all the features and plot a heatmap
    s <- extractFeatures(res)
    # do this a few times with less and more features
    for (max_num_features in c(50,100,500,1000)) {
      region_set <- c()
      for (i in 1:nmf_rank){
        if (!is.na(s[[i]])) {
          region_set <- c(region_set, s[[i]][1:min(max_num_features,length(s[[i]]))])
        }        
      }
      region_set <- unique(sort(region_set))
      selected_input <- result_comb_nmf[region_set,]
      write.table(rownames(selected_input),file=paste0(out_dir,"/NMF_regions_rank",nmf_rank,"_maxfeatures", max_num_features, ".tsv"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)     
      #pheatmap(t(selected_input),cluster_rows = TRUE,cluster_cols=FALSE,fontsize = 8, 
      #         fontsize_row=2,fontsize_col=6,show_colnames = TRUE,
      #         filename = paste0(out_dir,"/heatmap_NMF_few_regions_mean_meth_rank", nmf_rank, "_maxfeatures", max_num_features, ".png"))               
    }    
    
    # Follow with tSNE, for up to 5 perplexities
    for (perp in unique(floor(seq(8,max_perplexity,length.out=5)))) {
      clusterer <- paste0("NMF_rank",nmf_rank,"_tSNE_perp",perp)
      rfile <- paste0(out_dir,"results_", clusterer ,".txt")
      if (!file.exists(rfile)) {            
        print(paste0("    tSNE with perplexity ", perp))
        htsne <- Rtsne(t(coef(res)), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)
        #cl <- hdbscan(htsne$Y, minPts = min_dbscan_points)
        #vmeasure <- calc_vmeasure(cl, clusterer)
        png(file=paste0(out_dir,data_id,"_NMF_rank",nmf_rank,"_tSNE_perp",perp,".png"))
        plot(htsne$Y, col= as.factor(true$epigenotype_id), pch = as.numeric(factor(cluster_map$cluster)))
        #plot(htsne$Y, col=as.factor(c(cluster_map$cluster, rep(max(cluster_map$cluster)+1, max(cluster_map$cluster)))), pch=as.numeric(factor(true$epigenotype_id)))
        title(main=paste0(data_id, ", NMF (rank ",nmf_rank,")+tSNE (perp ",perp, ")"))
        dev.off()
      } else {
        print(paste0("Results already exist in file ", rfile))
      }            
    }
  }
}

true <- fread("/home/igor/Data/Epiclomal_results/true_files/SW_32_mm10_new.txt.gz")
best_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_SW_mm10_CGI_clust/0_0.95_10000/epi_region/462/"
input_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_sw_mm10_CGI_preproc/epiclomal_input/0_0.95_10000/"
out_dir <- "/home/igor/Data/Epiclomal_results/18_02_32_SW_mm10_CGI_clust/0_0.95_10000/custom_visual/"
data_id <- "32_Sw"

prepare_epiclomal(best_dir,input_dir,out_dir,data_id,true)


#true for epid top200
epid_200 = list.files("/home/igor/Data/mysql_data/epi_dermis_hg19/top200/")
epid_200 = gsub("_cov.tsv.gz", "", epid_200)

true_epid200 = data.table(cell_id = epid_200, epigenotype_id = c(rep(1, 100), rep(2,100)))


##for the datasets with numeric colnames
#result_comb <- as.data.table(result_comb, keep.rownames = TRUE)
#result_comb <- column_to_rownames(result_comb, var = "rn")
#new_colorder <- sort(as.numeric(colnames(result_comb))) %>% as.character()
#setcolorder(result_comb, new_colorder)
#colnames(result_comb) <- paste0("cell", colnames(result_comb))

setwd(out_dir)
#M3C umap
png(f = sprintf("UMAP_%s.pdf", data_id), width = 9, height = 6)
M3C::umap(result_comb,
          labels=as.factor(true$epigenotype_id),#c(cluster_map$cluster, rep(max(cluster_map$cluster)+1, max(cluster_map$cluster))),
          dotsize = 2,
          axistextsize = 10,
          legendtextsize = 15
)
dev.off()

#PCA
result_comb_pca <- result_comb
result_comb_pca[result_comb_pca==0] <- 0.001

res_prcomp <- prcomp(t(result_comb_pca), scale = FALSE)

pdf(sprintf("PCA_%s.pdf", data_id), width = 9, height = 6)

plot(res_prcomp$x, col=as.numeric(as.factor(true$epigenotype_id)))

fviz_pca_ind(res_prcomp,
             pointsize = 3,
             labelsize = 4,
             geom = c("point"),
             col.ind = as.factor(true$epigenotype_id),
             #col.ind = as.factor(c(cluster_map$cluster, rep("cluster centers", max(cluster_map$cluster)))), # color by groups
             #palette = c(rep("#000000", 7), custom_palette),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = FALSE,
             title = "PCA of 963 in-house cells"
             #title = sprintf("PCA of %s cells", data_id)
) 
dev.off()
