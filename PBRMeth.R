library(BPRMeth)

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
  cluster_cell <- cluster_profiles_vb(X = met_region$met, K = 1, model = "bernoulli", 
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
res.pca = PCA(w_main, scale.unit=TRUE, ncp=5, graph=T)

#UMAP
library(umap)
result_umap <- umap(w_data, config = umap.defaults, method = "umap-learn")
plot_umap(result_umap$layout, w_label)

plot_umap <- function(x, labels,
         main="UMAP visualization of the clustered weights (promoter 1kb, rbf = 4)",
         colors=c("#ff7f00", "#17becf"),
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

#create matrix with 395 cols (79x5)

infer_prof_f

