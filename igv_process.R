#for the first region

setwd("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/first_reg")
seg_files = list.files(".")
file_id <- gsub(".tsv", "", seg_files)
for (f in seq_along(seg_files)) {
  tmp <- fread(seg_files[f])
  tmp$ID <- "chr1:9647976-9648975"
  tmp <- tmp[,c(5,1:4)]
  colnames(tmp) <- c("ID","chrom", "loc.start", "loc.end", "seg.mean")
  tmp$loc.end = tmp$loc.end + 1
  fwrite(tmp, sprintf("%s.seg", file_id[f]), sep='\t')
  print(f)
}

#add 1 base to loc.end
seg_files = list.files("/Users/igorbulanov/Documents/Data/Epiclomal_results/15_12_300KC_allchr_clustering/5_0.99_10000/visual/cl_for_ucsc/first_reg")
seg_files <- gsub(".tsv", "", seg_files)

for (f in seq_along(seg_files)) {
  tmp <- fread(seg_files[f])
  fwrite(tmp, seg_files[f], sep='\t') 
}
