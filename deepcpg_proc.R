library(data.talbe)
library(tibble)
library(dplyr)

cov_files = list.files("/Users/igorbulanov/Documents/Data/KC_coverage/DeepCpG_cov_300KC_all_chr")

for (f in seq_along(cov_files)) {
  tmp <- fread(cov_files[f]) %>% .[, c("chr", "CpG_start", "meth_frac")]
  tmp$chr <- gsub("chr", "", tmp$chr)
  fwrite(tmp, sprintf("%s", cov_files[f]), sep='\t', col.names = FALSE)
  print(f)
}

