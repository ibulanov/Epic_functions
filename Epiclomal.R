library(data.table)
library(REpiclomal)

epi_load_test <- load_data(outdir = ".",
                      input_CpG_data_file = "/home/igor/Data/epiclomal/zuse1db_1k.tsv",
                      input_regions_file = "/home/igor/Data/zuse1db/epiclomal/anno_2k_chr1.tsv",
                      use_cache = FALSE
                      )

#generate test data
zuse1db <- fread( "/home/igor/Data/zuse1db/melissa/zuse1db_perm_3_1.tsv", sep='\t')

zuse1db <- split(zuse1db, zuse1db$expID)
zuse1db[c(6:72)] <- NULL

annos <- fread("/home/igor/Data/zuse1db/melissa/anno_files/anno_prom_1k_37k.tsv")
zuse1db_df <- do.call(rbind.data.frame, zuse1db)


#warning - argument is not numeric or logical: returning NA
#try to convert chrom values

zuse1db_df$chrom <- gsub("chr", "", zuse1db_df$chrom)
zuse1db_df$chrom <- as.numeric(zuse1db_df$chrom)
zuse1db_df <- zuse1db_df[complete.cases(zuse1db_df), ]

fwrite(zuse1db_df, "/home/igor/Data/zuse1db/epiclomal/zuse1db_10_test.tsv", sep='\t')

#try to convert chrom values to numeric


replace_chrom <- function(input) {
  
  input$chrom <- gsub("chr", "", input$chrom)
  input$chrom[input$chrom == "Un_random"] <- 30
  input$chrom[input$chrom == "1_random"] <- 31
  input$chrom[input$chrom == "3_random"] <- 33
  input$chrom[input$chrom == "4_random"] <- 34
  input$chrom[input$chrom == "5_random"] <- 35
  input$chrom[input$chrom == "7_random"] <- 37
  input$chrom[input$chrom == "8_random"] <- 38
  input$chrom[input$chrom == "9_random"] <- 39
  input$chrom[input$chrom == "X_random"] <- 40
  input$chrom[input$chrom == "Y_random"] <- 41
  input$chrom[input$chrom == "13_random"] <- 43
  input$chrom[input$chrom == "16_random"] <- 46
  input$chrom[input$chrom == "17_random"] <- 47
  input$chrom[input$chrom == "M"] <- 48
  input$chrom[input$chrom == "X"] <- 49
  input$chrom[input$chrom == "Y"] <- 50
  
  return(input)
  
}

zuse1db$chrom <- as.numeric(zuse1db$chrom)
fwrite(zuse1db, "/home/igor/Data/zuse1db/epiclomal/zuse1db_met_chr_conv.tsv")


#load_data

input_CpG_data_file <- "/home/igor/Data/zuse1db/epiclomal/zuse1db_10_test.tsv"
input_regions_file <- "/home/igor/Data/zuse1db/epiclomal/anno_chr_num.tsv"

function (outdir, input_CpG_data_file, input_regions_file, use_cache) 
{
  cached_data <- file.path(outdir, gsub(".tsv.gz", ".RDa.gz", 
                                        basename(input_CpG_data_file)))
  print(cached_data)
  if (file.exists(cached_data) & use_cache) {
    print("loading cached data")
    load(cached_data)
  }
  else {
    print("reading methylation data from TSV")
    tmp <- read.csv(input_CpG_data_file, sep = "\t", header = TRUE, 
                    check.names = FALSE)
    input_CpG_data <- as.matrix(tmp[, -1])
    rownames(input_CpG_data) <- tmp$cell_id
    rm(tmp)
    print("reading regions from TSV")
    tmp <- read.csv(input_regions_file, sep = "\t", header = TRUE, 
                    check.names = FALSE)
    input_regions <- as.matrix(tmp[, -1]) + 1
    colnames(input_regions) <- c("start", "end")
    rownames(input_regions) <- tmp$region_id
    rm(tmp)
    print("calculating mean_meth_matrix")
    mean_meth_matrix <- t(apply(input_CpG_data, 1, .extract_mean_meth_per_cell, 
                                region_coord = input_regions))
    if (use_cache) {
      print("Saving data to Rda file")
      save(input_CpG_data, input_regions, mean_meth_matrix, 
           file = cached_data, compress = "gzip")
    }
  }
  return(list(input_CpG_data = input_CpG_data, input_regions = input_regions, 
              mean_meth_matrix = mean_meth_matrix))
}


#delete Nth row

Nth_delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]

regions_cut <- Nth_delete(regions_cut, 2)
