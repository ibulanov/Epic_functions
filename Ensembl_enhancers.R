#9/02/21

#Downoad enhancers from Ensembl Regulatory Build for mm10 #3 run

library(biomaRt)
listMarts()
ensembl <- useMart("ensembl")

datasets <- listDatasets(ensembl)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

listEnsembl()

listEnsembl("GRCh=37")

regulation = useEnsembl(biomart="regulation", GRCh = 37)
head(listDatasets(regulation))

regulation = useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh=37)

head(listAttributes(regulation))

head(listFilters(regulation))

result = getBM(attributes = c('chromosome', 'entrezgene_id'),
               filters = 'affy_hg_u133_plus_2',
               values = affyids, 
               mart = ensembl)

#Another approach from ensembl ftp server

gm_enh <- fread("/home/igor/Data/ensembl/human_GRCh17/release_102/homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20191101.gff.gz")
fibr_enh <- fread("/home/igor/Data/ensembl/human_GRCh17/release_102/homo_sapiens.GRCh37.dermal_fibroblast.Regulatory_Build.regulatory_activity.20191101.gff.gz")

gm_fibr_enh <- rbind(gm_enh, fibr_enh)
gm_fibr_enh <- gm_fibr_enh[V3 == "enhancer"]

test = strsplit(gm_fibr_enh$V9,";")

library(magrittr)
gm_fibr_enh$V11 <- strsplit(gm_fibr_enh$V9, ";") %>% sapply(extract2, 1)

all_enh <- gm_fibr_enh[,c(1,4,5)]

act_enh <- gm_fibr_enh[V11 == "activity=ACTIVE"]
act_enh <- act_enh[,c(1,4,5)]

#create table of ERB
gm_all <- fread("/home/igor/Data/ensembl/human_GRCh17/release_102/homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20191101.gff.gz")
fibro_all <- fread("/home/igor/Data/ensembl/human_GRCh17/release_102/homo_sapiens.GRCh37.dermal_fibroblast.Regulatory_Build.regulatory_activity.20191101.gff.gz")

all_enh <- rbind(gm_all, fibro_all) %>% .[,"V3":= "enhancer"]

all_enh$V11 <- strsplit(all_enh$V9, ";") %>% sapply(extract2, 1)

#active enh - too much 
act_enh <- all_enh[V11 == "activity=ACTIVE"]
act_enh <- act_enh[,c(1,4,5)]
act_enh <- distinct(act_enh)

#active enh for 2 cell type separately
gm_active <- gm_all[V3 == "enhancer"]
gm_active$V11 <- strsplit(gm_active$V9, ";") %>% sapply(extract2, 1)
gm_active <- gm_active[V11 == "activity=ACTIVE"]

fibro_active <- fibro_all[V3 == "enhancer"]
fibro_active$V11 <- strsplit(fibro_active$V9, ";") %>% sapply(extract2, 1)
fibro_active <- fibro_active[V11 == "activity=ACTIVE"]

gm_fibro_active <- rbind(gm_active, fibro_active)
gm_fibro_active <- gm_fibro_active[,c(1,4,5)]
gm_fibro_active <- distinct(gm_fibro_active, .keep_all = TRUE)

#take another type of enhancers
gm_enh <- gm_all[V3 == "enhancer"]
fb_enh <- fibro_all[V3 == "enhancer"]

gm_enh$V11 <- strsplit(gm_enh$V9, ";") %>% sapply(extract2, 1)
fb_enh$V11 <- strsplit(fb_enh$V9, ";") %>% sapply(extract2, 1)

#take INACTIVE type of enhancers
gm_inact <- gm_enh[V11 == "activity=INACTIVE"] %>% .[, c("V1", "V4", "V5")]
fb_inact <- fb_enh[V11 == "activity=INACTIVE"] %>% .[, c("V1", "V4", "V5")]

enh_inact <- rbind(gm_inact, fb_inact)
enh_inact <- distinct(enh_inact)

setkey(enh_inact)
enh_inact$V1 <- paste0("chr", enh_inact$V1)
colnames(enh_inact) <- c("chr", "start", "end")
fwrite(enh_inact, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/gm_fibro_enh_inactive.tsv", sep='\t')

#take POISED/REPRESSEED type of enhancers
gm_type <- gm_enh[V11 == "activity=REPRESSED"] %>% .[, c("V1", "V4", "V5")]
fb_type <- fb_enh[V11 == "activity=REPRESSED"] %>% .[, c("V1", "V4", "V5")]

enh_type <- rbind(gm_type, fb_type)
enh_type <- distinct(enh_type)

setkey(enh_type)
enh_type$V1 <- paste0("chr", enh_type$V1)
colnames(enh_type) <- c("chr", "start", "end")
fwrite(enh_type, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/gm_fibro_enh_repressed.tsv", sep='\t')

#take another regulatory elements - all of them

all_CTCF_BS <- rbind(gm_all, fibro_all) %>% .[V3 == "CTCF_binding_site"]
all_OCR <- rbind(gm_all, fibro_all) %>% .[V3 == "open_chromatin_region"]
all_prom <- rbind(gm_all, fibro_all) %>% .[V3 == "promoter"]
all_prom_FR <- rbind(gm_all, fibro_all) %>% .[V3 == "promoter_flanking_region"]
all_TFBS <- rbind(gm_all, fibro_all) %>% .[V3 == "TF_binding_site"]

all_CTCF_BS$V11 <- strsplit(all_CTCF_BS$V9, ";") %>% sapply(extract2, 1)
all_OCR$V11 <- strsplit(all_OCR$V9, ";") %>% sapply(extract2, 1)
all_prom$V11 <- strsplit(all_prom$V9, ";") %>% sapply(extract2, 1)
all_prom_FR$V11 <- strsplit(all_prom_FR$V9, ";") %>% sapply(extract2, 1)
all_TFBS$V11 <- strsplit(all_TFBS$V9, ";") %>% sapply(extract2, 1)

#separate each RE by activity
all_RE_active <- all_TFBS[V11 == "activity=ACTIVE"] %>% .[, c("V1", "V4", "V5")] %>% setkey(.)
all_RE_active$V1 <- paste0("chr", all_RE_active$V1)
colnames(all_RE_active) <- c("chr", "start", "end")
all_RE_active <- distinct(all_RE_active)
fwrite(all_RE_active, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/TFBS/OCR_active.tsv", sep='\t')

all_RE_inactive <- all_TFBS[V11 == "activity=INACTIVE"] %>% .[, c("V1", "V4", "V5")]
all_RE_inactive$V1 <- paste0("chr", all_RE_inactive$V1)
colnames(all_RE_inactive) <- c("chr", "start", "end")
all_RE_inactive <- distinct(all_RE_inactive)
fwrite(all_RE_inactive, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/TFBS/CTCF_BS_inactive.tsv", sep='\t')

all_RE_poised <- all_TFBS[V11 == "activity=POISED"] %>% .[, c("V1", "V4", "V5")]
all_RE_poised$V1 <- paste0("chr", all_RE_poised$V1)
colnames(all_RE_poised) <- c("chr", "start", "end")
all_RE_poised <- distinct(all_RE_poised)
fwrite(all_RE_poised, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/TFBS/CTCF_BS_poised.tsv", sep='\t')

all_RE_repr <- all_TFBS[V11 == "activity=REPRESSED"] %>% .[, c("V1", "V4", "V5")]
all_RE_repr$V1 <- paste0("chr", all_RE_repr$V1)
colnames(all_RE_repr) <- c("chr", "start", "end")
all_RE_repr <- distinct(all_RE_repr)
fwrite(all_RE_repr, "/home/igor/Data/ensembl/human_GRCh17/release_102/new_ERB/TFBS/CTCF_BS_repr.tsv", sep='\t')

#lists


#to add data to count table
count_table <- data.table(ERB = rep(NA, 21), activity = rep(NA, 21), count = rep(NA, 21), sum_length = rep(NA, 21))
count_table$ERB <- c(rep("enhancer", 4), rep("CTCF_binding_site", 4), rep("open_chromatin_region", 1),
                     rep("promoter", 4), rep("promoter_flanking_region", 4),rep("TF_binding_site", 4))

count_table$activity[18:21] <- c("ACTIVE", "INACTIVE", "POISED", "REPRESSED")

#
count_table$count[18:21] <- c(1666, 12040, 411, 2137) 

count_table$sum_length[18:21] <- c(sum(all_RE_active$end - all_RE_active$start),
                                 sum(all_RE_inactive$end - all_RE_inactive$start),
                                 sum(all_RE_poised$end - all_RE_poised$start),
                                 sum(all_RE_repr$end - all_RE_repr$start))

