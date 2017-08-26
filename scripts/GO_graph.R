########
# Cody Markelz
# rjcmarkelz@gmail.com
# GPL License
########
library(qtl)

# need to rerun 
setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
load("QTL-genes_data.RData")
load("parent-DE-data.RData")
load("Parent-QTL.RData")
load("FT-traits.RData")
load("genome-db")
ls()
f# cluttered workspace so here are the key players
brassica_db <- makeTxDb(transcripts, splicings)
brassica_db #transcript database
brass_total # ril gene expression
FT_pheno # FT and yield traits
head(transcripts) # total gene list and locations
head(flr_genes) # known flowering time genes
head(cis_df) # non-aggregated cis 
head(cis_ag) 
dim(cis_ag) #9907
head(cistrans_df) # total cis and trans eQTL