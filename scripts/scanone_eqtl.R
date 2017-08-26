##########
# Cody Markelz
# markelz@gmail.com
# eQTL analysis pipeline for Brassica rapa
# GPL License 
##########

library(qtl)
library(qtlhot)
library(plyr)
setwd("/vol_c/brassica-eqtl-paper/input/")

brass_total <- read.cross("csvsr", genfile ="../../brassica_genetic_map_paper/output/snp_map_rqtl_Mbp_ref1.5.2_cross_output_gen.csv", 
	                       phefile="br_blues_RQTL.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brass_total)

class(brass_total)[1] <- "riself"
brass_total <- jittermap(brass_total)
brass_total
# individuals that were not genotyped. 

brass_total <- est.rf(brass_total)
plot.rf(brass_total) 

#about a minute
brass_total <- calc.errorlod(brass_total, error.prob=0.001)

system.time(scanone_imp_tot <- scanone(brass_total, pheno.col = 1:33353, 
	         method = "imp", use = "all.obs"))
# user  system elapsed 
# 699.787   8.045 705.035 

save.image(file = "scanone-eqtl.RData", version = NULL,
           ascii = FALSE, safe = TRUE)

set.seed(12345)
permtest <- scanone(brass_total, method = "imp", n.perm = 10000)
permtest
alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs <- summary(permtest, alphas)
lod.thrs # 1% cutoff 3.84, 5 to be more conservative

high1 <- highlod(scanone_imp_tot, lod.thr = 5, drop.lod = 0)
# A06x22222087 A06  70   294

#infile genomic coordinates of genes
setwd("/Users/Cody_2/git.repos/brassica_genome_db/raw_data")

transcripts <- read.table("transcripts_eqtl_start_stop_eqtl.csv", sep = ",", header = TRUE)
head(transcripts)
str(transcripts)
dim(transcripts)
# [1] 43463     6 
# includes scaffolds, need to figure this out
transcripts$tx_name <- as.character(transcripts$tx_name)
transcripts$tx_chrom <- as.character(transcripts$tx_chrom)

# it is just easier to use merge to get the type of dataframe that is useful for plotting
# and other applications for now
markers <- as.data.frame(rownames(scanone_imp_tot))
head(markers)
markers[1]
names(markers) <- "marker_name"
markers$index <- rownames(markers)

lods <- as.data.frame(high1$highlod)
head(lods)
dim(lods)
# [1] 38553     3


gene_names <- as.data.frame(high1$names)
head(gene_names)
names(gene_names) <- "tx_name"
gene_names$old <- rownames(gene_names)

gene_lods <- merge(gene_names, lods, by.x = "old", by.y = "phenos", all.y = TRUE)
head(gene_lods)
dim(gene_lods)

gene_lod_mark <- merge(gene_lods, markers, by.x = "row", by.y = "index", all.x = TRUE)
head(gene_lod_mark)
dim(gene_lod_mark)

cistrans_df <- merge(gene_lod_mark, transcripts, by.x = "tx_name", by.y = "tx_name", all.x = TRUE)
head(cistrans_df)
dim(cistrans_df)
tail(cistrans_df)

#load plyr library to do a quick string split
chr_pos <- ldply(strsplit(as.character(cistrans_df$marker_name), split = "x"))
head(chr_pos)
cistrans_df$qtl_chrom <- chr_pos$V1
cistrans_df$qtl_pos <- chr_pos$V2

head(cistrans_df)
str(cistrans_df)

# ignoring the scaffolds for a moment
cistrans_df$cis_trans <- ifelse(cistrans_df$tx_chrom == cistrans_df$qtl_chrom, paste("cis"), paste("trans"))
#take a quick look
head(cistrans_df)
cistrans_df[10000:10100,]
tail(cistrans_df)

# with LOD threshold of 5 and 0 lod support interval
setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
write.table(cistrans_df, "cis_trans_scanone.csv", sep = ",")
save.image(file = "scanone-eqtl.RData", version = NULL, ascii = FALSE, safe = TRUE)
# end



