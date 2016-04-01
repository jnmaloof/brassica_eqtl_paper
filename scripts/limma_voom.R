##########
# Cody Markelz
# markelz@gmail.com
# voom/limma analysis for brassica RNAseq data
# Image is limma-voom.RData
##########

# one time and will take a while to install all dependencies
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

# load libs
library(edgeR)
library(limma)

# go to data directory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")

# large file (112 MB) takes about 1 minute to get into memory
mapped_counts <- read.delim("RIL_v1.5_mapping.tsv", header = TRUE, sep = "\t")
dim(mapped_counts)
# [1] 43151   843

colnames(mapped_counts)

#replace all NA values with 0 
mapped_counts[is.na(mapped_counts)] <- 0
head(mapped_counts)
tail(mapped_counts)

# remove first row
mapped_counts <- mapped_counts[-1,]
head(mapped_counts)[,1:10]

# move the first column to row names
rownames(mapped_counts) <- mapped_counts[,1]
mapped_counts <- mapped_counts[,-1]
head(mapped_counts)
dim(mapped_counts)

# remove cr trt samples
names(mapped_counts)
mapped_counts <- mapped_counts[,grepl("*\\UN*",names(mapped_counts))]
head(mapped_counts)[1:10]
names(mapped_counts)

# use this cleaned table for analysis
write.table(mapped_counts, file = "RIL_mapped_counts_cleaned.csv", sep=",") 

#######
# design matrix for regression
#######
Br_RIL   <- factor(sub("(RIL_)(\\d+)(_)(\\w+)(_)(Rep)(\\d)+(.)+",
                       "\\2", colnames(mapped_counts)))
Br_RIL  # 122 levels
Br_RIL <- relevel(Br_RIL, ref = "103") # if needed, RIL is abitrary

design <- model.matrix(~ 0 + Br_RIL)
head(design)

######
# edgeR and Limma
######

# make necessary objects
brassica_DE <- DGEList(counts = mapped_counts, group = Br_RIL)
brassica_DE$samples
dim(brassica_DE)
# [1] 43150   434

# keep genes with at least 1 count per million in at least 50 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
# [1] 33353   434

brassica_DE <- calcNormFactors(brassica_DE)
system.time(brass_voom <- voom(brassica_DE, design, plot = FALSE))
system.time(geno_fit <-lmFit(brass_voom, design)) 
   # user  system elapsed 
# 167.291  27.120 193.555 

geno_fit <- eBayes(geno_fit)
toptable(geno_fit)
glmLRT(geno_fit)

# output plots, they are much to large to fit into memory
setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/output")
png(file = "Brassica_MDS.png", width = 1000 , height = 1000, res = 100)
plotMDS(brassica_DE)
dev.off()

#check out coeffs
brassica_BLUEs <- geno_fit$coeff
head(brassica_BLUEs)
str(brassica_BLUEs)
brassica_BLUE_df <- as.data.frame(geno_fit$coeff)
head(brassica_BLUE_df)

plot(brassica_BLUE_df)[1,]
hist(brassica_BLUE_df[,4])

brassica_BLUE_df_t <- as.data.frame(t(brassica_BLUE_df))
head(brassica_BLUE_df_t)[,1:10]

#export BLUES for eQTL mapping
setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
write.table(brassica_BLUE_df, file="brassica_blues.csv", sep=",",
           row.names = TRUE, col.names = TRUE) 

save.image(file = "limma-voom.RData")
# end

