###########
# Cody Markelz
# rjcmarkelz@gmail.com
# Modified March 17, 2016
###########

library(dplyr)
library(data.table)
library(ggplot2)

# load dataset
load('~/git.repos/brassica_eqtl_paper/input/scanone-eqtl.RData')
ls()

# choose between cis-peaks
# trans peaks go to scaffolds
# this is why I need the scaffold data from Mike
dim(cistrans_df)
head(cistrans_df)
tail(cistrans_df)
cistrans_df$qtl_pos <- as.numeric(cistrans_df$qtl_pos)

scaffolds <- cistrans_df[grepl("^Sc", cistrans_df$tx_chrom),]
dim(scaffolds)
head(scaffolds)

cis_df <- subset(cistrans_df, cis_trans == "cis")
dim(cis_df)
# [1] 26095    13
head(cis_df)
length(unique(cis_df$tx_name))
# [1] 8660
str(cis_df)
cis_df$qtl_pos <- as.numeric(cis_df$qtl_pos)

########
########
# figure out trans genes on same chromosome as gene is located
# if on same chrom AND within X bp difference then cis, else trans
########
########

plot(cis_df$lod)
dim(cis_df)
head(cis_df)

cis_df$gene_qtl_diff <- abs(cis_df$qtl_pos - cis_df$tx_start)
str(cis_df)
plot(cis_df$gene_qtl_diff)
hist(cis_df$gene_qtl_diff)
with(cis_df[cis_df$gene_qtl_diff < 100000,], plot(gene_qtl_diff))
with(cis_df[(cis_df$gene_qtl_diff > 100000) & (cis_df$qtl_chrom == "A10"),], plot(tx_start, gene_qtl_diff))

cis_plot <- ggplot(cis_df)
cis_plot <- cis_plot +  theme_bw() + geom_point(aes(x = gene_qtl_diff, y = lod), size = 1.5) +
                        facet_grid(qtl_chrom ~ . ) +
                        theme(text = element_text(size = 20))
cis_plot

?aggregate
cis_ag <- aggregate(lod ~ tx_name, data = cis_df, max)
dim(cis_ag)
head(cis_df)
cis_df$tx_start_Mbp <- cis_df$tx_start/1000000
cis_df$diff_Mbp <- cis_df$gene_qtl_diff/1000000

cis_df_ag <- merge(cis_ag, cis_df, by = c("tx_name", "lod"))
with(cis_df_ag[(cis_df_ag$gene_qtl_diff > 100000) & (cis_df_ag$qtl_chrom == "A10"),], plot(tx_start, gene_qtl_diff))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A03"),], plot(diff_Mbp, lod))



with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A08"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A01"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A02"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A03"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A04"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A05"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A06"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A07"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A08"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A09"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A10"),], plot(diff_Mbp, lod))


with(cis_df[cis_df$gene_qtl_diff < 100000,], plot(gene_qtl_diff))
with(cis_df[(cis_df$gene_qtl_diff > 100000) & (cis_df$qtl_chrom == "A10"),], plot(tx_start, gene_qtl_diff))


# end of A03
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A03"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A03"),], plot(tx_start_Mbp, diff_Mbp))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A09"),], plot(tx_start_Mbp, diff_Mbp))


with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A04"),], plot(diff_Mbp, lod))
with(cis_df[(cis_df$diff_Mbp > 0.5) & (cis_df$qtl_chrom == "A01"),], plot(tx_start_Mbp, diff_Mbp))

dim(cis_df)
head(cis_df)
length(unique(cis_df$tx_name))

dim(cis_df[duplicated(cis_df$tx_name),])
# 353
?duplicated
dup <- cis_df[duplicated(cis_df$tx_name),]
dim(dup)

dim(cis_df)
head(cis_df)
test <- aggregate(gene_qtl_diff ~ tx_name, data = cis_df, max)
dim(test)


cis_plot <- ggplot(cis_df)
cis_plot <- cis_plot +  theme_bw() + geom_point(aes(x = tx_start, y = lod), size = 1.5) +
                        facet_grid(qtl_chrom ~ . ) +
                        theme(text = element_text(size = 20))
cis_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_eqtl_plot.png", width = 10, height = 15)

#get large cis effect genes
# arbitrary cutoff
large_cis <- cis_df[cis_df$lod > 100,]
dim(large_cis)
large_cis
#79
setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(large_cis, "large_effect_cis.csv", sep = ",", col.names = TRUE, row.names = TRUE)

#######
# trans eQTL
#######
trans_df <- subset(cistrans_df, cis_trans == "trans")

dim(trans_df)
trans_df <- trans_df[!grepl("^Sc", trans_df$tx_chrom),]
dim(trans_df)
# [1] 11024    13

head(trans_df)
tail(trans_df)
with(trans_df[(cis_df$qtl_chrom == "A01"),], plot(qtl_pos, lod))


out <- subset(trans_df, trans_df$tx_chrom == trans_df$qtl_chrom)
dim(out) # double check trans

large_trans <- trans_df[trans_df$lod > 100,]
dim(large_trans)
large_trans
#20
setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(large_trans, "large_effect_trans.csv", sep = ",", col.names = TRUE, row.names = TRUE)

# we really do not have the resolution to care about how EXACT the qtl is next
# to the start site for the cis so we should just choose largest lod
# that is the statistical information that we have for the population

str(trans_df)
trans_ag <- aggregate(lod ~ tx_name, data = trans_df, max)
dim(trans_ag)
head(trans_ag)
trans_df <- merge(trans_ag, trans_df, by = c("tx_name", "lod"))
with(trans_df[(cis_df$qtl_chrom == "A01"),], plot(qtl_pos, lod))

dim(trans_df)

head(trans_df)
trans_df$qtl_pos <- as.numeric(trans_df$qtl_pos)

trans_plot <- ggplot(trans_df)
trans_plot <- trans_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
trans_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("trans_eqtl_plot.png", width = 10, height = 10)

head(trans_df)
head(cis_df)
str(cis_df)
str(trans_df)

########
#Trans Counts
########
head(trans_df) # post aggregation for largest eQTL value
str(trans_df)
trans_df$marker_name <- as.character(trans_df$marker_name)

#table will count each occurance
trans_counts <- as.data.frame(table(trans_df$marker_name))
head(trans_counts)
str(trans_counts)
trans_counts$chr <- sub("(A\\d+)(x)(\\d+)", "\\1", trans_counts$Var1)
trans_counts$pos <- as.numeric(sub("(A\\d+)(x)(\\d+)", "\\3", trans_counts$Var1))
trans_counts$pos <- trans_counts$pos/1000000

trans_plot <- ggplot(trans_counts)
trans_plot <- trans_plot +  theme_bw() + geom_point(aes(x = pos, y = Freq), size = 1) +
                        facet_grid(chr ~ . )
                        theme(text = element_text(size = 20))
trans_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("trans_hotspot.png", width = 10, height = 10)

# put data together
ct_merge <- rbind(cis_df, trans_df)
head(ct_merge)
dim(ct_merge)
str(ct_merge)

merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot +  theme_bw() + geom_point(aes(x = qtl_pos, y = lod, color = cis_trans), size = 2) +
                        facet_grid(qtl_chrom ~ . )
                        theme(text = element_text(size = 20))
merge_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_trans_eqtl_plot.png", width = 10, height = 10)

# cis trans plot
merge_plot <- ggplot(ct_merge)
merge_plot <- merge_plot + geom_point(aes(x = qtl_pos, y = tx_start, color = lod ), size = 1.5) +
                        scale_y_reverse() +
                        facet_grid(tx_chrom ~ qtl_chrom) + theme_bw() + 
                        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_blank())
merge_plot
?ggsave
ggsave("cis_diagonal.png", width = 10, height = 10, dpi = 300)

############
# flowering gene distribution
############

setwd('~/git.repos/brassica_eqtl_v1.5/data/')
br_flr <- read.delim("br_flowering_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_flr)
dim(br_flr)
colnames(br_flr)

str(flr_genes)
flr_genes <- br_flr[9]
colnames(flr_genes)[1] <- "gene_name"
head(flr_genes)

unique(br_flr$Gene.name.in.B..rapa)

trans_df$tx_name <- as.character(trans_df$tx_name)
cis_df$tx_name <- as.character(cis_df$tx_name)

trans_flr_df <- trans_df[trans_df$tx_name %in% flr_genes$gene_name,]
head(trans_flr_df)
dim(trans_flr_df)
trans_flr_df

trans_flr_plot <- ggplot(trans_df)
trans_flr_plot <- trans_flr_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_flr_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red', size = 1) +
                theme(text = element_text(size = 20))
trans_flr_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("trans_flr_plot.png", width = 10, height = 10)

cis_flr_df <- cis_df[cis_df$tx_name %in% flr_genes$gene_name,]
head(cis_flr_df)
dim(cis_flr_df)
sub <- subset(cis_flr_df, tx_chrom == "A06")
sub
flr_A10_eqtl <- subset(cis_flr_df, qtl_chrom = "A10")
head(flr_A10_eqtl)
dim(flr_A10_eqtl)
?sort
sort(flr_A10_eqtl, flr_A10_eqtl$lod)
flr_A10_eqtl
setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(flr_A10_eqtl, "flowering_eQTL.csv", sep = ",")

head(cis_df)

flr_can <- as.data.frame(subset(cis_df, cis_df$qtl_chr == "A10" & cis_df$qtl_pos > 12093873 & cis_df$qtl_pos < 13111909))
flr_can
dim(flr_can)

setwd('~/git.repos/brassica_eqtl_paper/input/')
write.table(flr_can, "flowering_eQTL_candidates.csv", sep = ",")


head(trans_df)
dim(cis_flr_df)

cis_flr_plot <- ggplot(cis_df)
cis_flr_plot <- cis_flr_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_flr_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_flr_plot
setwd('~/git.repos/brassica_eqtl_paper/output/')
ggsave("cis_flr_plot.png", width = 10, height = 10)

############
# parental analysis data
############
setwd('~/git.repos/brassica_parents/data/')

br_fruits <- read.delim("parental_fruit_field_DE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_fruits)
tail(br_fruits)
dim(br_fruits)
br_fruits$tx_name <- rownames(br_fruits)

cis_de_df <- cis_df[cis_df$tx_name %in% br_fruits$tx_name,]
head(cis_de_df)
dim(cis_de_df)
dim(cis_df)

trans_de_df <- trans_df[trans_df$tx_name %in% br_fruits$tx_name,]
head(trans_de_df)
dim(trans_de_df)
dim(trans_df)

cis_de_plot <- ggplot(cis_de_df)
cis_de_plot <- cis_de_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_de_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_de_plot

ggsave("cis_de_fruit_plot.png", width = 10, height = 10)

trans_de_plot <- ggplot(trans_de_df)
trans_de_plot <- trans_de_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_plot
ggsave("trans_de_fruit_plot.png", width = 10, height = 10)


br_leaf <- read.delim("parental_leaf_field_DE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(br_leaf)
dim(br_leaf)
br_leaf$tx_name <- rownames(br_leaf)


cis_de_leaf_df <- cis_df[cis_df$tx_name %in% br_leaf$tx_name,]
head(cis_de_leaf_df)
dim(cis_de_leaf_df)
dim(cis_df)

trans_de_leaf_df <- trans_df[trans_df$tx_name %in% br_leaf$tx_name,]
head(trans_de_leaf_df)
dim(trans_de_leaf_df)
dim(trans_df)

cis_de_leaf_plot <- ggplot(cis_df)
cis_de_leaf_plot <- cis_de_leaf_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = cis_de_leaf_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
cis_de_leaf_plot
ggsave("cis_de_leaf_plot.png", width = 10, height = 10)

trans_de_leaf_plot <- ggplot(trans_df)
trans_de_leaf_plot <- trans_de_leaf_plot + theme_bw() + geom_point(aes(x = qtl_pos, y = lod), size = 1.5) +
                facet_grid(qtl_chrom ~ . ) + 
                geom_segment(data = trans_de_leaf_df, aes(x = qtl_pos, xend = qtl_pos), y = 0 , yend = 100, color = 'red') +
                theme(text = element_text(size = 20))
trans_de_leaf_plot
ggsave("trans_de_leaf_plot.png", width = 10, height = 10)

setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
save.image(file = "un_eqtl_parent_field.RData", version = NULL, ascii = FALSE, safe = TRUE)

# end