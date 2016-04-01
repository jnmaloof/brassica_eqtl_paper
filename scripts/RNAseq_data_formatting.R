##########
# Cody Markelz
# markelz@gmail.com
# formatting for RQTL
# Image is RNAseq-data-formatting.RData
# GPL License- http://www.gnu.org/licenses/gpl-3.0.en.html
##########

setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
br_blues <- read.delim("brassica_blues.csv", 
	                   header = TRUE, row.names = 1, sep = ",")

#transpose
br_blues_t <- as.data.frame(t(br_blues))
head(br_blues_t)[,1:10]

#add ril numbers for RQTL
br_blues_t$id <- sub("(Br_group)(\\d+)(_)(UN)", "RIL_\\2", row.names(br_blues_t))
dim(br_blues_t)
head(br_blues_t)[,33350:33354]

#transpose back
#takes a bit to transpose back in this direction
br_blues_final <- as.data.frame(t(br_blues_t))

head(br_blues_final)[,1:10]
tail(br_blues_final)[,1:10]

names(br_blues_final) <- sub("(Br_)(RIL)(\\d+)", "\\2_\\3", names(br_blues_final))
names(br_blues_final)

#save output
setwd("/Users/Cody_2/git.repos/brassica_eqtl_paper/input")
write.table(br_blues_final, "br_blues_RQTL.csv", col.names= FALSE, 
	         row.names = TRUE, sep = ",")
save.image(file = "RNAseq-data-formatting.RData")
# end
