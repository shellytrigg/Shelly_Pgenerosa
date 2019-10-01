#load libraries
library(ggplot2)
library(heatmap3)
library(gplots)

#read in data

sep26_DMRs_allAmb <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20190926/amb_AllTimes_DMR250bp_MCmax25_cov5x_clst3_rms_results_collapsed.tsv", header = TRUE, sep = "\t")
sep26_DMRs_day10 <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20190926/day10_AllpH_DMR250bp_MCmax25_cov5x_clst3_rms_results_collapsed.tsv", header = TRUE, sep = "\t")
sep26_DMRs_day135 <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20190926/day135_AllpH_DMR250bp_MCmax25_cov5x_clst3_rms_results_collapsed.tsv", header = TRUE, sep = "\t")
sep26_DMRs_day145 <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20190926/day145_AllpH_DMR250bp_MCmax25_cov5x_clst3_rms_results_collapsed.tsv", header = TRUE, sep = "\t")

#make uniq ID column

sep26_DMRs_allAmb$ID <- paste(sep26_DMRs_allAmb$chr,":",sep26_DMRs_allAmb$start,"-",sep26_DMRs_allAmb$end, sep = "")
sep26_DMRs_allAmb$ID <- gsub("__.*__.*:",":",sep26_DMRs_allAmb$ID)

sep26_DMRs_day10$ID <- paste(sep26_DMRs_day10$chr,":",sep26_DMRs_day10$start,"-",sep26_DMRs_day10$end, sep = "")
sep26_DMRs_day10$ID <- gsub("__.*__.*:",":",sep26_DMRs_day10$ID)

sep26_DMRs_day135$ID <- paste(sep26_DMRs_day135$chr,":",sep26_DMRs_day135$start,"-",sep26_DMRs_day135$end, sep = "")
sep26_DMRs_day135$ID <- gsub("__.*__.*:",":",sep26_DMRs_day135$ID)

sep26_DMRs_day145$ID <- paste(sep26_DMRs_day145$chr,":",sep26_DMRs_day145$start,"-",sep26_DMRs_day145$end, sep = "")
sep26_DMRs_day145$ID <- gsub("__.*__.*:",":",sep26_DMRs_day145$ID)
#create matrix for heatmap

allAmb_m <- as.matrix(sep26_DMRs_allAmb[,7:22])
rownames(allAmb_m) <- sep26_DMRs_allAmb$ID

ColSideColors <- cbind(day = c(rep("#969696",4),rep("#737373",4),rep("#525252",4),rep("#252525",4)))
pdf("~/Documents/GitHub/Shelly_Pgenerosa/analyses/DMRs_heatmap/allAmb_DMR_heatmap.pdf")
heatmap.2(allAmb_m,cexRow = 0.25, cexCol = 0.8,ColSideColors = ColSideColors, Colv=NA, col = bluered, na.color = "black", density.info = "none", trace = "none", scale = "row")
dev.off()

#subset out the data and order it
day10_m <- as.matrix(sep26_DMRs_day10[,c(11:14,7:8,15:16,9:10,17:18)])
rownames(day10_m) <- sep26_DMRs_day10$ID

ColSideColors <- cbind(pH = c(rep("cyan",4),rep("plum2",4),rep("magenta",4)))
pdf("~/Documents/GitHub/Shelly_Pgenerosa/analyses/DMRs_heatmap/day10_DMR_heatmap.pdf")
heatmap.2(day10_m,cexRow = 0.25, cexCol = 0.8,ColSideColors = ColSideColors, Colv=NA, col = bluered, na.color = "black", density.info = "none", trace = "none", scale = "row")
dev.off()

#subset out the data and order it
day135_m <- as.matrix(sep26_DMRs_day135[,c(7:18)])
rownames(day135_m) <- sep26_DMRs_day135$ID

ColSideColors <- cbind(pH = c(rep("cyan",4),rep("plum2",4),rep("magenta",4)))
pdf("~/Documents/GitHub/Shelly_Pgenerosa/analyses/DMRs_heatmap/day135_DMR_heatmap.pdf")
heatmap.2(day135_m,cexRow = 0.25, cexCol = 0.8,ColSideColors = ColSideColors, Colv=NA, col = bluered, na.color = "black", density.info = "none", trace = "none", scale = "row")
dev.off()

#subset out the data and order it
day145_m <- as.matrix(sep26_DMRs_day145[,c(9:12,7:8,15:16,13:14,17:18,19:20,27:28,21:22,29:30,23:26)])
rownames(day145_m) <- sep26_DMRs_day145$ID

ColSideColors <- cbind(pH = c(rep("cyan",4),rep("green1",4),rep("green3",4),rep("mediumpurple1",4),rep("plum2",4),rep("magenta",4)))
pdf("~/Documents/GitHub/Shelly_Pgenerosa/analyses/DMRs_heatmap/day145_DMR_heatmap.pdf")
heatmap.2(day145_m,cexRow = 0.25, cexCol = 0.8,ColSideColors = ColSideColors, Colv=NA, col = bluered, na.color = "black", density.info = "none", trace = "none", scale = "row")
dev.off()
