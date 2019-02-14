#read in bed file with coverage coordinates 
EPI103.cov <- read.table("/Users/Shelly/EPI-103_S27_L005_dedup.W10000.bed", header = FALSE, sep = "\t")

#rename columns
colnames(EPI103.cov) <- c("chr_name","start", "stop", "features", "num.bases", "window.size", "frac.cov")

library(ggplot2)

ggplot(EPI103.cov[which(EPI103.cov$chr_name == "PGA_scaffold1__77_contigs__length_89643857"),], aes(x = start, y = frac.cov)) + geom_density()

#histogram of 
ggplot(EPI103.cov, aes(x = frac.cov)) + geom_histogram(bins = 20, color = "darkblue", fill = "lightblue") + stat_bin(bins = 20, aes(y=..count.., label=..count..), geom="text", vjust=-.1) + theme_bw() + xlab("fraction of coverage") + ylab("number of unique 10K windows") 



#total number of bases not covered
sum(EPI103.cov[which(EPI103.cov$frac.cov == 0),"window.size"])
[1] 833012981
#number of bases in windows covered (not all bases are actually covered but windows are)
sum(EPI103.cov[which(EPI103.cov$frac.cov != 0),"window.size"])
[1] 1372675707
#total number of bases covered
sum(EPI103.cov[which(EPI103.cov$num.bases != 0),"num.bases"])
[1] 558963312

558963312/833012981
67% of all bases are covered by reads
833012981/1372675707

#number of bases in each chromosome
xsum <- rowsum(EPI103.cov$window.size, EPI103.cov$chr_name)
#number of bases in each chromosome covered
ysum <- rowsum(EPI103.cov$num.bases, EPI103.cov$chr_name)

#merge rowsums files
zsum <- data.frame(cbind(xsum, ysum))
zsum$frac.cov <- zsum$X2/zsum$X1
zsum$X3 <- zsum$X1-zsum$X2
colnames(zsum) <- c("chr_length_bp", "chr_cov_bp", "frac.cov","chr_notcov_bp" )




#histogram of coverage x chromosome
ggplot(zsum, aes(x = frac.cov)) + geom_histogram(bins = 20, color = "darkblue", fill = "lightblue") + stat_bin(bins = 20, aes(y=..count.., label=..count..), geom="text", vjust=-.1) + theme_bw() + xlab("fraction of coverage") + ylab("chromosomes") 

#histogram of chromosome x size
ggplot(zsum, aes(x = chr_length_bp)) + geom_histogram(bins = 20, color = "darkblue", fill = "lightblue") + stat_bin(bins = 20, aes(y=..count.., label=..count..), geom="text", vjust=-.1) + theme_bw() + xlab("chromosome size (bp)") + ylab("number of chromosomes") 


#histogram of chromosome x size for chromosomes < 10000 bp
length(zsum$X1 < 10000)
#[1] 313649
ggplot(zsum[which(zsum$chr_length_bp < 10000),], aes(x = chr_length_bp)) + geom_histogram(bins = 20, color = "darkblue", fill = "lightblue") + stat_bin(bins = 20, aes(y=..count.., label=..count..), geom="text", vjust=-.1) + theme_bw() + xlab("chromosome size (bp)") + ylab("number of chromosomes") 

ggplot(zsum[which(zsum$chr_length_bp <=10000),], aes(x = frac.cov)) + geom_histogram()

ggplot(zsum[which(zsum$chr_length_bp >10000),], aes(x = frac.cov)) + geom_histogram()


#percent covereage of chromosomes < 10000 bases
ggplot(zsum[which(zsum$chr_length_bp < 10000),], aes(x = frac.cov)) + geom_histogram(bins = 20, color = "darkblue", fill = "lightblue") + stat_bin(bins = 20, aes(y=..count.., label=..count..), geom="text", vjust=-.1) + theme_bw() + xlab("fraction covered") + ylab("number of chromosomes") 
ggplot(zsum[which(zsum$chr_length_bp < 10000),], aes(x = chr_length_bp, y = frac.cov)) + geom_point()