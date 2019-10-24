
#manually removed the "#" in front the chromosome column name before reading in files

#read in ambient DMR table
amb_DMR <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20191024_10bp/amb_AllTimes_DMR250bp_MCmax10_cov5x_rms_results_collapsed.tsv", header = TRUE,sep = "\t",stringsAsFactors = FALSE)
#loop through table and keep only lines where up to one sample contains NA for % methylation

df <- data.frame() #create empty data frame to bind filtered rows into
for(i in (1:nrow(amb_DMR))){
  day0 <- amb_DMR[i,7:10] #define columns from the category Day 0
  day10 <- amb_DMR[i,11:14] #define columns from the category Day 10
  day135 <- amb_DMR[i,15:18] #define columns from the category Day 135
  day145 <- amb_DMR[i,19:22] #define columns from the category Day 145
  if(length(which(is.na(day0))) < 2 & length(which(is.na(day10))) < 2 & length(which(is.na(day135))) < 2 & length(which(is.na(day145))) < 2){
    df <- rbind(df,amb_DMR[i,]) #conditional statement: if less than 2 sameples/category have NA for % methylation bind the whole row to the new dataframe
  }
}
#write the output
write.table(df,"amb_AllTimes_DMR250bp_MCmax10_cov5x_rms_results_filtered.tsv", sep = "\t", row.names= FALSE, quote = FALSE)

#do the same for day 10

day10_DMR <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20191024_10bp/day10_AllpH_DMR250bp_MCmax10_cov5x_rms_results_collapsed.tsv", header = TRUE,sep = "\t",stringsAsFactors = FALSE)
#loop through table and keep only lines where up to one sample contains NA for % methylation

df <- data.frame() #create empty data frame to bind filtered rows into
for(i in (1:nrow(day10_DMR))){
  LowpH <- day10_DMR[i,c(7,8,15,16)] #define columns from the category low pH
  SuperLowpH <- day10_DMR[i,c(9,10,17,18)] #define columns from the category super low pH
  amb <- day10_DMR[i,c(11:14)] #define columns from the category ambient
  if(length(which(is.na(LowpH))) < 2 & length(which(is.na(SuperLowpH))) < 2 & length(which(is.na(amb))) < 2){
    df <- rbind(df,day10_DMR[i,]) #conditional statement: if less than 2 sameples/category have NA for % methylation bind the whole row to the new dataframe
  }
}
#write the output
write.table(df,"day10_AllpH_DMR250bp_MCmax10_cov5x_rms_results_filtered.tsv", sep = "\t", row.names= FALSE, quote = FALSE)

#do the same for day 135

day135_DMR <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20191024_10bp/day135_AllpH_DMR250bp_MCmax10_cov5x_rms_results_collapsed.tsv", header = TRUE,sep = "\t",stringsAsFactors = FALSE)
#loop through table and keep only lines where up to one sample contains NA for % methylation

df <- data.frame() #create empty data frame to bind filtered rows into
for(i in (1:nrow(day135_DMR))){
  amb <- day135_DMR[i,7:10] #define columns from the category low pH
  LowpH <- day135_DMR[i,11:14] #define columns from the category super low pH
  SuperLowpH <- day135_DMR[i,15:18] #define columns from the category ambient
  if(length(which(is.na(LowpH))) < 2 & length(which(is.na(SuperLowpH))) < 2 & length(which(is.na(amb))) < 2){
    df <- rbind(df,day135_DMR[i,]) #conditional statement: if less than 2 sameples/category have NA for % methylation bind the whole row to the new dataframe
  }
}
#write the output
write.table(df,"day135_AllpH_DMR250bp_MCmax10_cov5x_rms_results_filtered.tsv", sep = "\t", row.names= FALSE, quote = FALSE)

#do the same for day 145

day145_DMR <- read.table("/Volumes/web/metacarcinus/Pgenerosa/analyses/20191024_10bp/day145_AllpH_DMR250bp_MCmax10_cov5x_rms_results_collapsed.tsv", header = TRUE,sep = "\t",stringsAsFactors = FALSE)
#loop through table and keep only lines where up to one sample contains NA for % methylation

df <- data.frame() #create empty data frame to bind filtered rows into
for(i in (1:nrow(day145_DMR))){
  AmbAmb <- day145_DMR[i,9:12] #define columns from the category low pH
  LowAmb <- day145_DMR[i,c(7,8,15,16)] #define columns from the category super low pH
  SuperLowAmb <- day145_DMR[i,c(13,14,17,18)] #define columns from the category ambient
  AmbLow <- day145_DMR[i,c(19,20,27,28)] #define columns from the category ambient
  LowLow <- day145_DMR[i,c(21,22,29,30)] #define columns from the category ambient
  SuperLowLow <- day145_DMR[i,23:26] #define columns from the category ambient
  if(length(which(is.na(AmbAmb))) < 2 & length(which(is.na(LowAmb))) < 2 & length(which(is.na(SuperLowAmb))) < 2 & length(which(is.na(AmbLow))) < 2 & length(which(is.na(LowLow))) < 2 & length(which(is.na(SuperLowLow))) < 2){
    df <- rbind(df,day145_DMR[i,]) #conditional statement: if less than 2 sameples/category have NA for % methylation bind the whole row to the new dataframe
  }
}
#write the output
write.table(df,"day145_AllpH_DMR250bp_MCmax10_cov5x_rms_results_filtered.tsv", sep = "\t", row.names= FALSE, quote = FALSE)

