## make insertion plot diagram to compare insertions at TA sites in mbovis and mtb

# import insertion data (reads per insertion site)

mtb_temp<-read.table("dejesus_combined_ttr.wig", header=F, sep="\t",
                     comment.char="#", stringsAsFactors = F)
mtb_df<- mtb_temp[,1:15]
read_sums<- rowSums(mtb_df[,2:15])
read_means<- rowMeans(mtb_df[,2:15])
mtb_df<- cbind(mtb_df, read_sums, read_means)
head(mtb_df)

bovis_temp<-read.table("bovis_combined_TTR.wig", header=F, sep="\t",
                       comment.char="#", stringsAsFactors = F)
head(bovis_temp)
bovis_df <- bovis_temp[,1:4]
bovis_df <- cbind(bovis_df, rowSums(bovis_df[,2:4]), rowMeans(bovis_df[,2:4]))
head(bovis_df)
colnames(bovis_df)<-c("pos", "B", "C", "G", "sum", "mean")

max(bovis_df$sum)
#90726.9
max(bovis_df$mean)
#30242.3

max(mtb_df$read_sums)
#946019.6
max(mtb_df$read_means)
#67572.83