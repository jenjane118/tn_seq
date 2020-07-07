

sample_names<-read.table("names_wig.txt",stringsAsFactors = F)

#for (j in 1:nrow(sample_names)){
  #insertion_table<-read.csv(paste(sample_names[j,1],".wig",sep = ""))
  #non_norm<-as.data.frame(matrix(0, nrow = nrow(insertion_table), ncol = 3))
  
# If .txt tab file, use this
#B_data <- read.delim(file.choose("B_R1.wig"))  

## non normalised

B__table<-read.csv("B_R1.wig", sep = " ", skip=2)   # positions start on 3rd line
D_table<-read.csv("D_R1.wig", sep = " ", skip=2)
G_table<-read.csv("G_R1.wig", sep = " ", skip=2)

pos_vector<-B_table[,1]
pos_vector

non_norm<-as.data.frame(matrix(0, nrow = nrow(B_table), ncol =3), row.names <-pos_vector)
colnames(non_norm)<-c('B', 'D', 'G')
non_norm[,1]<-B_table[,2]
non_norm[,2]<-D_table[,2]
non_norm[,3]<-G_table[,2]

head(non_norm, 20)

cor.test(non_norm$B, non_norm$G, method = "kendall")

res_non <- cor(non_norm)
round(res_non, 2)

## TTR normalised
  
B_ttr_table<-read.csv("B_R1_ttr.wig", sep = " ", skip=2)   # positions start on 3rd line
D_ttr_table<-read.csv("D_R1_ttr.wig", sep = " ", skip=2)
G_ttr_table<-read.csv("G_R1_ttr.wig", sep = " ", skip=2)


pos_vector<-B_ttr_table[,1]

ttr_norm<-as.data.frame(matrix(0, nrow = nrow(B_ttr_table), ncol =3), row.names <-pos_vector)
colnames(ttr_norm)<-c('B', 'D', 'G')
ttr_norm[,1]<-B_ttr_table[,2]
ttr_norm[,2]<-D_ttr_table[,2]
ttr_norm[,3]<-G_ttr_table[,2]

head(ttr_norm, 20)

res_ttr <- cor(ttr_norm)
round(res_ttr, 2)


tn_stats_df<-read.csv("bovis_combined_TTR.wig", header=TRUE, sep = "\t", skip=6, row.names=NULL, col.names = c('position', 'B', 'D', 'G', 'gene'))
View(tn_stats_df[1:50,])

# df for sum of reads for each ta site
comb_names<-c("position", "sum reads", "gene")
combined_df<-as.data.frame(matrix(0, nrow = nrow(tn_stats_df), ncol = 3))
combined_df[,1]<-tn_stats_df[,1]
for (i in 1:nrow(combined_df)){
  combined_df[i,2]<-sum(tn_stats_df[i,2:4])
}
combined_df[,3]<-tn_stats_df[,5]
colnames(combined_df)<-comb_names
View(combined_df[0:200,])

# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_count <- combined_df[,2]
ins_count<-0
for (i in 1:length(tot_count)){
  if (tot_count[i] != 0.0){
    ins_count <- ins_count + 1
  }
}
ins_count
# 34178
length(tot_count)
# 73536
insertion_density <- ins_count/length(tot_count)
insertion_density
# 0.4647792

# max read count
max(tot_count)
# 100024.5  (uses normalisation factor?)

# mean of no insertions per gene
mean(tot_count)
# 489.7486


#non zero mean (mean of read counts at all non-zero sites)
# change all zeros to na
is.na(tot_count) <- tot_count==0.0
tot_count[1:100]
mean(tot_count, na.rm=TRUE) 
# 1053.723


# do this again with non_normal

tn_non_stats_df<-read.csv("bovis_combined_nonorm.wig", header=TRUE, sep = "\t", skip=6, row.names=NULL, col.names = c('position', 'B', 'D', 'G', 'gene'))
View(tn_non_stats_df[1:20,])

# df for sum of reads for each ta site
comb_names<-c("position", "sum reads", "gene")
combined_non_df<-as.data.frame(matrix(0, nrow = nrow(tn_non_stats_df), ncol = 3))
combined_non_df[,1]<-tn_non_stats_df[,1]
for (i in 1:nrow(combined_non_df)){
  combined_non_df[i,2]<-sum(tn_non_stats_df[i,2:4])
}
combined_non_df[,3]<-tn_non_stats_df[,5]
colnames(combined_non_df)<-comb_names
View(combined_non_df[0:200,])

# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_non_count <- combined_non_df[,2]
ins_non_count<-0
for (i in 1:length(tot_non_count)){
  if (tot_non_count[i] != 0.0){
    ins_non_count <- ins_non_count + 1
  }
}

length(tot_non_count)
# 73535   (why is this one less than ttr?)
ins_non_count
# 34178  (exactly same as normalised with TTR?)
# same number of sites inserted but number of counts different
max(tot_non_count)
# 5913

mean(tot_non_count)
# 29.54797

#non zero mean 
is.na(tot_non_count) <- tot_non_count==0.0
mean(tot_non_count, na.rm=TRUE) 
# 63.57335


# try with one of the dejesus libraries
dejesus_40_df<-read.csv("SRR4113440_1.wig", sep = " ", skip=2, header=TRUE, row.names=NULL, col.names=c('position', 'count'))
View(dejesus_40_df[1:50,])

# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_40 <- dejesus_40_df[,2]
ins_count_40<-0
for (i in 1:length(tot_40)){
  if (tot_40[i] != 0.0){
    ins_count_40 <- ins_count_40 + 1
  }
}

length(tot_40)
# 74603
ins_count_40
# 39247

#insertion density
id_40<-ins_count_40/length(tot_40)
id_40
# 0.526078

max(tot_40)
# 4121

mean(tot_40)
# 29.79875

#non zero mean 
is.na(tot_40) <- tot_40==0.0
mean(tot_40, na.rm=TRUE) 
# 56.64321


# use same methods but get data from bwa mapping (dong's script)
# count number of insertion sites from each sample ('Total No. IS')
# add together
# summary_stats just for last sample. need to read from csv?
tot_bovis_sites<-sum(summary_stats$`No.Insertion sites`)
tot_bovis_sites
