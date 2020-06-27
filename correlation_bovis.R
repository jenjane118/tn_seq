

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
