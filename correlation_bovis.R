sample_names<-read.table("names_wig.txt",stringsAsFactors = F)

#for (j in 1:nrow(sample_names)){
  #insertion_table<-read.csv(paste(sample_names[j,1],".wig",sep = ""))
  #non_norm<-as.data.frame(matrix(0, nrow = nrow(insertion_table), ncol = 3))
  


B_table<-read.csv("B_R1.wig", sep = " ", skip=2)   # positions start on 3rd line
D_table<-read.csv("D_R1.wig", sep = " ", skip=2)
G_table<-read.csv("G_R1.wig", sep = " ", skip=2)

pos_vector<-B_table[,1]
pos_vector
D_vector<-D_table[,1]



non_norm<-as.data.frame(matrix(0, nrow = nrow(B_table), ncol =3), row.names <-pos_vector)
colnames(non_norm)<-c('B', 'D', 'G')
non_norm[,1]<-B_table[,2]
non_norm[,2]<-D_table[,2]
non_norm[,3]<-G_table[,2]

View(non_norm)

cor.test(non_norm$B, non_norm$G, method = "kendall")
  
  
  