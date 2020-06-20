GFF<-read.table("short_bovis.gff",header = F,skip = 2,sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
#GFF<-GFF[1:(grep("##FASTA",GFF[,1])-1),]


features_table<-data.frame(matrix(NA,ncol=5,nrow=nrow(GFF)))
colnames(features_table)<-c("START","END","NAME","TYPE","PRODUCT")


for (i in 1:nrow(features_table)){
  if (GFF[i,5]<GFF[i,4]){
    features_table$START[i]<-GFF[i,5]
    features_table$END[i]<-GFF[i,4]
  } else {
    features_table$START[i]<-GFF[i,4]
    features_table$END[i]<-GFF[i,5]
  }
  features_table$TYPE[i]<-GFF[i,3]
  a<-unlist(strsplit(GFF[i,9],split = ";"))
  b<-a[grep(pattern = "gene=",x = a)]
  c<-a[grep(pattern = "product=",x = a)]
  d<-a[grep(pattern = "locus_tag=",x = a)]
  if (length(b)>0){
    features_table$NAME[i]<-substr(b,6,nchar(b))
  } else if (length(d)>0){
    features_table$NAME[i]<-substr(d,11,nchar(d))
  } else {
    features_table$NAME[i]<-a[1]
  }
  if (length(c)>0){
    features_table$PRODUCT[i]<-substr(c,9,nchar(c))
  }
}

View(features_table)

write.csv(features_table,"M_bovis_features.csv",row.names = F)
