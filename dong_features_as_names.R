#GFF<-read.table("Genbank-Mbovis-sequence.gff3",header = F,skip = 2,sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
GFF<-read.table("LT708304_updated_aug19.gff",header = F,skip = 2,sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
View(GFF)
#GFF<-GFF[1:(grep("##FASTA",GFF[,1])-1),]
# add row for gene 'name' (i.e. 'dnaA' and aa length)

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
  
  if (length(d)>0){
    features_table$NAME[i]<-substr(d,18,nchar(d))
  }
  
  else if (length(b)>0){
    features_table$NAME[i]<-substr(b,6,nchar(b))
    
  } else {
    features_table$NAME[i]<-a[1]
  }
  if (length(c)>0){
    features_table$PRODUCT[i]<-substr(c,9,nchar(c))
  }
  
  
}

View(features_table)
#write.csv(features_table,"M_bovis_features-locusAsName.csv",row.names = F)
write.csv(features_table,"M_bovis_features-locusAsName-aug19.csv",row.names = FALSE)
