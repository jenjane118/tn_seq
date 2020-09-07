## compare essential/non-essential calls between tradis and Sharon's lab

library(dplyr)
# results file for Sharon's group
custom_results<-read.table("results/BCD_bothBatch_TraDIS_summary.csv", sep=",",header=TRUE,stringsAsFactors=F, quote="\"")
head(custom_results)

# values for gamma fits of insertion index column in results file (custom_results)
ess_change<-0.0045
amb_change<-0.0072

bovis_custom <- select(custom_results, 1,3,8)
colnames(bovis_custom)<-c("Name", "Product", "Ins index")
bovis_custom$Name<-toupper(bovis_custom$Name)
head(bovis_custom)

for (i in 1:nrow(bovis_custom)){
  if (bovis_custom[i,3] <= ess_change){
      bovis_custom$Call[i] <- 'ES'
  }
  else if (bovis_custom[i,3] < amb_change & bovis_custom[i,3] > ess_change){
      bovis_custom$Call[i] <- 'AMB'
  }
  else{bovis_custom$Call[i] <-'NE'
  }
}

View(bovis_custom)
length(bovis_custom$Name)
essential_genes_custom <- bovis_custom[bovis_custom$Call == 'ES',1]

length(essential_genes_custom)
sort(essential_genes_custom)

# results file for Transit hmm

hmm_results<-read.delim("results/bovis_hmm_genes.txt", sep="\t", header=FALSE, stringsAsFactors=F, skip=6)
head(hmm_results)


# change to just gene and call--will affect downstream plots
bovis_hmm <- select(hmm_results, 2, 1, 11)
colnames(bovis_hmm) <- c("gene","ORF", "call")
# have to trim whitespace before gene to compare with mtb
bovis_hmm$gene<-trimws(bovis_hmm$gene)
View(bovis_hmm)
length(bovis_hmm$call)

#this calls two for each gene, have to ask for 'unique'
length(unique(bovis_hmm$ORF))

essential_genes_hmm <- unique(bovis_hmm[bovis_hmm$'call' == 'ES',1])
sort(essential_genes_hmm)
length(essential_genes_hmm)

gd_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GD',1])
length(gd_genes)
ga_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GA',1])
length(ga_genes)
non_ess_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'NE',1])
length(non_ess_genes)

# biotradis results

biotradis_ess_change <- 0.0018
biotradis_amb_change <- 0.0028


biotradis_res<-read.delim("results/joined_output.m_bovis_BDG.csv", sep="\t",header=TRUE,stringsAsFactors=F,quote = "")
head(biotradis_res)

biotradis <- select(biotradis_res, 1, 2, 8)
colnames(biotradis)<-c("Name", "Product", "Ins index")
for (i in 1:nrow(biotradis)){
  biotradis$Name[i] <- substr(biotradis[i,1], 8, nchar(biotradis[i,1]))
}
biotradis$Name<-toupper(biotradis$Name)
head(biotradis)
View(biotradis)

for (i in 1:nrow(biotradis)){
  if (biotradis[i,3] <= biotradis_ess_change){
    biotradis$Call[i] <- 'ES'
  }
  else if (biotradis[i,3] < biotradis_amb_change & biotradis[i,3] > biotradis_ess_change){
    biotradis$Call[i] <- 'AMB'
  }
  else{biotradis$Call[i] <-'NE'
  }
}

## compare with transit gumbel algorithm
gumbel_results<-read.delim("results/bovis_nonorm_edit.txt", sep="\t", header=FALSE, stringsAsFactors=F, skip=12)
head(gumbel_results)

bovis_gumbel <- select(gumbel_results, 1, 2, 9)
colnames(bovis_gumbel) <- c("ORF","gene","call")
bovis_gumbel$ORF<-toupper(bovis_gumbel$ORF)
View(bovis_gumbel)
length(bovis_gumbel$call)

#this calls two for each gene, have to ask for 'unique'
length(unique(bovis_gumbel$ORF))
#4101

essential_genes_gumbel <- unique(bovis_gumbel[bovis_gumbel$'call' == 'E',1])
sort(essential_genes_gumbel)
length(essential_genes_gumbel)
#change 'E' to 'ES'
#bovis_gumbel$call <- lapply(bovis_gumbel$call, gsub, pattern = "E", replacement = "ES", fixed = TRUE)
#this changes NE to NES

for (i in 1:nrow(bovis_gumbel)){
  if (bovis_gumbel$call[i] == 'E'){
    bovis_gumbel$call[i] = 'ES'
  }
}
View(bovis_gumbel)


uncertain_genes <- unique(bovis_gumbel[bovis_gumbel$'call' == 'U',1])
length(uncertain_genes)
short_genes <- unique(bovis_gumbel[bovis_gumbel$'call' == 'S',1])
length(short_genes)
non_ess_gumbel <- unique(bovis_gumbel[bovis_gumbel$'call' == 'NE',1])
length(non_ess_gumbel)


#####



essential_biotradis <-biotradis[biotradis$Call == 'ES',1]
length(essential_biotradis)
#531
length(essential_genes_hmm)
#499
length(essential_genes_custom)
#664
length(essential_genes_gumbel)
#362

common_ess<-NULL
for (i in 1:length(essential_genes_custom)){
  if (essential_genes_gumbel[i] %in% essential_genes_hmm){
    common_ess<-c(common_ess,essential_genes_gumbel[i])
  }
}
common_ess
length(common_ess)

bio_hmm_common <- NULL
for (i in 1:length(essential_biotradis)){
  if (essential_biotradis[i] %in% essential_genes_hmm){
    bio_hmm_common<-c(bio_hmm_common, essential_biotradis[i])
  }
}
length(bio_hmm_common)

# need to use annotation file and make chart, then fill in what each call is

pos_vector<-bovis_custom[,1]
pos_vector

compare_df<-as.data.frame(matrix(0, nrow = nrow(bovis_custom), ncol = 6))#, row.names <-pos_vector)
colnames(compare_df)<-c('ORF', 'gene', 'hmm_call', 'custom_call', 'biotradis_call', 'gumbel_call')
compare_df$ORF<-bovis_custom$Name

head(compare_df)

compare_df$custom_call<-with(bovis_custom, Call[match(compare_df$ORF, Name)])
head(compare_df)

compare_df$hmm_call<-with(bovis_hmm, call[match(compare_df$ORF, ORF)])
head(compare_df)

compare_df$biotradis_call<-with(biotradis, Call[match(compare_df$ORF, Name)])
head(compare_df)

compare_df$gumbel_call<-with(bovis_gumbel, call[match(compare_df$ORF, ORF)])
head(compare_df)

compare_df$gene<-with(bovis_hmm, gene[match(compare_df$ORF, ORF)])

View(compare_df)

write.table(essential_biotradis, file='results/biotradis_ess_list', quote=FALSE, sep = '', row.names=FALSE)
write.table(essential_genes_custom, file='results/custom_ess_list', quote=FALSE, sep = '', row.names=FALSE)
write.table(essential_genes_hmm, file='results/hmm_ess_list', quote=FALSE, sep = '\t', row.names=FALSE)
write.table(essential_genes_gumbel, file='results/gumbel_ess_list', quote=FALSE, sep='\t', row.names=FALSE)
write.table(compare_df, file = 'results/tnseq_compare.txt', quote = FALSE, sep = '\t', row.names=FALSE)

# install.packages('VennDiagram')
library(VennDiagram)

# compare essential genes with hmm results from mtb

mtb_genes<-read.delim("results/dejesus-Mtb_genes.txt", sep="\t", header=FALSE, stringsAsFactors=F, skip=4)
head(mtb_genes)

mtb_hmm <- select(mtb_genes, 2, 1, 11)
colnames(mtb_hmm) <- c("gene","ORF","call")
mtb_hmm$ORF<-toupper(mtb_hmm$ORF)
View(mtb_hmm)
length(mtb_hmm$call)

# change to just gene and call--will affect downstream plots
bovis_hmm <- select(hmm_results, 2, 1, 11)
colnames(bovis_hmm) <- c("gene","ORF", "call")
# have to trim whitespace before gene to compare with mtb
bovis_hmm$gene<-trimws(bovis_hmm$gene)
View(bovis_hmm)
length(bovis_hmm$call)


#this calls two for each gene, have to ask for 'unique'
# get rid of extra lines 
newbovis_hmm<-unique(bovis_hmm)
View(newbovis_hmm)
length(newbovis_hmm$gene)

# how many genes in common between mtb and mbovis using hmm analysis

# prints list of genes in both lists
newbovis_hmm$gene[newbovis_hmm$gene %in% mtb_hmm$gene]
common_genes <- newbovis_hmm$gene[newbovis_hmm$gene %in% mtb_hmm$gene]
common_genes <- sort(common_genes)
length(common_genes)

bovis_mtb <- as.data.frame(matrix(0, nrow = length(common_genes), ncol = 3))
colnames(bovis_mtb)<-c('gene', 'bovis_call', 'mtb_call')

for (i in 1:length(common_genes)){
  bovis_mtb$gene[i] <- common_genes[i]
  bovis_mtb$bovis_call[i] <- newbovis_hmm[which(newbovis_hmm$gene == common_genes[i]),3]
  bovis_mtb$mtb_call[i]   <- mtb_hmm[which(mtb_hmm$gene == common_genes[i]),3]
}

View(bovis_mtb)

write.table(bovis_mtb, file = 'results/comp_mtb_bovis.txt', quote = FALSE, sep = '\t', row.names=FALSE)

#number of genes same call
same_call <- NULL
for (i in 1:nrow(bovis_mtb)){
  if (bovis_mtb$bovis_call[i]==bovis_mtb$mtb_call[i]){
    same_call <- c(same_call, bovis_mtb$gene[i])
  }
}
length(bovis_mtb$gene)
#[1] 1484
length(same_call)
#[1] 980

#essential/ga in mtb but non-essential in mbovis
essentials<-c("ES", "GA")
non_ess<-c("NE", "GD")
ess_mtb <- NULL
for (i in 1:nrow(bovis_mtb)){
  if (bovis_mtb$mtb_call[i] %in% essentials & bovis_mtb$bovis_call[i] %in% non_ess){
    ess_mtb <- c(ess_mtb, bovis_mtb$gene[i])
  }
}
ess_mtb
length(ess_mtb)

#essential/ga in mtb called es/ga in mbovis
ess_both<-NULL
for (i in 1:nrow(bovis_mtb)){
  if (bovis_mtb$mtb_call[i] %in% essentials & bovis_mtb$bovis_call[i] %in% essentials){
    ess_both <- c(ess_both, bovis_mtb$gene[i])
  }
}
length(ess_both)
ess_both

# non-essential/gd in both
non_ess_both<-NULL
for (i in 1:nrow(bovis_mtb)){
  if (bovis_mtb$mtb_call[i] %in% non_ess & bovis_mtb$bovis_call[i] %in% non_ess){
    non_ess_both <- c(non_ess_both, bovis_mtb$gene[i])
  }
}
non_ess_both
length(non_ess_both)


essential_genes_hmm <- unique(bovis_hmm[bovis_hmm$'call' == 'ES',1])
sort(essential_genes_hmm)
length(essential_genes_hmm)

gd_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GD',1])
length(gd_genes)
ga_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GA',1])
length(ga_genes)
non_ess_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'NE',1])
length(non_ess_genes)

