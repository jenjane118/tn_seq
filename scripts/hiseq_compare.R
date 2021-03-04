# sum up results file for Transit hmm for new bovis hiseq set

library(dplyr)

hmm_results<-read.delim("hmm_bovis_hiseq_combo_genes.wig", 
                        sep="\t", header=FALSE, stringsAsFactors=F, 
                        comment.char = '#')
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
bovis_hmm<-unique(bovis_hmm)
nrow(bovis_hmm)
View(bovis_hmm)

essential_genes_hmm <- bovis_hmm[bovis_hmm$'call' == 'ES',1]
sort(essential_genes_hmm)
length(essential_genes_hmm)

gd_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GD',1])
length(gd_genes)
ga_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'GA',1])
length(ga_genes)
non_ess_genes <- unique(bovis_hmm[bovis_hmm$'call' == 'NE',1])
length(non_ess_genes)

# sum up results file for Transit hmm for new mtb hiseq set

hmm_results<-read.delim("mtb_hiseq_combined_genes.wig", 
                        sep="\t", header=FALSE, stringsAsFactors=F, 
                        comment.char = '#')
head(hmm_results)


# change to just gene and call--will affect downstream plots
mtb_hmm <- select(hmm_results, 2, 1, 11)
colnames(mtb_hmm) <- c("gene","ORF", "call")

nrow(mtb_hmm)
View(mtb_hmm)

#write to .csv for sharon
write.csv(mtb_hmm, "hiseq_mtb_transit.csv", quote=F)


install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#load up wig file insertions for bovis
b19<-read.delim("bovis_hiseq/tpp/bovis_hiseq_tpp_19.wig", header=T, comment.char = "#", sep=" ")
b20<-read.delim("bovis_hiseq/tpp/bovis_hiseq_tpp_20.wig", header=T, comment.char = "#", sep=" ")
b21<-read.delim("bovis_hiseq/tpp/bovis_hiseq_tpp_21.wig", header=T, comment.char = "#", sep=" ")

insertions<-b19[,2]
insertions_df<-cbind(insertions, b20[,2], b21[,2])
head(insertions_df)

totals<-rowSums(insertions_df[, c(1, 2, 3)])
head(totals)

# sum insertions into one df
positions<-b19[,1]
mb_ins_df<-data.frame(positions, totals)
#mb_ins_df<-cbind(positions, totals)
colnames(mb_ins_df)<-c("positions", "insertions")
head(mb_ins_df)
nrow(mb_ins_df)


max(mb_ins_df$insertions)
mb_ins_df[which.max(mb_ins_df[,2])]


plot(mb_ins_df$positions, mb_ins_df$insertions, type='l', main="insertion plot mbovis",
              xlab="TA position", ylab="number of insertions")

#ggplot(data = mb_ins_df) +
  #geom_bar(mapping = aes(x = positions, y = insertions), stat = "identity")

## this works better for getting pretty axes
ggplot(data=mb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  scale_x_continuous(name = positions, 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver())

## make bigger plot for export
## to make bigger plot for export
png(file="mbo_insertions.png", width=2048, height=1536)
ggplot(data=mb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title=element_text(size=16,face="bold")) +
  ylim(0, 10000) +
  scale_x_continuous(name = "position", 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver()) +
  ggtitle("Insertion plot M.bovis (hiseq)") + 
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20, hjust=0.5))
dev.off()


# do the same with tb
# load .wig files
b22<-read.delim("mtb_hiseq/TPP/mtb22_hiseq_tpp.wig", header=T, comment.char = "#", sep=" ")
b23<-read.delim("mtb_hiseq/TPP/mtb23_hiseq_tpp.wig", header=T, comment.char = "#", sep=" ")

ins_df<-data.frame(b22[,2], b23[,2])
totals<-rowSums(ins_df[, c(1,2)])
positions<-b22[,1]
mtb_ins_df<-data.frame(positions, totals)
colnames(mtb_ins_df)<-c("positions", "insertions")
head(mtb_ins_df)

ggplot(data=mtb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

## to make bigger plot for export
png(file="mtb_insertions.png", width=2048, height=1536)
ggplot(data=mtb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title=element_text(size=16,face="bold")) +
  ylim(0, 10000) +
  scale_x_continuous(name = "position", 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver()) +
  ggtitle("Insertion plot M.tb (hiseq)") + 
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20, hjust=0.5))
dev.off()
