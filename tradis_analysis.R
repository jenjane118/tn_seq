# Viewing tradis_summary.csv report for Biotradis analysis of mbovis data
# generated from Dong Xia's script: TraDIS_summary_Mbovis-LT708304.1.R


# read summary report for each sample (created per sample by script)
tradis_summary_B<-read.csv("B_S1_L001_R1_001_TraDIS_summary.csv",header=T)
colnames(tradis_summary_B)
View(tradis_summary_B)[0:20]

