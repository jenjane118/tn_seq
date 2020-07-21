# looking at transit analysis results

transit_data <- read.csv("bovis_TTR.txt", sep="\t", skip=11, header= FALSE)
head(transit_data)
colnames(transit_data)<-c("ORF","Name", "Descr", "no ins in ORF", "total TA in ORF", "len max run non-ins", "nt span for max run non-ins", "zbar", "call")
