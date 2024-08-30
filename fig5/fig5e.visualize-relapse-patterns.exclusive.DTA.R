#!/usr/bin/Rscript
# plot breakdown of relapse phylogenetic patterns
library(ggpubr)

# relapse pattern categorization was done manually, looking at the "AML_data.muts.chroms.txt" file

data=read.table("relapse-pattern-categories.txt",sep="\t",header=TRUE,row.names=1)
data$dtbin = paste(data$DNMT3A_Diagnosis,data$TET2_Diagnosis,data$ASXL1_Diagnosis, sep=":")

categories = as.data.frame.matrix(t(cbind(table(data$category),table(data$category)/nrow(data))))

bincats = as.data.frame.matrix(table(data$dtbin,data$category))
bincatspct = bincats/rowSums(bincats)*100

bincatspct[8,] = categories[2,]*100
rownames(bincatspct) = c("WT","ASXL1mut","TET2mut","TET2.ASXL1.mut","DNMT3Amut","DNMT3A.ASXL1.mut","DNMT3A.TET2.mut","all")

library(reshape2)
bincatspct$label = rownames(bincatspct)
meltdata = melt(bincatspct)

write.table(bincatspct,"relapse-categories.DNMT3A-TET2-ASXL1.bin.txt",sep="\t",row.names=FALSE)
write.table(bincats,"relapse-categories.DNMT3A-TET2-ASXL1.bin.cts.txt",sep="\t",row.names=TRUE)


library(ghibli)
meltdata = meltdata[meltdata$label %in% c("WT","DNMT3Amut","TET2mut","ASXL1mut"),]
meltdata$label = factor(meltdata$label,levels=c("WT","DNMT3Amut","TET2mut","ASXL1mut"))

meltdata$variable = factor(meltdata$variable,levels=c("relapse of same clone with no new mutations", "relapse of same clone with loss of mutations","relapse of same clone with new mutations","subclonal swap, both gain and loss of mutations"))

levels(meltdata$variable) = c("stable mutations","mutation loss","mutation gain","subclonal swap")

meltdata$variable = factor(meltdata$variable, levels=levels(meltdata$variable))


pdf("phylo-categories-DNMT3A-vs-TET2-vs-ASXL1.exclusive.barplot.pdf",height=4,width=9,useDingbats=FALSE)
ggbarplot(meltdata,x="label",y="value",group="variable",fill="variable", color=NA,position = position_dodge(.8),palette=ghibli_palettes$LaputaMedium[c(3,5,6,7)])
dev.off()
