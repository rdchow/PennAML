#!/usr/bin/Rscript
# plot breakdown of relapse phylogenetic patterns
library(ggpubr)
library(reshape2)
library(ghibli)

# relapse pattern categorization was done manually, looking at the "AML_data.muts.chroms.txt" file

data=read.table("refractory-pattern-categories.txt",sep="\t",header=TRUE,row.names=1)
categories = as.data.frame.matrix(t(cbind(table(data$category),table(data$category)/nrow(data)*100)))
meltdata = melt(categories[2,])
meltdata$variable = factor(meltdata$variable,levels=(c("relapse of same clone with no new mutations", "relapse of same clone with loss of mutations","relapse of same clone with new mutations","subclonal swap, both gain and loss of mutations")))
levels(meltdata$variable) = c("stable mutations","mutation loss","mutation gain","subclonal swap")
meltdata$variable = factor(meltdata$variable, levels=rev(levels(meltdata$variable)))

g1 = ggbarplot(meltdata,x="variable",y="value",group="variable",fill="variable",color=NA,orientation = "horiz",palette=rev(ghibli_palettes$LaputaMedium[c(3,5,6,7)])) 

data2=read.table("relapse-pattern-categories.txt",sep="\t",header=TRUE,row.names=1)
categories2 = as.data.frame.matrix(t(cbind(table(data2$category),table(data2$category)/nrow(data2)*100)))
meltdata2 = melt(categories2[2,])
meltdata2$variable = factor(meltdata2$variable,levels=(c("relapse of same clone with no new mutations", "relapse of same clone with loss of mutations","relapse of same clone with new mutations","subclonal swap, both gain and loss of mutations")))
levels(meltdata2$variable) = c("stable mutations","mutation loss","mutation gain","subclonal swap")
meltdata2$variable = factor(meltdata2$variable, levels=rev(levels(meltdata2$variable)))

g2 = ggbarplot(meltdata2,x="variable",y="value",group="variable",fill="variable",color=NA,orientation = "horiz",palette=rev(ghibli_palettes$LaputaMedium[c(3,5,6,7)])) 

pdf("phylo-categories-allSamples.barplot.REL1.REF1.pdf",height=3.5,width=7.5,useDingbats=FALSE)
ggarrange(g2,g1,ncol=2)
dev.off()