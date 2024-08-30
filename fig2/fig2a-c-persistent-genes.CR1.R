#!/usr/bin/Rscript
# barplot showing persistence of mutations at CR1 by gene cateogry


data=read.table("pt-level-mut-freqs.txt",sep="\t",header=TRUE)
data = data[data$gene != "TPMT",]
subdata = data[!is.na(data$CR1) & (data$diagnosis >0),]  # select variarnts present at diagnosis, in patients that were also sequenced at CR1
#subdata = data[!is.na(data$CR1),]  # select all patients/variants that had CR1 sequenced

subdata$CR1Bin = subdata$CR1
subdata[subdata$CR1Bin>0,"CR1Bin"] = 1

mydf = as.data.frame.matrix(table(subdata$gene,subdata$CR1Bin))
colnames(mydf) = c("wt","mut")

dna = c("ATM","TP53")
dtai = c("DNMT3A","TET2","ASXL1","IDH1","IDH2")
splice = c("SRSF2","U2AF1","ZRSR2")
prcrunx = c("BCOR","BCORL1","RUNX1","EZH2")
cohesin = c("SMC1A","RAD21","STAG2")
signaling = c("CSF1R","FLT3","NF1","KRAS","NRAS","BRAF","KIT","PTPN11","JAK2","CSF3R","CBL")
npm1 = "NPM1"


mydf2 = mydf[rownames(mydf) %in% c(dna,dtai,splice,prcrunx,cohesin,signaling,npm1),]
mydf2$category = ""
mydf2[rownames(mydf2) %in% dna,"category"] = "DNA damage"
mydf2[rownames(mydf2) %in% dtai,"category"] = "DTAI"
mydf2[rownames(mydf2) %in% splice,"category"] = "Splicing"
mydf2[rownames(mydf2) %in% prcrunx,"category"] = "PRC/RUNX"
mydf2[rownames(mydf2) %in% cohesin,"category"] = "Cohesin"
mydf2[rownames(mydf2) %in% signaling,"category"] = "Signaling"
mydf2[rownames(mydf2) %in% npm1,"category"] = "NPM1"


## aggregated by gene categories
aggdata = aggregate( .~category,data=mydf2,sum)

aggdata$mutpct = aggdata$mut/(aggdata$mut+aggdata$wt)*100
aggdata = aggdata[order(aggdata$mutpct,decreasing=FALSE),]

library(reshape2)
meltdf = melt(aggdata,id.vars=c("category","wt","mut"),measure.vars="mutpct")
meltdf$fraction = paste(meltdf$mut,"/",meltdf$mut+meltdf$wt,sep="")

library(ggplot2)
library(ggpubr)
library(rcartocolor)
meltdf2 = meltdf[meltdf$variable == "mutpct",]

g1= ggbarplot(meltdf2,x="category",y="value",orientation="horizontal",fill="category",color=NA,label=meltdf2$fraction)+ geom_hline(yintercept=c(50),linetype="dashed",color="gray44")+scale_y_continuous(expand=c(0,0),limits=c(0,100))+ scale_fill_manual(values=carto_pal(7,"TealGrn"))+rremove("legend")

### split by individual genes
# None retained for NRAS; removed
mydf3 = mydf2[rownames(mydf2) %in% c("DNMT3A","TET2","ASXL1","IDH2","IDH1","SRSF2","RUNX1","FLT3","NPM1","NRAS","TP53"),]
mydf3$mutpct = mydf3$mut/(mydf3$wt+mydf3$mut)*100
mydf3$gene = rownames(mydf3)
mydf3$gene = factor(mydf3$gene,levels=rev(c("DNMT3A","TET2","ASXL1","IDH2","IDH1","SRSF2","RUNX1","FLT3","NPM1","NRAS","TP53")))
mydf3$fraction = paste(mydf3$mut,"/",mydf3$mut+mydf3$wt,sep="")

g2= ggbarplot(mydf3,x="gene",y="mutpct",orientation="horizontal",fill="gene",color=NA,label=mydf3$fraction,palette=carto_pal(11,"SunsetDark"))+ geom_hline(yintercept=c(50),linetype="dashed",color="gray44")+scale_y_continuous(expand=c(0,0),limits=c(0,100)) + rremove("legend")

# then make the violin with the actual VAFs @CR1
subdata2 = subdata[subdata$gene %in% c("DNMT3A","TET2","ASXL1","IDH2","IDH1","SRSF2","RUNX1","FLT3","NPM1","NRAS","TP53"),]
subdata3 = subdata2[subdata2$CR1 >0,]
#subdata3 = subdata2

# add dummy data for NRAS so that it gets plotted
subdata3 = as.data.frame(rbind(subdata3,c("x","NRAS","x",0,0,rep("NA",6)), c("x","NRAS","x",0,0,rep("NA",6))))

library(colorspace)
mycols = carto_pal(11,"SunsetDark")
darkcols = darken(mycols,0.25)

subdata3$gene = factor(subdata3$gene,levels=rev(c("DNMT3A","TET2","ASXL1","IDH2","IDH1","SRSF2","RUNX1","FLT3","NPM1","NRAS","TP53")))
subdata3$CR1 = as.numeric(subdata3$CR1)

g3 = ggplot(subdata3,aes(x=gene,y=CR1)) + geom_hline(yintercept=c(50),linetype="dashed",color="gray44")+ geom_violin(aes(fill=gene),trim=TRUE,scale="width",color=NA,width=0.8) + geom_jitter(aes(color=gene),width=0.2,size=2.5) + scale_y_continuous(expand=c(0,0),limits=c(-1,100)) + scale_fill_manual(values=mycols) + scale_color_manual(values=darkcols) + theme_pubr() + coord_flip()  + rremove("legend")


pdf("barplot-pct-retained-CR1-muts.categorized.pdf",height=6,width=15,useDingbats=FALSE)
ggarrange(g1, g2,g3,ncol=3,align="v")
dev.off()


library(plotrix)
aggregate(subdata3, CR1 ~ gene,FUN = function(x) c(mean = mean(x), se = std.error(x)))
