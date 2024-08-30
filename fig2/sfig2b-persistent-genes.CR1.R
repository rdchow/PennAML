#!/usr/bin/Rscript
# barplot showing gain/loss of mutations at CR1 by gene cateogry

library(ggplot2)
library(ggpubr)
library(rcartocolor)
library(reshape2)


data=read.table("../fig2/pt-level-mut-freqs.txt",sep="\t",header=TRUE)
data = data[data$gene != "TPMT",]
subdata = data[!is.na(data$CR1) & (data$diagnosis >0 | data$CR1 > 0),]  # select all variants in patients that had CR1 sequenced

subdata$CR1Bin = subdata$CR1
subdata[subdata$CR1Bin>0,"CR1Bin"] = 1
subdata$DiagBin = subdata$diagnosis
subdata[subdata$DiagBin>0,"DiagBin"] = 1
subdata$comboBin = paste(subdata$CR1Bin,"_", subdata$DiagBin)

mydf = as.data.frame.matrix(table(subdata$gene,subdata$comboBin))
colnames(mydf) = c("Diag_only","CR1_only","shared")
mydf = mydf[(mydf$CR1_only + mydf$Diag_only + mydf$shared) >= 2,]

mydf3 = mydf
mydf3$cr1Onlypct = mydf3$CR1_only/(mydf3$Diag_only + mydf3$CR1_only + mydf3$shared)*100
mydf3$diagOnlypct = mydf3$Diag_only/(mydf3$Diag_only + mydf3$CR1_only + mydf3$shared)*100
mydf3$sharedpct = mydf3$shared/(mydf3$Diag_only + mydf3$CR1_only + mydf3$shared)*100
mydf3 = mydf3[order(mydf3$sharedpct,decreasing=TRUE),]

mydf3$gene = rownames(mydf3)
mydf3$gene = factor(mydf3$gene,levels=rev(mydf3$gene))
mydf3$fraction = paste(mydf3$shared,mydf3$Diag_only,mydf3$CR1_only,sep=";")

meltdata = melt(mydf3,id.vars=c("gene","fraction"),measure.vars=c("sharedpct","diagOnlypct","cr1Onlypct"))
meltdata$variable = factor(meltdata$variable,levels=rev(c("sharedpct","diagOnlypct","cr1Onlypct")))
meltdata$gene = factor(meltdata$gene,levels=rev(unique(meltdata$gene)))

#pdf("diagnosis-CR1.barplot.pdf",height=6,width=6,useDingbats=FALSE)
g1 = ggbarplot(meltdata,x="gene",y="value",orientation="horizontal",fill="variable",group="variable",color=NA,label=mydf3$fraction,palette=rev(c("#BDBDBD","#EDAD08","#48AEAE")))+scale_y_continuous(expand=c(0,0),limits=c(0,101)) + rremove("legend")

# then make the violin/boxplots of the delta VAFS (CR1 - diag) 
subdata$cr_diag = subdata$CR1 - subdata$diagnosis
subdata2 = subdata[subdata$gene %in% unique(meltdata$gene),]
subdata2$gene = factor(subdata2$gene,levels=levels(meltdata$gene))

g2 = ggviolin(subdata2,x="gene",y="cr_diag",orientation="horizontal",fill="gene",group="gene",scale="width") + scale_y_continuous(expand=c(0,0),limits=c(-100,100)) + rremove("legend") + geom_hline(yintercept = 0, linetype = 2, color="gray44")

library(colorspace)
mycols = carto_pal(36,"Temps")
darkcols = darken(mycols,0.45)

g2 = ggplot(subdata2,aes(x=gene,y=cr_diag)) + geom_hline(yintercept=c(0),linetype="dashed",color="gray44") + geom_violin(aes(fill=gene),trim=TRUE,scale="width",color=NA,width=0.8) + geom_jitter(aes(color=gene),width=0.2,size=2) +geom_boxplot(fill=NA,width=0.2,outlier.shape=NA)+ scale_y_continuous(expand=c(0,0),limits=c(-101,101)) + scale_fill_manual(values=mycols) + scale_color_manual(values=darkcols) + theme_pubr() + coord_flip()  + rremove("legend")

ggarrange(g1,g2,ncol=2,align="v")

pdf("barplot-diag-vs-CR1-muts.pdf",height=10,width=13,useDingbats=FALSE)
ggarrange(g1,g2,ncol=2,align="v",widths=c(0.7,1))
dev.off()