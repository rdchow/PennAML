#!/usr/bin/Rscript
# barplot showing persistence of mutations at CR1 by gene cateogry


library(ggpubr)
library(rcartocolor)
data=read.table("pt-level-mut-freqs.txt",sep="\t",header=TRUE)
data = data[data$gene != "TPMT",]

##########
# plot the diagnosis vs CR1 comparisons, gene level
# flt3, npm1, nras
dcdata = data[!is.na(data$diagnosis) & !is.na(data$CR1),]
dcdata1 = dcdata[dcdata$gene %in% c("FLT3","NPM1","NRAS"),]

g1 = ggscatter(dcdata1,x="diagnosis",y="CR1",fill="gene",color="gene",title="Diag vs CR1",shape=21,size=4,alpha=0.8) + xlim(0,100)+ylim(0,100) + geom_abline(slope=1,intercept=0,linetype="dashed",color="gray44") + theme(aspect.ratio=1) +scale_fill_manual(values=c("#81C770","#C070C7","#D9845B")) + scale_color_manual(values=c("#608558","#682B6D","#994015"))


# asxl1, dnmt3a, tet2
dcdata2 = dcdata[dcdata$gene %in% c("ASXL1","DNMT3A","TET2"),]

g2 = ggscatter(dcdata2,x="diagnosis",y="CR1",fill="gene",color="gene",title="Diag vs CR1",shape=21,size=4,alpha=0.8) + xlim(0,100)+ylim(0,100) + geom_abline(slope=1,intercept=0,linetype="dashed",color="gray44") + theme(aspect.ratio=1) +scale_fill_manual(values=c("#F16EA3","#D9C75B","#37A1DA")) + scale_color_manual(values=c("#80274B","#BC921B","#0D6C9F"))

# IDH1, IDH2
dcdata3 = dcdata[dcdata$gene %in% c("IDH1","IDH2"),]
g3 = ggscatter(dcdata3,x="diagnosis",y="CR1",fill="gene",color="gene",title="Diag vs CR1",shape=21,size=4,alpha=0.8) + xlim(0,100)+ylim(0,100) + geom_abline(slope=1,intercept=0,linetype="dashed",color="gray44") + theme(aspect.ratio=1) +scale_fill_manual(values=c("#3CBD9F","#5847B0")) + scale_color_manual(values=c("#126B56","#4F1EAE"))

pdf("diag-vs-CR1.selectGene.scatter.pdf",height=4,width=12,useDingbats=FALSE)
ggarrange(g3,g2,g1,ncol=3)
dev.off()

##########
# plot the diagnosis vs CR1 comparisons, gene level
# SRSF2, RUNX1, TP53
dcdata = data[!is.na(data$diagnosis) & !is.na(data$CR1) & (data$diagnosis >0 | data$CR1 > 0),]
dcdata1 = dcdata[dcdata$gene %in% c("SRSF2","RUNX1","TP53"),]

g4 = ggscatter(dcdata1,x="diagnosis",y="CR1",fill="gene",color="gene",title="Diag vs CR1",shape=21,size=4,alpha=0.8) + xlim(0,100)+ylim(0,100) + geom_abline(slope=1,intercept=0,linetype="dashed",color="gray44") + theme(aspect.ratio=1) +scale_fill_manual(values=c("#81C770","#C070C7","#D9845B")) + scale_color_manual(values=c("#608558","#682B6D","#994015"))

pdf("diag-vs-CR1.selectGene.SRSF2.RUNX1.TP53.scatter.pdf",height=4,width=4,useDingbats=FALSE)
g4
dev.off()
