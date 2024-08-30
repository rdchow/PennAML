#!/usr/bin/Rscript
library(ggpubr)
setwd("C:/Bowman/AML-longitudinal/final-version/fig4//fig4f-MDS-IWG")
data=read.table("mutations.txt",sep="\t",header=TRUE)  # from cbioportal
subdat = data[,3:9]
subdat[subdat == "WT"] = 0
subdat[subdat!= 0] = 1

dataDN = as.data.frame(as.matrix(sapply(subdat, as.numeric))) 

# mutually exclusive
datad = dataDN[dataDN$DNMT3A == 1 & dataDN$TET2 == 0 & dataDN$ASXL1 == 0 ,]
datat = dataDN[dataDN$DNMT3A == 0 & dataDN$TET2 == 1 & dataDN$ASXL1 == 0,]
datab = dataDN[dataDN$DNMT3A == 0 & dataDN$TET2 == 0 & dataDN$ASXL1 == 1,]

d1 = cbind(rownames(datad),datad$patient,datad$stage,"DNMT3Amut")
d2 = cbind(rownames(datat),datat$patient,datat$stage,"TET2mut")
d3 = cbind(rownames(datab),datab$patient,datab$stage,"ASXL1mut")
sampled2 = as.data.frame(rbind(d1,d2,d3))
write.table(sampled2,"diagnosis.any.DNMT3A-TET2-ASXL1.samples.MDS.txt",sep="\t",row.names=FALSE)

dmut = colSums(datad[,1:7])
dwt = nrow(datad) - dmut
dd = cbind(dmut,dwt)

tmut = colSums(datat[,1:7])
twt = nrow(datat) - tmut
tt = cbind(tmut,twt)

bmut = colSums(datab[,1:7])
bwt = nrow(datab) - bmut
bt = cbind(bmut,bwt)

all(rownames(dd)==rownames(tt))
all(rownames(dd)==rownames(bt))

fdata = cbind(dd,tt,bt)

##### 3 group barplots
fdata2 = as.data.frame(fdata[c("FLT3","NPM1","CBL","SRSF2"),])
fdata2$dpct = fdata2$dmut/(fdata2$dmut+fdata2$dwt)*100
fdata2$tpct = fdata2$tmut/(fdata2$tmut+fdata2$twt)*100
fdata2$bpct = fdata2$bmut/(fdata2$bmut+fdata2$bwt)*100

fdata3 = as.data.frame(t(fdata2))[c("dpct","tpct","bpct"),]
rownames(fdata3) = c("DNMT3Amut_only","TET2mut_only","ASXL1mutonly")
fdata3$catg = rownames(fdata3)
library(reshape2)
fdata4 = melt(fdata3)
fdata4$catg = factor(fdata4$catg,levels=c("DNMT3Amut_only","TET2mut_only","ASXL1mutonly"))

pdf("dnmt3a-vs-tet2-vs-ASXL1-MDS.mut-freqs.barplot.pdf",height=3,width=10,useDingbats=FALSE)
ggbarplot(fdata4,x="catg",y="value",fill="catg",color=NA,position = position_dodge(0.8),palette=c("#189ACD","#EE8C00","#4A9649"),facet.by="variable",ncol=4,scales="free_y")+ theme_pubr()+rremove("legend") 
dev.off()

# FLT3
fisher.test(rbind(c(8,331),c(8,587)))
fisher.test(rbind(c(8,331),c(15,537)))
fisher.test(rbind(c(8,587),c(15,537)))

#NPM1
fisher.test(rbind(c(24,315),c(1,594)))
fisher.test(rbind(c(24,315),c(4,548)))
fisher.test(rbind(c(1,594),c(4,548)))

#CBL
fisher.test(rbind(c(10,329),c(50,545)))
fisher.test(rbind(c(10,329),c(55,497)))
fisher.test(rbind(c(50,545),c(55,497)))

#SRSF2
fisher.test(rbind(c(17,322),c(166,429)))$p.value
fisher.test(rbind(c(17,322),c(147,405)))$p.value
fisher.test(rbind(c(166,429),c(147,405)))