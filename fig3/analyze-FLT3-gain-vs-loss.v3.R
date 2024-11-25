#!/usr/bin/Rscript
setwd("C:/Residency/Bowman/AML-longitudinal/version3/fig3")
library(ggpubr)
library(rcartocolor)
library(logistf)
library(dplyr)
library(reshape2)
library(stringr)

##### Part 1: FLT3 gain at REL1 #####

# first, get the binary mutation profiles at diagnosis
diagdata = read.table("../fig4/AML_data.muts.chroms.matrix.diagnosis.txt",sep="\t",header=TRUE,row.names=1)
rownames(diagdata) = str_split_i(rownames(diagdata),":",i=2)

# Establish which patients gained a FLT3 mutation at REL1
data=read.table("../fig2/pt-level-mut-freqs.txt",sep="\t",header=TRUE)
data = data[data$gene != "TPMT",]
dcdata = data[(!is.na(data$diagnosis) & !is.na(data$REL1)) & (data$diagnosis >0 | data$REL1 > 0),]

subdata = dcdata[dcdata$gene == "FLT3",]

# choose the variant with the largest delta VAF
subdata$abs_rel_diag = abs(subdata$REL1 - subdata$diagnosis)
library(dplyr)
subdata2 = subdata %>% group_by(patientID) %>% top_n(1, abs_rel_diag)

# binarize variant status changes
subdata2$subgain = 0
subdata2[subdata2$diagnosis == 0 & subdata2$REL1 > 0,"subgain"] = 1
subdata2$subloss = 0
subdata2[subdata2$diagnosis > 0 & subdata2$REL1 == 0,"subloss"] = 1
subdata2 = as.data.frame(subdata2)

# Now go back and filter/analyze the diagnosis mutation matrix, based on gain/loss of a FLT3 variant (with the largest delta VAF)
diagdata$subgain = 0
diagdata[rownames(diagdata) %in% subdata2[subdata2$subgain == 1,"patientID"], "subgain"] = 1 # 12
diagdata$subloss = 0
diagdata[rownames(diagdata) %in% subdata2[subdata2$subloss == 1,"patientID"], "subloss"] = 1 # 16

# add flt3i treatment info
info = read.table("../preprocess/sample-info.p6.ELNlab.uniqFilt.txt",sep="\t",header=TRUE)
info = info[info$simpleTimepoint %in% c("REL1"),]
# if the REL1 sample was treated with FLT3i and timeBinInhib is "post" -> then FLT3i=yes
info$REL1_FLT3i = 0
info[is.na(info)] = "."
info[info$simplestInhibitor == "FLT3i" & info$timeBinInhib == "post","REL1_FLT3i"] = 1 # 5 REL1 samples after FLT3 treatment

diagdata2 = diagdata[rownames(diagdata) %in% info$patientID,] # filter the diagnosis matrix to patients with matched REL1 genomics
diagdata2 = diagdata2[,colSums(diagdata2) >= 3]
info2 = info[info$patientID %in% rownames(diagdata2),]
info2 = info2[match(rownames(diagdata2),info2$patientID),]
all(info2$patientID == rownames(diagdata2))

diagdata2$REL1_FLT3i = info2$REL1_FLT3i

diagdata2 = diagdata2[,colnames(diagdata2) != "FLT3"]

# Univariate
mydf = matrix(nrow=ncol(diagdata2),ncol=15)
for (i in 1:ncol(diagdata2)){
    mygene = colnames(diagdata2)[i]

    if (! mygene %in% c("subgain","subloss")){

        mydat1 = as.data.frame(cbind(diagdata2[,i],diagdata2$subgain))
        colnames(mydat1) = c("testGene","subgain")
        mytab = table(mydat1$testGene,mydat1$subgain)
        mydf[i,1] = mygene
        mydf[i,2] = mytab[2,1]
        mydf[i,3] = mytab[2,2]
        mydf[i,4] = mytab[2,1] + mytab[2,2]
        #fish1 = fisher.test(as.data.frame.matrix(table(mydat1[,1],mydat1[,2])))
        logist1 = logistf(subgain ~ testGene,data=mydat1)
        mydf[i,5] = as.numeric(coef(logist1)[2])
        mydf[i,6] = confint(logist1)[2,1]
        mydf[i,7] = confint(logist1)[2,2]
        mydf[i,8] = as.numeric(logist1$prob[2])

        mydat2 = as.data.frame(cbind(diagdata2[,i],diagdata2$subloss))
        colnames(mydat2) = c("testGene","subloss")
        mytab = table(mydat2$testGene,mydat2$subloss)
        mydf[i,9] = mytab[2,1]
        mydf[i,10] = mytab[2,2]
        mydf[i,11] = mytab[2,1] + mytab[2,2]

        logist2 = logistf(subloss ~ testGene,data=mydat2)
        mydf[i,12] = as.numeric(coef(logist2)[2])
        mydf[i,13] = confint(logist2)[2,1]
        mydf[i,14] = confint(logist2)[2,2]
        mydf[i,15] = as.numeric(logist2$prob[2])
    }
}

colnames(mydf) = c("gene","testGeneMut.noFLT3gain","testGeneMut.FLT3gain","totalTestGeneMut","beta_FLT3gain","beta_FLT3gain_lowCI","beta_FLT3gain_hiCI","pval_FLT3gain", "testGeneMut.noFLT3loss","testGeneMut.FLT3loss","totalTestGeneMut","beta_FLT3loss","beta_FLT3loss_lowCI","beta_FLT3loss_hiCI","pval_FLT3loss")

write.table(mydf,"fixed.FLT3-binary-gain-vs-loss.associatedMuts.logistf.v2.txt",sep="\t",row.names=FALSE)


## mutlivariate
# Univariate
mydf = matrix(nrow=ncol(diagdata2)-1,ncol=9)
for (i in 1:(ncol(diagdata2)-1)){
    mygene = colnames(diagdata2)[i]

    if (! mygene %in% c("subgain","subloss")){

        mydat1 = as.data.frame(cbind(diagdata2[,i],diagdata2$subgain,diagdata2$REL1_FLT3i))
        colnames(mydat1) = c("testGene","subgain","zREL1_FLT3i")

        logist1 = logistf(subgain ~ testGene+zREL1_FLT3i,data=mydat1)
        mydf[i,1] = mygene
        mydf[i,2] = as.numeric(coef(logist1)[2])
        mydf[i,3] = confint(logist1)[2,1]
        mydf[i,4] = confint(logist1)[2,2]
        mydf[i,5] = as.numeric(logist1$prob[2])

        mydat2 = as.data.frame(cbind(diagdata2[,i],diagdata2$subloss,diagdata2$REL1_FLT3i))
        colnames(mydat2) = c("testGene","subloss","zREL1_FLT3i")
        logist2 = logistf(subloss ~ testGene+zREL1_FLT3i,data=mydat2)
        mydf[i,6] = as.numeric(coef(logist2)[2])
        mydf[i,7] = confint(logist2)[2,1]
        mydf[i,8] = confint(logist2)[2,2]
        mydf[i,9] = as.numeric(logist2$prob[2])
    }
}

colnames(mydf) = c("gene","beta_FLT3gain","beta_FLT3gain_lowCI","beta_FLT3gain_hiCI","pval_FLT3gain", "beta_FLT3loss","beta_FLT3loss_lowCI","beta_FLT3loss_hiCI","pval_FLT3loss")

write.table(mydf,"fixed.FLT3-binary-gain-vs-loss.associatedMuts.logistf.v2.FLT3iadjust.txt",sep="\t",row.names=FALSE)


####################################
## forest plots
setwd("C:/Residency/Bowman/AML-longitudinal/version3/fig3")
library(ggpubr)
library(rcartocolor)
library(colorspace)
library(dplyr)

### without adjustment for FLT3i
mydf = read.table("fixed.FLT3-binary-gain-vs-loss.associatedMuts.logistf.v2.txt",sep="\t",header=TRUE)
mydf2 = as.data.frame(mydf)
mydf2 = mydf2[!is.na(mydf2$gene),]
rownames(mydf2) = mydf2$gene
mydf2 = mydf2 %>% mutate_at(c(2:ncol(mydf2)), as.numeric)
mydf2$gene = factor(mydf2$gene,levels=rev(mydf2$gene))

#FLT3 gain
g1= ggplot(mydf2,aes(x=beta_FLT3gain,y=gene)) + geom_errorbar(aes(xmin=beta_FLT3gain_lowCI, xmax=beta_FLT3gain_hiCI),width=0.1,linewidth=0.4) + geom_point(size=3,color="#ED5D00")+ theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + scale_color_distiller(palette="BrBG") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3gain,digit=3)),nudge_y=0.2)

#FLT3 loss
#mydf2 = mydf2[order(mydf2$beta_FLT3loss,decreasing=TRUE),]
#mydf2$gene = factor(mydf2$gene,levels=rev(mydf2$gene))
g2= ggplot(mydf2,aes(x=beta_FLT3loss,y=gene)) +  geom_errorbar(aes(xmin=beta_FLT3loss_lowCI, xmax=beta_FLT3loss_hiCI),width=0.1,linewidth=0.4) + geom_point(size=3,color="#0095A4") + theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3loss,digit=3)),nudge_y=0.2)

pdf("FLT3-diag-vs-REL1.gain-vs-loss.logistf.forest.pdf",height=4,width=7,useDingbats=FALSE)
ggarrange(g1,g2,ncol=2)
dev.off()


### adjusted for FLT3i
mydf = read.table("fixed.FLT3-binary-gain-vs-loss.associatedMuts.logistf.v2.FLT3iadjust.txt",sep="\t",header=TRUE)
mydf2 = as.data.frame(mydf)
mydf2 = mydf2[!is.na(mydf2$gene),]
rownames(mydf2) = mydf2$gene
mydf2 = mydf2 %>% mutate_at(c(2:ncol(mydf2)), as.numeric)
mydf2$gene = factor(mydf2$gene,levels=rev(mydf2$gene))

#FLT3 gain
g1= ggplot(mydf2,aes(x=beta_FLT3gain,y=gene)) + geom_errorbar(aes(xmin=beta_FLT3gain_lowCI, xmax=beta_FLT3gain_hiCI),width=0.1,linewidth=0.4) + geom_point(size=3,color="#ED5D00")+ theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + scale_color_distiller(palette="BrBG") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3gain,digit=3)),nudge_y=0.2)

#FLT3 loss
g2= ggplot(mydf2,aes(x=beta_FLT3loss,y=gene)) +  geom_errorbar(aes(xmin=beta_FLT3loss_lowCI, xmax=beta_FLT3loss_hiCI),width=0.1,linewidth=0.4) + geom_point(size=3,color="#0095A4") + theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3loss,digit=3)),nudge_y=0.2)

pdf("FLT3-diag-vs-REL1.gain-vs-loss.logistf.FLT3i_adjust.forest.pdf",height=4,width=7,useDingbats=FALSE)
ggarrange(g1,g2,ncol=2)
dev.off()