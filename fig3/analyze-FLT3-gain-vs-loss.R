#!/usr/bin/Rscript


library(ggpubr)
library(rcartocolor)
library(logistf)
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
subdata2[subdata2$diagnosis == 0 & subdata2$REL1 > 0,"subgain"] = 1 #14
subdata2$subloss = 0
subdata2[subdata2$diagnosis > 0 & subdata2$REL1 == 0,"subloss"] = 1 #13
subdata2$subshared = 0
subdata2[subdata2$diagnosis > 0 & subdata2$REL1 > 0,"subshared"] = 1 #21

subdata2 = as.data.frame(subdata2)

library(reshape2)
mutdata = dcast(dcdata, patientID ~ gene,value.var="diagnosis")

mutdata2 = mutdata[mutdata$patientID %in% subdata2$patientID,]
subdata3 = subdata2[match(mutdata2$patientID,subdata2$patientID),]
mutdata2 = mutdata2[,colSums(mutdata2)>=3]
all(subdata3$patientID == mutdata2$patientID) 

mutdata2[mutdata2>1] = 1
mutdata2 = mutdata2[,colnames(mutdata2) != "FLT3"]

mydf = matrix(nrow=ncol(mutdata2)-1,ncol=13)
for (i in 2:ncol(mutdata2)){
    mygene = colnames(mutdata2)[i]

    mydat1 = as.data.frame(cbind(mutdata2[,i],subdata3$subgain))
    colnames(mydat1) = c("testGene","subgain")
    #fish1 = fisher.test(as.data.frame.matrix(table(mydat1[,1],mydat1[,2])))
    logist1 = logistf(subgain ~ testGene,data=mydat1)

    mydat2 = as.data.frame(cbind(mutdata2[,i],subdata3$subloss))
    colnames(mydat2) = c("testGene","subloss")
    #fish2 = fisher.test(as.data.frame.matrix(table(mydat2[,1],mydat2[,2])))
    logist2 = logistf(subloss ~ testGene,data=mydat2)

    mydat3 = as.data.frame(cbind(mutdata2[,i],subdata3$subshared))
    colnames(mydat3) = c("testGene","subshared")
    #fish3 = fisher.test(as.data.frame.matrix(table(mydat3[,1],mydat3[,2])))
    logist3 = logistf(subshared ~ testGene,data=mydat3)

    mydf[i-1,1] = mygene

    mydf[i-1,2] = as.numeric(coef(logist1)[2])
    mydf[i-1,3] = confint(logist1)[2,1]
    mydf[i-1,4] = confint(logist1)[2,2]
    mydf[i-1,5] = as.numeric(logist1$prob[2])

    mydf[i-1,6] = as.numeric(coef(logist2)[2])
    mydf[i-1,7] = confint(logist2)[2,1]
    mydf[i-1,8] = confint(logist2)[2,2]
    mydf[i-1,9] = as.numeric(logist2$prob[2])

    mydf[i-1,10] = as.numeric(coef(logist3)[2])
    mydf[i-1,11] = confint(logist3)[2,1]
    mydf[i-1,12] = confint(logist3)[2,2]
    mydf[i-1,13] = as.numeric(logist3$prob[2])
}

colnames(mydf) = c("gene","beta_FLT3gain","beta_FLT3gain_lowCI","beta_FLT3gain_hiCI","pval_FLT3gain","beta_FLT3loss","beta_FLT3loss_lowCI","beta_FLT3loss_hiCI","pval_FLT3loss","beta_FLT3shared","beta_FLT3shared_lowCI","beta_FLT3shared_hiCI","pval_FLT3shared")

write.table(mydf,"FLT3-binary-gain-vs-loss.associatedMuts.logistf.txt",sep="\t",row.names=FALSE)


## forest plots
mydf2 = as.data.frame(mydf)
rownames(mydf2) = mydf2$gene
mydf2 = mydf2 %>% mutate_at(c(2:ncol(mydf2)), as.numeric)
mydf2$gene = factor(mydf2$gene,levels=rev(mydf2$gene))
library(ggpubr)
library(rcartocolor)
library(colorspace)
library(dplyr)

#FLT3 gain
g1= ggplot(mydf2,aes(x=beta_FLT3gain,y=gene)) + geom_errorbar(aes(xmin=beta_FLT3gain_lowCI, xmax=beta_FLT3gain_hiCI),width=0.1,size=0.4) + geom_point(size=5,color="#ED5D00")+ theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + scale_color_distiller(palette="BrBG") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3gain,digit=3)),nudge_y=0.2)

#FLT3 loss
#mydf2 = mydf2[order(mydf2$beta_FLT3loss,decreasing=TRUE),]
#mydf2$gene = factor(mydf2$gene,levels=rev(mydf2$gene))
g2= ggplot(mydf2,aes(x=beta_FLT3loss,y=gene)) +  geom_errorbar(aes(xmin=beta_FLT3loss_lowCI, xmax=beta_FLT3loss_hiCI),width=0.1,size=0.4) + geom_point(size=5,color="#0095A4") + theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3loss,digit=3)),nudge_y=0.2)

#FLT3 shared
g3= ggplot(mydf2,aes(x=beta_FLT3shared,y=gene)) +  geom_errorbar(aes(xmin=beta_FLT3shared_lowCI, xmax=beta_FLT3shared_hiCI),width=0.1,size=0.4) + geom_point(size=5,color="#0095A4") + theme_pubr()+geom_vline(xintercept=0,linetype=2,color="gray44") +rremove("legend") + geom_text(aes(x=0.2,y=gene,label=signif(pval_FLT3shared,digit=3)),nudge_y=0.2)

pdf("FLT3-diag-vs-REL1.gain-vs-loss.logistf.forest.pdf",height=4,width=7,useDingbats=FALSE)
ggarrange(g1,g2,ncol=2)
dev.off()