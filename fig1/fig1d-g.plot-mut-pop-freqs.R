#!/usr/bin/Rscript
library(ggpubr)
library(reshape2)
library(viridis)
library(ggplot2)

# mut-frequencies.txt comes from running calc-mut-frequency.pl
data=read.table("mut-frequencies.txt",sep="\t",skip=1,header=TRUE)
data = data[data$gene != "TPMT",]
data[,3] = data[,3]*100
data[,5] = data[,5]*100
data[,7] = data[,7]*100
data[,9] = data[,9]*100

# Diagnosis vs CR1
# calculate fisher stats
#182 diagnosis samples, 117 CR1 samples, 76 Rel1 samples, 27 Ref1 samples, 19 CR2 samples, 13 Rel2 samples, 4 Ref2 samples

numdiag = 182
numcr1 = 119
numrel1 = 76
numref1 = 27

mypvals = ""
for (i in 1:nrow(data)){
    ddf = rbind(c(data$DiagnosisMutCt[i],numdiag-data$DiagnosisMutCt[i]), c(data$CR1mutCt[i],numcr1-data$CR1mutCt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}
data$fisher.diag_cr1 = as.numeric(mypvals)
#data$fisher.diag_cr1.bin =  as.numeric(data$fisher.diag_cr1<0.05)
data$fisher.diag_cr1.nlp = log10(data$fisher.diag_cr1)*-1
data$fisher.diag_cr1.nlp[data$fisher.diag_cr1.nlp <= log10(0.05)*-1] = 0

g1=ggscatter(data,x= "DiagnosisMutFreq", y="CR1mutFreq",fill="DiagnosisMutFreq",shape=21,color="black",repel=TRUE,alpha=0.75,size="fisher.diag_cr1.nlp",xlim=c(0,50),ylim=c(0,50),label="gene",label.select=unique(c(data[data$fisher.diag_cr1 < 0.05,"gene"][1:10],data[order(data$DiagnosisMutFreq,decreasing=TRUE),"gene"][1:8]))) + scale_size_continuous(range = c(2, 8),limits=c(0,13)) + scale_fill_viridis(option="rocket",direction=-1) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)


# CR1 vs Relapse 1
# calculate fisher stats
mypvals = ""
for (i in 1:nrow(data)){
    ddf = rbind(c(data$Rel1MutCt[i],numrel1-data$Rel1MutCt[i]), c(data$CR1mutCt[i],numcr1-data$CR1mutCt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}
data$fisher.cr1_rel1 = as.numeric(mypvals)
data$fisher.cr1_rel1.nlp = log10(data$fisher.cr1_rel1)*-1
data$fisher.cr1_rel1.nlp[data$fisher.cr1_rel1.nlp <= log10(0.05)*-1] = 0


g2=ggscatter(data,x= "CR1mutFreq", y="Rel1MutFreq",fill="Rel1MutFreq",shape=21,color="black",repel=TRUE,alpha=0.75,size="fisher.cr1_rel1.nlp",xlim=c(0,50),ylim=c(0,50),label="gene",label.select=unique(c(data[data$fisher.cr1_rel1 < 0.05,"gene"][1:10],data[order(data$Rel1MutFreq,decreasing=TRUE),"gene"][1:8]))) + scale_size_continuous(range = c(2, 8),limits=c(0,13)) + scale_fill_viridis(option="rocket",direction=-1) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)


# Diagnosis vs Relapse 1
# calculate fisher stats
mypvals = ""
for (i in 1:nrow(data)){
    ddf = rbind(c(data$DiagnosisMutCt[i],numdiag-data$DiagnosisMutCt[i]), c(data$Rel1MutCt[i],numrel1-data$Rel1MutCt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}
data$fisher.diag_rel1 = as.numeric(mypvals)
data$fisher.diag_rel1.nlp = log10(data$fisher.diag_rel1)*-1
data$fisher.diag_rel1.nlp[data$fisher.diag_rel1.nlp <= log10(0.05)*-1] = 0
g3=ggscatter(data,x= "DiagnosisMutFreq", y="Rel1MutFreq",fill="DiagnosisMutFreq",shape=21,color="black",repel=TRUE,alpha=0.75,size="fisher.diag_rel1.nlp",xlim=c(0,50),ylim=c(0,50),label="gene",label.select=unique(c(data[data$fisher.diag_rel1 < 0.05,"gene"][1:10],data[order(data$DiagnosisMutFreq,decreasing=TRUE),"gene"][1:8]))) + scale_size_continuous(range = c(2, 8),limits=c(0,13)) + scale_fill_viridis(option="rocket",direction=-1) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)



# REF1 vs Relapse 1
# calculate fisher stats

mypvals = ""
for (i in 1:nrow(data)){
    ddf = rbind(c(data$Ref1MutCt[i],numref1-data$Ref1MutCt[i]), c(data$Rel1MutCt[i],numrel1-data$Rel1MutCt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}
data$fisher.ref1_rel1 = as.numeric(mypvals)
data$fisher.ref1_rel1.nlp = log10(data$fisher.ref1_rel1)*-1
data$fisher.ref1_rel1.nlp[data$fisher.ref1_rel1.nlp <= log10(0.05)*-1] = 0
write.table(data,"fisher-mutFreqs.diagnosis-CR1-REL1-REF1.txt",sep="\t",row.names=FALSE)


g4=ggscatter(data,x= "Ref1MutFreq", y="Rel1MutFreq",fill="Ref1MutFreq",shape=21,color="black",repel=TRUE,alpha=0.75,size="fisher.ref1_rel1.nlp",xlim=c(0,50),ylim=c(0,50),label="gene",label.select=unique(c(data[data$fisher.ref1_rel1 < 0.05,"gene"][1:10],data[order(data$Ref1MutFreq,decreasing=TRUE),"gene"][1:8]))) + scale_size_continuous(range = c(2, 8),limits=c(0,13)) + scale_fill_viridis(option="rocket",direction=-1) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)


pdf("population-level-mutFreq-scatter.pdf",height=4,width=18,useDingbats=FALSE)
ggarrange(g1,g2,g3,g4,ncol=4)
dev.off()

