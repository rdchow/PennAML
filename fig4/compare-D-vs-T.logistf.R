#!/usr/bin/Rscript
# use TET2/ASXL1 vs DNMT3A definitions at time of diagnosis to guide subsequent labels
# NOT from each stage
library(ggplot2)
library(ggpubr)

# binary mutation matrix is from "matrix-muts.R"
data = read.table("AML_data.muts.chroms.matrix.txt",sep="\t",header=TRUE,row.names=1)
data = data[,colnames(data) != "TPMT"]

info = read.table("../fig2/uniqueStage-sample-names.txt",sep="\t")
colnames(info) = c("sampleID","patient","stage")


#####################################
# Diagnosis 

infoDN = info[info$stage=="dn",]
dataDN = data[rownames(data) %in% infoDN$sampleID,]
infoDN = infoDN[match(rownames(dataDN),infoDN$sampleID),]
all(rownames(dataDN) == infoDN$sampleID)

dataDN$patient = infoDN$patient
dataDN$stage = infoDN$stage

# define mutually exclusive groups
# can only have TET2 or DNMT3A mutations; mutually exclusive
dataDN2 = dataDN[dataDN$TET2 + dataDN$DNMT3A == 1,]
dataDN2$DTA = dataDN2$DNMT3A - dataDN2$TET2 # -1: TET2; +1: DNMT3A

library(logistf)
mydf = matrix(nrow = ncol(dataDN2)-3,ncol=5)
for (i in c(1:(ncol(dataDN2)-3))){
    mygene = colnames(dataDN2)[i]

    mydat1 = as.data.frame(cbind(dataDN2[,i],dataDN2$DTA))
    colnames(mydat1) = c("testGene","DTA")
    logist1 = logistf(testGene ~ DTA,data=mydat1)

    mydf[i,1] = mygene
    mydf[i,2] = as.numeric(coef(logist1)[2])
    mydf[i,3] = confint(logist1)[2,1]
    mydf[i,4] = confint(logist1)[2,2]
    mydf[i,5] = as.numeric(logist1$prob[2])
}
colnames(mydf) = c("gene","beta_DNMT3A_vs_TET2","lowCI","hiCI","pval")
mydf = as.data.frame(mydf)
mydf = mydf[! mydf$gene %in% c("DNMT3A","TET2"),]
mydf = mydf[order(mydf$pval),]

write.table(mydf,"DNMT3A-vs-TET2.exclusive.diagnosis.logistf.txt",sep="\t",row.names=FALSE)

####################
#### get mut frequencies in each group

# mutually exclusive
datad = dataDN2[dataDN2$DTA == 1,] # DNMT3A
datat = dataDN2[dataDN2$DTA == -1,] # TET2

d1 = cbind(rownames(datad),datad$patient,datad$stage,"DNMT3Amut")
d2 = cbind(rownames(datat),datat$patient,datat$stage,"TET2mut")
sampled2 = as.data.frame(rbind(d1,d2))
write.table(sampled2,"diagnosis.any.DNMT3A-TET2.samples.txt",sep="\t",row.names=FALSE)

dmut = colSums(datad[,1:(ncol(datad)-3)])
dwt = nrow(datad) - dmut
dd = cbind(dmut,dwt)

tmut = colSums(datat[,1:(ncol(datat)-3)])
twt = nrow(datat) - tmut
tt = cbind(tmut,twt)

all(rownames(dd)==rownames(tt))
fdata = as.data.frame(cbind(dd,tt))
fdata$d_mutfreq = fdata$dmut/(fdata$dmut+fdata$dwt)*100
fdata$t_mutfreq = fdata$tmut/(fdata$tmut+fdata$twt)*100

######################################
# dot plots
# diagnosis
library(ggplot2)
library(ggpubr)
library(viridis)
library(rcartocolor)

data = read.table("DNMT3A-vs-TET2.exclusive.diagnosis.logistf.txt",sep="\t",header=TRUE,row.names=1)
fdata = fdata[rownames(fdata) %in% rownames(data),]
fdata = fdata[match(rownames(data),rownames(fdata)),]
all(rownames(fdata) == rownames(data))

data$nlp = log10(data$pval)*-1
data[data$nlp > 3,"nlp"] = 3 # cap at 3
data$gene = rownames(data)
data$d_mutfreq = fdata$d_mutfreq
data$t_mutfreq = fdata$t_mutfreq

data$sumPct = data$d_mutfreq + data$t_mutfreq

pdf("diagnosis.DNMT3A-vs-TET2.mutfreqs.logistf.pdf",height=5,width=4,useDingbats=FALSE)
ggscatter(data,x="d_mutfreq", y="t_mutfreq",fill="sumPct",shape=21,color="black",repel=TRUE,alpha=0.75,size="nlp",xlim=c(0,70),ylim=c(0,40),label="gene", label.select=unique(c(data[data$pval < 0.05,"gene"][1:10],data[order(data$d_mutfreq,decreasing=TRUE),"gene"][1:5],data[order(data$t_mutfreq,decreasing=TRUE),"gene"][1:8]))) + 
geom_text(data=data[data$pval <0.05,],size=8,hjust=0.5,vjust = 0.75,color="black",label="*") + 
scale_size_continuous(range = c(2, 10),limits=c(0,3)) + scale_fill_gradientn(colors=carto_pal(7,"ag_GrnYl")) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)
dev.off()


#####################################################################################
#####################################
# REL1
data = read.table("AML_data.muts.chroms.matrix.txt",sep="\t",header=TRUE,row.names=1)
data = data[,colnames(data) != "TPMT"]
info = read.table("../fig2/uniqueStage-sample-names.txt",sep="\t")
colnames(info) = c("sampleID","patient","stage")

infoREL = info[info$stage=="REL1",]
dataREL = data[rownames(data) %in% infoREL$sampleID,]
infoREL = infoREL[match(rownames(dataREL),infoREL$sampleID),]
all(rownames(dataREL) == infoREL$sampleID)

dataREL$patient = infoREL$patient
dataREL$stage = infoREL$stage

# get labels from diagnosis
labels = read.table("diagnosis.any.DNMT3A-TET2.samples.txt",sep="\t",header=TRUE)

# filter dataREL to REL1 samples from patients with a DNMT3A or TET2 exclusive mutation at diagnosis (so DNMT3A / TET2 may actually be *lost* at time of REL1)
dataREL = dataREL[dataREL$patient %in% labels$V2,]
labelsREL = labels[labels$V2 %in% dataREL$patient,]
labelsREL = labelsREL[match(dataREL$patient,labelsREL$V2),]
all(dataREL$patient == labelsREL$V2)
dataREL$DTA = labelsREL$V4

dataREL[dataREL$DTA == "DNMT3Amut","DTA"] = 1
dataREL[dataREL$DTA == "TET2mut","DTA"] = -1

library(logistf)
mydf = matrix(nrow = ncol(dataREL)-3,ncol=5)
for (i in 1:(ncol(dataREL)-3)){
    mygene = colnames(dataREL)[i]
    
    mydat1 = as.data.frame(cbind(dataREL[,i],dataREL$DTA))
    colnames(mydat1) = c("testGene","DTA")
    mydat1$DTA = as.factor(mydat1$DTA)
    mydat1$testGene = as.factor(mydat1$testGene)
    logist1 = logistf(testGene ~ DTA,data=mydat1)

    mydf[i,1] = mygene
    mydf[i,2] = as.numeric(coef(logist1)[2])
    mydf[i,3] = confint(logist1)[2,1]
    mydf[i,4] = confint(logist1)[2,2]
    mydf[i,5] = as.numeric(logist1$prob[2])
}
colnames(mydf) = c("gene","beta_DNMT3A_vs_TET2","lowCI","hiCI","pval")
mydf = as.data.frame(mydf)
mydf = mydf[! mydf$gene %in% c("DNMT3A","TET2"),]
mydf = mydf[order(mydf$pval),]

write.table(mydf,"DNMT3A-vs-TET2.exclusive.REL1.logistf.txt",sep="\t",row.names=FALSE)

####################
#### get mut frequencies in each group

# mutually exclusive
datad = dataREL[dataREL$DTA == 1,]
datat = dataREL[dataREL$DTA == -1,]

d1 = cbind(rownames(datad),datad$patient,datad$stage,"DNMT3Amut")
d2 = cbind(rownames(datat),datat$patient,datat$stage,"TET2mut")
sampled2 = as.data.frame(rbind(d1,d2))

dmut = colSums(datad[,1:(ncol(datad)-3)])
dwt = nrow(datad) - dmut
dd = cbind(dmut,dwt)

tmut = colSums(datat[,1:(ncol(datat)-3)])
twt = nrow(datat) - tmut
tt = cbind(tmut,twt)

all(rownames(dd)==rownames(tt))
fdata = as.data.frame(cbind(dd,tt))
fdata$d_mutfreq = fdata$dmut/(fdata$dmut+fdata$dwt)*100
fdata$t_mutfreq = fdata$tmut/(fdata$tmut+fdata$twt)*100

######################################
# dot plots
# diagnosis
library(ggplot2)
library(ggpubr)
library(viridis)
library(rcartocolor)

data = read.table("DNMT3A-vs-TET2.exclusive.REL1.logistf.txt",sep="\t",header=TRUE,row.names=1)
fdata = fdata[rownames(fdata) %in% rownames(data),]
fdata = fdata[match(rownames(data),rownames(fdata)),]
all(rownames(fdata) == rownames(data))

data$nlp = log10(data$pval)*-1
data[data$nlp > 3,"nlp"] = 3 # cap at 3
data$gene = rownames(data)
data$d_mutfreq = fdata$d_mutfreq
data$t_mutfreq = fdata$t_mutfreq

data$sumPct = data$d_mutfreq + data$t_mutfreq

pdf("REL1.DNMT3A-vs-TET2.mutfreqs.logistf.pdf",height=5,width=4,useDingbats=FALSE)
ggscatter(data,x="d_mutfreq", y="t_mutfreq",fill="sumPct",shape=21,color="black",repel=TRUE,alpha=0.75,size="nlp",xlim=c(0,70),ylim=c(0,40),label="gene", label.select=unique(c(data[data$pval < 0.05,"gene"][1:10],data[order(data$d_mutfreq,decreasing=TRUE),"gene"][1:5],data[order(data$t_mutfreq,decreasing=TRUE),"gene"][1:8]))) + 
geom_text(data=data[data$pval <0.05,],size=8,hjust=0.5,vjust = 0.75,color="black",label="*") + 
scale_size_continuous(range = c(2, 10),limits=c(0,3)) + scale_fill_gradientn(colors=carto_pal(7,"ag_GrnYl")) + geom_abline(slope=1,intercept=0,linetype="dashed") + theme(aspect.ratio=1)
dev.off()
