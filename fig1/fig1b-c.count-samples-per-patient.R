#!/usr/bin/Rscript
library(ggpubr)
library(ggplot2)
library(reshape2)

data = read.table("../preprocess/sample-info.p6.txt",sep="\t",header=TRUE,row.names=1) # taking all samples from patients treated with anything; there may be multiple CR1, REL1 samples etc per patient
data = data[data$simpleTimepoint != "MDS?",]

mydf = as.data.frame.matrix(table(data$patientID,data$simpleTimepoint))
mydf2 = mydf # only count one sample for a given time point in a given patient
mydf2[mydf2 > 1] = 1
mmdf = as.data.frame(t(rbind(colSums(mydf2),colSums(mydf2))))
mmdf = mmdf[! rownames(mmdf) %in% "REC",]
mmdf$type = rownames(mmdf)
colnames(mmdf) = c("patientCt","patientCt2","type")
mmdf$type = factor(mmdf$type,levels=c("dn","CR1","REL1","REF1","CR2","REL2","REF2"))
mmdf$pct = round(mmdf$patientCt/nrow(data)*100, digits=1)

g1=ggbarplot(data=mmdf,x="type",y="pct",label=mmdf$patientCt,fill="type",color=NA,palette=c("#162252","#403369","#5B5992","#AE93BE","#B4DAE5","#D8AF39","#D86339")) +rremove("legend")+scale_y_continuous(expand=c(0,0))

####
mydf = as.data.frame.matrix(table(data$patientID,data$simpleTimepoint))
mydf$sum = rowSums(mydf)

metadf = as.data.frame.matrix(cbind(table(mydf$sum),table(mydf$sum)))
metadf$sampleCt = rownames(metadf)
colnames(metadf) = c("sum","sum2","sampleCt")
metadf$pct = round(metadf$sum/nrow(data)*100,digits=1)


g2=ggbarplot(data=metadf,x="sampleCt",y="pct",fill="sampleCt",color=NA,label=metadf$sum,palette=c("#274638","#2D715E","#43A57C","#58A449","#CEC917","#CE8D17"))+rremove("legend")+scale_y_continuous(expand=c(0,0))

pdf("num-sample-per-patient.allTreatments.pdf",height=6,width=6,useDingbats=FALSE)
ggarrange(g2,g1,ncol=1)
dev.off()

# plot of treatment categories
#!/usr/bin/Rscript
library(ggpubr)
library(ggplot2)
library(reshape2)

data = read.table("../preprocess/sample-info.p6.txt",sep="\t",header=TRUE,row.names=1) # taking all samples from patients treated with anything; there may be multiple CR1, REL1 samples etc per patient
diagdata = data[data$simpleTimepoint != "MDS?" & data$simpleTimepoint == "dn",]

mydf = as.data.frame(table(diagdata$simpleInitialTreatment))
colnames(mydf) = c("treatment","sampleCt")
mydf$pct = round(mydf$sampleCt/nrow(diagdata)*100,digits=1)
mydf = mydf[order(mydf$pct,decreasing=TRUE),] 
mydf$treatment = factor(mydf$treatment,levels=rev(mydf$treatment))
mydf$color = "gray44"
#mydf[mydf$treatment == "anthracycline_nucleosideAnalog","color"] = "red3"


pdf("treatmentTypes-barplot.pdf",height=6,width=8,useDingbats=FALSE)
ggbarplot(data=mydf,x="treatment",y="pct",fill="color",color=NA,label=mydf$sampleCt,palette=c("gray44"),orientation="horizontal")+rremove("legend")+scale_y_continuous(expand=c(0,0))
dev.off()




