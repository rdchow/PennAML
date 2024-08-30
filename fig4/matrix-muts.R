#!/usr/bin/Rscript
# convert mutation file to matrix form
library(reshape2)

samples = read.table("../fig2/uniqueStage-sample-names.txt",sep="\t",header=FALSE)

data=read.table("../preprocess/AML_data.muts.chroms.uniqFilt.dropDiscordantPanelGenes.txt",sep="\t",header=FALSE,row.names=NULL)
colnames(data) = c("sample","gene","variant","VAF","type","sampleType","patientID")
data = data[data$sample %in% samples$V1,]
#write.table(data,"../inputData/AML_data.muts.chroms.inductOnly.uniqFilt.txt",sep="\t",row.names=FALSE)

data$VAF[data$gene != ""] = "1"

ldata = dcast(data,sample~gene,value.var="VAF",length) # binary
ldata[is.na(ldata)] = 0
#ldata[ldata=="-Inf"] = 0
ldata2= ldata[,-2]
rownames(ldata2)=ldata[,1]
ldata2 = ldata2[,-1]
ldata2[ldata2>1]=1
ldata2 = cbind(rownames(ldata2),ldata2)
colnames(ldata2)[1] = "sampleID"
write.table(ldata2,"AML_data.muts.chroms.matrix.txt",sep="\t",row.names=FALSE)

data2 = data[data$sampleType == "dn.AML",]
colnames(data2) = c("sample","gene","variant","VAF","type","sampleType","patientID")
write.table(data2,"../preprocess/AML_data.muts.chroms.uniqFilt.diagnosis.txt",sep="\t",row.names=FALSE)

ldata = dcast(data2,sample+patientID~gene,value.var="VAF",length) # binary
ldata2= ldata[,-3]
rownames(ldata2)=paste(ldata[,1],":",ldata[,2],sep="")
ldata2 = ldata2[,c(-1,-2)]
ldata2[ldata2>1]=1
ldata2 = cbind(rownames(ldata2),ldata2)
colnames(ldata2)[1] = "sampleID"
write.table(ldata2,"AML_data.muts.chroms.matrix.diagnosis.txt",sep="\t",row.names=FALSE)


######
data3 = data[data$sampleType == "rel1.AML",]
colnames(data3) = c("sample","gene","variant","VAF","type","sampleType","patientID")
write.table(data3,"../preprocess/AML_data.muts.chroms.uniqFilt.REL1.txt",sep="\t",row.names=FALSE)

ldata = dcast(data3,sample+patientID~gene,value.var="VAF",length) # binary
ldata3= ldata[,-3]
rownames(ldata3)=paste(ldata[,1],":",ldata[,2],sep="")
ldata3 = ldata3[,c(-1,-2)]
ldata3[ldata3>1]=1
ldata3 = cbind(rownames(ldata3),ldata3)
colnames(ldata3)[1] = "sampleID"
write.table(ldata3,"AML_data.muts.chroms.matrix.relapse.txt",sep="\t",row.names=FALSE)

######
data3 = data[data$sampleType == "CR1.AML",]
colnames(data3) = c("sample","gene","variant","VAF","type","sampleType","patientID")
write.table(data3,"../preprocess/AML_data.muts.chroms.uniqFilt.CR1.txt",sep="\t",row.names=FALSE)

ldata = dcast(data3,sample+patientID~gene,value.var="VAF",length) # binary
ldata3= ldata[,-3]
rownames(ldata3)=paste(ldata[,1],":",ldata[,2],sep="")
ldata3 = ldata3[,c(-1,-2)]
ldata3[ldata3>1]=1
ldata3 = cbind(rownames(ldata3),ldata3)
colnames(ldata3)[1] = "sampleID"
write.table(ldata3,"AML_data.muts.chroms.matrix.CR1.txt",sep="\t",row.names=FALSE)