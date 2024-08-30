#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(viridis)
library(logistf)

# diagnosis
data = read.table("../AML_data.muts.chroms.matrix.diagnosis.txt",sep="\t",header=TRUE,row.names=1)
data = data[,colnames(data) != "TPMT"]
data = data[,colSums(data) >= 2] # at least 3 mutations

mydf = matrix(nrow = choose(ncol(data),2),ncol=7)
mygenes = colnames(data)
counter = 2 # for incrementing genes
counter2 = 1 # for incrementing rows in the dataframe
for (i in 1:length(mygenes)){
    geneA = mygenes[i]

    if (counter <= length(mygenes)){
        for (j in counter:length(mygenes)){
            if (geneA != mygenes[j]){
                geneB = mygenes[j]

                dat = as.data.frame(cbind(data[,i],data[,j]))
                colnames(dat) = c("geneA","geneB")

                logist1 = logistf(geneB ~ geneA,data=dat)

                mydf[counter2,1] = geneA
                mydf[counter2,2] = geneB
                mydf[counter2,3] = paste(geneA,geneB,sep="_")
                mydf[counter2,4] = as.numeric(coef(logist1)[2])
                mydf[counter2,5] = confint(logist1)[2,1]
                mydf[counter2,6] = confint(logist1)[2,2]
                mydf[counter2,7] = as.numeric(logist1$prob[2])
                counter2 = counter2 + 1

            }
        }
        counter = counter + 1
    }
}

mydf = as.data.frame(mydf)
colnames(mydf) = c("geneA","geneB","genePair","beta","loCI","hiCI","pvalue")
mydf[,4:7] = sapply(mydf[,4:7],as.numeric)
mydf = mydf[order(mydf$pvalue),]

write.table(mydf,"logistf.diagnosis.co-occurrence.txt",sep="\t",row.names=FALSE)

# CR1
data = read.table("../AML_data.muts.chroms.matrix.CR1.txt",sep="\t",header=TRUE,row.names=1)
data = data[,colnames(data) != "TPMT"]
data = data[,colSums(data) >= 2] # at least 3 mutations

mydf = matrix(nrow = choose(ncol(data),2),ncol=7)
mygenes = colnames(data)
counter = 2 # for incrementing genes
counter2 = 1 # for incrementing rows in the dataframe
for (i in 1:length(mygenes)){
    geneA = mygenes[i]

     if (counter <= length(mygenes)){
        for (j in counter:length(mygenes)){
            if (geneA != mygenes[j]){
                geneB = mygenes[j]

            dat = as.data.frame(cbind(data[,i],data[,j]))
            colnames(dat) = c("geneA","geneB")

            logist1 = logistf(geneB ~ geneA,data=dat)

            mydf[counter2,1] = geneA
            mydf[counter2,2] = geneB
            mydf[counter2,3] = paste(geneA,geneB,sep="_")
            mydf[counter2,4] = as.numeric(coef(logist1)[2])
            mydf[counter2,5] = confint(logist1)[2,1]
            mydf[counter2,6] = confint(logist1)[2,2]
            mydf[counter2,7] = as.numeric(logist1$prob[2])
            counter2 = counter2 + 1

            }
        }
        counter = counter + 1
     }
}

mydf = as.data.frame(mydf)
colnames(mydf) = c("geneA","geneB","genePair","beta","loCI","hiCI","pvalue")
mydf[,4:7] = sapply(mydf[,4:7],as.numeric)
mydf = mydf[order(mydf$pvalue),]

write.table(mydf,"logistf.CR1.co-occurrence.txt",sep="\t",row.names=FALSE)

# REL1
data = read.table("../AML_data.muts.chroms.matrix.relapse.txt",sep="\t",header=TRUE,row.names=1)
data = data[,colnames(data) != "TPMT"]
data = data[,colSums(data) >= 2] # at least 3 mutations

mydf = matrix(nrow = choose(ncol(data),2),ncol=7)
#mydf2 = matrix(nrow = choose(ncol(data),2),ncol=7)
mygenes = colnames(data)
counter = 2 # for incrementing genes
counter2 = 1 # for incrementing rows in the dataframe
for (i in 1:length(mygenes)){
    geneA = mygenes[i]

     if (counter <= length(mygenes)){
        for (j in counter:length(mygenes)){
            if (geneA != mygenes[j]){
                geneB = mygenes[j]

            dat = as.data.frame(cbind(data[,i],data[,j]))
            colnames(dat) = c("geneA","geneB")

            logist1 = logistf(geneB ~ geneA,data=dat)

            mydf[counter2,1] = geneA
            mydf[counter2,2] = geneB
            mydf[counter2,3] = paste(geneA,geneB,sep="_")
            mydf[counter2,4] = as.numeric(coef(logist1)[2])
            mydf[counter2,5] = confint(logist1)[2,1]
            mydf[counter2,6] = confint(logist1)[2,2]
            mydf[counter2,7] = as.numeric(logist1$prob[2])
            counter2 = counter2 + 1

            }
        }
        counter = counter + 1
     }
}

mydf = as.data.frame(mydf)
colnames(mydf) = c("geneA","geneB","genePair","beta","loCI","hiCI","pvalue")
mydf[,4:7] = sapply(mydf[,4:7],as.numeric)
mydf = mydf[order(mydf$pvalue),]
write.table(mydf,"logistf.REL1.co-occurrence.txt",sep="\t",row.names=FALSE)

#mydf2 = as.data.frame(mydf2)
#colnames(mydf2) = c("geneA","geneB","genePair","beta","loCI","hiCI","pvalue")
#mydf2[,4:7] = sapply(mydf2[,4:7],as.numeric)
#mydf2 = mydf2[order(mydf2$pvalue),]
#write.table(mydf2,"logit.REL1.co-occurrence.txt",sep="\t",row.names=FALSE)


#### plot select points
mypairs = c("ASXL1_FLT3","ASXL1_NPM1","DNMT3A_FLT3","DNMT3A_NPM1","FLT3_TET2","NPM1_TET2","ASXL1_CBL","ASXL1_SRSF2","CBL_DNMT3A","DNMT3A_SRSF2","CBL_TET2","SRSF2_TET2")

ddata = read.table("logistf.diagnosis.co-occurrence.txt",sep="\t",header=TRUE)
ddata = ddata[ddata$genePair %in% mypairs,]
ddata$stage = "diag"

cdata = read.table("logistf.CR1.co-occurrence.txt",sep="\t",header=TRUE)
cdata = cdata[cdata$genePair %in% mypairs,]
cdata$stage = "CR1"

rdata = read.table("logistf.REL1.co-occurrence.txt",sep="\t",header=TRUE)
rdata = rdata[rdata$genePair %in% mypairs,]
rdata$stage = "REL1"

finaldf = rbind(ddata,cdata,rdata)

library(ggpubr)
library(ggplot2)
finaldf$nlp = log10(finaldf$pvalue)*-1
finaldf[finaldf$nlp > 3,"nlp"] = 3

# flip columns

for (i in 1:nrow(finaldf)){
    if (finaldf[i,2] %in% c("ASXL1","DNMT3A","TET2")){
        holderA = finaldf[i,1]
        finaldf[i,1] = finaldf[i,2]
        finaldf[i,2] = holderA
    }
}
finaldf$geneA = factor(finaldf$geneA,levels=rev(c("DNMT3A","TET2","ASXL1")))
finaldf$geneB = factor(finaldf$geneB,levels=rev(c("SRSF2","CBL","NPM1","FLT3")))

finaldf[finaldf$beta > 2,"beta"] = 2
finaldf[finaldf$beta < -1.5,"beta"] = -1.5
finaldf$stage = factor(finaldf$stage,levels=(c("diag","CR1","REL1")))

finaldf = finaldf[finaldf$pvalue < 1,]

pdf("co-ocurrence-logistf.heatDot.pdf",height=2.75,width=5,useDingbats=FALSE)
ggplot(finaldf,aes(x=stage,y=geneA,group=geneB))+geom_point(aes(color=beta,size=nlp))+geom_text(data=finaldf[finaldf$pvalue<0.05,],aes(label="*"),size=8, hjust=0.5,vjust = 0.75,color="black") + facet_grid(. ~ geneB) + scale_color_gradientn(colors=c("#1F7452","#60C99F","white","#F28FC8","#7B257A"),values= scales::rescale(c(-1.5,-0.75,0,1,2)),limits=c(-1.5,2) ) + scale_size(limits=c(0,3),range=c(1,9)) + scale_y_discrete(position="left") + theme_bw()+ theme(legend.position="bottom")
dev.off()
