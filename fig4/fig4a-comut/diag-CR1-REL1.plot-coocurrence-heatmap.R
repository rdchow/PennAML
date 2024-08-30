#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(viridis)

zcutoff = 4
freqcutoff = 0

# filter to genes that are mutated in >3.5% of diagnosis or REL1 samples
mutfreqs = read.table("../../fig1/mut-frequencies.txt",skip = 1,header=TRUE)
mygenes = mutfreqs[(mutfreqs$DiagnosisMutFreq >= freqcutoff) | (mutfreqs$Rel1MutFreq >= freqcutoff) ,"gene"]
mygenes = mygenes[mygenes != "TPMT"]

mygenesd = mutfreqs[(mutfreqs$DiagnosisMutFreq >= freqcutoff) ,"gene"]
mygenesd = mygenesd[mygenesd != "TPMT"]

mygenesc = mutfreqs[(mutfreqs$CR1mutFreq >= freqcutoff) ,"gene"]
mygenesc = mygenesc[mygenesc != "TPMT"]

mygenesr = mutfreqs[(mutfreqs$Rel1MutFreq >= freqcutoff) ,"gene"]
mygenesr = mygenesr[mygenesr != "TPMT"]


### diagnosis ######
# mutual exclusivity count matrices taken from cbioportal, then recalculated here for more precise statistics
ddata=read.table("mutual-exclusivity-table.diagnosis.2.tsv",sep="\t",header=TRUE)
ddata = ddata[ddata$A %in% mygenes & ddata$B %in% mygenesd,]

ddata = ddata[,c(-7:-10)]
ddata$or = ""
ddata$pval = ""

for (i in 1:nrow(ddata)){
    mydat = ddata[i,]
    mydf = data.frame(rbind(c(as.numeric(mydat[3]),as.numeric(mydat[4])),c(as.numeric(mydat[5]),as.numeric(mydat[6]))))
    myfish = fisher.test(mydf)
    ddata[i,"or"] = myfish$estimate
    ddata[i,"pval"] = myfish$p.value
}
ddata$pval = as.numeric(ddata$pval)
ddata$or = as.numeric(ddata$or)
ddata$logor = log(ddata$or)
ddata$nlp = log10(ddata$pval)*-1
ddata$adjp = p.adjust(ddata$pval,method="BH")
write.table(ddata,"mutual-exclusivity-table.diagnosis.reanalyze.txt",sep="\t",row.names=FALSE)

# trim the co-occurrence mat range
ddata$logor[ddata$logor < -zcutoff] = -zcutoff # capped at +/- 2 for visualization
ddata$logor[ddata$logor > zcutoff] = zcutoff
ddata$logor = as.numeric(ddata$logor)
ddata$nlp = as.numeric(ddata$nlp)
ddata$nlp[ddata$nlp > 3] = 3
ddata$A = factor(ddata$A,levels=rev(sort(unique(ddata$A))))
ddata$B = factor(ddata$B,levels=rev(sort(unique(ddata$B))))
#ddata$nlp[ddata$pval > 0.05] = NA

### CR1 ########
###
cdata=read.table("mutual-exclusivity-table.CR1.2.tsv",sep="\t",header=TRUE)
cdata = cdata[cdata$A %in% mygenes & cdata$B %in% mygenesc,]
cdata = cdata[,c(-7:-10)]
cdata$or = ""
cdata$pval = ""
for (i in 1:nrow(cdata)){
    mydat = cdata[i,]
    mydf = data.frame(rbind(c(as.numeric(mydat[3]),as.numeric(mydat[4])),c(as.numeric(mydat[5]),as.numeric(mydat[6]))))
    if (complete.cases(mydat[3:6])){
        myfish = fisher.test(mydf)
        cdata[i,"or"] = myfish$estimate
        cdata[i,"pval"] = myfish$p.value
    }
    else {
        cdata[i,"or"] = NA
        cdata[i,"pval"] = NA
    }
}
cdata$pval = as.numeric(cdata$pval)
cdata$or = as.numeric(cdata$or)
cdata$logor = log(cdata$or)
cdata$nlp = log10(cdata$pval)*-1
cdata$adjp = p.adjust(cdata$pval,method="BH")
write.table(cdata,"mutual-exclusivity-table.CR1.reanalyze.txt",sep="\t",row.names=FALSE)

# trim the co-occurrence mat range
cdata$logor[cdata$logor < -zcutoff] = -zcutoff
cdata$logor[cdata$logor > zcutoff] = zcutoff
cdata$logor = as.numeric(cdata$logor)
cdata$nlp = as.numeric(cdata$nlp)
cdata$nlp[cdata$nlp > 3] = 3
cdata$A = factor(cdata$A,levels=rev(sort(unique(cdata$A))))
cdata$B = factor(cdata$B,levels=rev(sort(unique(cdata$B))))
#cdata$nlp[cdata$pval > 0.05] = NA

### relapse ########
###
rdata=read.table("mutual-exclusivity-table.REL1.2.tsv",sep="\t",header=TRUE)
rdata = rdata[rdata$A %in% mygenes & rdata$B %in% mygenesr,]
rdata = rdata[,c(-7:-10)]
rdata$or = ""
rdata$pval = ""
for (i in 1:nrow(rdata)){
    mydat = rdata[i,]
    mydf = data.frame(rbind(c(as.numeric(mydat[3]),as.numeric(mydat[4])),c(as.numeric(mydat[5]),as.numeric(mydat[6]))))
    if (complete.cases(mydat[3:6])){
        myfish = fisher.test(mydf)
        rdata[i,"or"] = myfish$estimate
        rdata[i,"pval"] = myfish$p.value
    }
    else {
        rdata[i,"or"] = NA
        rdata[i,"pval"] = NA
    }
}
rdata$pval = as.numeric(rdata$pval)
rdata$or = as.numeric(rdata$or)
rdata$logor = log(rdata$or)
rdata$nlp = log10(rdata$pval)*-1
rdata$adjp = p.adjust(rdata$pval,method="BH")
write.table(rdata,"mutual-exclusivity-table.REL1.reanalyze.txt",sep="\t",row.names=FALSE)

# trim the co-occurrence mat range
rdata$logor[rdata$logor < -zcutoff] = -zcutoff
rdata$logor[rdata$logor > zcutoff] = zcutoff
rdata$logor = as.numeric(rdata$logor)
rdata$nlp = as.numeric(rdata$nlp)
rdata$nlp[rdata$nlp > 3] = 3
rdata$A = factor(rdata$A,levels=rev(sort(unique(rdata$A))))
rdata$B = factor(rdata$B,levels=rev(sort(unique(rdata$B))))
#rdata$nlp[rdata$pval > 0.05] = NA


ddata$pair = paste(ddata$A, ddata$B,sep="_")
cdata$pair = paste(cdata$A, cdata$B,sep="_")
rdata$pair = paste(rdata$A, rdata$B,sep="_")

# subset to specific pairs of interest
mypairs = c("DNMT3A_FLT3","DNMT3A_NPM1","FLT3_TET2","NPM1_TET2","CBL_DNMT3A","DNMT3A_SRSF2","CBL_TET2","SRSF2_TET2")

ddata2 = ddata[ddata$pair %in% mypairs,]

#######

##### calculate the overlaps ####
sigddata = ddata[ddata$pval < 0.05 & !is.na(ddata$pval),]
sigcdata = cdata[cdata$pval < 0.05 & !is.na(cdata$pval),]
sigrdata = rdata[rdata$pval < 0.05 & !is.na(rdata$pval),]

drsharedpairs = sigddata$pair[sigddata$pair %in% sigrdata$pair] # signif in diag + rel1
dcsharedpairs = sigddata$pair[sigddata$pair %in% sigcdata$pair] # signif in diag + cr1
dcrsharedpairs = drsharedpairs[drsharedpairs %in% dcsharedpairs] # signif in diag + cr1 + rel1

allsigpairs = unique(c(sigddata$pair,sigcdata$pair,sigrdata$pair)) # signif in diag, cr1, or rel1
drsigpairs = unique(c(sigddata$pair,sigrdata$pair))# signif in diag or rel1

# filter the data matrices to the ones that were significant in diag or REL1
library(stringr)
ddata2 = ddata[ddata$pair %in% drsigpairs,]
cdata2 = cdata[cdata$pair %in% drsigpairs,]
rdata2 = rdata[rdata$pair %in% drsigpairs,]

all(ddata2$pair %in% cdata2$pair)
all(ddata2$pair %in% rdata2$pair)
all(cdata2$pair %in% rdata2$pair)

# refactor genes
ddata2$A = factor(ddata2$A,levels=rev(sort(unique(ddata2$A))))
ddata2$B = factor(ddata2$B,levels=rev(sort(unique(ddata2$B))))
cdata2$A = factor(cdata2$A,levels=rev(sort(unique(cdata2$A))))
cdata2$B = factor(cdata2$B,levels=rev(sort(unique(cdata2$B))))
rdata2$A = factor(rdata2$A,levels=rev(sort(unique(rdata2$A))))
rdata2$B = factor(rdata2$B,levels=rev(sort(unique(rdata2$B))))

# then make "NA" pairs that have pval == 1
ddata2[ddata2$pval == 1 & !is.na(ddata2$pval),"nlp"] = NA
cdata2[cdata2$pval == 1 & !is.na(cdata2$pval),"nlp"] = NA
rdata2[rdata2$pval == 1 & !is.na(rdata2$pval),"nlp"] = NA

ddata2[ddata2$pval == 1 & !is.na(ddata2$pval),"logor"] = NA
cdata2[cdata2$pval == 1 & !is.na(cdata2$pval),"logor"] = NA
rdata2[rdata2$pval == 1 & !is.na(rdata2$pval),"logor"] = NA

# plots
library(ggpubr)

g1 = ggplot(ddata2,aes(x=B,y=A)) + geom_point(data = ddata2[ddata2$pval < 0.05,], aes(size=nlp,fill=logor),pch=21,stroke=1.3,color="black") + scale_fill_distiller(palette="RdBu",limits=c(-zcutoff,zcutoff)) + theme_pubr() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) + scale_size(limits=c(0,3),range=c(1,5)) + grids(linetype = "dashed") + scale_x_discrete(position = "top",limits=rev(levels(ddata2$B)), drop=FALSE) + scale_y_discrete(limits=rev(levels(ddata2$A)),drop=FALSE) +
geom_tile(data=ddata2[ddata2$pair %in% drsharedpairs,],fill=NA,color="black",linewidth=0.8,width=1.1,height=1.1) + 
geom_tile(data=ddata2[ddata2$pair %in% dcsharedpairs,],fill=NA,color="purple3",linewidth=0.8,width=1.1,height=1.1) + 
geom_tile(data=ddata2[ddata2$pair %in% dcrsharedpairs,],fill=NA,color="darkgreen",linewidth=0.8,width=1.1,height=1.1) + coord_equal()


g2 = ggplot(cdata2,aes(x=B,y=A)) + geom_point(data = cdata2[cdata2$pval < 0.05,], aes(size=nlp,fill=logor),pch=21,stroke=1.3,color="black") + scale_fill_distiller(palette="RdBu",limits=c(-zcutoff,zcutoff)) + theme_pubr() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) + scale_size(limits=c(0,3),range=c(1,5)) + grids(linetype = "dashed") + scale_x_discrete(position = "top",limits=rev(levels(cdata2$B)), drop=FALSE) + scale_y_discrete(limits=rev(levels(cdata2$A)),drop=FALSE) +
geom_tile(data=cdata2[cdata2$pair %in% dcsharedpairs,],fill=NA,color="purple3",linewidth=0.8,width=1.1,height=1.1) + 
geom_tile(data=cdata2[cdata2$pair %in% dcrsharedpairs,],fill=NA,color="darkgreen",linewidth=0.8,width=1.1,height=1.1) + coord_equal()


g3 = ggplot(rdata2,aes(x=B,y=A)) + geom_point(data = rdata2[rdata2$pval < 0.05,], aes(size=nlp,fill=logor),pch=21,stroke=1.3,color="black") + scale_fill_distiller(palette="RdBu",limits=c(-zcutoff,zcutoff)) + theme_pubr() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) + scale_size(limits=c(0,3),range=c(1,5)) + grids(linetype = "dashed") + scale_x_discrete(position = "top",limits=rev(levels(rdata2$B)), drop=FALSE) + scale_y_discrete(limits=rev(levels(rdata2$A)),drop=FALSE) +
geom_tile(data=rdata2[rdata2$pair %in% drsharedpairs,],fill=NA,color="black",linewidth=0.8,width=1.1,height=1.1) + 
geom_tile(data=rdata2[rdata2$pair %in% dcrsharedpairs,],fill=NA,color="darkgreen",linewidth=0.8,width=1.1,height=1.1) + coord_equal()

library(gridExtra)
pdf("mutual-exclusivity-dotplot.diag.cr1.rel1.pdf",height=19,width=16,useDingbats=FALSE,onefile=TRUE)
grid.arrange(g1,g2,g3,ncol=2)
dev.off()