#!/usr/bin/Rscript
# compare chr8gain vs no-chr8gain mut freqs at diagnosis

data=read.table("../fig4/AML_data.muts.chroms.matrix.diagnosis.txt",sep="\t",header=TRUE,row.names=1)

library(stringr)
data$patient = str_split_i(rownames(data),":",2)

# keep the first entry for a given patient
library(dplyr)
dataf <- data %>%
  group_by(patient) %>%
  slice(1) %>%
  ungroup()
  
ctsumm <- dataf %>%
  group_by(chr8gain) %>%
  summarize(across(where(is.numeric), ~ sum(. == 1), .names = "{.col}")) %>%
  ungroup()

cts = t(as.data.frame(ctsumm))

colnames(cts) = c("no_chr8gain.mut","chr8gain.mut")
cts = cts[-1,]
cts = as.data.frame(cts)
cts[,1] = as.numeric(cts[,1])
cts[,2] = as.numeric(cts[,2])
cts$gene = rownames(cts)
cts$chr8gain.wt = table(dataf$chr8gain)[2] - cts$chr8gain.mut
cts$no_chr8gain.wt = table(dataf$chr8gain)[1] - cts$no_chr8gain.mut

mypvals = ""
for (i in 1:nrow(cts)){
    ddf = rbind(c(cts$chr8gain.mut[i],cts$chr8gain.wt [i]), c(cts$no_chr8gain.mut[i],cts$no_chr8gain.wt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}




summary <- dataf %>%
  select(chr8gain, everything()) %>%
  group_by(chr8gain) %>%
  summarize(across(where(is.numeric), ~ mean(. == 1) * 100, .names = "{.col}")) %>%
  ungroup()

mydf = t(as.data.frame(summary))

colnames(mydf) = c("no_chr8gain.mut","chr8gain.mut")
mydf = mydf[-1,]
mydf = as.data.frame(mydf)
mydf[,1] = as.numeric(mydf[,1])
mydf[,2] = as.numeric(mydf[,2])
mydf$gene = rownames(mydf)

mydf$pvalue = as.numeric(mypvals)
mydf$fisher.nlp = log10(mydf$pvalue)*-1

# keep genes mutated in >= 5% of the entire cohort
#mydf2 = mydf[rownames(mydf) %in% colnames(data[,colSums(data[,-ncol(data)])/nrow(data)>=0.05]),]
mydf2 = mydf[mydf$chr8gain.mut >= 10 | mydf$no_chr8gain.mut >= 10,]

library(ggpubr)
library(viridis)

mydf2$no_chr8gain.mutneg = mydf2$no_chr8gain.mut*-1
mydf2$gene = factor(mydf2$gene,levels=mydf2[order(mydf2$no_chr8gain.mut),"gene"])

p1=ggbarplot(mydf2,x="gene",y="chr8gain.mut",fill="dodgerblue3",orientation="horizontal",label=round(mydf2$pvalue,digits=2))+ylim(c(0,50))+geom_hline(yintercept = c(20,40), linetype = "dashed", color = "gray55")

p2=ggbarplot(mydf2,x="gene",y="no_chr8gain.mutneg",fill="goldenrod",orientation="horizontal")+ylim(c(-50,0))+geom_hline(yintercept = c(-20,-40), linetype = "dashed", color = "gray55")

pdf("diagnosis.chr8gain-vs-nochr8gain.pop-mutfreq.barplot.pdf",height=6,width=7.5,useDingbats=FALSE)
ggarrange(p2,p1,ncol=2)
dev.off()

######
####
# now relapse1
# compare chr8gain vs no-chr8gain mut freqs at relapse

data=read.table("../fig4/AML_data.muts.chroms.matrix.relapse.txt",sep="\t",header=TRUE,row.names=1)

library(stringr)
data$patient = str_split_i(rownames(data),":",2)

# keep the first entry for a given patient
library(dplyr)
dataf <- data %>%
  group_by(patient) %>%
  slice(1) %>%
  ungroup()
  
ctsumm <- dataf %>%
  group_by(chr8gain) %>%
  summarize(across(where(is.numeric), ~ sum(. == 1), .names = "{.col}")) %>%
  ungroup()

cts = t(as.data.frame(ctsumm))

colnames(cts) = c("no_chr8gain.mut","chr8gain.mut")
cts = cts[-1,]
cts = as.data.frame(cts)
cts[,1] = as.numeric(cts[,1])
cts[,2] = as.numeric(cts[,2])
cts$gene = rownames(cts)
cts$no_chr8gain.wt = table(dataf$chr8gain)[1] - cts$no_chr8gain.mut
cts$chr8gain.wt = table(dataf$chr8gain)[2] - cts$chr8gain.mut


mypvals = ""
for (i in 1:nrow(cts)){
    ddf = rbind(c(cts$chr8gain.mut[i],cts$chr8gain.wt [i]), c(cts$no_chr8gain.mut[i],cts$no_chr8gain.wt[i]))
    mypvals[i] = fisher.test(ddf)$p.value
}




summary <- dataf %>%
  select(chr8gain, everything()) %>%
  group_by(chr8gain) %>%
  summarize(across(where(is.numeric), ~ mean(. == 1) * 100, .names = "{.col}")) %>%
  ungroup()

mydf = t(as.data.frame(summary))

colnames(mydf) = c("no_chr8gain.mut","chr8gain.mut")
mydf = mydf[-1,]
mydf = as.data.frame(mydf)
mydf[,1] = as.numeric(mydf[,1])
mydf[,2] = as.numeric(mydf[,2])
mydf$gene = rownames(mydf)

mydf$pvalue = as.numeric(mypvals)
mydf$fisher.nlp = log10(mydf$pvalue)*-1

# keep genes mutated in >= 5% of the entire cohort
#mydf2 = mydf[rownames(mydf) %in% colnames(data[,colSums(data[,-ncol(data)])/nrow(data)>=0.05]),]
mydf2 = mydf[mydf$chr8gain.mut >= 10 | mydf$no_chr8gain.mut >= 10,]

library(ggpubr)
library(viridis)

mydf2$no_chr8gain.mutneg = mydf2$no_chr8gain.mut*-1
mydf2$gene = factor(mydf2$gene,levels=mydf2[order(mydf2$no_chr8gain.mut),"gene"])

p1=ggbarplot(mydf2,x="gene",y="chr8gain.mut",fill="dodgerblue3",orientation="horizontal",label=round(mydf2$pvalue,digits=2))+ylim(c(0,50))+geom_hline(yintercept = c(20,40), linetype = "dashed", color = "gray55")

p2=ggbarplot(mydf2,x="gene",y="no_chr8gain.mutneg",fill="goldenrod",orientation="horizontal")+ylim(c(-50,0))+geom_hline(yintercept = c(-20,-40), linetype = "dashed", color = "gray55")

pdf("relapse.chr8gain-vs-nochr8gain.pop-mutfreq.barplot.pdf",height=6,width=7.5,useDingbats=FALSE)
ggarrange(p2,p1,ncol=2)
dev.off()