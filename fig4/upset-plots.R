
#!/usr/bin/Rscript
# use TET2 vs DNMT3A definitions at time of diagnosis to guide subsequent labels
# NOT from each stage
library(ggplot2)
library(ggpubr)
library(UpSetR)

# binary mutation matrix is from "matrix-muts.R"
data = read.table("AML_data.muts.chroms.matrix.txt",sep="\t",header=TRUE,row.names=1)

info = read.table("../fig2/uniqueStage-sample-names.txt",sep="\t")
colnames(info) = c("sampleID","patient","stage")

infoDN = info[info$stage=="dn",]
dataDN = data[rownames(data) %in% infoDN$sampleID,]
infoDN = infoDN[match(rownames(dataDN),infoDN$sampleID),]
all(rownames(dataDN) == infoDN$sampleID)

dataDN$patient = infoDN$patient
dataDN$stage = infoDN$stage

pdf("upset-plot.DTA.diagnosis.pdf",height=5,width=8,useDingbats=FALSE)
upset(dataDN,sets=c("ASXL1","DNMT3A","TET2"),empty.intersections= TRUE,point.size=5, mb.ratio=c(0.65,0.35))
dev.off()