#!/usr/bin/Rscript

library(clevRvis)
library(viridis)
library(RColorBrewer)
library(data.table)

myfiles = read.table("dotlist.txt",sep="\t",header=FALSE)
mypatients = sapply(strsplit(myfiles[,1], "\\."), function(x) x[2])

for (mypatient in mypatients){
  tryCatch({
#mypatient = "100"
    mytable = read.table(paste("./output-res/pt.",mypatient,".soln1.clevRvis.txt",sep=""),sep="\t",header=TRUE)
    mytable = as.data.frame(t(mytable))
    if (ncol(mytable)==3){
        timepoints = c(0,50,100)
        if (mypatient == "132"){
            timelabs = c("Diagnosis","REL1","REF1")
        }
        else {
            timelabs = c("Diagnosis","CR1","REL1")
        }
        #colnames(mytable) = mytable[1,1:3]
    } else {
        timepoints = c(0,100)
        timelabs = c("Diagnosis","REL1")
        #colnames(mytable) = mytable[1,1:2]
    }

    mytable = mytable[-1,]
    cloneLabels = rownames(mytable)
    mymat = matrix(unlist(sapply(mytable,as.numeric))*100,ncol=ncol(mytable))

    parents = read.table(paste("./output-res/pt.",mypatient,".soln1.clevRvis.indices.txt",sep=""),sep="\t",header=FALSE)
    parents = as.numeric(parents$V2)

    # reorder the matrix
    testdf = as.data.frame(cbind(mymat,parents))
    rownames(testdf) = cloneLabels
    testdf$parentsName = 0
    for (i in 1:nrow(testdf)){
        if (testdf[i,"parents"] != 0){
            testdf[i,"parentsName"] = rownames(testdf)[parents[i]]
        }
    }

    testdf2 = testdf[order(testdf$V1,decreasing=TRUE),]
    testdf2$rownum = seq(1:nrow(testdf2))
    testdf2$cloneName = rownames(testdf2)
    for (i in 1:nrow(testdf2)){
        if (testdf2[i,"parentsName"] != 0){
            testdf2[i,"parents"] = testdf2[rownames(testdf2) == testdf2$parentsName[i],"rownum"]
        }
    }

    mytable2 = testdf2[,c(1:(ncol(testdf2)-4))]

    mymat2 = matrix(unlist(sapply(mytable2,as.numeric)),ncol=ncol(mytable2))
    parents2 = testdf2$parents

    mycols =  brewer.pal(nrow(mymat2),"YlGnBu")
    if (nrow(mymat2) < 3){ mycols = mycols[1:2]}

    seaObject_tp <- createSeaObject(mymat2, parents2, timepoints, cloneLabels=rownames(mytable2), timepointInterpolation = TRUE,col = mycols)

    pdf(paste("./clevRvis-plots/dolphinPlot.pt.",mypatient,".pdf",sep=""),height=4,width=8,useDingbats=FALSE)
    print(dolphinPlot(seaObject_tp, showLegend = TRUE, vlines = timepoints, vlab = timelabs, vlabSize = 2,ylab = 'Cancer cell fractions (CCFs)',markMeasuredTimepoints = timepoints, separateIndependentClones = TRUE))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#sharkPlot(seaObject_tp, showLegend = TRUE, main = 'Example Shark plot')