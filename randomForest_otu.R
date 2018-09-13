setwd("~/Desktop/BMD_turkeys/data/16s_amplicon/RDP16mothur_update/")
#load data and normalize it
library(randomForest)
library(dplyr)

#load data and normalize it, eventually the samples should be in the rows and the OTUs in the columns
#OTU table should be subsampled to equal number per sample.
data <- read.table("stability.trim.opti_mcc.0.03.subsample.txt", header = T)
data2 <- data[,c(-1,-3)]
#design file
design_file <- read.table("BMD_design.txt", header = T)
design_file_trt <- design_file[,-2:-6]
#taxonomy file
taxonomy <- read.table(file ="stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.abund.opti_mcc.0.03.cons.taxonomy", header =T)

OTU_abund <- as.data.frame(colSums(data[,-1:-2])/nrow(data))
taxonomy <- merge(taxonomy, OTU_abund, by.x = "OTU", by.y = 0)

data2_trt <- merge(design_file_trt, data2, by.x = "name", by.y = "Group")

#So now we have a data fram with column 1 as sample name, column 2 as trt, and the rest of the columns as OTUs

#Now we want data frames with just the two treatment groups we want.
data <- data2_trt
#data_sub <- droplevels(as.data.frame(subset(data, trt!="ther" )))
#data_ther <- droplevels(as.data.frame(subset(data, trt!="sub" )))
early_sub <- droplevels(as.data.frame(subset(data, trt_window=="early_ctrl"|trt_window=="early_sub" )))
early_ther <- droplevels(as.data.frame(subset(data, trt_window=="early_ctrl"|trt_window=="early_ther" )))
later_sub <- droplevels(as.data.frame(subset(data, trt_window=="later_ctrl"|trt_window=="later_sub" )))
later_ther <- droplevels(as.data.frame(subset(data, trt_window=="later_ctrl"|trt_window=="later_ther" )))
withdraw_sub <- droplevels(as.data.frame(subset(data, trt_window=="withdraw_ctrl"|trt_window=="withdraw_sub" )))
withdraw_ther <- droplevels(as.data.frame(subset(data, trt_window=="withdraw_ctrl"|trt_window=="withdraw_ther" )))


comparisons <- c("early_sub", "early_ther", "later_sub", "later_ther", "withdraw_sub", "withdraw_ther")

#Now run the randomForest algorythm. x = the otu table, y = the response vector, or treatment groups.
set.seed(1)
for(i in comparisons){
  print(i)
  rf <- randomForest(x = get(i)[, -1:-2], y = get(i)[, 2], importance=TRUE, proximity=TRUE, ntree=5000)
  
  
  png(paste0(i,"_OTU_rf_figure1.png"), height=1800, width=1800, res = 300)
  par(mfrow=c(2,1))
  par(pty="s")
  varImpPlot(rf, type=1, pch=19, col=1, cex=.5, main="")
  varImpPlot(rf, type=2, pch=19, col=1, cex=.5, main="")
  dev.off()
  
  fo <- file(paste0(i, "_OTU_rf.txt"), "w")
  imp <- importance(rf)
  imp <- merge(imp, taxonomy, by.x = 0, by.y = "OTU")
  imp_sort <- arrange(imp, desc(MeanDecreaseAccuracy))
  write.table(imp_sort, fo, sep="\t", row.names = F)
  flush(fo)
  close(fo)
}
