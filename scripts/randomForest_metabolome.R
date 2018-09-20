setwd(dir = "~/PATH/TO/BMD_TURKEY")
#load data and normalize it
library(randomForest)
library(dplyr)
library(openxlsx)

#load data and normalize it
data <- read.csv("data/metabolome_cc.csv");
data2 <- data[,c(-1,-3:-9)]
metabolome_key <- data[,c(2,3,1,4,6,8)]
data2_t <- t(data2[,-1])
column_names <- (data2[,1])
colnames(data2_t) <- column_names
turkey_info <- read.table("data/BMD_design_cc.txt", header = TRUE)
turkey_info <- filter(turkey_info, location =="cc")
str(turkey_info)
turkey_info <- droplevels(turkey_info)
turkey_info <- turkey_info[match(row.names(data2_t), turkey_info$turkey),]

data2_t <- merge(turkey_info, data2_t, by.x = "turkey", by.y = 0)
data <- data2_t
str(data)
days <- unique(turkey_info$day)
trts <- c("sub", "ther")
#d<-7
#t<-"sub"
rm(rf_summary)
for(d in days){
  print(d)
  day_ctrl <- filter(data, day==d) %>%
    filter(trt=="ctrl")
  for(t in trts){
    day_trt <- filter(data, day==d) %>%
      filter(trt==t)
    data_comp <- rbind(day_ctrl, day_trt)
    data_comp1 <- data_comp[,c(-1:-3, -5:-9)]
    data_comp1$group <- factor(data_comp1$group)
    row.names(data_comp1) <- data_comp$turkey
    rf_temp <- randomForest(x = data_comp1[, 2:ncol(data_comp1)], y = data_comp1[, 1], importance=TRUE, proximity=TRUE, ntree=5000)
    
    png(paste0("output/d", d, t, "_metabolome_rf.png"), height=1800, width=1800, res = 300)
    par(mfrow=c(2,1))
    par(pty="s")
    varImpPlot(rf_temp, type=1, pch=19, col=1, cex=.5, main="")
    varImpPlot(rf_temp, type=2, pch=19, col=1, cex=.5, main="")
    dev.off()
    
    fo <- file(paste0("output/d", d, t, "_metabolome_rf.txt"), "w")
    imp <- as.data.frame(importance(rf_temp))
    write.table(imp, fo, sep="\t")
    flush(fo)
    close(fo)
    
    colnames(imp)[1:2] <- c("ctrl", "trt")
    imp$biochemical <- row.names(imp)
    imp$comparison <- t
    imp$day <- d
    imp <- arrange(imp, desc(MeanDecreaseAccuracy))
    imp_filtered <- filter(imp, MeanDecreaseAccuracy > 6)

    if(!exists("rf_summary")){
      rf_summary <- imp_filtered
    } else {
      rf_summary <- rbind(rf_summary, imp_filtered)
    }
  }
}

str(rf_summary$biochemical)

molecules <- read.xlsx("data/molecules_key.xlsx")
str(molecules)
rf_summary_key <- merge(molecules, rf_summary, by.x = "BIOCHEMICAL", by.y = "biochemical")
rf_summary_key$SUPER.PATHWAY <- factor(rf_summary_key$SUPER.PATHWAY)
rf_summary_key$BIOCHEMICAL <- factor(rf_summary_key$BIOCHEMICAL, levels = rf_summary_key$BIOCHEMICAL[order((rf_summary_key$day))], ordered = TRUE)
str(rf_summary_key)

#ggplot(complete(rf_summary_key, BIOCHEMICAL, comparison, day), aes(x = reorder(BIOCHEMICAL, -(day)), y = MeanDecreaseAccuracy, fill = comparison)) +
my_colors <- c(
  '#1f78b4', '#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(complete(rf_summary_key, BIOCHEMICAL, comparison, day), aes(x = BIOCHEMICAL, y = MeanDecreaseAccuracy, fill = comparison)) +
  facet_grid(.~day) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  #theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(c(0,11)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)) +
  coord_flip()

ggsave("output/rf_summary.png", height = 7, width = 6)

