#!/usr/bin/env Rscript

require(data.table,quietly = TRUE)
require(optparse,quietly = TRUE)
require(ggplot2,quietly = TRUE)

message("------------")
usage=paste(getopt::get_Rscript_filename(), "file")
message(paste("Usage: ", usage))

# Setup -------------------------------------------------------------------
opt_parser = OptionParser(usage = usage);
args <- parse_args(opt_parser, positional_arguments = 1)
kpnn.output.folder <- args$args
if(!dir.exists(kpnn.output.folder)) message("Missing folder of KPNN output:", kpnn.output.folder, " please provide an existing folder")
print(paste("Processing KPNN output from folder:", kpnn.output.folder))
xRot <- function(){theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))}

# Read in Data ------------------------------------------------------------
for(fx in c("Memory.tsv", "tf_cost.csv", "tf_NumGradMeans.csv", "tf_NumGradTestError.txt")){
  if(!fx %in% list.files(kpnn.output.folder)) message("Missing file ", fx, " in ", kpnn.output.folder)
}
mem <- fread(paste0(kpnn.output.folder, "Memory.tsv"))
cost <- fread(paste0(kpnn.output.folder, "tf_cost.csv"))
nodeWeights <- fread(paste0(kpnn.output.folder, "tf_NumGradMeans.csv"))
test.error <- fread(paste0(kpnn.output.folder, "tf_NumGradTestError.txt"))


# Training curve ----------------------------------------------------------
difftime <- as.POSIXct(mem[nrow(mem)]$Time) - as.POSIXct(mem[1]$Time)
title <- paste("Test error on unseen data:", round(test.error$NumGradTestError,5), "\nTime:", round(as.numeric(difftime),5), attr(difftime, "units"))
pDat <- melt.data.table(cost, measure.vars = c("train", "validation"))
ggplot(pDat, aes(x=iteration, y=value, color=variable)) + geom_line() +
  scale_color_discrete(name="Data") + xlab("Iteration") + ylab("Loss") +
  theme_bw(12) +
  ggtitle(title)
ggsave(paste0(kpnn.output.folder, "TrainingCurve.pdf"), w=5,h=4)



# Node weights ------------------------------------------------------------
nodeWeights[,NodeWeight := abs(output)]
if("nodeX" %in% nodeWeights$Node) nodeWeights[,Node := gsub("nodeX", "A", gsub("b", "", Node))]
nodeWeights$Node <- factor(nodeWeights$Node, levels = nodeWeights[order(NodeWeight)]$Node)
pDat <- if(nrow(nodeWeights) <= 40) nodeWeights else nodeWeights[order(NodeWeight)][c(1:20,(nrow(nodeWeights)-20):nrow(nodeWeights))]
ggplot(pDat, aes(x=Node, y=NodeWeight)) + geom_point(color="red") + # geom_bar(stat="identity") + 
  theme_classic(12) + xRot() + ggtitle("Top and bottom 20 nodes shown")
ggsave(paste0(kpnn.output.folder, "NodeWeights.pdf"), w=8, h=4)

message("Successfully summarized KPNN run in folder ", kpnn.output.folder ,"\n")