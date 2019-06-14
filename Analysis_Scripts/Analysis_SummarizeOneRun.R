source(paste0(Sys.getenv("KPNN_CODEBASE"), "/init.R"))


# Setup -------------------------------------------------------------------
args=commandArgs(trailingOnly = TRUE)
kpnn.output.folder <- if(length(args) == 0) paste0(Sys.getenv("KPNN_OUTPUTS"), "/run_2/") else args[1]
if(!dir.exists(kpnn.output.folder)) message("Missing folder of KPNN output:", kpnn.output.folder, " please provide the correct folder")
out <- "Demo_OneNetwork_Analysis/"
dir.create(dirout(out))
print(paste("Processing KPNN output from folder:", kpnn.output.folder))

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
title <- paste("Test error on unseen data:", round(test.error$NumGradTestError,5), "\nTime:", as.numeric(difftime), attr(difftime, "units"))
pDat <- melt.data.table(cost, measure.vars = c("train", "validation"))
ggplot(pDat, aes(x=iteration, y=value, color=variable)) + geom_line() +
  scale_color_discrete(name="Data") + xlab("Iteration") + ylab("Loss") +
  theme_bw(12) +
  ggtitle(title)
ggsave(dirout(out, "TrainingCurve.pdf"), w=5,h=4)



# Node weights ------------------------------------------------------------
nodeWeights[,NodeWeight := abs(output)]
if("nodeX" %in% nodeWeights$Node) nodeWeights[,Node := gsub("nodeX", "A", gsub("b", "", Node))]
nodeWeights$Node <- factor(nodeWeights$Node, levels = nodeWeights[order(NodeWeight)]$Node)

ggplot(nodeWeights, aes(x=Node, y=NodeWeight)) + geom_bar(stat="identity") + theme_bw(12) + xRot() 
ggsave(dirout(out, "NodeWeights.pdf"), w=15, h=4)
