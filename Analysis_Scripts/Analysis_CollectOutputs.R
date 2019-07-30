#!/usr/bin/env Rscript

require(data.table,quietly = TRUE)
require(optparse,quietly = TRUE)
require(ggplot2,quietly = TRUE)

message("------------")
usage=paste(getopt::get_Rscript_filename(), "inputFolder outputFolder")
message(paste("Usage: ", usage))


# Setup -------------------------------------------------------------------
opt_parser = OptionParser(usage = usage);
args <- parse_args(opt_parser, positional_arguments = 2)
input.folder <- args$args[1]
# input.folder <- "~/projects_shared/pathway_learning/results_analysis/46_03_TCR_opt/"
output.folder <- args$args[2]
# output.folder <- "~/projects_shared/pathway_learning/results_analysis/46_03_TCR_opt/"
if(!dir.exists(input.folder)) message("Missing folder of trained KPNNs :", input.folder, " please provide an existing folder")
print(paste("Processing trained KPNNs from folder:", input.folder))
dir.create(output.folder,recursive = TRUE)


# FUNCTIONS ---------------------------------------------------------------
xRot <- function(){theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))}
write.tsv <- function(...){ write.table(..., sep="\t", row.names=FALSE, quote=FALSE);}
toMT <- function(dt, row, col, val){
  retDT <- dcast.data.table(dt, get(row) ~ get(col), value.var=val)
  retMT <- as.matrix(retDT[,-"row"])
  row.names(retMT) <- retDT$row
  return(retMT)
}


# Load settings files ----------------------------------------------------------
settings.files <- list.files(input.folder, recursive=T, full.names=T, pattern="tf_settings.csv")
settings <- lapply(settings.files, function(x){
  res <- fread(x, fill=TRUE)
  data.table(res, file=dirname(paste0(x)))
})
settings <- do.call(rbind, settings)
settings <- dcast.data.table(settings, file ~ V1, value.var = "V2")
settings$replicate <- paste0("Replicate", 1:nrow(settings))
for(i in 1:nrow(settings)){
  te.file <- paste0(settings[i]$file, "/tf_NumGradTestError.txt")
  if(file.exists(te.file)){
    te <- fread(te.file)
    settings[i, TestError := te$NumGradTestError]
  } else {
    settings[i, TestError := NA]
  }
}
write.tsv(settings, file = paste0(output.folder, "/", "MetaData_byRun.tsv"))


# Get iterations and worst step from cost files (training curves) ----------------------------------------------------------
i <- 1
for(i in 1:nrow(settings)){
  cost.data <- tryCatch(fread(paste0(settings[i]$file, "/tf_cost.csv")), error = function(e){message("Failed to read training curve (tf_cost.csv): ", e)})
  settings[i,final.iterations := cost.data[nrow(cost.data)]$iteration]
  settings[i,worstStep := min(cost.data$error[1:(nrow(cost.data)-1)] - cost.data$error[2:nrow(cost.data)])]
}
write.tsv(settings, file = paste0(output.folder, "/", "MetaData_byRun.tsv"))

# Get time for training from Memory file ----------------------------------------------------------
i <- 728
for(i in 1:nrow(settings)){
  mem.file <- tryCatch(fread(paste0(settings[i]$file, "/Memory.tsv")), error = function(e){message("Failed to read memory file: ", e)})
  mem.file[, TimePOSX := as.POSIXct(Time)]
  time.diff.sec <- as.numeric(mem.file[Step == "Trained Tensorflow"]$TimePOS - mem.file[Step == "Prepared Training"]$TimePOSX, units = "secs")
  if(length(time.diff.sec) == 0)time.diff.sec <- NA 
  settings[i,time.diff.seconds := time.diff.sec]
}
write.tsv(settings, file = paste0(output.folder, "/", "MetaData_byRun.tsv"))


# Load Node weights -------------------------------------------------------
nodeWeights <- data.table()
outputs <- c()
i <- 1
for(i in 1:nrow(settings)){
  te <- tryCatch(fread(paste0(settings[i]$file, "/tf_NumGradMeans.csv")), error = function(e){message("Failed to read Node Weights: ", e)})
  te.outputs <- tryCatch(fread(paste0(settings[i]$file, "/outputs.txt"), header = F), error = function(e){message("Failed to read outputs: ", e)})
  outputs <- unique(c(outputs, te.outputs$V1))
  te$file <- settings[i]$file
  te$replicate <- settings[i]$replicate
  nodeWeights <- rbind(nodeWeights, te, fill=TRUE)
}
for(cc in c("Node", outputs)){
  if(is.null(nodeWeights[[cc]])) nodeWeights[[cc]] <- NA
}
nodeWeights <- melt(nodeWeights, measure.vars = outputs, id.vars = c("Node", "replicate", "file"))
#nodeWeights <- nodeWeights[!is.na(value)]


# Output a matrix of node weights and a table of meta information -------------------------------------------------------
if(length(outputs) == 1){
  write.tsv(settings, file = paste0(output.folder, "/", "NodeWeights_MetaData.tsv"))
  write.csv(toMT(nodeWeights, row = "Node", col="replicate", val = "value"), file = paste0(output.folder, "/", "NodeWeights.csv"), quote=F)
} else {
  nodeWeights[,id := paste0(replicate, "_", variable)]
  write.csv(toMT(nodeWeights, row = "Node", col="id", val = "value"), file = paste0(output.folder, "/", "NodeWeights.csv"), quote=F)
  
  settings2 <- data.table()
  for(i in 1:nrow(settings)){
    xx <- unique(nodeWeights[replicate == settings[i]$replicate][,c("replicate", "variable", "id")])
    settings2 <- rbind(settings2, data.table(settings[i][,-"replicate"], replicate = xx$id, outputNode=xx$variable))
  }
  write.tsv(settings2, file = paste0(output.folder, "/", "NodeWeights_MetaData.tsv"))
}

message("\nSuccessfully collected output from folder: ", input.folder, " and stored it in folder: ", output.folder)
