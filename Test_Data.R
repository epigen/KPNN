#!/usr/bin/env Rscript

require(data.table,quietly = TRUE)
require(optparse,quietly = TRUE)

message("------------")
usage=paste(getopt::get_Rscript_filename(), "[options] OutFile")
message(paste("Usage: ", usage))

# Parse options
option_list = list(
  make_option(c("--nInputs"), type="integer", default=20,help="Number of inputs (genes)", metavar="integer"),
  make_option(c("--nOutputs"), type="integer", default=3,help="Number of outputs (class labels)", metavar="integer"),
  make_option(c("--nNodes"), type="integer", default=5,help="Number of nodes", metavar="integer"),
  make_option(c("--nSamples"), type="integer", default=100,help="Number of samples (cells)", metavar="integer"),
  make_option(c("--Seed"), type="integer", default=NA,help="Seed for random sampling", metavar="integer")
)

opt_parser = OptionParser(usage = usage, option_list=option_list);
args <- parse_args(opt_parser, positional_arguments = 1)
opt  <- args$options    # Parameters
fx <- paste0(args$args[1], "/")      # output location

# opt <- list()
# opt$nOutputs <- 2
# opt$nNodes <- 2
# opt$nInputs <- 20
# opt$nSamples <- 100
# opt$Seed <- NA
# fx <- "~/KPNN/KPNN/Inputs/"

if(!is.na(opt$Seed)) set.seed(opt$Seed)

# Name elements
outputs = paste0("output", 1:opt$nOutputs)
nodes = paste0("node", 1:opt$nNodes)
inputs = paste0("input", 1:opt$nInputs)
samples= paste0("sample", 1:opt$nSamples)

# Edgelist
edgelist <- data.table()
edgelist <- rbind(edgelist, data.table(parent = outputs, child = rep(nodes, each=opt$nOutputs)))
edgelist <- rbind(edgelist, data.table(parent = nodes, child = rep(paste0(inputs, "_gene"), each=opt$nNodes)))
edges.rm <- sample(which(grepl("node", edgelist$parent)), ceiling(sum(grepl("node", edgelist$parent))*0.5))
edgelist <- edgelist[-edges.rm]
write.table(edgelist, paste0(fx, "TEST_Edgelist.csv"), sep=",", quote = F, row.names=F)

# Classlabels
classlabels <- data.table(barcode=samples)
for(ox in outputs){
  classlabels[[ox]] <- 0
}
if(length(outputs) == 1){
  classlabels[[outputs]][1:floor(opt$nSamples/2)] <- 1
} else {
  os <- split(1:opt$nSamples, outputs)
  for(ox in outputs){
    classlabels[[ox]][os[[ox]]] <- 1
  }
}
write.table(classlabels, paste0(fx, "TEST_ClassLabels.csv"), sep=",", quote = F, row.names=F)

# Data
data <- matrix(abs(rnorm(opt$nSamples * opt$nInputs)), ncol = opt$nSamples, nrow = opt$nInputs)
row.names(data) <- inputs
colnames(data) <- samples
# Include a difference
i2o <- split(inputs, f = sample(outputs,size = opt$nInputs, replace = T))
ox <- outputs[1]
for(ox in outputs){
  data[i2o[[ox]],classlabels[get(ox)==1]$barcode] <- data[i2o[[ox]],classlabels[get(ox)==1]$barcode]+2
}
write.table(data, paste0(fx, "TEST_Data.csv"), sep=",", quote = F)

message("Successfully generated test data\n")