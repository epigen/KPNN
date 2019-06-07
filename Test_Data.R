source(paste0(Sys.getenv("KPNN_CODEBASE"), "/init.R"))

# list.files(Sys.getenv("KPNN_INPUTS"))

# Define size of network
nInputs=20
nOutputs=2
nNodes=2
nSamples=100

# Name elements
outputs = paste0("output", 1:nOutputs)
nodes = paste0("node", 1:nNodes)
inputs = paste0("input", 1:nInputs)
samples= paste0("sample", 1:nSamples)

# Data
data <- matrix(abs(rnorm(nSamples * nInputs)), ncol = nSamples, nrow = nInputs)
row.names(data) <- inputs
colnames(data) <- samples
write.table(data, paste0(Sys.getenv("KPNN_INPUTS"), "TEST_Data.csv"), sep=",", quote = F)

# Edgelist
edgelist <- data.table()
edgelist <- rbind(edgelist, data.table(parent = outputs, child = rep(nodes, each=nOutputs)))
edgelist <- rbind(edgelist, data.table(parent = nodes, child = rep(paste0(inputs, "_gene"), each=nNodes)))
edges.rm <- sample(which(grepl("node", edgelist$parent)), ceiling(sum(grepl("node", edgelist$parent))*0.5))
edgelist <- edgelist[-edges.rm]
write.table(edgelist, paste0(Sys.getenv("KPNN_INPUTS"), "TEST_Edgelist.csv"), sep=",", quote = F, row.names=F)

# Classlabels
classlabels <- data.table(barcode=samples)
for(ox in outputs){
  classlabels[[ox]] <- 0
}
if(length(outputs) == 1){
  classlabels[[outputs]][1:floor(nSamples/2)] <- 1
} else {
  os <- split(1:nSamples, outputs)
  for(ox in outputs){
    classlabels[[ox]][os[[ox]]] <- 1
  }
}
write.table(classlabels, paste0(Sys.getenv("KPNN_INPUTS"), "TEST_ClassLabels.csv"), sep=",", quote = F, row.names=F)