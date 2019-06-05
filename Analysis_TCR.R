source(paste0(Sys.getenv("KPNN_CODEBASE"), "/init.R"))


# SETUP -----------------------------------------------------------
inDir <- paste0("Analysis/")
out <- paste0(inDir, "TCR/")
dir.create(dirout(out))



# LOAD DATA  ---------------------------------------------------------------
ff <- list.files(paste0(Sys.getenv("KPNN_OUTPUTS"), c("/TCR/")), recursive=T, full.names=T, pattern="tf_NumGradMeans.csv")
fx <- ff[1]
numGradMeans.res <- lapply(ff, function(fx){
  fx <- gsub("\\/\\/", "/", fx)
  #if(!gsub(dirout(), "", dirname(fx)) %in% testErr$folder) return(NULL)
  x <- fread(fx)
  x$name <- dirname(gsub(dirout(), "", gsub("\\/+", "/", fx)))
  x
})
numGradMeans <- lapply(numGradMeans.res, function(x) {colnames(x) <- c("variable", "numGradMean", "name"); return(x)})
numGradMeans <- do.call(rbind, numGradMeans)
numGradMeans <- cbind(numGradMeans, do.call(rbind, strsplit(numGradMeans$name, "/")))
numGradMeans[, run := V4]
numGradMeans[, metarun := V3]
numGradMeans$V1 <- NULL
numGradMeans$V2 <- NULL
numGradMeans$V3 <- NULL
numGradMeans$V4 <- NULL
colnames(numGradMeans) <- c("Node", "Weight", "Folder", "Replicate", "Experiment")
write.tsv(numGradMeans, file=dirout(out, "NodeWeights.tsv"))



# SIM1 --------------------------------------------------------------------


