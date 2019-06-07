source(paste0(Sys.getenv("KPNN_CODEBASE"), "/init.R"))


# SETUP -----------------------------------------------------------
inDir <- paste0("Analysis/")
out <- paste0(inDir, "HCA/")
dir.create(dirout(out))



# LOAD DATA  ---------------------------------------------------------------
ff <- list.files(paste0(Sys.getenv("KPNN_OUTPUTS"), c("/HCA/")), recursive=T, full.names=T, pattern="tf_NumGradMeans.csv")
fx <- ff[1]
for(type.x in c("Controls", "Real")){
  if(type.x == "Controls") ff.use <- ff[grepl("\\/Control\\/", ff)]
  numGradMeans.res <- lapply(ff.use, function(fx){
    fx <- gsub("\\/\\/", "/", fx)
    #if(!gsub(dirout(), "", dirname(fx)) %in% testErr$folder) return(NULL)
    x <- fread(fx)
    x$name <- dirname(gsub(dirout(), "", gsub("\\/+", "/", fx)))
    x
  })
  numGradMeans <- do.call(rbind, numGradMeans.res)
  numGradMeans <- cbind(numGradMeans, do.call(rbind, strsplit(numGradMeans$name, "/")))
  numGradMeans[, run := V4]
  numGradMeans[, metarun := V3]
  numGradMeans$V1 <- NULL
  numGradMeans$V2 <- NULL
  numGradMeans$V3 <- NULL
  numGradMeans$V4 <- NULL
  colnames(numGradMeans) <- c("Node", "Weight", "Folder", "Replicate", "Experiment")
  write.tsv(numGradMeans, file=dirout(out, "NodeWeights_",type.x,".tsv"))
}
