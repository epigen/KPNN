source(paste0(Sys.getenv("KPNN_CODEBASE"), "/init.R"))

# FILE SETTINGS -----------------------------------------------------------
out <- "Analysis/"
dir.create(dirout(out))


# Num Grad Error ----------------------------------------------------------
numGradTestError <- data.table()

(res.folders <- paste0(Sys.getenv("KPNN_OUTPUTS"), c("HCA", "TCR", "SIM1", "SIM2")))
fold <- res.folders[2]
for(fold in res.folders){
  res.files <- list.files(fold, recursive=T, full.names=T, pattern="tf_NumGradTestError.txt")
  if(length(res.files) == 0) next
  err <- lapply(res.files, function(x) fread(x)[[1]])
  names(err) <- gsub(".+(run_\\d+).+", "\\1", res.files)
  res.files.folders <- gsub("tf_NumGradTestError.txt", "", gsub(dirout(), "", res.files))
  numGradTestError <- rbind(numGradTestError,
        data.table(error=unlist(err), do.call(rbind, strsplit(res.files.folders, "\\/")), folder=res.files.folders)
  )
}
numGradTestError$V1 <- NULL
numGradTestError$V2 <- NULL
colnames(numGradTestError) <- c("Error", "Dataset", "Experiment", "Replicate", "Folder")
write.tsv(numGradTestError, dirout(out, "TestError.tsv"))

ggplot(numGradTestError, aes(x=Experiment, y=Error)) + geom_jitter(width=0.1, height=0)  + theme_bw(12) + xRot() + facet_wrap(~Dataset, scales="free")
ggsave(dirout(out, "TestError.pdf"), w=10, h=5)

# Accuracy for multiple classes ---------------------------------------------------------
acc.file <- dirout(out, "Unbalanced_Accuracy.RData")
if(file.exists(acc.file)){load(acc.file)} else {res.dt <- data.table()}
for(fold in res.folders){
  res.files <- list.files(fold, recursive=T, full.names=T, pattern="tf_yHat_val.csv")
  if(length(res.files) == 0) next
  ff <- res.files[1]
  for(ff in res.files){
    ff.clean <- gsub("tf_yHat_val.csv", "", gsub(dirout(), "", ff))
    if(ff.clean %in% res.dt$run) next
    message(ff)
    yHat <- read.csv(ff)
    head(yHat)
    yTrue <- read.csv(gsub("yHat", "yTrue", ff))
    head(yTrue)
    x <- t(apply(yHat, 1, function(row) (row == max(row)) + 0)) * yTrue
    col <- colnames(yTrue)[1]
    for(col in colnames(yTrue)){
      res.dt <- rbind(res.dt, data.table(
        run=ff.clean, ct=col, accuracy=mean(x[which(yTrue[,col] == 1),col])
        ))
    }
  }
}
colnames(res.dt) <- c("Folder", "OutputNode", "Accuracy")
write.tsv(res.dt, dirout(out, "MultiClass_Accuracy.tsv"))