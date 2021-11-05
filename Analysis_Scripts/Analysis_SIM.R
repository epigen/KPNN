#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------
inDir <- ""
out <- ""


# Check folders -----------------------------------------------------------
if(inDir == "") stop("Need to define an input directory, where the results of Analysis_CollectOutputs.R are stored")
if(out == "") stop("Need to define an output directory")
dir.create(out,recursive = TRUE)


# FUNCTIONS ---------------------------------------------------------------
xRot <- function(){theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))}
write.tsv <- function(...){ write.table(..., sep="\t", row.names=FALSE, quote=FALSE);}
minMax <- function(x, na.rm=T){ (x-min(x, na.rm=na.rm))/(max(x, na.rm=na.rm) - min(x, na.rm=na.rm))}
toMT <- function(dt, row, col, val){
  retDT <- dcast.data.table(dt, get(row) ~ get(col), value.var=val)
  retMT <- as.matrix(retDT[,-"row"])
  row.names(retMT) <- retDT$row
  return(retMT)
}



# LOAD DATA  ---------------------------------------------------------------
ff <- list.files(paste0(inDir, c("/SIM2/", "/SIM1/")), recursive=T, full.names=T, pattern="tf_NumGradMeans.csv")
fx <- ff[1]
numGradMeans.res <- lapply(ff, function(fx){
  fx <- gsub("\\/\\/", "/", fx)
  x <- fread(fx)
  x$name <- dirname(gsub(paste0(".+(",inDir, ")"), "\\1", gsub("\\/+", "/", fx)))
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
write.tsv(numGradMeans, file=paste0(out, "NodeWeights.tsv"))


# Filter By error ---------------------------------------------------------
err <- fread(paste(inDir, "TestError.tsv"))
numGradMeans <- numGradMeans[Folder %in% gsub("^.(.+).$", "\\1", err[Error < 0.2]$Folder)]
numGradMeans[, NormWeights := minMax(abs(Weight)), by=c("Experiment")]

# SIM1 --------------------------------------------------------------------
pDat <- numGradMeans[grepl("SIM1", Experiment)]
pDat$Node <- gsub("nodeX", "A", gsub("b", "", pDat$Node))
pDat$Node <- factor(pDat$Node, levels=rev(unique(pDat[,mean(NormWeights), by="Node"][order(V1)]$Node)))
ggplot(pDat, aes(x=Node, y=NormWeights)) + geom_jitter(height=0, shape=1) + theme_bw(12) + xRot()
ggsave(paste0(out, "SIM1.pdf"), width=15, h=6)


# SIM1 --------------------------------------------------------------------
err <- fread(paste0(inDir, "TestError.tsv"))
pDat <- numGradMeans[grepl("SIM2", Experiment)]
pDat[,Experiment := gsub("SIM2", "", Experiment)]
pDat$Node <- gsub("node", "", gsub("b", "", pDat$Node))
pDat$Node <- factor(pDat$Node, levels=rev(unique(pDat[,mean(NormWeights), by="Node"][order(V1)]$Node)))
ggplot(pDat, aes(x=Node, y=NormWeights)) + geom_jitter(height=0, shape=1) + xRot() + facet_wrap(~Experiment, ncol=1)
ggsave(paste0(out, "SIM2_AllNodes.pdf"), width=15, h=12)

pDat2 <- dcast.data.table(pDat[Node %in% c("A", "B")], Replicate + Experiment ~ Node, value.var="NormWeights")
ggplot(pDat2, aes(x=A, y=B)) + geom_point(shape=1) + facet_wrap(~Experiment) + theme_bw(12) +
  geom_text(data=pDat2[,round(cor(A,B),3), by="Experiment"], aes(label=paste("R =",  V1)), x=0.85, y=0.95)
ggsave(paste0(out, "SIM2_AvsB.pdf"), width=12, h=6)
