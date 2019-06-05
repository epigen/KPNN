require(data.table)
require(ggplot2)


dirout <- function (...){
  paste0(Sys.getenv("KPNN_ROOT"), ...)
}


xRot <- function(){theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))}

write.tsv <- function(...){
  write.table(..., sep="\t", row.names=FALSE, quote=FALSE);
}

toMT <- function(dt, row, col, val){
  retDT <- dcast.data.table(dt, get(row) ~ get(col), value.var=val)
  retMT <- as.matrix(retDT[,-"row"])
  row.names(retMT) <- retDT$row
  return(retMT)
}

minMax <- function(x){ (x-min(x))/(max(x) - min(x))}
