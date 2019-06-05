require(project.init)
require(Matrix)


# FILE SETTINGS -----------------------------------------------------------

project.init2("pathway_learning") 
out <- paste0("46_Analysis2/", "HCA2/")
dir.create(dirout(out))

datasets <- c("46_05_HCA2")
ds <- datasets[1]

outDS <- paste0(out, ds, "/")
dir.create(dirout(outDS))


# PROCESS ONE DATASET -----------------------------------------------------

(load(dirout("46_Analysis2/Unbalanced_Accuracy.RData")))
acc <- res.dt
acc[,cnt := .N, by="run"]
acc[,cntSuccess := sum(accuracy > 0.90), by="run"]
acc <- acc[cnt == cntSuccess]
acc[, path := gsub("\\/\\/", "/", dirout(run))]
acc[, path := gsub("\\/$", "", path)]
acc <- cbind(acc, do.call(rbind, strsplit(acc$run, "/")))
acc[, length(unique(V5)), by="V4"]
lapply(with(acc, split(cntSuccess, factor(gsub("run_\\d+", "", run)))), quantile)



# READ NUMGRAD DATA -------------------------------------------------------

ff <- list.files(dirout(ds, "/"), recursive=T, full.names=T, pattern="tf_NumGradMeans.csv")
fx <- ff[1]
numGradMeans.res <- lapply(ff, function(fx){
  fx <- gsub("\\/\\/", "/", fx)
  if(!dirname(fx) %in% acc$path) return(NULL)
  x <- fread(fx)
  x$name <- dirname(gsub(dirout(), "", gsub("\\/+", "/", fx)))
  x
})


# Format data -------------------------------------------------------------
outputs <- fread(paste(dirname(fx), "outputs.txt", sep="/"), header=F)[[1]]
numGradMeans.res <- numGradMeans.res[!sapply(numGradMeans.res, is.null)]
numGradMeans.ctrl <- do.call(rbind, numGradMeans.res[sapply(numGradMeans.res, function(x) "output" %in% colnames(x))])
numGradMeans.ctrl <- melt(numGradMeans.ctrl, measure.vars="output")
numGradMeans.real <- do.call(rbind, numGradMeans.res[sapply(numGradMeans.res, function(x) all(outputs %in% colnames(x)))])
numGradMeans.real <- melt(numGradMeans.real, measure.vars=outputs)
numGradMeans <- rbind(numGradMeans.real, numGradMeans.ctrl)
numGradMeans <- cbind(numGradMeans, do.call(rbind, strsplit(numGradMeans$name, "/")))


# NORMALIZE ---------------------------------------------------------------
numGradMeans[,valueAbs := abs(value)]
numGradMeans[,maxAbs := max(valueAbs), by=c("variable", "V1", "V3", "V4")]
numGradMeans[,valueNorm01 := valueAbs / maxAbs]
stopifnot(max(numGradMeans$valueNorm01) == 1)
numGradMeans[,id := paste(V1, V3, V4, variable)]
qnMT <- toMT(dt=numGradMeans, row="Node",col="id",val="valueAbs")
qnMT <- limma::normalizeQuantiles(qnMT)
qnDT <- data.table(melt(qnMT))
colnames(qnDT) <- c("Node", "id", "valueQN")
numGradMeans$valueQN <- NULL
numGradMeans <- merge(numGradMeans, qnDT, by=c("Node", "id"), all.x=T)
ggplot(numGradMeans, aes(x=id, y=valueQN)) + geom_boxplot() + scale_y_log10() + theme_bw(12) + xRot()
ggsave(dirout(outDS, "DistributionQN.jpg"), width=20, height=6)


# SPLIT INTO REAL AND CONTROL ---------------------------------------------
numGradMeansControl <- numGradMeans[variable == "output"]
numGradMeansControl[,output := value]
numGradMeans <- numGradMeans[variable != "output"]
numGradMeans[,length(unique(V4)), by="V3"]
numGradMeansControl[,length(unique(V4)), by="V3"]


(experiments <- unique(numGradMeans$V3))


# LIMIT TO EQUAL NUMBERS --------------------------------------------------
numGradMeans <- numGradMeans[!grepl("4G4", V3)]
numGradMeans[,length(unique(V4)), by="V3"]
(minCount <- min(numGradMeans[][,length(unique(V4)), by="V3"]$V1))
for(exp in unique(numGradMeans$V3)){
  numGradMeans <- numGradMeans[V3 != exp | V4 %in% sort(unique(numGradMeans[V3 == exp]$V4))[1:minCount]]
}
numGradMeans[,length(unique(V4)), by="V3"]


# SAVE
save(numGradMeans, numGradMeansControl, file=dirout(outDS, "NumGradMeans.RData"))

ggplot(numGradMeans[, length(unique(V4)), by="V3"], aes(x=V3, y=V1)) + geom_bar(stat='identity') + theme_bw(16) + xRot()
ggsave(dirout(out, ds, "_Numbers.pdf"), height=5, width=4)


# Plots for each experiments --------------------------------------------------
expt <- experiments[3]
for(expt in experiments){
  load(dirout(outDS, "NumGradMeans.RData"))
  numGradMeans <- numGradMeans[V3 == expt]
  if(nrow(numGradMeans) == 0) next
  
  edges <- fread(dirout(numGradMeans$name[1], "/edges.tsv"))
  
  #Distributions
  ggplot(numGradMeans, aes(x=value, color=variable)) + stat_ecdf() + facet_wrap(~V4, ncol=5) + theme_bw(12)
  ggsave(dirout(outDS,  "Distribution_",expt,".pdf"), 
         height=min(ceiling(length(unique(numGradMeans$V4))/5) * 2 + 1, 29), 
         width=10 + 2)
  # range
  ggplot(numGradMeans[, .(range = max(value) - min(value), median=median(value)), by=c("V4", "variable")], 
         aes(x=variable, y=range, fill=variable)) + geom_violin() + theme_bw(16) + xRot() + geom_jitter(height=0)
  ggsave(dirout(outDS,  "DistributionRange_",expt,".pdf"),height=4, width=3)
  
  # CORRELATIONS
  numGradMeans[,id := paste(variable, V4)]
  numGradMT <- toMT(dt=numGradMeans, "Node", "id", val="valueAbs")
  str(numGradMT)
  numGradCor <- cor(numGradMT, method="pearson")
  diag(numGradCor) <- NA
  dimx <- ncol(numGradMT) * 0.2 + 3
  cleanDev()
  pdf(dirout(outDS,  "CorHM_",expt,".pdf"), width=dimx, height=dimx,onefile=F)
  pheatmap(numGradCor, fontsize_row=6, fontsize_col=6)
  dev.off()
  
  # CORRELATION WITHIN / BETWEEN
  x <- numGradCor
  x[upper.tri(x, diag=TRUE)] <- NA
  x <- melt(x,value.factor=F)
  x$Var1 <- as.character(x$Var1)
  x$Var2 <- as.character(x$Var2)
  x <- data.table(x)
  x <- x[!is.na(value)]
  x <- cbind(x, do.call(rbind, strsplit(x$Var1, " ")))
  x <- cbind(x, do.call(rbind, strsplit(x$Var2, " ")))
  colnames(x) <- make.unique(colnames(x))
  x[V1 != V1.1 & V2 == V2.1, type := "sameRun_diffCell"]
  x[V1 == V1.1 & V2 != V2.1, type := "sameCell_diffRun"]
  x[V1 != V1.1 & V2 != V2.1, type := "diff"]
  save(x, file=dirout(outDS, "Correlation_WithinBetween_", expt, ".RData"))
  ggplot(x[!is.na(type)], aes(x=value, color=type)) + geom_density() + theme_bw(12)
  ggsave(dirout(outDS, "Correlation_WithinBetween_", expt, ".pdf"), width=6, height=5)
  
  # VALUE HEATMAPS
  hm.col <- colorRampPalette(c("lightgrey",  "blue"))(50)
  
  # RECEPTORS
  cleanDev()
  receptors <- edges[parent %in% outputs]$child
  numGradMeans[,id2 := paste(V4, variable)]
  recDat <- toMT(numGradMeans[Node %in% receptors], row="Node",col="id2", val="valueAbs")
  gapx <- cumsum(table(gsub("(run_\\d+).+", "\\1", colnames(recDat))))
  pdf(dirout(outDS, "RECHM_P_", expt, ".pdf"), width=dimx, height=12,onefile=F)
  pheatmap(recDat[,sort(colnames(recDat))], cluster_cols=F, fontsize_row=6, fontsize_col=6, gaps_col=gapx, color=hm.col)
  dev.off()
  
  # TFs
  edges.proteins <- edges[!grepl("_gene", child)]
  tfs <- edges.proteins[!child %in% edges.proteins$parent]$child
  tfDat <- toMT(numGradMeans[Node %in% tfs], row="Node",col="id2", val="valueAbs")
  quantile(apply(tfDat, 1, max))
  gapx <- cumsum(table(gsub("(run_\\d+).+", "\\1", colnames(recDat))))
  pdf(dirout(outDS, "TFHM_P_", expt, ".pdf"), width=25, height=10,onefile=F)
  pheatmap(tfDat[apply(tfDat, 1, max) > 0.001,sort(colnames(recDat))], cluster_cols=F, fontsize_row=6, fontsize_col=6, gaps_col=gapx, color=hm.col)
  dev.off()
  
  # REACHABILITY
  g <- graph.edgelist(as.matrix(edges[,c("parent", "child")]))
  reach <- apply(distances(g, v=which(!grepl("_gene", V(g)$name)), to=which(grepl("_gene", V(g)$name)), mode="out") != Inf, 1, sum)
  reach <- data.table(reach=reach, Node=names(reach))
  reach$depth <- distances(g, v=which(V(g)$name %in% outputs[1]), to=which(V(g)$name %in% edges[!parent %in% outputs]$parent), mode="out")[1,][reach$Node]
  ggplot(merge(numGradMeans[,mean(valueAbs), by=c("Node", "variable")], reach, by="Node"), aes(x=reach, y=V1)) + geom_point(color="darkturquoise", alpha=0.7) +
    geom_text(aes(label=Node), alpha=0.8) + facet_grid(variable ~ .) + theme_bw(16) +
    ylab("Mean Node Importance") + xlab("Number of genes connected to")
  ggsave(dirout(outDS, "NodeVsReach_", expt, ".pdf"), width=12, height=min(25, length(unique(numGradMeans$variable)) * 3+1))

  # celltype vs celltype
  outputs <- unique(as.character(numGradMeans$variable))
  for(oi1 in 1:(length(outputs) -1)){
    for(oi2 in (oi1 + 1):length(outputs)){
      x1MT <- toMT(dt=numGradMeans[variable == outputs[oi1]][order(V4)], row="Node", col="V4", val="valueQN")
      x2MT <- toMT(dt=numGradMeans[variable == outputs[oi2]][order(V4)], row="Node", col="V4", val="valueQN")
      stopifnot(all(colnames(x1MT) == colnames(x2MT)))
      
      pDat <- data.table(
        Node=row.names(x1MT),
        p=sapply(row.names(x1MT), function(node) t.test(x1MT[node,], x2MT[node,], paired=T)$p.value),
        mean1=rowMeans(x1MT[row.names(x1MT),]),
        mean2=rowMeans(x2MT[row.names(x1MT),]))
      pDat[,diff := mean1 - mean2]
      write.tsv(pDat, dirout(outDS, paste("NodeVsNodes", expt, outputs[oi1], outputs[oi2], sep="_"), ".tsv"))
      
      ggplot(pDat, aes(x=mean1, y=mean2, color=-log10(p), size=-log10(p))) + geom_abline() + geom_point(alpha=0.7) +
        scale_color_gradient(high="darkturquoise", low="grey") +
        geom_text(data=pDat[p < 0.05], aes(label=Node), alpha=0.8, color="black") + theme_bw(16) + 
        xlab(outputs[oi1]) + ylab(outputs[oi2])
      ggsave(dirout(outDS, paste("NodeVsNodes", expt, outputs[oi1], outputs[oi2], sep="_"), ".pdf"), width=8, height=8)
      
      ggplot(pDat, aes(x=mean1, y=mean2, color=-log10(p), size=-log10(p))) + geom_abline() + geom_point(alpha=0.7) +
        scale_color_gradient(high="darkturquoise", low="grey") +
        geom_text(data=pDat[p < 0.05], aes(label=Node), alpha=0.8, color="black") + theme_bw(16) + 
        xlab(outputs[oi1]) + ylab(outputs[oi2]) + xlim(0, 0.1) + ylim(0, 0.1)
      ggsave(dirout(outDS, paste("NodeVsNodes", expt, outputs[oi1], outputs[oi2], "zoom", sep="_"), ".pdf"), width=8, height=8)
    }
  }
}


load(dirout(outDS, "NumGradMeans.RData"))


# COLLECT NODE STATS ------------------------------------------------------

node.stats <- data.table()
for(ds in datasets){
  ff <- list.files(dirout(out, ds, "/"), pattern="NodeVsNodes_.+\\.tsv", full.names=T)
  for(fx in ff){
    node.stats <- rbind(node.stats, data.table(fread(fx), file=gsub("NodeVsNodes_(.+)\\.tsv", "\\1", basename(fx)), dataset=ds))
  }
}
node.stats <- cbind(node.stats, do.call(rbind, strsplit(node.stats$file, "_")))
node.stats[,qval := p.adjust(p, method="BY")]
save(node.stats, file=dirout(outDS, "NodeStats.RData"))


stopifnot(length(table(node.stats[,.N, by="V1"]$N)) == 1)
ggplot(node.stats[dataset==ds][qval < 0.05], aes(x=V1)) + geom_bar() + facet_grid(paste(V2, "vs", V3) ~ .) + theme_bw(16) + xRot()
ggsave(dirout(outDS, ds, "_vsCell_NumberOfHits.pdf"), width=6, height=10)



# SUMMARY PLOTS ------------------------------------------------------------

node.stats[,log2FC := log2(mean1/mean2)]
node.stats2 <- node.stats
node.stats2$diff <- -node.stats$diff
node.stats2$log2FC <- -node.stats$log2FC
node.stats2$V2 <- node.stats$V3
node.stats2$V3 <- node.stats$V2
node.stats2$file <- paste0(node.stats$file, 2)
ns <- rbind(node.stats, node.stats2)
save(ns, file=dirout(outDS, "NodeStats_symmetric.RData"))

nsx <- ns[dataset == ds]
unique(nsx[qval < 0.01]$Node)

edges <- fread(list.files(dirout(nsx$dataset[1]), pattern="edges.tsv", recursive=T, full.names=T)[1])
nodes.outputs <- unique(edges[!parent %in% edges$child]$parent)
nodes.rec <- unique(edges[parent %in% nodes.outputs]$child)
nodes.non.rec <- unique(edges[!parent %in% nodes.rec & !parent %in% nodes.outputs]$parent)

ggplot(nsx[grepl("STAT", Node)], aes(x=V1, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=nsx[grepl("STAT", Node) & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(V2 ~ V3, scales="free", space="free") + theme_bw(16) + xRot() +
  #scale_alpha_manual(values=c(0.2, 1)) +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsCell_STATs.pdf"), width=10, height=10)

ggplot(nsx[Node %in% nodes.rec], aes(x=V1, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=nsx[Node %in% nodes.rec & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(V2 ~ V3, scales="free", space="free") + theme_bw(16) + xRot() +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsCell_RECs.pdf"), width=10, height=25)

nodes.non.rec.sig <- intersect(nodes.non.rec, unique(nsx[qval < 0.01 & abs(diff) > 0.001]$Node))
ggplot(nsx[Node %in% nodes.non.rec.sig], aes(x=V1, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=nsx[Node %in% nodes.non.rec.sig & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(V2 ~ V3, scales="free", space="free") + theme_bw(16) + xRot() +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsCell_SigNOdes.pdf"), width=10, height=length(nodes.non.rec.sig) * 0.6 + 2)





# COLLECT CORRELATIONS ------------------------------------------------------

corr <- data.table()
for(ds in datasets){
  ff <- list.files(dirout(out, ds, "/"), pattern="Correlation_WithinBetween_.+.RData", full.names=T)
  fx <- ff[1]
  for(fx in ff){
    load(fx)
    corr <- rbind(corr, data.table(x, file=gsub("Correlation_WithinBetween_(.+)\\.RData", "\\1", basename(fx)), dataset=ds))
  }
}
corr <- corr[dataset == ds]
corr[,dropout := gsub("(.G.).*", "\\1", file)]
corr[,subset := gsub("(.G.)(..).*", "\\2", file)]
corr[subset == "no", subset := "CB"]
save(corr, file=dirout(outDS, "Correlations.RData"))
ggplot(corr, aes(y=value, x=dropout)) + 
  geom_boxplot(aes(color=type), coef=1e6, position=position_dodge(1)) + 
  facet_grid(. ~ subset) + theme_bw(16) + xRot()
ggsave(dirout(out, ds, "_Correlations.pdf"), height=4, width=8)






# VS CONTROL --------------------------------------------------------------
(load(dirout(outDS, "NumGradMeans.RData")))
statVsControl <- data.table()
for(exp in unique(numGradMeans$V3)){
  message(exp)
  for(ct in unique(numGradMeans[V3 == exp]$variable)){
    for(nodex in unique(numGradMeans[V3 == exp & variable == ct]$Node)){
      vR <- numGradMeans[V3 == exp & variable == ct & Node == nodex]$valueQN
      vC <- numGradMeansControl[Node == nodex]$valueQN
      pDat <- data.table(
        Node=nodex,
        p=t.test(vR, vC)$p.value,
        mean1=mean(vR),
        mean2=mean(vC),
        experiment=exp,
        output=ct)
      pDat[,diff := mean1 - mean2]
      statVsControl <- rbind(statVsControl, pDat)
    }
  }
}
statVsControl[,qval := p.adjust(p, method="BY")]
save(statVsControl, file=dirout(outDS, "NodeStats_vsControl.RData"))

# Numbers
ggplot(statVsControl[qval < 0.05], aes(x=experiment)) + geom_bar() + theme_bw(16) + xRot()
ggsave(dirout(out, ds, "_vsControl_NumberOfHits.pdf"), width=6, height=4)

ggplot(statVsControl[grepl("STAT", Node)], aes(x=experiment, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=statVsControl[grepl("STAT", Node) & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(. ~ output, scales="free", space="free") + theme_bw(16) + xRot() +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsControl_STATs",".pdf"), width=10, height=5)

ggplot(statVsControl[Node %in% nodes.rec], aes(x=experiment, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=statVsControl[Node %in% nodes.rec & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(. ~ output, scales="free", space="free") + theme_bw(16) + xRot() +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsControl_RECs.pdf"), width=10, height=15)

nodes.non.rec.sig <- intersect(nodes.non.rec, unique(statVsControl[qval < 0.01 & abs(diff) > 0.01]$Node))
ggplot(statVsControl[Node %in% nodes.non.rec.sig], aes(x=experiment, y=Node, color=diff, size=-log10(p))) + 
  geom_point(data=statVsControl[Node %in% nodes.non.rec.sig & qval < 0.05], aes(size=-log10(p) + 2), color="black") + 
  geom_point() + facet_grid(. ~ output, scales="free", space="free") + theme_bw(16) + xRot() +
  scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, ds, "_vsControl_SigNOdes",".pdf"), width=10, height=length(nodes.non.rec.sig) * 0.3 + 2)

ggplot(statVsControl, aes(y=diff, x=paste(experiment, output))) + geom_boxplot()



# BM vs CB ----------------------------------------------------------------
(load(dirout(outDS, "NumGradMeans.RData")))
numGradMeans[V3 == "XGXnoDAdj", V3 := "XGXCBnoDAdj"]
(x <- gsub("BM", "..", gsub("CB", "..", unique(numGradMeans$V3))))
expx <- x[1]
statVsOrgans <- data.table()
for(expx in unique(x)){
  expx.x <- unique(sort(numGradMeans[grepl(expx, V3)]$V3))
  gprs <- gsub(".*(CB).*", "\\1", gsub(".*(BM).*", "\\1", expx.x))
  if(length(expx.x) != 2) next
  ct <- as.character(unique(numGradMeans[grepl(expx, V3)]$variable))[1]
  nodex <- numGradMeans$Node[1]
  for(ct in as.character(unique(numGradMeans[grepl(expx, V3)]$variable))){
    for(nodex in unique(numGradMeans[grepl(expx, V3) & variable == ct]$Node)){
      v1 <- numGradMeans[V3 == expx.x[1] & variable == ct & Node == nodex]$valueQN
      v2 <- numGradMeans[V3 == expx.x[2] & variable == ct & Node == nodex]$valueQN
      pDat <- data.table(
        Node=nodex,
        p=t.test(v1, v2)$p.value,
        mean1=mean(v1),
        mean2=mean(v2),
        experiment=expx,
        output=ct,
        comparison = paste(gprs,collapse=" ")
        )
      pDat[,diff := mean1 - mean2]
      statVsOrgans <- rbind(statVsOrgans, pDat)
    }
  }
}

statVsOrgans[,qval := p.adjust(p, method="BH")]
ggplot(statVsOrgans, aes(x=p)) + facet_grid(experiment ~ output) + geom_histogram()
statVsOrgans[order(qval)][1:20]
save(statVsOrgans, file=dirout(outDS, "NodeStats_vsOrgan.RData"))

ggplot(statVsOrgans[qval < 0.05], aes(x=experiment, color=diff > 0, fill=output)) + geom_bar(position="dodge")
ggsave(dirout(out, ds, "_vsOrgan_Numbers.pdf"), height=6, width=12)

table(statVsOrgans$comparison)
ggplot(statVsOrgans[Node %in% unique(statVsOrgans[order(qval)]$Node)[1:10]], aes(x=experiment, y=Node, color=diff, size=pmin(5, -log10(qval)))) + geom_point() +
  facet_grid(. ~ output) + theme_bw(16) + xRot() + scale_color_gradient2(high="indianred", low="navyblue") + ggtitle(statVsOrgans$comparison[1])
ggsave(dirout(out, ds, "_vsOrgan.pdf"), height=6, width=12)


# BM vs CB using LM ----------------------------------------------------------------
(load(dirout(outDS, "NumGradMeans.RData")))
numGradMeans[V3 == "XGXnoDAdj", V3 := "XGXCBnoDAdj"]
(x <- gsub("BM", "..", gsub("CB", "..", unique(numGradMeans$V3))))
statVsOrgans.LM <- data.table()
expx <- x[1]
for(expx in unique(x)){
  expx.x <- sort(unique(sort(numGradMeans[grepl(expx, V3)]$V3)))
  if(length(expx.x) != 2) next
  nodex <- numGradMeans$Node[1]
  for(nodex in unique(numGradMeans[grepl(expx, V3)]$Node)){
    datx <- numGradMeans[Node == nodex & V3 %in% expx.x]
    datx$V3 <- factor(datx$V3, levels=expx.x)
    fit <- lm(valueQN ~ variable + V3, data=datx)
    coef.V3 <- grep("CB", names(coef(fit)), value=T)
    pDat <- data.table(
      Node=nodex,
      p=coef(summary(fit))[coef.V3, 4],
      estimate=coef(summary(fit))[coef.V3, 1],
      mean1=mean(datx[grepl("BM", V3)]$valueQN),
      mean2=mean(datx[grepl("CB", V3)]$valueQN),
      experiment=expx
    )
    pDat[,BM_vs_CB := mean1 - mean2]
    statVsOrgans.LM <- rbind(statVsOrgans.LM, pDat)
  }
}

statVsOrgans.LM[,qval := p.adjust(p, method="BH")]
ggplot(statVsOrgans.LM, aes(x=p)) + facet_grid(experiment ~ .) + geom_histogram()
statVsOrgans.LM[order(qval)][1:20]
save(statVsOrgans.LM, file=dirout(outDS, "NodeStats_vsOrgan_LM.RData"))

ggplot(statVsOrgans.LM[qval < 0.05], aes(x=experiment, fill=BM_vs_CB > 0)) + geom_bar()
ggsave(dirout(out, ds, "_vsOrgan_LM_Numbers.pdf"), height=6, width=12)

statVsOrgans.LM[Node == "KLF4"]

ggplot(statVsOrgans.LM[Node %in% unique(statVsOrgans.LM[qval < 0.05]$Node)], 
       aes(x=experiment, y=Node, color= BM_vs_CB, size=pmin(5, -log10(qval)))) + geom_point() +
  theme_bw(12) + xRot() + scale_color_gradient2(high="indianred", low="navyblue")
ggsave(dirout(out, ds, "_vsOrgan_LM.pdf"), height=29, width=12)