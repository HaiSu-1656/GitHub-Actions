#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(DESeq2))
outdir <- "results/GSE89408"
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
input <- "GSE89408_GEO_count_matrix_rename.txt.gz"
stopifnot(file.exists(input))
dat <- read.delim(gzfile(input),row.names=1,check.names=FALSE)
stopifnot(nrow(dat)>20000L,ncol(dat)==218L,!anyDuplicated(rownames(dat)))
samples <- c(grep("^normal_tissue_",colnames(dat),value=TRUE),
             grep("^OA_tissue_",colnames(dat),value=TRUE))
stopifnot(length(samples)==50L,
          sum(startsWith(samples,"normal_tissue_"))==28L,
          sum(startsWith(samples,"OA_tissue_"))==22L)
mat <- as.matrix(dat[,samples,drop=FALSE])
storage.mode(mat) <- "numeric"
stopifnot(!anyNA(mat),all(is.finite(mat)),all(mat>=0))
fractional <- sum(abs(mat-round(mat))>1e-8)
fractional_genes <- sum(rowSums(abs(mat-round(mat))>1e-8)>0)
message("input genes=",nrow(mat)," selected samples=",ncol(mat),
        " fractional values=",fractional," fractional genes=",fractional_genes)
stopifnot("IRF7" %in% rownames(mat),sum(rownames(mat)=="IRF7")==1L)
cts <- round(mat)
storage.mode(cts) <- "integer"
group <- factor(ifelse(startsWith(samples,"OA_tissue_"),"OA","Healthy"),levels=c("Healthy","OA"))
meta <- data.frame(group=group,row.names=samples)
dds <- DESeqDataSetFromMatrix(countData=cts,colData=meta,design=~group)
# Keep at least 5 samples >= 10 counts; filtering is applied uniformly to all genes.
keep <- rowSums(counts(dds)>=10)>=5
message("Passing count filter: ",sum(keep))
stopifnot(keep["IRF7"])
dds <- dds[keep,]
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("group","OA","Healthy"),alpha=0.05,independentFiltering=TRUE)
all <- data.frame(gene=rownames(res),as.data.frame(res),check.names=FALSE)
write.csv(all,file.path(outdir,"GSE89408_OA_vs_Healthy_DESeq2_all_genes.csv"),row.names=FALSE,na="NA")
irf <- all[all$gene=="IRF7",,drop=FALSE]
stopifnot(nrow(irf)==1L)
irf$n_OA <- 22L
irf$n_Healthy <- 28L
irf$comparison <- "OA_vs_Healthy"
write.csv(irf,file.path(outdir,"GSE89408_IRF7_DESeq2_statistics.csv"),row.names=FALSE,na="NA")
norm <- counts(dds,normalized=TRUE)["IRF7",]
per <- data.frame(sample=samples,group=as.character(group),
                  raw_count=as.numeric(mat["IRF7",]),
                  rounded_count=as.integer(cts["IRF7",]),
                  size_factor=as.numeric(sizeFactors(dds)),
                  normalized_count=as.numeric(norm))
write.csv(per,file.path(outdir,"GSE89408_IRF7_50_patients.csv"),row.names=FALSE)
fmt <- function(v) if(length(v)==0L||is.na(v)) "NA" else formatC(v,format="g",digits=10)
summ <- c("GSE89408 human synovium | native R DESeq2 (rounded estimated counts)",
          paste("DESeq2 version",packageVersion("DESeq2")),
          paste("Genes in input",nrow(mat)),
          paste("Genes analyzed",nrow(dds)),
          paste("Non-integer values",fractional),
          paste("Genes with non-integer values",fractional_genes),
          "Note: RSEM estimated counts rounded before DESeq2. This is an approximation.",
          "Comparison: OA (n=22) vs Healthy (n=28)",
          paste("IRF7 baseMean",fmt(irf$baseMean)),
          paste("IRF7 log2FoldChange",fmt(irf$log2FoldChange)),
          paste("IRF7 lfcSE",fmt(irf$lfcSE)),
          paste("IRF7 Wald_stat",fmt(irf$stat)),
          paste("IRF7 pvalue",fmt(irf$pvalue)),
          paste("IRF7 padj",fmt(irf$padj)),
          paste("IRF7 fold change",fmt(2^irf$log2FoldChange)),
          paste("IRF7 Healthy normalized mean",fmt(mean(norm[group=="Healthy"]))),
          paste("IRF7 OA normalized mean",fmt(mean(norm[group=="OA"]))),
          paste("IRF7 Healthy normalized median",fmt(median(norm[group=="Healthy"]))),
          paste("IRF7 OA normalized median",fmt(median(norm[group=="OA"]))),
          paste("Genes with adjusted P",sum(!is.na(all$padj))))
writeLines(summ,file.path(outdir,"GSE89408_IRF7_summary.txt"))
cat(paste(summ,collapse="\n"),"\n")
png(file.path(outdir,"GSE89408_IRF7_DESeq2_expression.png"),width=2000,height=1500,res=220)
y <- log2(norm+1)
set.seed(89408)
x <- ifelse(group=="Healthy",1,2)
boxplot(y~factor(x,levels=1:2),outline=FALSE,
  at=c(1,2),names=c("Healthy (n=28)","OA (n=22)"),
  col=c("#d9e5ec","#edd9d3"),border="#6d6d6d",
  ylab="log2(DESeq2 normalized counts + 1)",
  main="GSE89408 | IRF7 | OA vs Healthy",cex.main=1.4,cex.lab=1.15)
points(jitter(x,amount=0.18),y,pch=21,cex=1.15,
 bg=ifelse(x==1,"#637d92","#b26a56"),col="white")
legend("topright",legend=c(
    paste("DESeq2 log2FC =",fmt(irf$log2FoldChange)),
    paste("Wald P =",fmt(irf$pvalue)),
    paste("BH-FDR =",fmt(irf$padj))),bty="n",cex=1.05)
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
