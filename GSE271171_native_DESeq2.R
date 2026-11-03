#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(DESeq2))
input <- "GSE271171_gene_count.txt.gz"
outdir <- "results"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
x <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
sample_names <- c("C1","C2","C3","I1","I2","I3")
stopifnot(all(c("gene_id","gene_name",sample_names) %in% names(x)))
stopifnot(!anyNA(x$gene_id), !anyDuplicated(x$gene_id))
counts <- as.matrix(x[,sample_names,drop=FALSE])
rownames(counts) <- x$gene_id
storage.mode(counts) <- "numeric"
stopifnot(!anyNA(counts),all(counts>=0),all(abs(counts-round(counts))<1e-8))
storage.mode(counts) <- "integer"
coldata <- data.frame(condition=factor(c(rep("Control",3),rep("IL1b",3)),levels=c("Control","IL1b")),row.names=sample_names)
dds <- DESeqDataSetFromMatrix(countData=counts,colData=coldata,design=~condition)
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("condition","IL1b","Control"),alpha=0.05,independentFiltering=TRUE)
tab <- as.data.frame(res)
tab$gene_id <- rownames(tab)
tab$gene_name <- x$gene_name[match(tab$gene_id,x$gene_id)]
tab <- tab[,c("gene_id","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
write.csv(tab,file.path(outdir,"GSE271171_DESeq2_all_genes.csv"),row.names=FALSE,na="NA")
id <- "ENSMUSG00000025498"
one <- tab[tab$gene_id==id,,drop=FALSE]
stopifnot(nrow(one)==1,one$gene_name=="Irf7")
one$comparison <- "IL1b_vs_Control"
one$n_Control <- 3L
one$n_IL1b <- 3L
write.csv(one,file.path(outdir,"GSE271171_Irf7_DESeq2_Wald_padj.csv"),row.names=FALSE,na="NA")
norm <- counts(dds,normalized=TRUE)
sam <- data.frame(sample=sample_names,condition=as.character(coldata$condition),raw_count=as.integer(counts[id,]),normalized_count=as.numeric(norm[id,]))
write.csv(sam,file.path(outdir,"GSE271171_Irf7_normalized_counts.csv"),row.names=FALSE)
f <- function(v) if (is.na(v)) "NA" else formatC(v,format="g",digits=8)
lines <- c("GSE271171 native DESeq2",paste("DESeq2 version:",as.character(packageVersion("DESeq2"))),"Comparison: IL-1beta vs Control (n=3+3)",
  paste("Irf7 baseMean:",f(one$baseMean)),paste("Irf7 log2FoldChange:",f(one$log2FoldChange)),paste("Irf7 Wald stat:",f(one$stat)),
  paste("Irf7 Wald pvalue:",f(one$pvalue)),paste("Irf7 BH padj:",f(one$padj)),paste("Genes with non-NA padj:",sum(!is.na(tab$padj))))
writeLines(lines,file.path(outdir,"Irf7_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")
png(file.path(outdir,"GSE271171_Irf7_DESeq2.png"),width=1800,height=1350,res=220)
y <- sam$normalized_count
plot(c(0.5,2.5),c(0,max(y)*1.35),type="n",xaxt="n",xlab="",ylab="DESeq2 normalized counts",main="GSE271171 | Irf7",cex.main=1.4)
axis(1,at=c(1,2),labels=c("Control","IL-1beta"))
set.seed(271171)
points(jitter(rep(1:2,each=3),amount=0.07),y,pch=21,bg=rep(c("#56667b","#B46755"),each=3),cex=1.5)
m <- tapply(y,sam$condition,mean); s <- tapply(y,sam$condition,sd)
segments(1:2,m-s,1:2,m+s,lwd=2); points(1:2,m,pch=18,cex=1.5)
legend("topright",legend=c(paste0("Wald P = ",f(one$pvalue)),paste0("BH padj = ",f(one$padj))),bty="n")
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
