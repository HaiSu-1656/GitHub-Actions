#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(DESeq2))
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221483/suppl/GSE221483_merged_gene_name_expression2.txt.gz"
input <- "GSE221483_merged_gene_name_expression2.txt.gz"
if (!file.exists(input)) {
  cat("Downloading GEO counts:", url, "\n")
  download.file(url, input, mode="wb", method="libcurl", quiet=FALSE)
}
outdir <- "results/GSE221483"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(input), file.info(input)$size > 10000)
cat("GSE221483 md5:", unname(tools::md5sum(input)), "\n")
x <- read.delim(gzfile(input), check.names=FALSE, stringsAsFactors=FALSE)
sample_names <- c("Veh-1","Veh-2","Veh-3","Veh-4","IL1-1","IL1-2","IL1-3","IL1-4")
stopifnot(all(c("name","length",sample_names) %in% names(x)))
stopifnot(!anyNA(x$name), !anyDuplicated(x$name))
count_mat <- as.matrix(x[,sample_names,drop=FALSE])
rownames(count_mat) <- x$name
storage.mode(count_mat) <- "numeric"
stopifnot(!anyNA(count_mat), all(is.finite(count_mat)), all(count_mat>=0), all(abs(count_mat-round(count_mat)) < 1e-8))
storage.mode(count_mat) <- "integer"
stopifnot(nrow(count_mat)==54658L, ncol(count_mat)==8L)
stopifnot("Irf7" %in% rownames(count_mat))
stopifnot(identical(as.integer(count_mat["Irf7",]), c(27L,23L,23L,26L,26L,25L,26L,31L)))
condition <- factor(c(rep("Vehicle",4),rep("IL1b",4)),levels=c("Vehicle","IL1b"))
coldata <- data.frame(condition=condition,row.names=sample_names)
dds <- DESeqDataSetFromMatrix(countData=count_mat,colData=coldata,design=~condition)
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("condition","IL1b","Vehicle"),alpha=0.05,independentFiltering=TRUE)
tab <- as.data.frame(res)
tab$gene <- rownames(tab)
tab <- tab[,c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
write.csv(tab,file.path(outdir,"GSE221483_DESeq2_all_genes.csv"),row.names=FALSE,na="NA")
one <- tab[tab$gene=="Irf7",,drop=FALSE]
stopifnot(nrow(one)==1L)
one$comparison <- "IL1b_vs_Vehicle"
one$n_Vehicle <- 4L; one$n_IL1b <- 4L
write.csv(one,file.path(outdir,"GSE221483_Irf7_DESeq2_Wald_padj.csv"),row.names=FALSE,na="NA")
sam <- data.frame(sample=sample_names,condition=as.character(condition),raw_count=as.integer(count_mat["Irf7",]),normalized_count=as.numeric(counts(dds,normalized=TRUE)["Irf7",]),size_factor=as.numeric(sizeFactors(dds)))
write.csv(sam,file.path(outdir,"GSE221483_Irf7_normalized_counts.csv"),row.names=FALSE)
f <- function(v) if(is.na(v)) "NA" else formatC(v,format="g",digits=9)
lines <- c("GSE221483 native DESeq2",
 paste("DESeq2 version:",as.character(packageVersion("DESeq2"))),
 "Comparison: IL-1beta versus Vehicle, 4 biological samples per group",
 paste("Irf7 baseMean:",f(one$baseMean)),
 paste("Irf7 log2FoldChange:",f(one$log2FoldChange)),
 paste("Irf7 Wald stat:",f(one$stat)),
 paste("Irf7 Wald pvalue:",f(one$pvalue)),
 paste("Irf7 BH padj:",f(one$padj)),
 paste("Genes with non-NA adjusted P:",sum(!is.na(tab$padj))))
writeLines(lines,file.path(outdir,"Irf7_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")
write.csv(data.frame(sample=sample_names,size_factor=sizeFactors(dds)),file.path(outdir,"GSE221483_DESeq2_size_factors.csv"),row.names=FALSE)
png(file.path(outdir,"GSE221483_Irf7_barplot.png"),width=1200,height=1100,res=220)
op <- par(mar=c(5.5,5.7,2,1.5),family="serif",las=1,cex.axis=1.2,cex.lab=1.2)
y <- sam$normalized_count
means <- c(mean(y[1:4]),mean(y[5:8]))
sds <- c(sd(y[1:4]),sd(y[5:8]))
ylim <- c(0,max(means+sds)*1.33)
bp <- barplot(means,names.arg=c("Vehicle","IL-1beta"),col=c("#BED2ED","#D1C2E7"),border=c("#24539A","#5A3895"),lwd=2,ylim=ylim,ylab="DESeq2 normalized counts",xlab="",axes=FALSE)
axis(2,las=1,lwd=1.6); axis(1,at=bp,labels=c("Vehicle","IL-1beta"),tick=FALSE,lwd=0)
box(bty="l",lwd=1.6)
arrows(bp,means-sds,bp,means+sds,angle=90,code=3,length=0.055,lwd=1.8,col=c("#24539A","#5A3895"))
j <- c(-0.055,-0.018,0.018,0.055)
points(rep(bp,each=4)+rep(j,2),y,pch=21,bg=rep(c("#5C88CA","#9270C3"),each=4),col=rep(c("#24539A","#5A3895"),each=4),cex=1)
segments(bp[1],ylim[2]*0.92,bp[2],ylim[2]*0.92,lwd=1.4)
lab <- if(is.na(one$padj)) "padj = NA" else if(one$padj<0.001) "***" else if(one$padj<0.01) "**" else if(one$padj<0.05) "*" else "ns"
text(mean(bp),ylim[2]*0.96,lab,cex=1.25,font=2)
par(op); dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
