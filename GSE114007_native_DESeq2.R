#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
outdir <- "results/GSE114007"
x <- read.delim(file.path(outdir,"all_counts.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
meta <- read.csv(file.path(outdir,"sample_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(x)==23710,nrow(meta)==38)
stopifnot(ncol(x)==39,identical(colnames(x)[-1],meta$sample),!anyDuplicated(x$symbol))
mat <- as.matrix(x[,-1]);rownames(mat) <- x$symbol;storage.mode(mat) <- "numeric"
stopifnot(!anyNA(mat),all(is.finite(mat)),all(mat>=0),all(abs(mat-round(mat))<1e-8))
storage.mode(mat) <- "integer"
rownames(meta) <- meta$sample
meta$condition <- factor(meta$condition,levels=c("Normal","OA"))
meta$platform <- factor(meta$platform,levels=c("HiSeq2000","NextSeq500"))
stopifnot(identical(colnames(mat),rownames(meta)))
stopifnot(identical(as.integer(table(meta$condition)),c(18L,20L)))
print(table(meta$platform,meta$condition))
dds <- DESeqDataSetFromMatrix(countData=mat,colData=meta,design=~platform+condition)
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("condition","OA","Normal"),alpha=0.05,independentFiltering=TRUE)
tab <- as.data.frame(res)
tab$gene <- rownames(tab)
tab <- tab[,c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
write.csv(tab,file.path(outdir,"GSE114007_DESeq2_platform_adjusted_all_genes.csv"),row.names=FALSE)
one <- tab[tab$gene=="IRF7",,drop=FALSE]
stopifnot(nrow(one)==1)
write.csv(one,file.path(outdir,"GSE114007_IRF7_DESeq2_Wald_padj.csv"),row.names=FALSE)
sam <- data.frame(sample=meta$sample,condition=meta$condition,platform=meta$platform,
                  raw_count=as.integer(mat["IRF7",]),normalized_count=as.numeric(counts(dds,normalized=TRUE)["IRF7",]),
                  size_factor=sizeFactors(dds))
write.csv(sam,file.path(outdir,"GSE114007_IRF7_per_sample.csv"),row.names=FALSE)
f <- function(z) if (is.na(z)) "NA" else formatC(z,format="g",digits=9)
lines <- c("GSE114007 native DESeq2","Model: ~ platform + condition",
paste("DESeq2 version:",as.character(packageVersion("DESeq2"))),
paste("IRF7 baseMean:",f(one$baseMean)),
paste("IRF7 log2FC OA_vs_Normal:",f(one$log2FoldChange)),
paste("IRF7 Wald stat:",f(one$stat)),
paste("IRF7 Wald pvalue:",f(one$pvalue)),
paste("IRF7 BH padj:",f(one$padj)),
paste("Genes with non-NA padj:",sum(!is.na(tab$padj))))
writeLines(lines,file.path(outdir,"IRF7_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")
grp <- factor(sam$condition,levels=c("Normal","OA"))
mean_sd <- function(z) c(mean(z),sd(z))
sapply(split(sam$normalized_count,grp),mean_sd)
palette_fill <- c("#BDD0EB","#CFC0E4")
palette_stroke <- c("#22529A","#58388F")
png(file.path(outdir,"GSE114007_IRF7_expression.png"),width=1900,height=1500,res=300)
par(family="serif",mar=c(5,6,2.2,1),mgp=c(3.2,0.8,0),las=1)
meanVals <- c(mean(sam$normalized_count[grp=="Normal"]),mean(sam$normalized_count[grp=="OA"]))
sdVals <- c(sd(sam$normalized_count[grp=="Normal"]),sd(sam$normalized_count[grp=="OA"]))
ym <- max(sam$normalized_count,meanVals+sdVals)*1.30
bp <- barplot(meanVals,col=palette_fill,border=palette_stroke,lwd=2.3,width=0.6,space=0.5,
              ylim=c(0,ym),axes=FALSE,names.arg=c("Normal","OA"),
              cex.names=1.4,ylab="IRF7 normalized counts",cex.lab=1.3)
axis(2,at=pretty(c(0,ym),n=6),las=1,cex.axis=1.15,lwd=1.8)
axis(1,at=bp,labels=FALSE,tick=FALSE)
box(bty="l",lwd=2)
arrows(bp,meanVals-sdVals,bp,meanVals+sdVals,angle=90,code=3,length=0.045,lwd=1.8,col=palette_stroke)
for(i in seq_along(bp)){
  z <- sam$normalized_count[grp==levels(grp)[i]]
  set.seed(114007+i)
  points(bp[i]+jitter(rep(0,length(z)),amount=0.085),z,
         pch=21,bg=adjustcolor(palette_stroke[i],alpha.f=0.6),col=palette_stroke[i],cex=0.95)
}
segments(bp[1],ym*0.91,bp[2],ym*0.91,lwd=1.7)
sig <- if(is.na(one$padj)) "padj = NA" else if(one$padj<0.001) "***" else if(one$padj<0.01) "**" else if(one$padj<0.05) "*" else "ns"
text(mean(bp),ym*0.95,sig,font=2,cex=1.4)
dev.off()
# Platform-stratified descriptive statistics and chart
z <- aggregate(normalized_count~platform+condition,data=sam,FUN=function(v)c(mean=mean(v),sd=sd(v),n=length(v)))
write.csv(data.frame(platform=z$platform,condition=z$condition,z$normalized_count),file.path(outdir,"GSE114007_platform_stratified_means.csv"),row.names=FALSE)
png(file.path(outdir,"GSE114007_IRF7_by_platform.png"),width=1800,height=1300,res=300)
par(family="serif",mar=c(4.5,5.5,2,1))
boxplot(normalized_count~interaction(platform,condition,sep="\n"),data=sam,
        col=rep(palette_fill,each=2),outline=FALSE,ylab="IRF7 normalized counts",
        xlab="",cex.axis=1.08,cex.lab=1.2)
set.seed(114007);stripchart(normalized_count~interaction(platform,condition,sep="\n"),data=sam,
vertical=TRUE,method="jitter",jitter=0.12,add=TRUE,pch=16,col="#4B5060",cex=.7)
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
