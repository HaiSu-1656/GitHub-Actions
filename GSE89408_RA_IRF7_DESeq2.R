#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(DESeq2))
outdir <- "results/GSE89408_RA"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
input <- "GSE89408_GEO_count_matrix_rename.txt.gz"
dat <- read.delim(gzfile(input),row.names=1,check.names=FALSE)
stopifnot(nrow(dat)>20000,ncol(dat)==218,!anyDuplicated(rownames(dat)),"IRF7" %in% rownames(dat))
samples <- c(grep("^normal_tissue_",colnames(dat),value=TRUE),
             grep("^OA_tissue_",colnames(dat),value=TRUE),
             grep("^RA_tissue_",colnames(dat),value=TRUE))
stopifnot(length(samples)==202,
 sum(startsWith(samples,"normal_tissue_"))==28,
 sum(startsWith(samples,"OA_tissue_"))==22,
 sum(startsWith(samples,"RA_tissue_"))==152)
mat <- as.matrix(dat[,samples,drop=FALSE])
storage.mode(mat)<-"numeric"
stopifnot(!anyNA(mat),all(is.finite(mat)),all(mat>=0))
fractional <- sum(abs(mat-round(mat))>1e-8)
message("Fractional estimates: ",fractional," out of ",length(mat))
cts <- round(mat);storage.mode(cts)<-"integer"
group <- ifelse(startsWith(samples,"normal_tissue_"),"Healthy",
                ifelse(startsWith(samples,"OA_tissue_"),"OA","RA"))
meta<-data.frame(group=factor(group,levels=c("Healthy","OA","RA")),row.names=samples)
dds<-DESeqDataSetFromMatrix(countData=cts,colData=meta,design=~group)
keep<-rowSums(counts(dds)>=10)>=5
stopifnot(keep["IRF7"])
dds<-dds[keep,]
dds<-DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
contrasts<-list(RA_vs_Healthy=c("group","RA","Healthy"),
                RA_vs_OA=c("group","RA","OA"),
                OA_vs_Healthy=c("group","OA","Healthy"))
fmt<-function(v)if(length(v)==0||is.na(v))"NA" else formatC(v,format="g",digits=10)
summary<-c("GSE89408 | IRF7 | R DESeq2 v1.52 style NB GLM",
           paste("DESeq2 version",packageVersion("DESeq2")),
           paste("Total samples",length(samples)),
           paste("Analysis groups: Healthy 28, OA 22, RA 152"),
           "RA early and established pooled; stage not inferred from arbitrary sample order",
           paste("Non-integer estimates rounded: ",fractional),
           paste("Genes analyzed",nrow(dds)))
all_irf <- NULL
for (key in names(contrasts)){
  res<-results(dds,contrast=contrasts[[key]],alpha=0.05,independentFiltering=TRUE)
  all<-data.frame(gene=rownames(res),as.data.frame(res),check.names=FALSE)
  write.csv(all,file.path(outdir,paste0("GSE89408_",key,"_DESeq2_all_genes.csv")),row.names=FALSE,na="NA")
  irf<-all[all$gene=="IRF7",,drop=FALSE]
  irf$comparison<-key
  all_irf<-rbind(all_irf,irf)
  summary<-c(summary,paste(key,"baseMean",fmt(irf$baseMean),
           "log2FoldChange",fmt(irf$log2FoldChange),
           "lfcSE",fmt(irf$lfcSE),
           "Wald_stat",fmt(irf$stat),
           "pvalue",fmt(irf$pvalue),
           "padj",fmt(irf$padj),
           "FC",fmt(2^irf$log2FoldChange)))
}
write.csv(all_irf,file.path(outdir,"GSE89408_IRF7_RA_DESeq2_statistics.csv"),row.names=FALSE,na="NA")
norm<-counts(dds,normalized=TRUE)["IRF7",]
per<-data.frame(sample=samples,group=group,raw_count=as.numeric(mat["IRF7",]),
 rounded_count=as.integer(cts["IRF7",]),
 size_factor=as.numeric(sizeFactors(dds)),
 normalized_count=as.numeric(norm))
write.csv(per,file.path(outdir,"GSE89408_IRF7_202_patients.csv"),row.names=FALSE)
for (g in c("Healthy","OA","RA")){
 v<-norm[group==g]
 summary<-c(summary,paste(g,"mean_normalized_counts",fmt(mean(v)),
      "median_normalized_counts",fmt(median(v))))
}
writeLines(summary,file.path(outdir,"GSE89408_IRF7_RA_summary.txt"))
cat(paste(summary,collapse="\n"),"\n")
png(file.path(outdir,"GSE89408_IRF7_RA_3groups.png"),width=2100,height=1550,res=230)
y<-log2(norm+1)
idx<-match(group,c("Healthy","OA","RA"))
set.seed(89408)
boxplot(y~factor(idx,levels=1:3),
 names=c("Healthy (n=28)","OA (n=22)","RA (n=152)"),
 col=c("#cfdee7","#eddbce","#ddcce0"),
 outline=FALSE,ylab="log2(DESeq2 normalized counts + 1)",
 main="GSE89408 | IRF7 in human synovium",cex.main=1.4,cex.lab=1.2)
points(jitter(idx,amount=.17),y,pch=21,
 bg=c("#507995","#ba7658","#8b668f")[idx],col="white",cex=.87)
legend("topright",legend=c(
  paste0("RA vs Healthy: log2FC ",fmt(all_irf$log2FoldChange[1]),
         ", P ",fmt(all_irf$pvalue[1])),
  paste0("RA vs OA: log2FC ",fmt(all_irf$log2FoldChange[2]),
         ", P ",fmt(all_irf$pvalue[2]))),bty="n")
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
