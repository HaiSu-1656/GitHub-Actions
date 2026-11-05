#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(DESeq2))

input <- "GSE254682_processed_matrix.csv.gz"
outdir <- "results/GSE254682"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(input))
data <- read.csv(gzfile(input), check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(data)==62710L,
          all(c("Gene_ID","Gene_Name","OA/Control_log2FoldChange",
                "OA/Control_pvalue","OA/Control_padj") %in% colnames(data)),
          !anyNA(data$Gene_ID), !anyDuplicated(data$Gene_ID))
cols <- c(paste0("Ctrl",1:7,"_RC"),paste0("OA",1:8,"_RC"))
stopifnot(all(cols %in% colnames(data)))
mat <- as.matrix(data[,cols,drop=FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- data$Gene_ID
stopifnot(!anyNA(mat),all(is.finite(mat)),all(mat>=0),
          all(abs(mat-round(mat))<1e-8),
          sum(rownames(mat)=="ENSG00000185507")==1L)
stopifnot(identical(data$Gene_Name[match("ENSG00000185507",data$Gene_ID)],"IRF7"))
storage.mode(mat) <- "integer"
meta <- data.frame(group=factor(c(rep("Control",7),rep("OA",8)),
                   levels=c("Control","OA")),row.names=cols)

dds <- DESeqDataSetFromMatrix(countData=mat,colData=meta,design=~group)
keep <- rowSums(counts(dds)>=10)>=3
stopifnot(keep["ENSG00000185507"])
dds <- dds[keep,]
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("group","OA","Control"),
               alpha=0.05,independentFiltering=TRUE)
result <- data.frame(Gene_ID=rownames(res), as.data.frame(res),check.names=FALSE)
result$Gene_Name <- data$Gene_Name[match(result$Gene_ID,data$Gene_ID)]
result <- result[,c("Gene_ID","Gene_Name","baseMean","log2FoldChange",
                   "lfcSE","stat","pvalue","padj")]
write.csv(result,file.path(outdir,"GSE254682_DESeq2_all_genes.csv"),
          row.names=FALSE,na="NA")
irf <- result[result$Gene_ID=="ENSG00000185507",,drop=FALSE]
stopifnot(nrow(irf)==1L)
irf$n_Control<-7L; irf$n_OA<-8L;irf$comparison<-"OA_vs_Control"
irf$FoldChange <- 2^irf$log2FoldChange
irf$log2FC_CI95_lower <- irf$log2FoldChange-1.96*irf$lfcSE
irf$log2FC_CI95_upper <- irf$log2FoldChange+1.96*irf$lfcSE
write.csv(irf,file.path(outdir,"GSE254682_IRF7_native_DESeq2_statistics.csv"),
          row.names=FALSE,na="NA")
source_row <- data[data$Gene_ID=="ENSG00000185507",,drop=FALSE]
source_vals<-data.frame(
  method=c("source_uploaded_processed_matrix","GitHub_original_DESeq2"),
  log2FoldChange=c(source_row[["OA/Control_log2FoldChange"]],irf$log2FoldChange),
  pvalue=c(source_row[["OA/Control_pvalue"]],irf$pvalue),
  padj=c(source_row[["OA/Control_padj"]],irf$padj))
write.csv(source_vals,file.path(outdir,"GSE254682_IRF7_source_vs_reanalysis.csv"),
          row.names=FALSE)
norm <- counts(dds,normalized=TRUE)["ENSG00000185507",]
sample_data <- data.frame(
  sample=sub("_RC$","",cols),
  GSM=c(paste0("GSM80487",29:35),paste0("GSM80487",36:43)),
  group=as.character(meta$group),
  raw_count=as.integer(mat["ENSG00000185507",]),
  source_normalized_count=as.numeric(source_row[1,sub("_RC$","_Norm",cols)]),
  DESeq2_normalized_count=as.numeric(norm),
  size_factor=as.numeric(sizeFactors(dds))
)
write.csv(sample_data,file.path(outdir,"GSE254682_IRF7_15_samples.csv"),
          row.names=FALSE)

fmt <- function(v) if(is.na(v)) "NA" else formatC(v,format="g",digits=11)
lines<-c("GSE254682 | 7 Control vs 8 OA | IRF7",
         paste("DESeq2_version",packageVersion("DESeq2")),
         paste("Input_Genes",nrow(mat)),
         paste("Analyzed_Genes",nrow(dds)),
         paste("All_RC_counts_integer",TRUE),
         paste("Group_design","~ group"),
         paste("Reference","Control"),
         paste("Native_IRF7_baseMean",fmt(irf$baseMean)),
         paste("Native_IRF7_log2FoldChange",fmt(irf$log2FoldChange)),
         paste("Native_IRF7_lfcSE",fmt(irf$lfcSE)),
         paste("Native_IRF7_Wald_stat",fmt(irf$stat)),
         paste("Native_IRF7_pvalue",fmt(irf$pvalue)),
         paste("Native_IRF7_padj",fmt(irf$padj)),
         paste("Native_IRF7_fold_change",fmt(irf$FoldChange)),
         paste("Native_IRF7_log2FC_CI95_lower",fmt(irf$log2FC_CI95_lower)),
         paste("Native_IRF7_log2FC_CI95_upper",fmt(irf$log2FC_CI95_upper)),
         paste("Uploaded_matrix_log2FC",fmt(source_row[["OA/Control_log2FoldChange"]])),
         paste("Uploaded_matrix_pvalue",fmt(source_row[["OA/Control_pvalue"]])),
         paste("Uploaded_matrix_padj",fmt(source_row[["OA/Control_padj"]])),
         paste("Control_mean_normalized",fmt(mean(norm[meta$group=="Control"]))),
         paste("OA_mean_normalized",fmt(mean(norm[meta$group=="OA"]))),
         paste("Control_median_normalized",fmt(median(norm[meta$group=="Control"]))),
         paste("OA_median_normalized",fmt(median(norm[meta$group=="OA"]))),
         paste("Genes_with_finite_padj",sum(!is.na(result$padj))))
writeLines(lines,file.path(outdir,"GSE254682_IRF7_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")

png(file.path(outdir,"GSE254682_IRF7_DESeq2_expression.png"),
    width=2000,height=1450,res=230)
par(family="serif")
y <- log2(norm+1)
idx <- as.integer(meta$group)
boxplot(y~factor(idx,levels=1:2),
        names=c("Control (n=7)","OA (n=8)"),
        col=c("#d6e4ef","#edd9d4"),
        outline=FALSE,boxwex=0.45,
        ylab="log2(DESeq2 normalized counts + 1)",
        main="GSE254682 | IRF7",cex.main=1.6,cex.lab=1.3,cex.axis=1.1)
set.seed(254682)
points(jitter(idx,amount=.11),y,pch=21,cex=1.5,
       bg=ifelse(idx==1,"#567996","#b96e59"),col="white")
legend("topright",bty="n",cex=0.95,
       legend=c(paste0("Wald P = ",fmt(irf$pvalue)),
                paste0("BH FDR = ",fmt(irf$padj)),
                paste0("log2FC = ",fmt(irf$log2FoldChange))))
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
