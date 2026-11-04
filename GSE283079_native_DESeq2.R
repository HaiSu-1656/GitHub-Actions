#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(DESeq2))

input <- "GSE283079_gene_count_RNAseq.txt.gz"
outdir <- "results/GSE283079"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(input))
raw <- read.delim(gzfile(input), check.names=FALSE, stringsAsFactors=FALSE,
                  row.names=1)
samples <- colnames(raw)
stopifnot(ncol(raw)==41L, !anyDuplicated(rownames(raw)),
          sum(startsWith(samples, "OA_syn_"))==36L,
          sum(startsWith(samples, "non-OA_syn_"))==5L)
group <- ifelse(startsWith(samples, "OA_syn_"),"OA",
                ifelse(startsWith(samples,"non-OA_syn_"),"NonOA",NA_character_))
stopifnot(!anyNA(group))
mat <- as.matrix(raw)
storage.mode(mat) <- "numeric"
stopifnot(!anyNA(mat), all(is.finite(mat)), all(mat >= 0))
fractional_values <- sum(abs(mat-round(mat)) > 1e-8)
fractional_genes <- sum(rowSums(abs(mat-round(mat)) > 1e-8)>0)
message("Rows: ",nrow(mat)," samples: ",ncol(mat),
        "; fractional values: ",fractional_values,
        "; genes with fractional count: ",fractional_genes)
# RSEM-style estimated gene counts are fractional; rounded counts are a
# documented approximation, NOT an exact reconstruction of a tximport
# length-corrected RSEM DESeq2 workflow.
cts <- round(mat)
storage.mode(cts) <- "integer"
meta <- data.frame(group=factor(group,levels=c("NonOA","OA")),
                   row.names=samples)
dds <- DESeqDataSetFromMatrix(countData=cts, colData=meta,design=~group)
keep <- rowSums(counts(dds)>=10)>=5
stopifnot(any(grepl("_IRF7$",rownames(dds)[keep])))
dds <- dds[keep,]
message("Genes retained: ",nrow(dds))
dds <- DESeq(dds,test="Wald",fitType="parametric",quiet=FALSE)
res <- results(dds,contrast=c("group","OA","NonOA"),
               alpha=0.05,independentFiltering=TRUE)
tab <- as.data.frame(res)
tab$gene <- rownames(tab)
tab <- tab[,c("gene","baseMean","log2FoldChange","lfcSE",
             "stat","pvalue","padj")]
write.csv(tab,file.path(outdir,"GSE283079_DESeq2_all_genes.csv"),
          row.names=FALSE,na="NA")
irf <- tab[grepl("_IRF7$",tab$gene),,drop=FALSE]
stopifnot(nrow(irf)==1)
irf$comparison <- "OA_vs_NonOA"
irf$n_OA <- 36L
irf$n_NonOA <- 5L
write.csv(irf,file.path(outdir,"GSE283079_IRF7_DESeq2_result.csv"),
          row.names=FALSE,na="NA")

irf_id <- irf$gene[1]
norm <- counts(dds,normalized=TRUE)[irf_id,]
expr <- data.frame(sample=samples,group=as.character(meta$group),
                   raw_count=as.numeric(mat[irf_id,]),
                   rounded_count=as.integer(cts[irf_id,]),
                   size_factor=as.numeric(sizeFactors(dds)),
                   DESeq2_normalized_count=as.numeric(norm))
write.csv(expr,file.path(outdir,"GSE283079_IRF7_sample_expression.csv"),
          row.names=FALSE)

fmt <- function(x) if(is.na(x)) "NA" else formatC(x,digits=9,format="g")
lines <- c("GSE283079 | IRF7 | R DESeq2 NB-GLM Wald test",
           paste("DESeq2 version",packageVersion("DESeq2")),
           "Comparison OA (n=36) vs NonOA (n=5)",
           paste("All genes",nrow(mat)),
           paste("Genes analyzed",nrow(dds)),
           paste("Fractional gene counts",fractional_values),
           paste("Fractional genes",fractional_genes),
           "Fractional RSEM-style gene counts rounded: approximation",
           paste("baseMean",fmt(irf$baseMean)),
           paste("log2FoldChange",fmt(irf$log2FoldChange)),
           paste("lfcSE",fmt(irf$lfcSE)),
           paste("Wald_stat",fmt(irf$stat)),
           paste("Wald_pvalue",fmt(irf$pvalue)),
           paste("BH_adjusted_pvalue",fmt(irf$padj)),
           paste("Genes with non-NA padj",sum(!is.na(tab$padj))))
writeLines(lines,file.path(outdir,"IRF7_DESeq2_summary.txt"))
cat(paste(lines,collapse="\n"),"\n")
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
png(file.path(outdir,"GSE283079_IRF7_DESeq2.png"),width=1800,
    height=1400,res=220)
y <- log2(expr$DESeq2_normalized_count+1)
x <- ifelse(expr$group=="NonOA",1,2)
plot(c(0.5,2.5),range(y)+c(-0.25,0.45),type="n",
     xaxt="n",xlab="",ylab="log2(DESeq2 normalized counts + 1)",
     main="GSE283079 | IRF7",cex.main=1.4)
axis(1,at=c(1,2),labels=c("Non-OA (n=5)","OA (n=36)"))
boxplot(y~factor(x,levels=1:2),add=TRUE,at=1:2,outline=FALSE,
        boxwex=0.40,col=c("#dde4ed","#f1ddd7"),border="gray40",
        axes=FALSE)
set.seed(283079)
points(jitter(x,amount=0.15),y,pch=21,
       bg=ifelse(x==1,"#63768c","#b26954"),
       cex=1.15,col="white")
legend("topright",legend=c(paste0("Wald P = ",fmt(irf$pvalue)),
                           paste0("BH padj = ",fmt(irf$padj)),
                           paste0("log2FC = ",fmt(irf$log2FoldChange))),
       bty="n")
dev.off()
