#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, timeout=1200)
suppressPackageStartupMessages({library(GEOquery); library(Biobase); library(limma)})
outdir <- "results/GSE55235"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
cat("Reading GSE55235 series matrix from GEO\n")
g <- GEOquery::getGEO("GSE55235",GSEMatrix=TRUE,getGPL=FALSE,destdir=".")
stopifnot(length(g)==1)
es <- g[[1]]; raw <- exprs(es); pd <- pData(es)
stopifnot(ncol(raw)==30, nrow(raw)>=20000)
titles <- as.character(pd$title)
typ <- ifelse(grepl("ND_[0-9]+",titles),"Normal",ifelse(grepl("OA_[0-9]+",titles),"OA",ifelse(grepl("RA_[0-9]+",titles),"RA","Unknown")))
cat("Groups:",paste(names(table(typ)),as.integer(table(typ)),collapse=", "),"\n")
stopifnot(identical(as.integer(table(factor(typ,levels=c("Normal","OA","RA")))),c(10L,10L,10L)))
meta <- data.frame(sample=colnames(raw),title=titles,group=typ,stringsAsFactors=FALSE)
write.csv(meta,file.path(outdir,"GSE55235_all_samples_metadata.csv"),row.names=FALSE)
keep <- typ %in% c("Normal","OA")
data <- raw[,keep,drop=FALSE]
groups <- factor(typ[keep],levels=c("Normal","OA"))
stopifnot(ncol(data)==20,all(is.finite(data)),all(!is.na(data)))
cat("Raw distribution:"); print(quantile(as.numeric(data),c(0,.25,.5,.75,.95,.99,1),na.rm=TRUE))
was_logged <- quantile(as.numeric(data),.99,na.rm=TRUE)<25
logmas <- if(was_logged) data else log2(pmax(data,1))
norm <- limma::normalizeBetweenArrays(logmas,method="quantile")
cat("Input already log-like:",was_logged,"; quantile-normalized expression range:",paste(round(range(norm),3),collapse=", "),"\n")
cat("Annotation GPL96\n")
gpl <- GEOquery::getGEO("GPL96",destdir=".",AnnotGPL=FALSE)
a <- GEOquery::Table(gpl)
cat("GPL96 columns:",paste(names(a),collapse="; "),"\n")
id_col <- intersect(c("ID","ID_REF"),names(a))[1]
symbol_col <- intersect(c("Gene Symbol","Gene symbol","Symbol","GENE_SYMBOL","gene_symbol"),names(a))[1]
stopifnot(!is.na(id_col),!is.na(symbol_col))
sym <- as.character(a[[symbol_col]])
# Exact token match even for multiple gene aliases split by ///.
is_irf7 <- grepl("(^|[^[:alnum:]_])IRF7([^[:alnum:]_]|$)",sym,perl=TRUE)
probes <- intersect(unique(as.character(a[[id_col]][is_irf7])),rownames(norm))
cat("IRF7 probe IDs:",paste(probes,collapse=", "),"\n")
stopifnot(length(probes)>0)
write.csv(a[as.character(a[[id_col]])%in%probes,,drop=FALSE],file.path(outdir,"GSE55235_IRF7_GPL96_annotation.csv"),row.names=FALSE)
design <- model.matrix(~groups)
fit <- limma::eBayes(limma::lmFit(norm,design))
tt <- limma::topTable(fit,coef="groupsOA",number=Inf,sort.by="none",adjust.method="BH")
tt$probe_id <- rownames(tt)
write.csv(tt,file.path(outdir,"GSE55235_limma_all_22283_probes.csv"),row.names=FALSE)
res <- tt[match(probes,tt$probe_id),,drop=FALSE]
res$mean_Normal <- rowMeans(norm[probes,groups=="Normal",drop=FALSE])
res$mean_OA <- rowMeans(norm[probes,groups=="OA",drop=FALSE])
res$mean_overall <- rowMeans(norm[probes,,drop=FALSE])
res <- res[order(-res$mean_overall),,drop=FALSE]
res$primary_probe <- seq_len(nrow(res))==1L
write.csv(res,file.path(outdir,"GSE55235_IRF7_all_probes_limma_BH.csv"),row.names=FALSE)
cat("Full IRF7 probes result:\n")
print(res)
pr <- res$probe_id[1]
sample_out <- data.frame(sample=colnames(norm),group=as.character(groups),probe_id=pr,raw_signal=as.numeric(data[pr,]),log2_MAS5=as.numeric(logmas[pr,]),normalized_log2_signal=as.numeric(norm[pr,]),stringsAsFactors=FALSE)
write.csv(sample_out,file.path(outdir,"GSE55235_IRF7_20_samples.csv"),row.names=FALSE)
# Sensitivity: log2 MAS5 without quantile normalization, same samples
fit2 <- limma::eBayes(limma::lmFit(logmas,design))
tt2 <- limma::topTable(fit2,coef="groupsOA",number=Inf,sort.by="none",adjust.method="BH")
tt2$probe_id<-rownames(tt2)
sens <- tt2[match(res$probe_id,tt2$probe_id),c("probe_id","logFC","P.Value","adj.P.Val"),drop=FALSE]
write.csv(sens,file.path(outdir,"GSE55235_IRF7_no_quantile_sensitivity.csv"),row.names=FALSE)
cat("Sensitivity (log2 MAS5 without between-array quantile normalization):\n");print(sens)
pn <- sample_out$normalized_log2_signal[groups=="Normal"]
po <- sample_out$normalized_log2_signal[groups=="OA"]
summary <- c("GSE55235 (GPL96): synovium, 10 OA and 10 Healthy (10 RA excluded)",
 paste("R version:",R.version.string,"limma:",as.character(packageVersion("limma"))),
 paste("Raw GEO matrix originally log-like:",was_logged),
 "Method: log2 MAS5 signal + quantile normalization + limma, BH adjusted across all probes",
 paste("Primary probe:",pr), paste("All IRF7 probes:",paste(probes,collapse="; ")),
 paste("Mean Normal:",format(mean(pn),digits=9),"SD:",format(sd(pn),digits=9)),
 paste("Mean OA:",format(mean(po),digits=9),"SD:",format(sd(po),digits=9)),
 paste("Primary logFC (OA-Normal):",format(res$logFC[1],digits=9)),
 paste("Primary moderated P:",format(res$P.Value[1],digits=9)),
 paste("Primary BH padj:",format(res$adj.P.Val[1],digits=9)),
 paste("Genes/probes with nonNA padj:",sum(!is.na(tt$adj.P.Val))),
 paste("GEO raw probe rows:",nrow(raw))
)
cat(paste(summary,collapse="\n"),"\n")
writeLines(summary,file.path(outdir,"GSE55235_IRF7_summary.txt"))
write.csv(data.frame(probe_id=pr,normal_n=length(pn),OA_n=length(po),normal_mean=mean(pn),normal_sd=sd(pn),OA_mean=mean(po),OA_sd=sd(po),logFC=res$logFC[1],pvalue=res$P.Value[1],padj=res$adj.P.Val[1]),file.path(outdir,"GSE55235_IRF7_summary.csv"),row.names=FALSE)
# Real-data publication plot
png(file.path(outdir,"GSE55235_IRF7_barplot.png"),width=1700,height=1500,res=300)
par(family="serif",mar=c(5.3,5.9,2,1.1),las=1)
mu <- c(mean(pn),mean(po)); sdv<-c(sd(pn),sd(po)); rge <- range(c(pn,po))
ylim <- c(min(0, floor(rge[1]-0.5)),ceiling(rge[2]+1.3))
bp <- barplot(mu,names.arg=c("Normal","OA"),col=c("#BFD2ED","#D4C5EA"),border=c("#24549D","#58388E"),lwd=2.3,ylim=ylim,axes=FALSE,ylab="IRF7 normalized expression",cex.lab=1.25,cex.names=1.5)
axis(2,las=1,lwd=1.8,cex.axis=1.15);box(bty="l",lwd=2)
arrows(bp,mu-sdv,bp,mu+sdv,angle=90,code=3,length=.06,lwd=1.7,col=c("#24549D","#58388E"))
for(i in seq_along(bp)) {
 v <- if(i==1)pn else po
 set.seed(55235+i)
 points(bp[i]+jitter(rep(0,length(v)),amount=.085),v,pch=21,bg=if(i==1)"#6D97D5" else "#9E85D2",col=if(i==1)"#24549D" else "#58388E",cex=.92)
}
pval<-res$adj.P.Val[1]; lab<-if(is.na(pval))"padj=NA" else if(pval<.001)"***" else if(pval<.01)"**" else if(pval<.05)"*" else "ns"
segments(bp[1],ylim[2]-0.6,bp[2],ylim[2]-0.6,lwd=1.7);text(mean(bp),ylim[2]-.25,lab,cex=1.4,font=2)
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
