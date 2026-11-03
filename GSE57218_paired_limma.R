#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,timeout=900)
suppressPackageStartupMessages({library(GEOquery);library(Biobase);library(limma)})
outdir <- "results/GSE57218"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
cat("Downloading GSE57218 series matrix and GPL6947 annotation...\n")
gse <- GEOquery::getGEO("GSE57218",GSEMatrix=TRUE,getGPL=FALSE,destdir=".",AnnotGPL=FALSE)
stopifnot(length(gse)>=1L)
eset <- gse[[1]]
emat <- exprs(eset); pdata <- pData(eset)
stopifnot(ncol(emat)==nrow(pdata),ncol(emat)==73L,nrow(emat)>40000)
title <- as.character(pdata$title)
matched <- regexec("^cartilage_RAAK_([0-9]+)_(OA|Preserved|Healthy)$",title)
pieces <- regmatches(title,matched)
bad <- which(lengths(pieces)!=3L)
if(length(bad))stop("Unrecognized sample titles: ",paste(title[bad],collapse="; "))
patient <- vapply(pieces,function(p)p[2],character(1))
type <- vapply(pieces,function(p)p[3],character(1))
meta <- data.frame(gsm=colnames(emat),title=title,patient=patient,type=type,stringsAsFactors=FALSE)
write.csv(meta,file.path(outdir,"GSE57218_all_73_samples_metadata.csv"),row.names=FALSE)
pairs <- sort(unique(patient[type=="OA"][patient[type=="OA"]%in%patient[type=="Preserved"]]))
cat("All sample types: ");print(table(type))
cat("Matched OA and preserved donors:",length(pairs),"\n")
stopifnot(length(pairs)==33,all(table(patient[type%in%c("OA","Preserved")])==2))
keep <- type%in%c("OA","Preserved") & patient%in%pairs
sub <- emat[,keep,drop=FALSE]
m <- meta[keep,,drop=FALSE]
stopifnot(ncol(sub)==66,nrow(m)==66)
m$patient <- factor(m$patient)
m$condition <- factor(m$type,levels=c("Preserved","OA"))
design <- model.matrix(~patient+condition,data=m)
stopifnot("conditionOA"%in%colnames(design),qr(design)$rank==ncol(design))
cat("Downloading GPL6947 probe annotation...\n")
gpl <- GEOquery::getGEO("GPL6947",destdir=".",AnnotGPL=FALSE)
ann <- GEOquery::Table(gpl)
write.csv(data.frame(annotation_column=names(ann)),file.path(outdir,"GPL6947_annotation_columns.csv"),row.names=FALSE)
cat("GPL annotation columns:",paste(names(ann),collapse="; "),"\n")
idcol <- intersect(c("ID","ID_REF","ProbeID","PROBE_ID"),names(ann))[1]
stopifnot(!is.na(idcol))
candidates <- intersect(c("SYMBOL","Symbol","symbol","ILMN_GENE","Gene Symbol","Gene.symbol","GENE_SYMBOL"),names(ann))
stopifnot(length(candidates)>0)
cat("Symbol columns considered:",paste(candidates,collapse=", "),"\n")
target <- rep(FALSE,nrow(ann))
for (sc in candidates) {
 a <- trimws(as.character(ann[[sc]]))
 # Match exact symbols only; avoid matching IRF7-like family members or free-text descriptions
 target <- target | a=="IRF7"
}
probeids <- unique(as.character(ann[[idcol]][target]))
probeids <- intersect(probeids,rownames(sub))
cat("GPL6947 IRF7 probe IDs:",paste(probeids,collapse=", "),"\n")
stopifnot(length(probeids)>=1)
ann_irf7 <- ann[as.character(ann[[idcol]])%in%probeids,,drop=FALSE]
write.csv(ann_irf7,file.path(outdir,"GSE57218_IRF7_probe_annotation.csv"),row.names=FALSE)
fit <- limma::eBayes(limma::lmFit(sub,design))
tt <- limma::topTable(fit,coef="conditionOA",number=Inf,sort.by="none",adjust.method="BH",confint=TRUE)
tt$probe_id <- rownames(tt)
write.csv(tt,file.path(outdir,"GSE57218_limma_all_probes.csv"),row.names=FALSE)
one <- tt[match(probeids,tt$probe_id),,drop=FALSE]
one$mean_Preserved <- vapply(probeids,function(probe)mean(sub[probe,m$condition=="Preserved"]),numeric(1))
one$mean_OA <- vapply(probeids,function(probe)mean(sub[probe,m$condition=="OA"]),numeric(1))
one$mean_all_samples <- rowMeans(sub[probeids,,drop=FALSE])
one <- one[order(-one$mean_all_samples),,drop=FALSE]
one$is_primary_probe <- seq_len(nrow(one))==1
write.csv(one,file.path(outdir,"GSE57218_IRF7_all_probes_paired_limma.csv"),row.names=FALSE)
primary <- one$probe_id[1]; stopifnot(!is.na(primary))
cat("PRIMARY IRF7 probe:",primary,"\n")
cat("IRF7 probe limma results:\n");print(one)
ind <- which(m$condition=="Preserved")
oa <- which(m$condition=="OA")
ids <- as.character(m$patient[ind]); 
w <- match(ids,as.character(m$patient[oa]))
stopifnot(length(w)==33,all(!is.na(w)))
pr <- as.numeric(sub[primary,ind])
ov <- as.numeric(sub[primary,oa[w]])
paired <- data.frame(patient=ids, preserved_gsm=m$gsm[ind],OA_gsm=m$gsm[oa[w]],preserved_expression=pr,OA_expression=ov,paired_delta=ov-pr,stringsAsFactors=FALSE)
write.csv(paired,file.path(outdir,"GSE57218_IRF7_33_paired_samples.csv"),row.names=FALSE)
delta <- ov-pr
tr <- t.test(ov,pr,paired=TRUE)
wil <- wilcox.test(ov,pr,paired=TRUE,exact=FALSE)
summary <- c("GSE57218 GPL6947, 73 total GSM samples, 33 matched pairs (n=66), 7 healthy excluded from pairwise test",
 paste("IRF7 probe IDs:",paste(probeids,collapse=", ")),
 paste("Primary probe:",primary),
 paste("Primary-probe mean preserved:",format(mean(pr),digits=9)),
 paste("Primary-probe mean OA:",format(mean(ov),digits=9)),
 paste("Mean paired difference OA-Preserved:",format(mean(delta),digits=9)),
 paste("Limma paired model logFC:",format(one$logFC[1],digits=9)),
 paste("Limma moderated P.value:",format(one$P.Value[1],digits=9)),
 paste("Limma BH adjusted P.value:",format(one$adj.P.Val[1],digits=9)),
 paste("Paired t-test P.value:",format(tr$p.value,digits=9)),
 paste("Wilcoxon signed-rank P.value:",format(wil$p.value,digits=9)),
 paste("Samples increasing IRF7:",sum(delta>0),"decreasing:",sum(delta<0),"zero:",sum(delta==0)),
 paste("Expression data range:",paste(range(sub,na.rm=TRUE),collapse=" to ")))
cat(paste(summary,collapse="\n"),"\n")
writeLines(summary,file.path(outdir,"GSE57218_IRF7_summary.txt"))
write.csv(data.frame(primary_probe=primary,normal_mean=mean(pr),OA_mean=mean(ov),mean_delta=mean(delta),
limma_pvalue=one$P.Value[1],limma_padj=one$adj.P.Val[1],paired_t_pvalue=tr$p.value,paired_wilcoxon_pvalue=wil$p.value),file.path(outdir,"GSE57218_IRF7_summary.csv"),row.names=FALSE)
# Exact-data figure, 33 paired line segments, publication serif style
png(file.path(outdir,"GSE57218_IRF7_paired.png"),width=1900,height=1550,res=300)
par(mar=c(5.8,6.1,2.4,1.1),family="serif",las=1)
yl <- range(c(pr,ov)); padx <- diff(yl)*0.07
plot(NA,xlim=c(.7,2.3),ylim=c(yl[1]-padx,yl[2]+padx*2.6),xaxt="n",xlab="",ylab="IRF7 expression (normalized array intensity)",cex.lab=1.25,cex.axis=1.05,bty="l",lwd=2)
axis(1,at=1:2,labels=c("Preserved","OA"),tick=FALSE,cex.axis=1.25)
segments(rep(1,33),pr,rep(2,33),ov,col="#B0B1BE",lwd=1.1)
points(rep(1,33),pr,pch=21,bg="#BDD0ED",col="#245598",cex=.85)
points(rep(2,33),ov,pch=21,bg="#CFBEE9",col="#58368D",cex=.85)
segments(1,mean(pr),2,mean(ov),lty=2,lwd=1.8,col="#2E2F34")
segments(.86,mean(pr),1.14,mean(pr),lwd=3.5,col="#245598")
segments(1.86,mean(ov),2.14,mean(ov),lwd=3.5,col="#58368D")
sig <- if(is.na(one$adj.P.Val[1])) "padj = NA" else if(one$adj.P.Val[1]<.001) "***" else if(one$adj.P.Val[1]<.01) "**" else if(one$adj.P.Val[1]<.05) "*" else "ns"
segments(1,yl[2]+padx*1.2,2,yl[2]+padx*1.2,lwd=1.7)
text(1.5,yl[2]+padx*1.7,sig,font=2,cex=1.35)
dev.off()
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
