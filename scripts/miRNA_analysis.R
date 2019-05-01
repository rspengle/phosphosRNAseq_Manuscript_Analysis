# Analyzes miRNA data from healthy control plasma samples. Used for EMBO manuscript Fig EV1. 

library(data.table)
library(DESeq2)
library(tidyr)
library(rtracklayer)
library(ggplot2)
#setDTthreads(18)
contam.dirs <- "./data/Plasma/healthyControlsAnalysis/biofragmenta-v1.5_20180220_103636/02-map_contaminates"
sample.info.rdata.file <- "./data/Plasma/healthyControlsAnalysis/ALL.SAMPLE.INFO"
load(sample.info.rdata.file)
# Load ULMC Control contaminate data, which also has miRNA and other useful RNA data ----
contam.files <- list.files(contam.dirs, pattern="*Processed.feature", full.names=TRUE, recursive=TRUE)

contam.files.dt <- rbindlist(lapply(contam.files, 
                                    FUN=function(x){
                                      dt <- fread(x, sep="\t", fill = TRUE)
                                      nm <- sub("_[ATGC]*_l001_l004_l007_Processed.feature", "", basename(x))
                                      dt[, File.Base.ID2:=nm]
                                      return(dt)
                                    }), fill=TRUE)

contam.files.dt[ALL.SAMPLE.INFO, `:=`(PNK=i.PNK, participant.ID=i.participant.ID, replicate.n=i.replicate.n), on="File.Base.ID2"]
contam.files.dt[, c("db", "mis.orient"):=tstrsplit(V1[1], split=".mis_", fixed=TRUE), by=V1]
contam.files.dt[, c("n.mm", "orient"):=tstrsplit(mis.orient[1], split="", type.convert = TRUE), by=mis.orient]

contam.files.dt.miR <- copy(contam.files.dt[db=="human_miRNA"])
contam.files.dt.miR.1mm <- contam.files.dt.miR[orient=="+" & n.mm<2, .(count=sum(V3)), by=.(PNK, participant.ID, replicate.n, V2)]
setnames(contam.files.dt.miR.1mm, "V2", "miR.ID")
contam.files.dt.miR.1mm[count>1, .(sample.count=uniqueN(participant.ID)), by=.(miR.ID, PNK)][sample.count==5, .N, by=PNK]

g1 <- ggplot(contam.files.dt.miR.1mm[, .(sum(count)), by=.(participant.ID, PNK)], aes(x=PNK, y=V1)) +
  geom_boxplot() + 
  ggbeeswarm::geom_quasirandom() + 
  scale_y_continuous(breaks=seq(0,5e6,500000)) + labs(y="Total miRNA aligned (count)", x=NULL) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size=8, color="black")); g1
fig.ex1c.dat <- contam.files.dt.miR.1mm[, .(total.miRNA.aligned.counts=sum(count)), by=.(participant.ID, PNK)]

save_plot(filename = "./output/figures/Supplement_miRNA_counts_PNK_VS_NoPNK.pdf", plot = g1, ncol = 1, nrow = 1, base_height = 3.5, base_width = 3.5)
save_plot(filename = "./output/figures/FigEX1A.pdf", plot = g1, ncol = 1, nrow = 1, base_height = 3.5, base_width = 3.5)
fwrite(x = fig.ex1c.dat, file="./output/figures/source_data/FigEX1A.csv", sep=",")

g2 <- ggplot(contam.files.dt.miR.1mm[count>1, .(sample.count=uniqueN(participant.ID)), by=.(miR.ID, PNK)][sample.count==5, .N, by=PNK], aes(x=PNK, y=N)) + geom_bar(stat="identity") + labs(y="# miRs detected in all samples", x=NULL) + scale_y_continuous(breaks=seq(0,300,50)) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size=8, color="black")); g2

save_plot(filename = "./output/figures/Supplement_n_miRNAs_Detected_PNK_VS_NoPNK.pdf", plot = g2, ncol = 1, nrow = 1, base_height = 3.5, base_width = 3.5)
save_plot(filename = "./output/figures/Fig_EX1B.pdf", plot = g2, ncol = 1, nrow = 1, base_height = 3.5, base_width = 3.5)
fig.ex1b.dat <- contam.files.dt.miR.1mm[, .(number.miRs.detected=sum(count)), by=.(participant.ID, PNK)]
fwrite(x = fig.ex1b.dat, file="./output/figures/source_data/FigEX1B.csv", sep=",")

setnames(contam.files.dt.miR.1mm, "V2", "miR.ID")
fwrite(x = contam.files.dt.miR.1mm, "./output/intermediate_files/miRNA_Analysis/miR_counts_dt_1mm.txt")
miR.counts <- data.frame(dcast.data.table(contam.files.dt.miR.1mm, miR.ID~participant.ID+PNK+replicate.n, value.var = "count", fun.aggregate = sum, fill=0), row.names=1)
miR.counts.f <- miR.counts[rowSums(miR.counts>0)>1,]
sample.info.df <- data.frame(contam.files.dt.miR.1mm[, .N, by=.(participant.ID, PNK, replicate.n)][, sample.ID:=paste(participant.ID, PNK, replicate.n, sep="_"), by=.(participant.ID, PNK, replicate.n)][], row.names="sample.ID")
sample.info.df <- sample.info.df[colnames(miR.counts),]
design <- model.matrix(~participant.ID+PNK, data=sample.info.df)
miR.deseq <- DESeqDataSetFromMatrix(countData=miR.counts.f, colData=DataFrame(sample.info.df), design=~participant.ID+PNK)
miR.deseq <- DESeq(miR.deseq, test="LRT", fitType="local", sfType="poscounts", reduced=~participant.ID)
miR.deseq.res <- results(miR.deseq, contrast = c("PNK", "PNK", "NoPNK"))
miR.deseq.res.dt <- data.table(data.frame(miR.deseq.res), keep.rownames = TRUE); miR.deseq.res.dt
miR.deseq.res.dt.sig <- miR.deseq.res.dt[!is.na(padj) & padj<0.1]

contam.files.dt.miR.1mm.sample <- contam.files.dt.miR.1mm[, .(tot.count.)]
sample.counts.miRs <- contam.files.dt.miR.1mm[count>=1, .(tot.count=sum(count), n.participants=uniqueN(participant.ID), min.count=min(count, na.rm = TRUE), max.count=max(count, na.rm=TRUE)), by=.(PNK, miR.ID)]

sample.counts.miRs.cast <- dcast.data.table(sample.counts.miRs, miR.ID~PNK, value.var = c("tot.count", "n.participants", "min.count", "max.count"), fun.aggregate = sum, fill=0)
sample.counts.miRs.cast[, `:=`(min.count_NoPNK=ifelse(n.participants_NoPNK<5, 0, min.count_NoPNK),
                               min.count_PNK=ifelse(n.participants_PNK<5, 0, min.count_PNK))]
sample.counts.miRs.cast[, .N, by=.(n.participants_NoPNK, n.participants_PNK)][(n.participants_NoPNK==5 & n.participants_PNK==0) | (n.participants_NoPNK==0 & n.participants_PNK==5)]


sample.counts.missing.in.one <- sample.counts.miRs.cast[(n.participants_NoPNK==5 & n.participants_PNK==0) | (n.participants_NoPNK==0 & n.participants_PNK==5)]
sample.counts.missing.in.one[, miR.ID.noArm:=sub("-[35]p$", "", miR.ID)]
setorderv(sample.counts.missing.in.one, c("tot.count_PNK", "tot.count_NoPNK"), c(-1, -1))

mir.gff <- rtracklayer::readGFF("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3", version=3)
miR.sRNAnalzyer.info <- fread("/data/LargeApplications/sRNAnalyzer/MainDBs/miRBase/unifiedName/hairpin_hsa_all_anno.description", sep="\t", header=FALSE)
miR.sRNAnalzyer.info[, miR.ID:=tolower(V1)]
miR.sRNAnalzyer.info[, miR.acc:=tstrsplit(V2, split=" ", fixed=TRUE, keep=1)]
aliases.fl <- "./data/miRNA_aliases.txt.gz"
aliases.dt <- fread(cmd=paste0("zcat ", aliases.fl), header = FALSE)
mir.gff.dt <- data.table(as.data.frame(mir.gff))
mir.gff.dt[, Name:=tolower(Name)]
setnames(miR.deseq.res.dt.sig, "rn", "miR.ID")
miR.deseq.res.dt.sig[, miR.ID:=tolower(miR.ID)]
miR.deseq.res.dt.sig[miR.sRNAnalzyer.info, info:=i.V2, on="miR.ID"]
miR.deseq.res.dt.sig[, miR.acc:=tstrsplit(info, split=" ", fixed=TRUE, keep=1)]
miR.deseq.res.dt.sig[, miR.id.noArm:=sub("-[35]p$", "", miR.ID)]
miR.deseq.res.dt.sig[mir.gff.dt, `:=`(seqid=i.seqid, start=i.start, end=i.end, type=i.type), on=c("miR.acc"="ID")]


repeat.file <- "./data/annotations/hg38_miscRNA_and_rmsk.bed"
repeats.dt <- fread(repeat.file)
setnames(repeats.dt, c("seqnames", "start", "end", "name", "score", "strand"))
setkey(repeats.dt, seqnames, start, end)
setkey(miR.deseq.res.dt.sig, seqid, start, end)
miR.deseq.with.repeats <- foverlaps(miR.deseq.res.dt.sig[!is.na(seqid)], repeats.dt, mult="first", nomatch=0)

miR.deseq.res.dt.sig[miR.deseq.with.repeats, repeat.name:=i.name, on="miR.acc"]
miR.deseq.res.dt.sig[, miR.missing.or.repeat:=ifelse(!is.na(repeat.name), "REPEAT", ifelse(is.na(seqid), "MISSING", "miR"))]
miR.deseq.res.dt.sig[, miR.change:=ifelse(log2FoldChange>0, "PNK_Higher", "NoPNK_Higher")]
#ggplot(miR.deseq.res.dt.sig[, .N, by=.(miR.missing.or.repeat, miR.change)], aes(x=miR.change, fill=miR.missing.or.repeat, y=N)) + geom_bar(stat="identity", pos="dodge")

sample.counts.miRs.cast[, miR.ID:=tolower(miR.ID)]
sample.counts.miRs.cast[miR.sRNAnalzyer.info, miR.acc:=i.miR.acc, on="miR.ID"]
sample.counts.miRs.cast[mir.gff.dt, `:=`(seqid=i.seqid, start=i.start, end=i.end, type=i.type), on=c("miR.acc"="ID")]
setkey(sample.counts.miRs.cast, seqid, start, end)
sample.counts.miRs.cast.with.repeats <-  foverlaps(sample.counts.miRs.cast[!is.na(seqid)], repeats.dt, mult="first", nomatch=0)
sample.counts.miRs.cast[sample.counts.miRs.cast.with.repeats, repeat.name:=i.name, on="miR.acc"]
sample.counts.miRs.cast[, miR.missing.or.repeat:=ifelse(!is.na(repeat.name), "REPEAT", ifelse(is.na(seqid), "MISSING", "miR"))]
sample.counts.miRs.cast[, miR.change:=ifelse(n.participants_PNK>n.participants_NoPNK, "PNK_Higher", ifelse(n.participants_PNK==n.participants_NoPNK, "BOTH", "NoPNK_Higher"))]

sample.counts.missing.in.one <- sample.counts.miRs.cast[(n.participants_NoPNK==5 & n.participants_PNK==0) | (n.participants_NoPNK==0 & n.participants_PNK==5)]

ggplot(sample.counts.missing.in.one[,  uniqueN(miR.ID), by=.(miR.missing.or.repeat, miR.change)], aes(x=miR.change, fill=miR.missing.or.repeat, y=V1)) + geom_bar(stat="identity", pos="dodge")

contam.files.dt.miR.1mm[, miR.ID:=tolower(miR.ID)]
contam.files.dt.miR.1mm[sample.counts.miRs.cast, miR.missing.or.repeat:=i.miR.missing.or.repeat, on="miR.ID"]
contam.files.dt.miR.1mm.sample.sum <- contam.files.dt.miR.1mm[, .(count=sum(count)), by=.(PNK, participant.ID, miR.missing.or.repeat, miR.ID)]
contam.files.dt.miR.1mm.sample.sum[, sample.tot:=sum(count), by=.(participant.ID, PNK)]

median.sample.tot <- contam.files.dt.miR.1mm.sample.sum[, .N, by=.(PNK, participant.ID, sample.tot)][, .(median.sample.tot=median(sample.tot)), by=PNK]
contam.files.dt.miR.1mm.sample.sum[median.sample.tot, cpm:=(i.median.sample.tot*count)/sample.tot, on="PNK"]
contam.files.dt.miR.1mm.sample.sum[median.sample.tot, cpm.million:=(10^6*count)/sample.tot, on="PNK"]

miR.order <- contam.files.dt.miR.1mm.sample.sum[, .(cpm=median(cpm)), by=.(miR.ID, miR.missing.or.repeat)][order(miR.missing.or.repeat, -cpm),]$miR.ID
contam.files.dt.miR.1mm.sample.sum[, miR.ID.f:=factor(miR.ID, levels=miR.order)]
g3 <- ggplot(contam.files.dt.miR.1mm.sample.sum[sample.counts.missing.in.one[, .N, by=miR.ID]$miR.ID, on='miR.ID'], aes(x=miR.ID.f, color=miR.missing.or.repeat, y=cpm)) + ggbeeswarm::geom_quasirandom(size=0.5) + facet_wrap(~PNK, scale="free_x", ncol=1) + theme(axis.text.x = element_text(angle=50, hjust=1, size = 6), legend.position = "top", text = element_text(size=6)) + scale_y_log10() + labs(x=NULL, y="Counts Per Median")
save_plot(filename = "./output/figures/miRs_detected_PNK_or_NoPNK_only.pdf", plot = g3, ncol = 1, nrow = 1, base_height = 4, base_width = 5)

g3 <- ggplot(contam.files.dt.miR.1mm.sample.sum[sample.counts.missing.in.one[, .N, by=miR.ID]$miR.ID, on='miR.ID'], aes(x=miR.ID.f, color=miR.missing.or.repeat, y=cpm)) + ggbeeswarm::geom_quasirandom(size=0.5) + facet_wrap(~PNK, scale="free_x", ncol=1) + theme(axis.text.x = element_text(angle=50, hjust=1, size = 6), axis.text.y = element_text(size = 6), legend.position = "none", text = element_text(size=6)) + scale_y_log10(breaks=2^seq(0,10,1)) + labs(x=NULL, y="Counts Per Median")
save_plot(filename = "./output/figures/miRs_detected_PNK_or_NoPNK_only_NoLegend.pdf", plot = g3, ncol = 1, nrow = 1, base_height = 5, base_width = 6)


g4 <- ggplot(contam.files.dt.miR.1mm.sample.sum[sample.counts.missing.in.one[, .N, by=miR.ID]$miR.ID, on='miR.ID'], aes(x=miR.ID.f, color=miR.missing.or.repeat, y=cpm.million+1)) + ggbeeswarm::geom_quasirandom(size=0.5) + facet_wrap(~PNK, scale="free_x", ncol=1) + theme(axis.text.x = element_text(angle=50, hjust=1, size = 6), axis.text.y = element_text(size = 6), legend.position = "none", text = element_text(size=6)) + scale_y_log10(breaks=2^seq(0,10,1)) + labs(x=NULL, y="CPM")
cowplot::save_plot(filename = "./output/figures/miRs_detected_PNK_or_NoPNK_only_NoLegend_permillion.pdf", plot = g4, ncol = 1, nrow = 1, base_height = 5, base_width = 6)
save_plot(filename = "./output/figures/FigEX1C.pdf", plot = g4, ncol = 1, nrow = 1, base_height = 5, base_width = 6)
fwrite(x = fig.ex1c.dat, file="./output/figures/source_data/FigEX1C.csv", sep=",")

# pull full sample size info from the read length files

read.len.fls <- list.files("./data/Plasma/healthyControlsAnalysis/biofragmenta-v1.5_20180220_103636/01-cutadapt", pattern = "*lengthCount.txt", full.names = TRUE)
read.len.dt <- rbindlist(lapply(read.len.fls, 
                                FUN=function(x){
                                  dt <- fread(x, sep="\t")
                                  return(dt)
                                }), fill=TRUE)
setnames(read.len.dt, c("File.Base.ID2", "len", "count"))
read.len.dt[, File.Base.ID2:=sub("_[ATGC]*_l001_l004_l007", "", File.Base.ID2)]
read.len.dt[ALL.SAMPLE.INFO, `:=`(PNK=i.PNK, participant.ID=i.participant.ID, replicate.n=i.replicate.n), on="File.Base.ID2"]
read.len.dt.sum.all <- read.len.dt[, .(count=sum(count)), by=.(participant.ID, PNK)]
read.len.dt.sum.all.f <- read.len.dt[len>=16, .(count=sum(count)), by=.(participant.ID, PNK)]