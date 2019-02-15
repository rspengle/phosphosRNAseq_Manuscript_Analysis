library(data.table)
library(ggplot2)
library(Biostrings)
library(GenomicAlignments)
library(ggExtra)
library(scales)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

# Load scripts 
source("./scripts/process_synthpool_biofrag_sRNAnalysis.R")

# Input files
data.dir <- "./data/SyntheticPool"
univec.info.fl = dir(data.dir, pattern="UniVec.description", full.names = TRUE, recursive = TRUE)
sample.info.fl=dir(data.dir, pattern="sample_info_synth_use.csv", full.names = TRUE, recursive = TRUE)
biofragmenta.outdir.no4Nrun <- dir(data.dir, pattern="biofragmenta-vsynthetic1.0_20180508_175741", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
pool.fa.fl <- dir(data.dir, pattern="U01_Combined.fa", full.names = TRUE, recursive = TRUE)
subpool.seqs <- readDNAStringSet(pool.fa.fl, use.names=TRUE)
subpool.seqs.dt <- data.table(data.frame(Aln.SeqID=names(subpool.seqs), sequence=as.character(subpool.seqs)))
subpool.seqs.dt[, c("Aln.SeqID.nm", "Aln.Subpool.ID", "AlnSeq.length", "Modification"):=tstrsplit(Aln.SeqID[1], split="|", fixed=TRUE), by=Aln.SeqID]
dt.list.no4N <- process_synthpool_biofrag_sRNAnalysis(biofragmenta.outdir = biofragmenta.outdir.no4Nrun, univec.info = univec.info.fl, sample.info.fl = sample.info.fl)

all.freature.dt.annot.summary <- copy(dt.list.no4N[["annot.summary"]])

# Sanity check to confirm orientation of alignment relative to subpool. ----
# should see mostly Subpool (+) alignments
ggplot(all.freature.dt.annot.summary, aes(x=Sample_Name, y=count, fill=BestMatch.db2)) +
  geom_bar(pos="stack", stat="identity", color="black") + 
  theme_bw() + 
  theme(axis.text.x=element_text(hjust=1, angle=50),
        legend.position = "top",
        legend.direction = "horizontal") + 
  facet_wrap(~library.method, scale="free_x")

# Get process steps ----
all.process.step.summary <- copy(dt.list.no4N[["read.count.summary"]])
all.sample.info <- copy(dt.list.no4N[["sample.info.dt"]])
setnames(all.process.step.summary, "Subpool(+)", "Subpool.Aligned")
all.process.step.summary[, Preprocessed.count.pct.input:=Preprocessed.count/Input.ReadCount]
all.process.step.summary[, Subpool.aligned.pct.Preprocessed:=Subpool.Aligned/Preprocessed.count]
all.process.step.summary[, Unaligned.pct.Preprocessed:=Unaligned.count/Preprocessed.count]
all.process.step.summary[all.sample.info, `:=`(enzyme=i.enzyme, heat=i.HEAT, Treatment=i.Treatment2), on="IndexNum"]

dt.out <- dt.list.no4N[["feature.count.dt"]][BestMatch.db2=="Subpool(+)"]
dt.out[, tot.reads]

# Read full star align instead of sRNAnalysis
all.star.aln.dt.TruSeq <- fread(paste0(biofragmenta.outdir.no4Nrun, "/02b-star_alignment_collapsed/ALL_outAligned.out.txt"))
all.star.aln.dt.TruSeq[, File.Base.ID:=sub("_L00[1-9]_R1_001.*", "", File.Base.ID[1]), by=File.Base.ID]
all.star.aln.dt.TruSeq <- all.star.aln.dt.TruSeq[all.process.step.summary[library.method=="TruSeq"]$File.Base.ID, , on="File.Base.ID"]
all.star.aln.dt.TruSeq[, lib.method:="TruSeq"]
all.star.aln.dt <- copy(all.star.aln.dt.TruSeq)
all.star.aln.dt[, IndexNum:=tstrsplit(File.Base.ID[1], split="_")[2], by=File.Base.ID]
all.star.aln.dt[Pool.SeqID=="UNALIGNED", Pool.ID:="UNALIGNED"]
setkey(all.star.aln.dt, IndexNum)
setkey(all.sample.info, "IndexNum")
all.star.aln.dt.annot <- merge(all.sample.info, all.star.aln.dt, by=intersect(colnames(all.sample.info), colnames(all.star.aln.dt)))
all.star.aln.dt.annot[, Idx.SeqID.count:=paste0(File.Base.ID, "-", Seq.ID, "-", count)]
setnames(all.star.aln.dt.annot, "Pool.ID", "Aln.Subpool.ID")
all.star.aln.dt.annot[, db.matches:=Subpool.ID==Aln.Subpool.ID]
all.star.aln.dt.annot[, cigar.width:=cigarWidthAlongReferenceSpace(CIGAR)]
all.star.aln.dt.annot[, matched.bases:=cigar.width-NM]
setkey(all.star.aln.dt.annot, Idx.SeqID.count)

# Require >= 90% alignment to the aligned sequence. Those that aren't or are too short become "unmapped"
all.star.aln.dt.annot[, min.mm.by.seq:=min(NM), by=.(Idx.SeqID.count)]
all.star.aln.dt.annot[, max.basesmatched.by.seq:=max(matched.bases), by=.(Idx.SeqID.count)]
all.star.aln.dt.annot[, percent.aligned:=matched.bases/Pool.SeqLen]
all.star.aln.dt.annot[, max.percent.aligned.by.seq:=max(percent.aligned), by=.(Idx.SeqID.count)]
all.star.aln.dt.annot.filt <- subset(all.star.aln.dt.annot, max.percent.aligned.by.seq>=0.9 & db.matches==TRUE & matched.bases>=15)

new.unaligned <- unique(all.star.aln.dt.annot[max.percent.aligned.by.seq<0.9 | matched.bases<15 | Aln.Subpool.ID=="UNALIGNED"], by="Idx.SeqID.count")
new.unaligned.f <- new.unaligned[!all.star.aln.dt.annot.filt$Idx.SeqID.count]
new.unaligned.f[, count:=sum(count), by=.(IndexNum)]
new.unaligned.f2 <- subset(new.unaligned.f, Aln.Subpool.ID=="UNALIGNED")

all.star.aln.dt.annot.filt <- rbind(all.star.aln.dt.annot.filt, new.unaligned.f2)


# Get the min mismatches for correct db, correct strand
all.star.aln.dt.annot.filt[, min.mm.correct.db:=min(NM), by=.(Idx.SeqID.count)]
all.star.aln.dt.annot.filt[, max.basesmatched.correct.db:=max(matched.bases), by=.(Idx.SeqID.count)]
all.star.aln.dt.annot.filt[, max.percent.aligned.correct.db:=max(percent.aligned), by=.(Idx.SeqID.count)]

# Now filter to keep only the 
all.star.aln.dt.annot.filt2 <- subset(all.star.aln.dt.annot.filt, (matched.bases==max.basesmatched.correct.db | Aln.Subpool.ID=="UNALIGNED"))
all.star.aln.dt.annot.filt2[, NH:=.N, by=.(Idx.SeqID.count)]
all.star.aln.dt.annot.filt2[, scaled.count:=count/NH]

all.star.aln.dt.annot.filt2[, `:=`(count.tot.feat=sum(count), scaled.count.tot.feat=sum(scaled.count)), by=.(IndexNum, Pool.SeqID)]
setorder(all.star.aln.dt.annot.filt2, -scaled.count)
all.star.aln.dt.annot.filt.feat.tot <- unique(all.star.aln.dt.annot.filt2, by=c("IndexNum", "Pool.SeqID"))


star.aln.tot.summary.aligned <- rbind(all.star.aln.dt.annot.filt2[Pool.SeqID!="UNALIGNED", .(Pool.SeqID=Pool.SeqID[1], count=sum(scaled.count)), by=1:12], all.star.aln.dt.annot.filt2[Pool.SeqID=="UNALIGNED", .(Pool.SeqID=Pool.SeqID[1], count=sum(scaled.count)), by=1:12])
star.aln.tot.summary.aligned[, tot:=sum(count), by=IndexNum]
star.aln.tot.summary.aligned[, perc.tot:=count/tot]
ggplot(star.aln.tot.summary.aligned[Pool.SeqID=="UNALIGNED"], aes(x=Sample_Name, y=perc.tot, fill=enzyme)) + 
  geom_bar(pos="stack", stat="identity") + 
  scale_y_continuous(labels=scales::percent_format()) +
  theme_bw() + 
  labs(y="Subpool Alignment (% Total Trimmed Input)", alpha="Heat Inactivation") +
  theme(axis.text.x=element_text(hjust=1, angle=50),
        legend.position = "top",
        legend.direction = "horizontal") + 
  facet_wrap(~library.method, scale="free_x")

# Star Unaligned %
ggplot(star.aln.tot.summary.aligned[Pool.SeqID!="UNALIGNED"], aes(x=Sample_Name, y=perc.tot, fill=enzyme)) + 
  geom_bar(pos="stack", stat="identity") + 
  scale_y_continuous(labels=scales::percent_format()) +
  theme_bw() + 
  labs(y="Subpool Alignment (% Total Trimmed Input)", alpha="Heat Inactivation") +
  theme(axis.text.x=element_text(hjust=1, angle=50),
        legend.position = "top",
        legend.direction = "horizontal") + 
  facet_wrap(~library.method, scale="free_x")


ggplot(star.aln.tot.summary.aligned[Pool.SeqID!="UNALIGNED"], aes(x=Sample_Name, y=count, fill=enzyme)) + 
  geom_bar(pos="stack", stat="identity") + 
  theme_bw() + 
  labs(y="Subpool Alignment (% Total Trimmed Input)", alpha="Heat Inactivation") +
  theme(axis.text.x=element_text(hjust=1, angle=50),
        legend.position = "top",
        legend.direction = "horizontal") + 
  facet_wrap(~library.method, scale="free_x")

all.star.aln.dt.annot.filt.feat.tot.TruSeq <- all.star.aln.dt.annot.filt.feat.tot[all.sample.info[ library.method=="TruSeq"]$File.Base.ID, on="File.Base.ID"]

all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast <- dcast.data.table(all.star.aln.dt.annot.filt.feat.tot.TruSeq, Pool.SeqID+Pool.SeqLen+Aln.Subpool.ID+Modification~Treatment, value.var="scaled.count", fun.aggregate = sum, fill=0)
setnames(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, c("PNKheat", "none"), c("PNKheat.count", "Untreated.count"))
subpool.seqs.dt.missing <- subpool.seqs.dt[Aln.Subpool.ID=="SynthRNA476v1" & !Aln.SeqID.nm%in%all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast$Pool.SeqID, .(Pool.SeqID=Aln.SeqID.nm, Pool.SeqLen=AlnSeq.length, Aln.Subpool.ID=Aln.Subpool.ID, Modification=Modification, PNKheat.count=0, Untreated.count=0)]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast <- rbind(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, subpool.seqs.dt.missing)
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, `:=`(PNK.heat.tot=sum(PNKheat.count), Untreated.tot=sum(Untreated.count))]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, `:=`(PNK.heat.tot.aligned=sum(ifelse(Aln.Subpool.ID=="UNALIGNED", 0, PNKheat.count)), Untreated.tot.aligned=sum(ifelse(Aln.Subpool.ID=="UNALIGNED", 0, Untreated.count)))]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, `:=`(PNK.cpm.tot=(10^6 * (PNKheat.count+1))/PNK.heat.tot,
                                                       Untreated.cpm.tot=(10^6 * (Untreated.count+1))/Untreated.tot,
                                                       PNK.cpm.totAligned=(10^6 * (PNKheat.count+1))/PNK.heat.tot.aligned,
                                                       Untreated.cpm.totAligned=(10^6 * (Untreated.count+1))/Untreated.tot.aligned)]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, `:=`(log.dif.PNK.Unt=log2(PNK.cpm.tot/Untreated.cpm.tot),
                                                       log.dif.PNK.Unt.aligned=log2(PNK.cpm.totAligned/Untreated.cpm.totAligned))]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, Pool.SeqLen:=as.numeric(Pool.SeqLen)]

recode.cast.vars.unt.aln.cpms <- dcast.data.table(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[Aln.Subpool.ID!="UNALIGNED"], Pool.SeqID+Pool.SeqLen+PNK.cpm.totAligned+Untreated.cpm.totAligned~Modification, fun.aggregate = uniqueN, fill=0, value.var = "Pool.SeqID")
recode.cast.vars.unt.aln.cpms[, mean.logcpm:=rowMeans(log2(.SD)), .SDcols=c("Untreated.cpm.totAligned", "PNK.cpm.totAligned")]
recode.cast.vars.unt.aln.cpms[, logDiff.PNK.minus.Untr:=log2(PNK.cpm.totAligned/Untreated.cpm.totAligned)]
recode.cast.vars.unt.aln.cpms[, c("PNK.cpm.totAligned", "Untreated.cpm.totAligned"):=NULL]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, size.group:=factor(ifelse(Pool.SeqLen%in%seq(15, 90, 15), Pool.SeqLen, "(17-25)"), levels=c("15", "(17-25)", as.character(seq(30, 90, 15))))]
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[, Modification:=factor(Modification, levels=c("5p-phosphorylation", "None", "3p-phosphorylation", "5p-phosphorylation_AND_2pO-methylation"))]

g <- ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[!is.na(Modification) & Modification%in%c("5p-phosphorylation", "None")], aes(y=PNK.cpm.totAligned, x=Untreated.cpm.totAligned, color=Modification)) + geom_point(alpha=0.7, size=0.25) + scale_x_log10(limits=c(10^-1,10^5), breaks=10^seq(-1,5), labels = trans_format("log10", math_format())) + scale_y_log10(breaks=10^seq(-1,5), limits=c(10^-1,10^5), labels = trans_format("log10", math_format())) + labs(x="CPM Untreated", y="CPM PNK-treated") + theme_bw() + theme(legend.position = "bottom", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black", size=6)); g
fig.fl="./output/figures/main/FIG1/FIG1C_PNK_VS_Untreated_Synth_NODensity.pdf"
ggsave(plot=g, filename = basename(fig.fl), width = 2.5, height = 2.5)

fig.fl="output/figures/main/FIG1/FIG1C_PNK_VS_Untreated_Synth_WithDensity.pdf"
ggsave(plot=ggExtra::ggMarginal(g, groupFill = TRUE, type="density"), filename = fig.fl, width = 2.5, height = 2.5)
#ggsave(plot=ggExtra::ggMarginal(g, groupFill = TRUE, type="violin"), filename = basename(fig.fl), width = 2.5, height = 2.5)
#ggExtra::ggMarginal(g, groupFill = TRUE, type="violin", draw_quantiles=c(0.25, 0.5, 0.75), pos="dodge")

# 2019-01-27 Add a similar plot for 3' phosph ----
g <- ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast[!is.na(Modification) & Modification%in%c("5p-phosphorylation", "3p-phosphorylation")], aes(y=PNK.cpm.totAligned, x=Untreated.cpm.totAligned, color=Modification)) + geom_point(alpha=0.7, size=0.25) + scale_x_log10(limits=c(10^-1,10^5), breaks=10^seq(-1,5), labels = trans_format("log10", math_format())) + scale_y_log10(breaks=10^seq(-1,5), limits=c(10^-1,10^5), labels = trans_format("log10", math_format())) + labs(x="CPM Untreated", y="CPM PNK-treated") + theme_bw() + theme(legend.position = "bottom", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black", size=6)); g
fig.fl="./output/figures/main/FIG1/FIG1D_PNK_VS_Untreated_Synth_3pPhosph_NODensity.pdf"
ggsave(plot=g, filename = basename(fig.fl), width = 2.5, height = 2.5)

fig.fl="output/figures/main/FIG1/FIG1D_PNK_VS_Untreated_Synth_3pPhosph_WithDensity.pdf"
ggsave(plot=ggExtra::ggMarginal(g, groupFill = TRUE, type="density"), filename =  basename(fig.fl), width = 2.5, height = 2.5)
#ggsave(plot=ggExtra::ggMarginal(g, groupFill = TRUE, type="violin"), filename = basename(fig.fl), width = 2.5, height = 2.5)
#ggExtra::ggMarginal(g, groupFill = TRUE, type="violin", draw_quantiles=c(0.25, 0.5, 0.75), pos="dodge")


all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt <- melt(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, id.vars = c("Pool.SeqID", "Modification"), measure.vars = c("PNK.cpm.totAligned", "Untreated.cpm.totAligned"), variable.name = "Treatment.meas", value.name = "CPM")
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt[, Treatment:=ifelse(Treatment.meas=="PNK.cpm.totAligned", "PNK", ifelse(Treatment.meas=="Untreated.cpm.totAligned", "Untreated", NA))]
g <- ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt[!is.na(Modification) & Modification%in%c("5p-phosphorylation", "3p-phosphorylation", "None")], aes(y=CPM, x=Modification, fill=Treatment)) + 
  geom_violin(draw_quantiles = seq(0.25, 0.75, 0.25), size=0.25) + scale_y_log10(limits=c(10^-1,10^5), breaks=10^seq(-1,5), labels = trans_format("log10", math_format())) + theme_bw() + theme(legend.position = "top", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black", size=6)); g
fig.fl="output/figures/main/FIG1/FIG1C_PNK_VS_Untreated_VIOLIN.pdf"
ggsave(plot=g, filename = basename(fig.fl), width = 3, height = 2.5)
fig.fl="output/figures/main/FIG1/FIG1C_PNK_VS_Untreated_BOXPLOT.pdf"
g <- ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt[!is.na(Modification) & Modification%in%c("5p-phosphorylation", "3p-phosphorylation", "None")], aes(y=CPM, x=Modification, pos=Treatment)) + 
  geom_boxplot(size=0.25, width=0.5, aes(fill=Treatment), outlier.colour = NA) + 
  geom_quasirandom(dodge.width=0.5, size=0.2, alpha=0.5, fill="black", color="black") + 
  scale_y_log10(limits=c(10^-1,10^5), breaks=10^seq(-1,5), labels = trans_format("log10", math_format())) + theme_bw() + theme(legend.position = "top", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black", size=6)); g
ggsave(plot=g, filename = basename(fig.fl), width = 3, height = 2.5)


compare_means(CPM~Treatment, data = all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt, group="Modification", ref.group = "Untreated", paired=TRUE, alternative="greater")
all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.remelt[, median(CPM), by=.(Treatment, Modification)]

# Summarize libs for pub ----
ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, aes(x=log2(PNKheat.count+1), color=Modification)) + stat_ecdf()
ggplot(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, aes(x=log2(Untreated.count+1), color=Modification)) + stat_ecdf()


fwrite(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, file = "20181125_PNK_VS_UNTREATED_TruSeq_SynthPool_counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.geo <- subset(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast, select=Hmisc::Cs(Pool.SeqID, Pool.SeqLen, Modification, PNKheat.count, Untreated.count, PNK.cpm.totAligned, Untreated.cpm.totAligned))

fwrite(all.star.aln.dt.annot.filt.feat.tot.TruSeq.cast.geo, file = "../2019_EMBO_GEO_SUBMISSION/SYNTHETIC/PNK_VS_UNTREATED_TruSeq_SynthPool_counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
