# response to reviwer 1 ----
fwrite(read.count.seq.by.sample, "./data/Plasma/healthyControlsAnalysis/readCountsSequenceBySamples.txt")
read.count.seq.by.sample <- fread("./data/Plasma/healthyControlsAnalysis/readCountsSequenceBySamples.txt")
setkey(read.count.seq.by.sample, Sample.ID, origSequence)
unique.samples <- read.count.seq.by.sample[, .N, by=Sample.ID]$Sample.ID
dss.list <- DNAStringSetList(lapply(unique.samples, 
                                    function(x){
                                      DNAStringSet(read.count.seq.by.sample[x, on="Sample.ID"]$origSequence)
                                    }))
mono.freq <- sapply(dss.list, simplify = FALSE, USE.NAMES = TRUE, letterFrequency, letters=c("A", "T", "G", "C"))
di.freq <- sapply(dss.list, simplify=FALSE, USE.NAMES=TRUE, dinucleotideFrequency)
freq.dt <- rbindlist(lapply(seq_along(dss.list),
                            FUN=function(i){
                              
                              dss <- dss.list[[i]]
                              mono <- letterFrequency(dss, as.prob=TRUE, letters=c("A", "T", "G", "C"))
                              di <- dinucleotideFrequency(dss, as.prob=TRUE)
                              freqs.dt <- data.table(cbind(mono, di))
                              
                              return(freqs.dt)
                              
                            }))

max.pos <- max(unlist(sapply(dss.list, width)))
pos.freq.dt <- rbindlist(lapply(seq_along(dss.list),
                                FUN=function(i){
                                  dss <- dss.list[[i]]
                                  freqs.dt <- data.table(t(sapply(1:max.pos, nucleotideFrequencyAt, x=dss)))
                                  freqs.dt[, pos:=.I]
                                  freqs.dt.rev <- data.table(t(sapply(1:max.pos, nucleotideFrequencyAt, x=reverse(dss))))
                                  freqs.dt.rev[, pos:=-1*.I]
                                  freqs.dt.out <- rbind(freqs.dt, freqs.dt.rev)
                                  return(freqs.dt.out)
                                  
                                }))
pos.freq.dt[, Sample.ID:=unlist(lapply(unique.samples, rep.int, times=pos.freq.dt[, uniqueN(pos)]))]
pos.freq.dt[, c("participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]


read.count.freqs <- data.table(read.count.seq.by.sample, freq.dt)
read.count.freqs[, seq.len:=nchar(origSequence)]
read.count.freqs[, `:=`(first.nucl=substr(origSequence, 1, 1),
                        first.dinucl=substr(origSequence, 1, 2),
                        last.nucl=substr(origSequence, seq.len, seq.len),
                        last.dinucl=substr(origSequence, seq.len-1, seq.len))]
sd.cols <- colnames(freq.dt)

read.count.freq.summary <- read.count.freqs[, .(colSums(.SD)), by=.(Sample.ID, tot.uniq.count.sample), .SDcols=sd.cols]
read.count.freq.summary$nucl <- sd.cols
read.count.freq.summary[, perc.tot:=V1/tot.uniq.count.sample]
read.count.freq.summary[, c("participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]
ggplot(read.count.freq.summary[nchar(nucl)==1], aes(x=nucl, y=perc.tot, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge")
ggplot(read.count.freq.summary[nchar(nucl)==2], aes(x=nucl, y=perc.tot, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge")


read.count.freq.summary.wt <- read.count.freqs[, .(colSums(.SD*sample.read.count)), by=.(Sample.ID, tot.uniq.count.sample), .SDcols=sd.cols]
read.count.freq.summary.wt$nucl <- sd.cols
read.count.freq.summary.wt[, c("participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]


read.count.freqs.endsummary <- read.count.freqs[, .(count=.N, weight.count=sum(sample.read.count)), by=.(participant.ID, PNK, tot.count.sample, tot.uniq.count.sample, first.nucl, first.dinucl, last.nucl, last.dinucl)]

read.count.freqs.endsummary.tot <- melt.data.table(read.count.freqs.endsummary, id.vars=c("participant.ID", "PNK", "tot.count.sample", "tot.uniq.count.sample", "count", "weight.count"), measure.vars = c("first.nucl", "last.nucl", "first.dinucl", "last.dinucl"), value.name = "nucl")[, .(count=sum(count), weight.count=sum(weight.count)), by=.(participant.ID, PNK, tot.count.sample, tot.uniq.count.sample, variable, nucl)]
read.count.freqs.endsummary.tot.f <- read.count.freqs.endsummary.tot[!grepl("N", nucl)]
read.count.freqs.endsummary.tot.f[, count.perc:=count/tot.uniq.count.sample]
read.count.freqs.endsummary.tot.f[, weight.count.perc:=weight.count/tot.count.sample]
ggplot(read.count.freq.summary[nchar(nucl)==1], aes(x=nucl, y=perc.tot, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge")
ggplot(read.count.freq.summary[nchar(nucl)==2], aes(x=nucl, y=perc.tot, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge")
ggplot(read.count.freqs.endsummary.tot.f, aes(x=nucl, y=weight.count.perc, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge") + facet_wrap(~variable, scale="free_x")
ggplot(read.count.freqs.endsummary.tot.f, aes(x=nucl, y=count.perc, color=PNK, pos=PNK)) + geom_boxplot(pos="dodge") + facet_wrap(~variable, scale="free_x")

pos.freq.dt.melt <- melt.data.table(pos.freq.dt, id.vars=c("pos", 'Sample.ID', "participant.ID", "PNK"), measure.vars = c("A", "C", "G", "T"), value.name="count", variable.name = "nucleotide")
pos.freq.dt.melt[, abs.pos:=abs(pos)]
pos.freq.dt.melt[, which.end:=ifelse(pos>0, "5p", "3p")]
pos.freq.dt.melt[read.count.freq.summary, tot.uniq.count.sample:=i.tot.uniq.count.sample, on="Sample.ID"]
pos.freq.dt.melt[, tot.uniq.count.at.pos:=sum(count), by=.(Sample.ID, pos)]
pos.freq.dt.melt[, perc.tot.at.pos:=count/tot.uniq.count.at.pos]
pos.freq.dt.melt[, perc.tot:=count/tot.uniq.count.sample]
ggplot(pos.freq.dt.melt[abs.pos<10], aes(x=abs.pos, color=PNK, y=perc.tot, group=PNK)) + 
  stat_summary(geom="line", fun.y="mean") + stat_summary(size=0.25, fun.data = "mean_sdl") +
  scale_x_continuous(breaks=seq(1,10,1)) + facet_grid(nucleotide~which.end)

pos.freq.dt.melt[, which.end:=factor(which.end, levels=c("5p", "3p"))]
g <- ggplot(pos.freq.dt.melt[abs.pos<=10], aes(x=abs.pos, color=nucleotide, y=perc.tot, group=nucleotide)) + geom_hline(linetype=2, yintercept = 0.25) +
  stat_summary(geom="line", fun.y="mean") + stat_summary(size=0.25, fun.data = "mean_sdl") +
  scale_x_continuous(breaks=seq(1,10,1)) + scale_y_continuous(limits = c(0.0, 0.5)) + facet_grid(PNK~which.end) + theme(legend.position = "top")
library(cowplot)
g5 <- g %+% pos.freq.dt.melt[abs.pos<=10 & which.end=="5p"]
g3 <- g %+% pos.freq.dt.melt[abs.pos<=10 & which.end=="3p"] + scale_x_reverse(breaks=seq(1,10,1)) 
plot_grid(plotlist = list(g5, g3), ncol=2)


## END response to reviwer 1