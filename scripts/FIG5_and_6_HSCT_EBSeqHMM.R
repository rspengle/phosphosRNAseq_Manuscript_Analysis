library(data.table)
library(EBSeqHMM)
library(ggplot2)
library(scales)
library(Hmisc)
library(ggbeeswarm)
library(DESeq2)
library(clue)
library(e1071)
library(WGCNA)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(GGally)
library(clustree)
options(stringsAsFactors = FALSE);
#enableWGCNAThreads(nThreads = 8)

setwd( "/data/ANALYSIS/BIOFRAGMENTA/2018_BIOFRAGMENTA_COMBINED_ANALYSIS")

ebseq.file.p07iter="P07_EBSeqHMMOut_iter.RData"
ebseq.file.p04iter="P04_EBSeqHMMOut_iter.RData"
load(file = ebseq.file.p07iter)
load(file = ebseq.file.p04iter)
load("p04_rdata_for_ebseqhmm.RData")
load("p07_rdata_for_ebseqhmm.RData")
load("BMT_lab_data_table_orignal.RData")
qnames.merged <- fread("gene_annotation_information_bmt.txt")
load(file = system.file("extdata", "combine-expression.rda", package = "TissueEnrich"))
tissue.specific.data <- rbind(data.table(source.data="ProteinAtlas", dataset$`Protein-Atlas`$tissueSpecificGenes),
                              data.table(source.data="GTEx", dataset$`GTEx-Combine`$tissueSpecificGenes))
setnames(tissue.specific.data, "Gene", "ensgid")
# read median cpm values from GtEx
gtex.mediantpm <- fread("/home/ryanspen/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", skip=2)
gtex.mediantpm.matr <- data.frame(gtex.mediantpm, row.names=1)[, 2:(ncol(gtex.mediantpm.matr)-1)]
row.names(gtex.mediantpm.matr) <- sub("\\.[0-9]*$", "", row.names(gtex.mediantpm.matr))



teEnrichment.multi <- function(query.list, background.c, tissue.groups.use=c("Tissue-Enriched", "Tissue-Enhanced", "Group-Enriched"), source.use=c("GTEx", "ProteinAtlas")){
  e <- new.env()
  load(file = system.file("extdata", "combine-expression.rda", package = "TissueEnrich"), envir = e)
  tissue.specific.data <- rbind(data.table(source.data="ProteinAtlas", e$dataset$`Protein-Atlas`$tissueSpecificGenes),
                                data.table(source.data="GTEx", e$dataset$`GTEx-Combine`$tissueSpecificGenes))
  setkey(tissue.specific.data, "Gene")
  tissue.specific.data.use <- copy(tissue.specific.data[tissue.groups.use, on="Group"])
  nTeGenesInTissue.dt <- tissue.specific.data.use[, .N, by=.(Gene, Tissue, source.data), on="Gene", nomatch=0][, .(nTeGenesInTissue=0), by=.(Tissue, source.data)]
  query.dt <- rbindlist(lapply(seq_along(query.list), 
                               function(x){
                                 q.x <- query.list[[x]]
                                 name.x <- names(query.list)[x]
                                 name.x <- ifelse(is.null(name.x), paste0("Query.", x), name.x)
                                 if(is.vector(q.x)){
                                   dtx <- data.table(NAME=q.x, query.set=name.x)
                                   dtx[, Gene:=tstrsplit(NAME, split=".", fixed=TRUE, keep=1)]
                                 } else{
                                   stop("Query must be a list of ensembl gene ids")
                                 }
                                 return(dtx)
                               }))
  
  
  query.dt.f <- query.dt[tissue.specific.data.use[, .N, by=Gene]$Gene, on="Gene", nomatch=0]
  query.dt.f[, nTotalInputGenes:=.N, by=.(query.set)]
  nTeGenesInTissue.dt.tmp <- tissue.specific.data.use[background.c, .N, by=.(Gene, Tissue, source.data), on="Gene", nomatch=0][, .(nTeGenesInTissue=.N), by=.(Tissue, source.data)]
  nTeGenesInTissue.dt[nTeGenesInTissue.dt.tmp, nTeGenesInTissue:=nTeGenesInTissue+i.nTeGenesInTissue, on=c("Tissue", "source.data")]
  query.dt.f.dt.tissue.summary <- rbindlist(lapply(names(query.list), 
                                                   function(x){
                                                     return(data.table(nTeGenesInTissue.dt, query.set=x, nTotalInputGenes=0L, overlapGenes=0L))
                                                   }))
  
  nTotalGenes.dt <- tissue.specific.data[background.c, .N, on="Gene", by=.(Gene, source.data), nomatch=0][, .(nTotalGenes=.N), by=source.data]
  
  query.dt.f.dt.tissue <- merge(query.dt.f, tissue.specific.data.use, allow.cartesian=TRUE, on="Gene", all.x=TRUE)
  
  query.dt.f.dt.tissue[query.dt.f.dt.tissue.summary, nTeGenesInTissue:=i.nTeGenesInTissue, on=c("Tissue", "source.data", "query.set")]
  
  query.dt.f.dt.tissue.summary.tmp <- query.dt.f.dt.tissue[!is.na(source.data), .(nTotalInputGenes=nTotalInputGenes[1],
                                                                                  overlapGenes=.N), by=.(query.set, source.data, Tissue)]
  
  query.dt.f.dt.tissue.summary[query.dt.f.dt.tissue.summary.tmp, `:=`(overlapGenes=i.overlapGenes), on=c("Tissue", "source.data", "query.set")]
  query.dt.f.dt.tissue.summary[query.dt.f.dt.tissue.summary.tmp, nTotalInputGenes:=i.nTotalInputGenes, on=c("source.data", "query.set")]
  query.dt.f.dt.tissue.summary[, nTotalInputGenes:=max(nTotalInputGenes), by=c("query.set")]
  
  query.dt.f.dt.tissue.summary[nTotalGenes.dt, nTotalGenes:=i.nTotalGenes, on="source.data"]
  query.dt.f.dt.tissue.summary[, pValue:=stats::phyper(overlapGenes-1, nTeGenesInTissue, 
                                                       nTotalGenes - nTeGenesInTissue, nTotalInputGenes, 
                                                       lower.tail = FALSE)]
  query.dt.f.dt.tissue.summary[, `:=`(fold.change=(overlapGenes/nTotalInputGenes)/(nTeGenesInTissue/nTotalGenes),
                                      Log10PValue=(-1*log10(pValue))
  )]
  
  query.dt.f.dt.tissue.summary[, query.source:=paste0(query.set, ".", source.data)]
  pvals <- split(query.dt.f.dt.tissue.summary$pValue, f=query.dt.f.dt.tissue.summary$query.source)
  p.adjs <- lapply(pvals, p.adjust, method="BH")
  query.dt.f.dt.tissue.summary[, pAdj:=p.adjust(pValue, method="BH"), by=.(query.set, source.data)]
  query.dt.f.dt.tissue.summary[, Log10PAdj:=(-1*log10(pAdj))]
  
  
  
  setnames(query.dt.f.dt.tissue.summary, c("overlapGenes", "source.data"), c("Tissue.Specific.Genes", "Database"))
  return(query.dt.f.dt.tissue.summary[])
}


wgcna.clust.f <- function(m){
  counts.matr <- t(m)
  powers = c(c(1:10), seq(from =10, to=30, by=1))
  sft=pickSoftThreshold(counts.matr, dataIsExpr = TRUE,powerVector = powers,verbose =5,networkType = "signed")
  softPower <- sft$powerEstimate
  
  adjacency = adjacency(counts.matr, power = softPower, type = "signed") #specify network type
  
  
  # Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
  #translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
  TOM = TOMsimilarity(adjacency,TOMDenom="mean", TOMType="signed") # specify network type
  dissTOM = 1-TOM
  
  # Generate Modules --------------------------------------------------------
  
  library(flashClust)
  # Generate a clustered gene tree
  geneTree = flashClust(as.dist(dissTOM), method="average")
  
  #This sets the minimum number of genes to cluster into a module
  minModuleSize = 10
  dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2,respectSmallClusters=TRUE, verbose=4, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
  dynamicColors= labels2colors(dynamicMods)
  MEList= moduleEigengenes(counts.matr, colors= dynamicColors,softPower = softPower)
  MEs= MEList$eigengenes
  MEDiss= 1-cor(MEs)
  METree= flashClust(as.dist(MEDiss), method= "average")
  
  
  
  merge = mergeCloseModules(counts.matr, dynamicColors,  verbose =4)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  module.dt <- data.frame(ensgid=colnames(counts.matr), module.group=merge$colors)
  return(list("modules"=module.dt, "WGCNA"=merge, "disTOM"=as.dist(dissTOM)))
}



get.transition.probs.f <- function(EBSeqHMMGeneOut, path.ids=NULL){
  SigPPLarge <- EBSeqHMMGeneOut$MgAllPP
  SigMAPLargeChar <- EBSeqHMMGeneOut$MgAllMAPChar
  SigMaxValLarge <- EBSeqHMMGeneOut$MgAllMaxVal
  AllPaths <- colnames(EBSeqHMMGeneOut$MgAllPP)
  WithUp <- grep("Up", AllPaths)
  WithDown <- grep("Down", AllPaths)
  UpAndDown <- union(WithUp, WithDown)
  AllEEPath <- AllPaths[-UpAndDown]
  SigPPLarge[, AllEEPath]
  
  cond.split <- limma::strsplit2(colnames(SigPPLarge), split="-", fixed=TRUE)
  SigPPLarge.f <- t(SigPPLarge[!is.na(SigPPLarge[,1]),])
  cond.sum.ee <- sapply(seq(ncol(cond.split)),
                        FUN=function(x){
                          x.cols <- cond.split[,x]=="EE"
                          up.cols <- cond.split[,x]=="Up"
                          dn.cols <- cond.split[,x]=="Down"
                          pEE <- colSums(SigPPLarge.f[x.cols,])
                          pUP <- colSums(SigPPLarge.f[up.cols,])
                          pDN <- colSums(SigPPLarge.f[dn.cols,])
                          dir.out <- ifelse(pUP>pDN, 1, -1)
                          dir.out[pEE>0.05] <- 0
                          return(dir.out)
                        })
  
  cond.sum.dt <- rbindlist(lapply(seq(ncol(cond.split)),
                                  FUN=function(x){
                                    x.cols <- cond.split[,x]=="EE"
                                    up.cols <- cond.split[,x]=="Up"
                                    dn.cols <- cond.split[,x]=="Down"
                                    pEE <- colSums(SigPPLarge.f[x.cols,])
                                    pUP <- colSums(SigPPLarge.f[up.cols,])
                                    pDN <- colSums(SigPPLarge.f[dn.cols,])
                                    rbind(data.table(direction="EE", p=pEE, ol.geneid=colnames(SigPPLarge.f), x=x),
                                          data.table(direction="Up", p=pUP, ol.geneid=colnames(SigPPLarge.f), x=x),
                                          data.table(direction="Down", p=pDN, ol.geneid=colnames(SigPPLarge.f), x=x))
                                  }))
  cond.sum.dt.cast <- dcast.data.table(cond.sum.dt, ol.geneid+x~direction, fun.aggregate = min, value.var="p")
    cond.sum.dt.cast[, dir.majority.call:=ifelse(Down>0.5, "Down", ifelse(Up>0.5, "Up", "EE"))]
  cond.sum.dt.cast[, dir.sig.call:=ifelse(EE>0.1, "EE", ifelse(Down>Up, "Down", ifelse(Down<Up, "Up", NA)))]
  cond.sum.dt.cast[, dir.sig.call2:=ifelse((1-EE)<0.05, "EE", ifelse((1-Down)<0.05, "Down", ifelse((1-Up)<0.05, "Up", "EE")))]
  cond.sum.dt.paths <- cond.sum.dt.cast[, .(sig.path=paste0(dir.sig.call, collapse="-"),
                                            sig2.path=paste0(dir.sig.call2, collapse="-"),
                                            maj.path=paste0(dir.majority.call, collapse="-")), by=ol.geneid]
  
  
    cond.sum.dt.cast[, ensgid:=sub("\\.[0-9]*$", "", ol.geneid)]
  if(is.null(paths07)){
    cond.sum.dt.cast[, transition:=paste0("t", x-1, "-t", x, "_", dir.sig.call)]
    cond.sum.dt.cast[, dir.pos:=paste0(transition, "_", dir.sig.call)]
    
  } else{
    cond.sum.dt.cast[, transition:=path.ids[x]]
    cond.sum.dt.cast[, dir.pos:=paste0(transition, "_", dir.sig.call)]
    
  }
  return(cond.sum.dt.cast[])
}


GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {		
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))		
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}

p07.expr.dt <- data.table(melt(counts(p07.dds, normalized=TRUE)))
p07.expr.dt[data.table(melt(log2(fpm(p07.dds, robust=TRUE)+1))), gmpr.log2.cpm:=i.value, on=c("Var1", "Var2")]
p04.expr.dt <- data.table(melt(counts(p04.dds, normalized=TRUE)))
p04.expr.dt[data.table(melt(log2(fpm(p04.dds, robust=TRUE)+1))), gmpr.log2.cpm:=i.value, on=c("Var1", "Var2")]

p04.counts.raw <- counts(p04.dds, normalized=FALSE)
p07.counts.raw <- counts(p07.dds, normalized=FALSE)


p04.sizes <- sizeFactors(p04.dds) * exp(mean(log(colSums(counts(p04.dds, normalized=FALSE)))))
p07.sizes <- sizeFactors(p07.dds) * exp(mean(log(colSums(counts(p07.dds, normalized=FALSE)))))
p04.counts.per.median <- log2(sweep(p04.counts.raw+1, 2, p04.sizes, "/") * mean(p04.sizes))
p04.counts.per.median.z <- t(apply(p04.counts.per.median, 1, function(x) (x-mean(x))/sd(x)))
p07.counts.per.median <- log2(sweep(p07.counts.raw+1, 2, p07.sizes, "/") * mean(p07.sizes))
p07.counts.per.median.z <- t(apply(p07.counts.per.median, 1, function(x) (x-mean(x))/sd(x)))


p04.counts <- counts(p04.dds, normalized=TRUE) * mean(colSums(counts(p04.dds, normalized=FALSE)))
row.names(p04.counts) <- sub("\\.[0-9]*$", "", row.names(p04.counts))
p07.counts <- counts(p07.dds, normalized=TRUE)
row.names(p07.counts) <- sub("\\.[0-9]*$", "", row.names(p07.counts))


#mbseq.dat07 <- MBCluster.Seq::RNASeq.Data(p07.counts.raw, Normalizer = log(as.numeric(sizeFactors(p07.dds))), Treatment=seq_along(colnames(p07.counts)), GeneID=row.names(p07.counts))
#mbseq.dat04 <- MBCluster.Seq::RNASeq.Data(p04.counts.raw, Normalizer = log(as.numeric(sizeFactors(p04.dds))), Treatment=seq_along(colnames(p04.counts)), GeneID=row.names(p04.counts))


#c0=KmeansPlus.RNASeq(mbseq.dat07, nK=100, print.steps = TRUE)
#cls=Cluster.RNASeq(data=mbseq.dat07,model="nbinom",centers=c0$centers,method="DA",iter.max = 1000)
#tr=Hybrid.Tree(data=mbseq.dat07,cluste=cls$cluster,model="nbinom")
#plotHybrid.Tree(merge=tr,cluster=cls$cluster,logFC=mbseq.dat07$logFC,tree.title=NULL)

#clsSA=Cluster.RNASeq(data=mbseq.dat07,model="nbinom",centers=c0$centers,method="SA",iter.max = 1000)
#trSA=Hybrid.Tree(data=mbseq.dat07,cluste=clsSA$cluster,model="nbinom")
#plotHybrid.Tree(merge=trSA,cluster=clsSA$cluster,logFC=mbseq.dat07$logFC,tree.title=NULL)

p04.rlog <- rlog(p04.dds, blind=TRUE)
p07.rlog <- rlog(p07.dds, blind=TRUE)
p07.rlog.z <- t(apply(assay(p07.rlog), 1, function(x) (x-mean(x))/sd(x)))
p04.rlog.z <- t(apply(assay(p04.rlog), 1, function(x) (x-mean(x))/sd(x)))
p04.rlog.nb <- rlog(p04.dds, blind=FALSE)
p07.rlog.nb <- rlog(p07.dds, blind=FALSE)
p07.rlog.nb.z <- t(apply(assay(p07.rlog.nb), 1, function(x) (x-mean(x))/sd(x)))
p04.rlog.nb.z <- t(apply(assay(p04.rlog.nb), 1, function(x) (x-mean(x))/sd(x)))

p07.expr.dt[data.table(melt(assay(p07.rlog))), gmpr.rlog.blind:=i.value, on=c("Var1", "Var2")]
p04.expr.dt[data.table(melt(assay(p04.rlog))), gmpr.rlog.blind:=i.value, on=c("Var1", "Var2")]
p07.expr.dt[data.table(melt(assay(p07.rlog.nb))), gmpr.rlog.notblind:=i.value, on=c("Var1", "Var2")]
p04.expr.dt[data.table(melt(assay(p04.rlog.nb))), gmpr.rlog.notblind:=i.value, on=c("Var1", "Var2")]


setnames(p07.expr.dt, c("Var1", "Var2", "value"), c("ol.geneid", "File.Base.ID2", "gmpr.counts"))
setnames(p04.expr.dt, c("Var1", "Var2", "value"), c("ol.geneid", "File.Base.ID2", "gmpr.counts"))
p07.expr.dt[, `:=`(participant.ID="P07",
                   Day.num=colData(p07.dds)[File.Base.ID2[1],]$Day.num,
                   BMT.Day=colData(p07.dds)[File.Base.ID2[1],]$BMT.Day), by=File.Base.ID2]

p04.expr.dt[, `:=`(participant.ID="P04",
                   Day.num=colData(p04.dds)[File.Base.ID2[1],]$Day.num,
                   BMT.Day=colData(p04.dds)[File.Base.ID2[1],]$BMT.Day), by=File.Base.ID2]

expr.dt <- rbind(p07.expr.dt, p04.expr.dt)
expr.dt[, gmpr.log2.count:=log2(gmpr.counts+1)]
expr.dt[, `:=`(gmpr.log2.cpm.z=(gmpr.log2.count-mean(gmpr.log2.count))/sd(gmpr.log2.count),
               gmpr.log2.count.z=(gmpr.log2.cpm-mean(gmpr.log2.cpm))/sd(gmpr.log2.cpm),
               gmpr.rlog.blind.z=(gmpr.rlog.blind-mean(gmpr.rlog.blind))/sd(gmpr.rlog.blind),
               gmpr.rlog.notblind.z=(gmpr.rlog.notblind-mean(gmpr.rlog.notblind))/sd(gmpr.rlog.notblind)), by=.(participant.ID, ol.geneid)]
expr.dt[, ensgid:=sub("\\.[0-9]*$", "", ol.geneid[1]), by=ol.geneid]

expr.dt.tissue <- merge(expr.dt, tissue.specific.data, by="ensgid", all.x=TRUE, allow.cartesian=TRUE)

GeneDECalls07 <- GetDECalls(EBSeqHMMGeneOut07.iter, FDR=.01)
row.names(GeneDECalls07) <- sub("\\.[0-9]*$", "", row.names(GeneDECalls07))
GeneDECalls04 <- GetDECalls(EBSeqHMMGeneOut04.iter, FDR=.01)
row.names(GeneDECalls04) <- sub("\\.[0-9]*$", "", row.names(GeneDECalls04))


allpp07 <- EBSeqHMMGeneOut07.iter$MgAllPP[!is.na(EBSeqHMMGeneOut07.iter$MgAllPP[,1]),]
all.cond <- limma::strsplit2(colnames(allpp07), split="-")

paths07 <- sub("-Day.", "-", paste0(time.facts[-length(time.facts)], "-", time.facts[-1]))
paths04 <- sub("-Day.", "-", paste0(time.facts04[-length(time.facts04)], "-", time.facts04[-1]))



p07.paths.probs.dt <- get.transition.probs.f(EBSeqHMMGeneOut07.iter, path.ids=paths07)
p04.paths.probs.dt <- get.transition.probs.f(EBSeqHMMGeneOut04.iter, path.ids=paths04)
p04.paths.probs.dt[, max.val:=apply(.SD, 1, max), .SDcols=c("Down", "EE", "Up")]
p07.paths.probs.dt[, max.val:=apply(.SD, 1, max), .SDcols=c("Down", "EE", "Up")]

p07.counts.sig <- p07.counts[row.names(p07.counts)%in%row.names(GeneDECalls07),]
p04.counts.sig <- p04.counts[row.names(p04.counts)%in%row.names(GeneDECalls04),]

p07.counts.sig.z <- t(apply(p07.counts.sig, 1, FUN=function(x) (x-mean(x))/sd(x)))


p07.modules.all <- wgcna.clust.f(p07.counts.sig)
p04.modules.all <- wgcna.clust.f(p04.counts.sig)

final.07 <- p07.modules.all$modules$module.group
names(final.07)  <- p07.modules.all$modules$ensgid
as.data.frame(p07.counts.sig.z) %>% mutate(Cluster = final.07) %>% group_by(Cluster) %>% summarise_all("mean") %>% kable() %>% kable_styling()


p07.MEs.dt <- data.table(melt(p07.modules.all$WGCNA$newMEs))
setnames(p07.MEs.dt, c("module.group", "eigengene.value"))
p07.MEs.dt[, module.group:=sub("ME", "", module.group)]
p04.MEs.dt <- data.table(melt(p04.modules.all$WGCNA$newMEs))
setnames(p04.MEs.dt, c("module.group", "eigengene.value"))
p04.MEs.dt[, module.group:=sub("ME", "", module.group)]
p04.MEs.dt[, File.Base.ID2:=colnames(p04.counts)]
p04.MEs.dt[, Day.num:=p04.dds$Day.num]
p07.MEs.dt[, File.Base.ID2:=colnames(p07.counts)]
p07.MEs.dt[, Day.num:=p07.dds$Day.num]


p07.modules <- p07.modules.all[[1]]
p04.modules <- p04.modules.all[[1]]
all.modules <- rbind(data.table(participant.ID="P07", p07.modules),data.table(participant.ID="P04", p04.modules))
all.modules[, participant.module:=paste0(participant.ID, ".", module.group)]
all.modules.list07 <- split(p07.modules$ensgid, f=p07.modules$module.group)
all.modules.list04 <- split(p04.modules$ensgid, f=p04.modules$module.group)
all.modules.teSig07 <- teEnrichment.multi(all.modules.list07, row.names(p07.counts))
all.modules.teSig04 <- teEnrichment.multi(all.modules.list04, row.names(p04.counts))
all.modules.teSig07[, participant.ID:="P07"]
all.modules.teSig04[, participant.ID:="P04"]
all.modules.te <- rbind(all.modules.teSig04, all.modules.teSig07)


p07.paths.probs.dt[p07.modules, module.group:=i.module.group, on=c("ensgid")]
p04.paths.probs.dt[p04.modules, module.group:=i.module.group, on=c("ensgid")]
p04.paths.probs.dt[, transition:=factor(transition, levels=p04.paths.probs.dt[, .N, by=.(transition, x)][order(x)]$transition)]
p07.paths.probs.dt[, transition:=factor(transition, levels=p07.paths.probs.dt[, .N, by=.(transition, x)][order(x)]$transition)]
p07.paths.probs.dt[, module.group.size:=uniqueN(ol.geneid), by=module.group]
p04.paths.probs.dt[, module.group.size:=uniqueN(ol.geneid), by=module.group]
p04.paths.probs.dt.module.sum <- p04.paths.probs.dt[, .(pp.gene.trans=colSums(.SD)), .SDcols=c("Down", "EE", "Up"), by=.(transition, ensgid, x, module.group, module.group.size)]
p04.paths.probs.dt.module.sum[, expr.direction:=c("Down", "EE", "Up")]
p07.paths.probs.dt.module.sum <- p07.paths.probs.dt[, .(pp.gene.trans=colSums(.SD)), .SDcols=c("Down", "EE", "Up"), by=.(transition, ensgid, x,  module.group, module.group.size)]
p07.paths.probs.dt.module.sum[, expr.direction:=c("Down", "EE", "Up")]

ggplot(p07.paths.probs.dt.module.sum[expr.direction!="EE"], aes(x=transition, fill=expr.direction, color=expr.direction, group=expr.direction, y=pp.gene.trans)) + 
  stat_summary(fun.data = "mean_sdl", fun.args=list(mult=1)) +
  stat_summary(fun.y = "mean", geom="line") + facet_wrap(~module.group, ncol=2)


p07.paths.probs.dt.tissue <- merge(p07.paths.probs.dt, tissue.specific.data, by="ensgid", all.x=TRUE)
p07.paths.probs.dt.tissue[, n.genes.transition:=uniqueN(ensgid), by=.(transition, Tissue)]
p04.paths.probs.dt.tissue <- merge(p04.paths.probs.dt, tissue.specific.data, by="ensgid", all.x=TRUE)
p04.paths.probs.dt.tissue[, n.genes.transition:=uniqueN(ensgid), by=.(transition, Tissue)]
p07.paths.probs.dt.tissue.sum <- p07.paths.probs.dt.tissue[, .(n.genes.call=uniqueN(ensgid)), by=.(dir.sig.call, x, transition, dir.pos, Tissue, n.genes.transition,tree.module)]
p07.paths.probs.dt.tissue.sum[, perc.genes.call:=n.genes.call/n.genes.transition]
p04.paths.probs.dt.tissue.sum <- p04.paths.probs.dt.tissue[, .(n.genes.call=uniqueN(ensgid)), by=.(dir.sig.call, x, transition, dir.pos, Tissue, n.genes.transition)]
p04.paths.probs.dt.tissue.sum[, perc.genes.call:=n.genes.call/n.genes.transition]




ggplot(p04.paths.probs.dt.tissue.sum[!is.na(Tissue) & dir.sig.call!="EE"], aes(x=x, fill=dir.sig.call, y=perc.genes.call)) + geom_bar(pos="stack", stat="identity") + facet_wrap(~Tissue)
p04.paths.probs.dt.tissue.cast <- data.frame(dcast.data.table(p04.paths.probs.dt.tissue, ensgid~transition, fun.aggregate = function(x){ifelse(x=="EE", 0, ifelse(x=="Up", 1, -1))}, value.var="dir.sig.call", fill=0), row.names=1)
p07.paths.probs.dt.tissue.cast <- data.frame(dcast.data.table(p07.paths.probs.dt.tissue, ensgid~transition, fun.aggregate = function(x){ifelse(x=="EE", 0, ifelse(x=="Up", 1, -1))}, value.var="dir.sig.call", fill=0), row.names=1)
p07.gene.annot.liver.bm <- data.frame(dcast.data.table(p07.paths.probs.dt.tissue[c("Liver", "Bone.Marrow"), .N,on="Tissue", by=.(ensgid, Tissue)], ensgid~Tissue, fun.aggregate = length), row.names=1)
p04.gene.annot.liver.bm <- data.frame(dcast.data.table(p04.paths.probs.dt.tissue[c("Liver", "Bone.Marrow"), .N,on="Tissue", by=.(ensgid, Tissue)], ensgid~Tissue, fun.aggregate = length), row.names=1)
p04.paths.probs.dt.tissue.cast.f <- p04.paths.probs.dt.tissue.cast[row.names(p04.paths.probs.dt.tissue.cast)%in%row.names(p04.gene.annot.liver.bm),p04.paths.probs.dt.tissue[, .N, by=.(x, transition)][order(x), make.names(transition)]]
p07.paths.probs.dt.tissue.cast.f <- p07.paths.probs.dt.tissue.cast[row.names(p07.paths.probs.dt.tissue.cast)%in%row.names(p07.gene.annot.liver.bm),p07.paths.probs.dt.tissue[, .N, by=.(x, transition)][order(x), make.names(transition)]]
p07.gene.annot.liver.bm <- p07.gene.annot.liver.bm[row.names(p07.paths.probs.dt.tissue.cast.f),]
p07.gene.annot.liver.bm$module.group <- data.frame(p07.modules, row.names=1)[row.names(p07.gene.annot.liver.bm), ]
p04.gene.annot.liver.bm <- p04.gene.annot.liver.bm[row.names(p04.paths.probs.dt.tissue.cast.f),]
p04.gene.annot.liver.bm$module.group <- data.frame(p04.modules, row.names=1)[row.names(p04.gene.annot.liver.bm), ]



gtex.mediantpm.matr.melt <- data.table(melt(gtex.mediantpm, id=c("gene_id", "Description")))
gtex.mediantpm.matr.melt[, tissue.simple:=tstrsplit(variable[1], split=" - ", keep=1), by=variable]
gtex.mediantpm.matr.melt.med <- gtex.mediantpm.matr.melt[, .(median.tmp=median(value)), by=.(gene_id, Description, tissue.simple)]
gtex.mediantpm.matr.melt.med[, median.tmp.z:=(median.tmp-mean(median.tmp))/sd(median.tmp), by=gene_id]
gtex.mediantpm.matr.melt.med[, ensgid:=sub("\\.[0-9]*$", "", gene_id[1]), by=gene_id]
gtex.mediantpm.matr.melt.med.f <- gtex.mediantpm.matr.melt.med[!is.na(median.tmp.z)]
gtex.mediantpm.z <- data.frame(dcast.data.table(gtex.mediantpm.matr.melt.med.f, ensgid~tissue.simple, value.var = "median.tmp.z", fun.aggregate = median, fill=NA), row.names=1)




gtex.mediantpm.matr.melt.med.f.07 <- merge(gtex.mediantpm.matr.melt.med.f, p07.paths.probs.dt, by="ensgid")
ggplot(gtex.mediantpm.matr.melt.med.f.07[dir.sig.call!="EE"], aes(x=tissue.simple, y=median.tmp.z)) + geom_violin() + facet_wrap(~dir.pos)


row.names(p04.rlog.z) <- sub("\\.[0-9]*$", "", row.names(p04.rlog.z))
row.names(p07.rlog.z) <- sub("\\.[0-9]*$", "", row.names(p07.rlog.z))
gtex.mediantpm.matr.p04.f <- gtex.mediantpm.z[row.names(gtex.mediantpm.z)%in%row.names(p04.rlog.z), ] %>% .[rowSums(.>0)>1,]
gtex.mediantpm.matr.p07.f <- gtex.mediantpm.z[row.names(gtex.mediantpm.z)%in%row.names(p07.rlog.z), ] %>% .[rowSums(.>0)>1,]

gtex.mediantpm.matr.p04.f.merge <- cbind(gtex.mediantpm.matr.p04.f, p04.rlog.z[row.names(gtex.mediantpm.matr.p04.f),])
gtex.mediantpm.matr.p07.f.merge <- cbind(gtex.mediantpm.matr.p07.f, p07.rlog.z[row.names(gtex.mediantpm.matr.p07.f),])
plot_correlation(gtex.mediantpm.matr.p07.f.merge)

p04.paths.probs.dt.tissue.cast.f2 <- p04.paths.probs.dt.tissue.cast[row.names(p04.paths.probs.dt.tissue.cast)%in%row.names(gtex.mediantpm.matr.p04.f),p04.paths.probs.dt.tissue[, .N, by=.(x, transition)][order(x), make.names(transition)]] %>% .[, apply(., 2, uniqueN)>1]
p07.paths.probs.dt.tissue.cast.f2 <- p07.paths.probs.dt.tissue.cast[row.names(p07.paths.probs.dt.tissue.cast)%in%row.names(gtex.mediantpm.matr.p07.f),p07.paths.probs.dt.tissue[, .N, by=.(x, transition)][order(x), make.names(transition)]] %>% .[, apply(., 2, uniqueN)>1]

pheatmap::pheatmap(gtex.mediantpm.matr.p04.f, annotation_row = p04.paths.probs.dt.tissue.cast.f2, scale="row")
pheatmap::pheatmap(p07.counts[row.names(p07.gene.annot.liver.bm),], cluster_cols = FALSE, annotation_row = p07.gene.annot.liver.bm, scale="row")

# NEW Figure 5 Liver and Bone Marrow P07 ----


library(cowplot)

pid <- "P07"
tissue <- "Bone.Marrow"
lab.dat.use <- BMT.lab.data.orig[Participant.ID==pid & Lab%in%c("WBC Count")]
wbc.lab.mean <- lab.dat.use[, mean(lab.value)]
wbc.lab.sd <- lab.dat.use[, sd(lab.value)]

expr.dt.tissue[all.modules, module.group:=i.module.group, on=c("participant.ID", "ensgid")]
clusts.ngenes <- expr.dt.tissue[!is.na(module.group), .(n.genes=uniqueN(ensgid)), by=.(Tissue, participant.ID,module.group)]
clusts.ngenes[, n.clusts:=seq(Tissue), by=.(participant.ID, Tissue)]
fdr.min <- 0.01
all.modules.te[, module.group:=query.set]
sig.clusts <- all.modules.te[participant.ID==pid & pAdj<fdr.min & Tissue==tissue, .(pAdj=min(pAdj)), by=.(module.group, Tissue)]
sig.clusts[, n.clusts:=order(pAdj)]
sig.clusts[clusts.ngenes[participant.ID==pid], n.genes:=i.n.genes, on=c("Tissue", "module.group")]

g1 <- ggplot( expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z, color=gmpr.rlog.blind)) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5) +
  scale_color_gradientn(colours=heat.colors(20)) +
  stat_summary(geom="line", color="black", fun.y = "mean", linetype=2, size=0.25, aes(group=module.group)) + 
  stat_summary(geom="errorbar",color="black", fun.data="mean_cl_boot", fun.args=list(conf.int=.98), width=1, size=0.5, aes(group=module.group)) +
  stat_summary(geom="point", color="black", fun.y="mean", size=0.1, color="black", aes(group=module.group)) + 
  theme_bw(base_size = 6) + 
  geom_text(data=sig.clusts, color="black", aes(x=0, y=(0.75*n.clusts), label=paste0(module.group, " = ", n.genes))) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, linetype=2, aes(y= (lab.value - wbc.lab.mean) / wbc.lab.sd , group=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - wbc.lab.mean) / wbc.lab.sd , group=Lab)) +
  labs(x="Days from BMT", y="Abundance (z-score)", color="P07 Cluster ID")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) +
  scale_y_continuous(sec.axis = sec_axis(~ (. * wbc.lab.sd) + wbc.lab.mean , breaks=seq(-2, 20, by=2), name = "Lab Value (K/uL)"),breaks=seq(-8, 8, 0.5)) +
  facet_wrap(~module.group, ncol=1); g1
gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_BONEMARROW_Enriched_CLUSTERSSeparate_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=30*2.54, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_BONEMARROW_Enriched_CLUSTERSSeparate_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=30*2.54, units="mm")


line.cols <- sig.clusts$module.group
names(line.cols) <- sig.clusts$module.group
g1 <- ggplot(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z,  color=module.group, group=module.group)) + 
  #scale_color_gradientn(colours=heat.colors(20), aes(color=gmpr.rlog.blind)) +
  stat_summary(geom="line", fun.y = "mean", size=0.25, position = position_dodge(width=2.5)) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=3, size=0.2, position = position_dodge(width=2.5)) +
  stat_summary(geom="point", fun.y="mean", size=0.1,position = position_dodge(width=2.5)) + 
  scale_color_manual(values=line.cols) +
  theme_bw(base_size = 6) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5) +
  geom_text(data=sig.clusts, color="black", size=1, aes(x=0, y=(0.75*n.clusts), label=paste0(module.group, " = ", n.genes))) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, linetype=2, aes(y= (lab.value - wbc.lab.mean) / (wbc.lab.sd*0.75) , group=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - wbc.lab.mean) / (wbc.lab.sd*0.75) , group=Lab)) +
  labs(x="Days from BMT", y="Abundance (z-score)", color="P07 Cluster ID")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) +
  scale_y_continuous(sec.axis = sec_axis(~ (. * (wbc.lab.sd*0.75)) + wbc.lab.mean ,  name = "Lab Value (K/uL)", breaks=seq(-4,12,2)),  breaks=seq(-6,6,by=0.5), limits=c(-2.5,2.5)) + 
  facet_wrap(~Tissue, ncol=1); g1
gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_BONEMARROW_Enriched_CLUSTERS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_BONEMARROW_Enriched_CLUSTERS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")




tissue="Liver"
lab.dat.use <- copy(BMT.lab.data.orig[Participant.ID=="P07" & Lab%in%c("SGOT (AST)", "SGPT (ALT)")])
lab.dat.use[, lab.value.norm:=(lab.value-mean(lab.value))/sd(lab.value), by=Lab]
sig.clusts <- all.modules.te[participant.ID==pid & pAdj<fdr.min & Tissue==tissue, .(pAdj=min(pAdj)), by=.(module.group, Tissue)]
sig.clusts[, n.clusts:=order(pAdj)]
sig.clusts[clusts.ngenes[participant.ID==pid], n.genes:=i.n.genes, on=c("Tissue", "module.group")]
g1 <- ggplot(expr.dt.tissue[Tissue=="Liver" & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z , group=module.group)) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.25, aes(color=gmpr.rlog.blind)) +
  scale_color_gradientn(colours=heat.colors(20)) +
  stat_summary(geom="line", fun.y = "mean", size=0.5) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=1, size=0.5) +
  stat_summary(geom="point", fun.y="mean", size=0.2) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, aes(y= (lab.value - 45.5) / 57.7 , group=Lab, linetype=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - 45.5) / 57.7 , group=Lab)) +
  geom_text(data=sig.clusts[Tissue=="Liver"], color="black", aes(x=50, y=2, label=paste0("n = ", n.genes))) + 
  theme_bw(base_size = 6) +
  labs( x="Days from BMT", y="Abundance (z-score)", color="black")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) + 
  scale_y_continuous(sec.axis = sec_axis(~ (. * 57.7) + 45.5 , breaks=seq(-100, 500, by=50), name = "Lab Value (IU/L)"),  breaks=seq(-6,6,by=0.5)) + facet_wrap(~module.group); g1

gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_LiverEnriched_CLUSTERS_WITHLABS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_LiverEnriched_CCLUSTERS_WITHLABS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")

# Plot P07 liver clusts with z-scores
g1 <- ggplot(expr.dt.tissue[Tissue=="Liver" & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z , group=module.group)) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5, aes(color=gmpr.rlog.blind)) +
  scale_color_gradientn(colours=heat.colors(20)) +
  stat_summary(geom="line", fun.y = "mean", size=0.25) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=1, size=0.5) +
  stat_summary(geom="point", fun.y="mean", size=0.2) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, aes(y= lab.value.norm , group=Lab, linetype=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=lab.value.norm , group=Lab)) +
  geom_text(data=lab.dat.use[value.call=="HI"], size=2, color="black", aes(y=lab.value.norm , label=lab.value, group=Lab)) +
  theme_bw(base_size = 6) +
  geom_text(data=sig.clusts[Tissue=="Liver"], aes(x=-7, y=2.5, label=paste0("n=", n.genes))) +
  labs( x="Days from BMT", y="Abundance (z-score)")  + 
  scale_x_continuous(breaks=expr.dt.tissue[Tissue=="Liver" & participant.ID==pid  & module.group%in%sig.clusts$module.group][, .N, by=Day.num]$Day.num) + 
  scale_y_continuous(breaks=seq(-4,4,0.5)) + facet_wrap(~module.group); g1

gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_LiverEnriched_CCLUSTERS_WITH_ZNORM_LABS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")

gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P07", "_LiverEnriched_CLUSTERS_WITH_ZNORM_LABS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")


####
# Figure 5 Liver and Bone Marrow P04 ---

pid <- "P04"
tissue <- "Bone.Marrow"
lab.dat.use.tmp <- BMT.lab.data.orig[Participant.ID==pid & Lab%in%c("WBC Count")]
wbc.lab.mean <- lab.dat.use.tmp[, mean(lab.value)]
wbc.lab.sd <- lab.dat.use.tmp[, sd(lab.value)]
lab.dat.use <- lab.dat.use.tmp[as.numeric(Day.num)%between%(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group, range(Day.num)+c(-1,1)])]
fdr.min <- 0.01

sig.clusts <- all.modules.te[participant.ID==pid & pAdj<fdr.min & Tissue==tissue, .(pAdj=min(pAdj)), by=.(module.group, Tissue)]
sig.clusts[, n.clusts:=rank(pAdj)]
sig.clusts[clusts.ngenes[participant.ID==pid], n.genes:=i.n.genes, on=c("Tissue", "module.group")]

g1 <- ggplot(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z,  group=module.group)) + 
  scale_color_gradientn(colours=heat.colors(20), aes(color=gmpr.rlog.blind)) +
  stat_summary(geom="line", fun.y = "mean", size=0.5) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=1, size=0.5) +
  stat_summary(geom="point", fun.y="mean", size=0.1) + 
  theme_bw(base_size = 6) + 
  ggbeeswarm::geom_quasirandom(size=0.2, alpha=0.5) +
  geom_text(data=sig.clusts, color="black", aes(x=0, y=(1.5*n.clusts), label=paste0(module.group, " = ", n.genes))) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, linetype=2, aes(y= (lab.value - wbc.lab.mean) / wbc.lab.sd , group=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - wbc.lab.mean) / wbc.lab.sd , group=Lab)) +
  labs(x="Days from BMT", y="Abundance (z-score)", color="P04 Cluster ID")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) +
  scale_y_continuous(sec.axis = sec_axis(~ (. * wbc.lab.sd) + wbc.lab.mean ,  name = "Lab Value (K/uL)"),  breaks=seq(-6,6,by=0.5)) + 
  facet_wrap(~module.group, ncol=1); g1
gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_BONEMARROW_Enriched_CLUSTERS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_BONEMARROW_Enriched_CLUSTERS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")

line.cols <- sig.clusts$module.group
names(line.cols) <- sig.clusts$module.group
g1 <- ggplot(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z,  color=module.group, group=module.group)) + 
  #scale_color_gradientn(colours=heat.colors(20), aes(color=gmpr.rlog.blind)) +
  stat_summary(geom="line", fun.y = "mean", size=0.25, position = position_dodge(width=2.5)) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=3, size=0.2, position = position_dodge(width=2.5)) +
  stat_summary(geom="point", fun.y="mean", size=0.1,position = position_dodge(width=2.5)) + 
  scale_color_manual(values=line.cols) +
  theme_bw(base_size = 6) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5) +
  geom_text(data=sig.clusts, color="black", size=1, aes(x=0, y=(0.75*n.clusts), label=paste0(module.group, " = ", n.genes))) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, linetype=2, aes(y= (lab.value - wbc.lab.mean) / (wbc.lab.sd*0.75) , group=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - wbc.lab.mean) / (wbc.lab.sd*0.75) , group=Lab)) +
  labs(x="Days from BMT", y="Abundance (z-score)", color="P04 Cluster ID")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) +
  scale_y_continuous(sec.axis = sec_axis(~ (. * (wbc.lab.sd*0.75)) + wbc.lab.mean ,  name = "Lab Value (K/uL)", breaks=seq(-4,12,2)),  breaks=seq(-6,6,by=0.5), limits=c(-2.5,2.5)) + 
  facet_wrap(~Tissue, ncol=1); g1
gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_BONEMARROW_Enriched_CLUSTERS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_BONEMARROW_Enriched_CLUSTERS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")


# P04 liver with lab value not z-transformed
tissue="Liver"
lab.dat.use <- BMT.lab.data.orig[Participant.ID==pid & Lab%in%c("SGOT (AST)", "SGPT (ALT)") & as.numeric(Day.num) <= expr.dt.tissue[participant.ID==pid, max(Day.num)+1]]
liver.lab.mean <- lab.dat.use[, mean(lab.value)]
liver.lab.sd <- lab.dat.use[, sd(lab.value)]

sig.clusts <- all.modules.te[participant.ID==pid & pAdj<fdr.min & Tissue==tissue, .(pAdj=min(pAdj)), by=.(module.group, Tissue)]
sig.clusts[, n.clusts:=rank(pAdj)]
sig.clusts[clusts.ngenes[participant.ID==pid], n.genes:=i.n.genes, on=c("Tissue", "module.group")]

g1 <- ggplot(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z , group=module.group)) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5, aes(color=gmpr.rlog.blind)) +
  scale_alpha_continuous() +
  stat_summary(geom="line", color=sig.clusts$module.group,fun.y = "mean", size=0.25, position = position_dodge(width=2.5)) + 
  stat_summary(geom="errorbar", color=sig.clusts$module.group,fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=2, size=0.2, position = position_dodge(width=2.5)) +
  stat_summary(geom="point", fun.y="mean", size=0.1,color=sig.clusts$module.group,position = position_dodge(width=2.5)) + 
  geom_line(data=lab.dat.use, color="black", size=0.2, aes(y= (lab.value - liver.lab.mean) / liver.lab.sd , group=Lab, linetype=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=(lab.value - liver.lab.mean) / liver.lab.sd , group=Lab)) +
  theme_bw(base_size = 6) +
  geom_text(data=sig.clusts[Tissue=="Liver"], size=1, aes(x=-7, y=2.5, label=paste0(module.group, " = ", n.genes))) +
  labs( x="Days from BMT", y="Abundance (z-score)")  + 
  scale_x_continuous(breaks=seq(-14, 70, 7)) + 
  scale_y_continuous(sec.axis = sec_axis(~ (. * liver.lab.sd) + liver.lab.mean , breaks=seq(-100, 500, by=25), name = "Lab Value (IU/L)"),  breaks=seq(-6,6,by=0.5)) + facet_wrap(~Tissue); g1

gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_LiverEnriched_CLUSTERS_WITHLABS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")
gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_LiverEnriched_CCLUSTERS_WITHLABS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")


# P04 liver with lab value z-transformed

lab.dat.use[, lab.value.norm:=(lab.value-mean(lab.value))/sd(lab.value), by=Lab]
g1 <- ggplot(expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group], aes(x=Day.num, y=gmpr.rlog.blind.z , group=module.group)) + 
  ggbeeswarm::geom_quasirandom(size=0.1, alpha=0.5, aes(color=gmpr.rlog.blind)) +
  scale_color_gradientn(colours=heat.colors(20)) +
  stat_summary(geom="line", fun.y = "mean", size=0.25, position = position_dodge(width=2.5)) + 
  stat_summary(geom="errorbar", fun.data="mean_cl_boot", fun.args=list(conf.int=0.98), width=2, size=0.2, position = position_dodge(width=2.5)) +
  stat_summary(geom="point", fun.y="mean", size=0.1,position = position_dodge(width=2.5)) +
  geom_line(data=lab.dat.use, color="black", size=0.2, aes(y= lab.value.norm , group=Lab, linetype=Lab)) +
  geom_point(data=lab.dat.use, color="black", size=0.2, aes(y=lab.value.norm , group=Lab)) +
  geom_text(data=lab.dat.use[value.call=="HI"], size=2, color="black", aes(y=lab.value.norm , label=lab.value, group=Lab)) +
  theme_bw(base_size = 6) +
  geom_text(data=sig.clusts[Tissue=="Liver"], aes(x=-7, y=2.5, label=paste0("n=", n.genes))) +
  labs( x="Days from BMT", y="Abundance (z-score)")  + 
  scale_x_continuous(breaks=expr.dt.tissue[Tissue==tissue & participant.ID==pid  & module.group%in%sig.clusts$module.group, .N, by=Day.num]$Day.num) + 
  scale_y_continuous(breaks=seq(-4,4,0.5)) + facet_wrap(~ module.group); g1

gL <- g1 + theme(legend.position = "none", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_LiverEnriched_CCLUSTERS_WITH_ZNORM_LABS_NOLEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")

gL <- g1 + theme(legend.position = "top", axis.text=element_text(size=6), panel.grid = element_blank()); gL
file.out <- paste0("./output/figures/TissueSpecificPlot_", "P04", "_LiverEnriched_CLUSTERS_WITH_ZNORM_LABS_LEGEND.pdf")
save_plot(filename = file.out, plot = gL, base_width = 75, base_height=45, units="mm")




# Plot heatmaps for sig genes ----
sig.clusters.all <- all.modules.te[pAdj<fdr.min, .(pAdj=min(pAdj)), by=.(module.group, Tissue, participant.ID)]

p04.rlog.m <- log2(p04.counts.sig+1)
p07.rlog.m <- log2(p07.counts.sig+1)
row.names(p07.rlog.m) <- sub("\\.[0-9]*$", "", row.names(p07.rlog.m))
row.names(p04.rlog.m) <- sub("\\.[0-9]*$", "", row.names(p04.rlog.m))

colors07mod <- expr.dt.tissue[participant.ID=="P07", .N, by=module.group]$module.group
names(colors07mod) <- colors07mod
colors04mod <- expr.dt.tissue[participant.ID=="P04", .N, by=module.group]$module.group
names(colors04mod) <- colors04mod
ann_colors07 = list(
  module.group =colors07mod
)
ann_colors04 = list(
  module.group =colors04mod
)
dis.m.04 <- as.matrix(p04.modules.all$disTOM)
dis.m.07 <- as.matrix(p07.modules.all$disTOM)
qnames.merged[, plot.id:=ifelse(ol.symbol=="", ensgid, ol.symbol)]
dimnames(dis.m.04) <- list(row.names(p04.counts.sig), row.names(p04.counts.sig))
dimnames(dis.m.07) <- list(row.names(p07.counts.sig), row.names(p07.counts.sig))

p04.rlog.m[p04.rlog.m==0] <- NA
p07.rlog.m[p07.rlog.m==0] <- NA
sig.clusts.all.expr <- merge(expr.dt.tissue, sig.clusters.all, by=c("module.group", "Tissue", "participant.ID"))
sig.clusts.all.expr.cast <- dcast.data.table(sig.clusts.all.expr, ensgid+module.group+participant.ID~Tissue, fun.aggregate = uniqueN, value.var="module.group")

expr.dt.liver <- copy(expr.dt.tissue[Tissue=="Liver" & !is.na(module.group)])
expr.dt.bm <- copy(expr.dt.tissue[Tissue=="Bone.Marrow" & !is.na(module.group)])

p04.gene.annot.liver <- data.frame(expr.dt.tissue[participant.ID=="P04" & Tissue=="Liver" & !is.na(module.group), .(mean.cpm=mean(gmpr.rlog.blind)), by=.(ensgid, module.group)], row.names=1)
p07.gene.annot.liver <- data.frame(expr.dt.tissue[participant.ID=="P07" & Tissue=="Liver" & !is.na(module.group), .(mean.cpm=mean(gmpr.rlog.blind)), by=.(ensgid, module.group)], row.names=1)
p04.gene.annot.bm <- data.frame(expr.dt.tissue[participant.ID=="P04" & Tissue=="Bone.Marrow" & !is.na(module.group), .(mean.cpm=mean(gmpr.rlog.blind)), by=.(ensgid, module.group)], row.names=1)
p07.gene.annot.bm <- data.frame(expr.dt.tissue[participant.ID=="P07" & Tissue=="Bone.Marrow" & !is.na(module.group), .(mean.cpm=mean(gmpr.rlog.blind)), by=.(ensgid, module.group)], row.names=1)


pheatmap::pheatmap(p04.rlog.m[row.names(p04.gene.annot.liver),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p04.gene.annot.liver, scale="row", annotation_colors=ann_colors04, clustering_distance_rows = as.dist(dis.m.04.f), labels_row = qnames.merged[row.names(p04.rlog.m[row.names(p04.gene.annot.liver),]), on="ensgid", mult="first"]$plot.id, labels_col=p04.dds$BMT.Day)

heatmap.dist.list <- list("P07.Liver" = dis.m.07[row.names(p07.gene.annot.liver),row.names(p07.gene.annot.liver)],
     "P04.Liver" =  dis.m.04[row.names(p04.gene.annot.liver),row.names(p04.gene.annot.liver)],
     "P04.Bone.Marrow" = dis.m.04[row.names(p04.gene.annot.bm),row.names(p04.gene.annot.bm)],
     "P07.Bone.Marrow" =  dis.m.07[row.names(p07.gene.annot.bm),row.names(p07.gene.annot.bm)])
# P07 Heatmaps
outfl.07 <- "./output/figures/NEW_FIG5_Heatmap_P07_Liver_NOLEGEND.pdf"
pheatmap::pheatmap(p07.rlog.m[row.names(p07.gene.annot.liver),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p07.gene.annot.liver, scale="row", annotation_colors=ann_colors07, clustering_distance_rows = as.dist(heatmap.dist.list[["P07.Liver"]]), labels_row = qnames.merged[row.names(p07.rlog.m[row.names(p07.gene.annot.liver),]), on="ensgid", mult="first"]$plot.id, labels_col=p07.dds$BMT.Day,
filename = outfl.07, legend = FALSE, annotation_legend = FALSE,fontsize = 5, fontsize_row=6, treeheight_row = 15, fontsize_col = 6, width = 3, height=3, main="P07: Liver")

outfl.07 <- "./output/figures/NEW_FIG5_Heatmap_P07_Liver_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(p07.rlog.m[row.names(p07.gene.annot.liver),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p07.gene.annot.liver, scale="row", annotation_colors=ann_colors07, clustering_distance_rows = as.dist(heatmap.dist.list[["P07.Liver"]]), labels_row = qnames.merged[row.names(p07.rlog.m[row.names(p07.gene.annot.liver),]), on="ensgid", mult="first"]$plot.id, labels_col=p07.dds$BMT.Day,
                   filename = outfl.07, legend = TRUE, annotation_legend = TRUE,fontsize = 5, fontsize_row=6, fontsize_col = 6, width = 3, height=2.5, main="P07: Liver")

outfl.07 <- "./output/figures/NEW_FIG5_Heatmap_P07_BoneMarrow_NOLEGEND.pdf"
pheatmap::pheatmap(p07.rlog.m[row.names(p07.gene.annot.bm),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p07.gene.annot.bm, scale="row", annotation_colors=ann_colors07, clustering_distance_rows = as.dist(heatmap.dist.list[["P07.Bone.Marrow"]]), labels_row = qnames.merged[row.names(p07.rlog.m[row.names(p07.gene.annot.bm),]), on="ensgid", mult="first"]$plot.id, labels_col=p07.dds$BMT.Day,
                   filename = outfl.07, legend = FALSE, annotation_legend = FALSE,fontsize = 5, fontsize_row=6, treeheight_row = 15, fontsize_col = 6, width = 3, height=3, main="P07: BoneMarrow")

outfl.07 <- "./output/figures/NEW_FIG5_Heatmap_P07_BoneMarrow_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(p07.rlog.m[row.names(p07.gene.annot.bm),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p07.gene.annot.bm, scale="row", annotation_colors=ann_colors07, clustering_distance_rows = as.dist(heatmap.dist.list[["P07.Bone.Marrow"]]), labels_row = qnames.merged[row.names(p07.rlog.m[row.names(p07.gene.annot.bm),]), on="ensgid", mult="first"]$plot.id, labels_col=p07.dds$BMT.Day,
                   filename = outfl.07, legend = TRUE, annotation_legend = TRUE,fontsize = 5, fontsize_row=6, fontsize_col = 6, width = 3, height=2.5, main="P07: BoneMarrow")


# P04 heamaps ----
outfl.04 <- "./output/figures/NEW_FIG5_Heatmap_P04_Liver_NOLEGEND.pdf"
pheatmap::pheatmap(p04.rlog.m[row.names(p04.gene.annot.liver),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p04.gene.annot.liver, scale="row", annotation_colors=ann_colors04, clustering_distance_rows = as.dist(heatmap.dist.list[["P04.Liver"]]), labels_row = qnames.merged[row.names(p04.rlog.m[row.names(p04.gene.annot.liver),]), on="ensgid", mult="first"]$plot.id, labels_col=p04.dds$BMT.Day,
                   filename = outfl.04, legend = FALSE, annotation_legend = FALSE,fontsize = 5, fontsize_row=6, treeheight_row = 15, fontsize_col = 6, width = 3, height=2.5, main="P04: Liver")

outfl.04 <- "./output/figures/NEW_FIG5_Heatmap_P04_Liver_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(p04.rlog.m[row.names(p04.gene.annot.liver),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p04.gene.annot.liver, scale="row", annotation_colors=ann_colors04, clustering_distance_rows = as.dist(heatmap.dist.list[["P04.Liver"]]), labels_row = qnames.merged[row.names(p04.rlog.m[row.names(p04.gene.annot.liver),]), on="ensgid", mult="first"]$plot.id, labels_col=p04.dds$BMT.Day,
                   filename = outfl.04, legend = TRUE, annotation_legend = TRUE,fontsize = 5, fontsize_row=6, fontsize_col = 6, width = 3, height=2.5, main="P04: Liver")

outfl.04 <- "./output/figures/NEW_FIG5_Heatmap_P04_BoneMarrow_NOLEGEND.pdf"
pheatmap::pheatmap(p04.rlog.m[row.names(p04.gene.annot.bm),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p04.gene.annot.bm, scale="row", annotation_colors=ann_colors04, clustering_distance_rows = as.dist(heatmap.dist.list[["P04.Bone.Marrow"]]), labels_row = qnames.merged[row.names(p04.rlog.m[row.names(p04.gene.annot.bm),]), on="ensgid", mult="first"]$plot.id, labels_col=p04.dds$BMT.Day,
                   filename = outfl.04, legend = FALSE, annotation_legend = FALSE,fontsize = 5, fontsize_row=6, treeheight_row = 15, fontsize_col = 6, width = 3, height=3, main="P04: BoneMarrow")

outfl.04 <- "./output/figures/NEW_FIG5_Heatmap_P04_BoneMarrow_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(p04.rlog.m[row.names(p04.gene.annot.bm),], na_col = "darkgrey",cluster_cols = FALSE, annotation_row = p04.gene.annot.bm, scale="row", annotation_colors=ann_colors04, clustering_distance_rows = as.dist(heatmap.dist.list[["P04.Bone.Marrow"]]), labels_row = qnames.merged[row.names(p04.rlog.m[row.names(p04.gene.annot.bm),]), on="ensgid", mult="first"]$plot.id, labels_col=p04.dds$BMT.Day,
                   filename = outfl.04, legend = TRUE, annotation_legend = TRUE,fontsize = 5, fontsize_row=6, fontsize_col =6, width = 3, height=2.5, main="P04: BoneMarrow")

# Plot multi-tissue cluster as heatmap ----

clusts.07.ngenes <- expr.dt.tissue[sig.clusts, on="module.group"]
clusts.07.ngenes[, Source.Set:=ifelse(GTEx>0 & ProteinAtlas>0, "BOTH", ifelse(GTEx>0, "GTEx", ifelse(ProteinAtlas>0, "ProteinAtlas", "NONE")))]
clusts.07.ngenes[, Tissue.Set:=paste0(Tissue, ".", Source.Set)]
tissue.row.dat <- dcast.data.table(clusts.07.ngenes[!is.na(Tissue)], ensgid~Tissue, value.var=c("module.group"), fun.aggregate = function(x){uniqueN(x)}, fill=0)
tissue.row.dat[, tissue.tot:=rowSums((.SD)), .SDcols=clusts.07.ngenes[!is.na(Tissue), .N, by=Tissue]$Tissue]
tissue.tots <- tissue.row.dat[, colSums((.SD)), .SDcols=clusts.07.ngenes[!is.na(Tissue), .N, by=Tissue]$Tissue]
tissue.tots <- tissue.tots[order(-tissue.tots)]
tissue.row.dat[, ensgid:=as.character(ensgid)]
tissue.row.dat[, mean.log.cpm:=rowMeans(rlog.cpm.m)[ensgid]] 
setorderv(tissue.row.dat, c("mean.log.cpm", "tissue.tot", names(tissue.tots)), rep.int(-1, length(tissue.tots)+2))
row.dat07 <- data.frame(tissue.row.dat, row.names=1)[, c("mean.log.cpm", names(tissue.tots))]
p07.gmp.log2.cpm.f.liverclust <- gmp.list.log2.for.mfuzz.nonstand[["P07.HiSeq"]][as.character(tissue.row.dat$ensgid),]
# FIG 5 heatmap liver multi tissue p07 ----
outfl.07 <- "./output/figures/FIG5_Heatmap_P07_Liver_MULTITISSUE_CLUSTERNOLEGEND.pdf"
p07.gmp.log2.cpm.f.liverclust2 <- p07.gmp.log2.cpm.f.liverclust
p07.gmp.log2.cpm.f.liverclust2[p07.gmp.log2.cpm.f.liverclust2==0] <- NA
pheatmap::pheatmap(t(p07.gmp.log2.cpm.f.liverclust2[row.names(row.dat07),]), na_col = "darkgrey",
                   annotation_col = row.dat07, 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   scale="column",
                   labels_col = qnames.merged[row.names(row.dat07), ifelse(ol.symbol=="", ensgid, ol.symbol), on="ensgid"], 
                   labels_row = as.character(colData(sceset.p07.06.04)[colnames(p07.gmp.log2.cpm.f.liverclust),]$BMT.Day),
                   filename = outfl.07, legend = FALSE, annotation_legend = FALSE, main =paste0("CLUSTER ",liver.multitissue.clusts), fontsize = 5, fontsize_row=6, fontsize_col = 4, width = 5, height=2.5)
outfl.07 <- "./output/figures/FIG5_Heatmap_P07_Liver_MULTITISSUE_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(p07.gmp.log2.cpm.f.liverclust2[row.names(row.dat07),],  na_col = "darkgrey",
                   annotation_row = row.dat07, 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   labels_row = qnames.merged[row.names(row.dat07), ifelse(ol.symbol=="", ensgid, ol.symbol), on="ensgid"], 
                   labels_col = as.character(colData(sceset.p07.06.04)[colnames(p07.gmp.log2.cpm.f.liverclust),]$BMT.Day),
                   filename = outfl.07, legend = TRUE, annotation_legend = TRUE, main = paste0("CLUSTER ",liver.multitissue.clusts), fontsize = 5, fontsize_row=4, fontsize_col = 6, width = 3, height=6)

# FIG 5 heatmap bmt multi tissue p07 ----

bonemarrow.clusts07 <- sig.clusts[Tissue=="Bone.Marrow", .N, by=query.set]$query.set
clusts.07.ngenes <- dt.use07.f[bonemarrow.clusts07, on="query.set"]
clusts.07.ngenes[, Source.Set:=ifelse(GTEx>0 & ProteinAtlas>0, "BOTH", ifelse(GTEx>0, "GTEx", ifelse(ProteinAtlas>0, "ProteinAtlas", "NONE")))]
clusts.07.ngenes[, Tissue.Set:=paste0(Tissue, ".", Source.Set)]
tissue.row.dat <- dcast.data.table(clusts.07.ngenes, ensgid+query.set+membership.value~Tissue.Set, value.var=c("Group"), fun.aggregate = function(x){unique(x)}, fill=0)
tissue.row.dat[, ensgid:=as.character(ensgid)]
p07.gmp.log2.cpm.f.bonemarrowclust <- gmp.list.log2.for.mfuzz.nonstand[["P07.HiSeq"]][as.character(tissue.row.dat$ensgid),]
p07.gmp.log2.cpm.f.bonemarrowclust[p07.gmp.log2.cpm.f.bonemarrowclust==0] <- NA
tissue.row.dat[, mean.log.cpm:=rowMeans(p07.gmp.log2.cpm.f.bonemarrowclust, na.rm = TRUE)[ensgid]] 
setorderv(tissue.row.dat, c("query.set", "membership.value"), c(1, -1))
row.dat07 <- data.frame(tissue.row.dat, row.names=1)[, c("mean.log.cpm", "membership.value", "Bone.Marrow.ProteinAtlas", "query.set")]
setnames(row.dat07, "Bone.Marrow.ProteinAtlas", "Bone.Marrow" )
outfl.07 <- "./output/figures/FIG5_Heatmap_P07_BONEMARROW_CLUSTERNOLEGEND.pdf"
pheatmap::pheatmap(t(p07.gmp.log2.cpm.f.bonemarrowclust[row.names(row.dat07),]), gaps_col = which(!duplicated(row.dat07$query.set)[-1]), na_col = "darkgrey",
                   annotation_col = row.dat07, 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   scale="column",
                   filename = outfl.07, 
                   labels_col = qnames.merged[row.names(row.dat07), ifelse(ol.symbol=="", ensgid, ol.symbol), on="ensgid"], 
                   labels_row = as.character(colData(sceset.p07.06.04)[colnames(p07.gmp.log2.cpm.f.bonemarrowclust),]$BMT.Day),
                   legend = FALSE, annotation_legend = FALSE, main =paste0("CLUSTER ",bonemarrow.clusts07), fontsize = 5, fontsize_row=6, fontsize_col = 4, width = 5, height=3)
outfl.07 <- "./output/figures/FIG5_Heatmap_P07_BONEMARROW_CLUSTERNOLEGEND_TALL.pdf"
pheatmap::pheatmap((p07.gmp.log2.cpm.f.bonemarrowclust[row.names(row.dat07),]), gaps_row = which(!duplicated(row.dat07$query.set)[-1]), na_col = "darkgrey",
                   annotation_row = row.dat07, 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   scale="row",
                   filename = outfl.07, 
                   labels_row = qnames.merged[row.names(row.dat07), ifelse(ol.symbol=="", ensgid, ol.symbol), on="ensgid"], 
                   labels_col = as.character(colData(sceset.p07.06.04)[colnames(p07.gmp.log2.cpm.f.liverclust),]$BMT.Day),
                   legend = FALSE, annotation_legend = FALSE,  fontsize = 5, fontsize_row=6, fontsize_col = 4, width = 2.5, height=3.5)
outfl.07 <- "./output/figures/FIG5_Heatmap_P07_BONEMARROW_CLUSTER_LEGEND.pdf"
pheatmap::pheatmap(t(p07.gmp.log2.cpm.f.bonemarrowclust[row.names(row.dat07),]), gaps_col = which(!duplicated(row.dat07$query.set)[-1]),
                   annotation_col = row.dat07, 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   scale="column",
                   filename = outfl.07, 
                   labels_col = qnames.merged[row.names(row.dat07), ifelse(ol.symbol=="", ensgid, ol.symbol), on="ensgid"], 
                   labels_row = as.character(colData(sceset.p07.06.04)[colnames(p07.gmp.log2.cpm.f.liverclust),]$BMT.Day),
                   legend = TRUE, annotation_legend = TRUE, main =paste0("CLUSTER ",bonemarrow.clusts07), fontsize = 5, fontsize_row=6, fontsize_col = 4, width = 5, height=3)



#
library(RCAS)
gtf.loc="/data/genomic_data/Homo_sapiens/UCSC/GRCh38/ANNOTATIONS/GENCODE/v27/gencode.v27.basic.annotation.gtf"
gtf.tmp <- paste0("/dev/shm/", basename(gtf.loc))
if(!file.exists(gtf.tmp)){
  system(paste0("cp ", gtf.loc, " ", gtf.tmp))
}
gtf <- importGtf(gtf.tmp)
txdbFeatures <- getTxdbFeaturesFromGRanges(gtf)
txdbFeatures.dt <- lapply(txdbFeatures, FUN=function(x){dt <- as.data.table(x); setkey(dt, seqnames, strand, start, end)})
gtf.dt <- as.data.table(gtf)
gtf.dt[, gene_id:=as.character(gene_id)]
gene.types.dt <- gtf.dt[, .N, by=gene_type]
keep.features <- c("protein_coding", "processed_transcript", "lincRNA", "3prime_overlapping_ncRNA", "antisense", "non_coding", 
                   "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncRNA", 
                   "macro_lncRNA")
mRNA.lincRNA.txids <- gtf.dt[gene_type%in%keep.features & transcript_type%in%keep.features, .N, by=transcript_id]
mRNA.lincRNA.geneids <- gtf.dt[gene_type%in%keep.features & transcript_type%in%keep.features, sub("\\.[0-9]*$", "", gene_id[1]), by=gene_id]$V1
gtf.mRNA.lincRNA <- gtf[gtf$transcript_id%in%mRNA.lincRNA.txids$transcript_id]
mRNA.lincRNA.txdbFeatures <- lapply(txdbFeatures, FUN=function(x) x[x$tx_name%in%gtf.mRNA.lincRNA$transcript_id])
mRNA.lincRNA.txdbFeatures.dt <- lapply(mRNA.lincRNA.txdbFeatures, FUN=function(x){dt <- as.data.table(x); setkey(dt, seqnames, start, end)})

setkey(tissue.specific.data, ensgid, source.data, Tissue)
tissue.specific.data.groupgene <- dcast.data.table(tissue.specific.data, ensgid~source.data, value.var=c("Tissue", "Group"), fun.aggregate = paste, collapse=",", fill=NA)

expr.dt.tissue[qnames.merged, `:=`(ol.symbol=i.ol.symbol, GENENAME=i.GENENAME, ENTREZID=i.ENTREZID), on="ensgid"]
setkey(expr.dt.tissue, ensgid)
expr.dt.tissue[, ol.geneid:=as.character(ol.geneid)]
p04.sig.genes <- intersect(row.names(GeneDECalls04), mRNA.lincRNA.geneids)
p07.sig.genes <- intersect(row.names(GeneDECalls07), mRNA.lincRNA.geneids)
expr.dt.tissue[, mean.log2.count.gene:=mean(gmpr.log2.count[gmpr.log2.count>0]), by=.(ensgid, participant.ID)]
expr.dt.tissue.tg <- merge(expr.dt.tissue, tissue.specific.data.groupgene, by="ensgid")
expr.dt.tissue.tg[!is.na(module.group), WGCNA.cluster:=paste0(participant.ID, ".", module.group)]
expr.dt.tissue.summary.list <- list("P04"=dcast.data.table(expr.dt.tissue.tg[participant.ID=="P04"][p04.sig.genes, on="ensgid", nomatch=0], 
                                                           ol.geneid+ensgid+ol.symbol+GENENAME+ENTREZID+Tissue_GTEx+Group_GTEx+Tissue_ProteinAtlas+Group_ProteinAtlas+WGCNA.cluster+mean.log2.count.gene~BMT.Day,
                                                           value.var="gmpr.counts", fun.aggregate = max, fill=0.0),
                                    "P07"=dcast.data.table(expr.dt.tissue.tg[participant.ID=="P07"][p07.sig.genes, on="ensgid", nomatch=0], 
                                                           ol.geneid+ensgid+ol.symbol+GENENAME+ENTREZID+Tissue_GTEx+Group_GTEx+Tissue_ProteinAtlas+Group_ProteinAtlas+WGCNA.cluster+mean.log2.count.gene~BMT.Day,
                                                           value.var="gmpr.counts", fun.aggregate = max, fill=0.0))

invisible(lapply(expr.dt.tissue.summary.list, function(dt){
  setorderv(dt, c("WGCNA.cluster", "mean.log2.count.gene"), c(1,-1)) 
}))

WriteXLS::WriteXLS(x=expr.dt.tissue.summary.list,
                   ExcelFileName="output/Table_S5_BMT_sig_genes_and_modules.xls", verbose=TRUE)
