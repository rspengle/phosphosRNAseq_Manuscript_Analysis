# Ran STAR vs hg38 using the full, unaligned data. Now we'll compare the different stages to see how the performance improved.
# The goal will be to demonstrate the improvements gained in filtering at the different stages.

# Load libraries ----
library(data.table)
library(rbamtools)
library(Biostrings)
library(Rsamtools)
library(GenomicAlignments)
library(ChIPseeker)
library(Matrix)
library(igraph)
library(csaw)
library(Hmisc)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(dplyr)
library(ggbeeswarm)
library(GenomicFeatures)
library(ggExtra)
library(annotatr)
library(biovizBase)
library(ca)
library(cowplot)
library(enrichplot)
library(FactoMineR)
library(ggbio)
library(ggpubr)
library(grid)
library(Gviz)
library(MASS)
library(RCAS)
library(vcd)

# Set global options ----
setDTthreads(20) # sets numer of threads to use with data.table. Change to number suitible for your machine

# Load raw data ----
base.dir="./data/Plasma/healthyControlsAnalysis/biofragmenta-v1.5_20180220_103636"
sample.info.rdata.file <- "ALL.SAMPLE.INFO"
alignment.output.dir <- paste0(base.dir, "/02b-star_alignment_all")
sample.basenames.file <- dir(alignment.output.dir, pattern = "new_file_bn.tmp", full.names=TRUE)
load(sample.info.rdata.file)

# Gather all bam files from the subdirectories in base.dir. Bam files named and listed by subdir
bamfiles.list <- sapply(dir(base.dir, full.names=TRUE, include.dirs=TRUE, recursive=FALSE), dir, pattern="*bam$", full.names=TRUE)
bamfiles.list <- bamfiles.list[sapply(bamfiles.list, length)>0]
names(bamfiles.list) <- basename(names(bamfiles.list))






# OLD ----



#bamfiles.list.dt <- rbindlist(lapply(seq_along(bamfiles.list), FUN=function(X){ lst <- bamfiles.list[[X]]; data.table(dirname=rep(basename(names(bamfiles.list)[X]), length(lst)), file=lst)}))

#Sys.setenv(TCL_LIBRARY="/home/ryanspen/miniconda3/envs/BioSandbox/lib/R/share/tcl8.6")

# Functions ----
flag <- scanBamFlag(isSecondaryAlignment = FALSE)
what.info <- c("qname","rname", "strand", "pos", "qwidth", "seq", "cigar", "seq")
tags <- c("NH", "HI", "AS", "NM")
param <- ScanBamParam(flag = flag, tag = tags, what = what.info, reverseComplement = TRUE)
bam.to.dt.f <- function(bam.file, bam.params=param){
  readsGA <- readGAlignments(bam.file, param=param)
  dt <- as.data.table(readsGA)
  dt <- dt[grepl("N")]
  return(dt)
}
param2 <- readParam(pe="none", forward=NULL)
#x <- strandedCounts(bamfiles.list$`05-filtered`, width=30, filter = 1, param=param2)
#x.ranges <- rowRanges(x)
#bin <- createControlRegions(rowRanges(x))
#filter.stat <- filterWindows(x,  bin, type="local")
#min.fc <- 3
#keep <- filter.stat$filter > log2(min.fc)
#summary(keep)


fread.frwrite.choice.f <- function(file.paths, r.obj.names, read.file=TRUE, ...){
  if(length(file.paths)!=length(r.obj.names)){
    stop("Must be one file path for each r object")
  }
  for(i in seq_along(file.paths)){
    fl <- file.paths[i]
    r.obj <- r.obj.names[i]
    if(read.file==TRUE){
      message(paste0("reading file", fl))
      dt <- fread(file=fl, sep="\t", header=TRUE)
      assign(x=r.obj, value = dt, envir = .GlobalEnv)
    } else if(read.file==FALSE){
      message(paste0("Writing ", r.obj, " to ", fl))
      fwrite(x = get(r.obj), file = fl, quote = FALSE, sep="\t", row.names=FALSE, col.names=TRUE)  
    } else{
      stop("read.file must be TRUE or FALSE")
    }
  }
}


sparseM.from.dt <- function(dt, i, j, x=1, make.binary=TRUE){
  if(is.numeric(i) | is.numeric(j)){
    stop("i and j must be valid column names in dt")
  }
  if(is.numeric(x)){
    if(length(x)>1){
      stop("x must be a column of dt or left blank")
    } else if(length(x)==1){
      x.val <- x
    }
    
  } else{
    x.val <- dt[[x]]
  }
  i.lvl <- dt[, .N, keyby=eval(i)][[1]]
  j.lvl <- dt[, .N, keyby=eval(j)][[1]]
  
  s <- sparseMatrix(i=as.numeric(factor(dt[[i]], levels=i.lvl)),
                    j=as.numeric(factor(dt[[j]], levels=j.lvl)),
                    x=x.val)
  
  if(make.binary==TRUE){
    s@x <- rep.int(1, nnzero(s))
  }
  dimnames(s) <- list(i.lvl, j.lvl)
  return(s)
}

to.igraph <- function(x, i, j){
  if(is(x, "sparseMatrix") | is(x, "matrix")){
    return(graph_from_incidence_matrix(x))
  } else if(is(x, "data.frame")){
    if(is.null(i) | is.null(j)){
      stop("Must supply i or j with data.table/frame")
    } else if(!is(x, "data.table")){
      s <- sparseM.from.dt(as.data.table(x), i, j)
      return(graph_from_incidence_matrix(s))
    } else{
      s <- sparseM.from.dt(as.data.table(x), i, j)
      return(graph_from_incidence_matrix(s))
    }
  }
}

# Summarizes overlaps as in RCAS, but using data.table foverlaps. 
# q.colnames and txdb.colnames must take the format shown (chrom, <strand>, start, end), with strand being optional. 
# If the strand column name is given and ignore.strand=TRUE, then the second name is dropped.
# If the strand column is not given and ignore.strand==FALSE, then a warning is given but it moves on.
summarizeQueryRegions.dt <- function(queryRegions, txdbFeatures, q.colnames=c("seqnames", "strand", "start", "end"), txdb.colnames=c("seqnames", "strand", "start", "end"), check.orientation=FALSE, ignore.strand=FALSE, val.cols=NULL, by.col=NULL, output="SUMMARY"){
  if(length(output)>1 | any(!output%in%c("SUMMARY", "FULL"))){
    stop("<output> must be a single choice of SUMMARY or FULL")
  }
  if(!length(q.colnames)%in%c(3,4) | !length(txdb.colnames)%in%c(3,4)){
    stop("Column names must be 3 or 4 charatrs long")
  }
  q.colnames.full <- q.colnames
  txdb.colnames.full <- txdb.colnames
  get.strand.colnames.f <- function(dt, col.ids, check=check.orientation){
    if(check!=TRUE){
      strand.col <- NA
    } else if(length(col.ids)==4){
      strand.col <- col.ids[2]
    } else{
      strand.col.all <- grep("strand", colnames(dt), ignore.case = TRUE, value = TRUE)
      if(length(strand.col.all)==1){
        strand.col <- strand.col.all
      } else if(length(strand.col.all)==0){
        strand.col <- NA
      } else{
        strand.col <- strand.col.all[which.min(nchar(strand.col.all))][1] 
        warning(paste0("multiple strand columns detected in query. Using ", strand.col, "."))
      }
    }
    return(strand.col)
  }
  
  if(ignore.strand==FALSE & (length(q.colnames)==3 | length(txdb.colnames)==3)){
    warning("Strand column names arent provided. Ignoring strand")
    ignore.strand=TRUE
  }
  if(ignore.strand==TRUE){
    if(length(q.colnames==4)){
      q.colnames <- q.colnames[-2]  
    }
    if(length(txdb.colnames)==4){
      txdb.colnames <- txdb.colnames[-2]  
    }
  }
  
  if(!identical(key(queryRegions), q.colnames)){
    setkeyv(queryRegions, q.colnames)
  }
  is.dt=0
  if(is.list(txdbFeatures) & !is.data.table(txdbFeatures)){
    for(i in seq_along(txdbFeatures)){
      dt.i <- txdbFeatures[[i]]
      if(!is.data.table(dt.i)){
        txdbFeatures[[i]] <- as.data.table(dt.i)
      } else if(!identical(key(dt.i), txdb.colnames)){
        setkeyv(txdbFeatures[[i]], txdb.colnames)
      }
    }
  } else if(is.data.table(txdbFeatures) & !identical(key(txdbFeatures), txdb.colnames)){
    setkeyv(txdbFeatures, txdb.colnames)
  } else{
    txdbFeatures <- as.data.table(txdbFeatures)
    setkeyv(txdbFeatures, txdb.colnames)
  }
  strand.col.q <- get.strand.colnames.f(queryRegions, q.colnames.full)
  if(is.na(strand.col.q) & check.orientation==TRUE){
    warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
    check.orientation <- FALSE
  }
  
  strand.col.txdb.all <- sapply(txdbFeatures, function(tx) get.strand.colnames.f(tx, txdb.colnames.full, check=check.orientation))
  if(any(is.na(strand.col.txdb.all)) | is.null(strand.col.txdb.all)){
    warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
    strand.col.txdb.all <- NA
    check.orientation <- FALSE
  }
  q.dt <- copy(queryRegions)
  q.dt[, uid:=.I]
  
  summarize <- function(i, by.cols, strand.q=strand.col.q, strand.tx.list=strand.col.txdb.all, check.orient=check.orientation){
    x <- txdbFeatures[[i]]
    if(check.orient==TRUE){
      if(length(strand.tx.list)==1){
        strand.tx <- strand.tx.list
      } else if(length(strand.tx.list)==0){
        strand.tx <- NA
      } else{
        strand.tx <- unlist(strand.tx.list[i])
      }
    } else{
      strand.tx <- NA
    }
    
    if((is.null(strand.tx) | is.na(strand.tx)) & check.orient==TRUE){
      warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
      check.orient <- FALSE
    }
    myOverlaps <- data.table::foverlaps(x, q.dt)[!is.na(uid)]
    if("same.orient"%in%by.cols){
      if(check.orient==TRUE){
        if(strand.tx==strand.q){
          strand.tx <- paste0("i.", strand.tx)
        }
        myOverlaps[, same.orient:=get(strand.col.q)==get(strand.tx)]  
      } else{
        myOverlaps[, same.orient:="NOT.CHECKED"]
      }
    }
    if(is.null(val.cols)){
      ol.dt <- myOverlaps[, .(tot=.N), keyby=by.cols]  
      
    } else{
      ol.dt <- myOverlaps[, lapply(.SD, function(y) y[1]), keyby=by.cols, .SDcols=val.cols]
    }
    ol.dt[, feat:=names(txdbFeatures)[i]]
    return(ol.dt)
  }
  
  if(is.null(by.col)){
    by.col1 <- "uid"
    by.col2 <- "feat"
  } else{
    by.col1 <- c("uid", by.col)
    by.col2 <- c("feat", by.col)
  }
  
  if(check.orientation==TRUE){
    by.col1 <- c(by.col1, "same.orient")
    by.col2 <- c(by.col2, "same.orient")
  }
  
  results <- rbindlist(lapply(seq_along(txdbFeatures), FUN = function(x) summarize(i=x, by.cols=by.col1, strand.q=strand.col.q, strand.tx.list=strand.col.txdb.all, check.orient=check.orientation)))
  if("same.orient"%in%colnames(results)){
    by.col3 <- c(by.col, "same.orient")
  } else{
    by.col3 <- by.col
  }
  if(is.null(val.cols)){
    any.feats <- results[, .N, by=by.col1][, .(tot=.N), by=by.col3]
    
    if("same.orient"%in%colnames(results) & ignore.strand==TRUE){
      any.feats.both <- results[, .N, by=by.col1][, .(tot=.N), by=by.col]
      any.feats.both[, same.orient:=NA]
      any.feats <- rbind(any.feats, any.feats.both)
    }
    any.feats[, feat:="AnyFeatures"]
    by.cols.tot <- intersect(colnames(q.dt), by.col)
    tots <- q.dt[, .(tot=.N), by=by.cols.tot]
    tots[, feat:="total"]
    results.sum <- rbind(results[, .(tot=.N), by=by.col2], any.feats, tots, fill=TRUE)
  } else{
    any.feats.tmp <- results[, lapply(.SD, function(x) x[1]), by=by.col1, .SDcols=val.cols]
    any.feats <- any.feats.tmp[, lapply(.SD, sum), by=by.col3, .SDcols=val.cols]
    if("same.orient"%in%colnames(results) & ignore.strand==TRUE){
      any.feats.both <- any.feats.tmp[, lapply(.SD, sum), by=by.col, .SDcols=val.cols]
      any.feats.both[, same.orient:=NA]
      any.feats <- rbind(any.feats, any.feats.both)
    }
    any.feats[, feat:="AnyFeatures"]
    by.cols.tot <- intersect(colnames(q.dt), by.col)
    tots <- q.dt[, lapply(.SD, sum), by=by.cols.tot, .SDcols=val.cols]
    tots[, feat:="total"]
    results.sum <- rbind(results[, lapply(.SD, sum), by=by.col2, .SDcols=val.cols], any.feats, tots, fill=TRUE)
  }
  
  return(results.sum)
}

# Same as above, but dont summarize counts
summarizeQueryRegions.dt.full <- function(queryRegions, txdbFeatures, q.colnames=c("seqnames", "strand", "start", "end"), txdb.colnames=c("seqnames", "strand", "start", "end"), check.orientation=FALSE, ignore.strand=FALSE, val.cols=NULL, by.col=NULL){
  if(!length(q.colnames)%in%c(3,4) | !length(txdb.colnames)%in%c(3,4)){
    stop("Column names must be 3 or 4 charatrs long")
  }
  q.colnames.full <- q.colnames
  txdb.colnames.full <- txdb.colnames
  get.strand.colnames.f <- function(dt, col.ids, check=check.orientation){
    if(check!=TRUE){
      strand.col <- NA
    } else if(length(col.ids)==4){
      strand.col <- col.ids[2]
    } else{
      strand.col.all <- grep("strand", colnames(dt), ignore.case = TRUE, value = TRUE)
      if(length(strand.col.all)==1){
        strand.col <- strand.col.all
      } else if(length(strand.col.all)==0){
        strand.col <- NA
      } else{
        strand.col <- strand.col.all[which.min(nchar(strand.col.all))][1] 
        warning(paste0("multiple strand columns detected in query. Using ", strand.col, "."))
      }
    }
    return(strand.col)
  }
  
  if(ignore.strand==FALSE & (length(q.colnames)==3 | length(txdb.colnames)==3)){
    warning("Strand column names arent provided. Ignoring strand")
    ignore.strand=TRUE
  }
  if(ignore.strand==TRUE){
    if(length(q.colnames==4)){
      q.colnames <- q.colnames[-2]  
    }
    if(length(txdb.colnames)==4){
      txdb.colnames <- txdb.colnames[-2]  
    }
  }
  
  if(!identical(key(queryRegions), q.colnames)){
    setkeyv(queryRegions, q.colnames)
  }
  is.dt=0
  if(is.list(txdbFeatures) & !is.data.table(txdbFeatures)){
    for(i in seq_along(txdbFeatures)){
      dt.i <- txdbFeatures[[i]]
      if(!is.data.table(dt.i)){
        txdbFeatures[[i]] <- as.data.table(dt.i)
      } else if(!identical(key(dt.i), txdb.colnames)){
        setkeyv(txdbFeatures[[i]], txdb.colnames)
      }
    }
  } else if(is.data.table(txdbFeatures) & !identical(key(txdbFeatures), txdb.colnames)){
    setkeyv(txdbFeatures, txdb.colnames)
  } else{
    txdbFeatures <- as.data.table(txdbFeatures)
    setkeyv(txdbFeatures, txdb.colnames)
  }
  strand.col.q <- get.strand.colnames.f(queryRegions, q.colnames.full)
  if(is.na(strand.col.q) & check.orientation==TRUE){
    warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
    check.orientation <- FALSE
  }
  
  strand.col.txdb.all <- sapply(txdbFeatures, function(tx) get.strand.colnames.f(tx, txdb.colnames.full, check=check.orientation))
  if(any(is.na(strand.col.txdb.all)) | is.null(strand.col.txdb.all)){
    warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
    strand.col.txdb.all <- NA
    check.orientation <- FALSE
  }
  q.dt <- copy(queryRegions)
  q.dt[, uid:=.I]
  
  summarize <- function(i, by.cols, strand.q=strand.col.q, strand.tx.list=strand.col.txdb.all, check.orient=check.orientation){
    x <- txdbFeatures[[i]]
    if(check.orient==TRUE){
      if(length(strand.tx.list)==1){
        strand.tx <- strand.tx.list
      } else if(length(strand.tx.list)==0){
        strand.tx <- NA
      } else{
        strand.tx <- unlist(strand.tx.list[i])
      }
    } else{
      strand.tx <- NA
    }
    
    if((is.null(strand.tx) | is.na(strand.tx)) & check.orient==TRUE){
      warning(paste0("Check orientation = TRUE but no strand info given. Setting check orientation to FALSE"))
      check.orient <- FALSE
    }
    myOverlaps <- data.table::foverlaps(x, q.dt)[!is.na(uid)]
    return(myOverlaps)
  }
  
  if(is.null(by.col)){
    by.col1 <- "uid"
    by.col2 <- "feat"
  } else{
    by.col1 <- c("uid", by.col)
    by.col2 <- c("feat", by.col)
  }
  
  if(check.orientation==TRUE){
    by.col1 <- c(by.col1, "same.orient")
    by.col2 <- c(by.col2, "same.orient")
  }
  
  results <- lapply(seq_along(txdbFeatures), FUN = function(x) summarize(i=x, by.cols=by.col1, strand.q=strand.col.q, strand.tx.list=strand.col.txdb.all, check.orient=check.orientation))
  names(results) <- names(txdbFeatures)
  return(results)
}

# Mutual information stats (from https://tm4ss.github.io/docs/Tutorial_5_Co-occurrence.html) ----
########## MI: log(k*kij / (ki * kj) ########
# m = co-occurence matrix
# k = total reads/documents from all observations
# i = row index
# j = col.index
calcMI.stats.f <- function(i, m, k){
  kj <- diag(m)
  ki <- diag(m)[i]
  kij <- m[i,]
  id.i <- row.names(m)[i]
  nj <- length(kj)
  mutualInformationSig <- log(k * kij / (ki * kj))
  #mutualInformationSig <- mutualInformationSig[order(mutualInformationSig, decreasing = TRUE)]
  
  ########## DICE: 2 X&Y / X + Y ##############
  dicesig <- 2 * kij / (ki + kj)
  #dicesig <- dicesig[order(dicesig, decreasing=TRUE)]
  
  ########## Log Likelihood ###################
  logsig <- 2 * ((k * log(k)) - (ki * log(ki)) - (kj * log(kj)) + (kij * log(kij)) 
                 + (k - ki - kj + kij) * log(k - ki - kj + kij) 
                 + (ki - kij) * log(ki - kij) + (kj - kij) * log(kj - kij) 
                 - (k - ki) * log(k - ki) - (k - kj) * log(k - kj))
  #logsig <- logsig[order(logsig, decreasing=T)]
  df <- data.table(Sample.ID.i=rep(id.i, nj), Sample.ID.j=names(kij), ki=rep(ki, nj), kj=kj, k=rep(k, nj), n.cooc=kij, mutualInformationSig=mutualInformationSig, dicesig=dicesig, logsig=logsig)
  
  return(df)
}


# This script needs to be run first. The script takes the collapsed, processed read files (in fasta format) from all the input files, and further collapses them across samples
# The resulting fasta file contains a header in the format: ><original sequence>-<number of samples where detected>-<total counts across samples>
#aln.script.file <- "./data/Evaluate_Filtering_Steps_HealthyControls/biofragmenta-v1.5_20180220_103636/02b-star_alignment_all/1_collapse_from_collapsed_and_map.sh"
#sys.call(aln.script.file)

# import processed fastq file from original run to get raw read counts
contam.dir <- paste0(base.dir, "/02-map_contaminates")
processed.fasta.files <- dir(contam.dir, pattern="*Processed.fa", full.names=TRUE, recursive = TRUE)

read.count.fa.dt <- rbindlist(lapply(processed.fasta.files, 
                                     FUN=function(fl){
                                       bn <- gsub("_l001_l004_l007_Processed.fa", "", basename(fl))
                                       dt <- fread(fl, sep="\t", header=FALSE)
                                       dt[, grp:=c(0,1)]
                                       setkey(dt, grp)
                                       dt2 <- data.table(read.ID=dt[grp==0]$V1, origSequence=dt[grp==1]$V1)
                                       dt2[, File.Base.ID2:=bn]
                                       return(dt2)
                                     }))
read.count.fa.dt[, sample.read.count:=tstrsplit(read.ID, split="-", fixed=TRUE, type.convert = TRUE, keep=2)]
read.count.fa.dt[, read.ID:=NULL]


# Bam file from unfiltered read alignments (using collapsed FASTA file)
unfilt.bam.file <- grep("02b-star_alignment_all", dir(base.dir, pattern="All_Files_Processed_Aligned.out.sorted.bam$", recursive=TRUE, full.names=TRUE), value = TRUE)
unfilt.bam.index.file <- paste0(unfilt.bam.file, ".bai")
if(length(unfilt.bam.file)==0 | length(unfilt.bam.index.file)==0){
  stop("BAM file with index must exist")
}


# Import bam alignments from full input fasta ----
# Remove secondary alignments and convert to data table
reads.dt1 <- bam.to.dt.f(unfilt.bam.file)
reads.dt1[, c("origSequence", "nSamples.detected", "tot.count"):=tstrsplit(qname, split="-", fixed=TRUE)]
reads.dt1.copy <- copy(reads.dt1)
uniq.read.ids.unfilt <- reads.dt1[, .N, by=.(qname, origSequence, nSamples.detected, tot.count)]


# Import bam alignments after first contaminate filtering ----
step.name <- "04-split_readgroups"
bam.files <- unlist(bamfiles.list[grep(step.name, names(bamfiles.list))], use.names = FALSE)
reads.dt2.list <- lapply(bam.files, 
                         FUN=function(fl){
                           dirnm <- dirname(fl)
                           bn <- gsub("_rRNA-filtered_aligned_genome.*", "", basename(fl))
                           ga <- readGAlignments(fl, param = param)
                           mcols(ga)$File.Base.ID2 <- bn
                           mcols(ga)$dirnm <- basename(dirnm)
                           return(ga)
                         })
reads.dt2.gr <- do.call(c, reads.dt2.list)
reads.dt2 <- as.data.table(reads.dt2.gr)
setnames(reads.dt2, "qname", "read.ID")


# Import bam alignments after second contaminate filtering ----
step.name <- "05-filtered"
bam.files <- unlist(bamfiles.list[grep(step.name, names(bamfiles.list))], use.names = FALSE)
reads.dt3.list <- lapply(bam.files, 
                         FUN=function(fl){
                           dirnm <- dirname(fl)
                           bn <- gsub("_rRNA-filtered_aligned_genome.*", "", basename(fl))
                           ga <- readGAlignments(fl, param = param)
                           mcols(ga)$File.Base.ID2 <- bn
                           mcols(ga)$dirnm <- basename(dirnm)
                           return(ga)
                         })
reads.dt3.gr <- do.call(c, reads.dt3.list)
reads.dt3 <- as.data.table(reads.dt3.gr)
setnames(reads.dt3, "qname", "read.ID")

reads.dt23 <- rbind(reads.dt2, reads.dt3)
rm(reads.dt2)
rm(reads.dt3)
rm(list = c("reads.dt2.gr", "reads.dt3.gr"))
setnames(reads.dt23, "seq", "origSequence")
setkey(reads.dt23, origSequence, dirnm, File.Base.ID2)

uniq.read.count <- reads.dt23[, .N, by=.(origSequence, dirnm, File.Base.ID2, read.ID)][, .(uniq.read.count=.N), by=.(origSequence, dirnm, File.Base.ID2)]
read.detected.in.filt.set <- dcast.data.table(uniq.read.count, origSequence~dirnm, fill=0, value.var = "uniq.read.count", fun.aggregate = sum)

setkey(read.detected.in.filt.set, origSequence)
setnames(read.detected.in.filt.set, make.names(colnames(read.detected.in.filt.set)))
setkey(reads.dt1, origSequence)

sel.cols <- colnames(read.detected.in.filt.set)[!colnames(read.detected.in.filt.set)%in%key(read.detected.in.filt.set)]
new.cols <- paste0("pass.contam.filt", seq_along(sel.cols))
reads.dt1[read.detected.in.filt.set, (sel.cols):=lapply(sel.cols, function(x) get(paste0("i.", x)))][, (new.cols):=lapply(.SD, function(y) y>0 & !is.na(y)), .SDcols=sel.cols][, (sel.cols):=NULL][]

setkey(read.count.fa.dt, File.Base.ID2)
old.filenames <- read.count.fa.dt[, .N, by=File.Base.ID2]$File.Base.ID2

THIS.SAMPLE.INFO <- ALL.SAMPLE.INFO[, old.filename:=old.filenames[pmatch(File.Base.ID2, old.filenames)]][!is.na(old.filename)]
setkey(THIS.SAMPLE.INFO, old.filename)

read.count.fa.dt[THIS.SAMPLE.INFO, File.Base.ID2:=i.File.Base.ID2, on=c(File.Base.ID2="old.filename")]


# Annotate Alignments ----
# subsets the combined data table by the individual file ids.
# assumes it's mostly a bed file
# Annotation function ----

annotate.sample.peakfiles <- function(peakfile.dt = peak.files.all,
                                      file.id = NULL,
                                      chr.col = 1,
                                      start.col = 2,
                                      end.col = 3,
                                      TxDb = txdb.gencode,
                                      txdb.name = "txdb.gencode",
                                      chr.convert = chr.gencode2ucsc,
                                      annot.priority = c("Exon",
                                                         "Intron",
                                                         "5UTR",
                                                         "3UTR",
                                                         "Promoter",
                                                         "Downstream",
                                                         "Intergenic"),
                                      annoDb = NULL) {
  anno <- annoDb
  if (is.null(file.id)) {
    dt <- peakfile.dt
  } else{
    dt <- peakfile.dt[file.id]
  }
  col.names <- colnames(dt)
  gr.cols <-
    c(
      col.names[chr.col],
      col.names[start.col],
      col.names[end.col],
      grep(
        "strand",
        colnames(peakfile.dt),
        ignore.case = TRUE,
        value = TRUE
      )[1]
    )
  mcol.names <- col.names[!col.names %in% gr.cols]
  gr <-
    dt[, GRanges(seqnames = get(gr.cols[1]),
                 ranges = IRanges(start = get(gr.cols[2]),  end = get(gr.cols[3])))]
  if (length(gr.cols) == 4) {
    strand(gr) <- dt[, Rle(get(gr.cols[4]))]
  }
  if (length(mcol.names) > 0) {
    mcols(gr) <-
      data.frame(subset(dt, select = mcol.names), stringsAsFactors = FALSE)
  }
  if (!all(seqlevels(gr) %in% seqlevels(TxDb))) {
    warning(paste0("Not all seqlevels found in txdb"), call. = TRUE)
    if (txdb.name == "TxDb.Hsapiens.UCSC.hg38.knownGene" &
        all(seqlevels(gr) %in% seqlevels(txdb.gencode))) {
      warning("Renaming Seqlevels to UCSC", call. = TRUE)
      gr2 <-
        renameSeqlevels(gr, chr.gencode2ucsc[seqlevels(gr), ]$UCSC)
      assign(x = "gr",
             value = gr2,
             pos = -1)
    } else if (txdb.name == "txdb.gencode" &
               all(seqlevels(gr) %in% seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene))) {
      warning("Renaming Seqlevels to Gencode", call. = TRUE)
      gr2 <-
        renameSeqlevels(gr, chr.ucsc2gencode[seqlevels(gr), ]$gencode)
      assign(x = "gr",
             value = gr2,
             pos = -1)
    } else{
      warning("Dropping levels to include only those in provided txdb")
      gr2 <- keepSeqlevels(gr, seqlevels(TxDb))
      assign(x = "gr",
             value = gr2,
             pos = -1)
    }
  }
  
  annot.out <-
    annotatePeak(
      gr,
      TxDb = TxDb,
      genomicAnnotationPriority = annot.priority,
      annoDb = annoDb,
      level = "transcript"
    )
  return(annot.out)
  return(gr)
}

tx.db.file <- "./data/TxDb.Hsapiens.Gencode27.GRCh38.primaryAssembly"
txdb.gencode <- loadDb(tx.db.file)
# For some reason, one read with a "N" gapped alignment showed up in the ouptut, even though it was explicitly disallowed by the STAR alignment parameters.
# No other file across the ~400 Biofragmenta files had this. Remove it here
read.count.fa.dt <- subset(read.count.fa.dt, origSequence!="AATCCACTGCCTCCAAAGAAACGATTGAACAGGAGAAGCAGCAGGCGAATCGT")
reads.dt1 <- subset(reads.dt1, origSequence!="AATCCACTGCCTCCAAAGAAACGATTGAACAGGAGAAGCAGCAGGCGAATCGT")


# Save big files ----
reads.dt1.fl <- "20180810_All_Unfilt_read_collapsed_align.txt"
read.count.fa.dt.fl <- "20180810_All_Unfilt_read_collapsed_counts_per_sample.txt"
#fwrite(reads.dt1, file = reads.dt1.fl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#fwrite(read.count.fa.dt, file = read.count.fa.dt.fl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# Start here if already processed ----
reads.dt1 <- fread(reads.dt1.fl, header=TRUE, sep="\t")
read.count.fa.dt <- fread(read.count.fa.dt.fl, header=TRUE, sep="\t")

if(!"ALL.SAMPLE.INFO"%in%ls()){
  load("ALL.SAMPLE.INFO")
}
file.ids.all <- read.count.fa.dt[, .N, by=File.Base.ID2]$File.Base.ID2
keep.cols <- c("File.Base.ID2", "Description", "participant.ID", "replicate.n", "PNK")

THIS.SAMPLE.INFO <- subset(ALL.SAMPLE.INFO, File.Base.ID2%in%file.ids.all, select=keep.cols)
THIS.SAMPLE.INFO[, Sample.ID:=paste0(participant.ID, ".", PNK)]
THIS.SAMPLE.INFO[, Sample.Repl.ID:=paste0(participant.ID, ".", PNK, ".", replicate.n)]
setkey(read.count.fa.dt, File.Base.ID2)
setkey(THIS.SAMPLE.INFO, File.Base.ID2)
read.count.fa.dt <- merge(read.count.fa.dt, THIS.SAMPLE.INFO)


reads.dt1[, thru.biofrag.stage:=1+rowSums(.SD), .SDcols=c("pass.contam.filt1", "pass.contam.filt2")]
reads.dt1.gr <- reads.dt1[, GRanges(seqnames=seqnames, strand=strand, ranges=IRanges(start, end))]
gr.reduce <- GenomicRanges::reduce(reads.dt1.gr, with.revmap=FALSE)
gr.reduce.dt <- data.table(as.data.frame(gr.reduce))
gr.reduce.dt[, width:=NULL]
gr.reduce.dt[, overlappingQuery:=paste(seqnames, start, end, strand, sep=":")]
mcols(gr.reduce) <- gr.reduce.dt$overlappingQuery

setkey(gr.reduce.dt, seqnames, strand, start, end)
setkey(reads.dt1, seqnames, strand, start, end)
reads.dt1 <- foverlaps(gr.reduce.dt, reads.dt1, verbose = TRUE)


#genes.gr.anno.red <- annotate.sample.peakfiles(peakfile.dt=gr.reduce.dt[seqnames%in%seqlevels(txdb.gencode)], txdb.name="txdb.gencode", TxDb=txdb.gencode)
#genes.gr.anno.dt <- as.data.table(genes.gr.anno.red)
#genes.gr.anno.dt.detail <- as.data.table(genes.gr.anno.red@detailGenomicAnnotation)
#genes.gr.anno.dt <- data.table(genes.gr.anno.dt, genes.gr.anno.dt.detail)
#genes.gr.anno.dt[, `:=`(seqnames=as.character(seqnames), strand=as.character(strand))]

#setkeyv(genes.gr.anno.dt, c("seqnames", "strand", "start", "end"))
#setkeyv(reads.dt1, c("seqnames", "strand", "start", "end"))
#reads.dt1.annot <- foverlaps(reads.dt1, genes.gr.anno.dt, type="within", verbose=TRUE)
#setkeyv(reads.dt1.annot, c("annotation", "annotation"))
#reads.dt1.annot[, genomic.feature:=tstrsplit(annotation[1], split=" (", fixed=TRUE, keep = 1), by="annotation"]

# Get gene/trx id for overlapping exon feature
#reads.dt1.annot[genomic.feature=="Downstream", genomic.feature:=annotation]
#setkey(reads.dt1.annot, annotation, genomic.feature)
#reads.dt1.annot[c("Exon", "Intron"), c("ol.txid", "ol.geneid"):=tstrsplit(annotation[1], split="\\(|/|,", keep=c(2,3)), by=annotation, on="genomic.feature"]
#reads.dt1.annot[!c("Exon", "Intron"), `:=`(ol.txid=transcriptId, ol.geneid=geneId)]

# Annotation V2 ----

#gtf.loc="/data/genomic_data/Homo_sapiens/UCSC/GRCh38/ANNOTATIONS/GENCODE/v27/gencode.v27.primary_assembly.annotation.gtf"
gtf.loc="/data/genomic_data/Homo_sapiens/UCSC/GRCh38/ANNOTATIONS/GENCODE/v27/gencode.v27.basic.annotation.gtf"
gtf.tmp <- paste0("/dev/shm/", basename(gtf.loc))
if(!file.exists(gtf.tmp)){
  system(paste0("cp ", gtf.loc, " ", gtf.tmp))
}
gtf <- importGtf(gtf.tmp)
txdbFeatures <- getTxdbFeaturesFromGRanges(gtf)
txdbFeatures.dt <- lapply(txdbFeatures, FUN=function(x){dt <- as.data.table(x); setkey(dt, seqnames, strand, start, end)})


queryRegions <- gr.reduce
queryRegions.opp <- queryRegions
strand(queryRegions.opp) <- ifelse(strand(queryRegions)=="+", 
                                   "-",
                                   ifelse(strand(queryRegions)=="-",
                                          "+",
                                          "*"))
gtf.dt <- as.data.table(gtf)
gene.types.dt <- gtf.dt[, .N, by=gene_type]
keep.features <- c("protein_coding", "processed_transcript", "lincRNA", "3prime_overlapping_ncRNA", "antisense", "non_coding", 
                   "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncRNA", 
                   "macro_lncRNA")
mRNA.lincRNA.txids <- gtf.dt[gene_type%in%keep.features & transcript_type%in%keep.features, .N, by=transcript_id]
gtf.mRNA.lincRNA <- gtf[gtf$transcript_id%in%mRNA.lincRNA.txids$transcript_id]
mRNA.lincRNA.txdbFeatures <- lapply(txdbFeatures, FUN=function(x) x[x$tx_name%in%gtf.mRNA.lincRNA$transcript_id])
mRNA.lincRNA.txdbFeatures.dt <- lapply(mRNA.lincRNA.txdbFeatures, FUN=function(x){dt <- as.data.table(x); setkey(dt, seqnames, start, end)})

transcript.type.dt <- gtf.dt[!is.na(transcript_id), .N, by=.(transcript_id, transcript_type)]
setkey(transcript.type.dt, transcript_id)
invisible(lapply(txdbFeatures.dt, FUN=function(x) x[transcript.type.dt, transcript_type:=i.transcript_type, on=c(tx_name="transcript_id")]))
invisible(lapply(seq_along(txdbFeatures.dt), FUN=function(x){ txdbFeatures.dt[[x]][, feature_type:=names(txdbFeatures.dt)[x], on=c(tx_name="transcript_id")]}))

keep.features <- c("exons", "introns", "promoters")
txdbFeatures.dt.c <- rbindlist(txdbFeatures.dt, fill=TRUE)
txdbFeatures.dt.filt <- txdbFeatures.dt.c[feature_type%in%keep.features]
txdbFeatures.dt.filt[, c("tx_id", "exon_id", "exon_name", "exon_rank", "cds_id", "cds_name"):=NULL]



setkeyv(gtf.dt, c("seqnames", "start", "end"))
setkeyv(gr.reduce.dt, c("seqnames", "start", "end"))
overlaps2 <- foverlaps(gtf.dt, gr.reduce.dt, verbose=TRUE)
overlaps2.filt <- overlaps2[!is.na(start) & !is.na(i.start)]
overlaps2.filt[, aln.orientation:=ifelse(strand==i.strand, "SENSE", "ANTI")]

#overlaps <- as.data.table(queryGff(queryRegions=queryRegions, gffData = gtf))
overlaps2.filt[, gene_type_orient:=paste0(gene_type, ".", aln.orientation)]


# Now get the merged ID in to the read count info per file
setkey(reads.dt1, origSequence)
setkey(read.count.fa.dt, origSequence, File.Base.ID2)


read.count.fa.dt[, tot.count.sample:=sum(sample.read.count), by=File.Base.ID2]
read.count.fa.dt[, tot.uniq.count.sample:=.N, by=File.Base.ID2]
sample.replicate.count.totals <- read.count.fa.dt[, .(tot.count.sample.repl=tot.count.sample[1], tot.uniq.count.sample.repl=tot.uniq.count.sample[1]), by=.(File.Base.ID2, Sample.Repl.ID, Sample.ID, PNK)]
sample.replicate.count.totals[, `:=`(tot.count.sample=sum(tot.count.sample.repl)), by=Sample.ID]

# Merge replicate data
read.count.seq.by.sample <- read.count.fa.dt[, .(sample.read.count=sum(sample.read.count)), by=.(origSequence, participant.ID, PNK, Sample.ID)]
read.count.seq.by.sample[, `:=`(tot.count.sample=sum(sample.read.count), tot.uniq.count.sample=.N), by=Sample.ID]

sample.count.totals <- read.count.seq.by.sample[, .(tot.count.sample=tot.count.sample[1], tot.uniq.count.sample=tot.uniq.count.sample[1]), by=.(Sample.ID, participant.ID, PNK)]
read.count.seq.by.sample[, cpm.read:=(10^6*sample.read.count)/tot.count.sample]
setkey(read.count.seq.by.sample, origSequence)
# add sample count info to reads.dt1
reads.dt1.sample.counts <- merge(reads.dt1, read.count.seq.by.sample, by="origSequence", allow.cartesian=TRUE)
reads.dt1.sample.counts[, is.unique.map:=NH==1]
setkey(reads.dt1.sample.counts, overlappingQuery, Sample.ID)


reads.dt1.sample.counts.thru2 <- reads.dt1.sample.counts[thru.biofrag.stage!=1]
reads.dt1.sample.counts.thru3 <- reads.dt1.sample.counts.thru2[thru.biofrag.stage!=2]

OQ.summary.f <- function(dt.x, keep.cols){
  
  old.cols <- colnames(dt.x)
  if(!identical(key(dt.x), c("overlappingQuery", "Sample.ID"))){
    setkey(dt.x, overlappingQuery, Sample.ID)
  }
  dt <- copy(dt.x)
  dt[, `:=`(sample.read.count.OQ=sum(sample.read.count),
            sample.cpm.read.OQ=sum(cpm.read),
            sample.num.uniq.reads.OQ=.N),
     by=.(overlappingQuery, Sample.ID)]
  dt[NH==1, `:=`(uniqmap.sample.read.count.OQ=sum(sample.read.count),
                 uniqmap.sample.cpm.read.OQ=sum(cpm.read),
                 uniqmap.sample.num.uniq.reads.OQ=.N),
     by=.(overlappingQuery, Sample.ID)]
  added.cols <- setdiff(colnames(dt), old.cols)
  umap.cols <- grep("uniqmap", added.cols, value=TRUE)
  dt[is.na(uniqmap.sample.num.uniq.reads.OQ), (umap.cols):=0]
  keep.cols <- unique(c(keep.cols, c("overlappingQuery", "Sample.ID"), added.cols))
  print(keep.cols)
  dt.filt <- subset(unique(dt, by=c("overlappingQuery", "Sample.ID")), select=keep.cols)
  
  return(dt.filt)
}

keep.cols.v <- c("seqnames", "strand", "start", "end", "overlappingQuery", "participant.ID",  "PNK", "Sample.ID")

reads.dt1.sample.counts.thru1.sum <- OQ.summary.f(reads.dt1.sample.counts, keep.cols.v)
reads.dt1.sample.counts.thru2.sum <- OQ.summary.f(reads.dt1.sample.counts.thru2, keep.cols.v)
reads.dt1.sample.counts.thru3.sum <- OQ.summary.f(reads.dt1.sample.counts.thru3, keep.cols.v)



# Save summary counts files
filt.stage.summary.coord <- c("reads.dt1.sample.counts.thru1.sum", "reads.dt1.sample.counts.thru2.sum", "reads.dt1.sample.counts.thru3.sum")
import.files=TRUE # set to false, if you need to save them
fread.frwrite.choice.f(file.paths = paste0(filt.stage.summary.coord, ".txt"), r.obj.names = filt.stage.summary.coord, read.file = import.files)

# Get some summary plots from annotatoins
setkey(overlaps2.filt, overlappingQuery)

value.columns <- Hmisc::Cs(sample.read.count.OQ, sample.cpm.read.OQ, sample.num.uniq.reads.OQ, uniqmap.sample.read.count.OQ, uniqmap.sample.cpm.read.OQ, uniqmap.sample.num.uniq.reads.OQ)
setkey(reads.dt1.sample.counts.thru1.sum, seqnames, start, end)
setkey(reads.dt1.sample.counts.thru2.sum, seqnames, start, end)
setkey(reads.dt1.sample.counts.thru3.sum, seqnames, start, end)
thru3.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru3.sum, txdbFeatures = mRNA.lincRNA.txdbFeatures.dt, by.col="Sample.ID", val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)
thru2.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru2.sum, txdbFeatures = mRNA.lincRNA.txdbFeatures.dt, by.col="Sample.ID", val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)
thru1.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru1.sum, txdbFeatures = mRNA.lincRNA.txdbFeatures.dt, by.col="Sample.ID", val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)
thru3.queryRegionsSum[, thru.biofrag.stage:=3]
thru2.queryRegionsSum[, thru.biofrag.stage:=2]
thru1.queryRegionsSum[, thru.biofrag.stage:=1]

queryRegionsSum <- rbind(thru1.queryRegionsSum, thru2.queryRegionsSum, thru3.queryRegionsSum)
queryRegionsSum.melt <- melt(queryRegionsSum, measure.vars = value.columns, variable.name = "measurement", value.name = "value")

queryRegionsSum.melt.tot <- queryRegionsSum.melt[feat=="total"]
queryRegionsSum.melt[queryRegionsSum.melt.tot, total:=i.value, on=c("Sample.ID", "thru.biofrag.stage", "measurement")]
queryRegionsSum.melt[, percent.total:=value/total]
queryRegionsSum.melt[, c("Participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]
queryRegionsSum.melt.f <- queryRegionsSum.melt[!is.na(same.orient) & feat!="total"]

feat.means <- queryRegionsSum.melt.f[measurement=="sample.cpm.read.OQ" & same.orient==TRUE, mean(value), by=feat]
setorder(feat.means, -V1)
queryRegionsSum.melt.f[, feat.f:=factor(feat, levels=feat.means$feat)]

queryRegionsSum.melt.f[, thru.biofrag.stage:=as.factor(thru.biofrag.stage)]
g <- ggplot2::ggplot(queryRegionsSum.melt.f, aes(x = feat.f, y = percent.total)) + 
  geom_boxplot(aes(fill = same.orient), pos="dodge") +  
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_grid(paste0("Step", thru.biofrag.stage)~PNK)

g %+% queryRegionsSum.melt.f[measurement=="sample.read.count.OQ"]


# Get protein-coding exon/link-RNA reads and then see the corrsponding ovrlap types ----
mRNA.linRNA.exon.features.dt <- copy(mRNA.lincRNA.txdbFeatures.dt[["exons"]])
all.exon.features.dt <- copy(txdbFeatures.dt[["exons"]])



setkey(reads.dt1.sample.counts, seqnames, strand, start, end)
setkey(mRNA.linRNA.exon.features.dt, seqnames, strand, start, end)
setnames(reads.dt1.sample.counts, c("i.start", "i.end"), c("olQuery.start", "olQuery.end"))
mRNA.linRNA.exon.ol.thru1 <- foverlaps(reads.dt1.sample.counts, mRNA.linRNA.exon.features.dt, nomatch = 0, verbose = TRUE)
setnames(mRNA.linRNA.exon.ol.thru1, c("i.start", "i.end"), c("read.start", "read.end"))



# Get the reads that overlapped mRNA/lincRNA exons, and pull their original alignments


setkey(mRNA.linRNA.exon.ol.thru1, "origSequence")
mRNA.linRNA.exon.ol.thru1.simple <- subset(mRNA.linRNA.exon.ol.thru1, select=c("origSequence", "tx_name", "gene_name", "exon_id"))
setkey(mRNA.linRNA.exon.ol.thru1.simple, "origSequence")
reads.dt1.sample.counts.all.aligns.from.exon.aligned <- reads.dt1.sample.counts[mRNA.linRNA.exon.ol.thru1[, .N, by=origSequence]$origSequence, on="origSequence"]
setkey(reads.dt1.sample.counts.all.aligns.from.exon.aligned, seqnames, strand, start, end)
all.alns.mRNA.linRNA.exon.ol.thru1 <- foverlaps(reads.dt1.sample.counts.all.aligns.from.exon.aligned, mRNA.linRNA.exon.features.dt, verbose = TRUE)
setnames(all.alns.mRNA.linRNA.exon.ol.thru1, c("i.start", "i.end"), c("read.start", "read.end"))

repeat.file <- "/data/genomic_data/biofragmenta_reference/Homo_sapiens/UCSC/GRCh38/depletion/hg38_miscRNA_and_rmsk.bed"

repeat.dt <- fread(repeat.file, header=FALSE, sep="\t")
# Convert to 1-based coords
setnames(repeat.dt, c("seqnames", "start", "end", "rep.name", "score", "strand"))
repeat.dt[, start:=start+1]
setkeyv(repeat.dt, c("seqnames", "start", "end"))


setkey(all.alns.mRNA.linRNA.exon.ol.thru1, seqnames, read.start, read.end)
all.alns.mRNA.linRNA.exon.ol.thru1.repeats <- foverlaps(all.alns.mRNA.linRNA.exon.ol.thru1, repeat.dt)
setkey(all.alns.mRNA.linRNA.exon.ol.thru1.repeats, rep.name)

all.alns.mRNA.linRNA.exon.ol.thru1.repeats[!is.na(rep.name), split.rep.num:=nchar(rep.name[1])-nchar(gsub("|", "", rep.name[1], fixed=TRUE))+1L, by=rep.name]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats[!is.na(rep.name), rep.short:=tstrsplit(rep.name[1], split="|", fixed=TRUE, keep=split.rep.num[1]), by=rep.name]
setkey(all.alns.mRNA.linRNA.exon.ol.thru1.repeats, rep.short)
all.alns.mRNA.linRNA.exon.ol.thru1.repeats[, has.repeat.group:=ifelse(is.na(rep.short), "NONE", ifelse(split.rep.num==3, "RMSK", "ENS"))]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats[!is.na(rep.name), split.rep.num:=nchar(rep.short[1])-nchar(gsub(":", "", rep.short[1], fixed=TRUE))+1L, by=rep.short]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats[!is.na(rep.short), c("repID", "repFamily", "repClass"):=tstrsplit(rep.short[1], split=":", fixed=TRUE, fill = TRUE), by=rep.name]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats[, has.repeat.group:=ifelse(is.na(rep.short), "NONE", ifelse(split.rep.num==3, "RMSK", "ENS"))]


keep.cols <- c(colnames(reads.dt1.sample.counts.all.aligns.from.exon.aligned), c("exon_id", "exon_rank", "tx_name", "gene_name", "has.repeat.group"))
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast <- dcast.data.table(subset(all.alns.mRNA.linRNA.exon.ol.thru1.repeats, select=keep.cols), ...~has.repeat.group, fill=0, fun.aggregate = length)
rep.cols <- c("ENS", "NONE", "RMSK")
rep.cols.new <- paste0("mm.", rep.cols)
# Gets count of multi-mapped locations where there might be an ensembl/repeat annotation
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[, (rep.cols.new):=lapply(.SD, sum), by=origSequence, .SDcols=rep.cols]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[, rep.map.category:="NONE"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[mm.ENS>0, rep.map.category:="ENS.TRANS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[mm.RMSK>0, rep.map.category:="RMSK.TRANS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[ENS>0, rep.map.category:="ENS.CIS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[RMSK>0, rep.map.category:="RMSK.CIS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[rep.map.category=="NONE", rep.map.category.simple:="NONE"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[rep.map.category=="ENS.CIS"|rep.map.category=="RMSK.CIS", rep.map.category.simple:="CIS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[rep.map.category=="ENS.TRANS"|rep.map.category=="RMSK.TRANS", rep.map.category.simple:="TRANS"]
all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[, rep.map.category:=factor(rep.map.category, levels=c("RMSK.CIS", "ENS.CIS", "RMSK.TRANS", "ENS.TRANS", "NONE"))]

all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[, simple.scale.frac:=1/NH]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant <- all.alns.mRNA.linRNA.exon.ol.thru1.repeats.cast[!is.na(tx_name), .N, by=.(participant.ID, PNK, Sample.ID, rep.map.category.simple, thru.biofrag.stage, tx_name, gene_name, is.unique.map, qname, cpm.read, sample.read.count, simple.scale.frac)][,.(sample.read.count=sum(sample.read.count)), by=.(participant.ID, PNK, Sample.ID, rep.map.category.simple, thru.biofrag.stage, tx_name, gene_name, is.unique.map)]


by.cols <- setdiff(colnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant), c("sample.read.count", "is.unique.map"))
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u <- all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant[, .(sample.read.count.allmap=sum(sample.read.count), sample.read.count.uniq=sum(ifelse(is.unique.map==TRUE, sample.read.count, 0L))), by=by.cols]

# make "TRANS" be CIS+TRANS
val.vars <-  c("sample.read.count.allmap", "sample.read.count.uniq")
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.trans <- all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u[rep.map.category.simple!="NONE", lapply(.SD, sum), by=setdiff(by.cols, "rep.map.category.simple"), .SDcols=val.vars]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.trans[, rep.map.category.simple:="CISorTRANS"]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.all <- all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u[, lapply(.SD, sum), by=setdiff(by.cols, "rep.map.category.simple"), .SDcols=val.vars]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.all[, rep.map.category.simple:="ALL"]


all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2 <- rbind(rbind(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u[rep.map.category.simple!="TRANS"], all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.trans), all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u.all)
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2[, rep.map.category.simple:=factor(rep.map.category.simple, levels=c("ALL", "CISorTRANS", "CIS", "NONE"))]

count.tots.thru1 <- reads.dt1.sample.counts[thru.biofrag.stage==1, .(tot.count.sample=tot.count.sample[1], tot.uniq.count.sample=tot.uniq.count.sample[1]), by=.(Sample.ID, thru.biofrag.stage)]
count.tots.thru2 <- reads.dt1.sample.counts[thru.biofrag.stage>1, .N, by=.(origSequence, Sample.ID, sample.read.count, thru.biofrag.stage)][, .(thru.biofrag.stage=min(thru.biofrag.stage), tot.count.sample=sum(sample.read.count), tot.uniq.count.sample=.N), by=.(Sample.ID)]
count.tots.thru3 <- reads.dt1.sample.counts[thru.biofrag.stage>2, .N, by=.(origSequence, Sample.ID, sample.read.count, thru.biofrag.stage)][, .(thru.biofrag.stage=min(thru.biofrag.stage), tot.count.sample=sum(sample.read.count), tot.uniq.count.sample=.N), by=.(Sample.ID)]
count.tots <- rbindlist(list(count.tots.thru1, count.tots.thru2, count.tots.thru3))

stage.split.gene.quant <- lapply(1:3, FUN=function(x){
  dt <- copy(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2[thru.biofrag.stage>=x])
  dt[, thru.biofrag.stage:=NULL]
  this.sample.tot <- count.tots[thru.biofrag.stage==x]
  dt[this.sample.tot, `:=`(tot.count.sample=i.tot.count.sample, tot.uniq.count.sample=i.tot.uniq.count.sample), on="Sample.ID"]
  dt.cast <- dcast.data.table(dt, ...~rep.map.category.simple, drop = c(TRUE, FALSE), value.var = val.vars, fill=0, fun.aggregate = sum)
  dt.cast[, percent.rep.all:=sample.read.count.allmap_CISorTRANS/sample.read.count.allmap_ALL]
  dt.cast[, percent.rep.cis:=sample.read.count.allmap_CIS/sample.read.count.allmap_ALL]
  dt.cast[, CPM.ALL:=(10^6)*sample.read.count.allmap_ALL/tot.count.sample]
  dt.cast.summary <- dt.cast[, .(thru.stage=x,
                                 median.cpm.ALL=median(CPM.ALL),
                                 wt.mean.cpm.ALL=weighted.mean(CPM.ALL, w = sample.read.count.allmap_ALL),
                                 median.percent.rep.all=median(percent.rep.all), 
                                 wt.mean.percent.rep.all=weighted.mean(percent.rep.all, w = sample.read.count.allmap_ALL),
                                 wt.mean.percent.rep.cis=weighted.mean(percent.rep.cis, w = sample.read.count.allmap_ALL),
                                 median.percent.rep.cis=median(percent.rep.cis)), by=.(PNK, tx_name, gene_name)]
  return(dt.cast.summary)
})



stage.split.gene.quant.all <- melt(rbindlist(stage.split.gene.quant), measure.vars = list("CISorTRANS"=c("median.percent.rep.all", "wt.mean.percent.rep.all"), "CIS"=c("median.percent.rep.cis", "wt.mean.percent.rep.cis")))
stage.split.gene.quant.all2 <- melt(rbindlist(stage.split.gene.quant), measure.vars = c("median.percent.rep.all", "wt.mean.percent.rep.all", "median.percent.rep.cis", "wt.mean.percent.rep.cis"))
stage.split.gene.quant.all2[c("median.percent.rep.all", "wt.mean.percent.rep.all"), rep.map.category:="CISorTRANS", on="variable"]
stage.split.gene.quant.all2[c("median.percent.rep.cis", "wt.mean.percent.rep.cis"), rep.map.category:="CIS", on="variable"]
stage.split.gene.quant.all2[c("median.percent.rep.all", "median.percent.rep.cis"), measure:="median", on="variable"]
stage.split.gene.quant.all2[c("wt.mean.percent.rep.all", "wt.mean.percent.rep.cis"), measure:="wt.mean", on="variable"]
stage.split.gene.quant.all <- copy(stage.split.gene.quant.all2)

setnames(stage.split.gene.quant.all, "value", "percent.mapping.repeats")

ggplot(stage.split.gene.quant.all, aes(y=median.cpm.ALL, color=PNK, x=percent.mapping.repeats)) + geom_point(alpha=0.3) + facet_grid(rep.map.category+PNK~paste0("Step",thru.stage)) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1)) + theme_bw() + theme(legend.position = "top")

# Keep only genes with >= 1 count in >= 2 samples
gene.expressed.count <- all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2[, .(sample.read.count.uniq=1*(max(sample.read.count.uniq)>0)), by=.(Sample.ID, PNK, thru.biofrag.stage, tx_name, gene_name)][, .(n.samples.expressed=sum(sample.read.count.uniq)), by=.(PNK, thru.biofrag.stage, tx_name, gene_name)][, .(n.samples.expressed=max(n.samples.expressed)), by=.(PNK, thru.biofrag.stage, gene_name)]
setnames(gene.expressed.count, "thru.biofrag.stage", "thru.stage")


stage.split.gene.quant.all[, max.cpm.gene:=max(median.cpm.ALL), by=.(PNK, thru.stage, gene_name, variable, rep.map.category)]
stage.split.gene.quant.all[, max.cpm.gene.wtmean:=max(wt.mean.cpm.ALL), by=.(PNK, thru.stage, gene_name, variable, rep.map.category)]

stage.split.gene.quant.all.gene <- stage.split.gene.quant.all[median.cpm.ALL==max.cpm.gene & rep.map.category=="CISorTRANS" & measure=="median", .(median.cpm.ALL=max(median.cpm.ALL), percent.mapping.repeats=max(percent.mapping.repeats)), by=.(PNK, gene_name, thru.stage, variable, rep.map.category)]
stage.split.gene.quant.all.gene[gene.expressed.count, n.samples.expressed:=i.n.samples.expressed, on=c("gene_name", "PNK", "thru.stage")]

stage.split.gene.quant.all.gene.mean <- stage.split.gene.quant.all[wt.mean.cpm.ALL==max.cpm.gene.wtmean & rep.map.category=="CISorTRANS" & measure=="wt.mean", .(wt.mean.cpm.ALL=max(wt.mean.cpm.ALL), percent.mapping.repeats=max(percent.mapping.repeats)), by=.(PNK, gene_name, thru.stage, variable, rep.map.category)]
stage.split.gene.quant.all.gene.mean[gene.expressed.count, n.samples.expressed:=i.n.samples.expressed, on=c("gene_name", "PNK", "thru.stage")]

stage.split.gene.quant.all.gene[, perc.repeat.bin:=cut(percent.mapping.repeats, include.lowest = TRUE, breaks = c(0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0), labels = c("0%", "1-25%", "25-50%", "50-75%", "75-99%", "100%"))]

stage.split.gene.quant.all.gene[, decile.expr:=factor(dplyr::ntile(-median.cpm.ALL, 10)), by=.(PNK, thru.stage, variable, rep.map.category)]

# FIG 2C Percent Repeats By Stage median ----
g <- ggplot(stage.split.gene.quant.all.gene[n.samples.expressed>1], aes(y=median.cpm.ALL, color=PNK, x=percent.mapping.repeats)) + 
  geom_point(alpha=0.4, size=0.1) + 
  facet_grid(paste0("Step",thru.stage)~PNK) + 
  scale_y_log10(limits=c(10^-1.85,10^5), breaks=10^seq(-1,5), labels = scales::trans_format("log10", scales::math_format())) + 
  scale_x_continuous(labels=scales::percent) + 
  theme_bw(base_size = 6) + labs(x="% Reads Aligned Repeats / sRNA", y="CPM") +
  theme(legend.position = "none", axis.ticks = element_line(size=0.25), panel.grid = element_blank(), text = element_text(color="black", size=6)) 
out.fl <- "./output/figures/main/FIG2C_PERC_REPEAT_VS_EXPR.pdf"
save_plot(out.fl, plot = g, base_width = 75, base_height = 80, units="mm")

ggplot(stage.split.gene.quant.all.gene[n.samples.expressed>1], aes(y=median.cpm.ALL, fill=PNK, x=perc.repeat.bin)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggbeeswarm::geom_quasirandom(alpha=0.2, size=0.1) + 
  facet_grid(paste0("Step",thru.stage)~PNK) + 
  scale_y_log10(limits=c(10^-2,10^5), breaks=10^seq(-2,5), labels = scales::trans_format("log10", scales::math_format())) + 
  theme_bw(base_size = 6) + 
  theme(legend.position = "top", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black")) 





# FIG 2C Percent Repeats By Stage weighted mean----
ggplot(stage.split.gene.quant.all.gene.mean[n.samples.expressed>1], aes(y=wt.mean.cpm.ALL, color=PNK, x=percent.mapping.repeats)) + 
  geom_point(alpha=0.2) + 
  facet_grid(paste0("Step",thru.stage)~PNK, scale="free_y") + 
  scale_y_log10(limits=c(10^-2,10^5), breaks=10^seq(-2,5), labels = scales::trans_format("log10", scales::math_format())) + 
  scale_x_continuous(labels=scales::percent) + 
  theme_bw() + 
  theme(legend.position = "top", axis.ticks = element_line(size=0.5), panel.grid = element_blank(), text = element_text(color="black", size=6)) 


ggplot(stage.split.gene.quant.all.gene, aes(y=median.cpm.ALL, color=PNK, x=percent.mapping.repeats)) + geom_point(alpha=0.3) + facet_grid(paste0("Step",thru.stage)~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent) + theme_bw() + theme(legend.position = "top")


stage.split.gene.quant.stage1 <- stage.split.gene.quant[[1]]
stage.split.gene.quant.stage1[, percent.rep.all:=sample.read.count.allmap_CISorTRANS/sample.read.count.allmap_ALL]
stage.split.gene.quant.stage1[, percent.rep.cis:=sample.read.count.allmap_CIS/sample.read.count.allmap_ALL]
#ggplot(stage.split.gene.quant.stage1, aes(x=10^6*sample.read.count.allmap_ALL/tot.count.sample, y=percent.rep.all)) + geom_point() + facet_wrap(~participant.ID) + scale_x_log10()
#ggplot(stage.split.gene.quant.stage1, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.cis)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))
#ggplot(stage.split.gene.quant.stage1, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.all)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))

stage.split.gene.quant.stage2 <- stage.split.gene.quant[[2]]
stage.split.gene.quant.stage2[, percent.rep.all:=sample.read.count.allmap_CISorTRANS/sample.read.count.allmap_ALL]
stage.split.gene.quant.stage2[, percent.rep.cis:=sample.read.count.allmap_CIS/sample.read.count.allmap_ALL]
#ggplot(stage.split.gene.quant.stage2, aes(x=10^6*sample.read.count.allmap_ALL/tot.count.sample, y=percent.rep.all)) + geom_point() + facet_wrap(~participant.ID) + scale_x_log10()
ggplot(stage.split.gene.quant.stage2, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.cis)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))
ggplot(stage.split.gene.quant.stage2, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.all)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))

stage.split.gene.quant.stage3 <- stage.split.gene.quant[[3]]
stage.split.gene.quant.stage3[, percent.rep.all:=sample.read.count.allmap_CISorTRANS/sample.read.count.allmap_ALL]
stage.split.gene.quant.stage3[, percent.rep.cis:=sample.read.count.allmap_CIS/sample.read.count.allmap_ALL]
ggplot(stage.split.gene.quant.stage3, aes(x=10^6*sample.read.count.allmap_ALL/tot.count.sample, y=percent.rep.all)) + geom_point() + facet_wrap(~participant.ID) + scale_x_log10()
ggplot(stage.split.gene.quant.stage3, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.cis)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))
ggplot(stage.split.gene.quant.stage3, aes(y=10^6*sample.read.count.allmap_ALL/tot.count.sample, x=percent.rep.all)) + geom_point() + facet_grid(participant.ID~PNK) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1))


# Look at strand specificity for individual features types. Combine the groups together ----


repeat.dt[, is.enst:=startsWith(rep.name, "ENST")]
rmsk.dt <- copy(repeat.dt[is.enst==FALSE])
rmsk.dt[, rep.short:=tstrsplit(rep.name, split="|", fixed=TRUE, keep=4)]
setkey(rmsk.dt, rep.short)
rmsk.dt[, c("repID", "repFamily", "repClass"):=tstrsplit(rep.short[1], split=":", fixed=TRUE, fill = TRUE), by=rep.short]
rmsk.dt[, c("score", "is.enst"):=NULL]
rmsk.dt[, width:=end-start]
setnames(rmsk.dt, c("rep.name", "repID", "repFamily", "repClass"), c("tx_name", "gene_name", "transcript_type", "feature_type"))
rmsk.dt[, feature.type.simple:="TE"]
rmsk.dt[feature_type%in%c("Simple_repeat", "Low_complexity"), feature.type.simple:="simple_RE"]
rmsk.dt[feature_type%in%c("tRNA", "rRNA", "Satellite", "snRNA"), feature.type.simple:=feature_type]
rmsk.dt[feature_type%in%c("srpRNA", "scRNA", "RNA"), feature.type.simple:=gene_name]
rmsk.dt[, annot.source:="RMSK"]
txdbFeatures.dt.filt[, feature.type.simple:=feature_type]
txdbFeatures.dt.filt[, annot.source:="GENCODE"]


txFeatures.rmsk <- rbind(txdbFeatures.dt.filt, subset(rmsk.dt, select=colnames(txdbFeatures.dt.filt)))
setkey(txFeatures.rmsk, seqnames, start, end)
setkey(reads.dt1.sample.counts, seqnames, start, end)
reads.dt1.annot.all <- foverlaps(reads.dt1.sample.counts, txFeatures.rmsk, nomatch = 0, verbose = TRUE)
reads.dt1.annot.all[, aln.orientation:=ifelse(strand==i.strand, "SENSE", "ANTI")]
reads.dt1.annot.all[, feature.type.simple:="TE"]
setkey(reads.dt1.annot.all, feature_type)
reads.dt1.annot.all[c("Simple_repeat", "Low_complexity"), feature.type.simple:="simple_RE"]
reads.dt1.annot.all[c("tRNA", "rRNA", "Satellite", "snRNA"), feature.type.simple:=feature_type]
reads.dt1.annot.all[c("srpRNA", "scRNA", "RNA"), feature.type.simple:=gene_name]
reads.dt1.annot.all[txdbFeatures.dt.filt[,.N, by=feature_type]$feature_type, feature.type.simple:=feature_type]

reads.dt1.annot.seqCount <- reads.dt1.annot.all[ , .N, by=.(origSequence, Sample.ID, PNK, sample.read.count)]
reads.dt1.annot.all.cast <- dcast.data.table(reads.dt1.annot.all, origSequence+aln.orientation+transcript_type+annot.source~feature.type.simple, fun.aggregate =length , fill=0L, verbose = TRUE)
val.cols <- setdiff(colnames(reads.dt1.annot.all.cast),  key(reads.dt1.annot.all.cast))
rmsk.val.cols <- setdiff(val.cols, txdbFeatures.dt.filt[,.N, by=feature_type]$feature_type)
gencode.val.cols <- txdbFeatures.dt.filt[,.N, by=feature_type]$feature_type
reads.dt1.annot.all.cast[, (val.cols):=lapply(.SD, as.double), .SDcols=val.cols]
setkey(reads.dt1.annot.all.cast, origSequence, aln.orientation, annot.source)
reads.dt1.annot.all.cast[, (rmsk.val.cols):=lapply(.SD, sum), by=.(origSequence, aln.orientation), .SDcols=rmsk.val.cols]
reads.dt1.annot.all.cast[, n.rmsk.feats:=rowSums(.SD>0), .SDcols=rmsk.val.cols]
reads.dt1.annot.all.cast[, n.gencode.feats:=rowSums(.SD>0), .SDcols=gencode.val.cols]
reads.dt1.annot.all.cast[, n.gencode.feat.types:=sum(n.gencode.feats>0), by=.(origSequence, aln.orientation)]


reads.dt1.annot.all.uniqmap <- copy(reads.dt1.annot.all[is.unique.map==TRUE])
reads.dt1.annot.all.cast.uniqmap <- reads.dt1.annot.all.cast[reads.dt1.annot.all.uniqmap[, .N, by=origSequence]$origSequence, on="origSequence"]
feat.cols <- colnames(reads.dt1.annot.all.cast.uniqmap)[which(reads.dt1.annot.all.cast.uniqmap[, unlist(lapply(.SD, is.numeric))])]

reads.dt1.annot.all.cast.uniqmap[, `:=`(n.rmsk.feats.bothOrient=sum(n.rmsk.feats),
                                        n.gencode.feats.bothOrient=sum(n.gencode.feats),
                                        n.gencode.feat.types.bothOrient=sum(n.gencode.feat.types)), keyby=origSequence]



reads.dt1.annot.all.cast.uniqmap[, n.feature.types.all:=apply(.SD, 1, max, na.rm=TRUE), .SDcols=c("n.rmsk.feats", "n.gencode.feats", "n.gencode.feat.types")]
reads.dt1.annot.all.cast.uniqmap[, n.feature.types.all:=sum(n.feature.types.all), by=origSequence]
reads.dt1.annot.all.cast.uniqmap.ufeat <- melt(reads.dt1.annot.all.cast.uniqmap[n.feature.types.all==1], measure.vars = val.cols, id.vars = c("origSequence", "aln.orientation", "transcript_type", "annot.source"))[value!=0]

reads.dt1.annot.seqCount.uniqmap <- reads.dt1.annot.seqCount[reads.dt1.annot.all.cast.uniqmap.ufeat[, .N, by=origSequence]$origSequence, on="origSequence"]
reads.dt1.annot.all.cast.uniqmap.ufeat.sample.info <- merge(reads.dt1.annot.seqCount.uniqmap, reads.dt1.annot.all.cast.uniqmap.ufeat, on="origSequence")


reads.dt1.annot.all.cast.uniqmap.ufeat.count <- reads.dt1.annot.all.cast.uniqmap.ufeat.sample.info[, .(sample.read.count.tot=sum(sample.read.count)), by=.(Sample.ID, PNK, transcript_type, annot.source, aln.orientation, variable)]
setnames(reads.dt1.annot.all.cast.uniqmap.ufeat.count, "variable", "feature.type.simple")
reads.dt1.annot.all.cast.uniqmap.ufeat.count[, feat.tot.sample:=sum(sample.read.count.tot), by=.(Sample.ID, PNK, feature.type.simple)]
reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast <- dcast.data.table(reads.dt1.annot.all.cast.uniqmap.ufeat.count, ...~aln.orientation, value.var="sample.read.count.tot", fun.aggregate = sum)
reads.dt1.annot.all.cast.uniqmap.ufeat.count[, percent.tot:=sample.read.count.tot/feat.tot.sample]

reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast[, percent.sense:=SENSE/(SENSE+ANTI)]
reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast[, feat.tot.sample:=SENSE+ANTI]
ggplot(reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast[feat.tot.sample>=10,], aes(x=feature.type.simple, y=percent.sense, color=log2(feat.tot.sample))) + geom_boxplot() + geom_jitter() + facet_grid(PNK~.) 

reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast.gencode <- reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast[annot.source=="GENCODE"]
ggplot(reads.dt1.annot.all.cast.uniqmap.ufeat.count.cast.gencode[feat.tot.sample>=10,], aes(x=transcript_type, y=SENSE, color=percent.sense-0.5)) + geom_boxplot(pos="dodge") + geom_jitter(pos=position_jitterdodge()) + coord_flip() + facet_grid(feature.type.simple~PNK) + scale_y_log10()



reads.dt1.annot.all.uniqmap[, transcript.type.orient:=paste0(transcript_type, ".", annot.source)]
reads.dt1.annot.all.uniqmap[, feature.type.simple.orient:=paste0(feature.type.simple, ".", annot.source)]

reads.dt1.annot.all.uniqmap[, same.strand:=strand==strand.1]
reads.dt1.annot.all.uniqmap.summ <- dcast.data.table(reads.dt1.annot.all.uniqmap[annot.source=="RMSK"], feature.type.simple~aln.orientation, fun.aggregate = sum, fill = 0, value.var = "sample.read.count")
ggplot(reads.dt1.annot.all.uniqmap.summ, aes(x=feature.type.simple, y=log2(SENSE/ANTI))) + geom_bar(stat="identity")


# disjoint 
setkey(reads.dt1.annot.all, transcript_type, feature.type.simple, aln.orientation)
gencode.feat.info.dt <- reads.dt1.annot.all[annot.source!="RMSK", .N, by=.(annot.source, transcript_type, feature.type.simple, aln.orientation)]
gencode.feat.info.dt[, transcript.type.orient:=paste(feature.type.simple, transcript_type, aln.orientation, sep="|")]

rmsk.feat.info.dt <- reads.dt1.annot.all[annot.source=="RMSK", .N, by=.(annot.source, feature_type, transcript_type, feature.type.simple, aln.orientation)]
rmsk.feat.info.dt[, transcript_type.simple:=tstrsplit(transcript_type, split="-")[1]]

rmsk.feat.info.orig.tmp <- rmsk.dt[, .N, by=.(transcript_type, feature_type, feature.type.simple)]
rmsk.feat.info.orig <- rbind(rmsk.feat.info.orig.tmp, rmsk.feat.info.orig.tmp)
rmsk.feat.info.orig[, aln.orientation:=c(rep.int("ANTI", nrow(rmsk.feat.info.orig.tmp)), rep.int("SENSE", nrow(rmsk.feat.info.orig.tmp)))]
rmsk.feat.info.orig[, tsplit:=tstrsplit(transcript_type, split="-", keep = 1)][, tsplit:=gsub("?", "", tsplit, fixed=TRUE)][]
rmsk.feat.info.orig[, tsplit.n:=.N, by=tsplit]
rmsk.feat.info.orig[tsplit.n==1 & tsplit!=transcript_type, tsplit:=transcript_type]
rmsk.feat.info.orig[, feature.type.simple.n:=.N, by=feature.type.simple]
rmsk.feat.info.orig[, transcript.type.orient:=paste(feature.type.simple, tsplit, aln.orientation, sep="|")]
rmsk.feat.info.orig[tsplit.n>=feature.type.simple.n, transcript.type.orient:=paste(tsplit, feature.type.simple, aln.orientation, sep="|")]
rmsk.feat.info.orig[, annot.source:="RMSK"]
feat.info.dt <- rbind(rmsk.feat.info.orig, gencode.feat.info.dt, fill=TRUE)
setkey(feat.info.dt, transcript_type, feature.type.simple, aln.orientation)

setkey(reads.dt1.annot.all.uniqmap, transcript_type, feature.type.simple, aln.orientation)
reads.dt1.annot.all.uniqmap[feat.info.dt, transcript.type.orient:=i.transcript.type.orient, on=key(reads.dt1.annot.all.uniqmap)]

reads.dt1.annot.all.uniqmap.summ.all <- reads.dt1.annot.all.uniqmap[,  .N, keyby=.(origSequence, annot.source, transcript.type.orient)]
s.m.read.feat.group.all.uniqmap <- sparseM.from.dt(reads.dt1.annot.all.uniqmap.summ.all, i="origSequence", j="transcript.type.orient", x=1, make.binary = TRUE)
s.m.read.feat.group.all.uniqmap.ig <- to.igraph(s.m.read.feat.group.all.uniqmap)
s.m.read.feat.group.all.uniqmap.cp <- crossprod(s.m.read.feat.group.all.uniqmap)


s.m.read.feat.group.all.uniqmap.ngc <- as(s.m.read.feat.group.all.uniqmap, "ngCMatrix")
s.im.x <- as(t(s.m.read.feat.group.all.uniqmap.ngc), "transactions")
ident.set.id.x <- .Call(arules:::R_pnindex, s.im.x@data, NULL, FALSE)
first.names <- s.m.read.feat.group.all.uniqmap[sapply(seq(max(ident.set.id.x)), match, table=ident.set.id.x),]


transcript.type.orient.m <- crossprod(sparseM.from.dt(reads.dt1.annot.all.uniqmap, i="origSequence", j="transcript.type.orient"))
feature.type.simple.orient.m <- crossprod(sparseM.from.dt(reads.dt1.annot.all.uniqmap, i="origSequence", j="feature.type.simple.orient"))


# Guitar analysis ----

gencode.biotypes.groups <- fread("./data/biotypes_annot.txt", header=FALSE)
setnames(gencode.biotypes.groups, c("gene.or.trx", "group", "biotype"))
tx.db.file <- "./data/TxDb.Hsapiens.Gencode27.GRCh38.basic"
txdb.gencode <- loadDb(tx.db.file)
#gc_txdb <- makeGuitarCoordsFromTxDb(txdb.gencode, minimalComponentLength = 20)

gencode.biotypes.groups <- fread("./data/biotypes_annot.txt", header=FALSE)
setnames(gencode.biotypes.groups, c("gene.or.trx", "group", "biotype"))


exons <- exonsBy(txdb.gencode, by="tx", use.names=TRUE)
count <- countOverlaps(exons, exons)


reads.dt1.sample.counts.uniq.gr.list <- GRangesList(lapply(split(reads.dt1.sample.counts.uniq, by="PNK"), FUN=function(x) x[, GRanges(seqnames, IRanges(start, end), strand = strand)]))

reads.dt1.sample.counts.uniq.gr <- reads.dt1.sample.counts.uniq[, GRanges(seqnames, IRanges(start, end), strand = strand, mcols=DataFrame(.SD))]
reads.dt1.sample.counts.uniq.gr <- GNCList(reads.dt1.sample.counts.uniq.gr)

#GuitarPlot(gfeatures = reads.dt1.sample.counts.uniq.gr, GuitarCoordsFromTxDb = gc_txdb)
tx.info <- gtf.dt[type=="transcript"]
setkey(tx.info, transcript_id)

tx <- transcripts(txdb.gencode, columns=c("tx_id", "tx_name", "gene_id"))


tx.info[, tx_row:=match(transcript_id, tx$tx_name)]
tx$transcript_type <- tx.info[tx$tx_name]$transcript_type
lncRNA.features <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncRNA", "antisense", "antisense_RNA", "non_coding", 
                     "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncRNA", 
                     "macro_lncRNA")

sRNA.features <- unlist(unique(sapply(c("miRNA", "misc_RNA", "snRNA", "snoRNA", "rRNA", "ribozyme", "Mt_tRNA", "Mt_rRNA", "scRNA", "scaRNA", "vaultRNA", "sRNA"), grep, x=tx.info[, .N, by=transcript_type]$transcript_type, ignore.case=TRUE, value=TRUE)))


exons <- exonsBy(txdb.gencode, by = c("tx"), use.names=TRUE)
mRNA.lncRNA.txids <- tx.info[!transcript_type%in%sRNA.features]
sRNA.txids <- tx.info[transcript_type%in%sRNA.features]
mRNA.lncRNA.tx <- tx[mRNA.lncRNA.txids$tx_row]
mRNA.lncRNA.exons <- exons[mRNA.lncRNA.tx$tx_name]
sRNA.exons <-  exons[sRNA.txids$transcript_id]
sRNA.exons.gr <- unlist(sRNA.exons, use.names = TRUE)


invisible(lapply(txdbFeatures.dt, FUN=function(x) x[gencode.biotypes.groups[gene.or.trx=="transcript_type"], transcript.type.group:=i.group, on=c(transcript_type="biotype")]))
keep.feats <- c("lncRNA", "protein_coding")
txdbFeatures.gr.mRNA <- lapply(txdbFeatures.dt, FUN=function(x){ x[transcript.type.group=="protein_coding"] })
names(txdbFeatures.gr.mRNA) <- paste0("mRNA.", names(txdbFeatures.dt))
txdbFeatures.gr.lncRNA <- lapply(txdbFeatures.dt, FUN=function(x){ x[transcript.type.group=="lncRNA"] })
names(txdbFeatures.gr.lncRNA) <- paste0("lncRNA.", names(txdbFeatures.dt))
txdbFeatures.gr.sRNA <- lapply(txdbFeatures.dt, FUN=function(x){ x[transcript.type.group=="sRNA"] })
names(txdbFeatures.gr.sRNA) <- paste0("sRNA.", names(txdbFeatures.dt))
txdbFeatures.gr.types <- c(txdbFeatures.gr.mRNA, txdbFeatures.gr.lncRNA)
txdbFeatures.gr.types <- txdbFeatures.gr.types[elementNROWS(txdbFeatures.gr.types)>0]
txdbFeatures.gr.types.grl <- lapply(seq_along(txdbFeatures.gr.types), 
                                    FUN=function(i){
                                      x <- txdbFeatures.gr.types[[i]]
                                      nm <- names(txdbFeatures.gr.types)[i]
                                      get.cols <- setdiff(colnames(x),  c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))
                                      setkey(x, tx_name, seqnames, start, end, strand)
                                      x[, feat.type.grp:=nm]
                                      x.gr <- x[, GRanges(seqnames, IRanges(start, end), strand=strand)]; 
                                      mcols(x.gr) <- x[, DataFrame(.SD), .SDcols=get.cols]
                                      split(x.gr, f=x.gr$tx_name)
                                      
                                      return(x.gr)
                                    })

txdbFeatures.gr.types.dt <- rbindlist(txdbFeatures.gr.types, fill=TRUE)
txdbFeatures.gr.sRNA.exon.dt <- txdbFeatures.gr.sRNA$sRNA.exons
setkey(txdbFeatures.gr.sRNA.exon.dt, seqnames, start, end)
setkey(txdbFeatures.gr.types.dt, seqnames, start, end)
setkey(rmsk.dt, seqnames, start, end)

txdbFeatures.gr.sRNA.exon.dt[, row.num.tmp:=.I]
txdbFeatures.gr.types.dt[, row.num.tmp:=.I]
rmsk.dt[, row.num.tmp:=.I]

setindex(txdbFeatures.gr.sRNA.exon.dt, row.num.tmp)
setindex(txdbFeatures.gr.types.dt, row.num.tmp)
setindex(rmsk.dt, row.num.tmp)

txdbFeatures.gr.types.grl.sRNAOL <- foverlaps(x = txdbFeatures.gr.types.dt, y = txdbFeatures.gr.sRNA.exon.dt, nomatch=0, which=TRUE)
txdbFeatures.gr.types.dt.sRNAOL <- copy(txdbFeatures.gr.types.dt[txdbFeatures.gr.types.grl.sRNAOL, on=c(row.num.tmp="xid")])
txdbFeatures.gr.types.dt.sRNAOL[txdbFeatures.gr.sRNA.exon.dt, `:=`(start=i.start, end=i.end, transcript_type=i.transcript_type, transcript.type.group=i.transcript.type.group, rep.strand=i.strand), on=c(yid="row.num.tmp")]

txdbFeatures.gr.types.grl.sRNAOL.dt <- txdbFeatures.gr.sRNA.exon.dt[txdbFeatures.gr.types.grl.sRNAOL$yid, .SD][, xid:=txdbFeatures.gr.types.grl.sRNAOL$xid][]
txdbFeatures.gr.types.grl.sRNAOL.dt[, c("feat.type.grp", "host_tx_name"):=txdbFeatures.gr.types.dt[xid, .(feat.type.grp, tx_name)]]

txdbFeatures.gr.types.grl.sRNAOL.dt2 <- txdbFeatures.gr.sRNA.exon.dt[txdbFeatures.gr.types.grl.sRNAOL$yid, .SD][, xid:=txdbFeatures.gr.types.grl.sRNAOL$xid][]


names(txdbFeatures.gr.types.grl) <- names(txdbFeatures.gr.types)

# Drop sRNAs before summarizing ----
setkey(reads.dt1.sample.counts.thru1.sum, seqnames, start, end)
setkey(reads.dt1.sample.counts.thru2.sum, seqnames, start, end)
setkey(reads.dt1.sample.counts.thru3.sum, seqnames, start, end)
setkey(txdbFeatures.gr.sRNA.exon.dt, seqnames, start, end)
gtf.full.loc="/data/genomic_data/Homo_sapiens/UCSC/GRCh38/ANNOTATIONS/GENCODE/v27/gencode.v27.primary_assembly.annotation.gtf"

gtf.full.tmp <- paste0("/dev/shm/", basename(gtf.full.loc))
if(!file.exists(gtf.full.tmp)){
  system(paste0("cp ", gtf.loc, " ", gtf.full.tmp))
}
gtf.full <- importGtf(gtf.full.tmp, keepStandardChr = FALSE, readFromRds = FALSE, saveObjectAsRds = TRUE, overwriteObjectAsRds = TRUE)
gtf.dt.full <- as.data.table(gtf.full)
sRNA.features <- unlist(unique(sapply(c("miRNA", "misc_RNA", "snRNA", "snoRNA", "rRNA", "ribozyme", "Mt_tRNA", "Mt_rRNA", "scRNA", "scaRNA", "vaultRNA", "sRNA"), grep, x=gtf.dt.full[, .N, by=transcript_type]$transcript_type, ignore.case=TRUE, value=TRUE)))

sRNA.exon.dt <- gtf.dt.full[sRNA.features, on="transcript_type"][type=="exon"]
sRNA.rmsk.dt <- rmsk.dt[feature.type.simple!="TE" & feature.type.simple!="simple_RE" & feature.type.simple!="Satellite"]
sRNA.dt <- rbind(sRNA.rmsk.dt, sRNA.exon.dt, fill=TRUE)

setkey(sRNA.dt, seqnames, start, end)
setkey(reads.dt1.sample.counts, seqnames, start, end)

reads.dt1.sample.counts[foverlaps(reads.dt1.sample.counts, sRNA.dt, nomatch = 0, which = TRUE)[, .N, by=xid]$xid, overlaps.sRNA:=TRUE][is.na(overlaps.sRNA), overlaps.sRNA:=FALSE]
reads.dt1.sample.counts[, overlaps.sRNA.indir:=max(overlaps.sRNA), by=qname]

reads.dt1.sample.counts.thru1.sum[foverlaps(reads.dt1.sample.counts.thru1.sum, sRNA.dt, nomatch = 0, which = TRUE)[, .N, by=xid]$xid, overlaps.sRNA:=TRUE][is.na(overlaps.sRNA), overlaps.sRNA:=FALSE]
reads.dt1.sample.counts.thru2.sum[foverlaps(reads.dt1.sample.counts.thru2.sum, sRNA.dt, nomatch = 0, which = TRUE)[, .N, by=xid]$xid, overlaps.sRNA:=TRUE][is.na(overlaps.sRNA), overlaps.sRNA:=FALSE]
reads.dt1.sample.counts.thru3.sum[foverlaps(reads.dt1.sample.counts.thru3.sum, sRNA.dt, nomatch = 0, which = TRUE)[, .N, by=xid]$xid, overlaps.sRNA:=TRUE][is.na(overlaps.sRNA), overlaps.sRNA:=FALSE]

setkey(txdbFeatures.gr.types.dt, seqnames, start, end)
if(txdbFeatures.gr.types.dt[, !all.equal(.I, row.num.tmp)]){
  txdbFeatures.gr.types.dt[, row.num.tmp:=.I]
}

simple.colnames <- c("seqnames", "start", "end", "strand", "tx_name", "gene_name", "feature_type", "transcript_type", "transcript.type.group", "feat.type.grp", "row.num.tmp")
txdbFeatures.gr.types.dt.simple <- subset(txdbFeatures.gr.types.dt, select=simple.colnames)
setindex(txdbFeatures.gr.types.dt.simple, row.num.tmp)



reads.dt1.sample.counts.thru3.sum.featol <- foverlaps(reads.dt1.sample.counts.thru3.sum, txdbFeatures.gr.types.dt.simple, nomatch = 0)
reads.dt1.sample.counts.thru2.sum.featol <- foverlaps(reads.dt1.sample.counts.thru2.sum, txdbFeatures.gr.types.dt.simple, nomatch = 0)
reads.dt1.sample.counts.thru1.sum.featol <- foverlaps(reads.dt1.sample.counts.thru1.sum, txdbFeatures.gr.types.dt.simple, nomatch = 0, verbose=TRUE)
reads.dt1.sample.counts.thru1.sum.featol[, same.orient:=strand==i.strand]
reads.dt1.sample.counts.thru2.sum.featol[, same.orient:=strand==i.strand]
reads.dt1.sample.counts.thru3.sum.featol[, same.orient:=strand==i.strand]
reads.dt1.sample.counts.thru1.sum.featol[, read.query.width:=i.end-i.start+1]
reads.dt1.sample.counts.thru1.sum.featol[, ol.width:=width(pintersect(IRanges(start, end), IRanges(i.start, i.end)))]
reads.dt1.sample.counts.thru1.sum.featol[, perc.ol:=ol.width/read.query.width]
reads.dt1.sample.counts.thru1.sum.featol[, max.perc.ol:=max(perc.ol, na.rm=TRUE), by=.(seqnames, i.start, i.end, i.strand, feat.type.grp, Sample.ID)]

reads.dt1.sample.counts.thru2.sum.featol[, read.query.width:=i.end-i.start+1]
reads.dt1.sample.counts.thru2.sum.featol[, ol.width:=width(pintersect(IRanges(start, end), IRanges(i.start, i.end)))]
reads.dt1.sample.counts.thru2.sum.featol[, perc.ol:=ol.width/read.query.width]
reads.dt1.sample.counts.thru2.sum.featol[, max.perc.ol:=max(perc.ol, na.rm=TRUE), by=.(seqnames, i.start, i.end, i.strand, feat.type.grp, Sample.ID)]

reads.dt1.sample.counts.thru3.sum.featol[, read.query.width:=i.end-i.start+1]
reads.dt1.sample.counts.thru3.sum.featol[, ol.width:=width(pintersect(IRanges(start, end), IRanges(i.start, i.end)))]
reads.dt1.sample.counts.thru3.sum.featol[, perc.ol:=ol.width/read.query.width]
reads.dt1.sample.counts.thru3.sum.featol[, max.perc.ol:=max(perc.ol, na.rm=TRUE), by=.(seqnames, i.start, i.end, i.strand, feat.type.grp, Sample.ID)]

setkey(reads.dt1.sample.counts, seqnames, strand, start, end)
setkey(txdbFeatures.gr.types.dt.simple, seqnames, strand, start, end)
txdbFeatures.gr.types.dt.simple.exons <- txdbFeatures.gr.types.dt.simple[c("lncRNA.exons", "mRNA.exons"), on="feat.type.grp"]
setkey(txdbFeatures.gr.types.dt.simple.exons, seqnames, strand, start, end)

reads.dt1.exon.ol <- foverlaps(reads.dt1.sample.counts, txdbFeatures.gr.types.dt.simple.exons, nomatch = 0, mult = "first", verbose=TRUE)
reads.dt1.exon.ol.nhtally <- reads.dt1.exon.ol[, .N, by=.(origSequence, NH, PNK, participant.ID, thru.biofrag.stage)][, .N, by=.(NH, PNK,participant.ID, thru.biofrag.stage)]
reads.dt1.exon.ol.nhtally[, tot:=sum(N), by=.(thru.biofrag.stage, PNK, participant.ID)]
reads.dt1.exon.ol.nhtally[, percent.tot:=N/tot]
ggplot(reads.dt1.exon.ol.nhtally, aes(x=NH, group=factor(thru.biofrag.stage), color=factor(thru.biofrag.stage), y=percent.tot)) + stat_summary() + geom_line(data=reads.dt1.exon.ol.nhtally[, .(percent.tot=mean(percent.tot)), by=.(NH, PNK, thru.biofrag.stage)]) + scale_x_continuous(limits=c(0,31), breaks=seq(1,31,2)) + theme_bw() + theme(legend.pos="top") + facet_wrap(~PNK, ncol=1)

g1 <- ggplot(reads.dt1.exon.ol.nhtally[NH==1], aes(x=PNK, fill=factor(thru.biofrag.stage), y=percent.tot)) + geom_boxplot(pos="dodge") + scale_y_continuous(labels=scales::percent, limit=c(0,1)) + theme_bw() + labs(y="Uniquely-Mapped Reads (% Total Aligned)", fill="Filtering Stage #")  + theme(legend.pos="top")
g2 <- ggplot(reads.dt1.exon.ol.nhtally[NH==1], aes(x=PNK, fill=factor(thru.biofrag.stage), y=N)) + geom_boxplot(pos="dodge") + scale_y_continuous() + theme_bw() + labs(y="Uniquely-Mapped Reads (% Total Aligned)", fill="Filtering Stage #")  + theme(legend.pos="top")
cowplot::plot_grid(plotlist = list(g1, g2), nrow = 2)
out.fl <- "./output/figures/main/FIG2D_Filtering_Stage_UniquleyMapped_withLegend.pdf"
save_plot(plot = g1, filename = out.fl)
# FIG 2D Filtering_Stage_ PERCENT UniquleyMapped ----
out.fl <- "./output/figures/main/FIG2D_Filtering_Stage_UniquleyMapped_NOLegend.pdf"
g1 <- ggplot(reads.dt1.exon.ol.nhtally[NH==1], aes(x=PNK, fill=factor(thru.biofrag.stage), y=percent.tot)) + 
  geom_boxplot(pos=position_dodge(width=0.75), width=0.5, outlier.colour = NA, alpha=0.7) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.75, size=0.2) + 
  scale_y_continuous(labels=scales::percent, limit=c(0,1)) + theme_bw(base_size = 6) + 
  labs(y="Uniquely-Mapped Exonic Reads (% Total Aligned)", x=NULL, fill="Filtering Stage #")  + 
  theme(legend.pos="none", panel.grid = element_blank(), text=element_text(colour="black"), axis.text=element_text(size=6)); g1
save_plot(plot = g1, filename = out.fl, base_height = 30, base_width = 50, units="mm")

# FIG 2E FINAL Filtering_Stage_ # UniquleyMapped and total----
reads.dt1.exon.ol.nhtally[, uniq.mapped:=ifelse(NH==1, TRUE, FALSE)]
reads.dt1.exon.ol.nhtally.uniq.vs.total <- reads.dt1.exon.ol.nhtally[, .(read.count=sum(N), percent.tot=sum(percent.tot)), by=.(PNK, participant.ID, thru.biofrag.stage, uniq.mapped)]
g1 <- ggplot(reads.dt1.exon.ol.nhtally[NH==1 & thru.biofrag.stage==3], aes(x=PNK, y=N)) + 
  geom_boxplot(pos=position_dodge(width=0.75), width=0.5, outlier.colour = NA, alpha=0.7) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.75, size=0.2) + 
  scale_y_log10(breaks=5) +
  theme_bw(base_size = 6) + 
  labs(y="Uniquely-Mapped Exonic Reads (#)", x=NULL, fill="Filtering Stage #")  +  
  theme(legend.pos="none", panel.grid = element_blank(), text=element_text(colour="black"), axis.text=element_text(size=6)); g1

g1 <- ggplot(reads.dt1.exon.ol.nhtally.uniq.vs.total[thru.biofrag.stage==3], aes(x=PNK, fill=uniq.mapped, y=read.count)) + 
  geom_boxplot(pos=position_dodge(width=0.75), width=0.5, outlier.colour = NA, alpha=0.7) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.75, size=0.2) + 
  scale_y_log10(breaks=5) +
  theme_bw(base_size = 6) + 
  labs(y="Uniquely-Mapped Exonic Reads (#)", x=NULL, fill="Filtering Stage #")  + 
  theme(legend.pos="none", panel.grid = element_blank(), text=element_text(colour="black"), axis.text=element_text(size=6)); g1


out.fl <- "./output/figures/main/FIG2E_Filtering_Stage3_ExonMapped_barplot_by_sample.pdf"
g1 <- ggplot(reads.dt1.exon.ol.nhtally.uniq.vs.total[thru.biofrag.stage==3], aes(x=participant.ID, fill=uniq.mapped, y=read.count)) + facet_wrap(~PNK) +
  geom_bar(pos="stack", stat="identity", size=0.25) + geom_text(aes(label=read.count), size=2, position = position_stack(vjust = 0.5)) +
  #scale_y_log10() +
  theme_bw(base_size = 6) + 
  labs(y="mRNA / lncRNA Exonic Reads (#)", x=NULL, fill="Uniquely-Mapped")  + 
  theme(legend.pos="top", panel.grid = element_blank(), axis.text=element_text(size=6), axis.text.x=element_text(hjust=1, angle=50), text=element_text(colour="black")); g1
save_plot(plot = g1, filename = out.fl, base_height = 75, base_width = 75, units="mm")
out.fl <- "./output/figures/main/FIG2E_Filtering_Stage3_ExonMapped_barplot_by_sample_NOLEGEND.pdf"
g1 <- ggplot(reads.dt1.exon.ol.nhtally.uniq.vs.total[thru.biofrag.stage==3], aes(x=participant.ID, fill=uniq.mapped, y=read.count)) + facet_wrap(~PNK) +
  geom_bar(pos="stack", stat="identity", size=0.25) + geom_text(aes(label=read.count), size=1.5, position = position_stack(vjust = 0.5)) +
  #scale_y_log10() +
  theme_bw(base_size = 6) + 
  labs(y="mRNA / lncRNA Exonic Reads (#)", x=NULL, fill="Uniquely-Mapped")  + 
  theme(legend.pos="none", panel.grid = element_blank(), axis.text.x=element_text(hjust=1, angle=50), text=element_text(colour="black"), axis.text=element_text(size=6)); g1
save_plot(plot = g1, filename = out.fl, base_height = 50, base_width = 75, units="mm")

#out.fl <- "./output/figures/main/FIG2E_Stage3_UniquleyMappedCount_NOLegend.pdf"
#save_plot(plot = g1, filename = out.fl, base_height = 50, base_width = 50, units="mm")

g1 <- ggplot(reads.dt1.exon.ol.nhtally[NH==1 & thru.biofrag.stage==3], aes(x=PNK, y=tot)) + 
  geom_boxplot(pos=position_dodge(width=0.75), width=0.5, outlier.colour = NA, alpha=0.7) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.75, size=0.2) + 
  scale_y_log10(breaks=5) +
  theme_bw(base_size = 6) + 
  labs(y="Total Reads (% Total Aligned)", x=NULL, fill="Filtering Stage #")  + 
  theme(legend.pos="none", panel.grid = element_blank(), text=element_text(colour="black")); g1



setorder(reads.dt1.exon.ol.nhtally, thru.biofrag.stage, NH)
reads.dt1.exon.ol.nhtally.sum <- reads.dt1.exon.ol.nhtally[, .(percent.tot=sum(percent.tot)), keyby=.(PNK, thru.biofrag.stage, NH)]
reads.dt1.exon.ol.nhtally.sum[, cumsum.perc.tot:=cumsum(percent.tot), by=.(PNK, thru.biofrag.stage)]
reads.dt1.exon.ol.nhtally.sum[, cumsum.perc.tot:=cumsum.perc.tot/sum(percent.tot), by=.(PNK, thru.biofrag.stage)]
ggplot(reads.dt1.exon.ol.nhtally.sum, aes(x=NH, group=factor(thru.biofrag.stage), color=factor(thru.biofrag.stage), y=cumsum.perc.tot)) + geom_point() + geom_line() + scale_x_continuous(limits=c(0,30), breaks=seq(1,31,2)) + theme_bw() + theme(legend.pos="top") +  facet_wrap(~PNK, ncol=1)
ggplot(reads.dt1.exon.ol.nhtally, aes(x=NH, group=factor(thru.biofrag.stage), color=factor(thru.biofrag.stage), y=cumsum.n)) + geom_point() + geom_line()  + scale_x_continuous(limits=c(0,30)) + theme_bw()

reads.dt1.sample.counts.thru1.sum.featol.typesum <- dcast.data.table(reads.dt1.sample.counts.thru1.sum.featol, seqnames+i.start+i.end+i.strand+feat.type.grp+Sample.ID~same.orient+feat.type.grp, value.var="perc.ol", fun.aggregate = max, na.rm=TRUE, fill=0)
reads.dt1.sample.counts.thru1.sum.featol.typesum.long <- reads.dt1.sample.counts.thru1.sum.featol[, .(perc.ol=max(perc.ol, na.rm=TRUE), sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames,i.start,i.end,i.strand,feat.type.grp,Sample.ID,same.orient)]

setkey(reads.dt1.sample.counts.thru1.sum.featol.typesum.long, seqnames, i.start, i.end, i.strand, Sample.ID)
reads.dt1.sample.counts.thru1.sum.featol.typesum.long[, grp:=.GRP, by=key(reads.dt1.sample.counts.thru1.sum.featol.typesum.long)]
reads.dt1.sample.counts.thru1.sum.featol.typesum.long[, feat.type.orient:=paste(feat.type.grp[1], same.orient[1], sep="|"), by=.(feat.type.grp, same.orient)]
feat.types <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[, .N, by=feat.type.orient]$feat.type.orient
feat.pair.summary <- rbindlist(lapply(feat.types, FUN=function(feat){
  sum.cols <- c("sample.read.count.OQ", "sample.num.uniq.reads.OQ")
  out.dt <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[reads.dt1.sample.counts.thru1.sum.featol.typesum.long[feat.type.orient==feat], lapply(.SD, sum), on="grp", .SDcols=sum.cols, by=.(feat.type.orient, Sample.ID)]
  out.dt[, feat.type.orient2:=feat]
  out.dt[, c("feat.type.grp", "same.orient"):=tstrsplit(feat.type.orient, split="|", fixed=TRUE)]
  out.dt[, c("feat.type.grp2", "same.orient2"):=tstrsplit(feat.type.orient2, split="|", fixed=TRUE)]
  return(out.dt)
}))

feat.pair.summary2 <- rbindlist(lapply(feat.types, FUN=function(feat){
  sum.cols <- c("sample.read.count.OQ", "sample.num.uniq.reads.OQ")
  #out.dt <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[reads.dt1.sample.counts.thru1.sum.featol.typesum.long[feat.type.orient==feat], lapply(.SD, sum), on="grp", .SDcols=sum.cols, by=.(feat.type.orient, Sample.ID)]
  # out.dt[, feat.type.orient2:=feat]
  feat.groups.dt <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[feat.type.orient==feat]
  group.featcount.dt <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[feat.groups.dt, .N, on="grp", by=.(Sample.ID, grp, feat.type.orient)][, .N, by=.(Sample.ID, grp)]
  out.dt.solo <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[group.featcount.dt[N==1], lapply(.SD, sum), on="grp", .SDcols=sum.cols, by=.(feat.type.orient, Sample.ID)]
  out.dt.solo[, is.unambig.annot.cis:=TRUE]
  group.featcount.dt.mult <- group.featcount.dt[N>1]
  out.dt.mult <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[group.featcount.dt.mult, lapply(.SD, sum), on="grp", .SDcols=sum.cols, by=.(feat.type.orient, Sample.ID)]
  out.dt.mult[, is.unambig.annot.cis:=FALSE]
  out.dt.sum <- rbind(out.dt.solo, out.dt.mult)
  out.dt.sum[, feat.type.orient2:=feat]
  out.dt.sum[, c("feat.type.grp", "same.orient"):=tstrsplit(feat.type.orient, split="|", fixed=TRUE)]
  out.dt.sum[, c("feat.type.grp2", "same.orient2"):=tstrsplit(feat.type.orient2, split="|", fixed=TRUE)]
  return(out.dt.sum)
}))


sum.cols <- c("sample.read.count.OQ", "sample.num.uniq.reads.OQ")


tot.cols <- paste0(sum.cols, ".bothOrient")

feat.pair.summary[, is.dup:=FALSE]
feat.pair.summary[feat.type.grp==feat.type.grp2 & same.orient!=same.orient2, is.dup:=duplicated(feat.type.grp), by=.(Sample.ID)]
dup.val <- feat.pair.summary[is.dup==TRUE, lapply(.SD, max, na.rm=TRUE), by=.(Sample.ID, feat.type.grp), .SD=sum.cols]
feat.pair.summary[feat.type.orient==feat.type.orient2, (tot.cols):=lapply(.SD, sum), by=.(Sample.ID, feat.type.grp), .SDcols=sum.cols]
feat.pair.summary[, (tot.cols):=lapply(.SD, max, na.rm=TRUE), by=.(Sample.ID, feat.type.grp), .SDcols=tot.cols]
feat.pair.summary[dup.val, sample.read.count.OQ.bothOrient:=sample.read.count.OQ.bothOrient-i.sample.read.count.OQ, on=.(Sample.ID, feat.type.grp)]
feat.pair.summary[dup.val, sample.num.uniq.reads.OQ.bothOrient:=sample.num.uniq.reads.OQ.bothOrient-i.sample.num.uniq.reads.OQ, on=.(Sample.ID, feat.type.grp)]

tot.cols <- paste0(sum.cols, ".bothOrient.thisPair")
feat.pair.summary[, (tot.cols):=lapply(.SD, sum), by=.(Sample.ID, feat.type.grp, feat.type.orient2), .SDcols=sum.cols]
feat.pair.summary[, perc.tot.thisPair:=sample.read.count.OQ/sample.read.count.OQ.bothOrient.thisPair]
feat.pair.summary[, perc.tot.bothOrient:=sample.read.count.OQ/sample.read.count.OQ.bothOrient]
g <- ggplot2::ggplot(feat.pair.summary, aes(x = feat.type.orient2, y = perc.tot.bothOrient)) + 
  geom_boxplot(aes(fill = same.orient), pos="dodge") +  
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_wrap(~feat.type.grp)
g %+% feat.pair.summary[feat.type.grp=="mRNA.introns"]


tot.cols <- paste0(sum.cols, ".bothOrient")
thru1.queryRegionsSum2 <- reads.dt1.sample.counts.thru1.sum.featol[, .(sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames, start, end, strand, feat.type.grp, participant.ID, PNK, Sample.ID, same.orient)][, lapply(.SD, sum), by=.(feat.type.grp, participant.ID, PNK, Sample.ID, same.orient), .SDcols=sum.cols]
thru2.queryRegionsSum2 <- reads.dt1.sample.counts.thru2.sum.featol[, .(sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames, start, end, strand, feat.type.grp, participant.ID, PNK, Sample.ID, same.orient)][, lapply(.SD, sum), by=.(feat.type.grp, participant.ID, PNK, Sample.ID, same.orient), .SDcols=sum.cols]
thru3.queryRegionsSum2 <- reads.dt1.sample.counts.thru3.sum.featol[, .(sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames, start, end, strand, feat.type.grp, participant.ID, PNK, Sample.ID, same.orient)][, lapply(.SD, sum), by=.(feat.type.grp, participant.ID, PNK, Sample.ID, same.orient), .SDcols=sum.cols]
thru3.queryRegionsSum2[, thru.biofrag.stage:=3]
thru2.queryRegionsSum2[, thru.biofrag.stage:=2]
thru1.queryRegionsSum2[, thru.biofrag.stage:=1]

queryRegionsSum2 <- rbind(thru1.queryRegionsSum2, thru2.queryRegionsSum2, thru3.queryRegionsSum2)

thru3.queryRegionsSum.lncRNA <- reads.dt1.sample.counts.thru3.sum.featol[transcript.type.group=="lncRNA", .(sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames, start, end, strand, feat.type.grp, transcript_type, participant.ID, PNK, Sample.ID, same.orient)][, lapply(.SD, sum), by=.(feat.type.grp, transcript_type, participant.ID, PNK, Sample.ID, same.orient), .SDcols=sum.cols]
thru3.queryRegionsSum.lncRNA[, thru.biofrag.stage:=3]
thru3.queryRegionsSum.lncRNA[, (tot.cols):=lapply(.SD, sum), by=.(feat.type.grp, transcript_type, participant.ID, PNK, Sample.ID, thru.biofrag.stage), .SDcols=sum.cols]
thru3.queryRegionsSum.lncRNA[, percent.tot.readcount:=sample.read.count.OQ/sample.read.count.OQ.bothOrient]
thru3.queryRegionsSum.lncRNA[, percent.tot.uniqread:=sample.num.uniq.reads.OQ/sample.num.uniq.reads.OQ.bothOrient]
thru3.queryRegionsSum.lncRNA[, nsamples.present:=sum(ifelse(sample.num.uniq.reads.OQ>=10, 1, 0)), by=.(feat.type.grp, transcript_type, participant.ID, thru.biofrag.stage)]
g <- ggplot2::ggplot(thru3.queryRegionsSum.lncRNA, aes(x = sub("lncRNA.", "", feat.type.grp), fill=PNK, y = percent.tot.uniqread)) + 
  geom_boxplot(pos="dodge", width=0.75, outlier.color = NA, aes(fill=PNK)) + geom_point(size=0.1, pos=position_jitterdodge(jitter.height = 0, jitter.width = 0.3, dodge.width = 0.75)) + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_wrap(~transcript_type) + geom_hline(yintercept = 0.5, linetype=2) +
  scale_y_continuous(labels=scales::percent_format(), breaks=seq(0,1,0.1), limits=c(0,1)) +
  labs(y="% Feature reads aligned in sense orientation", x=NULL); g %+% thru3.queryRegionsSum.lncRNA[same.orient==TRUE & nsamples.present>1]

# Co-occurrence calculations ----
sample.ids <- reads.dt1.sample.counts.thru1.sum.featol.typesum.long[, .N, by=Sample.ID]$Sample.ID
feat.ol.thru3.incidence.m.list <- rbindlist(lapply(sample.ids, 
                                                   FUN=function(x){
                                                     m <- sparseM.from.dt(reads.dt1.sample.counts.thru3.sum.featol.typesum.long[Sample.ID==x], i = "grp", j="feat.type.orient", make.binary = TRUE)
                                                     cp <- crossprod(m)
                                                     nr <- nrow(m)
                                                     dt1 <- rbindlist(lapply(seq_along(colnames(m)), calcMI.stats.f, m=cp, k=nr))
                                                     dt1[, Sample.ID:=x]
                                                     return(dt1[Sample.ID.i!=Sample.ID.j])
                                                   }))

feat.ol.thru3.incidence.m.list2 <- lapply(sample.ids, 
                                          FUN=function(x){
                                            sparseM.from.dt(reads.dt1.sample.counts.thru1.sum.featol.typesum.long[Sample.ID==x], i = "grp", j="feat.type.orient", make.binary = TRUE)
                                          })



features.all.l <- data.table(features.l1 = c("lncRNA.promoters", "mRNA.promoters", "lncRNA.exons", "lncRNA.introns", "mRNA.introns", "mRNA.fiveUTRs", "mRNA.cds", "mRNA.threeUTRs"),
                             features.l2 = c("lncRNA.promoters", "mRNA.promoters", "lncRNA.exons", "lncRNA.introns", "mRNA.introns", "mRNA.exons", "mRNA.exons", "mRNA.exons"),
                             features.l3 = c("lncRNA.promoters", "mRNA.promoters", "lncRNA.transcripts", "lncRNA.transcripts", "mRNA.transcripts", "mRNA.transcripts", "mRNA.transcripts", "mRNA.transcripts"))
txdbFeatures.gr.types.dt.simple.l1 <-  txdbFeatures.gr.types.dt.simple[unique(features.all.l$features.l2), on="feat.type.grp"]
setkey(txdbFeatures.gr.types.dt.simple.l1, seqnames, start, end)
setkey(reads.dt1.sample.counts, seqnames, start, end)
reads.dt1.sample.counts.thru1.featol <- foverlaps(reads.dt1.sample.counts, txdbFeatures.gr.types.dt.simple.l1, nomatch = 0, verbose=TRUE)
reads.dt1.sample.counts.thru3.featol <- foverlaps(reads.dt1.sample.counts[thru.biofrag.stage==3], txdbFeatures.gr.types.dt.simple.l1, nomatch = 0, verbose=TRUE)

reads.dt1.sample.counts.thru3.featol[, same.orient:=strand==i.strand]
reads.dt1.sample.counts.thru3.featol[, read.query.width:=i.end-i.start+1]
reads.dt1.sample.counts.thru3.featol[, ol.width:=width(pintersect(IRanges(start, end), IRanges(i.start, i.end)))]
reads.dt1.sample.counts.thru3.featol[, perc.ol:=ol.width/read.query.width]
reads.dt1.sample.counts.thru3.featol[, max.perc.ol:=max(perc.ol, na.rm=TRUE), by=.(seqnames, i.start, i.end, i.strand, feat.type.grp, Sample.ID)]
reads.dt1.sample.counts.thru3.featol.typesum.long <- reads.dt1.sample.counts.thru3.featol[, .(perc.ol=max(perc.ol, na.rm=TRUE), sample.read.count.OQ=sample.read.count[1]), by=.(seqnames,i.start,i.end,i.strand,feat.type.grp,Sample.ID, participant.ID, PNK, same.orient, qname, overlaps.sRNA.indir)]
reads.dt1.sample.counts.thru3.featol.typesum.long[, max.perc.ol.read:=max(perc.ol, na.rm=TRUE), by=.(qname, Sample.ID)]

reads.dt1.sample.counts.thru3.featol.typesum.long[, feat.type.orient:=paste(feat.type.grp[1], same.orient[1], sep="|"), by=.(feat.type.grp, same.orient)]
reads.dt1.sample.counts.thru3.featol.typesum.long[, feat.type.sample.orient:=paste(Sample.ID[1], feat.type.grp[1], same.orient[1], sep="|"), by=.(Sample.ID, feat.type.grp, same.orient)]
test.m.read <- sparseM.from.dt(reads.dt1.sample.counts.thru3.featol.typesum.long[perc.ol==max.perc.ol.read & Sample.ID==sample.ids[1]], i = "qname", j="feat.type.sample.orient", make.binary = TRUE)
test.m.read.all <- sparseM.from.dt(reads.dt1.sample.counts.thru3.featol.typesum.long[perc.ol==max.perc.ol.read], i = "qname", j="feat.type.sample.orient", make.binary = TRUE)

test.mfa.pnk.dt <- reads.dt1.sample.counts.thru3.featol.typesum.long[, .N, by=.(qname, participant.ID, PNK, feat.type.grp, same.orient)]
test.mfa.pnk.dt[, orient.count:=.N, by=.(qname, participant.ID, PNK, feat.type.grp)]
test.mfa.pnk.dt.sum <- dcast.data.table(test.mfa.pnk.dt, qname+participant.ID+PNK~feat.type.grp, value.var="orient.count", fun.aggregate = max, na.rm=TRUE, fill=0)
test.mfa.pnk.dt.sum.df <- data.frame(test.mfa.pnk.dt.sum, row.names=1)
test.mfa.pnk.dt.sum.mfa <- MFA(test.mfa.pnk.dt.sum, ind.sup = 1, type = length(4:9), group=4:9)

test.m.read.cp <- data.table(melt(as.matrix(crossprod(test.m.read))))
test.m.read.all.cp <- data.table(melt(as.matrix(crossprod(test.m.read.all))))

test.m.read.cp.matr <- as.matrix(crossprod(test.m.read))
test.m.read.cp.matr.all <- as.matrix(crossprod(test.m.read.all))
test.m.read.cp.matr.all.SENSE <- test.m.read.cp.matr.all[, grep("TRUE$", colnames(test.m.read.cp.matr.all))]
test.m.read.cp.matr.all.AS <- test.m.read.cp.matr.all[, grep("FALSE$", colnames(test.m.read.cp.matr.all))]


test.m.read.cp[, c("Sample.ID","feat.type.grp1", "same.orient1"):=tstrsplit(Var1, split="|", fixed=TRUE)]
test.m.read.cp[, c("feat.type.grp2", "same.orient2"):=tstrsplit(Var2, split="|", fixed=TRUE, keep=c(2,3))]
test.m.read.cp[, c("Participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]
test.m.pnk.xt <- xtabs(value~feat.type.grp1+same.orient1+feat.type.grp2+same.orient2+Participant.ID, data=test.m.read.cp[PNK=="PNK"])
test.m.pnk.st <- structable(~feat.type.grp1+same.orient1+feat.type.grp2+same.orient2, data = test.m.pnk.xt)

test.m.pnk.st <- structable(~feat.type.grp1+same.orient1+feat.type.grp2+same.orient2, data = test.m.pnk.xt)

test.m.pnk.xt.tot <- xtabs(value~Var2+same.orient1+feat.type.grp1, data=test.m.read.cp[PNK=="PNK"])
test.m.pnk.xt.tot <- xtabs(value~Var2+Var1, data=test.m.read.cp[PNK=="PNK"])

test.m.pnk.xt.tot.dt <- data.frame(dcast.data.table(data.table(as.data.frame(test.m.pnk.xt.tot)), Var2~feat.type.grp1+same.orient1, value.var="Freq", fill=0), row.names=1)
mfa <- MFA(test.m.pnk.xt.tot.dt, group=rep(2, ncol(test.m.pnk.xt.tot.dt)/2), type=rep("f", ncol(test.m.pnk.xt.tot.dt)/2))

plot(mfa, choix="freq", invisible="ind", habillage="group")

dcast.data.table(test.m.read.cp[PNK=="PNK"], feat.type.grp1~same.orient1, value.var="value", fun.aggregate = sum)

null.v <- rep(0, 12)
combn.matr <- do.call(rbind, lapply(1:12, function(x){ do.call(rbind, combinat::combn(x=1:12, m=x, fun=function(y){ x.v <- null.v; x.v[y]<-1; return(x.v)}, simplify=FALSE))}))
combn.matr.melt <- data.table(melt(combn.matr))
combn.matr.melt[, n.feat.types:=sum(value), by=Var1]
feat.type.count.read <- reads.dt1.sample.counts.thru3.featol.typesum.long[overlaps.sRNA.indir==0, .N, by=.(qname, feat.type.orient)][, n.feat.types.read:=.N, by=.(qname)][]
feat.type.lvls <- feat.type.count.read[, .N, keyby=feat.type.orient]$feat.type.orient
feat.type.count.read[, feat.type.orient:=factor(feat.type.orient, levels=feat.type.lvls)]
feat.type.count.read[, feat.type.orient.n:=as.numeric(feat.type.orient)]
setkey(feat.type.count.read, qname, feat.type.orient.n)

feat.type.count.read.tmp <- head(feat.type.count.read[n.feat.types.read>1], n=50)
combn.lvls <- lapply(seq_along(feat.type.lvls), combn, x=seq_along(feat.type.lvls))
feat.type.count.read[n.feat.types.read==1, combn.idx:=as.double(feat.type.orient.n)]
feat.type.count.read[n.feat.types.read!=1, combn.idx:=as.double(which.max(apply(combn.lvls[[n.feat.types.read[1]]], 2, function(x) sum(x==feat.type.orient.n)))), by=qname]
feat.type.count.read.combn.idx <- feat.type.count.read[, .(combn.idx.n=paste0(n.feat.types.read[1], "-", combn.idx[1])), keyby=qname]
setkey(reads.dt1.sample.counts.thru3.featol.typesum.long, qname)
reads.dt1.sample.counts.thru3.featol.typesum.long[feat.type.count.read.combn.idx, combn.idx.n:=i.combn.idx.n, on="qname"]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count <- reads.dt1.sample.counts.thru3.featol.typesum.long[, .N, by=.(qname, Sample.ID, participant.ID, PNK, combn.idx.n)][, .(combn.count=.N), by=.(Sample.ID, participant.ID, PNK, combn.idx.n)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count[, c("n.feat.types.read", "combn.idx"):=tstrsplit(combn.idx.n, split="-", type.convert = TRUE)]
combn.lvls.dt <- data.table(melt(combn.lvls))
setnames(combn.lvls.dt, c("idx.in.combn", "combn.idx", "feat.type.orient.n", "n.feat.types.read"))

reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls <- merge(reads.dt1.sample.counts.thru3.featol.typesum.combn.count, combn.lvls.dt, by=c("combn.idx", "n.feat.types.read"))
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.orient:=feat.type.lvls[feat.type.orient.n]]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, c("feat.type.grp", "same.orient"):=tstrsplit(feat.type.orient, split="|", fixed=TRUE)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, n.orients.in.combn:=.N, by=.(Sample.ID, participant.ID, PNK, n.feat.types.read, combn.idx, feat.type.grp)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.count.in.combn:=sum((combn.count/n.orients.in.combn)), by=.(Sample.ID, participant.ID, PNK, n.feat.types.read, feat.type.grp)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.count.all.combn:=sum(combn.count/n.orients.in.combn), by=.(Sample.ID, participant.ID, PNK, feat.type.grp)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.stranded.count.all.combn:=sum(combn.count), by=.(Sample.ID, participant.ID, PNK, feat.type.grp, same.orient)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.stranded.count.all.combn.pct:=feat.type.stranded.count.all.combn/feat.type.count.all.combn]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.stranded.count.in.combn.pct:=combn.count/feat.type.count.in.combn]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, read.orient:=ifelse(same.orient==TRUE, "S", "AS")]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, c("RNA.family", "feat.struct"):=tstrsplit(feat.type.grp, split=".", fixed=TRUE)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, RNA.family:=factor(RNA.family, levels=c("mRNA", "lncRNA"))]
ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[read.orient=="S" & n.feat.types.read==1], aes(x=toupper(feat.struct), y=feat.type.stranded.count.all.combn.pct)) + geom_boxplot(pos="dodge") + ggbeeswarm::geom_quasirandom() + geom_hline(yintercept = 0.5, linetype=2) + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks=seq(0,1,by=0.1)) + labs(x="Feature", y="% Sense-Aligned", title="Unambiguous Annotations") + facet_grid(RNA.family~PNK) + theme_bw(base_size = 7) +
  theme(legend.position = "top", 
        strip.text=element_text(),
        text = element_text(color="black")) 


# PLOT: Shows % of reads for each feature that align with N other feature types. ----
# The first positoins shows unambiguous annotaitons to the same feature type. 

reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, n.feat.types.read.simple:=factor(ifelse(n.feat.types.read>2, ">2", as.character(n.feat.types.read)), levels=c("1", "2", ">2"))]
g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK", .(combn.count=sum(combn.count)/feat.type.count.all.combn), 
                                                                     by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read, feat.type.count.all.combn)], 
       aes(x=factor(n.feat.types.read),
           y=combn.count, 
           fill=read.orient,
           pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.25) +
  facet_grid(toupper(feat.struct)~RNA.family) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x="# Annotated Features", y="Unique Reads (% Feature Total)", title="PNK (+)") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", 
        strip.text=element_text(),
        text = element_text(color="black"),
        legend.direction = "horizontal") 
ggsave(plot = g, filename = "PNK_Read_Percent_By_Num_Features.png", width = 6, height=5, units = "in")  

# SUPPL FIG 4A Gene S/AS Mapping ----
# Simplify to >2 
g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK", .(combn.count=sum(combn.count)/feat.type.count.all.combn), 
                                                                     by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read.simple, feat.type.count.all.combn)], 
       aes(x=n.feat.types.read.simple,
           y=combn.count, 
           fill=read.orient,
           pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.25) +
  facet_grid(RNA.family~toupper(feat.struct)) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x="# Annotated Features", y="Unique Reads (% Feature Total)") +
  theme_bw(base_size = 6) +
  theme(legend.position = "top", 
        panel.grid=element_blank(),
        text = element_text(color="black"),
        legend.direction = "horizontal") 
ggsave(plot = g , filename = "./output/figures/main/FIG4ASUPPLEMENT_PNK_Read_Percent_By_Num_Features_Simple.png", width = 7, height=5, units = "in")  

# FIG 4A Gene S/AS Mapping ----
g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK" & n.feat.types.read.simple==1, .(combn.count=sum(combn.count)/feat.type.count.all.combn), 
                                                                          by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read.simple, feat.type.count.all.combn)], 
            aes(x=toupper(feat.struct),
                y=combn.count, 
                fill=read.orient,
                pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), width=0.5, size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.2) +
  facet_wrap(~RNA.family) + 
  scale_y_continuous(labels = scales::percent) +
  labs(y="Unique Reads (% Feature Total)", x=NULL) +
  theme_bw(base_size = 6) 
g +  theme(legend.position = "top", 
           panel.grid=element_blank(),
           text = element_text(color="black"),
           legend.direction = "horizontal", axis.text=element_text(size=6)) 
ggsave(plot = g,filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Num_Features_Simple_LEGEND.pdf", width = 75, height=50, units = "mm")  
g +  theme(legend.position = "none", 
           panel.grid=element_blank(),
           text = element_text(color="black"),
           legend.direction = "horizontal", axis.text=element_text(size=6)) 

ggsave(plot = g,filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Num_Features_Simple_NOLEGEND.pdf", width = 80, height=50, units = "mm")  

# 2019-27-01 FIG 4A Gene S/AS Mapping Update -- show as % of uniquely-mapped only ----


reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[, feat.type.count.all.combn.at.numfeat.types:=sum(combn.count/n.orients.in.combn), by=.(Sample.ID, participant.ID, PNK, feat.type.grp, n.feat.types.read.simple)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq <- dcast.data.table(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK" & n.feat.types.read.simple==1], feat.type.grp + RNA.family + feat.struct + Sample.ID + PNK + n.feat.types.read.simple + feat.type.count.all.combn.at.numfeat.types ~ read.orient, value.var = "combn.count", fun.aggregate = sum)

reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq[, perc.S:=S/feat.type.count.all.combn.at.numfeat.types]

reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long <- reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK" & n.feat.types.read.simple==1, .(combn.count=sum(combn.count)), by=.(feat.type.grp , RNA.family , feat.struct , Sample.ID , PNK , n.feat.types.read.simple , feat.type.count.all.combn.at.numfeat.types, read.orient)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long[, tot.uniqannot.sample:=sum(combn.count), by=.(Sample.ID, RNA.family)]


setorder(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long, feat.type.grp, RNA.family, feat.struct, read.orient, Sample.ID)
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long[, perc.combn.count:=combn.count/feat.type.count.all.combn.at.numfeat.types ]
p <- ggboxplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long, x = "feat.struct", y = "perc.combn.count",
               color = "read.orient", palette = "npg",  facet.by = "RNA.family", outlier.color=NULL) + ggbeeswarm::geom_quasirandom(dodge.width = 0.8, size=0.2, aes(color=read.orient)) + scale_y_continuous(labels=scales::percent) + theme(text = element_text(size=6, colour="black")); p

p2 <- p + stat_compare_means(aes(group=read.orient, label=paste0("p = ", ..p.adj..)),  paired = FALSE, method.args = list(alternative = "greater",p.adjust.methods = "fdr"), size=1)
ggsave(plot = p2, filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Unambig_Num_Features_Simple_LEGEND.pdf", width = 75, height=50, units = "mm")  
ggsave(plot = p2 + theme(legend.position = "none"), filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Unambig_Num_Features_Simple_NOLEGEND.pdf", width = 75, height=50, units = "mm")  

p <- ggpaired(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.long, x = "read.orient", y = "perc.combn.count",
              color = "read.orient", palette = "jco", 
              line.color = "gray", line.size = 0.4,
              facet.by = "feat.type.grp", short.panel.labs = TRUE, scale="free_y") + scale_y_continuous(labels=scales::percent)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method.args = list(alternative = "less"), ref.group = "AS", size=6)


feat.table.s.as <- as.matrix(data.frame(dcast.data.table(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK" & n.feat.types.read.simple==1], feat.type.grp ~ read.orient, value.var = "combn.count", fun.aggregate = sum), row.names=1))
myChiSq <- chisq.test(feat.table.s.as)

reads.dt1.sample.counts.thru3.sum.featol.copy <- copy(reads.dt1.sample.counts.thru3.sum.featol)

trx.start.loc <- gtf.dt[type=="transcript", .(TSS=min(ifelse(strand=="+", start, end))), by=.(transcript_id)]
reads.dt1.sample.counts.thru3.sum.featol.copy[trx.start.loc, TSS:=i.TSS, on=c("tx_name"="transcript_id")]
n.feat.types.read <- reads.dt1.sample.counts.thru3.sum.featol.copy[, .N, by=.(feat.type.grp, seqnames, i.start, i.strand, i.end, participant.ID, PNK)][, .(n.feat.types=.N), by=.(seqnames, i.start, i.strand, i.end, participant.ID, PNK)]

reads.dt1.sample.counts.thru3.sum.featol.copy.promoter <- reads.dt1.sample.counts.thru3.sum.featol.copy[feature_type=="promoters"]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[n.feat.types.read, n.feat.types:=i.n.feat.types, on=.(seqnames, i.start, i.strand, i.end, participant.ID, PNK)]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[, abs.dist.to.tss:=abs((i.start+i.end)/2-TSS)]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[, min.abs.dist.to.tss:=min(abs.dist.to.tss), by=.(seqnames, i.start, i.strand, i.end, participant.ID, PNK)]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[, rel.dist.to.tss:=min.abs.dist.to.tss]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[(((i.start+i.end)/2-TSS > 0) & strand=="-") | (((i.start+i.end)/2-TSS < 0) & strand=="+"), rel.dist.to.tss:=min.abs.dist.to.tss*-1]

reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary <- reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[min.abs.dist.to.tss==abs.dist.to.tss & n.feat.types==1, .(sample.read.count.OQ=sum(sample.num.uniq.reads.OQ)), by=.(same.orient, feat.type.grp, Sample.ID, PNK, participant.ID, rel.dist.to.tss) ]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary[, tot:=sum(sample.read.count.OQ), by=.(feat.type.grp, Sample.ID)]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary[, sample.read.perc.OQ:=sample.read.count.OQ/tot, by=.(feat.type.grp, Sample.ID)]
ggplot(reads.dt1.sample.counts.thru3.sum.featol.copy.promoter[PNK=="PNK"], 
       aes(x=participant.ID, y=rel.dist.to.tss, fill=same.orient)) + geom_boxplot(position = "dodge") + facet_wrap(~feat.type.grp, ncol=1)
ggplot(reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary,
       aes(x=rel.dist.to.tss, y=sample.read.perc.OQ, group=same.orient, color=same.orient)) + geom_point() + stat_smooth() + facet_wrap(~feat.type.grp, ncol=1)

reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary[, pos.cut:=cut_width(rel.dist.to.tss, width = 50)]
reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary.ints <- reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary[, .(sample.read.perc.OQ.int=sum(sample.read.perc.OQ)), by=.(same.orient, feat.type.grp, Sample.ID, PNK, participant.ID, pos.cut)]
ggplot(reads.dt1.sample.counts.thru3.sum.featol.copy.promoter.summary.ints[PNK=="PNK"],
       aes(x=pos.cut, y=sample.read.perc.OQ.int, group=same.orient, color=same.orient)) + geom_point() + stat_smooth() + facet_wrap(~feat.type.grp, ncol=1) + theme(axis.text.x = element_text(hjust=1, angle=50))

overall.sense.anti.orient.counts <- reads.dt1.annot.all[, .N, by=.(participant.ID, PNK, aln.orientation, origSequence)][, .N, by=.(participant.ID, PNK, aln.orientation)]
overall.sense.anti.orient.counts[, tot:=sum(N), by=.(participant.ID, PNK)]
overall.sense.anti.orient.counts[, percent.orient:=N/tot, by=.(participant.ID, PNK)]

reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq.list <- lapply(split(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq[, .(AS=sum(AS), S=sum(S)), by=.(RNA.family, feat.struct)], f=reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvlsuniq[, .(AS=sum(AS), S=sum(S)), by=.(RNA.family, feat.struct)]$RNA.family), function(x){ as.matrix(data.frame(x[, c("feat.struct", "S", "AS")], row.names=1))})


g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="PNK" & n.feat.types.read.simple==1, .(combn.count=sum(combn.count)), 
                                                                          by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read.simple, feat.type.count.all.combn.at.numfeat.types)], 
            aes(x=toupper(feat.struct),
                y=combn.count, 
                fill=read.orient,
                pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), width=0.5, size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.2) +
  facet_wrap(~RNA.family) + 
  scale_y_continuous(labels = scales::percent) +
  labs(y="Read Count (% Unambiguously-Annotated)", x=NULL) +
  theme_bw(base_size = 6) 
g +  theme(legend.position = "top", 
           panel.grid=element_blank(),
           text = element_text(color="black"),
           legend.direction = "horizontal", axis.text=element_text(size=6)) 
ggsave(plot = g,filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Unambig_Num_Features_Simple_LEGEND.pdf", width = 75, height=50, units = "mm")  
g +  theme(legend.position = "none", 
           panel.grid=element_blank(),
           text = element_text(color="black"),
           legend.direction = "horizontal", axis.text=element_text(size=6)) 

ggsave(plot = g,filename = "./output/figures/main/FIG4A_PNK_Read_Percent_By_Unambig_Num_Features_Simple_NOLEGEND.pdf", width = 80, height=50, units = "mm")  



g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK!="PNK", .(combn.count=sum(combn.count)), 
                                                                     by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read, feat.type.count.all.combn)], 
       aes(x=factor(n.feat.types.read),
           y=combn.count/feat.type.count.all.combn, 
           fill=read.orient,
           pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.25) +
  facet_grid(toupper(feat.struct)~RNA.family) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x="# Annotated Features", y="Unique Reads (% Feature Total)", title="PNK (-)") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", 
        strip.text=element_text(),
        text = element_text(color="black"),
        legend.direction = "horizontal") 
ggsave(plot = g,filename = "NOPNK_Read_Percent_By_Num_Features.png", width = 6, height=5, units = "in")  


# Simplify to >2
g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[PNK=="NoPNK", .(combn.count=sum(combn.count)/feat.type.count.all.combn), 
                                                                     by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, read.orient, n.feat.types.read.simple, feat.type.count.all.combn)], 
       aes(x=n.feat.types.read.simple,
           y=combn.count, 
           fill=read.orient,
           pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.25) +
  facet_grid(RNA.family~toupper(feat.struct)) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x="# Annotated Features", y="Unique Reads (% Feature Total)", title="PNK (-)") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", 
        strip.text=element_text(),
        text = element_text(color="black"),
        legend.direction = "horizontal") 
ggsave(plot = g,filename = "NOPNK_Read_Percent_By_Num_Features_Simple.png", width = 7, height=5, units = "in")  



# Look at pairwise interactions
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.2feat <- reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls[n.feat.types.read==2]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.2feat[, tot.feat.count.pairs:=sum(combn.count), by=.(feat.type.grp, RNA.family, feat.struct, Sample.ID, PNK, n.feat.types.read)]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs <- reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.2feat[reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.2feat[, .(combn.idx.n, Sample.ID, feat.type.orient, feat.type.grp, read.orient)], on=c("combn.idx.n", "Sample.ID")][feat.type.orient!=i.feat.type.orient]
reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs[, i.feat.type.orient.rename:=ifelse(grepl("TRUE", i.feat.type.orient), 
                                                                                                        gsub("|TRUE", "(S)", i.feat.type.orient, fixed=TRUE),
                                                                                                        gsub("|FALSE", "(AS)", i.feat.type.orient, fixed=TRUE)), by=1:nrow(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs)]

g <- ggplot(reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs, 
            aes(x=i.feat.type.orient.rename,
                y=combn.count/feat.type.count.all.combn, 
                fill=read.orient,
                pos=read.orient)) + 
  geom_boxplot(pos=position_dodge(width=0.9), size=0.5, alpha=0.7, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(dodge.width=0.9, size=0.25) +
  facet_grid(toupper(feat.struct)~RNA.family , scale="free_y") + 
  scale_y_continuous(labels = scales::percent) +
  labs(x="Ambiguously Co-annotated Feature", y="Unique Reads (% Feature Total)") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", 
        axis.text.x=element_text(angle=50, hjust=1),
        strip.text=element_text(),
        text = element_text(color="black"),
        legend.direction = "horizontal"); g + labs(title="PNK (+)") %+% reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs[PNK=="PNK"]
ggsave(plot = g,filename = "PNK_Read_Percent_By_PairedAmbig_Features.png", width = 6, height=5, units = "in")  
g  + labs(title="PNK (-)") %+% reads.dt1.sample.counts.thru3.featol.typesum.combn.count.lvls.pairs[PNK!="PNK"]
ggsave(plot = g,filename = "NoPNK_Read_Percent_By_PairedAmbig_Features.png", width = 6, height=5, units = "in")  

test.m.read.cp.comp <- dcast.data.table(test.m.read.cp, feat.type.grp1+feat.type.grp2+same.orient2~same.orient1, value.var = "value", fun.aggregate = sum, fill=0)
setnames(test.m.read.cp.comp, c("FALSE", "TRUE"), c("FALSE.pair", "TRUE.pair"))
test.m.read.cp.comp[, `:=`(orig.TRUE=sum(ifelse(feat.type.grp1==feat.type.grp2, TRUE.pair, 0)), orig.FALSE=sum(ifelse(feat.type.grp1==feat.type.grp2, FALSE.pair, 0))), by=feat.type.grp1]
test.m.read.cp.comp.pairs <- test.m.read.cp.comp[feat.type.grp1!=feat.type.grp2]
test.m.read.cp.comp.pairs[, `:=`(droppair.TRUE=orig.TRUE-TRUE.pair, droppair.FALSE=orig.FALSE-FALSE.pair)]
test.m.read.cp.comp.pairs[, `:=`(pct.sence.orig=orig.TRUE/(orig.TRUE+orig.FALSE), pct.sense.droppair=droppair.TRUE/(droppair.TRUE+droppair.FALSE))]
ggplot(test.m.read.cp.comp.pairs, aes(x=feat.type.grp1, fill=same.orient2, y=pct.sence.orig-pct.sense.droppair)) + geom_bar(stat="identity", pos="dodge") + facet_wrap(~paste0(feat.type.grp2, "|", same.orient2)) + coord_flip()



test.m.read.all.cp[, c("Sample.ID","feat.type.grp1", "same.orient1"):=tstrsplit(Var1, split="|", fixed=TRUE)]
test.m.read.all.cp[, c("feat.type.grp2", "same.orient2"):=tstrsplit(Var2, split="|", fixed=TRUE, keep=c(2,3))]
test.m.read.all.cp[, c("participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]
test.m.read.all.cp.comp <- dcast.data.table(test.m.read.all.cp, participant.ID + PNK + feat.type.grp1+feat.type.grp2+same.orient2~same.orient1, value.var = "value", fun.aggregate = max, fill=0)
setnames(test.m.read.all.cp.comp, c("FALSE", "TRUE"), c("FALSE.pair", "TRUE.pair"))
test.m.read.all.cp.comp[, `:=`(orig.TRUE=sum(ifelse(feat.type.grp1==feat.type.grp2, TRUE.pair, 0)), orig.FALSE=sum(ifelse(feat.type.grp1==feat.type.grp2, FALSE.pair, 0))), by=.(feat.type.grp1, participant.ID, PNK)]
test.m.read.all.cp.comp[, `:=`(droppair.TRUE=orig.TRUE-TRUE.pair, droppair.FALSE=orig.FALSE-FALSE.pair)]
test.m.read.all.cp.comp[, `:=`(pct.sence.orig=orig.TRUE/(orig.TRUE+orig.FALSE), pct.sense.droppair=droppair.TRUE/(droppair.TRUE+droppair.FALSE))]
type.tots.all <- reads.dt1.sample.counts.thru3.featol.typesum.long[perc.ol==max.perc.ol.read, .N, by=.(qname, feat.type.grp, participant.ID, PNK)][, .(uReads.bothOrient=.N), by=.(feat.type.grp, participant.ID, PNK)]
type.tots.all.cast.orient <- dcast.data.table(reads.dt1.sample.counts.thru3.featol.typesum.long[perc.ol==max.perc.ol.read, .N, by=.(qname, feat.type.grp, participant.ID, PNK, same.orient, feat.type.orient)], feat.type.grp + participant.ID + PNK~same.orient, fill=0, fun.aggregate = length, value.var="qname")
setnames(type.tots.all.cast.orient, c("TRUE", "FALSE"), c("uReads.TRUE", "uReads.FALSE"))
type.tots.all[type.tots.all.cast.orient, `:=`(uReads.FALSE=i.uReads.FALSE, uReads.TRUE=i.uReads.TRUE), on=c("feat.type.grp", "participant.ID", "PNK")]

test.m.read.all.cp.comp.pairs <- test.m.read.all.cp.comp[type.tots.all, on=c("feat.type.grp1"="feat.type.grp", "participant.ID", "PNK")]
test.m.read.all.cp.comp.pairs[, tot.pair:=apply(.SD, 1, sum), .SDcols=c("FALSE.pair", "TRUE.pair")]
test.m.read.all.cp.comp.pairs[, feat.grp.strata:=paste0(feat.type.grp2, "|", same.orient2)]
test.m.read.all.cp.comp.pairs[, `:=`(FALSE.notpair=uReads.FALSE-FALSE.pair, TRUE.notpair=uReads.TRUE-TRUE.pair)]
ex.m <- melt.data.table(test.m.read.all.cp.comp.pairs, id.vars = c("participant.ID", "PNK", "feat.type.grp1", "feat.grp.strata"), measure.vars=c("FALSE.pair", "TRUE.pair", "TRUE.notpair", "FALSE.notpair"))
ex.m[, c("is.sense", "ol.feat.strata"):=tstrsplit(variable, split=".", fixed=TRUE)]
ex.m[, `:=`(read.orient=ifelse(is.sense=="TRUE", "S", "AS"), ol.feat.strata=ifelse(ol.feat.strata=="notpair", FALSE, TRUE))]
ex.m[, n.both.strand:=sum(value), by=.(participant.ID, PNK, feat.type.grp1, feat.grp.strata, ol.feat.strata)]
ex.m[, perc.both:=value/n.both.strand]
ex.m[, grp.matches.strata:=startsWith(feat.grp.strata, feat.type.grp1)]

ex.m.het.pairs <- ex.m[grp.matches.strata==FALSE]
ex.m.self.pairs <- ex.m[grp.matches.strata==TRUE]




ggplot(ex.m[read.orient=="S" & PNK=="PNK" & grp.matches.strata==FALSE & grepl("TRUE", feat.grp.strata)], aes(x=feat.grp.strata, y=perc.both, fill=ol.feat.strata)) + geom_boxplot(pos="dodge") + facet_wrap(~feat.type.grp1) + theme(axis.text=element_text(angle=50, hjust=1)) + geom_hline(yintercept = 0.5,  linetype=2)




reads.dt1.sample.counts.thru3.sum.featol.typesum.long <- reads.dt1.sample.counts.thru3.sum.featol[, .(perc.ol=max(perc.ol, na.rm=TRUE), sample.read.count.OQ=sample.read.count.OQ[1], sample.num.uniq.reads.OQ=sample.num.uniq.reads.OQ[1]), by=.(seqnames,i.start,i.end,i.strand,feat.type.grp,Sample.ID,same.orient)]

setkey(reads.dt1.sample.counts.thru3.sum.featol.typesum.long, seqnames, i.start, i.end, i.strand, Sample.ID)
reads.dt1.sample.counts.thru3.sum.featol.typesum.long[, grp:=.GRP, by=key(reads.dt1.sample.counts.thru3.sum.featol.typesum.long)]
reads.dt1.sample.counts.thru3.sum.featol.typesum.long[, feat.type.orient:=paste(feat.type.grp[1], same.orient[1], sep="|"), by=.(feat.type.grp, same.orient)]

reads.dt1.sample.counts.thru3.sum.featol.typesum.long.l1 <- reads.dt1.sample.counts.thru3.sum.featol.typesum.long[features.all.l$features.l2, on="feat.type.grp"]

type.tots <- dcast.data.table(reads.dt1.sample.counts.thru3.sum.featol.typesum.long.l1, Sample.ID+feat.type.grp+grp~same.orient, fill=0, value.var="sample.num.uniq.reads.OQ", fun.aggregate = sum)
type.tots[, row.minv:=apply(.SD, 1, min), .SDcols=c("FALSE", "TRUE")]
type.tots[, c("tot.OPP", "tot.SAME", "tot.BOTH"):=lapply(.SD, sum), by=.(Sample.ID, feat.type.grp), .SDcols=c("FALSE", "TRUE", "row.minv")]
type.tots.sum <- type.tots[, .(tot.OPP=tot.OPP[1], tot.SAME=tot.SAME[1], tot.BOTH=tot.BOTH[1]), by=.(Sample.ID, feat.type.grp)]
type.tots.sum[, tot.BOTH:=tot.OPP+tot.SAME-tot.BOTH]
type.tots.sum[, p.opp:=tot.OPP/tot.BOTH]
type.tots.sum[, p.same:=tot.SAME/tot.BOTH]
type.tots.sum.pairs <- type.tots.sum[type.tots.sum, on=.(Sample.ID), allow.cartesian=TRUE]
type.tots.sum.pairs[, `:=`(FALSE_FALSE.prob=p.opp*i.p.opp, FALSE_TRUE.prob=p.opp*i.p.same, TRUE_FALSE.prob=p.same*i.p.opp, TRUE_TRUE.prob=p.same*i.p.same)]

type.tots.thru3.pairs.obs <- dcast.data.table(reads.dt1.sample.counts.thru3.sum.featol.typesum.long.l1[reads.dt1.sample.counts.thru3.sum.featol.typesum.long.l1, on=c("Sample.ID", "grp"), allow.cartesian=TRUE], Sample.ID+feat.type.grp+i.feat.type.grp~same.orient+i.same.orient, value.var="sample.num.uniq.reads.OQ", fun.aggregate = sum)
type.tots.thru3.pairs.obs[, tot.PAIR:=apply(.SD, 1, sum), .SDcols=c("FALSE_FALSE", "FALSE_TRUE", "TRUE_FALSE", "TRUE_TRUE")]

type.tots.sum.pairs[type.tots.thru3.pairs.obs, `:=`(tot.PAIR.obs=i.tot.PAIR, FALSE_FALSE.obs=i.FALSE_FALSE, FALSE_TRUE.obs=i.FALSE_TRUE, TRUE_FALSE.obs=i.TRUE_FALSE, TRUE_TRUE.obs=i.TRUE_TRUE) , on=c("Sample.ID", "feat.type.grp", "i.feat.type.grp")]
expect.cols <- paste0(c("FALSE_FALSE", "FALSE_TRUE", "TRUE_FALSE", "TRUE_TRUE"), ".expected")
type.tots.sum.pairs[, (expect.cols):=lapply(.SD, function(x) x * tot.PAIR.obs), .SDcols=c("FALSE_FALSE.prob", "FALSE_TRUE.prob", "TRUE_FALSE.prob", "TRUE_TRUE.prob")]

val.cols <- unlist(lapply(c("FALSE_FALSE", "FALSE_TRUE", "TRUE_FALSE", "TRUE_TRUE"), function(x) grep(x, colnames(type.tots.sum.pairs), value=TRUE)), recursive=TRUE)
type.tots.sum.pairs.melt <- melt(type.tots.sum.pairs, measure.vars = val.cols)
type.tots.sum.pairs.melt[, c("orient.pair", "measure"):=tstrsplit(variable, split=".", fixed=TRUE)]
type.tots.sum.pairs.melt.c <- dcast.data.table(type.tots.sum.pairs.melt[,-c("variable"), with=FALSE], ...~measure, value.var="value", fill=0, fun.aggregate = sum)
type.tots.sum.pairs.melt.c[, enr:=(obs-expected)/sqrt(expected)]
type.tots.sum.pairs.melt.c[, c("Participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]


test.m <- sparseM.from.dt(reads.dt1.sample.counts.thru3.sum.featol.typesum.long.l1[Sample.ID==sample.ids[1]], i = "grp", j="feat.type.orient", x = "sample.num.uniq.reads.OQ" ,make.binary = FALSE)


test.m.cp <- data.table(melt(as.matrix(crossprod(test.m))))
test.m.cp[, c("feat.type.grp1", "same.orient1"):=tstrsplit(Var1, split="|", fixed=TRUE)]
test.m.cp[, c("feat.type.grp2", "same.orient2"):=tstrsplit(Var2, split="|", fixed=TRUE)]
test.m.cp[, feat.pair:=paste0(feat.type.grp1, "|", feat.type.grp2)]
test.m.cp.comp <- dcast.data.table(test.m.cp, feat.type.grp1+feat.type.grp2~same.orient1+same.orient2, value.var = "value", fun.aggregate = sum, fill=0)
test.m.cp.comp.pairs <- test.m.cp.comp[feat.type.grp1!=feat.type.grp2]
test.m.cp.comp.pairs[, tot.PAIR:=apply(.SD, 1, sum), .SDcols=c("FALSE_FALSE", "FALSE_TRUE", "TRUE_FALSE", "TRUE_TRUE")]




test.m.cp.unambig <- data.table(melt(as.matrix(crossprod(test.m[rowSums(test.m>0)==1,]))))[Var1==Var2]



test.m.xt <- test.m.cp[, xtabs(value~same.orient1+same.orient2+feat.type.grp1+feat.type.grp2, data = .SD)]
test.m.xt.list <- lapply(split(test.m.cp, by="feat.pair"), function(x) x[, xtabs(value~same.orient1+same.orient2, data = .SD)])

test.m.xt2 <- test.m.cp[, xtabs(value~feat.type.grp1+same.orient1+Var2, data = .SD)]
test.m.ft <- ftable(same.orient1+same.orient2~feat.type.grp1+feat.type.grp2, test.m.xt)
test.m.ft2 <- ftable(same.orient1~feat.type.grp1+Var2, test.m.xt2)

dt.x <- data.table(as.data.frame((test.m.xt)))
expected <- independence_table(test.m.xt)
dt.x$expected <- as.matrix(expected)


expected <- independence_table(test.m.xt2)
x <- (test.m.xt2-expected)/sqrt(expected)

st <- structable(same.orient1 ~ ., data = test.m.xt)



lm1 <- loglm(~same.orient1+same.orient2+feat.type.grp1+feat.type.grp2, data=test.m.xt)
lm2 <- loglm(~(feat.type.grp1*same.orient1)+(feat.type.grp2*same.orient2), data=test.m.xt)
mosaic(~feat.type.grp1 + same.orient1 + feat.type.grp2 + same.orient2, data=test.m.xt, expected=~feat.type.grp1:same.orient1 + feat.type.grp2:same.orient2 + same.orient1:same.orient2 + feat.type.grp1:feat.type.grp2, legend=TRUE, gp=shading_Friendly, split_vertical=c(FALSE, FALSE, TRUE, FALSE))
mosaic(~feat.type.grp1 + same.orient1 + feat.type.grp2 + same.orient2, data=test.m.xt, expected=~feat.type.grp1:same.orient1 + feat.type.grp2:same.orient2 + same.orient1:same.orient2 + feat.type.grp1:feat.type.grp2, legend=TRUE, gp=shading_Friendly)
plot(ca(test.m.xt2), mass = TRUE, contrib = "absolute", map ="rowgreen", arrows = c(FALSE, TRUE))


reads.dt1.sample.counts.thru1.sum.featol.typesum.long[reads.dt1.sample.counts.thru1.sum.featol.typesum.long[, .(grp, feat.type.grp, same.orient)], .N, by=.(Sample.ID, feat.type.grp, same.orient, i.feat.type.grp, i.same.orient), on="grp", allow.cartesian=TRUE][, .N, by=.(Sample.ID, feat.type.grp, same.orient, i.feat.type.grp, i.same.orient)]

# look at the lncRNA promoter reads
as.gene <- reads.dt1.sample.counts.thru3.sum.featol[feat.type.grp=="lncRNA.promoters" & same.orient==FALSE & transcript_type=="antisense_RNA", lapply(.SD, sum), by=.(tx_name, gene_name, transcript_type, participant.ID, PNK, Sample.ID), .SDcols=sum.cols][, lapply(.SD, mean), by=.(tx_name, gene_name, transcript_type, PNK), .SDcols=sum.cols][order(-sample.num.uniq.reads.OQ)]$gene_name[1]

gene.ranges <- reduce(txdbFeatures.dt.c[gene_name==as.gene, GRanges(seqnames, IRanges(start, end), strand=strand, seqinfo=seqinfo(txdb.gencode))])
as.gene.tx <- biovizBase::crunch(txdb.gencode, which=gene.ranges)
colnames(values(as.gene.tx))[4] <- "model"
lvls <- levels(as.gene.tx$model)
lvls[lvls=="gap"] <- "intron"
as.gene.tx$model <- factor(ifelse(as.gene.tx$model=="gap", "intron", as.character(as.gene.tx$model)), levels=lvls)

grl <- split(as.gene.tx, as.gene.tx$tx_id)
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
ggplot() + geom_alignment(grl, type="model")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gr.txdb <- crunch(txdb, which = gene.ranges)
grl <- split(gr.txdb, gr.txdb$tx_id)

exons.by.gene <- exonsBy(txdb.gencode, by="gene")
tx.by.gene <- transcriptsBy(txdb.gencode, by="gene")
exons.ol <- exons.by.gene[unique(as.gene.tx$gene_id)]
tx.ol <- tx.by.gene[unique(as.gene.tx$gene_id)]
options(ucscChromosomeNames=FALSE)
gtrack <- GenomeAxisTrack()

reads.dt1.gr.thru3 <- reads.dt1.annot.all[thru.biofrag.stage==3 & PNK=="PNK", .N, by=.(seqnames, i.start, i.end, i.strand, Sample.ID)][,GRanges(seqnames, IRanges(i.start, i.end), strand=i.strand, seqinfo = seqinfo(txdb.gencode))]

cov.plus <- coverage(reads.dt1.gr.thru3[strand(reads.dt1.gr.thru3)=="+"])
cov.minus <- coverage(reads.dt1.gr.thru3[strand(reads.dt1.gr.thru3)=="-"])
cov.both <- coverage(reads.dt1.gr.thru3)
red.query <- reduce(unlist(tx.ol), ignore.strand=TRUE)
names(red.query) <- "myTrx"
q.plus <- as(cov.plus, "GRanges")
q.minus <- as(cov.minus, "GRanges")
q.both <- as(cov.both, "GRanges")
cov.med <- as(runmed(cov.both[red.query], k = 101), "GRanges")
seqlevels(cov.med) <- "myTrx"
cov.med.gr <- mapFromTranscripts(x = cov.med, transcripts = red.query)
cov.med.gr$score <- cov.med[cov.med.gr$xHits]$score
cov.med.gr.filt  <- cov.med.gr[cov.med.gr$score/width(cov.med.gr)>1]


atrack <- AnnotationTrack(unlist(exons.ol), group=names(unlist(exons.ol)), name="Gene Model")
atrack.l <- split(atrack, f=names(unlist(exons.ol)))
gene.info.query <- gtf.dt[names(atrack.l), .N, by=.(gene_id, gene_name, transcript_type), on="gene_id"]
for(i in seq_along(atrack.l)){
  names(atrack.l[[i]]) <- gene.info.query[names(atrack.l)[i], .N, by=gene_name, on="gene_id"]$gene_name
}
dtrack1 <- DataTrack(q.plus[q.plus%over%red.query], name = "Coverage (+)")
dtrack2 <- DataTrack(q.minus[q.minus%over%red.query], name = "Coverage (-)")
plotTracks(unlist(list(gtrack, atrack.l, dtrack2), recursive=FALSE), from = min(start(cov.med.gr.filt)), to=max(end(cov.med.gr.filt)), type="polygon")
txdbFeatures.dt.c[gtf.dt, gene_id:=gene_id, on=c(tx_name="transcript_id")]


transcripts.gr <- txdbFeatures.dt.c[feature_type=="transcripts", GRanges(seqnames, IRanges(start, end), strand, .SD), .SDcols=setdiff(colnames(txdbFeatures.dt.c), c("seqnames", "start", "end", "strand", "width"))]

setkey(txdbFeatures.dt.c.trx, seqnames, start, end)
scheme <- getScheme()
scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

cov.fun.reads.thru3 <- function(as.gene, show.full.range="AUTO"){
  gene.ranges <- txdbFeatures.dt.c[gene_name==as.gene, GRanges(seqnames, IRanges(start, end), strand=strand, seqinfo=seqinfo(txdb.gencode), gene_name=gene_name, tx_name=tx_name, transcript_type=transcript_type, feature_type=feature_type, gene_id=gene_id)]
  promoter.ranges <- gene.ranges[gene.ranges$feature_type=="promoters"]
  gene.ranges <- reduce(gene.ranges, ignore.strand=FALSE)
  
  promoter.track <- AnnotationTrack(promoter.ranges, id=paste0(as.gene, " promoter"), genome="hg38", name="promoters", featureAnnotation="id", fill="darkred", fontsize=10, col="red", lty=2)
  #start(gene.ranges) <- start(gene.ranges)-100
  #end(gene.ranges) <- end(gene.ranges)+100
  #gene.ranges <- trim(gene.ranges)
  as.gene.tx <- transcripts.gr[unique(queryHits(findOverlaps(transcripts.gr, gene.ranges, ignore.strand=TRUE)))]
  txdbFeatures.dt.c.f <- txdbFeatures.dt.c[unique(as.character(as.gene.tx$tx_name)), on="tx_name"]
  mRNA.f.tmp <- txdbFeatures.dt.c.f[transcript_type=="protein_coding"][c("cds", "fiveUTRs", "threeUTRs", "exons"), , on="feature_type"]
  mRNA.f <- mRNA.f.tmp[feature_type!="exons"]
  mRNA.f.exons <- mRNA.f.tmp[feature_type=="exons"]
  setkey(mRNA.f, tx_name, strand, start, end)
  setkey(mRNA.f.exons, tx_name, strand, start, end)
  mRNA.exon.ol <- unique(rbind(foverlaps(mRNA.f, mRNA.f.exons, type="start", which=TRUE), foverlaps(mRNA.f, mRNA.f.exons, type="end", which=TRUE)))[, .(yid=min(yid, na.rm=TRUE), ycount=sum(ifelse(!is.na(yid), 1, 0))), by=xid]
  mRNA.f[mRNA.exon.ol$xid, exon_name:=mRNA.f.exons$exon_name[mRNA.exon.ol$yid]]
  if(any(mRNA.exon.ol$ycount!=1)){
    warning("Unambiguous assignment of exon ids to feature. picking first one")
  }
  
  other.f <- txdbFeatures.dt.c.f[transcript_type!="protein_coding"][c("exons"), , on="feature_type"]
  setnames(other.f,
           c("seqnames","start", "end", "width", "strand", "transcript_type", "gene_id", "exon_name", "tx_name", "gene_name"),
           c("chromosome","start", "end", "width", "strand", "feature", "gene", "exon", "transcript", "symbol"))
  setnames(mRNA.f,
           c("seqnames","start", "end", "width", "strand", "gene_id", "exon_name", "tx_name", "gene_name"),
           c("chromosome","start", "end", "width", "strand", "gene", "exon", "transcript", "symbol"))
  mRNA.feats.replace <- c("protein_coding", "utr3", "utr5")
  names(mRNA.feats.replace) <- c("cds", "fiveUTRs", "threeUTRs")
  mRNA.f[, feature:=mRNA.feats.replace[feature_type]]
  mRNA.f[, protein.coding:="protein_coding"]
  other.f[, protein.coding:="noncoding"]
  feats.out <- subset(rbind(mRNA.f, other.f, fill=TRUE), select=c("chromosome","start", "end", "width", "strand", "feature", "gene", "exon", "transcript", "symbol", "protein.coding"))
  
  gtrack <- GenomeAxisTrack()
  
  exons.ol <- exons.by.gene[unique(as.gene.tx$gene_id)]
  tx.ol <- tx.by.gene[unique(as.gene.tx$gene_id)]
  
  red.query <- reduce(unlist(tx.ol), ignore.strand=TRUE)
  
  if(length(red.query)>1){
    warning("red.query has more than one range. Reducing further to one maximal range.")
    if(length(unique(seqnames(red.query)))>1){
      stop("red.query must have only one seqlevel")
    }
    red.query <- GRanges(seqnames = seqnames(red.query)[1], IRanges(start=min(start(red.query)), end=max(end(red.query))))
  }
  
  names(red.query) <- "myTrx"
  atrack <- GeneRegionTrack(feats.out)
  atrack.l <- split(atrack, f=feats.out$protein.coding)
  displayPars(atrack.l[["noncoding"]]) <- list(col="red", lwd=2)
  names(atrack.l[["noncoding"]]) <- "noncoding"
  names(atrack.l[["protein_coding"]]) <- "protein_coding"
  #print(names(atrack.l))
  gene.info.query <- txdbFeatures.dt.c[names(atrack.l), .N, by=.(gene_id, gene_name, transcript_type), on="gene_id"]
  q.plus.x <- q.plus[q.plus%over%red.query]
  q.minus.x <- q.minus[q.minus%over%red.query]
  q.cov.plus <- cov.plus[red.query]
  q.cov.minus <- cov.minus[red.query]
  q.cov.diff <- as(q.cov.plus-q.cov.minus, "GRanges")
  #y.max <- max(score(c(q.plus.x, q.minus.x)))
  y.max <- max(abs(score(c(q.cov.diff))))
  #dtrack1 <- DataTrack(q.plus.x, name = "Coverage (+)", type="polygon", ylim=c(0, y.max))
  #dtrack2 <- DataTrack(q.minus.x, name = "Coverage (-)", type="polygon", ylim=c(0, y.max))
  dtrack <- DataTrack(q.cov.diff, name = "Coverage", type="polygon", ylim=c(-1*y.max, y.max))
  cov.med <- as(runmed(cov.both[red.query], k = 101), "GRanges")
  seqlevels(cov.med) <- "myTrx"
  cov.med.gr <- mapFromTranscripts(x = cov.med, transcripts = red.query)
  cov.med.gr$score <- cov.med[cov.med.gr$xHits]$score
  cov.med.gr.filt  <- which(cov.med.gr$score/width(cov.med.gr)>1)
  which.start <- max(1, min(cov.med.gr.filt)-1)
  which.end <- min(length(cov.med.gr), max(cov.med.gr.filt)+1)
  from.s1 <- start(cov.med.gr[which.start])
  to.e1 <- end(cov.med.gr[which.end])
  from.s <- start(red.query)
  to.e <- end(red.query)
  w <- to.e-from.s
  w.s <- abs(from.s-from.s1)
  w.e <- abs(to.e-to.e1)
  if(w.s/w > 0.2 & show.full.range=="AUTO"){
    from.s <- from.s1 - round(0.1*w)
  }
  if(w.e/w > 0.2 & show.full.range=="AUTO"){
    to.e <- to.e1 + round(0.1*w)
  }
  #track.list <- unlist(list(gtrack, atrack.l["protein_coding"], promoter.track, atrack.l["noncoding"], dtrack1, dtrack2), recursive=FALSE)
  track.list <- unlist(list(gtrack, atrack.l["protein_coding"], promoter.track, atrack.l["noncoding"], dtrack), recursive=FALSE)
  track.list2 <- list(track.list=track.list, from.s=from.s, to.e=to.e)
  return(track.list2)
  
}
as.promoter.genes <- reads.dt1.sample.counts.thru3.sum.featol[feat.type.grp=="lncRNA.promoters" & same.orient==FALSE & transcript_type=="antisense_RNA", lapply(.SD, sum), by=.(tx_name, gene_name, transcript_type, participant.ID, PNK, Sample.ID), .SDcols=sum.cols][, lapply(.SD, mean), by=.(gene_name, tx_name, transcript_type, PNK), .SDcols=sum.cols][order(-sample.num.uniq.reads.OQ)][PNK=="PNK" & sample.num.uniq.reads.OQ>5]

as.promoter.genes.gr <- txdbFeatures.dt.c[as.promoter.genes$tx_name, on="tx_name"][feature_type=="promoters", GRanges(seqnames, IRanges(start, end), strand=strand, tx_name, gene_name)]
as.promoter.genes.nonOverlapping <- as.promoter.genes.gr[as.data.table(findOverlaps(as.promoter.genes.gr, as.promoter.genes.gr))[, w.min:=apply(.SD, 1, min)][, .(wmin.v=min(w.min)), by=queryHits][, .N, by=wmin.v]$wmin.v]

ncols=2
ngenes=6
nrows <- ngenes%/%ncols
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrows, ncols)))
for(i in 1:ngenes){
  this.gene <- as.promoter.genes.nonOverlapping$gene_name[i]
  pushViewport(viewport(layout.pos.col = ((i - 1)%%ncols) + 1, layout.pos.row = (((i) - 1)%/%ncols) + 1))
  tl <- cov.fun.reads.thru3(as.gene = this.gene, show.full.range = FALSE)
  plotTracks(tl[[1]], from = tl[[2]], to=tl[[3]], antisense_RNA="darkred", transcriptAnnotation="symbol", collapseTranscripts="meta", add=TRUE)
  popViewport(1)
}


queryRegionsSum2[, (tot.cols):=lapply(.SD, sum), by=.(feat.type.grp, participant.ID, PNK, Sample.ID, thru.biofrag.stage), .SDcols=sum.cols]
queryRegionsSum2[, percent.tot.readcount:=sample.read.count.OQ/sample.read.count.OQ.bothOrient]
queryRegionsSum2[, percent.tot.uniqread:=sample.num.uniq.reads.OQ/sample.num.uniq.reads.OQ.bothOrient]
g <- ggplot2::ggplot(queryRegionsSum2, aes(x = thru.biofrag.stage, y = percent.tot.uniqread)) + 
  geom_boxplot(aes(fill = factor(thru.biofrag.stage)), pos="dodge") +  
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_grid(feat.type.grp~PNK); 
g %+% queryRegionsSum2[same.orient==TRUE]

g <- ggplot2::ggplot(queryRegionsSum2, aes(x = factor(thru.biofrag.stage), fill=PNK, y = percent.tot.uniqread)) + 
  geom_boxplot(pos="dodge", width=0.75, aes(fill=PNK)) + geom_point(size=0.1, pos=position_jitterdodge(jitter.height = 0, jitter.width = 0.3, dodge.width = 0.75)) + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_wrap(~feat.type.grp) + geom_hline(yintercept = 0.5, linetype=2) +
  scale_y_continuous(labels=scales::percent_format(), breaks=seq(0,1,0.1), limits=c(0,1)) +
  labs(y="% Feature reads aligned in sense orientation"); g %+% queryRegionsSum2[same.orient==TRUE]

thru3.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru3.sum[overlaps.sRNA==FALSE], txdbFeatures = txdbFeatures.gr.types, by.col=c("Sample.ID"), val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)

thru2.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru2.sum[overlaps.sRNA==FALSE], txdbFeatures = txdbFeatures.gr.types, by.col="Sample.ID", val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)
thru1.queryRegionsSum <- summarizeQueryRegions.dt(reads.dt1.sample.counts.thru1.sum[overlaps.sRNA==FALSE], txdbFeatures = txdbFeatures.gr.types, by.col="Sample.ID", val.cols = value.columns, ignore.strand = TRUE, check.orientation = TRUE)
thru3.queryRegionsSum[, thru.biofrag.stage:=3]
thru2.queryRegionsSum[, thru.biofrag.stage:=2]
thru1.queryRegionsSum[, thru.biofrag.stage:=1]

queryRegionsSum <- rbind(thru1.queryRegionsSum, thru2.queryRegionsSum, thru3.queryRegionsSum)
queryRegionsSum.melt <- melt(queryRegionsSum, measure.vars = value.columns, variable.name = "measurement", value.name = "value")

queryRegionsSum.melt.tot <- queryRegionsSum.melt[feat=="total"]
queryRegionsSum.melt[queryRegionsSum.melt.tot, total:=i.value, on=c("Sample.ID", "thru.biofrag.stage", "measurement")]
queryRegionsSum.melt[, percent.total:=value/total]
queryRegionsSum.melt[, c("Participant.ID", "PNK"):=tstrsplit(Sample.ID, split=".", fixed=TRUE)]
queryRegionsSum.melt.f <- queryRegionsSum.melt[!is.na(same.orient) & feat!="total"]

feat.means <- queryRegionsSum.melt.f[measurement=="sample.cpm.read.OQ" & same.orient==TRUE, mean(value), by=feat]
setorder(feat.means, -V1)
queryRegionsSum.melt.f[, feat.f:=factor(feat, levels=feat.means$feat)]

queryRegionsSum.melt.f[, thru.biofrag.stage:=as.factor(thru.biofrag.stage)]
g <- ggplot2::ggplot(queryRegionsSum.melt.f, aes(x = feat.f, y = percent.total)) + 
  geom_boxplot(aes(fill = same.orient), pos="dodge") +  
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90), legend.pos="top") + facet_grid(paste0("Step", thru.biofrag.stage)~PNK)

g %+% queryRegionsSum.melt.f[measurement=="sample.read.count.OQ"]



exons <- exonsBy(txdb.gencode, by = "tx", use.names = TRUE)
tx <- GenomicFeatures:::.tidy_transcripts(txdb.gencode)
tx$tx_type <- tx.info[tx$tx_name]$transcript_type
tx$gene_name <- tx.info[tx$tx_name]$gene_name
ex_by_tx <- GenomicFeatures:::.exons_by_txids(txdb.gencode, mcols(tx)$tx_id)
tx.f <- tx[mRNA.lncRNA.txids$tx_row]
ex_by_tx.f <- ex_by_tx[mRNA.lncRNA.txids$tx_row]
introns_by_tx <- psetdiff(tx.f, ex_by_tx.f)
ans <- unlist(introns_by_tx, use.names = FALSE)
idx <- rep(seq_along(tx.f), lengths(introns_by_tx))
mcols(ans) <- mcols(tx.f)[idx, , drop = FALSE]
mRNA.lncRNA.intron.parts <- GenomicFeatures:::.break_in_parts(ans, linked.to.single.gene.only = FALSE)

ans <- unlist(ex_by_tx.f, use.names = FALSE)
idx <- rep(seq_along(tx.f), lengths(ex_by_tx.f))
mcols(ans) <- cbind(mcols(tx.f)[idx, , drop = FALSE], mcols(ans))
mRNA.lncRNA.exon.parts <- GenomicFeatures:::.break_in_parts(ans, linked.to.single.gene.only = FALSE)

tx.f <- tx[sRNA.txids$tx_row]
ex_by_tx.f <- ex_by_tx[sRNA.txids$tx_row]
ans <- unlist(ex_by_tx.f, use.names = FALSE)
idx <- rep(seq_along(tx.f), lengths(ex_by_tx.f))
mcols(ans) <- cbind(mcols(tx.f)[idx, , drop = FALSE], mcols(ans))
sRNA.exon.parts <- GenomicFeatures:::.break_in_parts(ans, linked.to.single.gene.only = FALSE)

rmsk.gr <- rmsk.dt[, GRanges(seqnames, IRanges(start, end), strand=strand, mcols=DataFrame(.SD)), .SDcols=c("tx_name", "rep.short", "gene_name", "transcript_type", "feature_type", "feature.type.simple")]
rmsk.gr.dj <- disjoin(rmsk.gr, with.revmap=TRUE)
ans_mcols <- lapply(mcols(rmsk.gr), function(col) {
  col <- unique(extractList(col, rmsk.gr.dj$revmap))
  col[!is.na(col)]
})


mcols(rmsk.gr.dj) <- DataFrame(ans_mcols)
gencode.seqinfo.df <- read.table("/data/genomic_data/Homo_sapiens/UCSC/GRCh38/STARIndex/chrNameLength.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
gencode.seqinfo <- Seqinfo(seqnames=gencode.seqinfo.df$V1, seqlengths=gencode.seqinfo.df$V2, isCircular=grepl("chrm", gencode.seqinfo.df$V1, ignore.case=TRUE), genome="hg38")

hg38_custom_mRNA.lncRNA.exon.parts <- mRNA.lncRNA.exon.parts
hg38_custom_mRNA.lncRNA.intron.parts <- mRNA.lncRNA.intron.parts
hg38_custom_sRNA.exon.parts <- sRNA.exon.parts
rmsk.gr.dj <- keepSeqlevels(rmsk.gr.dj, seqlevels(gencode.seqinfo), pruning.mode = "coarse")
seqinfo(rmsk.gr.dj) <- gencode.seqinfo

reads.dt1.sample.counts.uniq <- reads.dt1.sample.counts[thru.biofrag.stage==3]


annots = c("hg38_basicgenes", 'hg38_custom_mRNA.lncRNA.exon.parts', 'hg38_custom_mRNA.lncRNA.intron.parts', 'hg38_custom_sRNA.exon.parts', 'hg38_custom_rmsk')
annotations = build_annotations(genome = 'hg38', annotations = annots)



read_annotations.rs <- function (con, name, genome = NA, format, extraCols = character(), 
                                 ...) 
{
  if (missing(name)) {
    name = con
  }
  if (is.na(genome)) {
    genome = "genome"
  }
  gr <- get(con)
  protected_extraCols = c("gene_id", "symbol", "tx_id")
  extraCols <- setdiff(c("seqnames", "start", "end", "strand", protected_extraCols), colnames(gr))
  
  missing_extraCols = base::setdiff(protected_extraCols, names(extraCols))
  if (any(missing_extraCols == "gene_id")) {
    GenomicRanges::mcols(gr)$gene_id = NA
  }
  if (any(missing_extraCols == "symbol")) {
    GenomicRanges::mcols(gr)$symbol = NA
  }
  if (any(missing_extraCols == "tx_id")) {
    GenomicRanges::mcols(gr)$tx_id = NA
  }
  GenomicRanges::mcols(gr)$id = paste0(name, ":", seq_along(gr))
  GenomicRanges::mcols(gr)$type = sprintf("%s_custom_%s", genome, 
                                          name)
  GenomicRanges::mcols(gr) = GenomicRanges::mcols(gr)[, c("id", 
                                                          "tx_id", "gene_id", "symbol", "type")]
  annotatr_cache$set(sprintf("%s_custom_%s", genome, name), gr)
}
read_annotations.rs("mRNA.lncRNA.exon.parts", genome = "hg38")

#overlaps.m <- sparseM.from.dt(overlaps2.filt, "overlappingQuery", "gene_type_orient")
#overlaps.ig <- to.igraph(overlaps.m)
#overlaps.ig.bip <- bipartite.projection(to.igraph(overlaps.m), which = TRUE)
#overlaps.m.cp <- crossprod(overlaps.m)
#overlaps.jacc <- text2vec::sim2(t(overlaps.m), method="jaccard", norm = "none")

#overlaps.m.ngc <- as(overlaps.m, "ngCMatrix")
#s.im.x <- as(overlaps.m.ngc, "transactions")
#ident.set.id.x <- .Call(arules:::R_pnindex, s.im.x@data, NULL, FALSE)



df <- overlaps2.filt[, .(n.seqs=sum(n.seqs)), by=.(gene_type, aln.orientation, overlappingQuery)][, .(n.seqs=sum(n.seqs)), by=.(gene_type, aln.orientation)]
setnames(df, c("feature", "orient", "count"))
df$percent <- round(df$count / length(reads.dt1.gr) * 100, 1)
df <- df[order(count, decreasing = TRUE)]
ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 0.5), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(reads.dt1.gr), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))

txdbFeatures <- getTxdbFeaturesFromGRanges(gtf)
summary <- summarizeQueryRegions(queryRegions = queryRegions, 
                                 txdbFeatures = txdbFeatures)
df <- data.frame(summary)
df$percent <- round((df$count / length(queryRegions)), 3) * 100
df$feature <- rownames(df)
ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 3), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))

#Region coverage ----
cvgList <- calculateCoverageProfileList(queryRegions = queryRegions, 
                                        targetRegionsList = txdbFeatures, 
                                        sampleN = 10000)

ggplot2::ggplot(cvgList, aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(fill = 'lightgreen', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + theme_bw(base_size = 14) +
  facet_wrap( ~ feature, ncol = 3) 

# import sRNAnalysis profile counts ----

contam.profile.files <- dir(contam.dir, pattern="*profile", recursive = TRUE, full.names = TRUE)
cat.call <- paste0("cat ", paste(contam.profile.files, collapse=" "))
contam.profile.dt <- fread(cat.call, sep="\t", header=FALSE)

# Rescue the file name labels using the length of the inptu files
file.sizes <- fread(sub("cat", "wc -l", cat.call), header=FALSE, sep=" ")
file.sizes <- file.sizes[1:length(contam.profile.files)]
file.sizes[, end.row:=cumsum(V1)]
file.sizes[, start.row:=data.table::shift(cs, 1, 0, "lag")+1]
file.sizes[, File.Base.ID2:=gsub("_l00[1-9].*", "", basename(V2))]

contam.profile.dt[file.sizes$start.row, File.Base.ID2:=file.sizes$File.Base.ID2]
contam.profile.dt[, grp:=cumsum(!is.na(File.Base.ID2))]
contam.profile.dt[, File.Base.ID2:=File.Base.ID2[1], by=grp]

# Compare with exceRpt ----
excerpt.file.cpm <- "./exceRpt_comparison/ULMC135_Biofragmenta_excerpt_output/exrna-mtewa1/ULMC_Control_PNK/exceRptPipeline_v4.6.2/ULMC135_Biofragmenta_Comparison-2018-9-19-21_19_25/postProcessedResults_v4.6.3/ULMC135_Biofragmenta_Comparison-2018-9-19-21%3A19%3A25_exceRpt_gencode_ReadsPerMillion.txt"
excerpt.file.counts <- "./exceRpt_comparison/ULMC135_Biofragmenta_excerpt_output/exrna-mtewa1/ULMC_Control_PNK/exceRptPipeline_v4.6.2/ULMC135_Biofragmenta_Comparison-2018-9-19-21_19_25/postProcessedResults_v4.6.3/ULMC135_Biofragmenta_Comparison-2018-9-19-21%3A19%3A25_exceRpt_gencode_ReadCounts.txt"

excerpt.cpm <- fread(excerpt.file.cpm)
excerpt.counts <- fread(excerpt.file.counts)
file.name.col <- grep("_fastq$", colnames(excerpt.counts), value=TRUE)

excerpt.counts.melt <- melt(excerpt.counts, id.vars = "V1", variable.name = "file.name", value.name = "COUNT")
excerpt.cpm.melt <- melt(excerpt.cpm, id.vars = "V1", variable.name = "file.name", value.name = "CPM")
excerpt.counts.melt[excerpt.cpm.melt, CPM:=i.CPM, on=c("V1", "file.name")]
excerpt.counts.melt[, tot.counts.excerpt:=(COUNT*10^6)/CPM]
excerpt.counts.melt[, tot.counts.excerpt:=max(tot.counts.excerpt, na.rm = TRUE), by=file.name ]
excerpt.counts.melt[, c("gene.symbol", "feat.type"):=tstrsplit(V1[1], split=":", fixed=TRUE), by=V1]
excerpt.counts.melt[, File.Base.ID2:=sub("sample_", "", file.name[1]), by=file.name]
base.id.excerpt <- excerpt.counts.melt[, .N, by=file.name][, File.Base.ID2:=ALL.SAMPLE.INFO$File.Base.ID2[which(startsWith(sub("sample_", "", file.name), ALL.SAMPLE.INFO$File.Base.ID2))], by=file.name][]
excerpt.counts.melt[base.id.excerpt, File.Base.ID2:=i.File.Base.ID2, on="file.name"]
excerpt.counts.melt[ALL.SAMPLE.INFO, `:=`(participant.ID=i.participant.ID, PNK=i.PNK), on="File.Base.ID2"]
excerpt.counts.melt[, Sample.ID:=paste0(participant.ID, ".", PNK)]

all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135 <- subset(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2[participant.ID=="ULMC135"], select=setdiff(colnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2),"sample.read.count.uniq"))
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.cast <- dcast.data.table(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135, ...~thru.biofrag.stage, fun.aggregate = sum, fill=0)

thru.cols <- paste0("thru.", 1:3)
thru.incols <- as.character(1:3)
setnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.cast, thru.incols, thru.cols)
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.cast[, thru.2:=thru.2+thru.3]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.cast[, thru.1:=thru.1+thru.2]

all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast <- dcast.data.table(melt(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.cast, 
                                                                                     variable.name = "thru.stage.f",
                                                                                     value.name="sample.read.count.allmap",
                                                                                     measure.vars=thru.cols
), ...~rep.map.category.simple,
value.var="sample.read.count.allmap", fill=0, fun.aggregate = c)


all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast[transcript.type.dt, feat.type:=i.transcript_type, on= c("tx_name"="transcript_id")]
setnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast, c("gene_name"), "gene.symbol")
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast[, cpm.rank.biofrag.all:=rank(-1*ALL, na.last = TRUE, ties.method = "min"), by=.(feat.type, Sample.ID, thru.stage.f)]

# Sum over replicates
excerpt.counts.melt.sum <- excerpt.counts.melt[, .(CPM=sum(CPM), COUNT=sum(COUNT), tot.counts.excerpt=sum(tot.counts.excerpt)), by=.(V1, gene.symbol, feat.type, participant.ID, PNK, Sample.ID)]
#all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast[excerpt.counts.melt.sum, `:=`(CPM.excerpt=i.CPM, COUNT.excerpt=i.COUNT, tot.counts.excerpt=i.tot.counts.excerpt), on=c("gene.symbol", "feat.type", "Sample.ID")]

setkey(excerpt.counts.melt.sum, gene.symbol, feat.type, Sample.ID)
setkey(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast, gene.symbol, feat.type, Sample.ID)
excerpt.counts.melt.sum[, feat.type.simple:="other"]
excerpt.counts.melt.sum[unique(c("protein_coding", "antisense", lncRNA.features, tolower(lncRNA.features))), feat.type.simple:="mRNA.or.lncRNA", on="feat.type"]

excerpt.counts.melt.sum[, cpm.rank.by.feat:=rank(-1*CPM, na.last = TRUE, ties.method = "min"), by=.(Sample.ID, feat.type)]
excerpt.counts.melt.sum.merge.repquant.all <- merge(excerpt.counts.melt.sum[feat.type=="protein_coding"] , all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.recast, allow.cartesian=TRUE, by=c("gene.symbol", "feat.type", "PNK", "participant.ID", "Sample.ID"))
excerpt.counts.melt.sum.merge.repquant.all[, percent.repeat:=CISorTRANS/ALL]
excerpt.counts.melt.sum.merge.repquant.all[, all.repeat:=CISorTRANS==ALL]

excerpt.counts.melt.sum.merge.repquant <- excerpt.counts.melt.sum.merge.repquant.all
excerpt.counts.melt.sum.merge.repquant.max.finals <- dcast.data.table(excerpt.counts.melt.sum.merge.repquant, gene.symbol+tx_name+feat.type+Sample.ID~thru.stage.f, value.var="ALL", fill=0)
excerpt.counts.melt.sum.merge.repquant.max.finals.use <- excerpt.counts.melt.sum.merge.repquant.max.finals[order(-thru.3, -thru.2, -thru.1), .(tx_name=tx_name[1]), by=.(gene.symbol, feat.type, Sample.ID)]

excerpt.counts.melt.sum.merge.repquant.f <- excerpt.counts.melt.sum.merge.repquant[excerpt.counts.melt.sum.merge.repquant.max.finals.use, on=colnames(excerpt.counts.melt.sum.merge.repquant.max.finals.use)]
excerpt.counts.melt.sum.merge.repquant.f[, CPM.rank.atStage:=rank(-1*ALL,na.last = TRUE, ties.method = "min"), by=.(Sample.ID, thru.stage.f)]

g <- ggplot(excerpt.counts.melt.sum.merge.repquant.f[, ], aes(x=CPM, y=ALL+1, color=all.repeat)) + geom_point() + facet_grid(PNK~thru.stage.f) + scale_y_log10() + scale_x_log10() + scale_color_manual(values=c("grey", "red")) + theme_bw(); g

excerpt.counts.melt.sum.merge.repquant.f[, ALL.p1:=ALL+1]
ggplot(excerpt.counts.melt.sum.merge.repquant.f[, ], aes(x=PNK, fill=thru.stage.f, y=ALL+1)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), pos="dodge") + scale_y_log10() + theme_bw() 
# FIG3B Excerpt counts top 50 genes detected ----
excerpt.counts.melt.sum.merge.repquant.f[, CPM.rank.atStage.group:=cut(CPM.rank.atStage, breaks = c(0,50,500,5000,50000))]
excerpt.counts.melt.sum.merge.repquant.f[, is.top.50:=cpm.rank.by.feat<=50]
g <- ggplot(excerpt.counts.melt.sum.merge.repquant.f[cpm.rank.by.feat<=50, ], aes(x=PNK, y=ALL.p1)) + 
  geom_boxplot(pos=position_dodge(width=0.8), width=0.5, size=0.5, outlier.colour = NA, aes(pos=thru.stage.f)) + 
  ggbeeswarm::geom_quasirandom(dodge.width = 0.8, alpha=0.7, size=0.2, aes( group=thru.stage.f, color=CPM.rank.atStage.group)) + scale_colour_ordinal(direction = -1) +
  scale_y_log10(breaks=10^seq(0,6,1), labels = scales::trans_format("log10", scales::math_format())) + 
  theme_bw(base_size = 6) + labs(x=NULL, y="Read Count", color="Gene Rank")
 
g1 <- g + theme(legend.position = "bottom", axis.text = element_text(size=6)) ; g1

save_plot(filename = "./output/figures/main/FIG3A_Top50_ExcerPt_though_filter_legend_protein_coding.pdf",  plot = g1, base_height = 75, base_width = 75, units="mm")

g1 <- g +theme(legend.position = "none", panel.grid = element_blank() , axis.text = element_text(size=6)) ; g1

save_plot(filename = "./output/figures/main/FIG3A_Top50_ExcerPt_though_filter_NOlegend_protein_coding.pdf",  plot = g1, base_height = 40, base_width = 50, units="mm")
p <- ggboxplot(excerpt.counts.melt.sum.merge.repquant.f[cpm.rank.by.feat<=50, ], x = "thru.stage.f", y = "ALL.p1", color = "thru.stage.f", palette = "npg", add = "jitter", facet.by = "PNK") + scale_y_log10(breaks=10^seq(-1,6,1)) + theme(text=element_text(size=6, colour="black"))
p2 <- p + stat_compare_means(aes(group = thru.stage.f,  y = CPM.rank.atStage), paired=TRUE, size=1, ref.group = "thru.1")
save_plot(filename = "./output/figures/main/FIG3A_Top50_ExcerPt_though_filter_legend_protein_coding_wilcox.pdf", plot = p2)

#PNK thru.stage.f (0,50] (50,500] (500,5e+03] (5e+03,5e+04]
#1: NoPNK       thru.1     29       10           7             2
#2: NoPNK       thru.2      7        2          39             0
#3: NoPNK       thru.3      7        3          38             0
#4:   PNK       thru.1     22       18           7             2
#5:   PNK       thru.2      8        8          18            15
#6:   PNK       thru.3      7        5          21            16


excerpt.counts.melt.sum.merge.repquant.all.max.finals <- dcast.data.table(excerpt.counts.melt.sum.merge.repquant.all, gene.symbol+tx_name+feat.type+Sample.ID~thru.stage.f, value.var="ALL", fill=0)
excerpt.counts.melt.sum.merge.repquant.all.max.finals.use <- excerpt.counts.melt.sum.merge.repquant.all.max.finals[order(-thru.3, -thru.2, -thru.1), .(tx_name=tx_name[1]), by=.(gene.symbol, feat.type, Sample.ID)]
excerpt.counts.melt.sum.merge.repquant.all.use.f <- excerpt.counts.melt.sum.merge.repquant.all[excerpt.counts.melt.sum.merge.repquant.all.max.finals.use, on=colnames(excerpt.counts.melt.sum.merge.repquant.max.finals.use)][COUNT>1 & thru.stage.f=="thru.1" & ALL>0]



g <- ggplot(excerpt.counts.melt.sum.merge.repquant.f[thru.stage.f=="thru.1" & cpm.rank.by.feat<=50], aes(x=CPM, color=PNK, y=percent.repeat)) + 
  geom_point(size=0.2) + 
  scale_x_log10() + 
  scale_y_continuous(labels = scales::percent_format()) + theme_bw(base_size = 6) + labs(x="ExceRpt CPM", y="% Reads Biofrag-labled sRNA / Repeat", color="Filtering Step"); g
g2 <- ggExtra::ggMarginal(g+ theme(legend.position = "none", axis.text = element_text(size=6), panel.grid = element_blank()), groupFill = TRUE, type="histogram", margins="y") ; g2
save_plot(filename = "./output/figures/main/FIG3B_Top50_ExcerPt_PercentRep2_proteinCoding_noLegend.pdf", plot=g2, base_height = 60, base_width = 75, units="mm")
g2 <- ggExtra::ggMarginal(g+ theme(legend.position = "top", axis.text = element_text(size=6), panel.grid = element_blank()), groupFill = TRUE, type="histogram", margins="y")
save_plot(filename = "./output/figures/main/FIG3B_Top50_ExcerPt_PercentRep2_proteinCoding_Legend.pdf", plot=g2, base_height = 60, base_width = 75, units="mm")
ggExtra::ggMarginal(g,  groupFill = TRUE, type="density")


ggplot(excerpt.counts.melt.sum.merge.repquant.f, aes(x=thru.stage.f, fill=ALL>0)) + geom_bar(pos="stack")  + scale_y_continuous() + theme_bw() + theme(legend.position = "bottom") + labs(x="Biofragmenta Filtering Stage", y="# mRNAs", fill="Biofragmenta Detected (COUNT>0)") + facet_wrap(~PNK) + labs("")


excerpt.counts.melt.sum.merge.repquant.caststg <- dcast.data.table(excerpt.counts.melt.sum.merge.repquant, gene.symbol+feat.type+PNK+participant.ID+Sample.ID+CPM~thru.stage, value.var="median.cpm.ALL", fun.aggregate = sum, fill=0L)
setnames(excerpt.counts.melt.sum.merge.repquant.caststg, c("1", "2", "3"), paste0("CPM.Stage", c("1", "2", "3")))

genes.gt.10.ulmc135 <- dcast.data.table(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135[rep.map.category.simple=="ALL"], PNK+participant.ID+tx_name+gene_name~thru.biofrag.stage, value.var="sample.read.count.allmap", fun.aggregate = sum, fill=0)
setnames(genes.gt.10.ulmc135, c("1", "2", "3"), paste0("Count.Stage", c("1", "2", "3")))
genes.gt.10.ulmc135[transcript.type.dt, feat.type:=i.transcript_type, on= c("tx_name"="transcript_id")]
genes.gt.10.ulmc135[, max.stage3:=max(Count.Stage3), by=.(PNK, participant.ID, gene_name, feat.type)]
genes.stage3dropped.ulmc135 <- genes.gt.10.ulmc135[Count.Stage1>=10 & max.stage3==0]

ggplot(excerpt.counts.melt.sum.merge.repquant, aes(x=cpm.rank.by.feat, y=percent.mapping.repeats)) + geom_point() + facet_grid(PNK~paste0("Step",thru.stage)) + theme_bw() + theme(legend.position = "top") + scale_color_gradient2(low = "blue", mid = "white", high = "red", breaks=seq(0, 1, 0.1), midpoint=0.5) + scale_y_continuous(labels=scales::percent, limits = c(0,1))
ggplot(excerpt.counts.melt.sum.merge.repquant[thru.stage==3], aes(x=CPM+1, y=median.cpm.ALL+1)) + geom_jitter() + facet_grid(PNK~paste0("Step",thru.stage)) + theme_bw() + theme(legend.position = "top") + scale_color_gradient2(low = "blue", mid = "white", high = "red", breaks=seq(0, 1, 0.1), midpoint=0.5) + scale_y_log10(breaks=10^seq(-1,5,1)) + scale_x_log10(breaks=10^seq(-1,5,1))


dropped.but.detected.genes <- merge(excerpt.counts.melt.sum[feat.type=="protein_coding", .(CPM=max(CPM), cpm.rank.by.feat=min(cpm.rank.by.feat)), by=.(gene.symbol, PNK)], genes.stage3dropped.ulmc135[feat.type=="protein_coding", .(N=.N, gene.symbol=gene_name), by=.(PNK, gene_name)], by=c("PNK", "gene.symbol"))
ggplot(dropped.but.detected.genes, aes(x=PNK, y=CPM)) + geom_violin() + theme_bw() + scale_y_log10(breaks=10^seq(-1,5,1)) 

excerpt.counts.melt.sum.merge.repquant.use <- excerpt.counts.melt.sum.merge.repquant[feat.type=="protein_coding"][genes.stage3dropped.ulmc135[feat.type=="protein_coding", .N, by=.(gene_name, PNK)], on=c("gene.symbol"="gene_name", "PNK")][!is.na(feat.type)]
excerpt.counts.melt.sum.merge.repquant.rank <- excerpt.counts.melt.sum.merge.repquant[, .(maxcpm=max(CPM)), by=.(gene.symbol, feat.type, PNK)][, rank.cpm:=rank(-1*maxcpm), by=PNK][]
excerpt.counts.melt.sum.merge.repquant.rank.top50 <- excerpt.counts.melt.sum.merge.repquant.rank[rank.cpm<=50]
excerpt.counts.melt.sum.merge.repquant[excerpt.counts.melt.sum.merge.repquant.rank.top50, on=c("gene.symbol", "feat.type", "PNK")]

ggplot(excerpt.counts.melt.sum.merge.repquant.use, aes(x=paste0("Step",thru.stage), y=CPM)) + geom_boxplot(outlier.color = NA) + geom_jitter() + facet_grid(PNK~.) + theme_bw() + theme(legend.position = "top") + scale_color_gradient2(low = "blue", mid = "white", high = "red", breaks=seq(0, 1, 0.1))


ggplot(excerpt.counts.melt.sum.merge.repquant[feat.type=="protein_coding" & thru.stage==3], aes(y=CPM, color= percent.mapping.repeats, x=median.cpm.ALL)) + geom_point(alpha=0.3) + facet_grid(PNK~paste0("Step",thru.stage)) + theme_bw() + theme(legend.position = "top") + scale_color_gradient2(low = "blue", mid = "white", high = "red", breaks=seq(0, 1, 0.1), midpoint=0.5)


ggplot(excerpt.counts.melt.sum.merge.repquant[feat.type=="protein_coding" & CPM>=1], aes(y=CPM, color=PNK, x=percent.mapping.repeats)) + geom_point(alpha=0.3) + facet_grid(PNK~paste0("Step",thru.stage)) + scale_y_log10() + scale_x_continuous(labels=scales::percent, limits = c(0,1)) + theme_bw() + theme(legend.position = "top")

all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135[transcript.type.dt, feat.type:=i.transcript_type, on= c("tx_name"="transcript_id")]
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all <- all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135[rep.map.category.simple=="ALL"]
setorderv(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all, c("Sample.ID", "tx_name", "thru.biofrag.stage"), c(1, 1, -1))
all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all[, sample.read.count.allmap.thrustage:=cumsum(sample.read.count.allmap), by=.(Sample.ID, tx_name)]

setnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all, "gene_name", "gene.symbol")
excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all <- merge(excerpt.counts.melt.sum, all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all, allow.cartesian=TRUE, by=intersect(colnames(excerpt.counts.melt.sum), colnames(all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.max)))
excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use <- excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all[feat.type=="protein_coding"][!is.na(feat.type)]
excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use[, cpm.rank.by.feat.biofrag:=rank(-sample.read.count.allmap), by=.(Sample.ID, thru.biofrag.stage)]

# percentile rank funciton ----
percentilerank<-function(x){
  rx<-rle(sort(x))
  smaller<-cumsum(c(0, rx$lengths))[seq(length(rx$lengths))]
  larger<-rev(cumsum(c(0, rev(rx$lengths))))[-1]
  rxpr<-smaller/(smaller+larger)
  rxpr[match(x, rx$values)]
}
excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use[, cpm.percentrank.by.feat.biofrag:=percentilerank(sample.read.count.allmap), by=.(Sample.ID, thru.biofrag.stage)]
excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use[, cpm.percentrank.by.feat.excerpt:=percentilerank(CPM), by=.(Sample.ID)]




g <- ggplot(excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use, aes(x=factor(thru.biofrag.stage, c("1", "2", "3")), y=sample.read.count.allmap.thrustage+1, group=gene.symbol))  + stat_summary(fun.y=median, geom="line") + stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max, aes(color=factor(thru.biofrag.stage))) + facet_grid(PNK~.) + theme_bw() + scale_y_log10()
g %+% excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use[cpm.rank.by.feat<=50]
g <- ggplot(excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use, aes(x=as.numeric(thru.biofrag.stage), y=sample.read.count.allmap.thrustage+1, group=tx_name))  + geom_line() + geom_point(aes(color=factor(thru.biofrag.stage, c("1", "2", "3")))) + facet_grid(PNK~.) + theme_bw() + scale_y_log10()
cpm.rank.max=50
dt.top <- excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use[cpm.rank.by.feat<=cpm.rank.max]
dt.use <- rbind(dt.top, 
                melt(dcast.data.table(dt.top,
                                      gene.symbol+feat.type+participant.ID+PNK+Sample.ID+V1+rep.map.category.simple+tx_name+tot.count.sample~thru.biofrag.stage, 
                                      fill=0,
                                      fun.aggregate = max, 
                                      value.var = "sample.read.count.allmap"), 
                     measure.vars = c("1", "2", "3"), variable.name = "thru.biofrag.stage", value.name="sample.read.count.allmap")[sample.read.count.allmap==0], fill=TRUE)
dt.use.rank <- rbind(dt.top, 
                     melt(dcast.data.table(dt.top,
                                           gene.symbol+feat.type+participant.ID+PNK+Sample.ID+V1+rep.map.category.simple+tx_name+tot.count.sample~thru.biofrag.stage, 
                                           fill=0.0,
                                           fun.aggregate = max, 
                                           value.var = "cpm.percentrank.by.feat.biofrag"), 
                          measure.vars = c("1", "2", "3"), variable.name = "thru.biofrag.stage", value.name="cpm.percentrank.by.feat.biofrag")[cpm.percentrank.by.feat.biofrag==0.0], fill=TRUE)
dt.top.cast <- dcast.data.table(dt.top,
                                gene.symbol+feat.type+participant.ID+PNK+Sample.ID+V1+rep.map.category.simple+tx_name+tot.count.sample~thru.biofrag.stage, 
                                fill=0.0,
                                fun.aggregate = max, 
                                value.var = c("cpm.percentrank.by.feat.biofrag", "cpm.percentrank.by.feat.excerpt"))
dt.top.use <- melt(dt.top.cast, 
                   measure.vars = list(grep("cpm.percentrank.by.feat.biofrag", colnames(dt.top.cast), value=TRUE),
                                       grep("cpm.percentrank.by.feat.excerpt", colnames(dt.top.cast), value=TRUE)),
                   value.name=c("cpm.percentrank.by.feat.biofrag", "cpm.percentrank.by.feat.excerpt"),
                   variable.name="thru.biofrag.stage")    
dt.top.use[, max.rank.gene.excerpt:=max(cpm.percentrank.by.feat.excerpt), by=.(gene.symbol, Sample.ID)]

g <- ggplot(dt.top.use, aes(x=cpm.percentrank.by.feat.excerpt, y=cpm.percentrank.by.feat.biofrag, group=gene.symbol))  + geom_jitter() + theme_bw() 
g2 <- g %+% dt.top.use[cpm.percentrank.by.feat.excerpt>0.0]



g <- ggplot(dt.use, aes(x=as.numeric(thru.biofrag.stage), y=sample.read.count.allmap.thrustage+1, group=gene.symbol))  + stat_summary(fun.y=median, geom="line") + stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max, aes(color=factor(thru.biofrag.stage))) + facet_grid(PNK~.) + theme_bw() + scale_y_log10()
g %+% dt.use
g <- ggplot(excerpt.counts.all.alns.mRNA.linRNA.exon.ol.thru1.gene.quant.u2.135.all.use, aes(x=as.numeric(thru.biofrag.stage), y=sample.read.count.allmap.thrustage.percentmax, group=tx_name))  + geom_line() + geom_point(aes(color=factor(thru.biofrag.stage, c("1", "2", "3")))) + facet_grid(PNK~.) + theme_bw()
g %+% dt.use
ggplot(dt.use, aes(x=as.numeric(thru.biofrag.stage), y=sample.read.count.allmap.thrustage+1, group=gene.symbol))  + stat_smooth() + facet_grid(PNK~.) + theme_bw() + scale_y_log10()


ggplot(reads.dt1.sample.counts.all.aligns.from.exon.aligned[thru.biofrag.stage==3, .N, by=.(origSequence, qwidth, Sample.ID, PNK, sample.read.count)][, .(count=sum(sample.read.count)), by=.(qwidth, Sample.ID, PNK)][ , percent.tot:=count/sum(count), by=Sample.ID][], aes(x=qwidth, y=percent.tot, group=PNK, color=PNK)) + stat_summary(size=0.2) + stat_smooth(span=0.25, fullrange = FALSE) + theme_bw() + scale_x_continuous(breaks=seq(0,80,5)) + scale_y_continuous(label=scales::percent, breaks = seq(-1, 1, by = 0.01)) + labs(y="Unique mRNA Exon-mapped reads\n(% total exon-mapped", x="Fragment Length (nt)") + theme(legend.position = "top")

# FIG 4B. Read length distribution for exonic reads ----
g <- ggplot(reads.dt1.sample.counts.all.aligns.from.exon.aligned[thru.biofrag.stage==3 & PNK=="PNK", .N, by=.(origSequence, qwidth, Sample.ID, PNK, sample.read.count)][, .(count=sum(sample.read.count)), by=.(qwidth, Sample.ID, PNK)][ , percent.tot:=count/sum(count), by=Sample.ID][], aes(x=qwidth, y=percent.tot)) + geom_point(size=0.1, alpha=0.7) + stat_smooth(span=0.25, size=0.5, fullrange = FALSE) + theme_bw(base_size = 6) + scale_x_continuous(breaks=seq(0,80,5)) + scale_y_continuous(label=scales::percent, breaks = seq(-1, 1, by = 0.01)) + labs(y="Unique mRNA Exon-mapped reads\n(% total exon-mapped", x="Fragment Length (nt)") + theme(legend.position = "top", panel.grid = element_blank(), axis.text=element_text(size=6))
save_plot(filename = "./output/figures/main/FIG4B_ReadLengthDistrPNK.pdf", plot = g, base_height = 50, base_width = 75, units="mm")

# Read length distribution by gene ----
pnk.genes.thru3.for.lengths <- reads.dt1.exon.ol[thru.biofrag.stage==3 & overlaps.sRNA.indir==0 & PNK=="PNK", .(sample.read.count=sample.read.count[1]), by=.(origSequence, participant.ID, gene_name, qwidth)]
pnk.genes.thru3.for.lengths.sum <- pnk.genes.thru3.for.lengths[, .(sample.gene.len.read.count=sum(sample.read.count)), by=.(participant.ID, gene_name, qwidth)]
pnk.genes.thru3.for.lengths.sum[, tot.gene.count.sample:=sum(sample.gene.len.read.count), by=.(participant.ID, gene_name)]
pnk.genes.thru3.for.lengths.sum[, gene.len.count.sample.perc.tot:=sample.gene.len.read.count/tot.gene.count.sample]


pnk.genes.thru3.for.lengths.all <- reads.dt1.exon.ol[thru.biofrag.stage==3 & overlaps.sRNA.indir==0 & PNK=="PNK", .(sample.read.count=sample.read.count[1]), by=.(origSequence, participant.ID, gene_name, qwidth)][, .(sample.read.count=sum(sample.read.count)), by=.(origSequence, gene_name, qwidth)]
pnk.genes.thru3.for.lengths.all.sum <- pnk.genes.thru3.for.lengths.all[, .(gene.len.read.count=sum(sample.read.count), uniq.gene.len.read.count=.N), by=.(gene_name, qwidth)]
pnk.genes.thru3.for.lengths.all.sum[, `:=`(tot.gene.count=sum(gene.len.read.count),
                                           tot.gene.unique.count=sum(uniq.gene.len.read.count)), by=.(gene_name)]
pnk.genes.thru3.for.lengths.all.sum[, gene.len.count.perc.tot:=gene.len.read.count/tot.gene.count]

pnk.genes.thru3.for.lengths.all.ranks <- pnk.genes.thru3.for.lengths.all.sum[, .(tot.gene.unique.count=max(tot.gene.unique.count)), by=gene_name][, tot.gene.unique.count.rank:=rank(-tot.gene.unique.count)][]

# Fig 4B-2 Top 50 genes read length dist ----
pnk.genes.thru3.for.lengths.all[pnk.genes.thru3.for.lengths.all.ranks, tot.gene.unique.count.rank:=i.tot.gene.unique.count.rank, on="gene_name"]
top.gene.lengths.plot <- pnk.genes.thru3.for.lengths.all[tot.gene.unique.count.rank<=50]
gene.order <- top.gene.lengths.plot[, .(mean.width=median(qwidth)), by=gene_name][order(-mean.width),]$gene_name
top.gene.lengths.plot[, gene_name.f:=factor(gene_name, levels=rev(gene.order))]
g <- ggplot(top.gene.lengths.plot, aes(x=gene_name.f, y=qwidth)) + 
  geom_boxplot(outlier.size = 0.5, outlier.stroke = 0, size=0.25) + 
  coord_flip() + 
  theme_bw(base_size = 6) + theme(panel.grid = element_blank(), axis.text=element_text(size=6)) +
  labs(x="Top 50 Expressed Genes", y="Read Length")
save_plot(filename = "./output/figures/main/FIG4B_ALT_ReadLengthDistrPNK_Top50Genes.pdf", plot = g, base_height = 100, base_width = 75, units="mm")


# Fig 4C Abundance VS number of individuals ----
pnk.genes.thru3.for.lengths.sum.tot <- pnk.genes.thru3.for.lengths.sum[, .N, by=.(participant.ID, gene_name, tot.gene.count.sample)][, .(tot.gene.count=sum(tot.gene.count.sample),
                                                                                                                                         n.individuals.expressed=.N), by=gene_name]

ggplot(pnk.genes.thru3.for.lengths.sum.tot, aes(x=factor(n.individuals.expressed), y=tot.gene.count)) + geom_boxplot() + scale_y_log10()
g <- ggplot(pnk.genes.thru3.for.lengths.sum.tot, aes(x=factor(n.individuals.expressed), y=tot.gene.count/n.individuals.expressed)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size=0.25, fill="grey") + 
  scale_y_log10() + labs(x="# Participants Detected", y="gene abundance (mean read counts)") + theme_bw(base_size = 6) + theme(panel.grid = element_blank(), axis.text=element_text(size=6))
save_plot(filename = "./output/figures/main/FIG4C_Gene_Abundance_VS_Individuals_Detected.pdf", plot = g, base_height = 50, base_width = 60, units="mm")



# Number of genes detected ----
reads.thru3.gene.sample.count <- reads.dt1.exon.ol[thru.biofrag.stage==3 & is.unique.map==TRUE & overlaps.sRNA.indir==0, .N, by=.(gene_name, feat.type.grp, origSequence, Sample.ID, participant.ID, PNK, tot.uniq.count.sample, tot.count.sample)]
thru3.gene.sample.count <- reads.thru3.gene.sample.count[, .(gene.unique.read.count.sample=.N), by=.(gene_name, feat.type.grp, Sample.ID,participant.ID, PNK, tot.uniq.count.sample, tot.count.sample)]
thru3.gene.sample.count[, `:=`(pnk.samples.detected=sum(ifelse(PNK=="PNK", 1, 0)), nopnk.samples.detected=sum(ifelse(PNK=="NoPNK", 1, 0)), n.samples.detected=.N), by=.(gene_name, feat.type.grp)]
gene.count.by.n.samples.detected <- thru3.gene.sample.count[, .N, by=.(gene_name, pnk.samples.detected, nopnk.samples.detected)][, xtabs(~pnk.samples.detected+nopnk.samples.detected)]
min.uniq.reads=2
thru3.gene.sample.count.f <- thru3.gene.sample.count[gene.unique.read.count.sample>=min.uniq.reads]
thru3.gene.sample.count.f[, `:=`(pnk.samples.detected=sum(ifelse(PNK=="PNK", 1, 0)), nopnk.samples.detected=sum(ifelse(PNK=="NoPNK", 1, 0)), n.samples.detected=.N), by=.(gene_name, feat.type.grp)]
gene.count.by.n.samples.detected.f <- thru3.gene.sample.count.f[, .N, by=.(gene_name, pnk.samples.detected, feat.type.grp, nopnk.samples.detected)][, xtabs(~pnk.samples.detected+nopnk.samples.detected+feat.type.grp)]



thru3.gene.sample.count[, median.perc.count:=median(gene.unique.read.count.sample/tot.uniq.count.sample), by=.(gene_name, feat.type.grp)]
gene.uniquemap.uniqread.count.tbl <- dcast.data.table(thru3.gene.sample.count, gene_name+feat.type.grp+n.samples.detected+pnk.samples.detected+nopnk.samples.detected+median.perc.count~Sample.ID, fun.aggregate = max, fill=0, value.var = "gene.unique.read.count.sample")[order(-n.samples.detected, -median.perc.count)]
fwrite(gene.uniquemap.uniqread.count.tbl, file = "unique_read_count_table_uniqmapped_ULMC_control.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


thru3.gene.sample.count.N <- thru3.gene.sample.count[, .N, by=.(gene_name, n.samples.detected, pnk.samples.detected, nopnk.samples.detected)][, .N, by=.(n.samples.detected, pnk.samples.detected, nopnk.samples.detected)]