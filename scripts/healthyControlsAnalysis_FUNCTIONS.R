# FUNCTIONS for healthyControlsAnalysis.R

# Functions ----
# FUNCTION:  percentile rank funciton ----
percentilerank<-function(x){
  rx<-rle(sort(x))
  smaller<-cumsum(c(0, rx$lengths))[seq(length(rx$lengths))]
  larger<-rev(cumsum(c(0, rev(rx$lengths))))[-1]
  rxpr<-smaller/(smaller+larger)
  rxpr[match(x, rx$values)]
}

# FUNCTION: Import bam file and convert to data.table ----
bam.to.dt.f <- function(bam.file, bam.params=param){
  readsGA <- readGAlignments(bam.file, param=param)
  dt <- as.data.table(readsGA)
  #dt <- dt[grepl("N")]
  return(dt)
}


# FUNCTION: Local evaluation of a det of expressions. Returns last call.
local2 <- function(expr) {
  env <- env(caller_env())
  eval_bare(enexpr(expr), env)
}

# FUNTION to optionally read or write file. Default sets read.file=TRUE so that pre-computed files 
# available on Github are used. Use read.file=FALSE to output a new file from the run
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


# FUNCTION to create a sparse matrix from data.table. 
# i and j must be a character string specifying the column names for the resulting nodes. 
# make.binary. If FALSE, the output matrix contains the number of i (row) j (col) instances in the input data frame/table. If TRUE, matrix outputs 1
#
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

# FUNCTION to create a bipartite graph (igraph package) from a data frame or matrix.
# x can be a data.table, data.frame or matrix.
# If x is a data.table or data.frame, i and j must be a character string specifying the column names for the resulting nodes. 
# Calls sparseM.from.dt, if data table/frame is provided
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

# FUNCTION:  Summarizes overlaps as in RCAS, but using data.table foverlaps. 
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

# FUNCTION:  Same as above, but dont summarize counts
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

# FUNCTION: Mutual information stats (from https://tm4ss.github.io/docs/Tutorial_5_Co-occurrence.html) ----
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