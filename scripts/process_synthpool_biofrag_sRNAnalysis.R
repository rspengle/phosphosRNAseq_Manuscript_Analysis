library(data.table)
# Function for analyzing synthetic pool samples 
process_synthpool_biofrag_sRNAnalysis <- function(biofragmenta.outdir,
                                                  library.run.type=c("TruSeq", "4N", "6N"),
                                                  univec.info,
                                                  sample.info.raw.fl,
                                                  library.type.info.fl=""){
  
  
  message(timestamp())
  tstart <- Sys.time()
  message(paste0(tstart, "STARTING sRNAnalysis: ", biofragmenta.outdir))
  sRNA.out.basedir <- paste0(biofragmenta.outdir, "/02-map_synthpool")
  
  if(library.type.info.fl==""){
    warning("No Library info type file is provided. Will import all libraries")
    filt.lib=FALSE
  } else{
    if(length(library.run.type)>1){
      warning("If library type info file is provided, then one library.run.type must also be selected")
      stop("Select one library run type: library.run.type=4N or TruSeq")
    }
    message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "library type"))
    filt.lib=TRUE
    library.type.info <- fread(library.type.info.fl)
    if(!any(colnames(library.type.info)=="lib.method")){
      stop("Library.type.info file must contain a column named lib.method, with 4N or TruSeq as labels")
    }
    if(!any(library.type.info$lib.method%in%library.run.type)){
      stop("The lib.method column in the library.type.info file must contain the lib.method value given.")
    }
    library.type.info <- subset(library.type.info, lib.method%in%library.run.type)
    message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  }
  
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "UNIVEC"))
  univec.dt <- fread(univec.info, header=FALSE, sep="\t")
  setnames(univec.dt, c("Aln.SeqID", "Seq.Info"))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "Sample Info"))
  sample.info.dt <- fread(sample.info.raw.fl)
  
  sample.info.dt[, Subpool.ID:=ifelse(Source=="U01FULL", 
                                      "SynthRNA476v1",
                                      ifelse(grepl("^Subpool", Source),
                                             sub("Subpool", "Pool", Source),
                                             Source))]
  setkey(sample.info.dt, IndexNum)
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  
  feature.count.fls <- dir(path = sRNA.out.basedir, pattern = "*.profile", full.names = TRUE, recursive = TRUE)
  unaligned.count.fls <- dir(path = sRNA.out.basedir, pattern = "*_unMatch.fa", full.names = TRUE, recursive = TRUE)
  processed.read.fls <- dir(path = sRNA.out.basedir, pattern = "*_Processed.fa", full.names = TRUE, recursive = TRUE)
  
  if(filt.lib==TRUE){
    sample.info.dt <- sample.info.dt[IndexNum%in%library.type.info$IndexNum]
    feature.count.fls <- feature.count.fls[sub("_L00[1-4]_R1_001_.*", "", basename(feature.count.fls))%in%library.type.info$File.Base.ID]
    unaligned.count.fls <- unaligned.count.fls[sub("_L00[1-4]_R1_001_.*", "", basename(unaligned.count.fls))%in%library.type.info$File.Base.ID]
    processed.read.fls <- processed.read.fls[sub("_L00[1-4]_R1_001_.*", "", basename(processed.read.fls))%in%library.type.info$File.Base.ID]
    if(!all.equal(dim(sample.info.dt)[1], length(feature.count.fls), length(unaligned.count.fls), length(processed.read.fls))){
      stop("NOT ALL FILTERED FILES ARE EQUAL SIZES")
    } else if(length(feature.count.fls)==0){
      stop("ALL FILES FILTERED OUT")
    }
  }
  
  
  all.base.names <- sub("_Processed.profile", "", basename(feature.count.fls))
  # From running "get_raw_read_counts.sh"
  input.read.count.fls <- paste0(all.base.names, "_input_readcount.txt")
  if(!all(file.exists(input.read.count.fls))){
    warning("SOME FILES DO NOT EXIST: ", input.read.count.fls[file.exists(input.read.count.fls)])
  }
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "Input Read Counts"))
  input.read.count.dt <- rbindlist(lapply(input.read.count.fls, fread, header=FALSE))
  setnames(input.read.count.dt, c("File.Base.ID", "Sample.ID", "IndexNum", "Input.ReadCount"))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  # Read in Feature Count alignments 
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "Profile Counts"))
  all.feature.dt <- rbindlist(lapply(seq_along(feature.count.fls),
                                     FUN=function(x){
                                       fl <- feature.count.fls[x]
                                       bn <- sub("_L001_R1_001_Processed.*", "", basename(fl))
                                       dt <- fread(fl, header=FALSE)
                                       dt[, File.Base.ID:=bn]
                                       dt[, V7:=NULL]
                                       return(dt)
                                     }))
  
  setnames(all.feature.dt, c("db.mismatch.strand", "SeqID.count", "Sequence", "Aln.SeqID", "PosAln.PosMM", "count",  "File.Base.ID"))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  message("###")
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis WRANGLING"))
  db.info <- all.feature.dt[, tstrsplit(db.mismatch.strand[1], split=".", fixed=TRUE), by=db.mismatch.strand]
  setnames(db.info, c("V1", "V2"), c("db", "mis.mm.strand"))
  db.info[, `:=`(mm=as.numeric(substr(mis.mm.strand, start = nchar(mis.mm.strand)-1, stop = nchar(mis.mm.strand)-1)),
                 strand=substr(mis.mm.strand, start = nchar(mis.mm.strand), stop = nchar(mis.mm.strand)))]
  db.info[, mis.mm.strand:=NULL]
  all.feature.dt[db.info, `:=`(db=i.db, strand=i.strand, mm=i.mm), on="db.mismatch.strand"]
  all.feature.dt[, IndexNum:=tstrsplit(File.Base.ID[1], split="_")[2], by=File.Base.ID]
  all.freature.dt.annot <- merge(sample.info.dt, all.feature.dt, by="IndexNum")
  
  all.freature.dt.annot[db!="UniVec", c("Aln.SeqID.nm", "Aln.Subpool.ID", "AlnSeq.length", "Modification"):=tstrsplit(Aln.SeqID[1], split="|", fixed=TRUE), by=Aln.SeqID]
  all.freature.dt.annot[univec.dt, `:=`(Aln.SeqID.nm=i.Seq.Info, Aln.Subpool.ID="Univec"), on="Aln.SeqID"]
  all.freature.dt.annot[, db.matches:=Aln.Subpool.ID==Subpool.ID]
  all.freature.dt.annot[, Idx.SeqID.count:=paste0(File.Base.ID, "-", SeqID.count)]
  all.freature.dt.annot.filt <- subset(all.freature.dt.annot, db.matches==TRUE | Aln.Subpool.ID=="Univec")
  setkey(all.freature.dt.annot.filt, Idx.SeqID.count)
  # Get the min mismatches for each sequence overall
  all.freature.dt.annot.filt[, min.mm.by.seq:=min(mm), by=.(Idx.SeqID.count)]
  # Get the min mismatches for correct db, correct strand
  all.freature.dt.annot.filt[, min.mm.by.seq.db.correct.strand:=min(ifelse(strand=="+" & db.matches==TRUE, mm, Inf)), by=.(Idx.SeqID.count)]
  # Get the min mismatches for correct db, incorrect strand
  all.freature.dt.annot.filt[, min.mm.by.seq.db.incorrect.strand:=min(ifelse(strand=="-" & db.matches==TRUE, mm, Inf)), by=.(Idx.SeqID.count)]
  # Get the min mismatches for Univec. Either Strand
  all.freature.dt.annot.filt[, min.mm.univec:=min(ifelse(Aln.Subpool.ID=="Univec", mm, Inf)), by=.(Idx.SeqID.count)]
  
  
  
  # This sets the "BestMatch.db" variable to the sample Subpool.ID if the + strand match to the correct database is the lowest. Otherwise, it's labeled "Univec" if univec is the lowest. 
  all.freature.dt.annot.filt[, BestMatch.db:=ifelse(min.mm.by.seq==min.mm.by.seq.db.correct.strand,
                                                    Subpool.ID,
                                                    ifelse(min.mm.by.seq==min.mm.univec, 
                                                           "Univec",
                                                           ifelse(min.mm.by.seq==min.mm.by.seq.db.incorrect.strand,
                                                                  "Subpool(-)",
                                                                  "OTHER")
                                                    ))]
  # Now filter to keep only the 
  all.freature.dt.annot.filt2 <- subset(all.freature.dt.annot.filt, mm==min.mm.by.seq & ( Aln.Subpool.ID==BestMatch.db | ( BestMatch.db=="Subpool(-)" & Aln.Subpool.ID == Subpool.ID &  strand=="-")  ))
  all.freature.dt.annot.filt2[, BestMatch.db2:=ifelse(Subpool.ID==BestMatch.db, "Subpool(+)", BestMatch.db)]
  all.freature.dt.annot.filt2[, NH:=.N, by=.(Idx.SeqID.count)]
  all.freature.dt.annot.filt2[, scaled.count:=count/NH]
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  message("####")
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis SUMMARY"))
  all.freature.dt.annot.summary.mm <- all.freature.dt.annot.filt2[, .(count=sum(scaled.count)), by=.(Sample_Name, Treatment, Subpool.ID, mm, BestMatch.db2)]
  all.freature.dt.annot.summary.mm[, tot:=sum(count), by=Sample_Name]
  all.freature.dt.annot.summary.mm[, perc.tot:=count/tot]
  
  all.freature.dt.annot.summary <- all.freature.dt.annot.filt2[, .(count=sum(scaled.count)), by=.(Sample_Name, IndexNum, Treatment, Subpool.ID, BestMatch.db2)]
  all.freature.dt.annot.summary[, tot:=sum(count), by=Sample_Name]
  all.freature.dt.annot.summary[, perc.tot:=count/tot]
  
  feat.count.dt <- all.freature.dt.annot.filt2[, .(tot.feat=sum(count), tot.feat.scaled=sum(scaled.count)), by=.(IndexNum, Aln.SeqID, Subpool.ID, BestMatch.db, Aln.SeqID.nm, Modification, AlnSeq.length, BestMatch.db2)]
  
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "Input Processed Read Counts"))
  # Read in counts from the Preprocessed reads
  all.processed <- rbindlist(lapply(seq_along(processed.read.fls),
                                    FUN=function(x){
                                      fl <- processed.read.fls[x]
                                      bn <- sub("_L001_R1_001_.*", "", basename(fl))
                                      dt <- fread(input = paste0("awk \'FNR%2==1\' ", fl), sep = "-")
                                      dt.sum <- dt[, .(File.Base.ID=bn, Preprocessed.count=sum(V2))]
                                      dt.sum[, c("Sample.ID", "IndexNum"):=tstrsplit(File.Base.ID, split="_")]
                                      return(dt.sum)
                                    }))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  # Read counts from reads not aligned to any database
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Import: ", "Input Unaligned Read Counts"))
  all.unaligned <- rbindlist(lapply(seq_along(unaligned.count.fls),
                                    FUN=function(x){
                                      fl <- unaligned.count.fls[x]
                                      bn <- sub("_L001_R1_001_.*", "", basename(fl))
                                      dt <- fread(input = paste0("awk \'FNR%2==1\' ", fl), sep = "-")
                                      dt.sum <- dt[, .(File.Base.ID=bn, Unaligned.count=sum(V2))]
                                      dt.sum[, c("Sample.ID", "IndexNum"):=tstrsplit(File.Base.ID, split="_")]
                                      return(dt.sum)
                                    }))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  read.count.summary <- merge(merge(input.read.count.dt, all.processed, by=c("File.Base.ID", "Sample.ID", "IndexNum")), all.unaligned, by=c("File.Base.ID", "Sample.ID", "IndexNum"))
  all.freature.dt.annot.summary.cast <- dcast.data.table(all.freature.dt.annot.summary[BestMatch.db2!="Subpool(-)"], IndexNum~BestMatch.db2, fun.aggregate = sum, value.var = "count")
  read.count.summary2 <- merge(read.count.summary, all.freature.dt.annot.summary.cast, by="IndexNum", all.x=TRUE)
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- START: ", "sRNAnalysis Finalizing Output"))
  dt.list <- list(input.files=list(profile.counts=feature.count.fls,
                                   unaligned.fa=unaligned.count.fls,
                                   processed.read.fa=processed.read.fls),
                  sample.info.dt=sample.info.dt[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")],
                  read.count.summary=read.count.summary2[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")],
                  annot.summary=all.freature.dt.annot.summary[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")],
                  annot.summary.mm=all.freature.dt.annot.summary.mm[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")],
                  feature.count.dt=feat.count.dt[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")],
                  all.feature.dt=all.freature.dt.annot.filt2[, library.method:=ifelse(filt.lib==TRUE, library.run.type, "NOT.PROVIDED")]
  )
  
  return(dt.list)
  message(paste0(Sys.time()-tstart, " -- ", biofragmenta.outdir, "-- DONE ----"))
  message(paste0("ALL DONE", biofragmenta.outdir, "--- TIME ELAPSED: ", Sys.time()-tstart))
}