filterBamPE <- function(
  ## INPUTS
  BAM.FILE=NULL,
  BAM.NAME=NULL,
  STRAND.MODE=2,
  BS.SPECIES=NULL,
  
  ## FILTERS
  STANDARD.CONTIGS.ONLY=TRUE, 
  LIBRARY.MAX.WIDTH=1000,
  QUANTILE.WIDTH.FILTER=TRUE,
  WIDTH.QUANTILE.PROBS=0.999,
  
  ## OPTIONS
  TAGS=c("NH","NM"),
  SIMPLE.CIGAR=TRUE,
  IS.PAIRED = TRUE, 
  IS.DUPLICATE = FALSE,
  IS.PROPER.PAIR = TRUE, 
  IS.UNMAPPED.QUERY = FALSE,
  IS.SECONDARY.ALIGNMENT = FALSE,
  
  ## OPTIONS 2
  YIELD.SIZE=NA
){
  
  ## AUTHOR: Pavol Genzor
  ## Use: Load PE bam file into R
  ## 10.25.21; Version 3; adding tags, correcting the lane info
  
  ## NOTE ON BAM INFO
  ## 10.25.21; Version 2; added Yield.size and width filter
  ## 10.22.21; Version 1; original
  ## 
  
  ## libraries
  suppressPackageStartupMessages({library("data.table");library("dplyr");library("Rsamtools");
    library("GenomicAlignments");library("BSgenome.Hsapiens.UCSC.hg38");
    library("BSgenome.Dmelanogaster.UCSC.dm6"); library("BSgenome.Mmusculus.UCSC.mm10")})
  
  ## check input
  if(is.null(BAM.FILE)) stop("Please provide full path to a .bam file !!!")
  if(is.null(BAM.NAME)) stop("Please provide .bam name !!!")
  if(is.null(BS.SPECIES)) stop("Please provide BSSPECIES name !!!")
  
  message("Starting")
  message(paste0(" file: ",BAM.NAME))
  PROGRESS.L <- list()
  PROGRESS.L[["ID"]] <- BAM.NAME
  
  if(isTRUE(STANDARD.CONTIGS.ONLY)){
    WHICH <- keepStandardChromosomes(seqinfo(eval(parse(text = BS.SPECIES)))) 
  } else { WHICH = ""}
  
  message("\tloading parameters")
  message(paste0(
    "\t\t TAGS: ", paste(TAGS,collapse = ", "),"\n",
    "\t\t STRAND.MODE: ",STRAND.MODE,"\n",
    "\t\t SIMPLE.CIGAR: ",SIMPLE.CIGAR,"\n",
    "\t\t IS.PAIRED: ",IS.PAIRED,"\n",
    "\t\t IS.PROPER.PAIR: ",IS.PROPER.PAIR,"\n",
    "\t\t IS.DUPLICATE: ",IS.DUPLICATE,"\n",
    "\t\t IS.SECONDARY.ALIGNMENT: ",IS.SECONDARY.ALIGNMENT,"\n"))
  
  PARAM <- ScanBamParam(which = WHICH, 
                        simpleCigar = SIMPLE.CIGAR,
                        tag = TAGS,
                        flag = scanBamFlag(isPaired = IS.PAIRED, 
                                           isDuplicate = IS.DUPLICATE,
                                           isProperPair = IS.PROPER.PAIR, 
                                           isUnmappedQuery = IS.UNMAPPED.QUERY,
                                           isSecondaryAlignment = IS.SECONDARY.ALIGNMENT))
  
  message("\tloading PE .bam file into GAlignments")
  if(is.na(YIELD.SIZE)){message("\t\t all reads")}
  else{message(paste0("\t\t YIELD.SIZE: ",YIELD.SIZE))}
  
  PGR <- readGAlignmentPairs(file = BamFile(file = BAM.FILE, 
                                            yieldSize = YIELD.SIZE),
                             use.names = TRUE, 
                             strandMode = STRAND.MODE, 
                             param = PARAM)
  PROGRESS.L[["INPUT"]] <- length(PGR)
  
  message("\textracting First read")
  FIRST.DT <- as.data.table(first(PGR))
  FIRST.DT[,"TAG":=paste("W",width,"NH",NH,"NM",NM,sep = ":")]

  message("\textracting Second read")
  SECOND.DT <- as.data.table(second(PGR))
  SECOND.DT[,"TAG":=paste("W",width,"NH",NH,"NM",NM,sep = ":")]
  
  message("\tcombining First_Second")
  FIRST.SECOND <- paste(FIRST.DT[["TAG"]],SECOND.DT[["TAG"]],sep=":")
  
  message("\tmaking GR")
  PGR.GR <- granges(PGR)
  in.gr <- length(PGR.GR)
  
  message("\tadding tag info")
  mcols(PGR.GR)[["TAG"]] <- FIRST.SECOND
  
  if(!is.na(LIBRARY.MAX.WIDTH)){
    message("\tfiltering by width")
    message("\t\tmax width: ",LIBRARY.MAX.WIDTH)
    PGR.GR <- PGR.GR[width(PGR.GR) <= LIBRARY.MAX.WIDTH]
    PROGRESS.L[["LIB_MAX_WIDTH"]] <- length(PGR.GR)}
  
  if(isTRUE(QUANTILE.WIDTH.FILTER)){
    message("\tfiltering by quantile")
    message(paste0("\t\tquantile: ",max(WIDTH.QUANTILE.PROBS)))
    WIDTH.QUANTILE <- quantile(x = width(PGR.GR), probs = WIDTH.QUANTILE.PROBS)
    message(paste0("\t\tquantile width: ",max(WIDTH.QUANTILE)))
    PGR.GR <- PGR.GR[width(PGR.GR) <= max(WIDTH.QUANTILE)] 
    PROGRESS.L[["QUANTILE_WIDTH"]] <- length(PGR.GR)} 
  
  else { message("\tno quantile filter") }
  
  message("\tconverting to DT")
  PGR.DT <- as.data.table(PGR.GR)
  PGR.DT[,"new_name" := paste(seqnames,start,end,strand,sep = ":")]
  PGR.DT <- setDT(PGR.DT)
  
  message("\tinterpreting TAGs")
  PGR.DT[,"W1":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 2, type.convert = TRUE)]
  PGR.DT[,"W2":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 8, type.convert = TRUE)]
  PGR.DT[,"NH1":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 4, type.convert = TRUE)]
  PGR.DT[,"NH2":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 10, type.convert = TRUE)]
  PGR.DT[,"NM1":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 6, type.convert = TRUE)]
  PGR.DT[,"NM2":= tstrsplit(TAG,split=":",fixed=TRUE, keep = 12, type.convert = TRUE)]

  ## NOTE: This is not a real multiplicity - rather tells if two fragments occupy same location
  message("\tmaking new DT")
  PGR.DT.MULT <- PGR.DT[,.N,by="new_name"]
  PGR.DT.MULT[,"chr" := tstrsplit(new_name,split=":",keep = 1)]
  PGR.DT.MULT[,"start" := tstrsplit(new_name,split=":",keep = 2)]
  PGR.DT.MULT[,"end" := tstrsplit(new_name,split=":",keep = 3)]
  PGR.DT.MULT[,"strand" := tstrsplit(new_name,split=":",keep = 4)]
  colnames(PGR.DT.MULT) <- gsub("N","MULT",colnames(PGR.DT.MULT))

  message("\tsummarizing TAGs")
  PGR.DT.W <- PGR.DT[,lapply(.SD,min),by="new_name",.SDcols = c("W1","W2")]
  PGR.DT.NH <- PGR.DT[,lapply(.SD,max),by="new_name",.SDcols = c("NH1","NH2")]
  PGR.DT.NM <- PGR.DT[,lapply(.SD,function(s){sum(s)/length(s)}),by="new_name",.SDcols = c("NM1","NM2")]
  
  message("\tcombining DTs")
  PGR.DT.M.W <- PGR.DT.MULT[PGR.DT.W,on=c("new_name")]
  PGR.DT.M.W.NH <- PGR.DT.M.W[PGR.DT.NH,on=c("new_name")]
  PGR.DT.M.W.NH.NM <- PGR.DT.M.W.NH[PGR.DT.NM,on=c("new_name")]
  PGR.DT.M.W.NH.NM[["new_name"]] <- NULL

  message("\tmaking new GR")
  PGR.NGR <- makeGRangesFromDataFrame(df = PGR.DT.M.W.NH.NM, keep.extra.columns = TRUE)
  
  message("\treporting")
  PROGRESS.L[["FINAL UNIQUE FRAGMENTS"]] <- length(PGR.NGR)
  PROGRESS.DT <- melt.data.table(as.data.table(PROGRESS.L), id.vars = "ID", value.name = "count", variable.name = "step")
  PROGRESS.DT[["percent"]] <- (PROGRESS.DT[["count"]]/PROGRESS.DT[step %in% "INPUT"][["count"]])*100
  write.csv(file=paste0(getwd(),"/",BAM.NAME,"_filterBamPE.tab"),x = PROGRESS.DT, quote = FALSE, row.names = FALSE)
  
  message("Done!")
  message("")
  print(PROGRESS.DT)
  return(PGR.NGR)
}
