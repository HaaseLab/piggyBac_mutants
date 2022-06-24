findPeaksInPERegions <- function(IN.GR=NULL,
                                 GR.NAME=NULL,
                                 MIN.READS=5,
                                 PPM.FACTOR=1000000){
  
  ## GOAL: find simple stranded peaks in Genomic Ranges
  ## Author: Pavol Genzor
  ## 04.20.22; Version 1
  
  ## Libraries
  suppressPackageStartupMessages({library(data.table); library(IRanges); 
    library(GenomicRanges); library(GenomicAlignments)})
  
  ## Input checking
  if(is.null(IN.GR)) stop("Please provide IN.GR!")
  if(is.null(GR.NAME)) stop("Please provide GR.NAME!")
  if(!"NH1" %in% colnames(mcols(IN.GR)) | !"NH2" %in% colnames(mcols(IN.GR))) stop("This is not PE data. You need NH1 and NH2 columns!")
  
  message("Looking in sample:")
  message(paste0("  ",GR.NAME))
  
  message("\tusing unique anchored")
  NGR <- IN.GR[mcols(IN.GR)[["NH1"]] %in% 1 | mcols(IN.GR)[["NH2"]] %in% 1]
  total_regions <- length(NGR)
  
  message("\tcalculating stranded coverage")
  COV.PLUS <- GenomicRanges::coverage(NGR[strand(NGR) %in% "+"])
  COV.MINUS <- GenomicRanges::coverage(NGR[strand(NGR) %in% "-"])
  
  message("\tfinding peaks")
  message(paste0("\t\tmin. read: ", MIN.READS))
  SLICE.PLUS <- IRanges::slice(x = COV.PLUS, lower = MIN.READS)
  SLICE.MINUS <- IRanges::slice(x = COV.MINUS, lower = MIN.READS)
  
  message("\tmaking range with strand")
  PLUS.GR <- GenomicRanges::reduce(as(SLICE.PLUS, "GRanges"))
  MINUS.GR <- GenomicRanges::reduce(as(SLICE.MINUS, "GRanges"))
  strand(PLUS.GR) <- "+";  strand(MINUS.GR) <- "-"
  
  # adding total sequences
  mcols(PLUS.GR)[["total_regions"]] <- total_regions
  mcols(MINUS.GR)[["total_regions"]] <- total_regions
  
  message("\tcounting sequences")
  mcols(PLUS.GR)[["uniq_region_count"]] <- GenomicRanges::countOverlaps(query = PLUS.GR, subject = NGR)
  mcols(MINUS.GR)[["uniq_region_count"]] <- GenomicRanges::countOverlaps(query = MINUS.GR, subject = NGR)
  
  message("\tadding peak info")
  message(paste0("\t\tmaximum peak height"))
  mcols(PLUS.GR)[["peak_height"]] <- unlist(rbindlist(lapply(SLICE.PLUS,function(x){as.data.table(max(x))})))
  mcols(MINUS.GR)[["peak_height"]] <- unlist(rbindlist(lapply(SLICE.MINUS,function(x){as.data.table(max(x))})))
  
  message(paste0("\t\traw area under the curve"))
  mcols(PLUS.GR)[["peak_auc"]] <- unlist(rbindlist(lapply(SLICE.PLUS,function(x){as.data.table(sum(x))})))
  mcols(MINUS.GR)[["peak_auc"]] <- unlist(rbindlist(lapply(SLICE.MINUS,function(x){as.data.table(sum(x))})))
  
  message(paste0("\t\tnormalized area under the curve"))
  mcols(PLUS.GR)[["peak_nauc"]] <- as.numeric(mcols(PLUS.GR)[["peak_auc"]] / (total_regions/PPM.FACTOR))
  mcols(MINUS.GR)[["peak_nauc"]] <- as.numeric(mcols(MINUS.GR)[["peak_auc"]] / (total_regions/PPM.FACTOR))
  
  message(paste0("\t\tpeak width"))
  mcols(PLUS.GR)[["peak_width"]] <- width(PLUS.GR); mcols(PLUS.GR)[["sample"]] <- GR.NAME
  mcols(MINUS.GR)[["peak_width"]] <- width(MINUS.GR); mcols(MINUS.GR)[["sample"]] <- GR.NAME
  
  message("Finished.")
  aPEAKGR <- c(PLUS.GR,MINUS.GR)
  aPEAKGR <- sort.GenomicRanges(aPEAKGR)
  return(aPEAKGR)
}