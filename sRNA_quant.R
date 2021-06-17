#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#################################################################################################
#################################################################################################
#################################################################################################
## This script will count reads from a .bam sequence alignment file against the C. elegans
## genome (version WS262). The reads will be counted as sense or antisense to features in the
## genome, termed biotypes (like miRNA, protein coding genes, transposon etc.).
##
## As input, the script requires a .bam file to be counted and two reference files for the 
## genome annotations: WS262.gene.whole.transposons.repeats.RData and
## WS262.gene.parts.transposons.repeats.RData.
##  
## The output of this script will be a list with two elements: (1) results, and (2) results.size.
## results is a list of tables containing the count numbers per gene per biotype.
## results.sizes is a summary table of the amount of reads counted per biotype.
##
## To run this script, just call from the terminal on your .bam file: 
##>Rscript sRNA_quant.R file.bam 
#################################################################################################
#################################################################################################
#################################################################################################

# load libraries
suppressPackageStartupMessages({
  library(gridExtra)
  library(bitops)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(seqinr)
  library(dplyr)
  require(reshape2)
  require(data.table)
})

prepare_biotype_reference <- function(reference, biotype, sense){
  reference.biotype <- reference[reference$gene_biotype == biotype]
  if (biotype %in% biotypes.antisense) {
    if (sense == "expected") {
      reference.biotype <- invertStrand(reference.biotype)
      reference.biotype$gene_biotype <- paste(biotype, "_AS", sep="")
    }
  } else {
    if (sense == "unexpected") {
      reference.biotype <- invertStrand(reference.biotype)
      reference.biotype$gene_biotype <- paste(biotype, "_AS", sep="")
    }
  }
  # update whole.ref for parts references
  return(reference.biotype)
}

count_reads <- function(ov, bam) {
  # fix names and merge to remove ov.parts edge cases
  colnames(ov)[2] <- "gene"
  colnames(ov.parts)[2] <- "feature"
  ov.merged <- merge(ov, ov.parts, by=c("gene", "queryHits", "reads"))
  # determine which reads map to more than one feature
  ov.merged[, shared := .N, by=reads]
  # determine distribution of unique reads
  ov.merged[, prior := sum(shared == 1), by=list(gene, feature)]
  ov.merged[, sum.prior := sum(prior), by=reads]
  # partition reads between genes
  ov.merged[, fraction := prior / sum.prior ]
  ov.merged[is.nan(fraction), fraction := 1 / shared]  # can this be optimized
  # add sequence properties of reads
  ov.merged[, seq := as.character(mcols(bam[queryHits])$seq)]
  ov.merged[, cigar := as.character(mcols(bam[queryHits])$cigar)]
  ov.merged[, pos := paste(as.character(mcols(bam[queryHits])$rname), as.character(mcols(bam[queryHits])$pos), as.numeric(mcols(bam[queryHits])$flag), sep = ":")]
  # count reads
  count <- ov.merged[,.(count=sum(fraction)),by=.(gene,feature,seq,cigar,pos)]
  return(count)
}

# read in references
load("WS262.gene.parts.transposons.repeats.RData")
load("WS262.gene.whole.transposons.repeats.RData")

# prepare biotype vectors
biotypes <- c(
  "transposable_element", "repeat_region", "rRNA", "snoRNA", "snRNA", "tRNA", "piRNA",
  "miRNA", "ncRNA", "protein_coding", "pseudogene", "lincRNA"
)
biotypes.antisense <-
  c("transposable_element", "repeat_region", "ncRNA", "protein_coding", "pseudogene")

# set aside room for results
results <- list()
results.sizes <- list()

bam_file <- args[1]

  print(bam_file)
  
  #read in entire BAM file
  param <- ScanBamParam(what = c("seq", "flag", "qname", "cigar", "rname", "pos"), reverseComplement=TRUE)
  bam <- readGAlignments(bam_file, param=param)
  name <- str_extract(bam_file, "(?<=/)([^/]+)(?=/)")
  
  total.count <- length(unique(mcols(bam)$qname))

  results[[name]] <- list()
  seen <- c()

  counted.count <- 0
  print(paste("biotype", "#_features", "counts"))
  for(sense in c("expected", "unexpected")) {
    for (biotype in biotypes) {
      # prepare biotype specifc gene.whole and find overlaps
      gene.whole.biotype <- prepare_biotype_reference(gene.whole, biotype, sense)
      gene.parts.biotype <- prepare_biotype_reference(gene.parts, biotype, sense)
      biotype <- as.character(gene.whole.biotype[1]$gene_biotype)
      ov <- findOverlaps(bam, gene.whole.biotype, type = "within")
      ov.parts <- findOverlaps(bam, gene.parts.biotype)
      ov.parts <- as.data.table(ov.parts)
      ov.parts[, reads := mcols(bam[ov.parts[, queryHits]])$qname]
      
      ov.parts$gene <- gene.whole[gene.parts.biotype[ov.parts$subjectHits]$whole.ref]$parts.ref
      ov.parts$subjectHits <- gene.parts.biotype[ov.parts$subjectHits]$type
      
      # little sanity check
      stopifnot(sum(gene.parts.biotype$gene_id != gene.whole.biotype[gene.whole[gene.parts.biotype$whole.ref]$parts.ref]$gene_id) == 0)
      
      # map alignments (queryHits) to reads
      ov <- as.data.table(ov)
      ov[, reads := mcols(bam[ov[, queryHits]])$qname]
      
      # remove previously seen alignments and update seen set
      ov <- ov[!(reads %in% seen), ]
      toDistribute <- length(unique(ov[, reads]))
      potentials <- unique(ov[, subjectHits])
      seen <- union(seen, unique(ov[, reads]))  # update seen
      
      # count reads
      counted <- 0
      if (nrow(ov)) {
        counts <- count_reads(ov, bam)
        counts[, gene := gene.whole.biotype[gene]$gene_id]
        results[[name]][[biotype]] <- counts
        counted <- counts[, sum(count)]
      }
      counted.count <- counted.count + counted
      
      print(paste(biotype, length(gene.whole.biotype), counted))
      stopifnot(abs(counted - toDistribute) < 0.1)  # allow for floating point error
    }
  }
  
  # what didn't map to a feature
  results[[name]][["no_feature"]] <- data.table(
    gene = "no_feature",
    seq = as.character(mcols(bam[!(mcols(bam)$qname %in% seen)])$seq),
    cigar = as.character(mcols(bam[!(mcols(bam)$qname %in% seen)])$cigar),
    pos = paste(mcols(bam[!(mcols(bam)$qname %in% seen)])$rname, mcols(bam[!(mcols(bam)$qname %in% seen)])$pos, mcols(bam[!(mcols(bam)$qname %in% seen)])$flag, sep = ":"),
    count = 1 / as.numeric(table(mcols(bam[!(mcols(bam)$qname %in% seen)])$qname)[mcols(bam[!(mcols(bam)$qname %in% seen)])$qname])
  )[,.(count=sum(count)),by=.(gene,seq,cigar,pos)]
  
  results.sizes[[name]] <- c(
    mapped=total.count,
    counted=counted.count
  )
  rm(bam)

name <- gsub(".*/([^/]+)/.*", "\\1", bam_file, perl=TRUE)
dir.create("results")
save(list=c("results", "results.sizes"), file=paste("results/", name, ".RData", sep=""))
