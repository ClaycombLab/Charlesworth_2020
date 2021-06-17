#################################################################################################
#################################################################################################
#################################################################################################
## This script will prepare reference files to be used with the sRNA_quant.R script. This script
## only needs to be run once to generate the reference files. The output will be two files:
## WS262.gene.whole.transposons.repeats.RData and WS262.gene.parts.transposons.repeats.RData.
## These files contain data related to genomic features in the C. elegams genome (version WS262).
## It uses wormbase annotation (WS262) for all genomic features except for miRNAs which uses
## miRBase annotations (version 22.1).
#################################################################################################
#################################################################################################
#################################################################################################

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicAlignments)
  library(data.table)
  library(BiocParallel)
  library(bitops)
})
# Prepare combined reference at gene level
# read in genes file
genes <- "c_elegans.PRJNA13758.WS262.canonical_geneset.gtf"
genes <- import(genes)
genes <- genes[genes$type == "gene" & grepl("WBGene", genes$gene_id) & !grepl("miRNA", genes$gene_biotype) & !grepl("antisense", genes$gene_biotype)] #get all genes excluding miRNAs
mcols(genes) <- mcols(genes)[, c("gene_id", "gene_biotype")]

# read in gfp 
gfp <- "./3xFLAG_loxP_GFP/3xFLAG_loxP_GFP.gtf"
gfp <- import(gfp)
gfp <- gfp[gfp$type == "gene"]
mcols(gfp) <- mcols(gfp)[, c("gene_id", "gene_biotype")]

# read in miRNA data from miRBase
miRNAs <- "cel_mirbase_22.1.gff3"
miRNAs <- import(miRNAs)
newStyle <- mapSeqlevels(seqlevels(miRNAs), "NCBI") #change chromosome names to NCBI convention
miRNAs <- renameSeqlevels(miRNAs, newStyle)

miRNAs <- miRNAs[miRNAs$type == "miRNA"]
mcols(miRNAs) <- mcols(miRNAs)[, c("Name", "type")]
colnames(mcols(miRNAs)) <- c("gene_id","gene_biotype")

# readin transposable elements from transposons
transposons <- "c_elegans.PRJNA13758.WS262.annotations.transposable_element.gff3"
transposons <- import(transposons)
transposons <- transposons[transposons$type == "transposable_element"]
transposons$gene_id <- transposons$Family
transposons$gene_id[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"] <- 
  transposons$Name[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"]
transposons$gene_biotype <- transposons$type
mcols(transposons) <- mcols(transposons)[, c("gene_id", "gene_biotype")]

# readin repeat elements from transposons
repeats <- "c_elegans.PRJNA13758.WS262.annotations.repeat.gff3"
repeats <- import(repeats)
repeats <- repeats[grep("(LTR)|(Ce\\d+)", repeats$Target)]
repeats$gene_id <- gsub("^(\\w+).*", "\\1", repeats$Target)
repeats$gene_biotype <- repeats$type
mcols(repeats) <- mcols(repeats)[, c("gene_id", "gene_biotype")]

gene.whole <- c(genes, gfp, miRNAs, transposons, repeats)
gene.whole$parts.ref <- data.table(gene_biotype=gene.whole$gene_biotype)[, count := seq(.N), by=gene_biotype]$count

save(gene.whole, file="WS262.gene.whole.transposons.repeats.RData")

# read in canonical geneset and collapse transcripts at gene level into
# disjoint set with type of each element set to the union of overlapping
# features
gtf_file <- "c_elegans.PRJNA13758.WS262.canonical_geneset.gtf"
gtf <- import(gtf_file)
gtf <- gtf[!grepl("miRNA", gtf$gene_biotype)]#exclude miRNAs
gtf <- gtf[!grepl("antisense", gtf$gene_biotype)]

gfp <- "./3xFLAG_loxP_GFP/3xFLAG_loxP_GFP.gtf"
gfp <- import(gfp)
gfp$protein_id <- NA

gtf <- c(gtf, gfp)

gtf.genes <- bplapply(unique(grep("(WBGene)|(GFP)", gtf$gene_id, value=TRUE)), function(gene){
  gene.parts <- gtf[
    gtf$type %in% c("exon", "three_prime_utr", "five_prime_utr") & 
    (grepl("WBGene", gtf$gene_id) | grepl("GFP", gtf$gene_id)) &
    gtf$gene_id == gene
  ]
  # convert string to binary encoding (introns=8, full_feature=16)
  gene.parts$type <- as.integer(c("five_prime_utr"=1, "exon"=2, "three_prime_utr"=4)[as.character(gene.parts$type)])
  gene.biotype <- gtf[gtf$type == "gene" & gtf$gene_id == gene]$gene_biotype
  parts.by.transcripts <- split(gene.parts, gene.parts$transcript_id)
  parts.by.transcripts.disjoint <- lapply(parts.by.transcripts, disjoin, with.revmap=TRUE)
  
  for (name in names(parts.by.transcripts.disjoint)) {
    parts.by.transcripts.disjoint[[name]]$type <- sapply(parts.by.transcripts.disjoint[[name]]$revmap, function(x){
      # NaN because UTRs should never overlap
      c(1, 2, 1, 4, NaN, 4, NaN)[Reduce(bitOr, parts.by.transcripts[[name]][x]$type)]
    })
    parts.by.transcripts.disjoint[[name]]$revmap <- NULL
    introns <- setdiff(range(parts.by.transcripts.disjoint[[name]]), parts.by.transcripts.disjoint[[name]])
    introns$type <- rep(as.integer(8), length(introns))
    parts.by.transcripts.disjoint[[name]] <- c(
      parts.by.transcripts.disjoint[[name]],
      introns
    )
  }
  
  parts.by.transcripts.disjoint <- unlist(GRangesList(parts.by.transcripts.disjoint))
  parts.disjoint <- disjoin(parts.by.transcripts.disjoint, with.revmap=TRUE)
  parts.disjoint$type <- sapply(parts.disjoint$revmap, function(x){
    Reduce(bitOr, parts.by.transcripts.disjoint[x]$type)
  })
  parts.disjoint$revmap <- NULL
  parts.disjoint$gene_id <- gene
  parts.disjoint$gene_biotype <- gene.biotype
  parts.disjoint
})
gtf.genes <- unlist(GRangesList(gtf.genes))
gtf.genes$whole.ref <- structure(1:(length(genes)+1), names=c(genes$gene_id, "GFP"))[gtf.genes$gene_id]
save(gtf.genes, file="WS262.gene.parts.RData")

# read in miRNAs
miRNAs <- "cel_mirbase_22.1.gff3"
miRNAs <- import(miRNAs)
newStyle <- mapSeqlevels(seqlevels(miRNAs), "NCBI") #change chromosome names to NCBI convention
miRNAs <- renameSeqlevels(miRNAs, newStyle)

miRNAs <- miRNAs[miRNAs$type == "miRNA"]
mcols(miRNAs) <- mcols(miRNAs)[, c("Name", "type")]
colnames(mcols(miRNAs)) <- c("gene_id","gene_biotype")

miRNAs$type <- 16
mcols(miRNAs) <- mcols(miRNAs)[, c("type", "gene_id", "gene_biotype")]
miRNAs$whole.ref <- length(genes) + 1 + 1:length(miRNAs)
save(miRNAs, file="mirbase_22.1.miRNAs.RData")

# read in transposable elements, use entire transposable element as desired feature
transposons <- "c_elegans.PRJNA13758.WS262.annotations.transposable_element.gff3"
transposons <- import(transposons)
transposons <- transposons[transposons$type == "transposable_element"]
transposons$gene_id <- transposons$Family
transposons$gene_id[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"] <- 
  transposons$Name[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"]
transposons$gene_biotype <- transposons$type
mcols(transposons) <- mcols(transposons)[, c("type", "gene_id", "gene_biotype")]
transposons$type <- as.integer(16)  # represents full feature in current encoding scheme
transposons$whole.ref <- length(genes) + length(miRNAs) + 1 + 1:length(transposons)
save(transposons, file="WS262.transposons.RData")

# readin repeat elements from gff3
repeats <- "c_elegans.PRJNA13758.WS262.annotations.repeat.gff3"
repeats <- import(repeats)
repeats <- repeats[grep("(LTR)|(Ce\\d+)", repeats$Target)]
repeats$gene_id <- gsub("^(\\w+).*", "\\1", repeats$Target)
repeats$gene_biotype <- repeats$type
mcols(repeats) <- mcols(repeats)[, c("type", "gene_id", "gene_biotype")]
repeats$type <- as.integer(16)  # represents full gene in current encoding scheme
repeats$whole.ref <- length(genes) + 1 + length(miRNAs)+ length(transposons) + 1:length(repeats)
save(repeats, file="WS262.repeats.RData")

# combine gene parts and transposon annotations and save
gene.parts <- c(gtf.genes, miRNAs, transposons, repeats)
stopifnot(sum(gene.parts$gene_id != gene.whole[gene.parts$whole.ref]$gene_id) == 0)
save(gene.parts, file="WS262.gene.parts.transposons.repeats.RData")

