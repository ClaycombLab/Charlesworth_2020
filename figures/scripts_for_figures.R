# Here you will find the functions created to produce the different plots in the paper:
#   - Clustering diagrams
#   - Metagene plots


##  ==========================================
##    CLUSTERING DIAGRAMS (Figs.4B, 4D, S8B)
##  ==========================================

# Load required packages
library(data.table)
library(pheatmap)
library(dendsort)

# Main function
make.heatmap <- function(all_enriched, extratext){
  
  ## make.heatmap() produces...
  ##    all_enriched:
  ##    extratext: (optional) A convenience parameter to distinguish output file names.
  
  for (btype in unique(all_enriched$biotype)){
    print(btype)
    dt = dcast.data.table(all_enriched, formula = biotype + gene ~ Argonaute, value.var = "avg_enrichment")
    dt[is.na(dt)] = 0
    dt = dt[biotype==btype]
    
    if (nrow(dt) < 2){
      next
    }

    mat = as.matrix(dt[,3:ncol(dt)])
    row.names(mat) = dt$gene
    mat[mat>0] = 1

    sort.myhclust = function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
    sorted_cols = sort.myhclust(hclust(dist(t(mat))))

    n_genes <- c()
    for (i in 3:ncol(dt)){
      n <- length(which(dt[, ..i] != 0))
      n_genes <- c(n_genes, n)
    }
    ago <- sorted_cols$labels
    col_labels <- paste(ago, n_genes)

    p = pheatmap(mat, 
                 main = btype, 
                 show_rownames = FALSE, 
                 legend = FALSE, 
                 color = c("#ffffff","#2295db"), 
                 angle_col = 90, 
                 fontsize_col = 12, 
                 cluster_cols = sorted_cols, 
                 cluster_rows = TRUE,
                 treeheight_row = 0,
                 labels_col = col_labels,
    )
    
    pdf(paste0(extratext, btype,"_heatmap.pdf"), width = 11, height = 8.5)
    print(p)
    dev.off()
  }
}

make.heatmap(all_reads, extratext = "allReads_")

##  =========================================
##    METAGENE PLOTS (Figs. 3E, S5E, & S7A)
##  =========================================

#load packages
suppressPackageStartupMessages({
  library(GenomicFeatures) 
  library(metagene2)
  library(gridExtra)
  library(ggplot2)
  library(data.table)
  library(ggforce)
  library(dplyr)
})

## This script is based on the metagene2 vignette (http://www.bioconductor.org/packages/devel/bioc/vignettes/metagene2/inst/doc/metagene2.html).
## Please refer to that document for specifics on metagene2 parameters.

## Main function
make.metagene <- function(paths_to_bams, master_table, txdb, replicates=c("separate", "combined"), type=c("gene-body", "exons-only", "transcript"), extraText=""){
  
  ## make.metagene() produces metagene plots for provided RNA-seq reads mapped to a TxDB object.
  ##    paths_to_bams: A data frame providing file paths for bam files to analyze
  ##    master_table: A data frame containing the targets of interest (enriched targets for each Argonaute)
  ##    txdb: A TxDB object (either gene body, exons, or transcripts)
  ##    replicates: Either "separate" or "combined" to tell the function how replicates should be handled in the final plot.
  ##    type: A string to signify the type of TxDB object.
  ##    extraText: (optional) A convenience parameter to distinguish output file names.
  
  ## For added security, in case input file is not properly labeled
  names(master_table) <- c("X", "biotype", "gene", "avg_enrichment", "sample")
  
  ## If only considering reads mapping to exons, regions should be stiched together.
  ## Otherwise, each region is kept separate.
  if (type == "gene-body"){
    region_mode = "separate"
  } else if (type == "exons-only"){
    region_mode = "stitch"
  } else if (type == "transcript"){
    region_mode = "separate"
  }
  
  ## paths_to_bams can contain multiple experiments to be analyzed at once.
  ## Here, we iterate over all samples and generate a metagene plot for each.
  for (i in 1:nrow(paths_to_bams)){
    
    ## Extact the name of the Argonaute for sample i
    argonaute = row.names(paths_to_bams)[i]
    
    ## Subset the provided txdb object to only include those genes enriched in the sRNA profile of the Argonaute of interest.
    enriched = genes[master_table[sample == argonaute & biotype == "protein_coding_AS"]$gene]
    
    ## If only considering reads mapping to exons, additional metadata is required
    if (type == "exons-only"){region_metadata=data.frame(md = rep(argonaute,length(enriched)))} else {region_metadata=NULL}
    
    ## Printing the current Argonaute and the number of enriched targets provides a good sanity check.
    print(paste0(argonaute, " ====> ", length(enriched), " Targets"))
    
    if (replicates == "separate"){
      
      ## Expand paths of bam files, separating by replicate.
      R1_bams = c(paste0(getwd(),"/STAR/",paths_to_bams[i,]$Input1,"/", paths_to_bams[i,]$Input1, ".bam"),
                  paste0(getwd(),"/STAR/",paths_to_bams[i,]$IP1,"/", paths_to_bams[i,]$IP1, ".bam"))
      R2_bams = c(paste0(getwd(),"/STAR/",paths_to_bams[i,]$Input2,"/", paths_to_bams[i,]$Input2, ".bam"),
                  paste0(getwd(),"/STAR/",paths_to_bams[i,]$IP2,"/", paths_to_bams[i,]$IP2, ".bam"))
      
      ## Specify experimental design for each replicate
      R1_design = data.frame(samples = R1_bams, "input" = c(1,0), "IP" = c(0,1))
      R2_design = data.frame(samples = R2_bams, "input" = c(1,0), "IP" = c(0,1))
      
      ## Build metagene object (replicate 1)
      R1_mg = metagene2$new(regions = enriched,
                            bam_files = R1_bams,
                            region_mode = region_mode,
                            region_metadata = region_metadata,
                            strand_specific = FALSE,
                            padding_size = 0,
                            paired_end = FALSE,
                            invert_strand = TRUE, verbose=T)
      R1_mg$group_coverages(design = R1_design, normalization = "RPM")
      R1_mg$bin_coverages(bin_count = 50)
      if (type == "exons-only"){R1_mg$split_coverages_by_regions(split_by = "md")} else {R1_mg$split_coverages_by_regions()}
      R1_mg$calculate_ci()
      
      ## Build metagene object (replicate 2)
      R2_mg = metagene2$new(regions = enriched,
                            bam_files = R2_bams,
                            region_mode = region_mode,
                            region_metadata = region_metadata,
                            strand_specific = FALSE,
                            padding_size = 0,
                            paired_end = FALSE,
                            invert_strand = TRUE, verbose=T)
      R2_mg$group_coverages(design = R2_design, normalization = "RPM")
      R2_mg$bin_coverages(bin_count = 50)
      if (type == "exons-only"){R2_mg$split_coverages_by_regions(split_by = "md")} else {R2_mg$split_coverages_by_regions()}
      R2_mg$calculate_ci()
      
      ## Generate the plots (one for each replicate) and save to file
      R1_plot = R1_mg$plot(group_by = "design", title = paste0(argonaute," Targets Metagene Analysis R1"))
      R2_plot = R2_mg$plot(group_by = "design", title = paste0(argonaute," Targets Metagene Analysis R2"))
      ggsave(paste0(argonaute, "_", extraText, "_", type, "_metagene-analysis-separate.pdf"), arrangeGrob(R1_plot, R2_plot),
             path = paste0(getwd(), "/../"), device = "pdf",
             width = 11, height = 8.5, units = "in")
      
    } else if (replicates == "combined"){
      
      ## Expand paths of bam files
      bams <- c(paste0(getwd(),"/STAR/",paths_to_bams[i,]$Input1,"/", paths_to_bams[i,]$Input1, ".bam"),
                paste0(getwd(),"/STAR/",paths_to_bams[i,]$IP1,"/", paths_to_bams[i,]$IP1, ".bam"),
                paste0(getwd(),"/STAR/",paths_to_bams[i,]$Input2,"/", paths_to_bams[i,]$Input2, ".bam"),
                paste0(getwd(),"/STAR/",paths_to_bams[i,]$IP2,"/", paths_to_bams[i,]$IP2, ".bam"))
      
      ## Specify experimental design
      design <- data.frame(Samples = bams,
                           Input = c(1,0,1,0),
                           IP = c(0,1,0,1))
      
      ## Build metagene object
      mg = metagene2$new(regions = enriched,
                         bam_files = bams,
                         region_mode = region_mode,
                         region_metadata = region_metadata,
                         strand_specific = FALSE,
                         padding_size = 0,
                         paired_end = FALSE,
                         invert_strand = TRUE, verbose=T)
      mg$group_coverages(design = design, normalization = "RPM")
      mg$bin_coverages(bin_count = 50)
      if (type == "exons-only"){mg$split_coverages_by_regions(split_by = "md")} else {mg$split_coverages_by_regions()}
      mg$calculate_ci()
      
      ## Generate the plot and save to file
      plot = mg$plot(group_by="design", title = paste0(argonaute," Targets Metagene Analysis"))
      ggsave(paste0(argonaute, "_", extraText, "_", type, "_metagene-analysis_combined.pdf"), plot=plot,
             path = paste0(getwd(), "/../"), device = "pdf",
             width = 11, height = 8.5, units = "in")
    }
  }
}

## build txdb objects
gtf_file <- "WS262_protein_coding.gtf"
txdb = makeTxDbFromGFF(gtf_file, format="gtf")
## Gene-body
genes = genes(txdb)
## Exons
exons = exonsBy(txdb, by="gene")
exons = GRangesList(lapply(exons, reduce))
## Transcripts
transcripts = transcriptsBy(txdb, by = "gene")
transcripts = GRangesList(lapply(transcripts, reduce)) ## Necessary???

## AGO IP enrichment data
paths_to_bams = read.csv(file = "paths_to_bams-all.csv", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

## Running the function!!
## Example:
##  - CSR-1[total], CSR-1a, and CSR-1b
##  - All sRNA read sizes.
##  - Within gene body.
##  - Replicates as separate plots
##  - extra text to distinguish the output file name.
make.metagene(paths_to_bams, master_table_all, genes, replicates="separate", type="gene-body", extraText="all-reads")


