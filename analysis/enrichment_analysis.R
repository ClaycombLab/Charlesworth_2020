# Here you will find the funtion used to compute the statistics for the VennPie diagrams.

##  ===========================================================
##    FISHER'S EXACT TEST (Figs. 3D, 4E, 5C, S5D, S6CDE, S8C)
##  ===========================================================

## Contingency Table Overview:
##                    | Enriched                | NOT Enriched                |
##                    | ----------------------- | --------------------------- | ------------------------
##  IN anno group     | Enriched & IN anno      | NOT Enriched & IN anno      | total genes IN anno
##  NOT IN anno group | Enriched & NOT IN anno  | NOT Enriched & NOT IN anno  | total genes NOT IN anno
##                    | ----------------------- | --------------------------- | ------------------------
##                    | total genes Enriched    | total genes NOT enriched    | total genes 

# load libraries (required by enrichment_analysis_*_reads.R)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(ggseqlogo)
  library(VennDiagram)
  library(stringr)
  library(Rmisc)
  library(pheatmap)
  library(dendsort)
  library(dendextend)
  library(tidyverse)
})

## LOCAL DEPENDECIES
##  The functions in this script take .Rds files as input and produce
##    the enrichment tables required by test.enrichment().
source(make_enrichment_tables.R)

## Read and format Ortiz data (germline)
ortiz2014 <- read.delim("Ortiz_2014/Ortiz2014_TableS1.txt")
ortiz_all <- ortiz2014[3:nrow(ortiz2014),c(8,11)]
colnames(ortiz_all) <- c("category", "ID")
ortiz_all$category <- gsub(" ", "_", ortiz_all$category)
ortiz_all$category <- mapply(paste0, "Ortiz2014_", ortiz_all$category)

## Read and format Gu data (glp-4)
gu2009 <- read.delim("Gu_2009/Gu_2009_glp-4_22Gs.txt")
gu_all <- gu2009[which(gu2009$WormBase.Gene.ID != "not found"),c(1,3)]
colnames(gu_all) <- c("category", "ID")
gu_all$category <-  gsub("glp-4_$", "glp-4_no_change", mapply(paste0, "glp-4_", gu_all$category))

## Read and format Haenni data (intestine)
haenni2012 <- read.delim("Haenni_2012/haenni_2012_intestine.txt")
haenni_all <- haenni2012[which(haenni2012$WormBase.Gene.ID != "not found"),c(1,2)]
colnames(haenni_all) <- c("category", "ID")
haenni_all$category <- gsub("Unsorted", "intestine_downregulated", haenni_all$category)
haenni_all$category <- gsub("Sorted", "intestine_upregulated", haenni_all$category)
haenni_all$category <- gsub("Other", "intestine_no_change", haenni_all$category)

ref_set <- list(Ortiz2014=ortiz_all,
                Gu2009=gu_all,
                Haenni2012=haenni_all)

input = list.files("./", pattern = ".Rds$")
names(input) = str_extract(input, "^[^(.*?)[:punct:]]+") #just takes first part of file name for later naming. Be careful if two files start with the same start (ie. CSR_56 and CSR_49) they will be overwritten.
reps_names = read.table("replicates.csv", header = TRUE, row.names = 1, sep = ",", stringsAsFactors = FALSE) #this table has 5 coulmns Argonaute, Input1, Input2, IP1, IP2 => these should hold the name of the Argonaute and the names of the .Rdata files from the counting script. If no reps, duplicate.
geneIDs_WS262 = as.data.table(read.csv(file = "../../geneIDs_WS262.csv", header = TRUE)) #a table to convert WBGene names to common gene names.
rpm_cutoff = 0
fc_cutoff = 0

test.enrichment <- function(input, ref_set, fn){
  
  ## test.enrichment() produces...
  ##    input:
  ##    ref_set: List of dataframes of references datasets to test against.
  ##                - Individual dataframes contain 2 columns: "category" (annotation) and "ID" (genes)
  
  
  rpm_cutoff = 0
  fc_cutoff = 0
  
  results <- sapply(input, fn, simplify = FALSE, USE.NAMES = TRUE)
  
  output <- matrix(nrow=0, ncol=5)
  for(i in 1:length(results)){
    
    sample <- names(results)[i]
    enriched_genes <- results[[i]][results[[i]][, IP1 >= 5 &
                                                  IP2 >= 5 &
                                                  IP1_enrichment >= 2 &
                                                  IP2_enrichment >= 2],1]$gene
    not_enriched_genes <- results[[i]][results[[i]][, !(IP1 >= 5 &
                                                          IP2 >= 5 &
                                                          IP1_enrichment >= 2 &
                                                          IP2_enrichment >= 2)],1]$gene
    
    
    for (dataset in 1:length(ref_set)){
      dset <- names(ref_set)[dataset]
      categories <- unique(ref_set[[dataset]]$category)
      
      for (category in 1:length(categories)){
        
        annotation <- categories[category]
        anno <- ref_set[[dataset]][which(ref_set[[dataset]]$category == categories[category]),2]
        ## NOTICE!
        ##  - adding "$category" in the above line fixed an issue with the antibody comparison but is not tested for published dataset comparisons!
        
        a = length(which(enriched_genes %in% anno))  ## no. enriched & IN anno (overlap enriched & anno)
        b = length(which(!(enriched_genes %in% anno)))  ## no.enriched & NOT IN anno (= enriched - a)
        c = length(which(not_enriched_genes %in% anno))   ## no. NOT enriched & IN anno (overlap enriched & anno)
        d = length(which(!(not_enriched_genes %in% anno)))  ## no. NOT enriched & NOT IN anno (= NOT enriched - c)
        counts <- (matrix(data = c(a, b, c, d), nrow = 2))
        rownames(counts) <- c("IN", "NOT_IN")
        colnames(counts) <- c("Enriched", "NOT_Enriched")
        
        enrichment <- fisher.test(counts, alternative = "greater") #enrichment
        depletion <- fisher.test(counts, alternative = "less") #depletion
        change <- fisher.test(counts, alternative = "two.sided") #different
        
        row <- c(sample, dset, annotation, enrichment$p.value, depletion$p.value)
        
        ## p-val threshold
        if (enrichment$p.value <= 0.01 | depletion$p.value <= 0.01){
          output <- rbind(output, row)
        }
      }
    }
  }
  output <- as.data.frame(output, row.names = FALSE)
  colnames(output) <- c("sample", "dataset", "annotation", "over-representation (p)", "under-representation (p)")
  return(output)
}

output_all <- test.enrichment2(input, ref_set, enrichment.all.reads)

## mRNA
sperm_proteome <- mRNA_comparisons[,1]
sperm_proteome <- sperm_proteome[sperm_proteome != ""]
sperm_proteome <- as.data.frame(cbind(rep("sperm_proteome", length(sperm_proteome)), as.character(sperm_proteome)),stringsAsFactors = F)
colnames(sperm_proteome) <- c("category", "ID")


TRAP_intestine <- mRNA_comparisons[,2]
TRAP_intestine <- TRAP_intestine[TRAP_intestine != ""]
TRAP_intestine <- as.data.frame(cbind(rep("TRAP_intestine", length(TRAP_intestine)), as.character(TRAP_intestine)),stringsAsFactors = F)
colnames(TRAP_intestine) <- c("category", "ID")

TRAP_muscle <- mRNA_comparisons[,3]
TRAP_muscle <- TRAP_muscle[TRAP_muscle != ""]
TRAP_muscle <- as.data.frame(cbind(rep("TRAP_muscle", length(TRAP_muscle)), as.character(TRAP_muscle)),stringsAsFactors = F)
colnames(TRAP_muscle) <- c("category", "ID")

TRAP_neuron <- mRNA_comparisons[,4]
TRAP_neuron <- TRAP_neuron[TRAP_neuron != ""]
TRAP_neuron <- as.data.frame(cbind(rep("TRAP_neuron", length(TRAP_neuron)), as.character(TRAP_neuron)),stringsAsFactors = F)
colnames(TRAP_neuron) <- c("category", "ID")

L4_EE <- mRNA_comparisons[,5]
L4_EE <- L4_EE[L4_EE != ""]
L4_EE <- as.data.frame(cbind(rep("L4_EE", length(L4_EE)), as.character(L4_EE)),stringsAsFactors = F)
colnames(L4_EE) <- c("category", "ID")

ref_set <- list(Ortiz2014=ortiz_all,
                Gu2009=gu_all,
                Haenni2012=haenni_all,
                sperm_proteome=sperm_proteome,
                TRAP_intestine=TRAP_intestine,
                TRAP_muscle=TRAP_muscle,
                TRAP_neuron=TRAP_neuron,
                L4_EE=L4_EE)

test.enrichment.mRNA <- function(input, ref_set){
  
  enriched_genes <- input[,1]
  enriched_genes <- enriched_genes[enriched_genes != ""]
  print(length(enriched_genes))
  not_enriched_genes <- input[,2]
  not_enriched_genes <- not_enriched_genes[not_enriched_genes != ""]
  print(length(not_enriched_genes))
  output <- matrix(nrow=0, ncol=4)
  for (dataset in 1:length(ref_set)){
    dset <- names(ref_set)[dataset]
    categories <- unique(ref_set[[dataset]]$category)
    
    for (category in 1:length(categories)){
      
      annotation <- categories[category]
      anno <- ref_set[[dataset]][which(ref_set[[dataset]] == categories[category]),2]
      
      a = length(which(enriched_genes %in% anno))  ## no. enriched & IN anno (overlap enriched & anno)
      b = length(which(!(enriched_genes %in% anno)))  ## no.enriched & NOT IN anno (= enriched - a)
      c = length(which(not_enriched_genes %in% anno))   ## no. NOT enriched & IN anno (overlap enriched & anno)
      d = length(which(!(not_enriched_genes %in% anno)))  ## no. NOT enriched & NOT IN anno (= NOT enriched - c)
      counts <- (matrix(data = c(a, b, c, d), nrow = 2))
      rownames(counts) <- c("IN", "NOT_IN")
      colnames(counts) <- c("Enriched", "NOT_Enriched")
      #print(counts)
      enrichment <- fisher.test(counts, alternative = "greater") #enrichment
      depletion <- fisher.test(counts, alternative = "less") #depletion
      change <- fisher.test(counts, alternative = "two.sided") #different
      
      row <- c(dset, annotation, enrichment$p.value, depletion$p.value)
      
      if (enrichment$p.value <= 0.01 | depletion$p.value <= 0.01){
        output <- rbind(output, row)
      }
    }
  }
  output <- as.data.frame(output, row.names = FALSE)
  colnames(output) <- c("dataset", "annotation", "over-representation (p)", "under-representation (p)")
  return(output)
}

output_mRNA_enriched_2 <- test.enrichment.mRNA(mRNA_enriched, ref_set)

