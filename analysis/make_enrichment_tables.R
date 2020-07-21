#
#   Input: results_tables.Rds from enrichment_analysis_input_results_table.R, geneIDs_WS262.csv, & replicates1.csv
#   Output: table of enrichment values

# load libraries
# suppressPackageStartupMessages({
#   library(data.table)
#   library(dplyr)
#   library(reshape2)
#   library(ggplot2)
#   library(ggseqlogo)
#   library(VennDiagram)
#   library(stringr)
#   library(Rmisc)
#   library(pheatmap)
#   library(dendsort)
#   library(dendextend)
#   library(tidyverse)
# })

##  ===================================
##    ENRICHMENT ANALYSIS - ALL READS
##  ===================================

enrichment.all.reads <- function(x){
  
  results.table = readRDS(x)
  print(head(results.table))  
  
  ago = row.names(reps_names[names(x),])
  current_row = reps_names[names(x),]
  
  print(current_row)
  print(ago)
  
  enrichment_table = results.table[,.(RPM=sum(RPM)),by=.(gene,biotype,sample)]
  enrichment_table = as.data.table(dcast(enrichment_table, gene + biotype ~ sample, value.var = "RPM"))
  
  print(head(enrichment_table))
  
  ### THIS PART OF THE CODE IS SUPER MESSY BUT IT WORKS
  #converts all col names to Input1 IP1 etc
  ##NOTE if you only have two reps - comment out sections for column 5&6
  a = colnames(enrichment_table[,3]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,3]))
  }
  ###
  a = colnames(enrichment_table[,4]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,4]))
  }
  # ###
  a = colnames(enrichment_table[,5]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,5]))
  }
  ###
  a = colnames(enrichment_table[,6]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,6]))
  }
  
  print(enrichment_table)
  
  enrichment_table[,Argonaute:=ago]
  if (length(unique(results.table$sample)) == 2){
    enrichment_table[,Input2:=enrichment_table$Input1]
    enrichment_table[,IP2:=enrichment_table$IP1]
    print(head(enrichment_table))
  }
  
  enrichment_table[,c("Input1","Input2","IP1","IP2")] = enrichment_table[,c("Input1","Input2","IP1","IP2")] + 0.01
  enrichment_table[is.na(enrichment_table)] = 0.01
  enrichment_table[,avg_input := (Input1+Input2)/2]
  enrichment_table[,avg_IP := (IP1+IP2)/2]
  enrichment_table[,IP1_enrichment:=IP1/Input1]
  enrichment_table[,IP2_enrichment:=IP2/Input2]
  enrichment_table[,avg_enrichment := (IP1_enrichment+IP2_enrichment)/2]
  
  enriched = enrichment_table[enrichment_table[, IP1 >= rpm_cutoff &
                                                 IP2 >= rpm_cutoff &
                                                 IP1_enrichment >= fc_cutoff &
                                                 IP2_enrichment >= fc_cutoff]]
  
}


##  ===================================
##    ENRICHMENT ANALYSIS - 22G READS
##  ===================================

enrichment.22Gs <- function(x){
  
  results.table = readRDS(x)
  print(head(results.table))  
  
  ago = row.names(reps_names[names(x),])
  current_row = reps_names[names(x),]
  
  enrichment_table = results.table[matched_read_length==21 | matched_read_length==22 | matched_read_length==23]
  enrichment_table = enrichment_table[,.(RPM=sum(RPM)),by=.(gene,biotype,sample)]
  enrichment_table = as.data.table(dcast(enrichment_table, gene + biotype ~ sample, value.var = "RPM"))
  
  ### THIS PART OF THE CODE IS SUPER MESSY BUT IT WORKS
  #converts all col names to Input1 IP1 etc
  a = colnames(enrichment_table[,3]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,3]))
  }
  ###
  a = colnames(enrichment_table[,4]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,4]))
  }
  ###
  a = colnames(enrichment_table[,5]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,5]))
  } 
  ###
  a = colnames(enrichment_table[,6]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,6]))
  }
  
  enrichment_table[,Argonaute:=ago]
  if (length(unique(results.table$sample)) == 2){
    enrichment_table[,Input2:=enrichment_table$Input1]
    enrichment_table[,IP2:=enrichment_table$IP1]
  }
  
  enrichment_table[,c("Input1","Input2","IP1","IP2")] = enrichment_table[,c("Input1","Input2","IP1","IP2")] + 0.01
  enrichment_table[is.na(enrichment_table)] = 0.01
  enrichment_table[,avg_input := (Input1+Input2)/2]
  enrichment_table[,avg_IP := (IP1+IP2)/2]
  enrichment_table[,IP1_enrichment:=IP1/Input1]
  enrichment_table[,IP2_enrichment:=IP2/Input2]
  enrichment_table[,avg_enrichment := (IP1_enrichment+IP2_enrichment)/2]
  
  enriched = enrichment_table[enrichment_table[, IP1 >= rpm_cutoff &
                                                 IP2 >= rpm_cutoff &
                                                 IP1_enrichment >= fc_cutoff &
                                                 IP2_enrichment >= fc_cutoff]]
}


##  ===================================
##    ENRICHMENT ANALYSIS - 26G READS
##  ===================================

enrichment.26Gs <- function(x){
  
  results.table = readRDS(x)
  print(head(results.table))  
  
  ago = row.names(reps_names[names(x),])
  current_row = reps_names[names(x),]
  
  enrichment_table = results.table[matched_read_length==25 | matched_read_length==26 | matched_read_length==27]
  enrichment_table = enrichment_table[,.(RPM=sum(RPM)),by=.(gene,biotype,sample)]
  enrichment_table = as.data.table(dcast(enrichment_table, gene + biotype ~ sample, value.var = "RPM"))
  
  ### THIS PART OF THE CODE IS SUPER MESSY BUT IT WORKS
  #converts all col names to Input1 IP1 etc
  a = colnames(enrichment_table[,3]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,3]))
  }
  a = colnames(enrichment_table[,3]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,3]))
  }
  ###
  a = colnames(enrichment_table[,4]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,4]))
  }
  a = colnames(enrichment_table[,4]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,4]))
  }
  ###
  a = colnames(enrichment_table[,5]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,5]))
  }
  a = colnames(enrichment_table[,5]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,5]))
  } 
  ###
  a = colnames(enrichment_table[,6]) %in% reps_names$Input1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$Input2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "Input2" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP1
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP1" = colnames(enrichment_table[,6]))
  }
  a = colnames(enrichment_table[,6]) %in% reps_names$IP2
  if (a ==TRUE){
    enrichment_table<- dplyr::rename(enrichment_table, "IP2" = colnames(enrichment_table[,6]))
  } 

  enrichment_table[,Argonaute:=ago]
  if (length(unique(results.table$sample)) == 2){
    enrichment_table[,Input2:=enrichment_table$Input1]
    enrichment_table[,IP2:=enrichment_table$IP1]
  }
  
  enrichment_table[,c("Input1","Input2","IP1","IP2")] = enrichment_table[,c("Input1","Input2","IP1","IP2")] + 0.01
  enrichment_table[is.na(enrichment_table)] = 0.01
  enrichment_table[,avg_input := (Input1+Input2)/2]
  enrichment_table[,avg_IP := (IP1+IP2)/2]
  enrichment_table[,IP1_enrichment:=IP1/Input1]
  enrichment_table[,IP2_enrichment:=IP2/Input2]
  enrichment_table[,avg_enrichment := (IP1_enrichment+IP2_enrichment)/2]
  
  enriched = enrichment_table[enrichment_table[, IP1 >= rpm_cutoff &
                                                 IP2 >= rpm_cutoff &
                                                 IP1_enrichment >= fc_cutoff &
                                                 IP2_enrichment >= fc_cutoff]]
  
}
