# Two isoforms of the essential <I>C. elegans</I> Argonaute CSR-1 differentially regulate sperm and oocyte fertility
Amanda G. Charlesworth<sup>1,#</sup>, Uri Seroussi<sup>1,#</sup>, Nicolas J. Lehrbach<sup>2,3,#</sup>,  Mathias S. Renaud<sup>1</sup>, Adam E. Sundby<sup>1</sup>, Ruxandra I. Molnar<sup>1</sup>, Robert X. Lao<sup>1</sup>, Alexandra R. Willis<sup>1</sup>, Jenna R. Woock<sup>1</sup>, Matthew J. Aber<sup>2,3</sup>, Annette J. Diao<sup>1</sup>, Aaron W. Reinke<sup>1</sup>, Gary Ruvkun<sup>2,3</sup>, Julie M. Claycomb<sup>1*</sup>

<sup>1</sup>Department of Molecular Genetics, University of Toronto Toronto, ON M6G2Z7 Canada 

<sup>2</sup>Department of Molecular Biology, Massachusetts General Hospital, Boston, MA 02114, USA 

<sup>3</sup>Department of Genetics, Harvard Medical School, Boston, MA 02115, USA

<sup>#</sup>Equal contribution

<sup>*</sup>Correspondence: julie.claycomb@utoronto.ca

### Read the published paper [here](https://academic.oup.com/nar/article/49/15/8836/6331683?login=true).
### All data are available from GEO accession [GSE154678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154678).

## Summary of custom scripts
### sRNA_quant.R

This script will count reads from a .bam sequence alignment file against the C. elegans
genome (version WS262). The reads will be counted as sense or antisense to features in the
genome, termed biotypes (like miRNA, protein coding genes, transposon etc.).

As input, the script requires a .bam file to be counted and two custom reference files for the 
genome annotations (see below for generating refrence files): WS262.gene.whole.transposons.repeats.RData and
WS262.gene.parts.transposons.repeats.RData.
  
The output of this script will be a list with two elements: (1) results, and (2) results.size.
results is a list of tables containing the count numbers per gene per biotype.
results.sizes is a summary table of the amount of reads counted per biotype.

To run this script, just call from the terminal on your .bam file: 
`>Rscript sRNA_quant.R file.bam`

### c_elegans_reference_prep.R

This script will prepare reference files to be used with the sRNA_quant.R script. This script
only needs to be run once to generate the reference files. The output will be two files:
WS262.gene.whole.transposons.repeats.RData and WS262.gene.parts.transposons.repeats.RData.
These files contain data related to genomic features in the C. elegans genome (version WS262).
It uses wormbase annotation (WS262) for all genomic features except for miRNAs which uses
miRBase annotations (version 22.1). Contact us for the specific input data used in this paper.
