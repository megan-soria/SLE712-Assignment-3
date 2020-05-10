library("seqinr")
library("R.utils")
library("rBLAST")
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")


# Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. 
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=FALSE)

# Use the makeblast() function to create a blast database. 
rBLAST::makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl", "-parse_seqids")

#How many sequences are present in the E.coli set?

## makeblastdb() console result
# Building a new DB, current time: 05/10/2020 23:30:48
# New DB name:   /mnt/student11/projects/SLE712-Assignment-3/Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa
# New DB title:  Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 4140 sequences in 0.330607 seconds.

# Download the sample fasta sequences and read them in as above. 
# For your allocated sequence, determine the length (in bp) and the proportion of GC bases.
sample_fasta <- seqinr::read.fasta("part2_ecoliGene_sample.fa")
seqinr::getLength(sample_fasta$`1`)
# 615
seqinr::GC(sample_fasta$`1`)
# 0.5560976

# You will be provided with R functions to create BLAST databases and perform blast searches. 
# Use blast to identify what E. coli gene your sequence matches best. Show a table of the top 
# 3 hits including percent identity, E-value and bit scores.

myblastn_tab # provided R function

results <- myblastn_tab(myseq = sample_fasta, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
top3_hits <- results[1:3,]
top3_hits

# You will be provided with a function that enables you to make a set number of point mutations 
# to your sequence of interest. Run the function and write an R code to check the number of 
# mismatches between the original and mutated sequence.


# Using the provided functions for mutating and BLASTing a sequence, determine the number 
# and proportion of sites that need to be altered to prevent the BLAST search from matching the 
# gene of origin. Because the mutation is random, you may need to run this test multiple times 
# to get a reliable answer.

# Provide a chart or table that shows how the increasing proportion of mutated bases reduces 
# the ability for BLAST to match the gene of origin. Summarise the results in 1 to 2 sentences.
