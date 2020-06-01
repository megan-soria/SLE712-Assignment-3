library("seqinr")
library("R.utils")
library("rBLAST")
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")


# 1. Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. 
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-47/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", 
              destfile = "Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
R.utils::gunzip("Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=FALSE)

#    Use the makeblast() function to create a blast database. 
rBLAST::makeblastdb("Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl", "-parse_seqids")

#    How many sequences are present in the E.coli set?

## makeblastdb() console result
# Building a new DB, current time: 05/10/2020 23:30:48
# New DB name:   /mnt/student11/projects/SLE712-Assignment-3/Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa
# New DB title:  Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 4140 sequences in 0.330607 seconds.



# 2. Download the sample fasta sequences and read them in as above. 
#    For your allocated sequence, determine the length (in bp) and the proportion of GC bases.
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "Data/sample.fa")
sample_fasta <- seqinr::read.fasta("Data/sample.fa")

# Subset sequence 11 
seq11 <- sample_fasta[[11]]

# Sequence length in bp
seqinr::getLength(seq11)
# 1497 bp

# Proportion of GC bases
seqinr::GC(seq11)
# 0.5744823



# 3. You will be provided with R functions to create BLAST databases and perform blast searches. 
#    Use blast to identify what E. coli gene your sequence matches best. Show a table of the top 
#    3 hits including percent identity, E-value and bit scores.

myblastn_tab # provided R function

ecoli_seq <- seqinr::read.fasta("Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")

results <- myblastn_tab(myseq = seq11, db = "Data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
results

top3_hits <- results[1:3,]
top3_hits

#       qseqid   sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# 1        11 AAC76604    100   1497        0       0      1 1497      1 1497      0     2878
# NA       NA     <NA>     NA     NA       NA      NA     NA   NA     NA   NA     NA       NA
# NA.1     NA     <NA>     NA     NA       NA      NA     NA   NA     NA   NA     NA       NA


# 4. You will be provided with a function that enables you to make a set number of point mutations 
#    to your sequence of interest. Run the function and write an R code to check the number of 
#    mismatches between the original and mutated sequence.

mutator # Provided R function
# create a mutated copy with 100 substitutions
seq11_mut <- mutator(myseq=seq11,100)

# now create a pairwise alignment
seq11_mut_ <- DNAString(c2s(seq11_mut))
seq11_ <- DNAString(c2s(seq11))
aln <- Biostrings::pairwiseAlignment(seq11_,seq11_mut_)
pid(aln)
# 94.92318
nmismatch(aln)
# 76

# 5. Using the provided functions for mutating and BLASTing a sequence, determine the number 
#    and proportion of sites that need to be altered to prevent the BLAST search from matching the 
#    gene of origin. Because the mutation is random, you may need to run this test multiple times 
#    to get a reliable answer.

# Write a fasta file and make a blast db from the traget sequence, seq11

write.fasta(seq11, names= "seq11", file.out = "Data/seq11.fa")
makeblastdb(file = "Data/seq11.fa", dbtype = "nucl")


# blast_tester is a function to test the limits of a BLAST search with the following inputs

# init_mut      initial number of mutations
# mut_incr      number of mutations added per iteration
# target_seq    sequence of interest
# blast_db      blast database you are comparing your sequence to
#               (in this case, blast_db is made from the original sequence, seq11)

blast_tester <- function(init_mut, mut_incr, target_seq, blast_db){
        # number of mutations
        mut <- init_mut  
        seq_mut <- mutator(myseq=target_seq,mut)
        results <- myblastn_tab(myseq = seq_mut, db = blast_db)
        
        # Save BLAST search results in a dataframe
        results_table <- as.data.frame(results)
        # Insert new columnn for number of mutations  
        results_table$num_mut <- mut
        
        # Keep mutating until BLAST search returns null
        while (!is.null(results)){
                # Number of added mutations per iteration
                mut = mut + mut_incr
                seq_mut <- mutator(myseq=target_seq,mut)
                results <- myblastn_tab(myseq = seq_mut, db = blast_db)
                if (is.null(results)){ # Do not append the null search result to the table
                        results_table 
                # Append search results and mutations if it is not empty
                } else (results_table[nrow(results_table) + 1,] <- c(results,mut)) 
        }
        return(results_table)
}


# Test the limits of BLAST search with different initial mutations 
# and increments using the blast_tester function

test1 <- blast_tester(1,1,seq11,"Data/seq11.fa")
test2 <- blast_tester(1,5,seq11,"Data/seq11.fa")
test3 <- blast_tester(1,10,seq11,"Data/seq11.fa") 
test4 <- blast_tester(2,2,seq11,"Data/seq11.fa")
test5 <- blast_tester(2,5,seq11,"Data/seq11.fa")

# Merge all test results in one table

all_tests <- rbind(test1, test2, test3, test4, test5)

# 6. Provide a chart or table that shows how the increasing proportion of mutated bases reduces 
#    the ability for BLAST to match the gene of origin. Summarise the results in 1 to 2 sentences.

x <- results_table$pident/100
y <- results_table$mismatch/1497

plot(results_table$mismatch,x)

plot(y,x)