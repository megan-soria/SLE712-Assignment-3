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

# How many sequences are present in the E.coli set? 4140 sequences


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


# blast_lim is a function that tests the maximum number of mutations that can still return
# a BLAST search match when compared to the original sequence. It takes an initial number of 
# mutations used to mutate the original sequence, makes a BLAST search, and repeats this process in 
# defined increments until the search returns NULL. It stores each iteration in a table with the 
# last row as the highest number of mutations that returned a match. The following are the inputs:

# init_mut      initial number of mutations
# mut_incr      number of mutations added per iteration


blast_lim <- function(init_mut, mut_incr){
        # number of mutations
        mut <- init_mut  
        seq_mut <- mutator(myseq=seq11,mut)
        results <- myblastn_tab(myseq = seq_mut, db = "Data/seq11.fa")
        
        # Save BLAST search results in a dataframe
        results_table <- as.data.frame(results)
        # Insert new columnn for number of mutations  
        results_table$num_mut <- mut
        
        # Keep mutating until BLAST search returns null
        while (!is.null(results)){
                # Number of added mutations per iteration
                mut = mut + mut_incr
                seq_mut <- mutator(myseq=seq11,mut)
                results <- myblastn_tab(myseq = seq_mut, db = "Data/seq11.fa")
                if (is.null(results)){ # Do not append the null search result to the table
                        results_table 
                # Append search results and mutations if it is not empty
                } else (results_table[nrow(results_table) + 1,] <- c(results,mut)) 
        }
        return(results_table)
}


# Test the limits of BLAST search with different initial mutations 
# and increments using the blast_tester function

test1 <- blast_lim(1,1)
test2 <- blast_lim(1,10)
test3 <- blast_lim(1,20) 
test4 <- blast_lim(1,30)
test5 <- blast_lim(2,50) 


# Merge all test results in one table and take the top 10 highest number of mutations

all_tests <- rbind(test1, test2, test3, test4, test5)
max_mut <- all_tests[order(-all_tests$num_mut),]
head(max_mut, 10)

# blast_tester mutates a sequence in a defined number of places ("mut").
# If a BLAST search against the original sequence returns a match, the function returns a 1. 
# If the search result is NULL, it returns a 0. Input:

# mut   number of mutations to be applied in the sequence

blast_tester <- function(mut){
        seq_mut <- mutator(myseq=seq11,mut)
        results <- myblastn_tab(myseq = seq_mut, db = "Data/seq11.fa")
        if (!is.null(results)){
                return(1)
        } else (return(0))
}


# Since the mutations are random, the BLAST search results changes in each run.
# The following code uses the results from the blast_lim function and replicates the run 
# of the blast_tester function 100 times to get a mean value and a better grasp of the
# BLAST search behavior


# Create an empty data frame for blast_tester function results
blast_test_res <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("num_of_mut", "Mean_blast_res"))))

i <- 0

while (i < 800){
        mean_blast_res <- mean(replicate(100,blast_tester(i)))
        blast_test_res[nrow(blast_test_res) + 1,] <- c(i,mean_blast_res)
        i <- i + 50
}
blast_test_res

# 6. Provide a chart or table that shows how the increasing proportion of mutated bases reduces 
#    the ability for BLAST to match the gene of origin. Summarise the results in 1 to 2 sentences.

blast_test_res$sites <- blast_test_res$num_of_mut

blast_test_res$prop <- blast_test_res$Mean_blast_res

blast_test_res$random <- blast_test_res$num_of_mut/1497

library("ggplot2")
png(filename = "Data/part2_blast1.png")
ggplot(blast_test_res, aes(sites,prop)) +
        geom_line(color = "violet", size = 1 ) + 
        geom_point(shape=19,color="turquoise3", size= 3) + 
        theme_bw() + labs(title="100 iterations, using sequence number 11", x="Number of sites randomised",y="Proportion of successful BLASTs")
dev.off()

png(filename = "Data/part2_blast2.png")
ggplot(blast_test_res, aes(random,prop)) + 
        geom_line(color = "violet", size = 1 ) + 
        geom_point(shape=19,color="turquoise3", size= 3) + 
        theme_bw()+ labs(title="100 iterations, using sequence number 11", x="Number of sites randomised",y="Proportion of successful BLASTs")
dev.off()
