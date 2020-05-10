
# 1. Read in the file, making the gene accession numbers the row names. 
#    Show a table of values for the first six genes.

gene_expression_data <- read.csv("~/projects/SLE712-Assignment-3/Data/part1_gene_expression.tsv", 
                                 sep = '\t', row.names = "GeneID")
head(gene_expression_data, n=6)

# 2. Make a new column which is the mean of the other columns. 
#    Show a table of values for the first six genes.

gene_expression_data$Mean <- rowMeans(gene_expression_data[,1:2])
head(gene_expression_data, n=6)

# 3. List the 10 genes with the highest mean expression
gene_ordered <- gene_expression_data[order(-gene_expression_data$Mean),]
head(gene_ordered, 10)

# 4. Determine the number of genes with a mean <10
lessthan10 <- gene_expression_data[(gene_expression_data[,3]<10),]
nrow(lessthan10)

# 5. Make a histogram plot of the mean values in png format and paste it into your report.
hist(gene_expression_data$Mean)

