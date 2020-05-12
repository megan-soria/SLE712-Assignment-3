# 1. Read in the file, making the gene accession numbers the row names. 
#    Show a table of values for the first six genes.

gene_expression_data <- read.csv("~/projects/SLE712-Assignment-3/Data/part1_gene_expression.tsv", 
                                 sep = '\t', row.names = "GeneID")
head(gene_expression_data, n=6)

# 2. Make a new column which is the mean of the other columns. 
#    Show a table of values for the first six genes.

gene_expression_data$Mean <- rowMeans(gene_expression_data[,1:2])
head(gene_expression_data, 6)

# 3. List the 10 genes with the highest mean expression
gene_ordered <- gene_expression_data[order(-gene_expression_data$Mean),]
head(gene_ordered, 20)

# 4. Determine the number of genes with a mean <10
lessthan10 <- gene_expression_data[(gene_expression_data[,3]<10),]
nrow(lessthan10) # 43124 genes

# 5. Make a histogram plot of the mean values in png format and paste it into your report.
hist_data <- gene_expression_data[,3]

hist(hist_data, xlim = c(0, 100), ylim = c(0, 100000))


# 6. Import this csv file into an R object.
# What are the column names?

growth_data <- read.csv("~/projects/SLE712-Assignment-3/Data/part1_growth_data.csv",
                        header = TRUE, stringsAsFactors = FALSE)

column <- colnames(growth_data)
column

# 7. Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites.
# For northeast
ne <- subset(growth_data, Site=="northeast")
head(ne)

mean_end2 <- mean(ne$Circumf_2004_cm)
mean_end1 <- mean(ne$Circumf_2019_cm)

sd(ne$Circumf_2004_cm)
sd(ne$Circumf_2019_cm)

# For southwest
sw <- subset(growth_data, Site == "southwest")
head(sw)

mean_start2 <- mean(sw$Circumf_2004_cm)
mean_end2 <- mean(sw$Circumf_2019_cm)

sd(sw$Circumf_2004_cm)
sd(sw$Circumf_2019_cm)

# 8. Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(ne$Circumf_2004_cm, ne$Circumf_2019_cm, sw$Circumf_2004_cm, sw$Circumf_2019_cm,
        names = c("NE 2004", "NE 2019", "SW 2004","SW 2019"),
        ylab= "Circumfrence (cm)" , xlab = "Sites and year" ,
        main = "Growth at two different sites during 2004 and 2019")

# 9. Calculate the mean growth over the past 10 years at each site.
