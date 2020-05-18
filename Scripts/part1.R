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
head(gene_ordered, 10)

# 4. Determine the number of genes with a mean <10
lessthan10 <- gene_expression_data[(gene_expression_data[,3]<10),]
nrow(lessthan10) # 43124 genes

# 5. Make a histogram plot of the mean values in png format and paste it into your report.
png(filename = "Data/part1_histogram" )
hist(gene_expression_data$Mean, breaks = 5, xlab = "Mean", main = "Histogram of Gene Expression Data Mean")
dev.off()

# 6. Import this csv file into an R object.
#    What are the column names?

growth_data <- read.csv("~/projects/SLE712-Assignment-3/Data/part1_growth_data.csv",
                        header = TRUE, stringsAsFactors = FALSE)

column <- colnames(growth_data)
column

# 7. Calculate the mean and standard deviation of tree circumference at the start and end of the study 
#    at both sites.

# For Northeast
ne <- subset(growth_data, Site=="northeast")
head(ne)

# Mean for Northeast data
mean_end2 <- mean(ne$Circumf_2004_cm)
mean_end1 <- mean(ne$Circumf_2019_cm)
mean_end1
mean_end2

# Standard Deviation for Northeast data
sd(ne$Circumf_2004_cm)
sd(ne$Circumf_2019_cm)

# For Southwest
sw <- subset(growth_data, Site == "southwest")
head(sw)

# Mean for Southwest data
mean_start2 <- mean(sw$Circumf_2004_cm)
mean_end2 <- mean(sw$Circumf_2019_cm)
mean_start2
mean_end2

# Standard Deviation for Southwest data
sd(sw$Circumf_2004_cm)
sd(sw$Circumf_2019_cm)

# 8. Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(ne$Circumf_2004_cm, ne$Circumf_2019_cm, sw$Circumf_2004_cm, sw$Circumf_2019_cm,
        names = c("NE 2004", "NE 2019", "SW 2004","SW 2019"),
        ylab= "Circumfrence (cm)" , xlab = "Sites and years" ,
        main = "Growth at two different sites during 2004 and 2019", col= "green")

# 9. Calculate the mean growth over the past 10 years at each site.

# Mean growth for Northeast data
ne$Circumf_2019_cm - ne$Circumf_2009_cm
ne$growth <- (ne$Circumf_2019_cm - ne$Circumf_2009_cm)
mean_growth_ne <- mean(ne$growth)
mean_growth_ne

# Mean growth for Southwest data
sw$Circumf_2019_cm - sw$Circumf_2009_cm
sw$growth <- (sw$Circumf_2019_cm - sw$Circumf_2009_cm)
mean_growth_sw<- mean(sw$growth)
mean_growth_sw

# 10. Use the t.test and wilcox.test functions to estimate the p-value that the 10 year growth is different at the two sites.

# t test
t_test <- t.test(ne$growth,sw$growth)
t_test

# wilcox.test
wilcox_test <- wilcox.test(ne$growth, sw$growth)
wilcox_test
