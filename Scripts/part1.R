# 1. Read in the file, making the gene accession numbers the row names. 
#    Show a table of values for the first six genes.

gene_expression_data <- read.csv("Data/part1_gene_expression.tsv", 
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
mean_lessthan10 <- gene_expression_data[(gene_expression_data[,3]<10),]
nrow(mean_lessthan10) # 43124 genes

# 5. Make a histogram plot of the mean values in png format and paste it into your report.
png(filename = "Data/part1_histogram.png")
hist(gene_expression_data$Mean, xlab = "Mean", main = "Histogram of Gene Expression Data Mean")
dev.off()

# The histogram for the Mean column produces only one bar graph near Mean = zero. 
# This tells us that the number of genes with mean=0 & mean<500 far exceeds the number of genes with, 
# say, mean =< 500. The frequency is skewed to the right and does not give us a proper picture of the data. 
# We can solve this by subsetting the Mean column and taking out the extremes. 

# To inspect the data, we can use the ftable() function to give us the numerical data of the
# mean frequencies

ftable(gene_expression_data$Mean)

# -----Console Result-----
#     0   0.5     1   1.5     2   2.5     3   3.5     4   4.5     5
# 29694  4269  2130  1392   976   766   586   517   388   367   298

mean(gene_expression_data$Mean) # column mean = 354.5272

# From the result above, we notice that our observation was correct and the count for mean = 0 skews the data.
# In the script below, the mean of the column was set as the maximum (354.5272 ~ 360) and only mean values greater than 2 is 
# set as minimum. By removing the extremes, the histogram can now reveal how the frequency of the mean steeply
# drops after about the first 10 bars.

png(filename = "Data/part1_histogram_adjusted.png")
with(gene_expression_data, hist
     (Mean[Mean>2 & Mean<360], breaks=seq(2,360,by=1), 
                                xlab = "Mean", main = "Histogram of Gene Expression Data Mean>2 & Mean<360"))
dev.off()

# We can also recommend using other data visualization tools to represent the frequencies of the means
# like a scatterplot to accurately capture the data.

# 6. Import this csv file into an R object.
#    What are the column names?

growth_data <- read.csv("Data/part1_growth_data.csv",
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
png(filename = "Data/part1_boxplot.png")
boxplot(ne$Circumf_2004_cm, ne$Circumf_2019_cm, sw$Circumf_2004_cm, sw$Circumf_2019_cm,
        names = c("NE 2004", "NE 2019", "SW 2004","SW 2019"),
        ylab= "Circumfrence (cm)" , xlab = "Sites and years" ,
        main = "Growth at two different sites during 2004 and 2019", col= "green")
dev.off()

# 9. Calculate the mean growth over the past 10 years at each site.

# Mean growth for Northeast data

ne$growth <- (ne$Circumf_2019_cm - ne$Circumf_2009_cm)
mean_growth_ne <- mean(ne$growth)
mean_growth_ne

# Mean growth for Southwest data

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
