# In this analysis, the pickrell1 rna-seq dataset from a study on RNA sequencing by Joseph Pickrell et al. The dataset will be analyzed based on stratification by gender, and an analysis of the relationship between mean and variance for expression will also be conducted. The distribution of the mean and variance for expression data will be compared to the Poisson distribution. 

# We can begin by loading the dataset and some necessary libraries.

# Source biocLite, which is used to install Bioconductor packages in R
source("https://bioconductor.org/biocLite.R")

# Install and load the edgeR package
  
# biocLite("edgeR")
library(edgeR)

# Install and load the tweeDEseqCountD package

# biocLite("tweeDEseqCountData")
library(tweeDEseqCountData)

# Load the rna-seq data from the pickrell1 dataset
data(pickrell1)

# Create object 'Counts' to hold read count information
Counts <- exprs(pickrell1.eset)

# Define and extract the variable of interest, which is gender
Gender <- pickrell1.eset$gender

# Prepare and extract annotation information that will contain genes names to be used later as part of the DGEList below
data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]

# Remove annotation ensembl object
rm(annotEnsembl63)

# Create DGEList object and add the read counts and annotations
y <- DGEList(counts=Counts, genes=annot[rownames(Counts),])

# Define genes with high enough expression levels. Here we define it to be genes with more than 1 cpm (count per million reads) in at least 20 samples
isexpr <- rowSums(cpm(y)>1) >= 20

# Define whether genes have annotated functions
hasannot <- rowSums(is.na(y$genes))==0

# Pick only genes with high enough expression AND having annotated functions
y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]

# Check libsize for the two genders using a boxplot 
boxplot(split(log10(y$samples$lib.size),Gender),main="Gender Based Library Size Comparison for Pickrell1", ylab="log10 library size")

#From the library size comparison, we can see that the sizes are roughly the same for males and females, except that females have somewhat significantly greater number of smaller library sizes. 

# We want to compare mean-var relationship at gene-level relative to Poisson
# First let's get the mean and coefficient of variation values 
gene.mean <- apply(y$count,1,mean)
CV.sq <- apply(y$count,1,var)/gene.mean^2

# Plot log of mean and coefficient of variation for the data
plot(log(gene.mean), log(CV.sq),cex=0.3,col="grey",ylim=c(-5,3), main="Mean-Variance Relationship for Pickrell1 Dataset Compared to Poisson")

# Add linear regression line for the mean-var relationship
abline(lm(log(gene.mean)~log(CV.sq)),col='red')

# Add predicted regression line under Poisson 
abline(c(0,-1),col='blue')

# The red line indicates the regression line for the actual dataset before any normalization, and the blue line indicates the expected regression slope under the Poisson distribution. Therefore, we can see that the relationship between the mean and the variance for this dataset follows very closely to the expected relationship under the Poisson distribution.

# Now we can compare the mean-var relationship to Poisson the same way after normalization of the expression counts. 

# First calculate relative NORMALIZATION factors (using default method)
y <- calcNormFactors(y)

# Obtain normalized counts for expression
sf <- colSums(y$count); sf <- sf/mean(sf)
eff.sf <- sf*y$samples$norm.factors
norm.y <- y$count %*% diag(1/eff.sf)

# Again compare mean-var relationship at gene-level relative to Poisson
gene.mean <- apply(norm.y,1,mean)
CV.sq <- apply(norm.y,1,var)/gene.mean^2
plot(log(gene.mean), log(CV.sq),cex=0.3,col="grey",ylim=c(-5,3), main="Normalized Mean-Variance Relationship for Pickrell1 Dataset Compared to Poisson")

# Add linear regression line for the mean-var relationship
abline(lm(log(gene.mean)~log(CV.sq)), col='red')

# Add predicted regression line under Poisson 
abline(c(0,-1),col='blue')

# There is quite a difference in the slope of the actual regression line for the mean-var relationship after normalization of the read counts. The slope (from the red regression line), seems to increase as compared to before normalization, and is therefore not as well represented by the Poisson model anymore. 

# From this analysis, we can conclude that before normalization, the mean-variance relationship for the Pickrell1 rna-seq dataset was roughly approximated by the Poisson distribution. However, after normalization of the counts, we can see that the mean-var relationship became skewed. The y-intercept for the mean-var relationship regression line increased as compared to before, indicating that overall variance may have increased after normalization.

# Further analysis would involve a more thorough investigation of the relationship, such as quantitatively comparing the slopes to determine the difference. Also other commonly used distributions can be compared to the rna-seq data, such as the negative binomial distribution. The data also looked to be severely affected by outlier, so another aspect of analysis could involve removing outliers, or normalizing their values using some sort of 'fencing' technique. 
