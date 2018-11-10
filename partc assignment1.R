# Source biocLite, which is used to install Bioconductor packages in R
source("https://bioconductor.org/biocLite.R")
# Install the edgeR package
biocLite("edgeR")
# Load the edgeR package
library(edgeR)
# Source biocLite, which is used to install Bioconductor packages in R
source("https://bioconductor.org/biocLite.R")
# Install the tweeDEseqCountData package
biocLite("tweeDEseqCountData")
# Load the tweeDEseqCountData package
library(tweeDEseqCountData)

#Load the data from the pickrell1 dataset
data(pickrell1)

# Create object 'Counts' to hold read count information
Counts <- exprs(pickrell1.eset)

# Define the variable of interest
Gender <- pickrell1.eset$gender

# Prepare and extract annotation information that will contain genes names to be used later as part of the DGEList below.
data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]

# Remove object
rm(annotEnsembl63)

#Create DGEList object and add the read counts and annotations
y <- DGEList(counts=Counts, genes=annot[rownames(Counts),])

# define genes with high enough expression levels. Here we define it to be genes with more than 1 cpm (count per million reads) in at least 20 samples
isexpr <- rowSums(cpm(y)>1) >= 20

# define whether genes have annotated functions
hasannot <- rowSums(is.na(y$genes))==0

# pick only genes with high enough expression AND has annotated functions
y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]

# Check libsize for the two genders using a boxplot 
boxplot(split(log10(y$samples$lib.size),Gender),main="log10 library size")

# Compare mean-var relationship at gene-level relative to Poisson
gene.mean <- apply(y$count,1,mean)
CV.sq <- apply(y$count,1,var)/gene.mean^2

# Plot CV
plot(log(gene.mean), log(CV.sq),cex=0.3,col="grey",ylim=c(-5,3))

# Add predicted CV.sq under Poisson (do you know why under Poisson model the slope is -1?
abline(c(0,-1))

# First calculate relative NORMALIZATION factors (using default method)
y <- calcNormFactors(y)

# Obtain normalized count
sf <- colSums(y$count); sf <- sf/mean(sf)
eff.sf <- sf*y$samples$norm.factors
norm.y <- y$count %*% diag(1/eff.sf)

# Compare mean-var relationship at gene-level relative to Poisson
gene.mean <- apply(norm.y,1,mean)
CV.sq <- apply(norm.y,1,var)/gene.mean^2
plot(log(gene.mean), log(CV.sq),cex=0.3,col="grey",ylim=c(-5,3))

# Add predicted CV.sq under Poisson 
abline(c(0,-1))

##Any conclusions?
