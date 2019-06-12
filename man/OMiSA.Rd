\name{OMiSA}
\alias{OMiSA}
\title{
Optimal Microbiome-based Survival Analysis (OMiSA)
}
\description{
OMiSA is a non-parametric method which tests the association between a microbial group (e.g., community, taxon) composition and a survival (time-to-event) response on human health or disease with or without covariate adjustments (e.g., age, sex).
}
\usage{
OMiSA(obstime, delta, X, total.reads = NULL, tree, cov = NULL, 
pow = c(1/4,1/3,1/2,1), g.unif.alpha = c(0.5), n.perm = 5000)
}

\arguments{
  \item{obstime}{
A numeric vector for the observed times.
}
  \item{delta}{
A numeric vector for the event/censoring indicators (1: event, 0: censoring).
}
  \item{X}{
A matrix for the OTU table (1. Elements are counts. 2. Rows are samples and columns are OTUs. 3. Monotone/singletone OTUs must be removed.).
}
  \item{total.reads}{
A numeric vector for the total reads per sample in the entire community. If you survey the entire community, you do not need to specify this. If you test a microbial taxon, you need to specify this (see the examples below). Default is NULL for the entire community.
}
  \item{tree}{
A rooted phylogenetic tree.
}
  \item{cov}{
A data.frame (or vector) for covariate adjustment(s) (Rows are samples and columns are covariate variables).
}
  \item{pow}{
A set of the candidate gamma values. Default is c(1/4, 1/3, 1/2, 1).
}
  \item{g.unif.alpha}{
A set of the candidate alpha parameter value(s) for the generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5).
}
  \item{n.perm}{
The number of permutations. Default is 5000.  
}
}

\value{
pvs.misaln: The estimated p-values for individual MiSALN tests

pvs.mirkats: The estimated p-values for individual MiRKAT-S tests 

p.omisaln: The estimated p-value for OMiSALN

p.omirkats: The estimated p-value for OMiRKAT-S 

p.omisa: The estimated p-value for OMiSA 
}
\references{
Koh H, Livanos AE, Blaser MJ, Li H. (2018) A highly adaptive microbiome-based association test for survival traits. BMC Genomics 19, 210

Plantinga A, Zhan X, Zhao N, Chen J, Jenq RR, Wu MC. (2017) MiRKAT-S: a community-level test of association between the microbiota and survival times. Microbiome 5, 17
}
\author{
Hyunwook Koh
}

\examples{
library(dirmult) 
library(phyloseq) 
library(robustbase)
library(robCompositions) 
library(BiasedUrn)
library(CompQuadForm)
library(GUniFrac) 
library(ecodist) 
library(survival)
library(OMiSA)

### An example simulated microbiome data with survival outcomes

data(MiSurv.Data) 

otu.tab <- otu_table(MiSurv.Data)
tax.tab <- tax_table(MiSurv.Data)
tree <- phy_tree(MiSurv.Data)
obstime <- as.numeric(unlist(sample_data(MiSurv.Data)[,1]))
delta <- as.numeric(unlist(sample_data(MiSurv.Data)[,2]))
x1 <- as.numeric(unlist(sample_data(MiSurv.Data)[,3]))
x2 <- as.numeric(unlist(sample_data(MiSurv.Data)[,4]))
covs <- as.data.frame(cbind(x1, x2))
covs[,2] <- as.factor(covs[,2])

### Example 1. To test the entire community (e.g., kingdom)

set.seed(100)

OMiSA(obstime, delta, otu.tab, total.reads=NULL, tree, cov=covs)

### Example 2. To test the higher-level taxon, p__Firmicutes

# Important notice: Create n total reads counts "in the entire community" and 
# specify it into the element, "total.reads", of the function, OMiSA, as below.
total.reads <- rowSums(otu.tab)

ind.Firmicutes <- which(tax.tab[,2] == "p__Firmicutes")
otu.tab.Firmicutes <- otu.tab[,ind.Firmicutes]
tree.Firmicutes <- prune_taxa(colnames(otu.tab.Firmicutes), tree)

set.seed(100)

OMiSA(obstime, delta, otu.tab.Firmicutes, total.reads=total.reads, 
tree.Firmicutes, cov=covs)

### Example 3. To test the higher-level taxon, p__Bacteroidetes

# Important notice: Create n total reads counts "in the entire community" and 
# specify it into the element, "total.reads", of the function, OMiSA, as below.
total.reads <- rowSums(otu.tab)

ind.Bacteroidetes <- which(tax.tab[,2] == "p__Bacteroidetes")
otu.tab.Bacteroidetes <- otu.tab[,ind.Bacteroidetes]
tree.Bacteroidetes <- prune_taxa(colnames(otu.tab.Bacteroidetes), tree)

set.seed(100)

OMiSA(obstime, delta, otu.tab.Bacteroidetes, total.reads=total.reads, 
tree.Bacteroidetes, cov=covs)

}