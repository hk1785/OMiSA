# R package: OMiSA

Version: 1.5

Date: 2019-6-11

Title: Optimal Microbiome-based Survival Analysis (OMiSA)

Author: Hyunwook Koh, Alexandra E. Livanos, Martin J. Blaser, Huilin Li

Maintainer: Hyunwook Koh hk1785@nyu.edu

Description: This software package provides facilities for the non-parametric methods 1) Optimal Microbiome-based Survival Analysis (OMiSA), 2) Optimal Microbiome-based Survival Analysis using Linear and Non-linear bases of OTUs (OMiSALN) and 3) Optimal Microbiome Regression-based Kernel Association Test for Survival traits (OMiRKAT-S). OMiSA, OMiSALN and OMiRKAT-S test the association between a microbial group (e.g., community, taxon) composition and a survival (time-to-event) response on human health/disease with or without covariate adjustments (e.g., age, sex).

NeedsCompilation: No

Depends: R(>= 3.2.3)

Imports: BiasedUrn, CompQuadForm, dirmult, ecodist, GUniFrac, phyloseq, robCompositions, robustbase, survival

License: GPL-2

## Reference

* Koh H, Livanos AE, Blaser MJ, Li H. (2018) A highly adaptive microbiome-based association test for survival traits. _BMC Genomics_ 19:210.

* DOI: https://doi.org/10.1186/s12864-018-4599-8

## Installation
```
library(devtools)
install_github("hk1785/OMiSA", force=T)
```

## Data format
```
library(phyloseq)
```
* URL: https://joey711.github.io/phyloseq/

## Prerequites

BiasedUrn
```
install.packages("BiasedUrn")
```
CompQuadForm
```
install.packages("CompQuadForm")
```
dirmult
```
install.packages("dirmult")
```
ecodist
```
install.packages("ecodist")
```
GUniFrac
```
install.packages("GUniFrac")
```
phyloseq
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```
robCompositions
```
install.packages("robCompositions")
```
robustbase
```
install.packages("robustbase")
```
survival
```
install.packages("survival")
```

## Installation

```
library(devtools)
install_github("hk1785/OMiSA", force=T)
```

## Data format

```
library(phyloseq)
URL: https://joey711.github.io/phyloseq/
```

---------------------------------------------------------------------------------------------------------------------------------------

# Manual
This R package includes two core functions, OMiSA, OMiSALN and OMiRKAT. Please find the details below.

## :mag: OMiSA

### Description
OMiSA is a non-parametric method which tests the association between a microbial group (e.g., community, taxon) composition and a survival (time-to-event) response on human health or disease with or without covariate adjustments (e.g., age, sex).

### Usage
```
OMiSA(obstime, delta, X, total.reads = NULL, tree, cov = NULL, pow = c(1/4,1/3,1/2,1), g.unif.alpha = c(0.5), n.perm = 5000)
```

### Arguments
* _obstime_ - A numeric vector for the observed times.
* _delta_ - A numeric vector for the event/censoring indicators (1: event, 0: censoring).
* _X_ - A matrix for the OTU table (1. Elements are counts. 2. Rows are samples and columns are OTUs. 3. Monotone/singletone OTUs must be removed.).
* _total.reads_ - A numeric vector for the total reads per sample in the entire community. If you survey the entire community, you do not need to specify this. If you test a microbial taxon, you need to specify this (see the examples below). Default is NULL for the entire community.
* _tree_ - A rooted phylogenetic tree.
* _cov_ - A data.frame (or vector) for covariate adjustment(s) (Rows are samples and columns are covariate variables).
* _pow_ - A set of the candidate gamma values. Default is c(1/4, 1/3, 1/2, 1).
* _g.unif.alpha_ - A set of the candidate alpha parameter value(s) for the generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5).
* _n.perm_ - The number of permutations. Default is 5000.  

### Values
_$pvs.misaln_ - The estimated p-values for individual MiSALN tests

_$pvs.mirkats_ - The estimated p-values for individual MiRKAT-S tests 

_$p.omisaln_ - The estimated p-value for OMiSALN

_$p.omirkats_ - The estimated p-value for OMiRKAT-S 

_$p.omisa_ - The estimated p-value for OMiSA 

### References
* Koh H, Livanos AE, Blaser MJ, Li H. (2018) A highly adaptive microbiome-based association test for survival traits. BMC Genomics 19, 210
* Plantinga A, Zhan X, Zhao N, Chen J, Jenq RR, Wu MC. (2017) MiRKAT-S: a community-level test of association between the microbiota and survival times. Microbiome 5, 17

### Example
Import requisite R packages
```
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
```
Import example microbiome data
```
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
```

Example 1. To test the entire community (e.g., kingdom)
```
set.seed(100)
OMiSA(obstime, delta, otu.tab, total.reads=NULL, tree, cov=covs)
```

Example 2. To test the higher-level taxon, p__Firmicutes
```
total.reads <- rowSums(otu.tab)
ind.Firmicutes <- which(tax.tab[,2] == "p__Firmicutes")
otu.tab.Firmicutes <- otu.tab[,ind.Firmicutes]
tree.Firmicutes <- prune_taxa(colnames(otu.tab.Firmicutes), tree)

set.seed(100)
OMiSA(obstime, delta, otu.tab.Firmicutes, total.reads=total.reads, 
tree.Firmicutes, cov=covs)
```

Example 3. To test the higher-level taxon, p__Bacteroidetes
```
total.reads <- rowSums(otu.tab)
ind.Bacteroidetes <- which(tax.tab[,2] == "p__Bacteroidetes")
otu.tab.Bacteroidetes <- otu.tab[,ind.Bacteroidetes]
tree.Bacteroidetes <- prune_taxa(colnames(otu.tab.Bacteroidetes), tree)

set.seed(100)
OMiSA(obstime, delta, otu.tab.Bacteroidetes, total.reads=total.reads, 
tree.Bacteroidetes, cov=covs)
```

## :mag: OMiSALN

### Description
OMiSALN is a non-parametric method which tests the association between a microbial group (e.g., community, taxon) composition and a survival (time-to-event) response on human health or disease with or without covariate adjustments (e.g., age, sex).

### Usage
```
OMiSALN(obstime, delta, X, total.reads = NULL, cov = NULL, pow = c(1/4, 1/3, 1/2, 1), n.perm = 5000)
```

### Arguments
* _obstime_ - A numeric vector for the observed times.
* _delta_ - A numeric vector for the event/censoring indicators (1: event, 0: censoring).
* _X_ - A matrix for the OTU table (1. Elements are counts. 2. Rows are samples and columns are OTUs. 3. Monotone/singletone OTUs must be removed.).
* _total.reads_ - A numeric vector for the total reads per sample in the entire community. If you survey the entire community, you do not need to specify this. If you test a microbial taxon, you need to specify this (see the examples below). Default is NULL for the entire community.
* _tree_ - A rooted phylogenetic tree.
* _cov_ - A data.frame (or vector) for covariate adjustment(s) (Rows are samples and columns are covariate variables).
* _pow_ - A set of the candidate gamma values. Default is c(1/4, 1/3, 1/2, 1).
* _n.perm_ - The number of permutations. Default is 5000.  

### Values
_$pvs.misaln_ - The estimated p-values for individual MiSALN tests

_$p.omisaln_ - The estimated p-value for OMiSALN

### References
* Koh H, Livanos AE, Blaser MJ, Li H. (2018) A highly adaptive microbiome-based association test for survival traits. BMC Genomics 19, 210

### Example
Import requisite R packages
```
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
```
Import example microbiome data
```
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
```

Example 1. To test the entire community (e.g., kingdom)
```
set.seed(100)
OMiSALN(obstime, delta, otu.tab, total.reads=NULL, cov=covs)
```

Example 2. To test the higher-level taxon, p__Firmicutes
```
total.reads <- rowSums(otu.tab)
ind.Firmicutes <- which(tax.tab[,2] == "p__Firmicutes")
otu.tab.Firmicutes <- otu.tab[,ind.Firmicutes]
tree.Firmicutes <- prune_taxa(colnames(otu.tab.Firmicutes), tree)

set.seed(100)
OMiSALN(obstime, delta, otu.tab.Firmicutes, total.reads=total.reads, cov=covs)
```

Example 3. To test the higher-level taxon, p__Bacteroidetes
```
total.reads <- rowSums(otu.tab)
ind.Bacteroidetes <- which(tax.tab[,2] == "p__Bacteroidetes")
otu.tab.Bacteroidetes <- otu.tab[,ind.Bacteroidetes]
tree.Bacteroidetes <- prune_taxa(colnames(otu.tab.Bacteroidetes), tree)

set.seed(100)
OMiSALN(obstime, delta, otu.tab.Bacteroidetes, total.reads=total.reads, cov=covs)
```

## :mag: OMiRKATS

### Description
OMiRKATS is a non-parametric method which tests the association between a microbial group (e.g., community, taxon) composition and a survival (time-to-event) response on human health or disease with or without covariate adjustments (e.g., age, sex).

### Usage
```
OMiRKATS(obstime, delta, X, total.reads = NULL, tree, cov = NULL, g.unif.alpha = c(0.5), n.perm = 5000)
```

### Arguments
* _obstime_ - A numeric vector for the observed times.
* _delta_ - A numeric vector for the event/censoring indicators (1: event, 0: censoring).
* _X_ - A matrix for the OTU table (1. Elements are counts. 2. Rows are samples and columns are OTUs. 3. Monotone/singletone OTUs must be removed.).
* _total.reads_ - A numeric vector for the total reads per sample in the entire community. If you survey the entire community, you do not need to specify this. If you test a microbial taxon, you need to specify this (see the examples below). Default is NULL for the entire community.
* _tree_ - A rooted phylogenetic tree.
* _cov_ - A data.frame (or vector) for covariate adjustment(s) (Rows are samples and columns are covariate variables).
* _g.unif.alpha_ - A set of the candidate alpha parameter value(s) for the generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5).
* _n.perm_ - The number of permutations. Default is 5000.  

### Values
_$pvs.mirkats_ - The estimated p-values for individual MiRKAT-S tests

_$p.omirkats_ - The estimated p-value for OMiRKAT-S

### References
* Koh H, Livanos AE, Blaser MJ, Li H. (2018) A highly adaptive microbiome-based association test for survival traits. BMC Genomics 19, 210
* Plantinga A, Zhan X, Zhao N, Chen J, Jenq RR, Wu MC. (2017) MiRKAT-S: a community-level test of association between the microbiota and survival times. Microbiome 5, 17

### Example
Import requisite R packages
```
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
```
Import example microbiome data
```
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
```

Example 1. To test the entire community (e.g., kingdom)
```
set.seed(100)
OMiRKATS(obstime, delta, otu.tab, total.reads=NULL, tree, cov=covs)
```

Example 2. To test the higher-level taxon, p__Firmicutes
```
total.reads <- rowSums(otu.tab)
ind.Firmicutes <- which(tax.tab[,2] == "p__Firmicutes")
otu.tab.Firmicutes <- otu.tab[,ind.Firmicutes]
tree.Firmicutes <- prune_taxa(colnames(otu.tab.Firmicutes), tree)

set.seed(100)
OMiRKATS(obstime, delta, otu.tab.Firmicutes, total.reads=total.reads, tree.Firmicutes, cov=covs)
```

Example 3. To test the higher-level taxon, p__Bacteroidetes
```
total.reads <- rowSums(otu.tab)
ind.Bacteroidetes <- which(tax.tab[,2] == "p__Bacteroidetes")
otu.tab.Bacteroidetes <- otu.tab[,ind.Bacteroidetes]
tree.Bacteroidetes <- prune_taxa(colnames(otu.tab.Bacteroidetes), tree)

set.seed(100)
OMiRKATS(obstime, delta, otu.tab.Bacteroidetes, total.reads=total.reads, tree.Bacteroidetes, cov=covs)
```
