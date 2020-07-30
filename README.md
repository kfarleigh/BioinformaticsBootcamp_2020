# Exploring genome-environment associations with Redundnacy Analysis

## Purpose
To use Redundnacy analysis (RDA) to explore genome-environment associations and detect loci under selection. 

## Overview
RDA is a multivariate ordination technique that can analyze many loci and environmental predictors simaltaneously. This procedure determines how loci covary in response to a multivariate environment and is capable of detecting weak multilocus selection ([Rellstab et al., 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13322); [Forester et al., 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=ghJOZ0-A-a0AAAAA%3Ajp7tQPDyymFyemaC7poslehIB1Eli3WhPAqea9nMhasQs4NPFK0rg6jf8V498iG8ixedk2R6optMm4U)). RDA is a two step procedure where multivarite linear regression is used to analyze genetic and environmental data, creating a matrix of fitted values. Then PCA is performed on the fitted values to produce canonical axes. We then use SNP loadings on these axes to identify candidates that may be under selection. Finally, we determine the correlation between each predictor and candidate.

![RDA_Concept](https://user-images.githubusercontent.com/54188848/88552471-22bbf500-cff2-11ea-9042-72b145f037b4.JPG)

### Things to consider before you start your analysis
1. What is your sampling design? 
RDA can accommodate both individual and population based sampling designs, though your preperation before RDA will be different depending on your sampling. Individual based sampling will use allele counts, whereas population sampling will use allele frequencies within demes. 
2. Does linkage disequlibrium affect your dataset? 
