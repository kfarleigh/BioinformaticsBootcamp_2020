# Exploring genome-environment associations with Redundnacy Analysis

## Purpose
To use Redundnacy analysis (RDA) to explore genome-environment associations and identify candidate loci under selection. This script is modified from ([Forester et al., 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=ghJOZ0-A-a0AAAAA%3Ajp7tQPDyymFyemaC7poslehIB1Eli3WhPAqea9nMhasQs4NPFK0rg6jf8V498iG8ixedk2R6optMm4U)), and the original vignette can be found [here].(https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

## Overview
RDA is a multivariate ordination technique that can analyze many loci and environmental predictors simaltaneously. This procedure determines how loci covary in response to a multivariate environment and is capable of detecting weak multilocus selection ([Rellstab et al., 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13322); [Forester et al., 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=ghJOZ0-A-a0AAAAA%3Ajp7tQPDyymFyemaC7poslehIB1Eli3WhPAqea9nMhasQs4NPFK0rg6jf8V498iG8ixedk2R6optMm4U)). RDA is a two step procedure where multivarite linear regression is used to analyze genetic and environmental data, creating a matrix of fitted values. Then PCA is performed on the fitted values to produce canonical axes. We then use SNP loadings on these axes to identify candidates that may be under selection. Finally, we determine the correlation between each predictor and candidate.

![RDA_Concept](https://user-images.githubusercontent.com/54188848/88552471-22bbf500-cff2-11ea-9042-72b145f037b4.JPG)

### Things to consider before you start your analysis
**1.** What is your sampling design? 

RDA can accommodate both individual and population based sampling designs, though your preperation before RDA will be different depending on your sampling. Individual based sampling will use allele counts, whereas population sampling will use allele frequencies within demes. 

**2.** Does linkage disequlibrium affect your dataset? 

Correlation among loci (e.g. linkage disequilibrium) can result in spurious signal of correlation. You can address this by taking a single snp per locus as implemented in STACKS [populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php)(Rochette et al., 2019) module or you can estimate amounts of linkage disequilibrium yourself with a program like [plink](http://zzz.bwh.harvard.edu/plink/ld.shtml)(Chang et al., 2015) or custom scripts. 
## Load in Data and Packages
We will use genomic data from 64 desert horned lizards (*Phrynosoma platyrhinos*) sampled across western deserts of North America (Great Basin, Mojave, Sonoran). We will use a dataset of 7,500 SNPs and are interested in understanding how individuals may be locally adapted to conditions across their range. For our purposes we will use an individual based RDA, using allele counts for each individual rather than frequencies per deme. The data is availabe in the repository (see data folder).

```
# Load or install packages if necessary
install.packages(c("psych","vegan", "dplyr"), dependencies=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LEA")

library(pysch)
library(vegan)
library(dplyr)
library(LEA)
```

## Analysis 
First we will inspect the data, individuals are represented by rows and SNPs by columns. We also check to make sure that there is no missing data.
```
Genotypes<-read.lfmm("Platyrhinos_subset.lfmm") # Read in genetic data
EnvVars<-read.csv("Platyrhinos_EnvData.csv") # Read in environmental variables

# Inspect the data, look for any missing data
dim(Genotypes)
sum(is.na(Genotypes))
```

### Do the same for the environmental data
```
str(EnvVars)
```

We will change our samples to characters and set the rownames of our genetic data
```
EnvVars$Samples<-as.character(EnvVars$Samples)
Samples<-EnvVars$Samples
row.names(Genotypes)<- Samples
row.names(Genotypes)
identical(rownames(Genotypes), EnvVars$Samples)
```

### Check for correlation between our environmental variables
```
cor.plot(EnvVars [,5:23], xlas =  2, cex.axis = 0.75)
```
Zoom in on the figure, notice that we do have some variables that are correlated (greater than |0.7|). Next we will select 4 variables to perform our analysis with. 
```
# Pull out the selected variables and rename the columns so they are easier to read
SelectedVars<-select(EnvVars,5,8,22,23)
colnames(SelectedVars) <- c("AMT","TS","PWQ","PCQ")

# Check again to make sure our variables are not correlated
pairs.panels(SelectedVars, scale = TRUE)
```
![Pairs_plot](https://user-images.githubusercontent.com/54188848/89056683-6b451c80-d32a-11ea-804d-c14751333e66.jpeg)

We don't have any relationships with a value above 0.7, so we are good to go.

### Run the RDA
Now we run the RDA on our dataset. It should only take a minute maximum, though if you are using large datsetes (tens of thousands of SNPs), the run can take 30 minutes depending on the computer. For our RDA, the data does not include any factors, which means we can use the shortand formula (using ~)to run the analysis. If you did have factors in your dataset you have to write out the formula if you want to test for significance. 
```
# Run the RDA and examine it
HL_rda<- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)
HL_rda
```
Looking at our rda, we have as many constrained RDA axes as we did predictors. The proportion of variance explained by the environmental predictors is provided under the proportion column for the constrained axes.  Next we will calculate an adjusted R-squared. 

```
Rsquared<-RsquareAdj(HL_rda)
Rsquared
```
Our RDA axes explain about 21% of the variation; this is realtively high, as we expect a majority of our SNPs will not show a relationship with environmental predictors. 

Now let's take a look at how much of the variance the individual RDA axes explain
```
summary(eigenvals(HL_rda, model = "constrained"))
```
We see that the first two RDA axes account for most of the variation. If you use RDA in your own research these numbers are often reported. Knowing that variacne explained by our model and individual axes, we will check for the significance of both. These will take a few minutes to run.
```
# test the significance of the rda
Fullsig <- anova.cca(HL_rda, parallel=getOption("mc.cores"))
Fullsig

# find which constrained axes are significant
Axissig <- anova.cca(HL_rda, by="axis", parallel=getOption("mc.cores"))
Axissig
```
Pay attention to the signficance codes produced by each of these. The full significance anova represents our global RDA. The axis signficance represents our RDA axes. 
Since all of our constrained axes are signficant, we will use all of them to identify candidate SNPs under selection. 

But first, one more check to make sure collinearity is not an issue. We will use variance inflation factors (vif) to check. We want values under 10.
```
vif.cca(HL_rda)
```
Since all of our values are below 10 we can move on to some plots.

```
# Make the populations a factor for plotting 
EnvVars$Population = as.factor(EnvVars$Population)
Pop<- EnvVars$Population
bg <- c("#228B22","#ffff33", "#32CD32", "#a6cee3","#33a02c","#e31a1c")
# pull out populations and get some colors

###############
# axes 1 & 2
plot(HL_rda, type="n", scaling=3)
points(HL_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[EnvVars$Population]) #displays individuals
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25)  #displays predictors
legend("bottomright", legend=levels(Pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend

###############
# axes 1 & 3

plot(HL_rda, type="n", scaling=3, choices = c(1,3))
points(HL_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices = c(1,3)) #displays SNPs
points(HL_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[Pop], choices = c(1,3)) #displays individuals
text(HL_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices = c(1,3))  #displays predictors
legend("bottomright", legend=levels(Pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
```
![RDA_12](https://user-images.githubusercontent.com/54188848/89062402-5ec5c180-d334-11ea-986d-af1e03c947c4.jpeg)

We use this plot to evaluate how our individuals (colored dots) relate to the enivornmental predictors. For example, we see that individuals from the Southern cluster have a postitive relationship with annual mean temeperature, whereas individuals from the Great Basin (East & West) have positive relationships with precipitation.

### Test for outliers
To test for outliers, we will look for SNPs in the tails of the loadings. We will use a 3 standard deviation cutoff. This can be modified depending on your goals and tolerance for true postives vs false positves. For example, if you want to be very conservative use a standard deviation of 3.5, but keep in mind that this may increase the false negative rate. Contrastingly you could use a cutoff of 2.5 standard deviations if you were less concerned with false positives.
```
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# writes a function identify outliers

load.rda <- scores(HL_rda, choices=c(1:4), display="species") 
# extract SNP loadings from the significant axes

hist(load.rda[,1], main="Loadings on Axis 1", xlab = NULL)
hist(load.rda[,2], main="Loadings on Axis 2", xlab = NULL)
hist(load.rda[,3], main="Loadings on Axis 3", xlab = NULL)
hist(load.rda[,4], main="Loadings on Axis 4", xlab = NULL)
# plot histograms of the loadings
``` 

![Hist_load](https://user-images.githubusercontent.com/54188848/89063048-7d788800-d335-11ea-97d7-5c4177e9a328.jpeg)

We expect relatively normal distributions since most SNPs will not show a relationship with environmental predictors. 

Now we apply our function to each axis and see how many outliers we have.
```
# apply the function you wrote to each significant axis
# the number by itself indicates how many standard deviations to go
# choosing higher numbers i.e. 3 looks for loci under strong selection, lesser i.e. 2 will be more lenient
cand1 <- outliers(load.rda[,1],3)
cand2 <- outliers(load.rda[,2],3)
cand3 <- outliers(load.rda[,3],3)
cand4 <- outliers(load.rda[,4],3)

# lets see how many outliers we have
ncand <- length(cand1) + length(cand2) +length(cand3) + length(cand4)
ncand
```

We will combine our outliers with the environmental data to determine which predictor is most associated with each candidate.
```
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(1,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(1,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(1,times=length(cand4)), names(cand4), unname(cand4))
colnames(cand1) <- colnames(cand2) <- colnames(cand3)<- colnames(cand4)<- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3, cand4)
cand$snp <- as.character(cand$snp)

# Now we combine this data with our environmental variables

foo <- matrix(nrow=(ncand), ncol=4)
colnames(foo) <- c("AMT","TS","PWQ","PCQ")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- Genotypes[,nam]
  foo[i,] <- apply(SelectedVars,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
```

Next we investigate the candidates, and look for any duplicates
```
# look for any canidate SNPs
length(cand$snp[duplicated(cand$snp)])

# check for duplicates on axis 1
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2])

# check for duplicates on axis 2
table(foo[foo[,1]==2,2])

# remove any duplicate SNPs
cand <- cand[!duplicated(cand$snp),]

# Will tell us which predictor each SNP is associated with
for (i in 1:length(cand$snp)) {
      bar <- cand[i,]
      cand[i,8] <- names(which.max(abs(bar[4:7]))) # makes a column with predictor
      cand[i,9] <- max(abs(bar[4:7])) } #makes a column with correlation

# Adds names to the columns we just made
colnames(cand)[8] <- "predictor"
colnames(cand)[9] <- "correlation"

# Lets see which variables our data are most associated with
table(cand$predictor) 

# Write out the candidates
write.csv(cand, "CanidateSNPs.csv")
```
Finally, we will plot the SNPs
```
sel <- cand$snp
env <- cand$predictor
# take out the SNPs and predictors

env[env=="AMT"] <- '#e31a1c'
env[env=="TS"] <- '#ffa500'
env[env=="PWQ"] <- '#a6cee3'
env[env=="PCQ"] <-  '#2b94db'

# add colors for plotting 

col.pred <- rownames(HL_rda$CCA$v)
# get SNP Names

for (i in 1:length(sel)) {          
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

# color code the canidate SNPs

col.pred[grep("V",col.pred)] <- '#f1eef6'
# add a color to the non-canidate SNPs

empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#e31a1c', '#ffa500','#a6cee3','#2b94db')
# more visualizaton prep


#####################
##### Plot SNPs #####
#####################

plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25)  #displays predictors

legend("topleft", legend= c("AMT", "TS", "PWQ", "PCQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend

##################
# axis 1 & 3
plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices = c(1,3))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices = c(1,3))
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices = c(1,3))
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25, choices = c(1,3))  #displays predictors

legend("bottomright", legend= c("AMT", "TS", "PWQ", "PCQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend

##################
# axis 1 & 4
plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices = c(1,4))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices = c(1,4))
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices = c(1,4))
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25, choices = c(1,4))  #displays predictors

legend("bottomright", legend= c("AMT", "TS", "PWQ", "PCQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend
```

![SNP_plot](https://user-images.githubusercontent.com/54188848/89064209-92eeb180-d337-11ea-9e63-5262a511f150.jpeg)

Again, we interpret this plot like the individual plot. The colored SNPs are our candidates, while the gray SNPs are neutral. Note that some of the candidates do not load strongly onto RDA axis 1 or 2. We expect that those candidates load strongly on another of the RDA axes. 

## References

Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience, 4(1), s13742-015.

Forester, B. R., Lasky, J. R., Wagner, H. H., & Urban, D. L. (2018). Comparing methods for detecting multilocus adaptation with multivariate genotype–environment associations. Molecular Ecology, 27(9), 2215-2233.

Rellstab, C., Gugerli, F., Eckert, A. J., Hancock, A. M., & Holderegger, R. (2015). A practical guide to environmental association analysis in landscape genomics. Molecular Ecology, 24(17), 4348-4370.

Rochette, N. C., Rivera‐Colón, A. G., & Catchen, J. M. (2019). Stacks 2: Analytical methods for paired‐end sequencing improve RADseq‐based population genomics. Molecular ecology, 28(21), 4737-4754.
