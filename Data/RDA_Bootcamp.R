install.packages(c("psych","vegan", "dplyr"), dependencies=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



library(psych)
library(vegan)
library(LEA)

###############################################
##### Bring in files and check formatting #####
###############################################

Genotypes<-read.lfmm("Platyrhinos_subset.lfmm") # Read in genetic data
EnvVars<-read.csv("Platyrhinos_EnvData.csv") # Read in environmental variables

# Inspect the data, look for any missing data
dim(Genotypes)
sum(is.na(Genotypes))

# Look at the structure of your environmental data and change sample names to characters
str(EnvVars)
EnvVars$Samples<-as.character(EnvVars$Samples)


Samples<-EnvVars$Samples
row.names(Genotypes)<- Samples
row.names(Genotypes)
identical(rownames(Genotypes), EnvVars$Samples)

################################################################
##### Check for Correlation between Environmental Variables#####
################################################################

# Visualize correlation between variables
corPlot(EnvVars [,5:23], xlas =  2, cex.axis = 0.75)

# Pull out the selected variables and rename the columns so they are easier to read
SelectedVars<-select(EnvVars,5,8,22,23)
colnames(SelectedVars) <- c("AMT","TS","PWQ","PCQ")

# Check again to make sure our variables are not correlated
pairs.panels(SelectedVars, scale = TRUE)

########################
##### Run your RDA #####
########################

# Run the RDA and examine it
HL_rda<- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)
HL_rda

Rsquared<-RsquareAdj(HL_rda)
Rsquared

# provides a summary of the eigenvalues for each axis
summary(eigenvals(HL_rda, model = "constrained"))

# test the significance of the rda
Fullsig <- anova.cca(HL_rda, parallel=getOption("mc.cores"))
Fullsig


# find which constrained axes are significant
Axissig <- anova.cca(HL_rda, by="axis", parallel=getOption("mc.cores"))
Axissig

# checks the variance inflation factors
# the lower the better
vif.cca(HL_rda)


#############################
# lets make some cool plots

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


################################
##### Testing for Outliers #####
################################

# writes a function identify outliers
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


# extract SNP loadings from the significant axes
load.rda <- scores(HL_rda, choices=c(1:4), display="species") 

# plot histograms of the loadings
hist(load.rda[,1], main="Loadings on Axis 1", xlab = NULL)
hist(load.rda[,2], main="Loadings on Axis 2", xlab = NULL)
hist(load.rda[,3], main="Loadings on Axis 3", xlab = NULL)
hist(load.rda[,4], main="Loadings on Axis 4", xlab = NULL)

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


###################################
##### Organize into Dataframe #####
###################################

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

#################################
##### Investigate Canidates #####
#################################

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

###############################
##### Let's Plot the SNPs #####
###############################

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


# !!!!!!!!!!!!!!!!!!!!!!! Your Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

