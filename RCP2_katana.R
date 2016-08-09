#library(Matrix, lib.loc="../R/x86_64-unknown-linux-gnu-library/3.1/")
#library(Matrix, lib.loc="../R/x86_64-pc-linux-gnu-library/3.2/")
library(methods)
library(Matrix)
library(glmnet)
library(RCPmod)

# data prep ---------------------------------------------------------------
## load data that is prepped from "RCP2_dataPrep.R"
## alternatively, load or create a single object that is site by species/covariates

# covariates.species = read.csv("covariates_species.csv", 
#                               header=TRUE, stringsAsFactors=FALSE)

# loading an .RData file, since uncompressed .csv is a little unwieldly in this case
load("covariates_species_RCP2.RData")

# parameratise model ------------------------------------------------------
arr.job = TRUE # is this an array job or a single job
if (arr.job) {
  # get job number from array job env. var.
  job = Sys.getenv("PBS_ARRAYID")
  job = as.numeric(job)
}

# nominal nRCP if not array job
# nRCP = 3

# should the models be fit with species formula?
# should be one of "species", "nospecies", "both"
# should be passed from the PBS batch job
if ("both" %in% commandArgs()) {species.model = "both"}
if ("species" %in% commandArgs()) {species.model = "species"}
if ("nospecies" %in% commandArgs()) {species.model = "nospecies"}


# use all communities?
subset.data = FALSE # specify T/F to subset data
subset.size = 1000 # specifcy random subset size

# conditional processing to subset data to subset - ensure sitename column is correctly specified
if (subset.data) {
  # specify the  subset to use in the analysis
  set.seed(subset.size)
  sample.sites = sample(covariates.species$SiteNo, subset.size)
  covariates.species = covariates.species[covariates.species$SiteNo %in% sample.sites,]
  print(paste0("Successfully subsetted [",subset.size,"] random sites"))
} else {
  print("No subsetting performed")
}

# record site order
site.names = covariates.species$SiteNo

# define nRCP (number of communities/RCPs)
if (arr.job) {
#   # make the first nRCP to test = 2 (can't solve nRCP=1)
#   nRCP = job+1
  # temp hack code to run multiple starts
  starts = rep(c(2:20),50)
  nRCP = starts[job]
}

# define where abundance data starts in covariates.species, test with names(covariates.species)[n.abund]
n.abund = 19

# choose species to model with
# create a list of species with occurance > n to use for modelling
species.n = 50 # change this to desired minimum occurance
species.count = data.frame(count=sort(colSums(covariates.species[,n.abund:ncol(covariates.species)]), decreasing=T))
species.count$species = row.names(species.count)
model.species.vector = species.count$species[species.count$count>species.n]
model.species.string = paste0("cbind(", paste(model.species.vector, collapse=","),")")

# choose which covariates to use
model.covariates.string = "tempmtcp.1+tempmtcp.2+precipseas.1+precipseas.2+precipann.1+precipann.2+tempmtwp.1+tempmtwp.2+rough500.1+rough500.2+bd200.1+bd200.2+ph200.1+ph200.2"
model.covariates.vector = c("tempmtcp","precipseas","precipann","tempmtwp","rough500","bd200","ph200")

covar.data = covariates.species[,model.covariates.vector]
# calculate quadratic polynomial cols
covar.data = data.frame(poly(covar.data$tempmtcp, 2),
                        poly(covar.data$precipseas, 2),
                        poly(covar.data$precipann, 2),
                        poly(covar.data$tempmtwp, 2),
                        poly(covar.data$rough500, 2),
                        poly(covar.data$bd200, 2),
                        poly(covar.data$ph200, 2))
names(covar.data) = c("tempmtcp.1","tempmtcp.2",
                      "precipseas.1","precipseas.2",
                      "precipann.1","precipann.2",
                      "tempmtwp.1","tempmtwp.2",
                      "rough500.1","rough500.2",
                      "bd200.1","bd200.2",
                      "ph200.1","ph200.2")


## convert categorical variables to factors?

# generate model data
model.data = data.frame(covariates.species[,model.species.vector], covar.data)

# define model form
RCP.form = paste0(model.species.string,"~","1","+",model.covariates.string)

# add column containing factor for form.spp
# e.g. model.data$observer = as.factor(covariates.species$Observers)
model.data$score.method = as.factor(covariates.species$Species.score.method)
model.data$date.int = scale(as.integer(covariates.species$Date), center=T, scale=T)

# define species form
# e.g. "~Observer"
species.form = "~score.method+date.int"

# clear unused variables
rm(covariates.species, species.count)
gc()



# fit mixture models ------------------------------------------------------
my.cont = list(maxit=3000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)

if (species.model %in% c("both", "species")) {
  tic = proc.time()
  fit.regi = regimix(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=nRCP, 
                     dist="Bernoulli", control=my.cont, inits="noPreClust", titbits=TRUE)
  toc = proc.time()
  
  # write model fit stats
  modelStats=list(sites=site.names, covariates=model.covariates.vector, species=model.species.vector, 
                  SppMin=species.n, SppN=fit.regi$S, nRCP=fit.regi$nRCP, runtime=round((toc-tic)[3]/60),
                  AIC=fit.regi$AIC, BIC=fit.regi$BIC, postProbs=fit.regi$postProbs, logl=fit.regi$logl, 
                  coefs=fit.regi$coefs, species.form=species.form, penalties=unlist(my.cont), conv=fit.regi$conv)
  save(modelStats, file=paste0("results/spec/RegimixStats.n",fit.regi$n,
                               ".rcp",fit.regi$nRCP,".s",fit.regi$S,round(fit.regi$logl),".RData"))
}

if (species.model %in% c("both", "nospecies")) {
  tic = proc.time()
  fit.regi = regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=nRCP, 
                     dist="Bernoulli", control=my.cont, inits="noPreClust", titbits=TRUE)
  toc = proc.time()
  
  # write model fit stats
  modelStats=list(sites=site.names, covariates=model.covariates.vector, species=model.species.vector, 
                  SppMin=species.n, SppN=fit.regi$S, nRCP=fit.regi$nRCP, runtime=round((toc-tic)[3]/60),
                  AIC=fit.regi$AIC, BIC=fit.regi$BIC, postProbs=fit.regi$postProbs, logl=fit.regi$logl, 
                  coefs=fit.regi$coefs, species.form=species.form, penalties=unlist(my.cont), conv=fit.regi$conv)
  save(modelStats, file=paste0("results/nospec/RegimixStats.n",fit.regi$n,
                               ".rcp",fit.regi$nRCP,".s",fit.regi$S,round(fit.regi$logl),".RData"))
}