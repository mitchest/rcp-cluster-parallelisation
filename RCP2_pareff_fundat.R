## -> need to run main analysis script to do the load (*_katana.R)
load("covariates_species_RCP2.RData")
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
# # calculate quadratic polynomial cols
# covar.data = data.frame(poly(covar.data$tempmtcp, 2),
#                         poly(covar.data$precipseas, 2),
#                         poly(covar.data$precipann, 2),
#                         poly(covar.data$tempmtwp, 2),
#                         poly(covar.data$rough500, 2),
#                         poly(covar.data$bd200, 2),
#                         poly(covar.data$ph200, 2))
# names(covar.data) = c("tempmtcp.1","tempmtcp.2",
#                       "precipseas.1","precipseas.2",
#                       "precipann.1","precipann.2",
#                       "tempmtwp.1","tempmtwp.2",
#                       "rough500.1","rough500.2",
#                       "bd200.1","bd200.2",
#                       "ph200.1","ph200.2")
# ## convert categorical variables to factors?
# # generate model data
# model.data = data.frame(covariates.species[,model.species.vector], covar.data)
# # define model form
# RCP.form = paste0(model.species.string,"~","1","+",model.covariates.string)

# # add column containing factor for form.spp
# # e.g. model.data$observer = as.factor(covariates.species$Observers)
# covar.data$score.method = as.factor(covariates.species$Species.score.method)
# covar.data$date.int = scale(as.integer(covariates.species$Date), center=T, scale=T)

rm(covariates.species, model.covariates.vector, model.covariates.vector, model.covariates.string, n.abund, species.count, model.species.string, model.species.vector)

# functions defs ----------------------------------------------------------

colMedian <- function(x) {
  meds <- numeric(ncol(x))
  for (i in 1:ncol(x)){
    meds[i] <- median(x[,i])
  }
  return(meds)
}

colMin <- function(x) {
  mins <- numeric(ncol(x))
  for (i in 1:ncol(x)){
    mins[i] <- min(x[,i])
  }
  return(mins)
}

ParPlotIndiv.regimix = function(fm, covar.data, newdata, variable, RCPs) {
  var.range = seq(min(covar.data[,variable]),max(covar.data[,variable]),length.out=1000)
  newdata[,grep(variable, colnames(newdata))] = predict(poly(covar.data[,variable],2),newdata=var.range)
  # do prediction
  preds = predict.regimix(fm, newdata=newdata, nboot=0)
  # plot it
  for (i in RCPs) {
    matplot(var.range, preds[,i], type='l', ylab=paste0("RCP",i), 
            main=paste0("Partial effect of ",variable), xlab=variable)
  }
}

ParPlotMany.regimix = function(fm, covar.data, newdata, variable, RCPs) {
  var.range = seq(min(covar.data[,variable]),max(covar.data[,variable]),length.out=1000)
  newdata[,grep(variable, colnames(newdata))] = predict(poly(covar.data[,variable],2),newdata=var.range)
  # do prediction
  preds = predict.regimix(fm, newdata=newdata, nboot=0)
  # plot it
  matplot(var.range, preds[,RCPs], type='l', lty=rep(1:3,length.out=length(RCPs)), col=RCPs, lwd=3,  
          ylab=paste0("Probability"), main=paste0("Partial effect of ",variable), xlab=variable)
  legend("topleft", legend=paste0("RCP",RCPs), lty=rep(1:3,length.out=length(RCPs)),
         col=RCPs, lwd=3, horiz=FALSE)
}