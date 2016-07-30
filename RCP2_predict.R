#library(raster)
library(RCPmod)
library(Cairo)
library(foreign)
library(data.table)



setwd("A:/1_UNSW/floristic/RCP2")

# file/model code
#model.code = "n4813.s499" # working vcov


# load coefs from best fitted model refit ---------------------------------

# files = list.files(path="results/largeData_NoSpeciesModel/0.0001_10/",
#                    full.names=T,pattern="*.rcp12.*")
# for (i in 1:length(files)) {
#   load(files[i])
#   print(files[i])
#   print(round(min(colSums(modelStats$postProbs)),4))
# }


## --> run anlaysis script to load species/covar data

# load model WITH species params 
load("results/V2_species/RegimixStats.n4715.rcp6.s488-228351.RData")
my.cont = list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
nRCP = 6
params = unlist(modelStats$coefs)
fit.regi.sp = regimix(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=nRCP,
                      dist="Bernoulli", control=my.cont, inits=params, titbits=TRUE)
save(fit.regi.sp, file="predict/fit.regi.sp.RData")

# load model WITHOUT species params, at same number of RCPs as best model WITH species params...
load("results/V2_NoSpecies/RegimixStats.n4715.rcp6.s488-239712.RData")
my.cont = list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
nRCP = 6
params = unlist(modelStats$coefs)
fit.regi.nospsp = regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=nRCP,
                      dist="Bernoulli", control=my.cont, inits=params, titbits=TRUE)
save(fit.regi.nospsp, file="predict/fit.regi.nosp-at-sp.RData")

# load model WITHOUT species params
load("results/V2_NoSpecies/RegimixStats.n4715.rcp8.s488-234104.RData")
my.cont = list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
nRCP = 8
params = unlist(modelStats$coefs)
fit.regi.nosp = regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=nRCP,
                        dist="Bernoulli", control=my.cont, inits=params, titbits=TRUE)
save(fit.regi.nosp, file="predict/fit.regi.nosp.RData")



# make  prediction to env/geo space ---------------------------------------

# load covariates
covars = read.dbf("covars/NSW_2km.dbf", as.is=T)
# # convert categorical variables to factors
# covars$soil_gsg = as.factor(covars$soil_gsg)

# check which vars needed
attr(fit.regi.sp$titbits$X, "dimnames")[[2]]

# make covar data and location data
locs = covars[,c("lat","long")]

# get poly() coefficients from model fit and predict to newdata
model.covars = covariates.species[,model.covariates.vector]

covars = data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=covars$tempmtcp),
                    predict(poly(model.covars$precipseas, 2), newdata=covars$precipseas),
                    predict(poly(model.covars$precipann, 2), newdata=covars$precipann),
                    predict(poly(model.covars$tempmtwp, 2), newdata=covars$tempmtwp),
                    predict(poly(model.covars$rough500, 2), newdata=covars$rough500),
                    predict(poly(model.covars$bd200, 2), newdata=covars$bd200),
                    predict(poly(model.covars$ph200, 2), newdata=covars$ph200))
names(covars) = c("tempmtcp.1","tempmtcp.2",
                  "precipseas.1","precipseas.2",
                  "precipann.1","precipann.2",
                  "tempmtwp.1","tempmtwp.2",
                  "rough500.1","rough500.2",
                  "bd200.1","bd200.2",
                  "ph200.1","ph200.2")

# need to choose the species model levels
covars$score.method = model.data$score.method[1]

scale.att = attributes(scale(as.integer(covariates.species$Date), center=T, scale=T))
covars$date.int = as.numeric(scale(as.integer(as.Date("2015-01-01")),
                                   center=scale.att$`scaled:center`, scale=scale.att$`scaled:scale`))



# combine regiboot objects ------------------------------------------------

# the regiboot call should be farmed out to a cluster if it's for a larger model...

nboot = 5 # this is the number of bootstrap iterations done when processing is farmed out 

boot.nosp = list.files(path="A:/1_UNSW/floristic/RCP2/predict/bootstrap/",
                       pattern="\\.nosp\\.", full.names=T)
boot.nosp.rbind = list()
#for (i in 1:length(boot.nosp)) {
for (i in 1:length(boot.nosp)) {
  load(boot.nosp[i])
  boot.nosp.rbind[[i]] = data.frame(matrix(fit.regiboot.nosp, nrow=nboot))
  #names(boot.nosp.rbind[[i]]) = attributes(fit.regiboot.nosp)[[2]][[2]]
}
fit.regiboot.nosp = as.matrix(rbindlist(boot.nosp.rbind))
class(fit.regiboot.nosp) = "regiboot"
save(fit.regiboot.nosp, file="fit.regiboot.nosp.RData")

boot.nospsp = list.files(path="A:/1_UNSW/floristic/RCP2/predict/bootstrap/",
                       pattern="\\.nosp-at-sp\\.", full.names=T)
boot.nospsp.rbind = list()
#for (i in 1:length(boot.nosp)) {
for (i in 1:length(boot.nospsp)) {
  load(boot.nospsp[i])
  boot.nospsp.rbind[[i]] = data.frame(matrix(fit.regiboot.nospsp, nrow=nboot))
  #names(boot.nosp.rbind[[i]]) = attributes(fit.regiboot.nosp)[[2]][[2]]
}
fit.regiboot.nospsp = as.matrix(rbindlist(boot.nospsp.rbind))
class(fit.regiboot.nospsp) = "regiboot"
save(fit.regiboot.nospsp, file="fit.regiboot.nosp-at-sp.RData")

boot.sp = list.files(path="A:/1_UNSW/floristic/RCP2/predict/bootstrap/",
                       pattern="\\.sp\\.", full.names=T)
boot.sp.rbind = list()
#for (i in 1:length(boot.nosp)) {
for (i in 1:length(boot.nosp)) {
  load(boot.sp[i])
  boot.sp.rbind[[i]] = data.frame(matrix(fit.regiboot.sp, nrow=nboot))
  #names(boot.sp.rbind[[i]]) = attributes(fit.regiboot.sp)[[2]][[2]]
}
fit.regiboot.sp = as.matrix(rbindlist(boot.sp.rbind))
class(fit.regiboot.sp) = "regiboot"
save(fit.regiboot.sp, file="fit.regiboot.sp.RData")



# predict out to covariates -----------------------------------------------

#################
## point preds ##
#################
# predict for newdata - just point predictions
load("predict/fit.regi.sp.RData")
predicted.sp = predict.regimix(object=fit.regi.sp, newdata=covars, nboot=0)
predicted.sp.estimate = data.frame(locs, predicted.sp)
save(predicted.sp.estimate, file="predict/predicted.sp.estimate.RData")

load("predict/fit.regi.nosp-at-sp.RData")
predicted.nospsp = predict.regimix(object=fit.regi.nospsp, newdata=covars, nboot=0)
predicted.nospsp.estimate = data.frame(locs, predicted.nospsp)
save(predicted.nospsp.estimate, file="predict/predicted.nosp-at-sp.estimate.RData")

load("predict/fit.regi.nosp.RData")
predicted.nosp = predict.regimix(object=fit.regi.nosp, newdata=covars, nboot=0)
predicted.nosp.estimate = data.frame(locs, predicted.nosp)
save(predicted.nosp.estimate, file="predict/predicted.nosp.estimate.RData")

# dummy ci's until boots done
predicted.sp.low = data.frame(locs, matrix(nrow=nrow(predicted.sp), ncol=ncol(predicted.nosp)))
predicted.sp.up = data.frame(locs, matrix(nrow=nrow(predicted.sp), ncol=ncol(predicted.nosp)))
predicted.nosp.low = data.frame(locs, matrix(nrow=nrow(predicted.nosp), ncol=ncol(predicted.nosp)))
predicted.nosp.up = data.frame(locs, matrix(nrow=nrow(predicted.nosp), ncol=ncol(predicted.nosp)))

###################
## bootstrapping ##
###################

# boot strap sample for parameter variance (this nees to be farmed out to multiple cores for large models)
#fit.regiboot.sp = regiboot(fit.regi.sp, nboot=100)
#fit.regiboot.nosp = regiboot(fit.regi.nosp, nboot=100)

# first look at variance on model data
fit.regivar.sp = predict.regimix(object=fit.regi.sp, object2=fit.regiboot.sp)
fit.regivar.nosp = predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp)

# predict for newdata with variances (split up newdata for memory issues)
predicted.sp.list = list()
predicted.nosp.list = list()
predicted.nospsp.list = list()
idx = round(nrow(covars)/2)
for (i in 1:2) {
  if (i==1) {
#     predicted.sp.list[[i]] = predict.regimix(object=fit.regi.sp, object2=fit.regiboot.sp, 
#                                              newdata=covars[1:idx,])
#     gc()
#     predicted.nosp.list[[i]] = predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp, 
#                                                newdata=covars[1:idx,])
#     gc()
    predicted.nospsp.list[[i]] = predict.regimix(object=fit.regi.nospsp, object2=fit.regiboot.nospsp, 
                                               newdata=covars[1:idx,])
    gc()
  }
  if (i==2) {
#     predicted.sp.list[[i]] = predict.regimix(object=fit.regi.sp, object2=fit.regiboot.sp, 
#                                              newdata=covars[(idx+1):nrow(covars),])
#     gc()
#     predicted.nosp.list[[i]] = predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp, 
#                                                newdata=covars[(idx+1):nrow(covars),])
#     gc()
    predicted.nospsp.list[[i]] = predict.regimix(object=fit.regi.nospsp, object2=fit.regiboot.nospsp, 
                                               newdata=covars[(idx+1):nrow(covars),])
    gc()
  }
}


# add location data back
# predicted.ptPreds = data.frame(locs, predicted$ptPreds)
# predicted.bootPred = data.frame(locs, predicted$bootPreds)
# predicted.SE = data.frame(locs, predicted$bootSEs)
predicted.sp.btpreds = data.frame(locs, rbind(predicted.sp.list[[1]]$bootPreds[,],
                                          predicted.sp.list[[2]]$bootPreds[,]))
predicted.sp.low = data.frame(locs, rbind(predicted.sp.list[[1]]$bootCIs[,,1],
                                            predicted.sp.list[[2]]$bootCIs[,,1]))
predicted.sp.up = data.frame(locs, rbind(predicted.sp.list[[1]]$bootCIs[,,2],
                                           predicted.sp.list[[2]]$bootCIs[,,2]))

predicted.nosp.btpreds = data.frame(locs, rbind(predicted.nosp.list[[1]]$bootPreds[,],
                                              predicted.nosp.list[[2]]$bootPreds[,]))
predicted.nosp.low = data.frame(locs, rbind(predicted.nosp.list[[1]]$bootCIs[,,1],
                                            predicted.nosp.list[[2]]$bootCIs[,,1]))
predicted.nosp.up = data.frame(locs, rbind(predicted.nosp.list[[1]]$bootCIs[,,2],
                                           predicted.nosp.list[[2]]$bootCIs[,,2]))

predicted.nospsp.btpreds = data.frame(locs, rbind(predicted.nospsp.list[[1]]$bootPreds[,],
                                                predicted.nospsp.list[[2]]$bootPreds[,]))
predicted.nospsp.low = data.frame(locs, rbind(predicted.nospsp.list[[1]]$bootCIs[,,1],
                                            predicted.nospsp.list[[2]]$bootCIs[,,1]))
predicted.nospsp.up = data.frame(locs, rbind(predicted.nospsp.list[[1]]$bootCIs[,,2],
                                           predicted.nospsp.list[[2]]$bootCIs[,,2]))

save(predicted.sp.low, file="predict/predicted.sp.low.RData")
save(predicted.sp.up, file="predict/predicted.sp.up.RData")
save(predicted.sp.btpreds, file="predict/predicted.sp.btpreds.RData")
save(predicted.nosp.low, file="predict/predicted.nosp.low.RData")
save(predicted.nosp.up, file="predict/predicted.nosp.up.RData")
save(predicted.nosp.btpreds, file="predict/predicted.nosp.btpreds.RData")
save(predicted.nospsp.low, file="predict/predicted.nosp-at-sp.low.RData")
save(predicted.nospsp.up, file="predict/predicted.nosp-at-sp.up.RData")
save(predicted.nospsp.btpreds, file="predict/predicted.nosp-at-sp.btpreds.RData")

# save off
# model.code_ = gsub(model.code, pattern="\\.", replacement="_")
# 
# write.csv(predicted.sp.estimate, file=paste0("predict/RCP_sp_ptpreds_",model.code_,".csv"), row.names=F)
# write.csv(predicted.sp.bootPred, file=paste0("predict/RCP_sp_bootpreds_",model.code_,".csv"), row.names=F)
# write.csv(predicted.sp.SE, file=paste0("predict/RCP_sp_SEs_",model.code_,".csv"), row.names=F)
# write.csv(predicted.sp.low, file=paste0("predict/RCP_sp_lowerCI_",model.code_,".csv"), row.names=F)
# write.csv(predicted.sp.up, file=paste0("predict/RCP_sp_upperCI_",model.code_,".csv"), row.names=F)
# 
# write.csv(predicted.nosp.estimate, file=paste0("predict/RCP_nosp_ptpreds_",model.code_,".csv"), row.names=F)
# write.csv(predicted.nosp.bootPred, file=paste0("predict/RCP_nosp_bootpreds_",model.code_,".csv"), row.names=F)
# write.csv(predicted.nosp.SE, file=paste0("predict/RCP_nosp_SEs_",model.code_,".csv"), row.names=F)
# write.csv(predicted.nosp.low, file=paste0("predict/RCP_nosp_lowerCI_",model.code_,".csv"), row.names=F)
# write.csv(predicted.nosp.up, file=paste0("predict/RCP_nosp_upperCI_",model.code_,".csv"), row.names=F)
# 
# 
# # plot residuals
# plot.regimix(x=fit.regi.sp)
# plot.regimix(x=fit.regi.nosp)

# save for plotting
save(predicted.nosp.up, predicted.nosp.low, predicted.nosp.estimate, predicted.nosp.btpreds, file="PredsWithoutSpeciesTerms.RData")
save(predicted.nospsp.up, predicted.nospsp.low, predicted.nospsp.estimate, predicted.nospsp.btpreds, file="PredsWithoutSpeciesTerms-at-sp.RData")
save(predicted.sp.up, predicted.sp.low, predicted.sp.estimate, predicted.sp.btpreds, file="PredsWithSpeciesTerms.RData")




# plot rasters ------------------------------------------------------------

library(raster)

load("predict/predicted.sp.estimate.RData")
load("predict/predicted.sp.up.RData")
load("predict/predicted.sp.low.RData")
load("predict/predicted.nosp.estimate.RData")
load("predict/predicted.nosp.up.RData")
load("predict/predicted.nosp.low.RData")
load("predict/predicted.nosp-at-sp.estimate.RData")
load("predict/predicted.nosp-at-sp.up.RData")
load("predict/predicted.nosp-at-sp.low.RData")

# plot points preds and 95% CIs for each RCP

## species model
#CairoWin()
CairoPDF(file="NSWVeg_SpeciesModel_preds.pdf", height=7, width=6)
par(mar=c(0.5,4,1,0.5), oma=c(2,4,2,0), mfrow=c(6,3))
#ext = extent(-35,-30,144,150.5)
#colour = c("#ffffff","#f5f5f5","#fcae91","#fb6a4a","#de2d26","#a50f15", "#a50f15")
colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
breaks = c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

for (i in 3:8) {
  if (i==3) {
    plot(rasterFromXYZ(predicted.sp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2), main="Lower CI")
    plot(rasterFromXYZ(predicted.sp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Point Prediction")
    plot(rasterFromXYZ(predicted.sp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Upper CI")
  } else if (i>3 & i<8) {
    plot(rasterFromXYZ(predicted.sp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.sp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
    plot(rasterFromXYZ(predicted.sp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
  } else {
    plot(rasterFromXYZ(predicted.sp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.sp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
    plot(rasterFromXYZ(predicted.sp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
  }
}
dev.off()

## no species model
#CairoWin()
CairoPDF(file="NSWVeg_NoSpeciesModel_preds.pdf", height=8, width=5.5)
par(mar=c(0.5,4,1,0.5), oma=c(2,4,2,0), mfrow=c(8,3))
#ext = extent(-35,-30,144,150.5)
#colour = c("#ffffff","#f5f5f5","#fcae91","#fb6a4a","#de2d26","#a50f15", "#a50f15")
colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
breaks = c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

for (i in 3:10) {
  if (i==3) {
    plot(rasterFromXYZ(predicted.nosp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2), main="Lower CI")
    plot(rasterFromXYZ(predicted.nosp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Point Prediction")
    plot(rasterFromXYZ(predicted.nosp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Upper CI")
  } else if (i>3 & i<10) {
    plot(rasterFromXYZ(predicted.nosp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.nosp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
    plot(rasterFromXYZ(predicted.nosp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
  } else {
    plot(rasterFromXYZ(predicted.nosp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.nosp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
    plot(rasterFromXYZ(predicted.nosp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
  }
}
dev.off()

# make one plot for the colour ramp
CairoWin()
colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
plot(rasterFromXYZ(cbind(1:100,1:100,runif(100,0,1))), breaks=breaks, col=colour, legend.only=T)


## no species model at nRCP=best model WITH speices terms
#CairoWin()
CairoPDF(file="NSWVeg_NoSpeciesModel-at-sp_preds.pdf", height=7, width=6)
par(mar=c(0.5,4,1,0.5), oma=c(2,4,2,0), mfrow=c(6,3))
#ext = extent(-35,-30,144,150.5)
#colour = c("#ffffff","#f5f5f5","#fcae91","#fb6a4a","#de2d26","#a50f15", "#a50f15")
colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
breaks = c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

for (i in 3:10) {
  if (i==3) {
    plot(rasterFromXYZ(predicted.nospsp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2), main="Lower CI")
    plot(rasterFromXYZ(predicted.nospsp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Point Prediction")
    plot(rasterFromXYZ(predicted.nospsp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F, main="Upper CI")
  } else if (i>3 & i<10) {
    plot(rasterFromXYZ(predicted.nospsp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, xaxt="n", ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.nospsp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
    plot(rasterFromXYZ(predicted.nospsp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, axes=F)
  } else {
    plot(rasterFromXYZ(predicted.nospsp.low[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, ylab=paste0("RCP",i-2))
    plot(rasterFromXYZ(predicted.nospsp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
    plot(rasterFromXYZ(predicted.nospsp.up[,c(2,1,i)]), breaks=breaks, col=colour, legend=F, yaxt="n")
  }
}
dev.off()

# make one plot for the colour ramp
CairoWin()
colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
plot(rasterFromXYZ(cbind(1:100,1:100,runif(100,0,1))), breaks=breaks, col=colour, legend.only=T)



# # plot just point preds
# CairoWin()
# #CairoPDF(file="SpeciesModel_RCP5_ptPreds.pdf", height=6,  width=9)
# par(mar=c(0.5,0.5,1,0.5), oma=c(2,2,2,0), mfrow=c(2,3))
# colour = c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d","#000000")
# breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# for (i in 3:7) {
#   plot(rasterFromXYZ(predicted.sp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, 
#        legend=F, axes=F, main=paste0("RCP ", i-2))
# }
# plot(rasterFromXYZ(cbind(1:100,1:100,runif(100,0,1))), breaks=breaks, col=colour, legend.only=T)
# #dev.off()
# 
# CairoWin()
# #CairoPDF(file="NoSpeciesModel_RCP12_ptPreds.pdf", height=12,  width=9)
# par(mar=c(0.5,0.5,1,0.5), oma=c(2,2,2,0), mfrow=c(3,3))
# for (i in 3:9) {
#   plot(rasterFromXYZ(predicted.nosp.estimate[,c(2,1,i)]), breaks=breaks, col=colour, 
#        legend=F, axes=F, main=paste0("RCP ", i-2))
#   if (i == 14) { plot(rasterFromXYZ(cbind(1:100,1:100,runif(100,0,1))), breaks=breaks, col=colour, legend.only=T, add=T) }
# }
# #dev.off()


# plot a hard clustered map
## function that returns hard cluster membership
HardClust = function(ptPreds){
  postProbs = data.frame(round(ptPreds[,c(-1,-2)], 5))
  postProbs.rowMax = apply(postProbs,1,max)
  RCPclusters = integer(nrow(postProbs)) # pre-allocate vector
  for (i in 1:nrow(postProbs)){
    RCPclusters[i] = which(postProbs[i,]==postProbs.rowMax[i])
  }
  return(RCPclusters)
}

ptPreds.hard = HardClust(predicted.nosp.estimate)
ptPreds.hard = data.frame(predicted.nosp.estimate[,c(2,1)],ptPreds.hard)

#SEs = predicted[,c(1,2,21:29)]

CairoWin()
#par(mar=c(1,1,1,1), mfrow=c(1,2))
colour = sample(rainbow(n=8))
plot(rasterFromXYZ(ptPreds.hard), legend=F, col=colour)
legend(x='bottomright', legend = c("RCP 1","RCP 2","RCP 3","RCP 4","RCP 5","RCP 6","RCP 7","RCP 8"),
       fill = colour)

