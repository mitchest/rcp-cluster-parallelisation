library(ggplot2)
library(Cairo)
library(dplyr)
library(RCPmod)



# load data and plot ------------------------------------------------------

#controls = "p0.05_t100_random"
regimix.results = function(path, pattern, controls="", remove.misfit=T, plot.only=F) {
  controls = controls

  files = list.files(path=paste0(path,controls),full.names=T,pattern=pattern)
  nRCP = as.numeric(rep(NA,length(files)))
  BIC = as.numeric(rep(NA,length(files)))
  AIC = as.numeric(rep(NA,length(files)))
  runtime = as.numeric(rep(NA,length(files)))
  minPostProb = as.numeric(rep(NA,length(files)))
  maxPostProb = as.numeric(rep(NA,length(files)))
  logl = as.numeric(rep(NA,length(files)))
  conv = as.numeric(rep(NA,length(files)))

  for (i in 1:length(files)) {
    load(files[i])
    nRCP[i] = modelStats$nRCP
    BIC[i] = modelStats$BIC
    AIC[i] = modelStats$AIC
    runtime[i] = modelStats$runtime
    logl[i] = modelStats$logl
    minPostProb[i] = min(colSums(modelStats$postProbs))
    maxPostProb[i] = max(colSums(modelStats$postProbs))
    conv[i] = modelStats$conv
    rm(modelStats)
  }

  nRCP.plot = data.frame(nRCP, BIC, AIC, runtime, logl, minPP=round(minPostProb, 3), maxPP=round(maxPostProb, 3), conv)
  nRCP.plot = nRCP.plot[!is.na(nRCP),]

  if (remove.misfit) {
    nRCP.plot = nRCP.plot[nRCP.plot$minPP>0,]
    nRCP.plot = nRCP.plot[nRCP.plot$conv==0,]
  }
  if (plot.only) {
    return(nRCP.plot[,c("nRCP","BIC","logl")])
  } else {
    return(nRCP.plot)
  }
}

plot.regimix.results = function(regimix.results, model.type) {
  print(plot(regimix.results$AIC~regimix.results$nRCP, main=paste0(model.type)))
  print(abline(v=regimix.results$nRCP[regimix.results$AIC==min(regimix.results$AIC)], lty=2))
  print(plot(regimix.results$BIC~regimix.results$nRCP, main=paste0("BIC")))
  print(abline(v=regimix.results$nRCP[regimix.results$BIC==min(regimix.results$BIC)], lty=2))
  print(plot(regimix.results$runtime/60~regimix.results$nRCP, main="runtime (hours)"))
  print(plot(log(regimix.results$minPP)~regimix.results$nRCP, main="log(min(colSums(postProbs))))"))
  print(plot(log(regimix.results$maxPP)~regimix.results$nRCP, main="log(max(colSums(postProbs)))"))
}


# paper results
results.NoSpeciesModel_.0001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/V3_NoSpecies/",
                                               "*.n4715.*.s488*", remove.misfit=T)
results.SpeciesModel_.0001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/V3_species/",
                                             "*.n4715.*.s488*", remove.misfit=T)

for (i in 1:20) {print(paste0(i,": ",sum(results.SpeciesModel_.0001$nRCP==i)))}
for (i in 1:20) {print(paste0(i,": ",sum(results.NoSpeciesModel_.0001$nRCP==i)))}

save(results.SpeciesModel_.0001, results.NoSpeciesModel_.0001, file="plotcodedata/modelresults.RData")

load("plotcodedata/modelresults.RData")
results.NoSpeciesModel_.0001 <- results.NoSpeciesModel_.0001[c("nRCP","BIC")]
results.SpeciesModel_.0001 <- results.SpeciesModel_.0001[c("nRCP","BIC")]
save(results.SpeciesModel_.0001, results.NoSpeciesModel_.0001, file="plotcodedata/BICresults.RData")

# paper figure
BICmin.sp = numeric(19)
BICmin.nosp = numeric(19)
for (i in 2:20) {
  BICmin.sp[i-1] = min(results.SpeciesModel_.0001$BIC[results.SpeciesModel_.0001$nRCP==i])
  BICmin.nosp[i-1] = min(results.NoSpeciesModel_.0001$BIC[results.NoSpeciesModel_.0001$nRCP==i])
}
#CairoWin()
#CairoPDF(file="NSWVegBIC.pdf", height=9, width=10)
#par(mfrow=c(1,1))
plot(1~1, type='n', ylab="BIC", xlab="nRCP", main="BIC vs. nRCP", ylim=c(492000, 555000), xlim=c(1,30))
points(results.NoSpeciesModel_.0001$BIC ~ c(results.NoSpeciesModel_.0001$nRCP-0.2), pch=16, col="grey", cex=0.5)
points(results.SpeciesModel_.0001$BIC ~ c(results.SpeciesModel_.0001$nRCP+0.2), pch=16, col="coral", cex=0.5)
points(BICmin.nosp ~ c(2:20), pch=16, type='b', col="black")
points(BICmin.sp ~ c(2:20), pch=16, type='b', col="red")
legend("bottomright", legend=c("No species model \n(min BIC)\n", "Species dependence \nmodel (min BIC)",
                               "No species model BIC", "Species dependence \nmodel BIC"),
       pch=c(16,16,16,16), col=c("black","red","grey","coral"))
#dev.off()


# regimix.results.plot = function(){
#   # spend time on this once we figure out a plotting style
# }

# # explore
# results.NoSpeciesModel_.01 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.01_10/",
#                                              "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.NoSpeciesModel_.001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.001_10/",
#                                               "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.NoSpeciesModel_.0001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.0001_10/",
#                                                "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
#
# results.NoSpeciesModel_.1_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.1_1/",
#                                               "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.NoSpeciesModel_.01_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.01_1/",
#                                                "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.NoSpeciesModel_.001_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.001_1/",
#                                                 "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.NoSpeciesModel_.0001_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_NoSpeciesModel/0.0001_1/",
#                                                  "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
#
# results.SpeciesModel_.01 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.01_10/",
#                                            "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.SpeciesModel_.001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.001_10/",
#                                             "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.SpeciesModel_.0001 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.0001_10/",
#                                              "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
#
# results.SpeciesModel_.1_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.1_1/",
#                                             "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.SpeciesModel_.01_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.01_1/",
#                                              "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.SpeciesModel_.001_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.001_1/",
#                                               "*.n4813.*.s499*", remove.misfit=F, plot.only=F)
# results.SpeciesModel_.0001_1 = regimix.results("A:/1_UNSW/floristic/RCP2/results/largeData_SpeciesModel/method_date_0.0001_1/",
#                                                "*.n4813.*.s499*", remove.misfit=F, plot.only=F)

# # results plots
# CairoPDF(file="n4813.s499.largeData.results.pdf", height=10, width=4)
# par(mfcol=c(5,2))
#
# # no species model
# # plot.regimix.results(results.NoSpeciesModel_.01, "No Spp. model \nk=0.01,tau=10")
# # plot.regimix.results(results.NoSpeciesModel_.001, "No Spp. model \nk=0.001,tau=10")
# plot.regimix.results(results.NoSpeciesModel_.0001, "No Spp. model \nk=0.0001,tau=gamma=10")
#
# # plot.regimix.results(results.NoSpeciesModel_.1_1, "No Spp. model \nk=0.1,tau/gamma=1")
# # plot.regimix.results(results.NoSpeciesModel_.01_1, "No Spp. model \nk=0.01,tau/gamma=1")
# # plot.regimix.results(results.NoSpeciesModel_.001_1, "No Spp. model \nk=0.001,tau/gamma=1")
# # plot.regimix.results(results.NoSpeciesModel_.0001_1, "No Spp. model \nk=0.0001,tau/gamma=1")
#
# # species model
# # plot.regimix.results(results.SpeciesModel_.01, "~method+date \nk=0.01,tau=10")
# # plot.regimix.results(results.SpeciesModel_.001, "~method+date \nk=0.001,tau=10")
# plot.regimix.results(results.SpeciesModel_.0001, "~method+date \nk=0.0001,tau=gamma=10")
#
# # plot.regimix.results(results.SpeciesModel_.1_1, "~method+date \nk=0.1,tau/gamma=1")
# # plot.regimix.results(results.SpeciesModel_.01_1, "~method+date \nk=0.01,tau/gamma=1")
# # plot.regimix.results(results.SpeciesModel_.001_1, "~method+date \nk=0.001,tau/gamma=1")
# # plot.regimix.results(results.SpeciesModel_.0001_1, "~method+date \nk=0.0001,tau/gamma=1")
#
# dev.off()



# # coef plots
# par(mfrow=c(2,3))
# plot(modelStats$coefs$tau, ylab="taus", main="~method+date \nk=0.001,tau=1")
# plot(modelStats$coefs$beta, ylab="betas", main="~method+date \nk=0.001,tau=1")
# plot(modelStats$coefs$gamma, ylab="gammas", main="~method+date \nk=0.001,tau=1")
#
# plot(modelStats$coefs$tau, ylab="taus", main="~method+date \nk=0.01,tau=10")
# plot(modelStats$coefs$beta, ylab="betas", main="~method+date \nk=0.01,tau=10")
# plot(modelStats$coefs$gamma, ylab="gammas", main="~method+date \nk=0.01,tau=10")


# par(mfrow=c(2,2))
# plot(results.NoSpeciesModel$AIC~results.NoSpeciesModel$nRCP, main="AIC")
# points(results.SpeciesModel$AIC~results.SpeciesModel$nRCP, col="red")
# plot(results.NoSpeciesModel$BIC~results.NoSpeciesModel$nRCP, main=paste0("BIC"))
# points(results.SpeciesModel$BIC~results.SpeciesModel$nRCP, col="red")
# plot(results.NoSpeciesModel$runtime/60~results.NoSpeciesModel$nRCP, main="runtime (hours)")
# points(results.SpeciesModel$runtime/60~results.SpeciesModel$nRCP, col="red")
# plot(results.NoSpeciesModel$minPP~results.NoSpeciesModel$nRCP, main="min(colSums(postProbs))")
# points(results.SpeciesModel$minPP~results.SpeciesModel$nRCP, col="red")

# # paper figure
# CairoWin()
# par(mfrow=c(1,2), cex.lab=1.5, oma=c(3,3,0,0)+0.1, mar=c(5,4,0,0)+0.1)
# plot(BIC~nRCP, data=nRCP.plot, ylab="BIC", xlab="number of RCPs", yaxt='n')
# rect(xleft=c(4.5,4.5), xright=c(14.5,14.5), ybottom=c(580000,580000), ytop=c(670000,670000), lty=2, lwd=2)
# plot(BIC~nRCP, data=nRCP.plot[nRCP.plot$nRCP>4 & nRCP.plot$nRCP<15 & nRCP.plot$BIC<800000,],
#      ylab="", xlab="number of RCPs", yaxt='n')
# points(BIC~nRCP, data=nRCP.plot[nRCP.plot$BIC==min(nRCP.plot$BIC),], pch=8, cex=4)



# load specific model for diagnostics -------------------------------------
load("covariates_species_RCP2.RData")
nospec.rcp7 = get(load("results/V2_NoSpecies/RegimixStats.n4715.rcp8.s488-234104.RData"))
spec.rcp5 = get(load("results/V2_species/RegimixStats.n4715.rcp6.s488-228351.RData"))

gammas = spec.rcp5$coefs$gamma
gamma.labs = c(levels(as.factor(covariates.species$Species.score.method))[-1], "Date")
names(gammas) = rep(gamma.labs,each=488)

#CairoWin()
CairoPDF(file="NSWVegRCP5_gammas.pdf", height=6, width=10)
par(mfrow=c(2,3), cex.lab=2.5)
for (i in gamma.labs) {
  hist(gammas[names(gammas)==i], ylab="", xlab=expression(gamma), main=i)
}
dev.off()


# check posterior probabilities -------------------------------------------
# colSums(postProbs)
# plot(colSums(fit.regi$postProbs),
#      ylim=c(0,900), pch=16,
#      xlab="RCP number", ylab="expected number sites")



# confusion matrix --------------------------------------------------------
library(foreign)
load("covariates_species_RCP2.RData")
load("predict/fit.regi.sp.RData")
load("predict/fit.regi.nosp.RData")

keithmap.locs = read.dbf("results/ConfusionMatrix/keithmap_locs.dbf", as.is=T)
keithmap.locs = keithmap.locs[keithmap.locs$SiteNo %in% covariates.species$SiteNo,]
keithmap.locs$RASTERVALU[keithmap.locs$RASTERVALU==-9999] = 1
keithmap.attr = read.dbf("results/ConfusionMatrix/keithmap_attributes.dbf", as.is=T)
classes = inner_join(x=keithmap.locs, y=keithmap.attr, by=c("RASTERVALU"="VALUE"))

# species model
# get RCP posterior probabilities
postProbs = fit.regi.sp$postProbs
# create binary matrix for veg classification
vegCom = model.matrix(~0+FORMATIONN, data=classes) # choose which veg classes to compare i.e. PCT/class/form
# build matrix of expected shared sites
sharedSites = t(postProbs) %*% vegCom
sharedSites = round(sharedSites, 1)
# clean up
rm(postProbs, vegCom)
# check out which classes have many matches with RCPs
sharedSites.df = data.frame(sharedSites)
# write to file
write.csv(sharedSites.df, file="results/ConfusionMatrix/species_sharedsites.csv", row.names=F)

# no species model
postProbs = fit.regi.nosp$postProbs
vegCom = model.matrix(~0+FORMATIONN, data=classes) # choose which veg classes to compare i.e. PCT/class/form
sharedSites = t(postProbs) %*% vegCom
sharedSites = round(sharedSites, 1)
rm(postProbs, vegCom)
sharedSites.df = data.frame(sharedSites)
write.csv(sharedSites.df, file="results/ConfusionMatrix/nospecies_sharedsites.csv", row.names=F)




# stability plots ---------------------------------------------------------

# modified functions so plooting can be done later - run on cluster
# stability.regimix.data <- function (model, oosSizeRange = NULL, times = model$n, mc.cores = 1,
#                                     quiet = FALSE) {
#   if (is.null(oosSizeRange))
#     oosSizeRange <- round(seq(from = 1, to = model$n%/%5,
#                               length = 10))
#   if (any(oosSizeRange < 1))
#     stop("Silly number of RCPs. Specified range is: ", oosSizeRange,
#          " and they should all be >= 1")
#   disty <- matrix(NA, nrow = length(oosSizeRange), ncol = model$nRCP)
#   predlogls <- array(NA, dim = c(length(oosSizeRange), model$n,
#                                  times))
#   for (ii in oosSizeRange) {
#     tmp <- cooks.distance.regimix(model, oosSize = ii, times = times,
#                           mc.cores = mc.cores, quiet = quiet)
#     disty[oosSizeRange == ii, ] <- colMeans(abs(tmp$cooksD))
#     predlogls[oosSizeRange == ii, , ] <- tmp$predLogL
#   }
#   # return stuff needed for plotting
#   return(list(oosSizeRange=oosSizeRange,
#               times=times,
#               disty=disty,
#               predlogls=predlogls,
#               model=list(n=model$n,nRCP=model$nRCP,logl.sites=model$logl.sites)))
# }

stability.regimix.plot <- function (stability.regimix.data) {
  # declare all the plotting variables
  oosSizeRange = stability.regimix.data$oosSizeRange
  times = stability.regimix.data$times
  disty = stability.regimix.data$disty
  predlogls = stability.regimix.data$predlogls
  model = stability.regimix.data$model
  # plotting
  par(mfrow = c(1, 2))
  matplot(c(0, oosSizeRange), rbind(0, disty), type = "b",
          ylab = "Distance from Full Model Predictions", xlab = "Number of Obs Removed",
          main = "Stability of Group Predictions", col = 1:model$nRCP,
          pch = as.character(1:model$nRCP), lty = 1)
  legend("center", bty = "n", lty = 1, pch = as.character(1:model$nRCP),
         col = 1:model$nRCP, legend = paste("RCP ", 1:model$nRCP,
                                            sep = ""))
  plot(rep(oosSizeRange, each = prod(dim(predlogls[1, , ]))),
       predlogls, pch = 20, ylab = "Pred LogL (OOS)", xlab = "Number of Obs Removed",
       main = "Stability of Pred Logl", xlim = c(0, max(oosSizeRange)),
       type = "n")
  # all data
  rbPal <- colorRampPalette(c('powderblue','blue'))
  histo <- hist(model$logl.sites, breaks=100, plot=F)
  breaks <- histo$breaks; counts <- histo$counts; cols <- rbPal(max(counts))[counts]
  points(rep(0, length(breaks)), breaks, pch = 20, col=cols)
  # leave outs
  rbPal <- colorRampPalette(c('peachpuff','red'))
  for (ii in oosSizeRange) {
    # points(rep(ii, prod(dim(predlogls[1, , ]))), predlogls[oosSizeRange ==
    #                                                          ii, , ], pch = 20)
    histo <- hist(predlogls[oosSizeRange==ii,,], breaks=100, plot=F)
    breaks <- histo$breaks; counts <- histo$counts; cols <- rbPal(max(counts))[counts]
    points(rep(ii,length(breaks)), breaks, pch=20, col=cols)
  }
  legend("bottom",title="Density",legend=c("zero","low","med","high"), pch=21, col="black", pt.bg=rbPal(4), horiz=T)
  lines(c(0, oosSizeRange), c(mean(model$logl.sites), apply(predlogls,
                                                            1, mean, na.rm = TRUE)), lwd = 2, col = "red")
}


# run stability plotting functions
setwd("A:/1_UNSW/floristic/RCP2/results/")

# stability.regimix run - in n chunks and times times for id model - on cluster
combine.registab = function(n, times, id) {
  registab.data = lapply(as.list(1:n), FUN=function(x) {get(load(paste0("stability/stability.",id,".",x,".RData")))})
  disty.data = lapply(registab.data, FUN=function(x) {x$disty})
  disty = Reduce("+", disty.data) / n
  predlogls = array(NA,c(length(registab.data[[1]]$oosSizeRange),
                         registab.data[[1]]$model$n,
                         registab.data[[1]]$time * length(registab.data)))
  for (i in 1:n) {predlogls[,,(i*times-(times-1)):(i*times)]=registab.data[[i]]$predlogls}
  registab.out = registab.data[[1]]
  registab.out$disty = disty
  registab.out$predlogls = predlogls
  return(registab.out)
}


stability.nosp = combine.registab(7, 100, "nosp")
CairoPDF(file="stability/stability.nosp.pdf", width=12, height=7)
stability.regimix.plot(stability.nosp)
dev.off()

stability.sp = combine.registab(7, 100, "sp")
CairoPDF(file="stability/stability.sp.pdf", width=12, height=7)
stability.regimix.plot(stability.sp)
dev.off()

