library(RCPmod)
library(Cairo)

# load fitted models
setwd("A:/1_UNSW/floristic/RCP2/")
load("predict/fit.regi.nosp.RData")

# load covariate data and funs
# not this source includes data prep/load, so needs to be updated if modelling data/analsysi changes
source("RCP2_pareff_fundat.R")


# partial effects data prep -----------------------------------------------

# make a new data matrix
new.dat <- matrix(NA, ncol=ncol(covar.data), nrow=1000)

# get the col means
X.means <- colMeans(covar.data)
# get the col medians
X.meds = colMedian(covar.data)
# get col mins
X.mins <- colMin(covar.data)

# add meds or means to the new data matrix
for (i in 1:length(X.meds)) {
  new.dat[,i] <- X.meds[i]
}

# give it the names used for fitting
new.dat <- data.frame(new.dat)
names(new.dat) <- names(covar.data)

## once you have categorical variables
# # hack the catgorical X's into shape - need a better solution for this
# # choose the group to be the prediciton level
# sort(table(covar.data$XXXXXX))
# new.dat$XXXXXX <- as.factor(covar.data$XXXXXX)[1]
# # ## or test the XXXXXX varaible
# # new.dat = new.dat[1:length(unique(covar.data$XXXXXX)),]
# # new.dat$XXXXXX = as.factor(unique(covar.data$XXXXXX))
# # new.dat.XXXXXX = new.dat

# make the polynomials - ensure names match original call
new.dat <- data.frame(predict(poly(covar.data$tempmtcp, 2), newdata=new.dat$tempmtcp),
                    predict(poly(covar.data$precipseas, 2), newdata=new.dat$precipseas),
                    predict(poly(covar.data$precipann, 2), newdata=new.dat$precipann),
                    predict(poly(covar.data$tempmtwp, 2), newdata=new.dat$tempmtwp),
                    predict(poly(covar.data$rough500, 2), newdata=new.dat$rough500),
                    predict(poly(covar.data$bd200, 2), newdata=new.dat$bd200),
                    predict(poly(covar.data$ph200, 2), newdata=new.dat$ph200))
names(new.dat) <- c("tempmtcp.1","tempmtcp.2",
                  "precipseas.1","precipseas.2",
                  "precipann.1","precipann.2",
                  "tempmtwp.1","tempmtwp.2",
                  "rough500.1","rough500.2",
                  "bd200.1","bd200.2",
                  "ph200.1","ph200.2")



# partial effects plots ---------------------------------------------------

CairoPDF(file = "RCP2_11RCP_pareff.pdf", width = 10, height = 10)
par(mfcol=c(1,1))
for (i in names(covar.data)) {
  ParPlotMany.regimix(fit.regi.nosp, covar.data, new.dat, i, 1:11)
}
dev.off()




# plot actual predicted values --------------------------------------------

CairoPDF(file = "RCP2_11RCP_predvals.pdf", width = 10, height = 14)
par(mfcol=c(11,7))
for (i in 1:11) {
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$tempmtcp, ylab = i, xlab = names(covar.data)[1])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$precipseas, ylab = i, xlab = names(covar.data)[2])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$precipann, ylab = i, xlab = names(covar.data)[3])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$tempmtwp, ylab = i, xlab = names(covar.data)[4])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$rough500, ylab = i, xlab = names(covar.data)[5])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$bd200, ylab = i, xlab = names(covar.data)[6])
  plot(fit.regi.nosp$postProbs[,i] ~ covar.data$ph200, ylab = i, xlab = names(covar.data)[7])
}
dev.off()







