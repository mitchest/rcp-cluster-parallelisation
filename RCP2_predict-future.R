library(RCPmod)
library(Cairo)
library(foreign)
library(data.table)
library(dplyr)
library(raster)
library(rgdal)

setwd("A:/1_UNSW/floristic/RCP2")



# load projections --------------------------------------------------------

# # load covariates
# future_covars <- read.dbf("covars/NSW_2km_future.dbf", as.is=T)
# 
# # scale precipseas back to percentage
# precipse_cols <- grep(x = names(future_covars), pattern = "*precipse*")[-1]
# future_covars[,precipse_cols] <- future_covars[,precipse_cols] * 100
# 
# # # check which vars needed
# # load("predict/fit.regi.nosp.RData")
# # attr(fit.regi.nosp$titbits$X, "dimnames")[[2]]
# 
# # make covar data and location data
# locs <- future_covars[,c("lat","long")]
# 
# # load data from main analysis script
# model.covars = covariates.species[,model.covariates.vector]
# 
# covar_names <- c("tempmtcp.1","tempmtcp.2",
#                  "precipseas.1","precipseas.2",
#                  "precipann.1","precipann.2",
#                  "tempmtwp.1","tempmtwp.2",
#                  "rough500.1","rough500.2",
#                  "bd200.1","bd200.2",
#                  "ph200.1","ph200.2")
# 
# # calculate polynomials
# covars15 = data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=future_covars$tempmtcp15),
#                     predict(poly(model.covars$precipseas, 2), newdata=future_covars$precipse15),
#                     predict(poly(model.covars$precipann, 2), newdata=future_covars$precipan15),
#                     predict(poly(model.covars$tempmtwp, 2), newdata=future_covars$tempmtwp15),
#                     predict(poly(model.covars$rough500, 2), newdata=future_covars$rough500),
#                     predict(poly(model.covars$bd200, 2), newdata=future_covars$bd200),
#                     predict(poly(model.covars$ph200, 2), newdata=future_covars$ph200))
# names(covars15) = covar_names
# 
# covars35 = data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=future_covars$tempmtcp35),
#                       predict(poly(model.covars$precipseas, 2), newdata=future_covars$precipse35),
#                       predict(poly(model.covars$precipann, 2), newdata=future_covars$precipan35),
#                       predict(poly(model.covars$tempmtwp, 2), newdata=future_covars$tempmtwp35),
#                       predict(poly(model.covars$rough500, 2), newdata=future_covars$rough500),
#                       predict(poly(model.covars$bd200, 2), newdata=future_covars$bd200),
#                       predict(poly(model.covars$ph200, 2), newdata=future_covars$ph200))
# names(covars35) = covar_names
# 
# covars55 = data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=future_covars$tempmtcp55),
#                       predict(poly(model.covars$precipseas, 2), newdata=future_covars$precipse55),
#                       predict(poly(model.covars$precipann, 2), newdata=future_covars$precipan55),
#                       predict(poly(model.covars$tempmtwp, 2), newdata=future_covars$tempmtwp55),
#                       predict(poly(model.covars$rough500, 2), newdata=future_covars$rough500),
#                       predict(poly(model.covars$bd200, 2), newdata=future_covars$bd200),
#                       predict(poly(model.covars$ph200, 2), newdata=future_covars$ph200))
# names(covars55) = covar_names
# 
# covars75 = data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=future_covars$tempmtcp75),
#                       predict(poly(model.covars$precipseas, 2), newdata=future_covars$precipse75),
#                       predict(poly(model.covars$precipann, 2), newdata=future_covars$precipan75),
#                       predict(poly(model.covars$tempmtwp, 2), newdata=future_covars$tempmtwp75),
#                       predict(poly(model.covars$rough500, 2), newdata=future_covars$rough500),
#                       predict(poly(model.covars$bd200, 2), newdata=future_covars$bd200),
#                       predict(poly(model.covars$ph200, 2), newdata=future_covars$ph200))
# names(covars75) = covar_names
# 
# future_covars_list <- list(covars15=covars15, covars35=covars35, covars55=covars55, covars75=covars75)
# save(future_covars_list, locs, file="future_projections/future_covars_list.RData")
load("future_projections/future_covars_list.RData")



# make future predictions -------------------------------------------------

# load the model and bootstrap sample distribuion for CI estimates
load("predict/fit.regiboot.nosp.RData")
load("predict/fit.regi.nosp.x.RData")

predicted.nosp.estimate <- list()
predicted.nosp.btpreds <- list()
predicted.nosp.low <- list()
predicted.nosp.up <- list()

for (ii in 2:length(future_covars_list)) {
  
  covars <- future_covars_list[[ii]]
  covar_year <- names(future_covars_list)[ii]
  
  predicted.nosp <- predict.regimix(object=fit.regi.nosp, newdata=covars, nboot=0)
  predicted.nosp.estimate[[covar_year]] <- data.frame(locs, predicted.nosp)
  print(round(colSums(predicted.nosp)))
  
  # splitting this just for memory issues
  predicted.nosp.list <- list()
  idx <- round(nrow(future_covars_list[[1]])/3)
  for (i in 1:3) {
    if (i==1) {
      predicted.nosp.list[[i]] <- predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp,
                                                 newdata=covars[1:idx,])
      gc()
    }
    if (i==2) {
      predicted.nosp.list[[i]] <- predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp,
                                                 newdata=covars[(idx+1):(idx*2),])
      gc()
    }
    if (i==3) {
      predicted.nosp.list[[i]] <- predict.regimix(object=fit.regi.nosp, object2=fit.regiboot.nosp,
                                                  newdata=covars[(idx*2+1):nrow(covars),])
      gc()
    }
  }
  # adding location data back for GIS plotting
  predicted.nosp.btpreds[[covar_year]] <- data.frame(locs, rbind(predicted.nosp.list[[1]]$bootPreds[,],
                                                                 predicted.nosp.list[[2]]$bootPreds[,],
                                                                 predicted.nosp.list[[3]]$bootPreds[,]))
  predicted.nosp.low[[covar_year]] <- data.frame(locs, rbind(predicted.nosp.list[[1]]$bootCIs[,,1],
                                                             predicted.nosp.list[[2]]$bootCIs[,,1],
                                                             predicted.nosp.list[[3]]$bootCIs[,,1]))
  predicted.nosp.up[[covar_year]] <- data.frame(locs, rbind(predicted.nosp.list[[1]]$bootCIs[,,2],
                                                            predicted.nosp.list[[2]]$bootCIs[,,2],
                                                            predicted.nosp.list[[3]]$bootCIs[,,2]))
  
  rm(covars, predicted.nosp, predicted.nosp.list)
  gc()
}


# save the results
save(predicted.nosp.estimate, file="future_projections/nosp_future_estimates.RData")
save(predicted.nosp.btpreds, file="future_projections/nosp_future_btpreds.RData")
save(predicted.nosp.low, file="future_projections/nosp_future_lowerCI.RData")
save(predicted.nosp.up, file="future_projections/nosp_future_upperCI.RData")

# # write them out for import into a GIS
# for (i in 1:length(predicted.nosp.estimate)) {
#   covar_year <- names(predicted.nosp.estimate)[i]
#   
#   estimates <- predicted.nosp.estimate[[i]]
#   write.csv(estimates, paste0("future_projections/predicted_nosp_estimate_",covar_year,".csv"), row.names = F)
#   
#   btpreds <- predicted.nosp.btpreds[[i]]
#   write.csv(btpreds, paste0("future_projections/predicted_nosp_btpreds_",covar_year,".csv"), row.names = F)
# 
#   lowers <- predicted.nosp.low[[i]]
#   write.csv(lowers, paste0("future_projections/predicted_nosp_low_",covar_year,".csv"), row.names = F)
#   
#   uppers <- predicted.nosp.up[[i]]
#   write.csv(uppers, paste0("future_projections/predicted_nosp_up_",covar_year,".csv"), row.names = F)
#   
# }



# export shapefile --------------------------------------------------------

for (i in 1:length(predicted.nosp.estimate)) {
  covar_year <- names(predicted.nosp.estimate)[i]
  
  estimate <- predicted.nosp.estimate[[i]]
  coordinates(estimate) <- ~long+lat
  proj4string(estimate) <- CRS("+init=epsg:4326")
  writeOGR(obj = estimate, driver = "ESRI Shapefile",
           dsn = "future_projections", layer = paste0("estimate_",covar_year))
  
  btpreds <- predicted.nosp.btpreds[[i]]
  coordinates(btpreds) <- ~long+lat
  proj4string(btpreds) <- CRS("+init=epsg:4326")
  writeOGR(obj = btpreds, driver = "ESRI Shapefile",
           dsn = "future_projections", layer = paste0("btpreds_",covar_year))
  
  low <- predicted.nosp.low[[i]]
  coordinates(low) <- ~long+lat
  proj4string(low) <- CRS("+init=epsg:4326")
  writeOGR(obj = low, driver = "ESRI Shapefile",
           dsn = "future_projections", layer = paste0("low_",covar_year))
  
  up <- predicted.nosp.up[[i]]
  coordinates(up) <- ~long+lat
  proj4string(up) <- CRS("+init=epsg:4326")
  writeOGR(obj = up, driver = "ESRI Shapefile",
           dsn = "future_projections", layer = paste0("up_",covar_year))
  
}











