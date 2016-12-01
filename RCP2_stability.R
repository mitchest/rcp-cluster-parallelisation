# stability plots ---------------------------------------------------------

library(RCPmod)
library(data.table)

# load, run, dump - stability data function
load("fit.regi.nosp.RData")
stability.nosp = stability.regimix(model=fit.regi.nosp, oosSizeRange=c(750), times=10, mc.cores=4, doPlot=F)
save(stability.nosp, file="stability.nosp.RData")
rm(fit.regi.nosp,stability.nosp)
gc()

load("fit.regi.sp.RData")
stability.sp = stability.regimix(model=fit.regi.sp, oosSizeRange=c(750), times=10, mc.cores=4, doPlot=F)
save(stability.sp, file="stability.sp.RData")
rm(fit.regi.sp,stability.sp)
gc()


# stitch separate runs
stablist_nosp <- lapply(X = c("1-20-50",100,300,750), FUN = function(x) {get(load(paste0("results/stability/stability.nosp.",x,".RData")))})
stablist_sp <- lapply(X = c("1-20-50",100,300,750), FUN = function(x) {get(load(paste0("results/stability/stability.sp.",x,".RData")))})

registab_sp <- stablist_sp[[1]]
registab_sp$oosSizeRange <- c(1,20,50,100,300,750)
registab_sp$disty <- as.matrix(rbindlist(lapply(X = stablist_sp, FUN = function(x){as.data.frame(x$disty)})))
registab_sp$predlogls <- array(NA, dim = c(6,4715,10))
registab_sp$predlogls[1:3,,] <- stablist_sp[[1]]$predlogls
registab_sp$predlogls[4,,] <- stablist_sp[[2]]$predlogls
registab_sp$predlogls[5,,] <- stablist_sp[[3]]$predlogls
registab_sp$predlogls[6,,] <- stablist_sp[[4]]$predlogls

registab_nosp <- stablist_nosp[[1]]
registab_nosp$oosSizeRange <- c(1,20,50,100,300,750)
registab_nosp$disty <- as.matrix(rbindlist(lapply(X = stablist_nosp, FUN = function(x){as.data.frame(x$disty)})))
registab_nosp$predlogls <- array(NA, dim = c(6,4715,10))
registab_nosp$predlogls[1:3,,] <- stablist_nosp[[1]]$predlogls
registab_nosp$predlogls[4,,] <- stablist_nosp[[2]]$predlogls
registab_nosp$predlogls[5,,] <- stablist_nosp[[3]]$predlogls
registab_nosp$predlogls[6,,] <- stablist_nosp[[4]]$predlogls


# test
plot.registab(registab_sp)
plot.registab(registab_nosp)

save(registab_sp, file="registab_sp.RData")
save(registab_nosp, file="registab_nosp.RData")



# modified functions so plooting can be done later
# cooks.distance.regimix.katana <- function (model, ..., oosSize = 1, times = model$n, mc.cores = 1, 
#           quiet = FALSE) {
#   if (oosSize > model$n%/%2) 
#     stop("Out of sample is more than half the size of the data! This is almost certainly an error.  Please set `oosSize' to something smaller.")
#   if (is.null(model$titbits)) 
#     stop("Model doesn't contain all information required for cross validation.  Please supply model with titbits (from titbits=TRUE in regimix call)")
#   if (!quiet) 
#     pb <- txtProgressBar(min = 1, max = times, style = 3, 
#                          char = "><(('> ")
#   funny <- function(x) {
#     if (!quiet) 
#       setTxtProgressBar(pb, x)
#     if (oosSize != 1 | times != model$n) 
#       OOBag <- sample(1:model$n, oosSize, replace = FALSE)
#     else OOBag <- x
#     inBag <- (1:model$n)[!(1:model$n) %in% OOBag]
#     new.wts <- model$titbits$wts
#     new.wts[OOBag] <- 0
#     control <- model$titbits$control
#     control$quiet <- TRUE
#     control$trace <- 0
#     tmpmodel <- regimix.fit(outcomes = model$titbits$Y, W = model$titbits$W, 
#                             X = model$titbits$X, offy = model$titbits$offset, 
#                             wts = new.wts, disty = model$titbits$disty, nRCP = model$nRCP, 
#                             power = model$titbits$power, inits = unlist(model$coef), 
#                             control = control, n = model$n, S = model$S, p.x = model$p.x, 
#                             p.w = model$p.w)
#     OOSppPreds <- matrix(NA, nrow = tmpmodel$n, ncol = tmpmodel$S)
#     for (ss in 1:tmpmodel$S) OOSppPreds[OOBag, ss] <- rowSums(tmpmodel$mus[OOBag, 
#                                                                            ss, ] * tmpmodel$pis[OOBag, , drop = FALSE])
#     newPis <- tmpmodel$pis
#     r.negi <- model$pis - newPis
#     r.negi[OOBag, ] <- NA
#     r.negi <- colMeans(r.negi, na.rm = TRUE)
#     alpha.score <- as.numeric(rep(NA, model$S))
#     tau.score <- as.numeric(matrix(NA, ncol = model$S, nrow = model$nRCP - 
#                                      1))
#     beta.score <- as.numeric(matrix(NA, ncol = ncol(model$titbits$X), 
#                                     nrow = model$nRCP - 1))
#     if (model$p.w > 0) {
#       gamma.score <- as.numeric(matrix(NA, nrow = model$S, 
#                                        ncol = model$p.w))
#       gamma <- tmpmodel$coef$gamma
#       W <- model$titbits$W
#     }
#     else gamma.score <- W <- gamma <- -999999
#     if (model$titbits$disty %in% 3:5) {
#       disp.score <- as.numeric(rep(NA, model$S))
#       disp <- coef(model)$logDisp
#     }
#     else disp.score <- -999999
#     scoreContri <- -999999
#     logCondDens <- as.numeric(matrix(NA, nrow = model$n, 
#                                      ncol = model$nRCP))
#     logls <- as.numeric(rep(NA, model$n))
#     conv <- as.integer(0)
#     tmplogl <- .Call("RCP_C", as.numeric(model$titbits$Y), 
#                      as.numeric(model$titbits$X), as.numeric(model$titbits$W), 
#                      as.numeric(model$titbits$offset), as.numeric(model$titbits$wts), 
#                      as.integer(model$S), as.integer(model$nRCP), as.integer(model$p.x), 
#                      as.integer(model$p.w), as.integer(model$n), as.integer(model$titbits$disty), 
#                      as.numeric(tmpmodel$coef$alpha), as.numeric(tmpmodel$coef$tau), 
#                      as.numeric(tmpmodel$coef$beta), as.numeric(gamma), 
#                      as.numeric(tmpmodel$coef$disp), as.numeric(model$titbits$power), 
#                      as.numeric(model$titbits$control$penalty), as.numeric(model$titbits$control$penalty.tau), 
#                      as.numeric(model$titbits$control$penalty.gamma), 
#                      as.numeric(model$titbits$control$penalty.disp[1]), 
#                      as.numeric(model$titbits$control$penalty.disp[2]), 
#                      alpha.score, tau.score, beta.score, gamma.score, 
#                      disp.score, scoreContri, as.numeric(tmpmodel$pis), 
#                      as.numeric(tmpmodel$mus), logCondDens, logls, as.integer(model$titbits$control$maxit), 
#                      as.integer(model$titbits$control$trace), as.integer(model$titbits$control$nreport), 
#                      as.numeric(model$titbits$control$abstol), as.numeric(model$titbits$control$reltol), 
#                      as.integer(conv), as.integer(FALSE), as.integer(TRUE), 
#                      as.integer(FALSE), as.integer(FALSE), as.integer(FALSE), 
#                      PACKAGE = "RCPmod")
#     ret.logl <- rep(NA, model$n)
#     ret.logl[OOBag] <- logls[OOBag]
#     return(list(OOSppPreds = OOSppPreds, cooksDist = r.negi, 
#                 predLogL = ret.logl))
#   }
#   if (!quiet & mc.cores > 1 & Sys.info()["sysname"] != "Windows") 
#     message("Progress bar may not be monotonic due to the vaguaries of parallelisation")
#   tmp <- parallel::mclapply(1:times, funny, mc.cores = mc.cores, mc.preschedule = FALSE)
#   if (!quiet) 
#     message("")
#   cooksD <- t(sapply(tmp, function(x) x$cooksDist))
#   OOpreds <- array(NA, dim = c(model$n, model$S, times), dimnames = list(rownames(model$titbits$X), 
#                                                                          colnames(model$titbits$Y), paste("CVset", 1:times, sep = "")))
#   for (bb in 1:times) OOpreds[, , bb] <- tmp[[bb]]$OOSppPreds
#   logls <- sapply(tmp, function(x) x$predLogL)
#   colnames(logls) <- rownames(cooksD) <- paste("OOS", 1:times, 
#                                                sep = "_")
#   ret <- list(Y = model$titbits$Y, CV = OOpreds, cooksD = cooksD, 
#               predLogL = logls)
#   class(ret) <- "regiCooksD"
#   return(ret)
# }


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

# stability.regimix.plot <- function (stability.regimix.data) {
#   # declare all the plotting variables
#   oosSizeRange = stability.regimix.data$oosSizeRange
#   time = stability.regimix.data$times
#   disty = stability.regimix.data$disty
#   predlogls = stability.regimix.data$predlogls
#   model = stability.regimix.data$model
#   # plotting
#   par(mfrow = c(1, 2))
#   matplot(c(0, oosSizeRange), rbind(0, disty), type = "b", 
#           ylab = "Distance from Full Model Predictions", xlab = "Number of Obs Removed", 
#           main = "Stability of Group Predictions", col = 1:model$nRCP, 
#           pch = as.character(1:model$nRCP), lty = 1)
#   legend("topleft", bty = "n", lty = 1, pch = as.character(1:model$nRCP), 
#          col = 1:model$nRCP, legend = paste("RCP ", 1:model$nRCP, 
#                                             sep = ""))
#   plot(rep(oosSizeRange, each = prod(dim(predlogls[1, , ]))), 
#        predlogls, pch = 20, ylab = "Pred LogL (OOS)", xlab = "Number of Obs Removed", 
#        main = "Stability of Pred Logl", xlim = c(0, max(oosSizeRange)), 
#        type = "n")
#   points(rep(0, model$n), model$logl.sites, col = "blue", pch = 20)
#   for (ii in oosSizeRange) {
#     points(rep(ii, prod(dim(predlogls[1, , ]))), predlogls[oosSizeRange == 
#                                                              ii, , ], pch = 20)
#   }
#   lines(c(0, oosSizeRange), c(mean(model$logl.sites), apply(predlogls, 
#                                                             1, mean, na.rm = TRUE)), lwd = 2, col = "red")
# }