#library(Matrix, lib.loc="../../R/x86_64-unknown-linux-gnu-library/3.1/")
library(RCPmod)

job = Sys.getenv("PBS_ARRAYID")

if ("species" %in% commandArgs()) {species.model = "species"}
if ("nospecies" %in% commandArgs()) {species.model = "nospecies"}
if ("nosp-at-sp" %in% commandArgs()) {species.model = "nosp-at-sp"}

if (species.model=="species") {
  load("fit.regi.sp.RData")
  fit.regiboot.sp = regiboot(fit.regi.sp, nboot=5)
  save(fit.regiboot.sp, file=paste0("fit.regiboot.sp.",job,".RData"))
}

if (species.model=="nospecies") {
  load("fit.regi.nosp.RData")
  fit.regiboot.nosp = regiboot(fit.regi.nosp, nboot=5)
  save(fit.regiboot.nosp, file=paste0("fit.regiboot.nosp.",job,".RData"))
}

if (species.model=="nosp-at-sp") {
  load("fit.regi.nosp-at-sp.RData")
  fit.regiboot.nospsp = regiboot(fit.regi.nospsp, nboot=5)
  save(fit.regiboot.nospsp, file=paste0("fit.regiboot.nosp-at-sp.",job,".RData"))
}