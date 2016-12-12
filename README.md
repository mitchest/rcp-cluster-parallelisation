# rcp-cluster-parallelisation  

# README #

This repo provides code related to the following two papers:  

  + Foster, Hill and Lyons (accepted) "Ecological Grouping of Survey Sites when Sampling Artefacts are Present". Journal of the Royal Statistical Society: Series C (Applied Statistics). DOI: xx.xxxxx/xxx

  + Lyons, Foster and Keith (in review) "Data driven, simultaneous vegetation classification and mapping at large spatial scales". Journal Name. DOI: xx.xxxxx/xxx  

The purpose of this repo is not necessarily to reproduce results form those papers, but to provide guidence for how analysis with the R package {RCPmod} can be parallelised on a (relatively generic) computing cluster. You would want (have) to do this when your analysis reaches 1000's of observations and 100's of species. The data is provided here as an .RData file for convenience. For some non-cluster work flows for fitting RCP models, see the repo https://github.com/mitchest/rcp-survey-artifacts  

### Running code ###

* The analysis is contained within the .R scripts (obviously)  
* The core fitting code is in "RCP2_katana.R", the other scripts are labelled more sensibly in terms of what they do =)
* The architechture specific paralellisation is found in the .pbs file (written for a Tourque cluster)  
* I am currently working on getting the hard-coded directories out of there, for the meantime, sorry

### Contribution guidelines ###

* The best route for feed back on {RCPmod} is to lookup the function authors:  
* https://cran.r-project.org/web/packages/RCPmod/index.html  
* Otherwise, feel free to suggest edits here  

### Contact ###

* Mitchell Lyons (mitchell.lyons@gmail.com / mitchell.lyons@unsw.edu.au)
