# rcp-cluster-parallelisation  

# README #

This repo provides code related to the following two papers:
  Foster, Hill and Lyons (201x) "Ecological Grouping of Survey Sites when Sampling Artefacts are Present". Journal Name: XX(X), XXX–XXX.  
  DOI: 10.xxxxx/xxx

  Lyons, Foster and Keith (201x) "Simultaneous vegetation classification and mapping using a statistical model". Journal Name: XX(X), XXX–XXX.  
  DOI: 10.xxxxx/xxx

The purpose of this repo is not necessarily to reproduce results form those papers, but to provide guidence for how analysis with the R package {RCPmod} can be parallelised on a (relatively generic) computing cluster. The data itself can be found at this repo:  
https://github.com/mitchest/rcp-survey-artifacts  

### Running code ###

* The analysis is contained within the .R scripts (obviously)  
* And the architechture specific paralellisation is found in the .pbs file (written for a Tourque cluster)  
* I am currently workign on getting the hard-coded directories out of there, for the meantime, sorry

### Contribution guidelines ###

* The best route for feed back on {RCPmod} is to lookup the function authors:  
* https://cran.r-project.org/web/packages/RCPmod/index.html  
* Otherwise, feel free to suggest edits here  

### Contact ###

* Mitchell Lyons (mitchell.lyons@gmail.com / mitchell.lyons@unsw.edu.au)
