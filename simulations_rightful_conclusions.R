## A function to run a single batch of simulations with a 3-parameter set
## and to return proportion of rigthful conclusions, proportion of S-type
## errors and mean M-type errors  
source("functions.R")
library(parallel)

## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4

################################################################################
## t-test simulations
################################################################################
## Hypercube sampling (see introductory vignette of pse package)
## Variables to be sampled: standard effect (as t-values), standandard error, sample size
factors <- c("D", "SE", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=3), list(min=0.1, max=1), list(min=5, max=100))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.ttest")))
}
## Hypercube sample of nsamp combinations of parameters
## Start cluster (set the file mpd.hosts with the desired number of
## cores and processes, see pse:machinefile)
## cl <- makePSOCKcluster(machinefile("../mpd.hosts"))
cl <- makePSOCKcluster(machinefile("mpd.hosts"))
## Export needed functions
##clusterEvalQ(cl, source("functions.R"))
clusterExport(cl, c("nrepl", "AICc", "wp.ttest", "sim.averages"))
## Run the simulations 
myLHS <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
t.results <- cbind(get.data(myLHS), get.results(myLHS))
names(t.results)[4:12] <- c("p.NHT.right", "p.AIC.right", "p.mismatch","mean.NHT.M", "mean.AIC.M",
                            "p.NHT.S",  "p.AIC.S", "mean.pvalue", "mean.wH0")
save.image()
## Same simulation with IT criterium deltaAIC >2
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.ttest", delta2=TRUE)))
}
myLHS.d2 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
t.results.d2 <- cbind(get.data(myLHS.d2), get.results(myLHS.d2))
names(t.results.d2)[4:12] <- c("p.NHT.right", "p.AIC.right", "p.mismatch","mean.NHT.M", "mean.AIC.M",
                            "p.NHT.S",  "p.AIC.S", "mean.pvalue", "mean.wH0")
save.image()
## Stop cluster
stopCluster(cl)
## Save results in separated binary
save(t.results, t.results.d2, file="tresults.RData")


################################################################################
## correlation simulations
################################################################################
## Hypercube sampling (see introductory vignette of pse package)
## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("D", "SE", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=15), list(min=0.1, max=1), list(min=5, max=100))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.cortest")))
}
## Hypercube sample of nsamp combinations of parameters
## Start cluster (set the file mpd.hosts with the desired number of
## cores and processes, see pse:machinefile)
## cl <- makePSOCKcluster(machinefile("../mpd.hosts"))
cl <- makePSOCKcluster(machinefile("mpd.hosts"))
## Export needed functions
##clusterEvalQ(cl, source("functions.R"))
clusterExport(cl, c("nrepl", "AICc", "wp.cortest", "sim.averages", "rmvnorm", "dmvnorm"))
## Run the simulations 
myLHS2 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
cor.results <- cbind(get.data(myLHS2),  get.results(myLHS2))
names(cor.results)[4:12] <- c("p.NHT.right", "p.AIC.right", "p.mismatch","mean.NHT.M", "mean.AIC.M",
                              "p.NHT.S",  "p.AIC.S", "mean.pvalue", "mean.wH0")
## Adds Pearson correlation coeff values to the results (which is a function of tval and N)
cor.results$rpears <- with(get.data(myLHS2), D/sqrt(N-2+D^2))
save.image()
## Run the simulations with citerium deltaAIC>2
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.cortest", delta2=TRUE)))
}
myLHS2.d2 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
cor.results.d2 <- cbind(get.data(myLHS2.d2),  get.results(myLHS2.d2))
names(cor.results,d2)[4:12] <- c("p.NHT.right", "p.AIC.right", "p.mismatch","mean.NHT.M", "mean.AIC.M",
                              "p.NHT.S",  "p.AIC.S", "mean.pvalue", "mean.wH0")
## Adds Pearson correlation coeff values to the results (which is a function of tval and N)
cor.results.d2$rpears <- with(get.data(myLHS2.d2), D/sqrt(N-2+D^2))
save.image()
stopCluster(cl)
## Save results in a separated binary
save(cor.results, cor.results.d2, file="corresults.RData")
