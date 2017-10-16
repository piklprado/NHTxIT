source("functions.R")
library(parallel)

## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4

################################################################################
## t-test simulations
################################################################################
## With H0 FALSE
## Hypercube sampling (see introductory vignette of pse package)
## Variables to be sampled: standard effect (as t-values), standandard deviation, sample size
factors <- c("D", "SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=0.1, max=8), list(min=5, max=150))
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
names(t.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
## Wiht H0 true ##
## Hypercube sampling (see introductory vignette of pse package)
## Variables to be sampled: standard effect (as t-values), standandard deviation, sample size
factors <- c("SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=5, max=150))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages.null, my.data[,1], my.data[,2],
                  MoreArgs = list(nrep = nrepl, function.name="wp.ttest")))
}
## Export needed functions
clusterExport(cl, c("sim.averages.null"))
## Run the simulations 
myLHS0 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
t.results <- cbind(get.data(myLHS), get.results(myLHS))
names(t.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                            "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
t.results0 <- cbind(get.data(myLHS0), get.results(myLHS0))
names(t.results0)[3:15] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
save(t.results, t.results0, file="tresults.RData")
stopCluster(cl)

save.image()

################################################################################
## correlation simulations
################################################################################
## For Null hypothesis false ##
## Hypercube sampling (see introductory vignette of pse package)
## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("D", "SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=0.1, max=8), list(min=10, max=300))
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
names(cor.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
## Adds Pearson correlation coeff values to the results (which is a function of tval and N)
cor.results$rpears <- with(get.data(myLHS2), D/sqrt(N-2+D^2))

## With H0 True #
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=10, max=300))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages.null, my.data[,1], my.data[,2],
                  MoreArgs = list(nrep = nrepl, function.name="wp.cortest")))
}

clusterExport(cl, c("sim.averages.null"))
## Run the simulations 
myLHS20 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
cor.results0 <- cbind(get.data(myLHS20),  get.results(myLHS20))
names(cor.results0)[3:15] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
##
save(cor.results, cor.results0, file="corresults.RData")
stopCluster(cl)
save.image()

################################################################################
## ANOVA simulations
################################################################################
## Hypercube sampling (see introductory vignette of pse package)
## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("D", "SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=0.1, max=8), list(min=10, max=100))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.anova")))
}
## Hypercube sample of nsamp combinations of parameters
## Start cluster (set the file mpd.hosts with the desired number of
## cores and processes, see pse:machinefile)
## cl <- makePSOCKcluster(machinefile("../mpd.hosts"))
cl <- makePSOCKcluster(machinefile("mpd.hosts"))
## Export needed functions
##clusterEvalQ(cl, source("functions.R"))
clusterExport(cl, c("nrepl", "AICc", "wp.anova", "sim.averages"))
## Run the simulations 
myLHS3 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
anova.results <- cbind(get.data(myLHS3),  get.results(myLHS3))
names(anova.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                                "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")

## H0 is true ##
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=10, max=100))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages.null, my.data[,1], my.data[,2],
                  MoreArgs = list(nrep = nrepl, function.name="wp.anova")))
}
##
clusterExport(cl, c("sim.averages.null"))
## Run the simulations 
myLHS30 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
anova.results0 <- cbind(get.data(myLHS30),  get.results(myLHS30))
names(anova.results0)[3:15] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                                "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
save.image()
save(anova.results, anova.results0, file="anovaresults.RData")
stopCluster(cl)

################################################################################
## Linear regression simulations
################################################################################
## 1. Without colinearity
## Hypercube sampling (see introductory vignette of pse package)
## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("D", "SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=0.1, max=8), list(min=10, max=300))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.lm")))
}
## Hypercube sample of nsamp combinations of parameters
## Start cluster (set the file mpd.hosts with the desired number of
## cores and processes, see pse:machinefile)
## cl <- makePSOCKcluster(machinefile("../mpd.hosts"))
cl <- makePSOCKcluster(machinefile("mpd.hosts"))
## Export needed functions
##clusterEvalQ(cl, source("functions.R"))
clusterExport(cl, c("nrepl", "AICc", "wp.lm", "sim.averages", "rmvnorm"))
## Run the simulations 
myLHS4 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
lm.results <- cbind(get.data(myLHS4),  get.results(myLHS4))
names(lm.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
## 2. With colinearity
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages, my.data[,1], my.data[,2], my.data[,3],
                  MoreArgs = list(nrep = nrepl, function.name="wp.lm", rpears=0.25)))
}
## Hypercube sample of nsamp combinations of parameters
## Run the simulations 
myLHS5 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
lm.colin.results <- cbind(get.data(myLHS5),  get.results(myLHS5))
names(lm.colin.results)[4:16] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                                   "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")

## H0 is true ##
## 1. Without colinearity
## Hypercube sampling (see introductory vignette of pse package)
## Number of hypercube samples
nsampH <- 2000
## Number of replicate simulations for each parameter sample in the hypercube
nrepl <- 1e4
## Variables to be sampled: standard effect (as tvalues), standandard error of marginal distributions, sample size
factors <- c("SD", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=0.1, max=8), list(min=10, max=300))
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages.null, my.data[,1], my.data[,2],
                  MoreArgs = list(nrep = nrepl, function.name="wp.lm")))
}
##
clusterExport(cl, c("sim.averages.null"))
## Run the simulations 
myLHS40 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
lm.results0 <- cbind(get.data(myLHS40),  get.results(myLHS40))
names(lm.results0)[3:15] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                              "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
## 2. With colinearity
## Accessory function for a single run with a combination of parameters
modelRun <- function (my.data) {
    return(mapply(sim.averages.null, my.data[,1], my.data[,2],
                  MoreArgs = list(nrep = nrepl, function.name="wp.lm", rpears=0.5)))
}
## Hypercube sample of nsamp combinations of parameters
## Run the simulations 
myLHS50 <- LHS(modelRun, factors, nsampH, q, q.arg, cl=cl)
## dataframe with parameter values and results
lm.colin.results0 <- cbind(get.data(myLHS50),  get.results(myLHS50))
names(lm.colin.results0)[3:15] <- c("p.NHT.right", "p.AIC.right", "p.AIC.right.2", "p.mismatch", "p.mismatch.2", "mean.NHT.M",
                                   "mean.AIC.M", "mean.AIC.M.2", "p.NHT.S",  "p.AIC.S", "p.AIC.S.2", "mean.pvalue", "mean.wH0")
##
save(lm.results, lm.colin.results, lm.results0, lm.colin.results0, file="lmresults.RData")
save.image()
stopCluster(cl)

