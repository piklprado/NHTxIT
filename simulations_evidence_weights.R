################################################################################
## Relationship between p-values an evidence weights in a t-test
################################################################################

## Hypercube sampling (see introductory vignette of pse package)
## Variables to be sampled: arguments
factors <- c("t.val", "sd1", "N")
## Discrete uniform random deviates function, to sample discrete
## sample sizes
qdunif<-function(p, min, max) floor(qunif(p, min, max))
## Distributions to sample the variables; all uniform in this case
q <- c("qunif", "qunif", "qdunif")
## Arguments of each distribution (min and max values of the uniform distributions)
q.arg <- list(list(min=1.28, max=3), list(min=1, max=5), list(min=10, max=50))
## To run the sampling simulation many times
modelRun <- function (my.data) {
 return(mapply(wp.ttest, my.data[,1], my.data[,2], my.data[,3]))
}

## Hypercube sample of 1000 combinations of parameters
myLHS <- LHS(modelRun, factors, 2000, q, q.arg)
## Relationship between p-value and weight of evidence for each
## parameter combination in the hypercube
tresults <- as.data.frame(get.results(myLHS)) #extract results from the LHS object
names(tresults) <- c("p","w")
pdf("../ecology/figures/fig3.pdf")
par(cex.axis=1.5, cex.lab=1.75, las=1, mar = c(5, 6, 4, 2), mgp = c(4,1,0), lwd=2)
plot(w~p, data=tresults, xlab="p-value", ylab="Akaike weight", cex=0.25)
dev.off()
