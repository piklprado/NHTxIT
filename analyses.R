################################################################################
## t-test
################################################################################
load("tresults.RData")
## proportion of rightfull conclusions
par(mfrow=c(1,2))
plot(p.NHT.right ~ D, data=t.results, cex=0.25)
plot(p.AIC.right ~ D, data=t.results, cex=0.25)
par(mfrow=c(1,1))
plot(p.NHT.right ~ D, data=t.results, cex=0.25)
points(p.AIC.right ~ D, data=t.results, cex=0.25, col="blue")
plot(p.NHT.right ~ p.AIC.right, data=t.results, cex=0.25)
abline(0,1)
##Odds-ratio of rightfull conclusions NHT x IT
plot(I(p.NHT.right/p.AIC.right) ~ D, data=t.results, cex=0.25)
## Proportion of conclusion mismatches
plot(p.mismatch ~ D, data=t.results, cex=0.25)
## Odd ratio rightfull x mismatches
plot(I(p.NHT.right/p.AIC.right) ~ I(p.mismatch/(1-p.mismatch)), data=t.results, cex=0.25)
## Type-M errors
plot(I(mean.NHT.M/mean.AIC.M) ~ D, data=t.results, cex=0.25)
abline(h=1, lty=2)
plot(mean.NHT.M ~ D, data=t.results,
     ylim=range(c(mean.NHT.M,mean.AIC.M)), cex=0.25, xlim=c(0,5))
points(mean.AIC.M ~ D, data=t.results, col="blue", cex=0.25)
## P-valus x w
plot(mean.pvalue~mean.wH0, data=t.results, cex=0.25)

################################################################################
## t-test
################################################################################
load("corresults.RData")
## proportion of rightfull conclusions
par(mfrow=c(1,2))
plot(p.NHT.right ~ D, data=cor.results, cex=0.25)
plot(p.AIC.right ~ D, data=cor.results, cex=0.25)
par(mfrow=c(1,1))
plot(p.NHT.right ~ D, data=cor.results, cex=0.25)
points(p.AIC.right ~ D, data=cor.results, cex=0.25, col="blue")
plot(p.NHT.right ~ p.AIC.right, data=cor.results, cex=0.25)
abline(0,1)
###
plot(p.NHT.right ~ D, data=cor.results,
     xlab="Standardized effect", ylab="P rightful conclusion",
     cex=0.25)
points(p.AIC.right ~ D, data=cor.results, col="blue", cex=0.25)
legend("topleft", c("NHT", "IT"), col=c("black","blue"), pch=1, bty="n")
## Investigando a regiao de maior de dicrepancia: 0 < t < 5
f1 <- function(tval, N) tval / sqrt(N-2+tval^2)
f1 <- Vectorize(f1)
x.tval <- seq(0.1,5.5, length=200)
N.tval <- 10:100
Y <- outer(x.tval,N.tval, f1)
contour(x.tval,N.tval,Y)
##Odds-ratio of rightfull conclusions NHT x IT
plot(I(p.NHT.right/p.AIC.right) ~ D, data=cor.results, cex=0.25)
## Proportion of conclusion mismatches
plot(p.mismatch ~ D, data=cor.results, cex=0.25)
## Odd ratio rightfull x mismatches
plot(I(p.NHT.right/p.AIC.right) ~ I(p.mismatch/(1-p.mismatch)), data=cor.results, cex=0.25)
## Type-M errors
plot(I(mean.NHT.M/mean.AIC.M) ~ D, data=cor.results, cex=0.25)
abline(h=1, lty=2)
plot(mean.NHT.M ~ D, data=cor.results,
     ylim=range(c(mean.NHT.M,mean.AIC.M)), cex=0.25, xlim=c(0,5))
points(mean.AIC.M ~ D, data=cor.results, col="blue", cex=0.25)
## P-valus x w
plot(mean.pvalue~mean.wH0, data=cor.results, cex=0.25)

