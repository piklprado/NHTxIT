## t-test ##
## P.rightfull
plot(p.NHT.right ~ D, data=t.results,
     xlab="Standardized effect", ylab="P rightful conclusion",
     cex=0.25, ylim=range(c(p.NHT.right,p.AIC.right)))
points(p.AIC.right ~ D, data=t.results, col="blue", cex=0.25)
legend("topleft", c("NHT", "IT"), col=c("black","blue"), pch=1, bty="n")
## Odds ratio P.rightfull
plot(I(p.NHT.right/p.AIC.right) ~ D, data=t.results,
     xlab="Standardized effect", ylab="P right NHT/ P right IT",
     cex=0.25)
abline(h=1, lty=2, col="blue")
## P mismatch
plot(p.mismatch ~ D, data=t.results,
     xlab="Standardized effect", ylab="P mismatch",
     cex=0.25)
## M-Errors
plot(mean.NHT.M ~ D, data=t.results,
     xlab="Standardized effect", ylab="Mean type-M error",
     cex=0.25)
points(mean.AIC.M ~ D, data=t.results, col="blue", cex=0.25)
## M-errors ratio
## P mismatch
plot(I(mean.NHT.M/mean.AIC.M) ~ D, data=t.results,
     xlab="Standardized effect", ylab="NHT/IT type-M error ratio",
     cex=0.25)
abline(h=1, lty=2)

## Correlation ##
## P.rightfull
plot(p.NHT.right ~ D, data=cor.results,
     xlab="Standardized effect", ylab="P rightful conclusion",
     cex=0.25, ylim=range(c(p.NHT.right,p.AIC.right)))
points(p.AIC.right ~ D, data=cor.results, col="blue", cex=0.25)
legend("topleft", c("NHT", "IT"), col=c("black","blue"), pch=1, bty="n")
## With correlations instead of standardize t effects
plot(p.NHT.right ~ rpears, data=cor.results,
     xlab="Standardized effect", ylab="P rightful conclusion",
     cex=0.25, ylim=range(c(p.NHT.right,p.AIC.right)))
points(p.AIC.right ~ rpears, data=cor.results, col="blue", cex=0.25)
legend("topleft", c("NHT", "IT"), col=c("black","blue"), pch=1, bty="n")
## Odds ratio P.rightfull
plot(I(p.NHT.right/p.AIC.right) ~ D, data=cor.results,
     xlab="Standardized effect", ylab="P right NHT/ P right IT",
     cex=0.25)
abline(h=1, lty=2, col="blue")

## P mismatch
plot(p.mismatch ~ D, data=cor.results,
     xlab="Standardized effect", ylab="P mismatch",
     cex=0.25)
## M-Errors
plot(mean.NHT.M ~ D, data=cor.results,
     xlab="Standardized effect", ylab="Mean type-M error",
     cex=0.25)
points(mean.AIC.M ~ D, data=cor.results, col="blue", cex=0.25)
## M-errors ratio
## P mismatch
plot(I(mean.NHT.M/mean.AIC.M) ~ D, data=cor.results,
     xlab="Standardized effect", ylab="NHT/IT type-M error ratio",
     cex=0.25)
abline(h=1, lty=2)
