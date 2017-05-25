library(ggplot2)
################################################################################
## t-test
################################################################################
load("tresults.RData")
## proportion of rightfull conclusions
par(mfrow=c(1,2))
plot(p.NHT.right ~ D, data=t.results, cex=0.25)
plot(p.AIC.right ~ D, data=t.results, cex=0.25)
par(mfrow=c(1,1))
plot(p.NHT.right ~ D, data=t.results, cex=0.25, ylim=range(c(p.NHT.right,p.AIC.right)),
     xlab="Effect size", ylab="P rightfull conclusions")
points(p.AIC.right ~ D, data=t.results, cex=0.25, col="blue")
points(p.AIC.right ~ D, data=t.results.d2, cex=0.25, col="red")
##
plot(p.NHT.right ~ p.AIC.right, data=t.results, cex=0.25)
abline(0,1)
##Odds-ratio of rightfull conclusions NHT x IT
plot(I(p.NHT.right/p.AIC.right) ~ D, data=t.results, cex=0.25)
## Proportion of conclusion mismatches
plot(p.mismatch ~ D, data=t.results, cex=0.25)
plot(p.mismatch ~ D, data=t.results.d2, cex=0.25)
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
## correlation
################################################################################
load("corresults.RData")
## proportion of rightfull conclusions
par(mfrow=c(1,2))
plot(p.NHT.right ~ D, data=cor.results, cex=0.25)
plot(p.AIC.right ~ D, data=cor.results, cex=0.25)
par(mfrow=c(1,1))
plot(p.NHT.right ~ D, data=cor.results, cex=0.25)
points(p.AIC.right ~ D, data=cor.results, cex=0.25, col="blue")
plot(p.NHT.right ~ rpears, data=cor.results, cex=0.25)
points(p.AIC.right ~ rpears, data=cor.results, cex=0.25, col="blue")
## transformção do r em variavel z
plot(p.NHT.right ~ I(atanh(rpears)*sqrt(N-3)), data=cor.results, cex=0.25, xlab="z") 
points(p.AIC.right ~ I(atanh(rpears)*sqrt(N-3)), data=cor.results, cex=0.25, col="blue")
plot(p.NHT.right ~ p.AIC.right, data=cor.results, cex=0.25)
abline(0,1)
    ggplot(cor.results , aes(D, p.AIC.right)) +
    geom_point(aes(colour=N))
plot(p.AIC.right ~ I(atanh(rpears)*sqrt(N-3)), data=cor.results, cex=0.25) # z-scores

###
pairs(cor.results[,1:5])
###
plot(p.NHT.right ~ D, data=cor.results,
     xlab="Standardized effect", ylab="P rightful conclusion",
     cex=0.25, ylim=range(c(p.NHT.right,p.AIC.right)))
points(p.AIC.right ~ D, data=cor.results, col="blue", cex=0.25)
points(p.AIC.right ~ D, data=cor.results.d2, col="red", cex=0.25)
legend("topleft", c("NHT", "IT", "IT, d>2"), col=c("black","blue", "red"), pch=1, bty="n")
## Investigando a regiao de maior de dicrepancia: 0 < t < 5
f1 <- function(tval, N) tval / sqrt(N-2+tval^2)
f1 <- Vectorize(f1)
x.tval <- seq(0.1,5.5, length=200)
N.tval <- 10:100
Y <- outer(x.tval,N.tval, f1)
contour(x.tval,N.tval,Y)
## Ou o contrario, verificando a isolinha de tval > 5
curve(5/sqrt(x-2+5^2), 5,500)
## ou zval = 6
curve(tanh(6/sqrt(x-3)), 5,500)
##Odds-ratio of rightfull conclusions NHT x IT
plot(I(p.NHT.right/p.AIC.right) ~ D, data=cor.results, cex=0.25)
## Proportion of conclusion mismatches
plot(p.mismatch ~ D, data=cor.results, cex=0.25)
plot(p.mismatch ~ rpears, data=cor.results, cex=0.25)
plot(p.mismatch ~ I(atanh(rpears)*sqrt(N-3)), data=cor.results, cex=0.25) ## z-score
plot(p.mismatch ~ D, data=cor.results.d2, cex=0.25)
plot(p.mismatch ~ N, data=cor.results.d2, cex=0.25)
## Odd ratio rightfull x mismatches
plot(I(p.NHT.right/p.AIC.right) ~ I(p.mismatch/(1-p.mismatch)), data=cor.results, cex=0.25)

## Type-M errors
plot(I(mean.NHT.M/mean.AIC.M) ~ D, data=cor.results, cex=0.25)
abline(h=1, lty=2)
plot(mean.NHT.M ~ D, data=cor.results,
     ylim=range(c(mean.NHT.M,mean.AIC.M)), cex=0.25)
points(mean.AIC.M ~ D, data=cor.results, col="blue", cex=0.25)
## P-valus x w
plot(mean.pvalue~mean.wH0, data=cor.results, cex=0.25)

## A model on logits
cor.m1 <- glm(cbind(p.AIC.right*1e4, (1-p.AIC.right)*1e4) ~ D + N + SD , data=cor.results, family=binomial)
cor.lm1 <- lm(log((p.AIC.right+1/1e4)/(1-p.AIC.right+1/1e4)) ~ D + N + SD , data=cor.results)
cor.lm2 <- lm(log((p.AIC.right+1/1e4)/(1-p.AIC.right+1/1e4)) ~ D * N + SD , data=cor.results)
summary(cor.m1)
avPlots(cor.m1)
avPlots(cor.lm1)
avPlots(cor.lm2)

