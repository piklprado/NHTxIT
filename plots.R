load("abacus/tresults.RData")
load("abacus/corresults.RData")
load("abacus/anovaresults.RData")
load("abacus/lmresults.RData")
source("functions.R")

### Graficos exploratorios, com os dois criterios de AIC

Os resultados das 

### Prob rightfull conclusions

par(mfrow = c(2,3))
p1b(t.results, main = "t-test")
p1b(cor.results, legend=FALSE, main = "Correlation")
p1b(anova.results, legend=FALSE, main = "Anova")
p1b(lm.results, legend=FALSE, main = "Linear regression, no collinearity")
p1b(lm.colin.results, legend=FALSE, main = "Linear regression, with collinearity")
par(mfrow = c(1,1))


## Conclusion mismatch

par(mfrow = c(2,3))
p2b(t.results, pos.leg = "topright", main = "t-test")
p2b(cor.results, pos.leg = "topright", legend=FALSE, main = "Correlation")
p2b(anova.results, pos.leg = "topright", legend=FALSE, main = "Anova")
p2b(lm.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, no collinearity")
p2b(lm.colin.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, with collinearity")
par(mfrow = c(1,1))


## Mean M-error

par(mfrow = c(2,3))
p3b(t.results, pos.leg = "topright", main = "t-test")
p3b(cor.results, pos.leg = "topright", legend=FALSE, main = "Correlation")
p3b(anova.results, pos.leg = "topright", legend=FALSE, main = "Anova")
p3b(lm.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, no collinearity")
p3b(lm.colin.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, with collinearity")
par(mfrow = c(1,1))


## prob of S-error

par(mfrow = c(2,3))
p4b(t.results, pos.leg = "topright", main = "t-test")
p4b(cor.results, pos.leg = "topright", legend=FALSE, main = "Correlation")
p4b(anova.results, pos.leg = "topright", legend=FALSE, main = "Anova")
p4b(lm.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, no collinearity")
p4b(lm.colin.results, pos.leg = "topright", legend=FALSE, main = "Linear regression, with collinearity")
par(mfrow = c(1,1))


## P(H0) x Akaike weight
## All simulations
plot(mean.pvalue ~ mean.wH0, data= t.results,
         xlab = "Akaike weight for H0",
     ylab = "p for H0", subset=mean.pvalue <0.1, cex=0.4,
     xlim = c(0,
              max(sapply(list(t.results, cor.results, anova.results, lm.results, lm.colin.results),
                         function(x) max(x$mean.wH0[x$mean.pvalue <0.1]))))
     )
points(mean.pvalue ~ mean.wH0, data= cor.results,
       subset=mean.pvalue <0.1, cex=0.4, col = "red")
points(mean.pvalue ~ mean.wH0, data= anova.results,
       subset=mean.pvalue <0.1, cex=0.4, col = "blue")
points(mean.pvalue ~ mean.wH0, data= lm.results,
       subset=mean.pvalue <0.1, cex=0.4, col = "orange")
legend("topright", c("t-test", "Correlation", "ANOVA", "Linear regression"),
       col = c("black", "red", "blue", "orange"), pch=19, bty="n")


