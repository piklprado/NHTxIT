## Functions used to produce the plots of the paper
## 'Towards a pragmatic use of statistics in ecology, by Leonardo Castilho and Paulo Inacio Prado'
## For any question please use the issue pages at GitHub: https://github.com/piklprado/NHTxIT/issues

source("functions.R")
source("plot.functions.R")
## Load simulation results.
## The R binaries with the simulation results are available under request to the authors,
## but you can reproduce the results running the simmulations again
##(script 'simulation_rightfull_conclusions.R')
## Please correct file paths accordingly
load("simulations/tresults.RData")
load("simulations/corresults.RData")
load("simulations/anovaresults.RData")
load("simulations/lmresults.RData")

## Plot Parameters (from http://shinyapps.org/apps/RGraphCompendium/)
#par(cex.main = 1.5, mar = c(5, 6, 4, 4) + 0.1, mgp = c(3.5, 1, 0), 
#    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 95,
               maxColorValue = 255)
greycol2 <- rgb(0, 0, 0, alpha = 0.7)
grey1 <- grey.colors(10, alpha=0.3)[1]
grey2 <- grey.colors(10, alpha=0.3)[10]

## Predicted relationship for NHT Gaussian designs
## (Gelman & Carlin Perspectives on Psychological Science 2014, Vol. 9(6) 641 –651)
D.seq <- c(seq(0,1,.01),seq(1,10,.1),100)
pGauss <- as.data.frame(t(sapply(D.seq, retrodesign, s=1, alpha = 0.05)))
##
pGauss2 <- as.data.frame(t(sapply(D.seq, retrodesign2, s=1, alpha = 0.05)))


## Probabilities of rightfull conclusions when H0 is false
pdf("../overleaf/figures/prightH1.pdf", width=8, height=8)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(2,2),
    oma=c(3,3,0,0))
p1b(t.results, main = "t-test", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, ylim =c(0,1))
legend("bottomright",c("NHT","IT"), pch=19, col=c("black","darkgrey"), bty="n", cex=1.2)
p1b(cor.results, main = "Correlation", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, ylim =c(0,1))
p1b(anova.results, main = "ANOVA", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, ylim =c(0,1))
p1b(lm.results, main = "Linear regression", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, ylim =c(0,1))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Proportion of rightfull conclusions", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()

## M-errors
pdf("../overleaf/figures/M-error.pdf", width=8, height=8)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(2,2),
    oma=c(3,3,0,0))
p3b(t.results, main = "t-test", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, xlim=c(0.1,2))
legend("topright",c("NHT","IT"), pch=19, col=c("black","darkgrey"), bty="n", cex=1.2)
p3b(cor.results, main = "Correlation", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,2))
p3b(anova.results, main = "ANOVA", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,2))
p3b(lm.results, main = "Linear regression", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,2))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Mean M-error rate", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()

## S-errors
pdf("../overleaf/figures/S-error.pdf", width=8, height=8)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(2,2),
    oma=c(3,3,0,0))
p4b(t.results, main = "t-test", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25, xlim=c(0.1,1.5))
legend("topright",c("NHT","IT"), pch=19, col=c("black","darkgrey"), bty="n", cex=1.2)
p4b(cor.results, main = "Correlation", AIC2=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,1.5))
p4b(anova.results, main = "ANOVA", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,1.5))
p4b(lm.results, main = "Linear regression", AIC1=FALSE, colours=c("black","grey","grey"),
    legend=FALSE, xlab="", ylab="", cex=0.25,  xlim=c(0.1,1.5))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Proportion of S-error", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()

## P x Akaike weights
pdf("../overleaf/figures/pXweights.pdf", width=8, height=8)
par(cex.main = 1.5, lwd=2,
    mar = c(5.5, 5.5, 4, 2) + 0.1,
    mgp = c(3.5, 0.5, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1)
plot(mean.pvalue ~ mean.wH0, data= t.results,
     pch = 3,
     xlab = "Akaike weight for H0",
     ylab = "p for H0", subset=mean.pvalue <0.1, cex=0.4,
     xlim = c(0,
              max(sapply(list(t.results, cor.results, anova.results, lm.results, lm.colin.results),
                         function(x) max(x$mean.wH0[x$mean.pvalue <0.1]))))
     )
points(mean.pvalue ~ mean.wH0, data= cor.results,
       subset=mean.pvalue <0.1, cex=0.4, pch=2, col="darkgrey")
points(mean.pvalue ~ mean.wH0, data= anova.results,
       subset=mean.pvalue <0.1, cex=0.4, pch=1)
points(mean.pvalue ~ mean.wH0, data= lm.results,
       subset=mean.pvalue <0.1, cex=0.4, pch=4, col="darkgrey")
legend("bottomright", c("ANOVA", "Linear regression", "t-test", "Correlation"),
       col = c("black", "darkgrey", "black", "darkgrey"),
       pch=c(1,4,3,2), bty="n")
dev.off()



### Supplementary figures ###

## Critical value of t and its IT correspondent 
## Minimum value of the likelihhod-ratio to reject H0 expressed in units of t (expression 9.5.3 from Edawrds, 1972)
st <- function(n, m) sqrt((n-1)*(exp(2*m/n)-1))
## Critical t for a certain significance value, from the t distribution
cp <- function(n, alpha) qt(p= 1-alpha/2, df = n-1)
## Fig S1
pdf("../overleaf/figures/figS1.pdf")
par(cex.axis=1.5, cex.lab=1.75, las=1, mar = c(5, 6, 4, 2), mgp = c(4,1,0), lwd=1.5)
curve(cp(x, alpha = 0.05), 3, 30, xlab= "sample size", ylab="| t | value to reject H0", lty=2, lwd=2)
curve(st(x, m=2), add=TRUE, lwd = 2)
legend("topright", c("t-test", "IT test"), lty=2:1, bty="n", cex=1.5)
dev.off()


## Relationship between t-value and evidence weights ##
x <- seq(-4,5.8,by=0.01)
y1 <- dt(x, df=20)
y2 <- dt(x, df=20, ncp = 1)
a <- c(0.7,2.5)
pdf("../overleaf/figures/weights_and_ttest.pdf", width=8, height=6.5)
par(cex= 1.25, lwd=2)
#par(cex.main = 1.5, mar = c(5, 6, 4, 4) + 0.1, mgp = c(3.5, 1, 0), 
#    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, cex=2)
plot(x, y1, type="n", axes=FALSE, xlab="", ylab="")
polygon(c(a[1],x[x>a[1]]), c(0,y1[x>a[1]]), col=grey1, lty=0)
polygon(c(a[2],x[x>a[2]]), c(0,y1[x>a[2]]), col=grey2, lty=0)
segments(x0=c(a[1], a[1], a[2], a[2]),
         y0=c(dt(a[1], df=20), dt(a[1], df=20, ncp=1), dt(a[2], df=20), dt(a[2], df=20, ncp=1)),
         x1=c(min(x)*.6, min(x)*.6, max(x)*.85, max(x)*.85), lty=2,
         lwd=1.5,
         col=c("darkgrey", "black", "darkgrey", "black"))
lines(x, y2, lwd=3)
lines(x, y1, lwd=3, col="darkgrey")
text(label=c("a", "b"), x= a, y=c(-0.01, -0.01))
arrows(x0=min(x)*.6, y0=dt(a[1], df=20), y1=dt(a[1], df=20, ncp=1),
       angle=90, code=3, length=0.05, lwd=1.5)
text(x=min(x)*.6, y=c(dt(a[1], df=20), dt(a[1], df=20, ncp=1)),
     labels=c("H0(a)", "H1(a)"), pos=2)
arrows(x0=max(x)*.85, y0=dt(a[2], df=20), y1=dt(a[2], df=20, ncp=1),
       angle=90, code=3, length=0.05, lwd=1.5)
text(x=max(x)*.85, y=c(dt(a[2], df=20), dt(a[2], df=20, ncp=1)),
     labels=c("H0(b)", "H1(b)"), pos=4)
dev.off()


## P rightfull conclusions for both IT criteria
pdf("../overleaf/figures/Sup_prightH1.pdf", width=12, height=4.5)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(1,3),
    oma=c(3,3,0,0))
p1b(anova.results, main = "ANOVA", 
    legend=FALSE, xlab="", ylab="", cex=0.25,
    colours=rainbow(3, alpha=0.2), ylim =c(0,1))
legend("bottomright",c("NHT","IT", "IT, adjusted"), pch=19,
       col=rainbow(3, alpha=0.5), bty="n", cex=1.2)
p1b(lm.results, main = "Linear regression", 
    legend=FALSE, xlab="", ylab="", cex=0.25,
    colours=rainbow(3, alpha=0.2), ylim =c(0,1))
p1b(lm.colin.results, main = "Linear regression, collinearity", 
    legend=FALSE, xlab="", ylab="", cex=0.25,
    colours=rainbow(3, alpha=0.2), ylim =c(0,1))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Prop. rightfull conclusions", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()


## M-Errors
pdf("../overleaf/figures/Sup_M-error.pdf", width=12, height=4.5)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(1,3),
    oma=c(3,3,0,0))
p3b(anova.results, main = "ANOVA", 
    legend=FALSE, xlab="", ylab="", cex=0.25, xlim =c(0.1,3),
    colours=rainbow(3, alpha=0.2))
legend("topright",c("NHT","IT", "IT, adjusted"), pch=19,
       col=rainbow(3, alpha=0.5), bty="n", cex=1.2)
p3b(lm.results, main = "Linear regression", 
    legend=FALSE, xlab="", ylab="", cex=0.25, xlim =c(0.1,3),
    colours=rainbow(3, alpha=0.2))
p3b(lm.colin.results, main = "Linear regression, collinearity",
    legend=FALSE, xlab="", ylab="", cex=0.25, xlim =c(0.1,3),
    colours=rainbow(3, alpha=0.2))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Mean M-error rate", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()

## S-Errors
pdf("../overleaf/figures/Sup_S-error.pdf", width=12, height=4.5)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1, mfrow=c(1,3),
    oma=c(3,3,0,0))
p4b(anova.results, main = "ANOVA", 
    legend=FALSE, xlab="", ylab="", cex=0.5, xlim =c(0.1,1.5),
    colours=rainbow(3, alpha=0.5))
legend("topright",c("NHT","IT", "IT, adjusted"), pch=19,
       col=rainbow(3, alpha=0.5), bty="n", cex=1.2)
p4b(lm.results, main = "Linear regression", 
    legend=FALSE, xlab="", ylab="", cex=0.5, xlim =c(0.1,1.5),
    colours=rainbow(3, alpha=0.5))
p4b(lm.colin.results, main = "Linear regression, collinearity",
    legend=FALSE, xlab="", ylab="", cex=0.5, xlim =c(0.1,1.5),
    colours=rainbow(3, alpha=0.5))
mtext("Effect Size", side=1, outer=TRUE, cex=2, line=0)
mtext("Proportion of S-error", side=2, outer=TRUE, cex=2, line=0.5, las=0)
dev.off()
