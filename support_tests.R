## Exemplos seção 9.5 de Edwards
## Teste t de uma amostra de que a amostra vem de uma distribuição co média mu
## Valor de t que dá um suporte (likelihood ratio) m (expressão 9.5.3)
st <- function(n, m) sqrt((n-1)*(exp(2*m/n)-1))
## Valor do t critico para um certo alfa
cp <- function(n, alfa) qt(p= 1-alfa, df = n-1)
## Fig 1
pdf("../ecology/figures/fig1.pdf")
par(cex.axis=1.5, cex.lab=1.75, las=1, mar = c(5, 6, 4, 2), mgp = c(4,1,0), lwd=1.5)
curve(cp(x, alfa = 0.025), 3, 30, xlab= "sample size", ylab="| t | value to reject H0", lty=2, lwd=2)
curve(st(x, m=2), add=TRUE, lwd = 2)
legend("topright", c("t-test", "IT test"), lty=2:1, bty="n", cex=1.5)
dev.off()

## Fig 3: duas distribuições de t para um t-value = 2.5 df=10
x <- seq(-4,5,by=0.01)
y1 <- dt(x, df=20)
y2 <- dt(x, df=20, ncp = 1)
a <- c(0.7,2.5)

pdf("../ecology/figures/fig4.pdf", width=8, height=6.5)
par(cex.axis=1.5, cex.lab=1.75, las=1, mar = c(5, 6, 4, 5), mgp = c(4,1,0), lwd=2)
plot(x, y1, type="n", axes=FALSE, xlab="", ylab="")
polygon(c(a[1],x[x>a[1]]), c(0,y1[x>a[1]]), col=grey.colors(10, alpha=0.5)[10], lty=0)
polygon(c(a[2],x[x>a[2]]), c(0,y1[x>a[2]]), col=grey.colors(10, alpha=0.5)[1], lty=0)
#segments(a, y0 = c(0,0), x1=a, y1 = dt(a, df=20), lwd=5, col="grey")
lines(x, y2, lwd=2)
lines(x, y1, lwd=2, col="darkgrey")
##segments(a, y0 = c(0,0), x1=a, y1 = dt(a, ncp=1, df=20), lwd=1.5)
axis(1, at=a, labels=c("a", "b"))
axis(2, at=c(dt(a[1], df=20),dt(a[1], df=20, ncp=1)), labels=c("H0(a)", "H1(a)"))
axis(4, at=c(dt(a[2], df=20),dt(a[2], df=20, ncp=1)), labels=c("H0(b)", "H1(b)"))
segments(x0=c(a[1], a[1], a[2], a[2]),
         y0=c(dt(a[1], df=20), dt(a[1], df=20, ncp=1), dt(a[2], df=20), dt(a[2], df=20, ncp=1)),
         x1=c(min(x)-1, min(x)-1, max(x)+.5, max(x))+.5, lty=2, col=c("darkgrey", "black", "darkgrey", "black"))
box()
dev.off()
