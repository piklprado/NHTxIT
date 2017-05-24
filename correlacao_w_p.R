library(bbmle)
library(pse)
library(mvtnorm)

## Relationship between p-values an evidence weights in a t-test

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
tresults <- as.data.frame(get.results(myLHS)) #extract results from
                                        #the LHS object
names(tresults) <- c("p","w")
pdf("../ecology/figures/fig3.pdf")
par(cex.axis=1.5, cex.lab=1.75, las=1, mar = c(5, 6, 4, 2), mgp = c(4,1,0), lwd=2)
plot(w~p, data=tresults, xlab="p-value", ylab="Akaike weight", cex=0.25)
dev.off()

#' p-value and evidence weight from correlation tests between two samples
#'
#' This function simulates samples from bivariate Gaussian distribution
#' with equal marginal variances. The hypothesis that there is a
#' correlation between the value of the samples the samples come from
#' the  is evaluated with a significance corrleation test and also comparing
#' the AIC of the fit of the data to a bicariate Gaussian with the
#' covariance parameter free or set to zero.
#'
#' @param rpears real number, the true Pearson correlation between the two
#'     variables sampled 
#' @param sd1 real positive, marginal of the standard deviations of the Gaussian
#'     populations from which the samples are drawn.
#' @param N integer positive, size of the two samples.
#' @return a vector with of size 2: the p-value of a correlation
#'     t-test (\code{cor.test} for the
#'     null hypothesis that the two samples come from independent
#'     distributions (zero corrleation); and the Akaike evidence weight of the null
#'     hypothesis compared to the alternative model that the
#'     samples come from a bivariate Gaussian with some degree of correlation.
#' @details This function simulates samples of the same size \code{N}
#'     drawn from a bivariate Gaussian distribution with the same
#'     marginal standard
#'     deviation for both samples \code{sd1} and corrleation given by \code{corr}. The mean of both
#'     distributions are set to zero.
#'     The null hypothesis that the samples come from
#'     independent distributions (or from a bivariate normal with no
#'     correlation, wich is the same) is then evaluated by a bicaudal
#'     t-test (see \code{cor.test}}
#'     and by model selection AIC corrected for small samples
#'     (AICc). The two models compared are the fits of a bivariate
#'     Gaussian to the data with a free parameter for correlation of
#'     with correlation fixed to zero. The evidence wieght is then
#'     calculated from the AIC difference of the two models.
wp.cortest <- function(rpears, sd1, N){
    cov1 <- rpears*sd1^2
    sigma <- matrix(c(sd1^2, cov1, cov1, sd1^2), nrow=2)
    x <- rmvnorm(n = N, sigma=sigma)
    pt <- cor.test(x[,1], x[,2])$p.value
    cov3 <- cov2 <- cov(x)
    cov3[2,1] <- 0
    cov3[1,2] <- 0
    mean1 <- apply(x, 2, mean)
    m1 <- sum(dmvnorm(x , mean = mean1, sigma = cov3, log=TRUE))
    m2 <- sum(dmvnorm(x , mean = mean1, sigma = cov2, log=TRUE))
    AICc.m1 <- -2*m1 + 8 + 40/(N - 3)
    AICc.m2 <- -2*m2 + 12 + 84/(N - 5)
    AICs <- c(AICc.m1, AICc.m2)
    dAICs <- AICs - min(AICs)
    w <- exp(-dAICs/2)
    wi <- w/sum(w)
    return(c(pt, wi[1]))
}

## A mesma funcao que retorna o p e a diferenca entre AIC do modelo nulo em relacao ao alternativo
wp.cortest.aic <- function(rpears, sd1, N){
    cov1 <- rpears*sd1^2
    sigma <- matrix(c(sd1^2, cov1, cov1, sd1^2), nrow=2)
    x <- rmvnorm(n = N, sigma=sigma)
    pt <- cor.test(x[,1], x[,2])$p.value
    ## Matriz de covariancia nas amostras: tem a estimativa de maxima verossimilhanca
    cov3 <- cov2 <- cov(x)
    cov3[2,1] <- 0
    cov3[1,2] <- 0
    mean1 <- apply(x, 2, mean)
    ## Modelo nulo: ajusta um normal bivariada com correlacao = 0
    m1 <- sum(dmvnorm(x , mean = mean1, sigma = cov3, log=TRUE))
    ## Modelo alternativo: ajusta normal biv com parametro de correlacao livre
    m2 <- sum(dmvnorm(x , mean = mean1, sigma = cov2, log=TRUE))
    AICc.m1 <- -2*m1 + 8 + 40/(N - 3)
    AICc.m2 <- -2*m2 + 12 + 84/(N - 5)
    ## Diferenca entre AIC modelo nulo em relacao ao alternativo
    delta <- AICc.m1 - AICc.m2 
    return(c(pt, delta))
}

## Um exemplo de 1000 simulacoes do coeficiente de regressao

## Repete um loop que armazena cada resultado em uma linha da matriz
## Cria um data.frame com todas as combinacoes possivels dos parametros da simulacao
rdf <- expand.grid(r=c(0, 0.5, 0.7, 0.8), sd1= 1, N = c(10,20,30))
## cria uma matriz para os resultados
r1 <- matrix(nrow=1000*nrow(rdf), ncol=2)
for(i in 1:nrow(r1)){
    r1[i, ] <- wp.cortest.aic(rpears = 0.8, sd1 = 1, N = 30)
}

## primeiras linhas da matriz de resultados: primeira coluna é p, segunda deltaAICc
head(r1)

## Hypercube sampling (see introductory vignette of pse package)
## To run the sampling simulation many times
modelcor <- function (my.data) {
 return(mapply(wp.corest, my.data[,1], my.data[,2], my.data[,3]))
}

## Hypercube sample of 1000 combinations of parameters
myLHS2 <- LHS(modelcor,
              factors = c("rpears", "sd1", "N"),
              N = 2000,
              q = c("qunif", "qunif", "qdunif"),
              q.arg = list(list(min=0, max=1),
                           list(min=1, max=5),
                           list(min=10, max=100))
              )
## Relationship between p-value and weight of evidence for each
## parameter combination in the hypercube
tresults2 <- as.data.frame(get.results(myLHS2)) #extract results from the LHS object
names(tresults2) <- c("p","w")
plot(w~p, xlab="p", ylab="w", data=tresults2, subset=p<=0.1, cex=0.2)
## Leo, aqui o padrão não é tão simples como o para o teste t. Até
## merece ser mais explorado, mas para os propósitos deste paper
## usaria o caso simples do t para mostrar a relação entre p-value e
## peso de evidência. De qq forma a função wp.cortest serve para vc
## refazer as comparações do deisgn de correlacao, como sugeri no texto.
