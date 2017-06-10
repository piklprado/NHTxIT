library(bbmle)
library(pse)
library(mvtnorm)

#' NHT and IT  comparison for t-test of simulated samples
#'
#' This function simulates samples from two Gaussian distributions
#' with equal variance. The hypothesis that the samples come from the
#' same distribution is evaluated with a t-test and also comparing
#' the AIC of linear models with or without a group label as
#' predictors.
#'
#' @param t.val real number, the true difference or effect size
#'     expressesed as a t-value (difference between true means divided
#'     by the standard error of the difference). Ignored if \code{mu2}
#'     is provided
#' @param sd1 real positive, standard deviation of the Gaussian
#'     populations from which the samples are drawn.
#' @param N integer positive, size of the two samples.
#' @param mu2 real number, the mean of one of the Gaussian populations
#'     from which the samples are drawn.  The mean of the other
#'     population is set to zero.
#' @param delta2 logical; should null hypothesis be rejected only if it has a value of 
#'     deltaAIC larger than two? If FALSE the null
#'     hypothesis is rejected if it has a larger deltaAIC value than the
#'     alternative, even if deltaAIC<2.
#' @return a vector of length 7: the p-value of a two-tailed t-test
#'     for the null hypothesis that the two samples come from the same
#'     distributions; the Akaike evidence weight of the null
#'     hypothesis compared to the alternative model that each sample
#'     comes from a distribution with a different mean; the difference
#'     between Akaike Information Criteria of the null to the
#'     alternative hypothesis; the ratio between the estimated
#'     difference among groups and the true difference, if the
#'     difference when significant is in the opposite direction (M and
#'     S-errorr, Gelman & Carlin (2014), if the null hypothesis tests
#'     reached to the right conclusion (yes=1, no=0), if the model
#'     selection based on AICs reached to the right conclusion.
#' @details This function simulates samples of the same size \code{N}
#'     drawn from Gaussian distributions with the same standard
#'     deviation \code{sd1}. The mean of the two distributions differ
#'     according to the true mean of one group or to the t-value
#'     indicated by the argument \code{t.val}, given that the mean of
#'     one group is zero.  The null hypothesis that the samples come
#'     from the same distribution is then evaluated by a bicaudal
#'     t-test and by model selection AIC corrected for small samples
#'     (AICc). The two models compared are that the expected value of
#'     all values are the same or that the the expected values of each
#'     sample are different. The evidence weight and delta-AIC of the
#'     null hypothesis is then calculated from the AIC difference of
#'     the two models. Also, the M-error (Gelman & Carlin 2014) is
#'     calculated dividing the estimated difference among groups by
#'     the true difference.
#' @references Gelman, A., & Carlin, J. (2014). Beyond power
#'     calculations assessing type s (sign) and type m (magnitude)
#'     errors. Perspectives on Psychological Science, 9(6), 641-651.
wp.ttest <- function(tval, sd1, N, mu2, delta2=FALSE){
    # sd1 <- se1*sqrt(2*N)
    if(missing(mu2))
        delta <- tval*sd1*sqrt(2/N)
    if(missing(tval))
        delta <-  mu2
    x1 <- rnorm(N, sd = sd1)
    x2 <- rnorm(N, mean = delta, sd = sd1)
    pvalue <- t.test(x1, x2, var.equal = TRUE)$p.value
    y <- factor(rep(c("am1","am2"), each = N))
    m1 <- lm(c(x1,x2) ~ 1)
    m2 <- lm(c(x1,x2)~ y)
    AICs <- c(AICc(m1, nobs=2*N), AICc(m2, nobs=2*N))
    dAICs <- AICs - min(AICs)
    if(delta==0){
        S.error <- NA
        M.error <- NA
        nht.right <- pvalue >= 0.05
        if(delta2)
            aic.right <- dAICs[2] > 2
        else
            aic.right <- dAICs[1] < dAICs[2] 
    }
    else{
        M.error <- coef(m2)[2]/delta
        S.error <- M.error < 0
        nht.right <- pvalue < 0.05
        if(delta2)
            aic.right <- dAICs[1] > 2
        else
            aic.right <- dAICs[1] > dAICs[2] 
    }
    w <- exp(-dAICs/2)
    wi <- w/sum(w)
    results <- c(pvalue, wi[1], dAICs[1], M.error, S.error, nht.right,
                 aic.right)
    names(results) <- c("p.value", "Akaike.weight.H0", "deltaAIC.H0",
                        "M.error", "S.error", "NHT.right",
                        "AIC.right")
    return(results)
}


#' NHT and IT comparison for correlation tests between two samples
#'
#' This function simulates samples from bivariate Gaussian distribution
#' with equal marginal variances. The hypothesis that there is a
#' correlation between the two variables is evaluated with a
#' significance correlation test and also comparing 
#' the AIC of the fit of the data to a bivariate Gaussian with the
#' covariance parameter free or set to zero.
#' @param tval real number, the true value of  correlation coefficient
#' expressed as a t-value (t = r sqrt((n-2)/(1-r^2))).
#' @param sd1 real positive, standard
#'     deviations of the marginal Gaussian
#'     populations from which the samples are drawn.
#' @param N integer positive, size of the two samples.
#' @param rpears real number, the true Pearson correlation between the two
#'     variables sampled. If missing calculated from \code{tval}.
#' @param delta2 logical; should null hypothesis be rejected only if it has a value of 
#'     deltaAIC larger than two? If FALSE the null
#'     hypothesis is rejected if it has a larger deltaAIC value than the
#'     alternative, even if deltaAIC<2.
#' @return a vector with of size 7: the p-value of a correlation
#'     t-test (\code{cor.test} for the
#'     null hypothesis that the two samples come from independent
#'     distributions (zero correlation); and the Akaike evidence weight of the null
#'     hypothesis compared to the alternative model that the samples come from a bivariate Gaussian with some degree of 
#'     correlation; the difference between Akaike Information Criteria of the null
#'     to the alternative hypothesis; and the ratio between the
#'     estimated correlation and the true value if
#'     the difference when significant is in the opposite direction (S
#'     and M type errors, Gelman & Carlin, 2014) if the null hypothesis tests
#'     reached to the right conclusion (yes=1,
#'     no=0), if the model selection based on AICs reached to the right conclusion.
#' @details This function simulates samples of the same size \code{N}
#'     drawn from a bivariate Gaussian distribution with the same
#'     marginal standard
#'     deviation for both samples (calculated from \code{se1} and
#'     \code{N}) and correlation given by \code{corr}. The mean of both
#'     distributions are set to zero.
#'     The null hypothesis that the samples come from
#'     independent distributions (or from a bivariate normal with no
#'     correlation, wich is the same) is then evaluated by a bicaudal
#'     t-test (see \code{cor.test}}
#'     and by model selection AIC corrected for small samples
#'     (AICc). The two models compared are the fits of a bivariate
#'     Gaussian to the data with a free parameter for correlation of
#'     with correlation fixed to zero. The evidence weight is then
#'     calculated from the AIC difference of the two models.
#' @references Gelman, A., & Carlin, J. (2014). Beyond power
#'     calculations assessing type s (sign) and type m (magnitude)
#'     errors. Perspectives on Psychological Science, 9(6), 641-651. 
wp.cortest <- function(tval, sd1, N, rpears, delta2=FALSE){
    if(missing(rpears))
        rpears <- tval/sqrt(N - 2 + tval^2)
    ## sd1 <- se1*sqrt(N)
    cov1 <- rpears*sd1^2
    sigma <- matrix(c(sd1^2, cov1, cov1, sd1^2), nrow=2)
    x <- rmvnorm(n = N, sigma=sigma)
    p.value <- cor.test(x[,1], x[,2])$p.value
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
    if(rpears == 0){
        M.error <-  NA
        S.error <-  NA
        nht.right <- p.value >= 0.05
        if(delta2)
            aic.right <- dAICs[2] > 2
        else
            aic.right <- dAICs[1] < dAICs[2] 
    }
    else{
        M.error <-  cor(x[,1],x[,2])/rpears
        S.error <- M.error < 0
        nht.right <- p.value < 0.05
        if(delta2)
            aic.right <- dAICs[1] > 2
        else
            aic.right <- dAICs[1] > dAICs[2]
        }
    w <- exp(-dAICs/2)
    wi <- w/sum(w)
    results <- c(p.value, wi[1], dAICs[1], M.error, S.error,
                 nht.right, aic.right)
    names(results) <- c("p.value", "Akaike.weight.H0", "deltaAIC.H0",
                        "M.error", "S.error", "NHT.right",
                        "AIC.right")
    return(results)
}

#' NHT and IT  comparison for F-test of three simulated samples
#'
#' This function simulates samples from three Gaussian distributions
#' with equal variance. The hypothesis that the samples come from the
#' same distribution is evaluated with a F-test. If the null
#' hypothesis is rejected (p<0.05) the differences among samples is
#' then tested with a
#' post-hoc Tukey HSD test. The differences among groups is also tested comparing
#' the AIC of linear models with or without a group label as
#' predictors.
#'
#' @param std.diff real number, the true difference or effect size
#'     expresseed as a t-value (difference between true mean of tretament one and treatment 2 divided
#'     by the standard error of the difference). Ignored if \code{mu1}
#'     is provided
#' @param mu1 real number, the mean of one of the Gaussian populations from which the samples are drawn.
#' The mean of the other two populations is set to zero.
#' @param se1 real positive, standard error of the effect, from which
#'     is calculated 
#'     the standard deviation of the Gaussian
#'     populations from which the samples are drawn.
#' @param N integer positive, size of the three samples.
#' @param delta2 logical; should null hypothesis be rejected only if it has a
#'     deltaAIC larger than two? If FALSE the null
#'     hypothesis is rejected if it has larger deltaAIC than
#'     the alternative hypothesis.
#' @return a vector of length 7: the p-value of a F-test for the
#'     null hypothesis that the three samples come from the same
#'     distribution; the Akaike evidence weight of the null
#'     hypothesis compared to the alternative models that at least one
#'     sample comes from a distribution with a different mean;
#'     the difference between Akaike Information Criteria of the null
#'     to the most plausible hypothesis; the ratio between the
#'     estimated difference among group 1 and the other two and the true difference if
#'     the difference when significant is in the opposite direction (S
#'     and M type errors, Gelman & Carlin, 2014), if the null hypothesis tests
#'     (F-test + Tukey HSD) reached to the right conclusion (yes=1,
#'     no=0), if the model selection based on AICs reached to the right conclusion.
#' @details This function simulates samples of the same size \code{N}
#'     drawn from Gaussian distributions with the same standard
#'     deviation \code{sd1}. The mean of one of the distributions distributions differ
#'     according to the standardized effect (or alternatively the true mean) of one group indicated by the argument 
#'     \code{mu1}, given that the mean of the other two groups is zero.
#'     The null hypothesis that the samples come from
#'     the same distribution is then evaluated by a F-test and a Tukey
#'     HSD Tukey if necessary 
#'     and by model selection AIC corrected for small samples
#'     (AICc). The  models compared are that the expected value of
#'     all values are the same, that the expected values
#'     of each sample are different, or that one group is different of
#'     the other two. The evidence weight and delta-AIC
#'     of the null hypothesis is then
#'     calculated from the AIC difference between the models. Also, the
#'     M-error (Gelman & Carlin 2014) is calculated dividing the estimated
#'     difference among group one and the other two by the true difference.
#' @references Gelman, A., & Carlin, J. (2014). Beyond power
#'     calculations assessing type s (sign) and type m (magnitude)
#'     errors. Perspectives on Psychological Science, 9(6), 641-651. 
wp.anova <- function(std.diff, sd1, N, mu1, delta2=FALSE){
    ## sd1 <- se1*sqrt(N*3)
    if(missing(mu1))
        delta <- std.diff*sd1*sqrt(3/N)
    if(missing(std.diff))
        delta <-  mu1
    x1 <- rnorm(N, sd=sd1)
    x2 <- rnorm(N, sd=sd1)
    x3 <- rnorm(N, mean = delta, sd=sd1)
    x <- c(x1,x2,x3)
    y <- factor(rep(c("am1","am2", "am3"), each=N))
    y3 <- factor(rep(c("am12","am12", "am3"), each=N))
    y2 <- factor(rep(c("am13","am2", "am13"), each=N))
    y1 <- factor(rep(c("am1","am23", "am23"), each=N))
    aov123 <- aov(x ~y)
    pF <- summary(aov123)[[1]][1,5]
    if(mu1 != 0) {
        if(pF < 0.05) {
            phoc123 <- TukeyHSD(aov123)
            aov.right <- phoc123$y[,4][1]>0.05 & all(phoc123$y[,4][-1]<0.05)
        }
        else
            aov.right <- FALSE
    }
    if(mu1 == 0) {
        aov.right <- pF >= 0.05
    }
    m0 <- lm( x ~ 1)
    m123 <- lm( x ~ y)
    m1 <-lm( x ~ y1)
    m2 <-lm( x ~ y2)
    m3 <-lm( x ~ y3)
    AICs <- c(AICc(m0, nobs=3*N),
              AICc(m1, nobs=3*N),
              AICc(m2, nobs=3*N),
              AICc(m3, nobs=3*N),
              AICc(m123, nobs=3*N))
    dAICs <- AICs - min(AICs)
    if(mu1 != 0){
        M.error <- coef(m3)[2]/mu1
        S.error <- M.error < 0
        if(delta2)
            aic.right <- dAICs[4]==0 & all(dAICs[-4]>2)
        ## agr1 <- (dAICs[4]==0 & all(dAICs[-4]>2))
        else
            aic.right <- dAICs[4]==0
    }
    if(mu1 == 0){
        M.error <-  NA
        S.error <- NA
        if(delta2)
            aic.right <- dAICs[1]==0 & all(dAICs[-1]>2)
        else                                   
            aic.right <- dAICs[1] < 2
        ##agr0 <-  (pF >= 0.05 & dAICs[1] < 2)
    }
    w <- exp(-dAICs/2)
    wi <- w/sum(w)
    results <- c(pF, wi[1], dAICs[1], M.error, S.error, aov.right, aic.right)
    names(results) <- c("p.value", "Akaike.weight.H0",
                        "deltaAIC.H0", "M.error", "S.error", "NHT.right", "AIC.right")
    return(results)
}

#' NHT and IT  comparison for regression test allowing collinearity
#'
#' This function simulates samples from a Gaussian distribution whose
#' mean is a linear function of a continuous predictor variable drawn from a standard
#' Gaussian distribution. A second variable that is used as another
#' predictor is also drawn, wich can have any correlation to the first
#' one, but does not affect the mean of the response variable.
#' The significance of the effect of each predictor is evaluated with
#' Wald marginal tests. The effects are also tested comparing
#' the AIC of linear models with each or both predictors. Four linear
#' regressions are then fit (with no effect -- null hypothesis --,
#' single effect of each predictor and additive effects of both
#' predictors). The
#' effects of each predictor are the checked by a null hypothesis test
#' (marginal Wald test from function \code{summary} and by model
#' selection using AIC corrected for small samples. For each test it
#' is considered a rightful conclusion if only the correct predictor
#' is spoted as a true effect, or if the null hypothesis is identified
#' in the absence of effects (beta = 0)
#'
#' @param std.beta real number, standardized effect (aka standardized
#'     slope coefficient) of the first predictor
#'     on the response variable (aka dependent variable). Ignored if
#'     bet is provided.
#' @param beta real number, effect (aka slope) of the first predictor
#'     on the response variable (aka dependent variable).
#' @param sd1 real positive, standard deviation of the effect, which
#'     is the standard deviation of the Gaussian of the response variable.
#' @param N integer positive, size of the sample.
#' @param rpears real number, Pearson correlation coefficient between
#'     the two variables to be tested as predictors (aka independent
#'     variables). Values of  \code{rpears} different from zero
#'     creates collinearity between the two predictor variables.
#' @param delta2 logical; should null hypothesis be rejected only if it has a
#'     deltaAIC larger than two? If FALSE the null
#'     hypothesis is rejected if it has larger deltaAIC than
#'     the alternative hypothesis.
#' @return a vector with of size 7: the p-value of a F-test for the
#'     null hypothesis of no effect; the Akaike evidence weight of the null
#'     hypothesis compared to the alternative models of effect of each
#'     predictor or both;
#'     the difference between Akaike Information Criteria of the null
#'     to the most plausible hypothesis; the ratio between the
#'     estimated effect of the first predictor and the true effect, if
#'     the effect when significant is in the opposite direction (S
#'     and M type errors, Gelman & Carlin, 2014), if the null hypothesis test
#'     reached to the right conclusion (yes=1, 
#'     no=0), if the model selection based on AICs reached to the right conclusion.
#' @details This function samples two predictor variables from a
#'     bivariate Gaussian with correlation term set by the user. The
#'     marginal distributions of both predictors are standard
#'     Gaussians (zero mean and unity standard error). A response
#'     variable is then simulated from a Gaussian distribution that
#'     has the mean value proportional to the first predictor. The
#'     proportionality factor is the effect, or slope , of the
#'     predictor and is set with argument \code{beta} or \{std.beta}. The standard
#'     deviation of the response is set by the argument \code{sd1}. If
#'     there is an effect of the predictor (beta!=0)
#'     it is also calculated the ratio between the estimated and true
#'     effect value (type M error, Gelman & Carlin 2014).
#' 
#' @references Gelman, A., & Carlin, J. (2014). Beyond power
#'     calculations assessing type s (sign) and type m (magnitude)
#'     errors. Perspectives on Psychological Science, 9(6), 641-651. 
wp.lm <- function(std.beta, sd1, N, beta, rpears, delta2=FALSE){
    if(missing(beta)) beta <- std.beta*sd1/sqrt(N)
    # sd1 <- se1*sqrt(N)
    sigma <- matrix(c(1, rpears, rpears, 1), nrow=2)
    x <- rmvnorm(n = N, sigma=sigma)
    y <- rnorm(length(x[,1]), mean = beta*x[,1], sd = sd1)
    m0 <- lm(y ~ 1)
    m1 <- lm( y ~ x[,1])
    m2 <- lm( y ~ x[,2])
    m12 <- lm( y ~ x[,1] + x[,2])
    F <- summary(m12)$fstatistic
    pF <- 1 - pf(F[1], F[2], F[3])
    p.vals <- summary(m12)$coefficients[2:3,4]
    AICs <- c(AICc(m0, nobs=N),
        AICc(m1, nobs=N),
        AICc(m2, nobs=N),
        AICc(m12, nobs=N))
    dAICs <- AICs - min(AICs)
    w <- exp(-dAICs/2)
    wi <- w/sum(w)
    if(beta==0){
        M.error <- NA
        S.error <- NA
        nht.right <- p.vals[1]>0.05 & p.vals[2]>0.05
        if(delta2)
            aic.right <- dAICs[1]==0 & all(dAICs[-1]>2)
        else
            aic.right <- dAICs[1]==0
    }
    else{
        M.error <-  coef(m12)[2]/beta
        S.error <- M.error < 0
        nht.right <- p.vals[1]<0.05 & p.vals[2]>0.05
        if(delta2)
            aic.right <- dAICs[2]==0 & all(dAICs[-2]>2)
        else
            aic.right <- dAICs[2]==0
        }
    results <- c(pF, wi[1], dAICs[1], M.error, S.error, nht.right, aic.right)
    names(results) <- c("p.value.F", "Akaike.weight.H0", "deltaAIC.H0",
                        "M.error", "S.error", "NHT.right", "AIC.right")
    return(results)
}


## A function to run replicates of any function above for
## a given combination of parameters standard effect, standard deviations
## and sample sizes and the simulation function to run
sim.averages <- function(effect, st.dev, sample.size, nrep,
                         function.name, delta2=FALSE, ...){
    f1 <- get(function.name)
    raw <- replicate(nrep, mapply(f1, effect, st.dev, sample.size,
                                  delta2=delta2, ...),
                     simplify=TRUE)
    results <- c(sum(raw[6,])/nrep, # Prop rightful conclusions NHT
                 sum(raw[7,])/nrep, # Prop rightfull conclusions IT
                 sum(raw[6,]!=raw[7,])/nrep,# Prop mismatches in conclusions IT x NHT
                 mean(raw[4,raw[6,]==1], na.rm=TRUE), # Mean M-error NHT
                 mean(raw[4,raw[7,]==1], na.rm=TRUE), # Mean M-error IT
                 mean(raw[5,raw[6,]==1], na.rm=TRUE), # Mean S-error NHT
                 mean(raw[5,raw[7,]==1], na.rm=TRUE), # Mean S-error IT
                 mean(raw[1,], na.rm=TRUE), # Mean p-value
                 mean(raw[2,], na.rm=TRUE) # Mean evidence weight for H0
                 ) 
    names(results) <- c("p.NHT.right", "p.AIC.right", "p.mismatch","mean.NHT.M",
                        "mean.AIC.M", "p.NHT.S",  "p.AIC.S", "mean.pvalue", "mean.wH0")
    return(results)
}

 
