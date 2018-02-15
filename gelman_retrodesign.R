## Not used ##
## Function to estimate M and S errors for Gaussian
## Gelman & Carlin Perspectives on Psychological Science 2014, Vol. 9(6) 641 –651
retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
    z <- qt(1-alpha/2, df)
    p.hi <- 1 - pt(z-A/s, df)
    p.lo <- pt(-z-A/s, df)
    power <- p.hi + p.lo
    typeS <- p.lo/power
    estimate <- A + s*rt(n.sims,df)
    significant <- abs(estimate) > s*z
    exaggeration <- mean(abs(estimate)[significant])/A
    return(c(power=power, typeS=typeS, exaggeration=exaggeration)) # changed from a list to a named vector
}

## Modified to include t-value equivalent to a log0likelhood=2
## Function to estimate M and S errors for Gaussian
## Gelman & Carlin Perspectives on Psychological Science 2014, Vol. 9(6) 641 –651
retrodesign2 <- function(A, s = 1 , alpha=.05, L = 2, df=1e6, n.sims=10000){
    n <- df+1
    z1 <- qt(1-alpha/2, df)
    z2 <- sqrt(df*(exp(2*L/n)-1))
    p.hi1 <- 1 - pt(z1-A/s, df)
    p.lo1 <- pt(-z1-A/s, df)
    power1 <- p.hi1 + p.lo1
    typeS1 <- p.lo1/power1
    estimate1 <- A + s*rt(n.sims,df)
    significant1 <- abs(estimate1) > s*z1
    exaggeration1 <- mean(abs(estimate1)[significant1])/A
    p.hi2 <- 1 - pt(z2-A/s, df)
    p.lo2 <- pt(-z2-A/s, df)
    power2 <- p.hi2 + p.lo2
    typeS2 <- p.lo2/power2
    estimate2 <- A + s*rt(n.sims,df)
    significant2 <- abs(estimate2) > s*z2
    exaggeration2 <- mean(abs(estimate2)[significant2])/A
    return(c(power=power1, Lpower=power2, typeS=typeS1, LtypeS=typeS2,
             exaggeration=exaggeration1, Lexaggeration=exaggeration2)) # changed from a list to a named vector
}

## Same using simulated data drawn from t-distribution
retrodesign3 <- function(A, df= 1e6 , s=1, alpha = 0.05, L = 2, n.sims=1e4){
    n <- df + 1
    z1 <- qt(1 - alpha/2, df)
    z2 <- sqrt(df*(exp(2*L/n)-1))
    estimate <- rt(n.sims, df = df, ncp = A/s)
    p.hi1 <- 1 - pt(z1, df = df, ncp = A/s)
    p.lo1 <- pt(-z1, df= df, ncp=A/s)
    power1 <- p.hi1 + p.lo1
    typeS1 <- p.lo1/power1
    #p.hi2 <- 1 - pt(z2, df = df, ncp = A/s)
    #p.lo2 <- pt(-z2, df= df, ncp=A/s)
    #power2 <- p.hi2 + p.lo2
    #typeS2 <- p.lo2/power2
    significant1 <- abs(estimate) > z1
    significant2 <- -2*dt(estimate, df=df, ncp = estimate, log=TRUE)+4 <  -2*dt(estimate, df=df, log=TRUE)+2
    exaggeration1 <- mean(abs(estimate*s)[significant1])/A
    exaggeration2 <- mean(abs(estimate*s)[significant2])/A
    power2 <- mean(significant2)
    typeS2 <- mean(estimate[significant2]/A < 0)
    return(c(power=power1, Lpower=power2, typeS=typeS1, LtypeS=typeS2,
             exaggeration=exaggeration1, Lexaggeration=exaggeration2))
    }
