## Plot functions ###

## Proportion of rightfull conclusions x effect size
p1b <- function(data, legend = TRUE, pos.leg="bottomright",
                AIC1 = TRUE, AIC2 = TRUE, H0= FALSE,
                colours=c("black","red","blue"),
                ylab="P rightfull conclusions",
                xlab="Effect size", ylim, ...){
    dots <- list(...)
    if(H0){
        Y <- data$SD/sqrt(data$N)
        if(xlab!="") xlab="Standard error"
    }
    else{
        Y <- data$D
        if(xlab!="") xlab="Effect size"
    }
    if(!AIC2&missing(ylim))
        ylim <- with(data, range(c(p.NHT.right,p.AIC.right)))
    else if(!AIC1&missing(ylim))
        ylim <- with(data, range(c(p.NHT.right,p.AIC.right.2)))
    else if(missing(ylim))
        ylim <- with(data,
                     range(c(p.NHT.right,p.AIC.right, p.AIC.right.2)))
    plot(p.NHT.right ~ Y, data=data,
         ylim= ylim,
         xlab=xlab, ylab=ylab, col=colours[1], ...)
    if(AIC1)
        points(p.AIC.right ~ Y, data=data, col=colours[2], ...)
    if(AIC2)
        points(p.AIC.right.2 ~ Y, data=data, col=colours[3], ...)
    if(legend)
        legend(pos.leg,
               c("NHT", "IT crit. 1", "IT crit. 2"),
       pch=19, col=colours[c(1,3,2)], bty="n")
}

## Proportion of mismatching conclusions x Effect size
p2b <- function(data, legend = TRUE,
                pos.leg="bottomright", AIC2 = TRUE, ...){
    plot(p.mismatch ~ D, data=data, cex=0.25,
         ylim=range(c(p.NHT.right,p.AIC.right, p.AIC.right.2)),
         xlab="Effect size", ylab="P mismatching conclusions",
         col = "blue", ...)
    if(AIC2)
        points(p.mismatch.2 ~ D, data=data, cex=0.25, col="red")
    if(legend)
        legend(pos.leg,
               c(
                 expression(paste("NHT x ", Delta,"AIC = 0")),
                 expression(paste("NHT x ",Delta,"AIC = 0, all other ",
                                  Delta, "AIC > 2"))),
       pch=19, col=c("blue", "red"), bty="n")
}

## Mean M-error x effect size
p3b <- function(data, legend = TRUE, pos.leg="bottomright",
                AIC1 = TRUE, AIC2 = TRUE,
                colours=c("black","red","blue"),
                ylab="Mean M-error",
                xlab="Effect size", ylim,
                pred.line=FALSE, x.pred, y.pred, ...){
    if(!AIC2&missing(ylim))
        ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M)))
    else if(!AIC1&missing(ylim))
        ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M.2)))
    else if(missing(ylim))
        ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M,
                                   mean.AIC.M.2))) 
    if(pred.line){
        plot(x.pred, y.pred,
         ylim=ylim,
         xlab=xlab, ylab=ylab, col= "red", type="l", lwd = 4, ...)
        points(mean.NHT.M ~ D, data=data, col= colours[1], ...)
    }
    else
        plot(mean.NHT.M ~ D, data=data,
             ylim=ylim,
             xlab=xlab, ylab=ylab,
             col= colours[1], ...)
    if(AIC1) 
        points(mean.AIC.M ~ D, data=data, col=colours[2], ...)
    if(AIC2)
        points(mean.AIC.M.2 ~ D, data=data, col=colours[3], ...)
    if(legend)
        legend(pos.leg,
               c("NHT", "IT crit. 1", "IT crit. 2"),
       pch=19, col=colours[c(1,3,2)], bty="n")
}

## Proportion of S-error x effect size
p4b <- function(data, legend = TRUE, pos.leg="bottomright",
                AIC1 = TRUE, AIC2 = TRUE,
                colours=c("black","red","blue"),
                ylab="Mean M-error",
                xlab="Effect size", ylim, ...){
    if(!AIC2&missing(ylim))
        ylim <- with(data, range(c(p.NHT.S,p.AIC.S)))
    else if(!AIC1&missing(ylim))
        ylim <- with(data, range(c(p.NHT.S,p.AIC.S.2)))
    else if(missing(ylim))
        ylim <- with(data, range(c(p.NHT.S,p.AIC.S, p.AIC.S.2)))
    plot(p.NHT.S ~ D, data=data,
         ylim=ylim,
         xlab=xlab, ylab=ylab, col= colours[1], ...)
    if(AIC1) 
        points(p.AIC.S ~ D, data=data, col=colours[2], ...)
    if(AIC2)
        points(p.AIC.S.2 ~ D, data=data, col=colours[3], ...)
    if(legend)
        legend(pos.leg,
               c("NHT", "IT crit. 1", "IT crit. 2"),
       pch=19, col=colours[c(1,3,2)], bty="n")
}

## Mean M-error x proportion rightfull conclusions
p5b <- function(data, legend = TRUE, pos.leg="bottomright",
                AIC1 = TRUE, AIC2 = TRUE,
                colours=c("black","red","blue"),
                ylab="Mean M-error",
                xlab="Effect size", ylim, xlim, ...){
    if(!AIC2){
        if(missing(ylim))
            ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M)))
        if(missing(xlim))
            xlim <- with(data, range(c(p.NHT.right,p.AIC.right)))
    }
    else if(!AIC1){
        if(missing(ylim))
            ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M.2)))
        if(missing(xlim))
        xlim <- with(data, range(c(p.NHT.right,p.AIC.right.2)))
        }
    else if(missing(ylim))
        ylim <- with(data, range(c(mean.NHT.M,mean.AIC.M,
                                   mean.AIC.M.2))) 
    else if(missing(xlim))
        xlim <- with(data, range(c(p.NHT.right,p.AIC.right,
                                   p.AIC.right.2))) 
    plot(mean.NHT.M ~ p.NHT.right, data=data,
         ylim=ylim,
         xlab=xlab, ylab=ylab, col= colours[1], ...)
    if(AIC1) 
        points(mean.AIC.M ~ p.AIC.right, data=data, col=colours[2],
               ...) 
    if(AIC2)
        points(mean.AIC.M.2 ~ p.AIC.right.2, data=data,
               col=colours[3], ...) 
    if(legend)
        legend(pos.leg,
               c("NHT", "IT crit. 1", "IT crit. 2"),
       pch=19, col=colours[c(1,3,2)], bty="n")
}

## Proportion of S-error x proportion rightfull conclusions
p6b <- function(data, legend = TRUE, pos.leg="bottomright",
                AIC1 = TRUE, AIC2 = TRUE,
                colours=c("black","blue","red"),
                ylab="Mean M-error",
                xlab="Effect size", ylim, xlim, ...){
    if(!AIC2){
        if(missing(ylim))
            ylim <- with(data, range(c(p.NHT.S,p.AIC.S)))
        if(missing(xlim))
            xlim <- with(data, range(c(p.NHT.right,p.AIC.right)))
    }
    else if(!AIC1){
        if(missing(ylim))
            ylim <- with(data, range(c(p.NHT.S,p.AIC.S.2)))
        if(missing(xlim))
        xlim <- with(data, range(c(p.NHT.right,p.AIC.right.2)))
        }
    else {
            if(missing(ylim))
                ylim <- with(data, range(c(p.NHT.S,p.AIC.S,
                                           p.AIC.S.2))) 
            if(missing(xlim))
                xlim <- with(data, range(c(p.NHT.right,p.AIC.right,
                                           p.AIC.right.2))) 
            }
    plot(p.NHT.S ~ p.NHT.right, data=data,
         ylim=ylim, xlim=xlim,
         xlab=xlab, ylab=ylab, col= colours[1], ...)
    if(AIC1) 
        points(p.AIC.S ~ p.AIC.right, data=data, col=colours[2], ...)
    if(AIC2)
        points(p.AIC.S.2 ~ p.AIC.right.2, data=data, col=colours[3],
               ...) 
    if(legend)
        legend(pos.leg,
               c("NHT", "IT crit. 1", "IT crit. 2"),
       pch=19, col=colours[c(1,3,2)], bty="n")
}

