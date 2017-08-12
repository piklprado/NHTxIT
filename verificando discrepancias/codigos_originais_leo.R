
## Teste-t (NHT)
R=10000
p=numeric(R)
for(i in 1:R){
    amostra1=sample(x1,size=N, replace=T)
    amostra2=sample(x2,size=N,replace=T)
    p[i]=t.test(amostra1,amostra2,var.equal=TRUE)$p.value
}
## Teste-t(IT)
R=10000
media=numeric(R)
RSSn=numeric(R)
RSSa=numeric(R)
aiccn=numeric(R)
aicca=numeric(R)
Delta_aicc=numeric(R)
for(i in 1:R){
    amostra1=sample(x1,size=N,replace=T)
    amostra2=sample(x2,size=N,replace=T)
    media[i]=((mean(amostra1)+mean(amostra2))/2)
    RSSn[i]=sum((amostra1-media[i])^2)+sum((amostra2-media[i])^2)
    RSSa[i]=sum((amostra1-(mean(amostra1)))^2)+sum((amostra2-(mean(amostra2)))^2)
    aiccn[i]=N*log(RSSn[i]/N)+4+(12/(N-3))
    aicca[i]=N*log(RSSa[i]/N)+6+(24/(N-4))
    Delta_aicc[i]=aiccn[i]-aicca[i]
}

## ANOVA (IT)

R=10000
media1=numeric(R)
media2=numeric(R)
media3=numeric(R)
media1.2=numeric(R)
media1.3=numeric(R)
media2.3=numeric(R)
media1.2.3=numeric(R)
RSSn=numeric(R)
RSS12.3=numeric(R)
RSS1.23=numeric(R)
RSS13.2=numeric(R)
RSS1.2.3=numeric(R)
aiccn=numeric(R)
aicc12.3=numeric(R)
aicc1.23=numeric(R)
aicc13.2=numeric(R)
aicc1.2.3=numeric(R)
for(i in 1:R){
    amostra1=sample(x1,size=N,replace=T)
    amostra2=sample(x2,size=N,replace=T)
    amostra3=sample(x3,size=N,replace=T)
    media1[i]=mean(amostra1)
    media2[i]=mean(amostra2)
    media3[i]=mean(amostra3)
    media1.2[i]=(mean(amostra1)+mean(amostra2))/2
    media1.3[i]=(mean(amostra1)+mean(amostra3))/2
    media2.3[i]=(mean(amostra2)+mean(amostra3))/2
    media1.2.3[i]=(mean(amostra1)+mean(amostra2)+mean(amostra3))/3
    RSSn[i]=sum((amostra1-media1.2.3[i])^2,(amostra2-media1.2.3[i])^2,(amostra3-media1.2.3[i])^2)
    RSS12.3[i]=sum((amostra1-media1.2[i])^2,(amostra2-media1.2[i])^2,(amostra3-media3[i])^2)
    RSS1.23[i]=sum((amostra1-media1[i])^2,(amostra2-media2.3[i])^2,(amostra3-media2.3[i])^2)
    RSS13.2[i]=sum((amostra1-media1.3[i])^2,(amostra2-media2[i])^2,(amostra3-media1.3[i])^2)
    RSS1.2.3[i]=sum((amostra1-media1[i])^2,(amostra2-media2[i])^2,(amostra3-media3[i])^2)
    aiccn[i]=N*log(RSSn[i]/N)+4+(12/(N-3))
    aicc12.3[i]=N*log(RSS12.3[i]/N)+6+(24/(N-4))
    aicc1.23[i]=N*log(RSS1.23[i]/N)+6+(24/(N-4))
    aicc13.2[i]=N*log(RSS13.2[i]/N)+6+(24/(N-4))
    aicc1.2.3[i]=N*log(RSS1.2.3[i]/N)+8+(40/(N-5))
}
## ANOVA(NHT)
R=10000
SSb=numeric(R)
SSw=numeric(R)
dfw=numeric(R)
MSb=numeric(R)
MSw=numeric(R)
F=numeric(R)
for(i in 1:R){
    amostra1=sample(x1,size=N,replace=T)
    amostra2=sample(x2,size=N,replace=T)
    amostra3=sample(x3,size=N,replace=T)
    SSb[i]=(((sum(amostra1)^2)/N)+((sum(amostra2)^2)/N)+((sum(amostra3)^2)/N))-((sum(amostra1,amostra2,amostra3)^2)/(3*N))
    SSw[i]=sum((amostra1-mean(amostra1))^2,(amostra2-mean(amostra2))^2,(amostra3-mean(amostra3))^2)
    dfw[i]=(3*N)-3
    MSb[i]=SSb[i]/2
    MSw[i]=SSw[i]/dfw[i]
    F[i]=MSb[i]/MSw[i]
}
## Tukey
R=10000
media1=numeric(R)
media2=numeric(R)
media3=numeric(R)
RSS=numeric(R)
MS=numeric(R)
Q12=numeric(R)
Q13=numeric(R)
Q23=numeric(R)
df=numeric(R)
for(i in 1:R){
    amostra1=sample(x1,size=N,replace=T)
    amostra2=sample(x2,size=N,replace=T)
    amostra3=sample(x3,size=N,replace=T)
    media1[i]=mean(amostra1)
    media2[i]=mean(amostra2)
    media3[i]=mean(amostra3)
    RSS[i]=sum((amostra1-media1[i])^2,(amostra2-media2[i])^2,(amostra3-media3[i])^2)
    df[i]=(3*N)-3
    MS[i]=RSS[i]/df[i]
    Q12[i]=(media1[i]-media2[i])/sqrt(MS[i]/N)
    Q13[i]=(media1[i]-media3[i])/sqrt(MS[i]/N)
    Q23[i]=(media2[i]-media3[i])/sqrt(MS[i]/N)
}


## Correlação (NHT)
R=10000
p=numeric(R)
for(i in 1:R){
    n=sample(1:10000,size=N,replace=F)
    amostra1=x[n]
    amostra2=y[n]
    p[i]=cor.test(amostra1,amostra2)$p.value
}
## Correlação (IT)
R=10000
delta_aicc=numeric(R)
for(i in 1:R){
    n=sample(1:10000,size=N,replace=F)
    amostra1=x[n]
    amostra2=y[n]
    modn<-lm(amostra1~1)
    moda<-lm(amostra1~amostra2)
    delta_aicc[i]<-AICc(modn)-AICc(moda)
}
##GLM (NHT)
R=10000
N=N
py_x1=numeric(R)
py_x2=numeric(R)
for(i in 1:R){
    n=sample(1:10000,size=N,replace=F)
    amostray=y[n]
    amostrax1=x1[n]
    amostrax2=x2[n]
    py_x1[i]=summary(glm(amostray~amostrax1+amostrax2))$coef[2,"Pr(>|t|)"]
    py_x2[i]=summary(glm(amostray~amostrax1+amostrax2))$coef[3,"Pr(>|t|)"]
}
R=10000
N=N
AICc_modn=numeric(R)
AICc_moda1=numeric(R)
AICc_moda2=numeric(R)
AICc_moda3=numeric(R)
for(i in 1:R){
    n=sample(1:10000,size=N,replace=F)
    amostray=y[n]
    amostrax1=x1[n]
    amostrax2=x2[n]
    modn<-glm(amostray~1)
    moda1<-glm(amostray~amostrax1)
    moda2<-glm(amostray~amostrax2)
    moda3<-glm(amostray~amostrax1+amostrax2)
    AICc_modn[i]=AICc(modn)
    AICc_moda1[i]=AICc(moda1)
    AICc_moda2[i]=AICc(moda2)
    AICc_moda3[i]=AICc(moda3)
}
