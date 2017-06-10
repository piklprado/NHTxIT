## Na versao friendly review foi usado o critério do menor dletaAIC para identificar
## a hipotes com mais suporte no teste t. A figura desta versao mostra
## que a proporcao de conclusoes corretas é um pouco maior para NHT,
## no caso de amostras pequenas e efeito pequeno. No entanto os novos
## codigos do PI mostram o contrario. Um resultado similar ao anterior
## foi obtido com o criterio deltaAIC>2 para a identificar a hipotese
## correta

## Verificando com os meus codigos. A funcao para simular testes t tem
## um argumento delta2 que controla qual criterio é usado, veja no
## arquivo com as funcoes que criei e estou usando para fazer os
## ultimos graficos que te enviei
source("functions.R")
set.seed(42) ## para fazer os mesmos sorteior nod dois casos
results1 <- matrix(ncol=7, nrow=1000)
## Com o criterio do maior deltaAIC: IT tem mais poder
for(i in 1:1000)
    results1[i,] <- wp.ttest(sd1 = 1, N = 10, mu2 = 0.5, delta2 = FALSE)
## Proporcoes de conclusoes corretas para NHT e IT, respectivamente:
apply(results1[,6:7], 2, sum)/1000

## Com o criterio do maior deltaAIC: NHT tem mais poder
set.seed(42) ## para fazer os mesmos sorteio nos dois casos
results2 <-  matrix(ncol=7, nrow=1000)
for(i in 1:1000)
    results2[i,] <- wp.ttest(sd1 = 1, N = 10, mu2 = 0.5, delta2 = TRUE)
## Proporcoes de conclusoes corretas para NHT e IT, respectivamente:
apply(results2[,6:7], 2, sum)/1000

## Uma possibilidade:
## a funcao acima calcula o AICcs diretamente com a funcao interna do
## R chamada AICc. Nas versoes anteriores Leo calculou o AIC pela expressao a
## partir do residual sum of squares de cada modelo
## Entao vamos fazer passo a passo esses calculos e comparar:

## Duas amostras de tamanho 10, de normais com desvio-padrao 1 e
## medias que diferem em 0,5
x1 <- rnorm(n=10)
x2 <- rnorm(n=10, mean=0.5)

## AIC a partir dos residual sum of squares (Anderson, quadro na
## página 62): N*log(ML) + 2*K*(N/(n-k-1))
## Onde ML é  estimador de maxima verossimilhanca da variancia, que é a soma dos quadrados
## dividido pelo n de observacoes. Uma funcao pra fazer este estimador
mle.var <- function(x) sum((x - mean(x))^2 / length(x))
## Para H0 o ML é
X <- c(x1,x2) # apenas junta as duas amostras em um mesmo vetor pra facilitar
(ML.h0 <- sum( (X - mean(X))^2/ length(X) ))
## E o AICc calculado com ele de acordo com Anderson é
(AIC1.h0 <- length(X) * log(ML.h0) + 4 * (length(X)/(length(X) - 3)))
## Para HA o ML é a soma dos desvios quadráticos em relacao à medias
## dos grupos, divididos pelo n total de observacoes
(ML.ha <- (sum((x1 - mean(x1))^2) + sum((x2 - mean(x2))^2))/length(X))
## E o AIC será
(AIC1.hA <- length(X) * log(ML.ha) + 6 * (length(X)/(length(X) - 4)))

## Calculando os AICs a partir da funcao AICc
## Para isso é preciso construir dois modelos lineares:
## Um em que a media nao depende da identidade da amostra e outro em
## que a media depende
## Uma fator para indicar a identidade da amostra
y <-  factor( rep(letters[1:2], c(length(x1),length(x2))))
## Modelo que corresponde a H0
m1 <-  lm(X ~ 1)
## Modelo que corresponde a HA
m2 <-  lm(X ~ y)
## Os AICs
AIC2.h0 <- AICc(m1, nobs=length(X))
AIC2.hA <- AICc(m2, nobs = length(X))
## E agora comparando os delta-AICs obtidos pelos dois metodos
## O primeiro valor é o deltaAIC de NHT e o segundo de IT
## Pelo metodo do Anderson
c(AIC1.h0 , AIC1.hA) - min(c(AIC1.h0 , AIC1.hA))
## Pelo metodo de modelos
c(AIC2.h0 , AIC2.hA) - min(c(AIC2.h0, AIC2.hA))

## Delta-AIC identicos como esperado, entao nao eh isso. Mas vale
## conferir como foram caculados os AICs nos codigos do Leo.

## Por fim, vou repetir varias vezes as simualcoes acima num loop e
## com o calculo do Anderson, para verificar se dá o mesmo resultado
nrep <- 1000
results3 <- matrix(ncol=3, nrow=nrep, dimnames=list(NULL,
                                                    c("NHT","IT menor AIC", "IT dAIC>2")))
se.seed(42)
for(i in 1:nrep){
    ## sorteio
    x1 <- rnorm(n=10)
    x2 <- rnorm(n=10, mean=0.5)
    ## pvalor do teste t
    pvalue <- t.test(x1, x2, var.equal=TRUE)$p.value
    ## Calculos do AIC pelo metodo do Anderson
    ## Para H0 o ML é
    X <- c(x1,x2) # apenas junta as duas amostras em um mesmo vetor pra facilitar
    ML.h0 <- sum( (X - mean(X))^2/ length(X) )
    ## E o AICc calculado com ele de acordo com Anderson é
    AIC1.h0 <- length(X) * log(ML.h0) + 4 * (length(X)/(length(X) - 3))
    ## Para HA o ML é a soma dos desvios quadráticos em relacao à medias
    ## dos grupos, divididos pelo n total de observacoes
    ML.ha <- (sum((x1 - mean(x1))^2) + sum((x2 - mean(x2))^2))/length(X)
    ## E o AIC será
    AIC1.hA <- length(X) * log(ML.ha) + 6 * (length(X)/(length(X) -4))
    ## Delta AICs
    dAIC1 <- c(AIC1.h0 , AIC1.hA) - min(c(AIC1.h0 , AIC1.hA))
    ## Comparacao dos metodos
    nht.right <- pvalue < 0.05
    it.right1 <- dAIC1[1] > dAIC1[2]
    it.right2 <- dAIC1[1] > 2
    results3[i,] <- c( nht.right, it.right1, it.right2)
}

## Perguntando quantos casos em cada criterio
## Os resultados sao coerentes com os obtidos pela minha funcao:
## Se usamos o criterio de deltaAIC>2 NHT e IT acertam por volta de
## 20%, e NHT é superior. Usando o criterio de menor AIC IT passa a
## frente, com por volta de 30% de acerto
apply(results3, 2, sum)/nrow(results3)

