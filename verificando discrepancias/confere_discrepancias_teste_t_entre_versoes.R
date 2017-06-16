## Na versao friendly review foi usado o critério do menor dletaAIC para identificar
## a hipotes com mais suporte no teste t. A figura desta versao mostra
## que a proporcao de conclusoes corretas é um pouco maior para NHT,
## no caso de amostras pequenas e efeito pequeno. No entanto os novos
## codigos do PI mostram o contrario. Um resultado similar ao anterior
## foi obtido com o criterio deltaAIC>2 para a identificar a hipotese
## correta

################################################################################
## PARTE 1: usando funcoes criadas pelo PI
################################################################################
## Verificando com os meus codigos. A funcao para simular testes t tem
## um argumento delta2 que controla qual criterio é usado, veja no
## arquivo com as funcoes que criei e estou usando para fazer os
## ultimos graficos que te enviei
source("../functions.R")
set.seed(42) ## para fazer os mesmos sorteior nod dois casos
results1 <- matrix(ncol=7, nrow=1000)
## Com o criterio do maior deltaAIC: IT tem mais poder
for(i in 1:1000)
    results1[i,] <- wp.ttest(sd1 = 1, N = 10, mu2 = 1, delta2 = FALSE)
## Proporcoes de conclusoes corretas para NHT e IT, respectivamente:
apply(results1[,6:7], 2, sum)/1000

## Com o criterio do maior deltaAIC: NHT tem mais poder
set.seed(42) ## para fazer os mesmos sorteio nos dois casos
results2 <-  matrix(ncol=7, nrow=1000)
for(i in 1:1000)
    results2[i,] <- wp.ttest(sd1 = 1, N = 10, mu2 = 1, delta2 = TRUE)
## Proporcoes de conclusoes corretas para NHT e IT, respectivamente:
apply(results2[,6:7], 2, sum)/1000

################################################################################
## Parte 2: conferindo se o calculo de deltaAIC das funcoes do PI
## batem com os calculos de deltaAIC pela expressao do Anderson,
## eu foi a usada pelo Leo nas versoes anteriores.
## A conslusao eh que as duas batem, como esperado.
################################################################################
## Uma possibilidade:
## a funcao acima calcula o AICcs diretamente com a funcao interna do
## R chamada AICc. Nas versoes anteriores Leo calculou o AIC pela expressao a
## partir do residual sum of squares de cada modelo
## Entao vamos fazer passo a passo esses calculos e comparar:

## Duas amostras de tamanho 10, de normais com desvio-padrao 1 e
## medias que diferem em 1
x1 <- rnorm(n=10)
x2 <- rnorm(n=10, mean=1)

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
## CONCLUSAO: ##
## Delta-AIC identicos como esperado, entao nao eh isso. Mas vale
## conferir como foram caculados os AICs nos codigos do Leo.

################################################################################
## PARTE 3: repete a simulacao tomando amostras de uma normal no R, calculando o AICc
## pelo criterio do Anderson
## CONCLUSAo: os resultados sao os mesmos dos obtidos com as funcoes do PI
################################################################################
## Um loop que repete os passos da Parte 2:
## N de repeticoes
nrep <- 1000
## Uma matriz para guardar os resultados
results3 <- matrix(ncol=3, nrow=nrep, dimnames=list(NULL,
                                                    c("NHT","IT menor AIC", "IT dAIC>2")))
## As simulacoes
set.seed(42)
for(i in 1:nrep){
    ## sorteio
    x1 <- rnorm(n=10)
    x2 <- rnorm(n=10, mean=1)
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
## 50%, e NHT é superior. Usando o criterio de menor AIC IT passa a
## frente, com por volta de 70% de acerto.
apply(results3, 2, sum)/nrow(results3)
## Esses valores sao os mesmos obtidos usando a minha função wp.ttest (PARTE 1)
## Como reiniciei o gerador de numeros aleatorios, todas as simulacoes estao usando as mesmas amostras,
## Por isso os valores sao exatamente iguais.
c(apply(results1[,6:7], 2, sum)/1000, apply(results2[,6:7], 2, sum)[-1]/1000)

################################################################################
## PARTE 4: repete as simulacoes da parte 3, mas com uma populacao gerada
## pelo Crystal Ball
## CONCLUSAO: os resultados continuam os mesmos.
################################################################################
## Arquivo txt com as populacoes geradas pelo Crystall Ball, enviado pelo Leo 15/06/17
## Nome original do arquivo era "N,1,1.txt"
CB <- read.table("crystal_ball_N_1_1.txt", header=TRUE)
## populacao de 10 mil observacoes gerada com diferenca entre medias=1 e desvio-padrao=1
summary(CB)
## Medias e sd populacionais ligeiramente diferentes do teorico, como era de se esperar
## Mas nao acredito que isso faca diferenca
apply(CB, 2, sd)
apply(CB, 2, mean)

## Repetindo o mesmo teste, com criterio do Anderson, para amostras de tamanho dez tomadas dessas
## Populacoes do Crystall Ball
## Matriz para guardar os resultados
results4 <- matrix(ncol=3, nrow=nrep, dimnames=list(NULL,
                                                    c("NHT","IT menor AIC", "IT dAIC>2")))
## Um loop identico ao anterior, exceto que toma amostras das populacoes
## ao invés de sortear no R 
set.seed(42) ## seta de novo o gerador de n aleatorios para reprodutibilidade
for(i in 1:nrep){
    ## sorteio: esta é a unica parte do codigo da seção anterior que muda,
    ## agora amostra com reposicao das populacoes criadas pelo Crystal Ball
    x1 <- sample(CB$x1, 10, replace=TRUE)
    x2 <- sample(CB$x2, 10, replace=TRUE)
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
    results4[i,] <- c( nht.right, it.right1, it.right2)
}
## Resulta no mesmo padrão: 
apply(results4, 2, sum)/nrow(results4)
## Comparando com as simulacoes da Parte 3:
apply(results3, 2, sum)/nrow(results3)
