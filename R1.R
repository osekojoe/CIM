### Q1

library(glm2)
library(lmtest)

data(crabs)
names(crabs)

plot(crabs$Width,crabs$Satellites)
head(crabs)

crabs <- na.omit(crabs)

########## 1.1
satfullfit <- glm(Satellites ~ Width + GoodSpine,
                  family=poisson(link="log"), data=crabs)
summary(satfullfit)


# LRT 
satredfit <- glm(Satellites ~ Width,
                 family=poisson(link="log"), data=crabs)

LRT.obs <- -2*(logLik(satredfit) - logLik(satfullfit))
LRT.obs

pval <- pchisq(LRT, df = 1, lower.tail = FALSE)
pval
(LRTsat <- lrtest(satfullfit, satredfit))

## -------------------------------------------------------------------
# parametric

n <- length(crabs)
index <- c(1:n)

set.seed(2025)

B <- 1000
LRT.nb <- c(1:B)

for (i in 1:B) {
  crabs.b <- crabs[sample(nrow(crabs), replace = TRUE), ]
  
  satful.b <- glm(Satellites ~ Width + GoodSpine, family=poisson, data=crabs.b)
  satred.b <- glm(Satellites ~ Width, family=poisson, data=crabs.b)
  
  LRT.nb[i] <- -2*(logLik(satred.b) - logLik(satful.b))
}

hist(LRT.nb, nclass=50, xlim=c(0,30))
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)

(pval.nonp <- mean(LRT.b >= LRT.obs))

## -------------------------------------------------------------------
# non-parametric

Satellites <- crabs$Satellites
Width <- crabs$Width
GoodSpine <- crabs$GoodSpine

set.seed(2025)
B <- 1000  # Number of bootstrap samples
LRT.pb <- c(1:B)

for (i in 1:B) {
  Satellites.b <- rpois(n = nrow(crabs), lambda = exp(predict(satredfit)))
  satful.b <- glm(Satellites.b ~ Width + GoodSpine, family=poisson)
  satred.b <- glm(Satellites.b ~ Width, family=poisson)
  
  LRT.pb[i] <- -2*(logLik(satred.b) - logLik(satful.b))
  #LRT.pb[i] <- lrtest(satredfit, satfullfit)
}

hist(LRT.pb, nclass=50, xlim=c(-1,13))
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)

(pval.p <- mean(LRT.pb >= LRT.obs))

## -------------------------------------------------------------------
# permutation
set.seed(2025)
B <- 1000  # Number of bootstrap samples
LRT.permb <- c(1:B)

for (i in 1:B) {
  permuted_data <- crabs
  permuted_data$Satellites <- sample(crabs$Satellites, replace = FALSE)
  
  satful.b <- glm(Satellites ~ Width + GoodSpine, family=poisson, data=permuted_data)
  satred.b <- glm(Satellites ~ Width, family=poisson, data=permuted_data)
  
  LRT.permb[i] <- -2*(logLik(satred.b) - logLik(satful.b))
}

hist(LRT.permb, nclass=50, xlim=c(0,40))
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)

(pval.p <- mean(LRT.permb >= LRT.obs))

