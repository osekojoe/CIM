### Q1

library(glm2)
library(lmtest)

data(crabs)
names(crabs)

plot(crabs$Width,crabs$Satellites)
head(crabs)

crabs <- na.omit(crabs)
crabs$GoodSpine <- as.factor(crabs$GoodSpine)


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
# non-parametric

n <- length(crabs)
index <- c(1:n)

set.seed(2025)

B <- 1000
LRT.nb <- c(1:B)

for (i in 1:B) {
  #crabs.b <- crabs[sample(nrow(crabs), replace = TRUE), ]
  #sample Satellites and Width together
  indices <- sample(1:n, replace = TRUE)
  Satellites.b <- crabs$Satellites[indices]
  Width.b <- crabs$Width[indices]
  #sample Goodspine seperately
  GoodSpine.b <- crabs$GoodSpine[sample(1:n, replace = TRUE)]
  
  crabs.b <- data.frame(Satellites = Satellites.b,
                        Width = Width.b,
                        GoodSpine = GoodSpine.b)
  
  satful.b <- glm(Satellites ~ Width + GoodSpine, family=poisson, data=crabs.b)
  satred.b <- glm(Satellites ~ Width, family=poisson, data=crabs.b)
  
  
  LRT.nb[i] <- -2*(logLik(satred.b) - logLik(satful.b))
}

png("pictures/histnonparam1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(LRT.nb, nclass=50, xlim=c(0,30), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)
dev.off()

(pval.nonp <- mean(LRT.nb >= LRT.obs)) # 0.759

## -------------------------------------------------------------------
# parametric

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

png("pictures/histparam1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(LRT.pb, nclass=50, xlim=c(-1,13), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)
dev.off()

(pval.p <- mean(LRT.pb >= LRT.obs)) # 0.556

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

png("pictures/histperm1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(LRT.permb, nclass=50, xlim=c(0,40), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs,LRT.obs) ,c(0,500),col=2)
dev.off()

(pval.p <- mean(LRT.permb >= LRT.obs)) # 0.747

