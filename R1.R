### Q1

library(glm2)
library(lmtest)

data(crabs)
names(crabs)
attach(crabs)

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
(LRTsat <- lrtest(satredfit, satfullfit))

## -------------------------------------------------------------------
# non-parametric

n <- length(crabs)
index <- c(1:n)

set.seed(2025)

B <- 10000
LRT.nb <- c(1:B)

for (i in 1:B) {
  # Sample Satellites and Width together 
  crabs.b <- crabs[sample(nrow(crabs), replace = TRUE), c("Satellites", "Width")]
  # sample GoodSpine separately
  crabs.b$GoodSpine <- sample(crabs$GoodSpine, size = nrow(crabs), replace = TRUE)
  
  satful.b <- glm(Satellites ~ Width + GoodSpine, family = poisson, data = crabs.b)
  satred.b <- glm(Satellites ~ Width, family = poisson, data = crabs.b)
  
  # Calculate LRT statistic
  LRT.nb[i] <- -2 * (logLik(satred.b) - logLik(satful.b))
}

png("pictures/histnonparam1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(LRT.nb, probability=T, nclass=50, xlim=c(0,40), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs,LRT.obs) ,c(0,5300),col=2)
x <- seq(0, 40, by = 0.01)  # Range for the chi-squared curve
lines(x, dchisq(x, df), col = "blue", lwd = 2, lty = 1)
dev.off()


df <- 1 
x <- seq(0, 10, by = 0.01)

# Plot the theoretical  chi-squared distribution
png("pictures/theoretical1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
plot(x, dchisq(x, df), type = "l", lwd = 2, col = "blue",
     xlab = "LRT Statistic",
     ylab = "Density",
     main = "Theoretical Chisq Distribution")
dev.off()

#(pval.nonp <- mean(LRT.nb >= LRT.obs)) # 0.759
pval.nonp <- (1 + sum(LRT.nb > LRT.obs)) / (B+1) # 0.7442256
pval.nonp
## -------------------------------------------------------------------
# parametric

Satellites <- crabs$Satellites
Width <- crabs$Width
GoodSpine <- crabs$GoodSpine

(hat.mu<-mean(Satellites))

set.seed(2025)
B <- 10000  # Number of bootstrap samples
LRT.pb <- c(1:B)

for (i in 1:B) {
  Satellites.b <- rpois(n = nrow(crabs), lambda = hat.mu)
  satful.b <- glm(Satellites.b ~ Width + GoodSpine, family=poisson)
  satred.b <- glm(Satellites.b ~ Width, family=poisson)
  
  LRT.pb[i] <- -2*(logLik(satred.b) - logLik(satful.b))
  #LRT.pb[i] <- lrtest(satredfit, satfullfit)
}

png("pictures/histparam1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(LRT.pb, nclass=50, xlim=c(-1,20), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs, LRT.obs) ,c(0,5300),col=2)
dev.off()

(pval.p <- (1 + sum(LRT.pb > LRT.obs)) / (B+1)) # 0.5624376


png("pictures/histcombined1.2.png", width = 18, 
    height = 10, units = "cm", res = 300)
par(mfrow=c(1,2))
hist(LRT.nb, probability=T, nclass=50, xlim=c(0,40), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs,LRT.obs) ,c(0,5300),col=2)
lines(x, dchisq(x, df), col = "blue", lwd = 2, lty = 1)

hist(LRT.pb, probability=T, nclass=50, xlim=c(-1,20), xlab="LRT.b", main = "Histogram of LRT.b")
lines(c(LRT.obs, LRT.obs) ,c(0,5300),col=2)
lines(x, dchisq(x, df), col = "blue", lwd = 2, lty = 1)
dev.off()
par(mfrow=c(1,1))

## -------------------------------------------------------------------
# permutation
set.seed(2025)
B <- 10000  # Number of bootstrap samples
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
lines(c(LRT.obs,LRT.obs) ,c(0,5300),col=2)
dev.off()

(pval.p <- (1 + sum(LRT.permb > LRT.obs)) / (B+1))  # 0.7472527

