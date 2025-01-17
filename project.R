
## Project 2
##------------------------------------------------------------------------------
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

(LRTsat <- lrtest(satfullfit, satredfit))
# no evidence to reject the null hypothesis. Adding GoodSpine doesnt improve
# the model fit.

##------------------------------------------------------------------------------
########## 1.2 
### parametric bootstrap
set.seed(3542)
(hat.mu <- mean(crabs$Satellites))

n <- length(crabs[,1])
B <- 10000
beta0.b <- beta1.b <- beta2.b <- c(1:B)
for (i in 1:B)
{
  Satellites.b <- rpois(n, hat.mu)
  Width.b <- Width
  GoodSpine.b <- GoodSpine
  fit.glm.b <- glm(Satellites.b ~ Width.b + GoodSpine.b, family=poisson)
  beta0.b[i]<-summary(fit.glm.b)$coeff[1,1]
  beta1.b[i]<-summary(fit.glm.b)$coeff[2,1]
  beta2.b[i] <- summary(fit.glm.b)$coeff[3,1]
}

(beta2.obs<-summary(satfullfit)$coeff[3,1])

par(mfrow=c(1,1))
hist(beta2.b, nclass=50)
lines(c(beta2.obs, beta2.obs),c(0,1000), col=2)

# monte carlo pvalue
(1+sum(abs(beta2.b)>abs(beta2.obs)))/(B+1)

##----------------------------

### non-parametric bootstrap
Satellites <- crabs$Satellites
Width <- crabs$Width
GoodSpine <- crabs$GoodSpine

set.seed(3542)
n <- length(crabs[,1])
B <- 10000
beta0.b <- beta1.b <- beta2.b <-c(1:B)
for (i in 1:B) {
  Satellites.b <- sample(Satellites, size = n, replace = TRUE)
  Width.b <- sample(Width, size = n, replace = TRUE)
  GoodSpine.b <- sample(GoodSpine, size = n, replace = TRUE)
  fit.glm.b <- glm(Satellites.b ~ Width.b + GoodSpine.b, family=poisson(link="log"))
  beta0.b[i] <- summary(fit.glm.b)$coeff[1,1]
  beta1.b[i] <- summary(fit.glm.b)$coeff[2,1]
  beta2.b[i] <- summary(fit.glm.b)$coeff[3,1]
}

# observed slope
(beta2.obs<-summary(satfullfit)$coeff[3,1])

# The distribution of the test statistic under the null
par(mfrow=c(1,1))
hist(beta2.b,nclass=50)
lines(c(beta2.obs, beta2.obs), c(0,1000), col=2)

# The Monte carlo p value.. 0.7546245, greater than 0.05
(1 + sum(abs(beta2.b) > abs(beta2.obs)))/(B+1)

##------------------------------------------------------------------------------
### permutation test

Satellites <- crabs$Satellites
Width <- crabs$Width
GoodSpine <- crabs$GoodSpine

n <- length(crabs[,1])
B <- 10000
beta0.b <- beta1.b <- beta2.b <-c(1:B)
for (i in 1:B) {
  Satellites.b <- sample(Satellites, size = n, replace = FALSE)
  Width.b <- sample(Width, size = n, replace = FALSE)
  GoodSpine.b <- sample(GoodSpine, size = n, replace = FALSE)
  fit.glm.b <- glm(Satellites.b ~ Width.b + GoodSpine.b, family=poisson(link="log"))
  beta0.b[i] <- summary(fit.glm.b)$coeff[1,1]
  beta1.b[i] <- summary(fit.glm.b)$coeff[2,1]
  beta2.b[i] <- summary(fit.glm.b)$coeff[3,1]
}

# observed slope
(beta2.obs<-summary(satfullfit)$coeff[3,1])

# The distribution of the test statistic under the null
par(mfrow=c(1,1))
hist(beta2.b,nclass=50)
lines(c(beta2.obs, beta2.obs), c(0,1000), col=2)

# The Monte carlo p value.. 0.7564244, greater than 0.05
(1 + sum(abs(beta2.b) > abs(beta2.obs)))/(B+1)



################################################################################
Q2
####################################

library(dplyr)

attach(sleep)

extra<- c(0.7,-1.6,-0.2,-1.2,-0.1,3.4,3.7,0.8,0.0,2.0,
         1.9,0.8,1.1,0.1,-0.1,4.4,5.5,1.6,4.6,4.3)
group<-c(rep(1,10),rep(2,10))
ID <- c(1:20)
(sleepdf <- data.frame(extra,group,ID))

# classical two-samples t-test for two independent samples.
#t.test(extra ~ group, var.equal = TRUE) 

(group1 <- sleepdf %>%
  filter(group==1) %>%
  select(extra))

(group2 <- sleepdf %>%
    filter(group==2) %>%
    select(extra))

group1 <- group1$extra
group2 <- group2$extra

## classical test
mean(group1)
mean(group2)

boxplot(group1, group2)

t.test(group1, group2, var.equal = TRUE) 


### Non parametric bootstrap test using the two-samples t-test statistic
# as a test statistic.

(t.obs <- t.test(group1, group2)$statistic) # extract test statistic

z <- c(group1, group2)

m <- length(group1)
n <- length(group2)
mn <- m + n

B <- 2500

t.boot <- c(1:B)
for (b in 1:B)
{
  z.b <- sample(z, size=mn, replace=T)
  x.b <- z.b[1:n]
  y.b <- z.b[(n+1):mn]
  t.boot[b]<-t.test(x.b,y.b)$statistic
}

hist(t.boot,nclass=50,probability=T)
lines(c(t.obs,t.obs),c(0,1),lwd=3,col=2)

# bootstrap p value
(Pmc <- (1+sum(t.boot > t.obs))/(B+1))  # 0.9584166


### 2.3 -----------------------------------------------------------------------

# define test statistic tm
group1 <- sleepdf$extra[sleepdf$group == 1]
group2 <- sleepdf$extra[sleepdf$group == 2]
M1 <- median(group1)
M2 <- median(group2)

SM_x <- sum(abs(group1 - M1)) + sum(abs(group2 - M2))
(tm.obs <- (M1 - M2) / SM_x)

(sd.group1 <- sd(group1))
(sd.group2 <- sd(group2))

# Parametric bootstrap
set.seed(1423) 
B <- 10000  
tm.boot <- c(1:B)
(pooled.mean <- mean(c(group1, group2)))


for (b in 1:B) {
  # Simulate new data under the null (common distribution)
  group1.boot <- rnorm(n, pooled.mean, sd.group1)
  group2.boot <- rnorm(n, pooled.mean, sd.group2)
  M1.boot <- median(group1.boot)
  M2.boot <- median(group2.boot)
  SM.boot <- sum(abs(group1.boot - M1.boot)) + sum(abs(group2.boot - M2.boot))
  tm.boot[b] <- (M1.boot - M2.boot) / SM.boot
}

hist(tm.boot,nclass=50,probability=T)
lines(c(tm.obs,tm.obs),c(0,10),lwd=3,col=2)

(pmc <- (1+sum(abs(tm.boot) >= abs(tm.obs)))/(B+1)) 


## 2.4 -----------------------------------------------------------------

par(mfrow = c(1, 2))

hist(t.boot,nclass=50,probability=T)
lines(c(t.obs,t.obs),c(0,1),lwd=3,col=2)

hist(tm.boot,nclass=50,probability=T)
lines(c(tm.obs,tm.obs),c(0,10),lwd=3,col=2)

par(mfrow = c(1, 1))


#############################################################################
## Question 3
##------------------------------------------------
ID<-c(1:10)
x1<-c(0.8,-1.23,1.25,-0.28,-0.03,0.61,1.43,-0.54,-0.35,-1.60)
x2<-c(0.64,-1.69,1.47,-0.14,-0.18,0.43,1.61,-0.31,-0.38,-1.82)
(q3.df <- data.frame(ID,x1,x2))

# 3.1 --------------------------------------
# Calculate the means
(mean_x1 <- mean(q3.df$x1))
(mean_x2 <- mean(q3.df$x2))

(theta.hat <- mean_x1 / mean_x2)

### 3.2 ------------------------------------




