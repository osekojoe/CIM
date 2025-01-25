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

png("pictures/boxplot2.png", width = 18, 
    height = 10, units = "cm", res = 300)
boxplot(group1, group2)
dev.off()

t.test(group1, group2, var.equal = TRUE) 

### ------2.2 -------------------------------------------------------------------
### Non parametric bootstrap test using the two-samples t-test statistic
# as a test statistic.

(t.obs <- t.test(group1, group2)$statistic) # extract test statistic

z <- c(group1, group2)
z

m <- length(group1)
n <- length(group2)
mn <- m + n

set.seed(2025)
B <- 10000

t.boot <- c(1:B)
for (b in 1:B)
{
  z.b <- sample(z, size=mn, replace=T)
  x.b <- z.b[1:n]
  y.b <- z.b[(n+1):mn]
  t.boot[b]<-t.test(x.b,y.b)$statistic
}

png("pictures/hist2.1.png", width = 18, 
    height = 10, units = "cm", res = 300)
hist(t.boot,nclass=50,probability=T)
lines(c(t.obs,t.obs),c(0,1),lwd=3,col=2)
dev.off()

# bootstrap p value
(Pmc <- (1+sum(abs(t.boot) > abs(t.obs)))/(B+1))  # 0.5593441


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

(pmc <- (1+sum(abs(tm.boot) >= abs(tm.obs)))/(B+1))  # 0.2143786


## 2.4 -----------------------------------------------------------------


png("pictures/hist2.3.png", width = 18, 
    height = 10, units = "cm", res = 300)

par(mfrow = c(1, 2))

hist(t.boot,nclass=50,probability=T, main = "Non-parametric bootstrap")
lines(c(t.obs,t.obs),c(0,1),lwd=3,col=2)

hist(tm.boot,nclass=50,probability=T, main = "Parametric bootstrap")
lines(c(tm.obs,tm.obs),c(0,10),lwd=3,col=2)

dev.off()

par(mfrow = c(1, 1))