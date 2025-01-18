
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

(theta.obs <- mean_x1 / mean_x2)


### 3.2 ----------------------------------------------------
## Non-parameteric bootstrap 
# 12.3 Example: the patch data

# B_values <- c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000)


# Function to compute theta hat
#set.seed(235)
theta_hat <- function(data) {
  mean(data$x1) / mean(data$x2)
}
# Bootstrap function
bootstrap_se <- function(data, B) {
  thetas <- numeric(B)
  for(i in 1:B) {
    boot_data <- data[sample(nrow(data), replace = TRUE), ]
    thetas[i] <- theta_hat(boot_data)
  }
  # Remove NA values or Inf if any before calculating SD
  thetas <- thetas[!is.na(thetas) & !is.infinite(thetas)]
  if(length(thetas) < 2) return(NA)  # If there aren't enough valid samples
  return(sd(thetas))
}
# Different values of B
B_values <- c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000)
se_bootstrap <- sapply(B_values, function(B) bootstrap_se(q3.df, B))
# Print results
for(i in seq_along(B_values)) {
  cat("Bootstrap SE with B =", B_values[i], ":", se_bootstrap[i], "\n")
}


##--------Jacknife ----------------
(n <- length(x1))
m.jack <- c(1:n)

for(i in 1:n) {
  cat(i)
  x1.jack <- x1[ - c(i)]
  x2.jack <- x2[ - c(i)]
  m.jack[i] <- mean(x2.jack)/mean(x1.jack)
}

(n - 1) * (mean(m.jack) - theta.obs)


# Jackknife function
jackknife_se <- function(x1, x2) {
  n <- length(x1)
  jackknife_thetas <- numeric(n)
  
  for (i in 1:n) {
    x1_jack <- x1[-i]
    x2_jack <- x2[-i]
    jackknife_thetas[i] <- mean(x1_jack) / mean(x2_jack)
  }
  
  # Jackknife standard error
  se <- sqrt((n - 1) * var(jackknife_thetas))
  list(se = se, estimates = jackknife_thetas)
}

# Jackknife SE
jackknife_results <- jackknife_se(x1, x2)
jackknife_se <- jackknife_results$se
jackknife_se


### 3.2 -------------###########################################

library(plotrix)
# Define the dataset
ID <- c(1:10)
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)

# Calculate the ratio statistic
theta_hat <- function(x1, x2) {
  mean(x1) / mean(x2)
}

theta_hat(x1, x2)

# Bootstrap standard error function
bootstrap_se <- function(x1, x2, B) {
  set.seed(3647)
  n <- length(x1)
  theta_boot <- numeric(B)
  for (b in 1:B) {
    # Resample with replacement
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    
    # Check for mean(x2) != 0
    if (mean(x2_sample) != 0) {
      theta_boot[b] <- theta_hat(x1_sample, x2_sample)
    } else {
      theta_boot[b] <- NA  # Assign NA if mean(x2) = 0
    }
  }
  sd(na.omit(theta_boot))/sqrt(n)
  #std.error(na.omit(theta_boot))  # Exclude NAs from sd calculation
}

#B=7500
#bootstrap_se(x1, x2, B)

bootstrap_B <- c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000)
(bootstrap_se_values <- sapply(bootstrap_B, function(B) bootstrap_se(x1, x2, B)))

# Combine results
(results <- data.frame(B = bootstrap_B, Bootstrap_SE = bootstrap_se_values))

# Plot Bootstrap SE vs B
plot(results$B, results$Bootstrap_SE, type = "b", log = "x",
     xlab = "B", ylab = "Boot. Std. Error",
     main = "Bootstrap SE vs B")
grid()

# Jackknife standard error function
jackknife_se <- function(x1, x2) {
  n <- length(x1)
  theta_jack <- numeric(n)
  for (i in 1:n) {
    # Leave one observation out
    theta_jack[i] <- theta_hat(x1[-i], x2[-i])
  }
  (n - 1) * sqrt(var(theta_jack) / n)
}

(jack.se <- jackknife_se(x1, x2))

### 3.3 ----------
# bootstrap confidence interval
bootstrap_ci <- function(x1, x2, B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(x1)
  theta_boot <- numeric(B)
  
  for (b in 1:B) {
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    
    if (mean(x2_sample) != 0) {
      theta_boot[b] <- mean(x1_sample) / mean(x2_sample)
    } else {
      theta_boot[b] <- NA 
    }
  }
  
  theta_boot <- na.omit(theta_boot)
  
  # 95% CI using the percentile method
  ci_lower <- quantile(theta_boot, probs = 0.025)
  ci_upper <- quantile(theta_boot, probs = 0.975)
  
  list(CI = c(ci_lower, ci_upper), Bootstrap_Estimates = theta_boot)
}

result <- bootstrap_ci(x1, x2, B = 10000, seed = 123)

# Print the confidence interval
print(result$CI)



#######===3.3 ======================================
# Confidence interval using bootstrap estimates
# Data
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)

# Function to calculate bootstrap estimates and CI (without histogram)
bootstrap_ci <- function(x1, x2, B) {
  set.seed(2523)  # Set seed if provided
  n <- length(x1)
  theta_boot <- numeric(B)
  
  for (b in 1:B) {
    # Resample with replacement
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    
    # Check for mean(x2) != 0
    if (mean(x2_sample) != 0) {
      theta_boot[b] <- mean(x1_sample) / mean(x2_sample)
    } else {
      theta_boot[b] <- NA  # Assign NA if mean(x2) = 0
    }
  }
  
  # Remove NA values
  theta_boot <- na.omit(theta_boot)
  
  ci <- quantile(theta_boot, probs = c(0.025, 0.975))
  
  return(list(CI = ci, Bootstrap_Estimates = theta_boot))
}

result <- bootstrap_ci(x1, x2, B = 10000)
print(result$CI)

### following prof code
n <- length(x1)
B = 10000
theta_boot<-c(1:B)

for (b in 1:B) {
  # Resample with replacement
  indices <- sample(1:n, size = n, replace = TRUE)
  x1_sample <- x1[indices]
  x2_sample <- x2[indices]
  theta_boot[b] <- mean(x1_sample) / mean(x2_sample)
}

t.up<-quantile(theta_boot,probs=c(0.95))
t.lo<-quantile(theta_boot,probs=c(0.05))
c(t.up,t.lo)
hist(theta_boot,nclass=10,col=0,probability=TRUE)
lines(c(t.up,t.up),c(0,20),col=2)
lines(c(t.lo,t.lo),c(0,10),col=2)




########=================================

#### 3.4 ------------------------------------------------------------------
# Hypothesis test using bootstrap
bootstrap_test <- function(x1, x2, B, theta_null = 1) {
  set.seed(5332)
  n <- length(x1)
  bootstrap_thetas <- numeric(B)
  
  for (b in 1:B) {
    indices <- sample(1:n, n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    bootstrap_thetas[b] <- mean(x1_sample) / mean(x2_sample)
  }
  
  # Compute p-value for one-sided test
  p_value <- mean(bootstrap_thetas >= theta_null)
  p_value
}

# Test with B = 1000 and H0: theta = 1
bootstrap_p_value <- bootstrap_test(x1, x2, B = 1000, theta_null = 1)
bootstrap_p_value


###########################################################################
## Q4

### ----- 4.1 
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)

(d <- x1 - x2)

(mean_d <- mean(d))
(sd_d <- sd(d))

n <- length(d)

df <- n - 1

t_value <- qt(0.975, df) # critical value
se_d <- sd_d / sqrt(n) # Standard error of the mean difference

# Calculate the confidence interval
(lower_ci <- mean_d - t_value * se_d)
(upper_ci <- mean_d + t_value * se_d)

### - 4.2 -------------------------------------------------------

B <- 10000
n <- length(x1)

(observed_diff <- mean(x1) - mean(x2))

# Function to calculate bootstrap mean difference
bootstrap_diff <- function(x1, x2, B) {
  n <- length(x1)
  bootstrap_diffs <- numeric(B)
  
  for (b in 1:B) {
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    bootstrap_diffs[b] <- mean(x1_sample) - mean(x2_sample)
  }
  
  return(bootstrap_diffs)
}

bootstrap_diffs <- bootstrap_diff(x1, x2, B)

# Percentile method
(percentile_ci <- quantile(bootstrap_diffs, probs = c(0.025, 0.975)))

##-----------------------------------
# Bootstrap t method
quantile(bootstrap_diffs,probs=c(0.05,0.95))
qt(0.95, 9)
qt(0.05, 9)

hist(bootstrap_diffs,probability=T,nclass=20,ylim=c(0,6),xlim=c(-0.5,0.5))
xx<-seq(from=-2,to=2,length=1000)
dx2<-dt(xx,6)
lines(xx,dx2,col=2)
lines(rep(qt(0.95,9),2),c(0,0.6),col=2)
calpha <- quantile(bootstrap_diffs, probs=c(0.95))
lines(rep(calpha,2),c(0,6),col=4)

##-----------------------------------

# BCa method
library(boot)
bca_ci <- boot.ci(boot.out = boot(data = data.frame(x1, x2), statistic = function(data, indices) {
  sample_x1 <- data$x1[indices]
  sample_x2 <- data$x2[indices]
  return(mean(sample_x1) - mean(sample_x2))
}, R = B), type = "bca")
bca_ci

## 4.3 

(p_value <- mean(abs(bootstrap_diffs) >= abs(observed_diff)))
