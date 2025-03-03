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


# Classical method for standard error
classical_se <- function(x1, x2) {
  n <- length(x1)
  mean_x1 <- mean(x1)
  mean_x2 <- mean(x2)
  
  var_x1 <- var(x1)
  var_x2 <- var(x2)
  
  se <- sqrt((1 / mean_x2)^2 * (var_x1 / n) + 
               (-mean_x1 / mean_x2^2)^2 * (var_x2 / n))
  se
}

# Calculate classical SE
classical_se(x1, x2)



### 3.2 ----------------------------------------------------
## Non-parameteric bootstrap 
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
  set.seed(2025)
  n <- length(x1)
  theta_boot <- numeric(B)
  for (b in 1:B) {
    # Resample with replacement
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    
    # Check for mean(x2) != 0
    if (mean(x2_sample) >= 0.0001) {
      theta_boot[b] <- theta_hat(x1_sample, x2_sample)
    } else {
      theta_boot[b] <- NA 
    }
  }
  #sd(na.omit(theta_boot))/sqrt(n)
  sqrt(var(na.omit(theta_boot)))
}

bootstrap_B <- c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000)
(bootstrap_se_values <- sapply(bootstrap_B, function(B) bootstrap_se(x1, x2, B)))

# Combine results
(results <- data.frame(B = bootstrap_B, Bootstrap_SE = bootstrap_se_values))

# Plot Bootstrap SE vs B
plot(results$B, results$Bootstrap_SE, type = "b", log = "x",
     xlab = "B", ylab = "Boot. Std. Error",
     main = "Bootstrap SE vs B")
grid()


##--------Jacknife ----------------------------

(n <- length(x1))
m.jack <- c(1:n)

for(i in 1:n) {
  cat(i)
  x1.jack <- x1[ - c(i)]
  x2.jack <- x2[ - c(i)]
  m.jack[i] <- mean(x2.jack)/mean(x1.jack)
}
mean(m.jack)
# for bias
(n - 1) * (mean(m.jack) - theta.obs) # 5.819458
# for standard error
(se.jack <- sqrt( ((n-1)/n)*sum((m.jack-mean(m.jack))^2) )) # 3.400245


### 3.3 ---------------------------------------------------------------
# Confidence interval using bootstrap estimates
# Data
ID <- c(1:10)
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)
q3df <- data.frame(ID,x1,x2, d=(x1-x2), theta=(x1/x2))

n <- length(x1)
set.seed(2025)
B = 10000
theta_boot<-c(1:B)

for (b in 1:B) {
  # Resample with replacement
  indices <- sample(1:n, size = n, replace = TRUE)
  x1_sample <- x1[indices]
  x2_sample <- x2[indices]
  theta_boot[b] <- mean(x1_sample) / mean(x2_sample)
}

t.up<-quantile(theta_boot,probs=c(0.975))
t.lo<-quantile(theta_boot,probs=c(0.025))
c(t.up,t.lo)   # [-1.476429 , 3.276681]

hist(theta_boot,nclass=50,col=0,probability=TRUE)
lines(c(t.up,t.up),c(0,20),col=2)
lines(c(t.lo,t.lo),c(0,10),col=2)


#### 3.4 ------------------------------------------------------------------
# Hypothesis test using bootstrap
# Data
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)

# Observed theta
theta.obs <- mean(x1) / mean(x2)

# Transform x2 under H0 to ensure theta = 1
x2.null <- x2 * (mean(x1) / mean(x2))
x2.null

#mean(x1) / mean(x2.null)

# Bootstrap procedure
n <- length(x1)
B = 10000
theta.boot <- c(1:B)

set.seed(2025)

for (b in 1:B) {
  indices <- sample(1:n, size = n, replace = TRUE)
  x1.boot <- x1[indices]
  x2.boot <- x2.null[indices]
  theta.boot[b] <- mean(x1.boot) / mean(x2.boot)
}

#(pmc <- (1+sum(abs(theta.boot) > abs(theta.obs)))/(B+1))
(pmc <- (1+sum(theta.boot >= theta.obs))/(B+1))

