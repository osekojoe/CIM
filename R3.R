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
  #sqrt(var(theta_boot))
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

