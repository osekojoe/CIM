### ----- 4.1 
ID <- c(1:10)
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)
q3df <- data.frame(ID,x1,x2, d=(x1-x2))

(d <- x1 - x2)
(mean_d <- mean(d))
(sd_d <- sd(d))

n <- length(d)
df <- n - 1

(t_value <- qt(0.975, df)) # critical value
(se_d <- sd_d / sqrt(n)) # Standard error of the mean difference

# Calculate the confidence interval
(lower_ci <- mean_d - t_value * se_d)
(upper_ci <- mean_d + t_value * se_d)

c(t.hat-se.t.hat*C.alpha,t.hat+se.t.hat*C.alpha)

### - 4.2 -------------------------------------------------------
# Percentile method ------
(mud.obs <- mean(x1) - mean(x2))

set.seed(2025)
B <- 10000
n <- length(x1)
mud.boot <- c(1:B)

for (b in 1:B) {
  indices <- sample(1:n, size = n, replace = TRUE)
  x1.boot <- x1[indices]
  x2.boot <- x2[indices]
  mud.boot[b] <- mean(x1.boot) - mean(x2.boot)
}


(CI.P <- quantile(mud.boot, probs = c(0.025, 0.975)))

lo <- quantile(mud.boot,probs=c(0.05,0.95))[1]
up <- quantile(mud.boot,probs=c(0.05,0.95))[2]
hist(mud.boot, probability=T,nclass=50)
lines(c(lo,lo),c(0,5),col=4)
lines(c(up,up),c(0,5),col=4)

# bootstrap t method --------

z <- q3df$d
t.hat <- mean(q3df$d)
se.t.hat<-sqrt(var(z)/10)

set.seed(2025)
B <- 10000
n <- length(x1)
mud.boot <- c(1:B)

for (b in 1:B) {
  x.boot <- sample(z, size = length(z), replace = TRUE)
  se.boot <- sqrt(var(x.boot) / length(z))
  mud.boot[b] <- (mean(x.boot) - t.hat) / se.boot
}

# Quantiles for bootstrap replicates
quantile(mud.boot, probs=c(0.05,0.95))

# bootstrap t interval
up <- quantile(mud.boot,probs=c(0.95))
lo <- quantile(mud.boot,probs=c(0.05))
c(t.hat + se.t.hat*lo, t.hat + se.t.hat*up)

t.test(z, conf.level=0.95)

# BCa method ----------------------------------------
library(boot)
bca_ci <- boot.ci(boot.out = boot(data = data.frame(x1, x2), statistic = function(data, indices) {
  sample_x1 <- data$x1[indices]
  sample_x2 <- data$x2[indices]
  return(mean(sample_x1) - mean(sample_x2))
}, R = B), type = "bca")
bca_ci

## 4.3 
(p_value <- mean(abs(bootstrap_diffs) >= abs(observed_diff)))

