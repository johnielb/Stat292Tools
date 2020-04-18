summ_bin <- function(size, prob) {
  cat(paste0("\nX ~ Bin(", size, ", ", prob, ")\n> E(X) = ", size*prob,
             "\n> Var(X) = ", size*prob*(1-prob), "\n> sd(x) = ", sqrt(size*prob*(1-prob)), "\n"))
}

mass_bin <- function(x, size, prob) {
  cat(paste0("\nX ~ Bin(", size, ", ", prob, ")\n> P(X = ", x, ") = ", 
               dbinom(x, size, prob), "\n"))
}

cuml_bin <- function(x, size, prob, left.tail=TRUE, exclusive=FALSE) {
  str <- paste0("\nX ~ Bin(", size, ", ", prob, ")\n> P(X ")
  if (left.tail) {
    if (exclusive) {
      x <- x-1
    }
    str <- paste0(str, "<= ", x, ") = ",
                  pbinom(x, size, prob))
  } else {
    if (!exclusive) {
      x <- x-1
    }
    str <- paste0(str, "> ", x, ") = ", 
                  pbinom(x, size, prob, lower=F))
  }
  cat(paste0(str, "\n"))
}

range_bin <- function(lower_bound, upper_bound, size, prob, norm.approx=FALSE,
                      exclusive.lower=FALSE, exclusive.upper=FALSE) {
  mean <- size*prob
  var <- size*prob*(1-prob)
  str <- paste0("\nX ~ Bin(", size, ", ", prob, ")")
  warning <- ""
  
  if (norm.approx) {
    if (mean < 5 || size*(1-prob) < 5) warning <- "WARNING: Sample size may not be large enough.\n"
      
    str <- paste0(str, " is approximated as Y ~ N(", mean, ", ", var, ")\n")
    final_lower_bound <- lower_bound
    final_upper_bound <- upper_bound
    if (exclusive.lower) {
      final_lower_bound <- final_lower_bound + 0.5
      initial_inequality <- paste0(lower_bound, " < X")
    } else {
      final_lower_bound <- final_lower_bound - 0.5
      initial_inequality <- paste0(lower_bound, " <= X")
    }
    if (exclusive.upper) {
      final_upper_bound <- final_upper_bound - 0.5
      initial_inequality <- paste0(initial_inequality, " < ", upper_bound)
    } else {
      final_upper_bound <- final_upper_bound + 0.5
      initial_inequality <- paste0(initial_inequality, " <= ", upper_bound)
    }
    str <- paste0(str, "> P(", initial_inequality, ") = P(", final_lower_bound, " < Y < ", final_upper_bound, ") = ",
                  pnorm(final_upper_bound, mean, sqrt(var)) - pnorm(final_lower_bound, mean, sqrt(var)))
  } else {
    str <- paste0(str,"\n")
    final_lower_bound <- lower_bound
    final_upper_bound <- upper_bound
    if (exclusive.lower) {
      initial_inequality <- paste0(lower_bound, " < X")
    } else {
      final_lower_bound <- final_lower_bound - 1
      initial_inequality <- paste0(lower_bound, " <= X")
    }
    if (exclusive.upper) {
      final_upper_bound <- final_upper_bound - 1
      initial_inequality <- paste0(initial_inequality, " < ", upper_bound)
    } else {
      initial_inequality <- paste0(initial_inequality, " <= ", upper_bound)
    }
    str <- paste0(str, "> P(", initial_inequality, ") = P(X <= ", final_upper_bound, ") - P(X <= ", final_lower_bound,  ") = ",
                  pbinom(final_upper_bound, size, prob) - pbinom(final_lower_bound, size, prob))
  }
  cat(paste0(warning, str, "\n"))
}

range_bin(3, 11, 20, 0.25, norm.approx=T, exclusive.lower=T, exclusive.upper=T)
range_bin(3, 11, 20, 0.25, norm.approx=T, exclusive.lower=F, exclusive.upper=F)
range_bin(3, 11, 20, 0.25, norm.approx=F, exclusive.lower=F, exclusive.upper=T)
range_bin(3, 11, 20, 0.25, norm.approx=F, exclusive.lower=T, exclusive.upper=T)

summ_bin(4, 12, 0.2)
mass_bin(4, 12, 0.2)
cuml_bin(4, 12, 0.2, left.tail=T, exclusive=F)
cuml_bin(4, 12, 0.2, left.tail=T, exclusive=T)
cuml_bin(4, 12, 0.2, left.tail=F, exclusive=F)
cuml_bin(4, 12, 0.2, left.tail=F, exclusive=T)

# Assignment 1, question 2b
summ_bin(40, 0.2)
# Assignment 1, question 2c
mass_bin(20, 40, 0.2)
# Assignment 1, question 2d
cuml_bin(10, 40, 0.2, exclusive = T)
range_bin(-Inf, 10, 40, 0.2, exclusive.upper = T)

z_from_confidence_level <- function(confidence.level) {
  significance_level <- 1-confidence.level
  qnorm(1-significance_level/2)
}

proportion_ci <- function(y, n, z=NA, confidence.level=NA, ac=FALSE) {
  if (is.na(z)) { # if no z-value specified, calculate z-value from confidence level
    if (is.na(confidence.level)) {
      cat("No z-value or confidence level specified.")
      return(NA)
    }
    z <- z_from_confidence_level(confidence.level)
  }
  if (is.na(confidence.level)) { # if no confidence level specified, calculate from z-value
    confidence.level <- pnorm(z)
  }
  str <- paste0(confidence.level*100, "% Clopper-Pearson confidence interval: ")
  if (ac) {
    str <- paste0(confidence.level*100, "% Agresti-Coull adjusted confidence interval: ")
    y <- y + 2
    n <- n + 4
  }
  p_hat <- y/n
  error <- z * sqrt((p_hat*(1-p_hat))/n)
  interval <- c(p_hat-error, p_hat+error)
  cat(paste0(str, "(", interval[1], ", ", interval[2], ")\n"))
}

# Slide 47, discrete distributions
proportion_ci(22, 300, confidence.level = 0.95)
proportion_ci(22, 300, z = 1.645)
proportion_ci(22, 300, confidence.level = 0.95, ac=T)
# Assignment 1, question 3b
proportion_ci(32, 172, z=1.96)
proportion_ci(32, 172, z=1.96, ac=T)

min_sample_size <- function(margin_of_error, confidence.level, p=0.5) {
  z <- z_from_confidence_level(confidence.level)
  cat(paste0("> Sample size should be greater or equal to ", (z/margin_of_error)^2*p*(1-p)),"\n")
}

# selection from slide 53, discrete distributions
min_sample_size(0.1, confidence.level=0.95)
min_sample_size(0.09, confidence.level=0.99)
min_sample_size(0.03, confidence.level=0.95, p=0.65)

proportion_z_test <- function(y1, n1, y2, n2, sig.level=0.05) {
  p1 <- y1/n1
  p2 <- y2/n2
  p_hat <- (y1+y2)/(n1+n2)
  p_hat_variance <- sqrt(p_hat*(1-p_hat)*(1/n1+1/n2))
  z <- abs((p1-p2)/p_hat_variance)
  p_val <- 2 * (1 - pnorm(z))
  result <- "\n> Can't reject null hypothesis, insufficient evidence to suggest the two proportions are different.\n"
  if (p_val < sig.level) {
    result <- "\n> Can reject null hypothesis, sufficient evidence to suggest the two proportions are different.\n"
  }
  cat(paste0("Variance of common proportion: ", p_hat_variance, "\nTest statistic: ", z, "\np-value of test: ", p_val, result))
}

# Slide 58, discrete distributions
proportion_z_test(3113, 24117, 8391, 35829)
# Assignment 1, question 3c
proportion_z_test(32, 172, 32, 163)

summ_poi <- function(x, lambda) {
  cat(paste0("\nX ~ Poi(", lambda, ")\n> E(X) = ", lambda,
             "\n> Var(X) = ", lambda, "\n> sd(x) = ", sqrt(lambda), "\n"))
}

mass_poi <- function(x, lambda) {
  cat(paste0("\nX ~ Poi(", lambda, ")\n> P(X = ", x, ") = ", 
             dpois(x, lambda), "\n"))
}

cuml_poi <- function(x, lambda, left.tail=TRUE, exclusive=FALSE) {
  str <- paste0("\nX ~ Poi(", lambda, ")\n> P(X ")
  if (left.tail) {
    if (exclusive) {
      str <- paste0(str, "< ", x, ") = ",
                    ppois(x-1, lambda))
    } else {
      str <- paste0(str, "<= ", x, ") = ",
                    ppois(x, lambda))
    }
  } else {
    if (exclusive) {
      str <- paste0(str, "> ", x, ") = ", 
                    ppois(x, lambda, lower=F))
    } else {
      str <- paste0(str, ">= ", x, ") = ", 
                    ppois(x-1, lambda, lower=F))
    }
  }
  cat(paste0(str, "\n"))
}

summ_poi(16, 12)
mass_poi(16, 12)
cuml_poi(16, 12, left.tail=T, exclusive=F)
cuml_poi(16, 12, left.tail=T, exclusive=T)
cuml_poi(16, 12, left.tail=F, exclusive=F)
cuml_poi(16, 12, left.tail=F, exclusive=T)

chisq_gof_pois_test <- function(n, f, sig.level=0.05) {
  lambda_hat <- sum(n*f)/sum(f)
  pois <- dpois(n, lambda_hat)
  pois[length(pois)] <- 1-sum(pois[1:length(pois)-1])
  f_hat_orig <- sum(f) * pois
  # carry over only those with an expected f > 5
  f_hat <- f_hat_orig[f_hat_orig>5]
  f_new <- f[f_hat_orig>5]
  if (!identical(f_hat_orig, f_hat)) { # regroup all groups less than 5, assuming all on right tail
    new_last_f_hat <- sum(f) * (1-sum(pois[f_hat_orig>5])) # total times Poisson probability not carried over
    new_last_f <- sum(f[f_hat_orig<=5]) # all observations not carried over
    if (new_last_f_hat < 5) { # only need to do this once as f_hat[length(f_hat)] is guaranteed to be >5
      new_last_f_hat <- new_last_f_hat + sum(f) * pois[length(f_hat)]
      new_last_f <- new_last_f + f[length(f_hat)]
    }
    # add these to the end
    f_hat[length(f_hat)+1] <- new_last_f_hat 
    f_new[length(f_new)+1] <- new_last_f
  }
  mse <- (f_new-f_hat)^2/f_hat
  df <- length(f_hat)-1-1
  p_val <- pchisq(sum(mse), df=df, lower.tail=F)
  result <- "\n> Can't reject null hypothesis, insufficient evidence to suggest distribution is not Poisson.\n"
  if (p_val < sig.level) {
    result <- "\n> Can reject null hypothesis, sufficient evidence to suggest distribution is not Poisson.\n"
  }
  cat(paste0("Test statistic: ", sum(mse), " with ", df, " degrees of freedom.\np-value of test: ", p_val, result))
}

# Slide 71, discrete distributions
chisq_gof_pois_test(0:10, c(114,25,15,10,6,5,2,1,1,0,1))
# Assignment 2, 1c
chisq_gof_pois_test(0:5, c(19,3,1,4,7,7))
