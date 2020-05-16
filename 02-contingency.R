library(tidyr)

z_from_confidence_level <- function(confidence.level) {
  significance_level <- 1-confidence.level
  qnorm(1-significance_level/2)
}

odds_ratio <- function(n11, n12, n21, n22, z=NA, confidence.level=NA) {
  theta_hat <- (n11*n22)/(n12*n21)
  if (is.na(z)) { # if no z-value specified, calculate z-value from confidence level
    if (is.na(confidence.level)) {
      cat("No z-value or confidence level specified.")
      return(NA)
    }
    z <- z_from_confidence_level(confidence.level)
  }
  if (is.na(confidence.level)) { # if no confidence level specified, calculate from z-value
    confidence.level <- 2*pnorm(z)-1
  }
  se <- z * sqrt(1/n11+1/n12+1/n21+1/n22)
  interval <- c(exp(log(theta_hat) - se), exp(log(theta_hat) + se))
  cat(paste0("Odds ratio estimate = ", theta_hat, "\n", confidence.level*100, "% confidence interval for odds ratio: (", interval[1], ", ", interval[2], ")\n\n"))
}

expected_fs <- function(tbl, check=F) {
  row_marginals <- apply(tbl, 1, sum)
  col_marginals <- apply(tbl, 2, sum)
  total <- sum(row_marginals)
  found <- F
  f_hat <- matrix(ncol=length(col_marginals), nrow=length(row_marginals))
  for (i in 1:length(row_marginals)) {
    for (j in 1:length(col_marginals)) {
      f_hat[i,j] <- row_marginals[[i]]*col_marginals[[j]]/total;
      if (f_hat[i,j] < 5) {
        found <- T
        if (check) cat(paste0("n",i,j," = ",tbl[[i]][[j]], " has an expected frequency of ", f_hat, ", smaller than 5.\n\n"))
      }
    }
  }
  if (check && !found) cat("All expected frequencies are at least 5.\n\n")
  if (!check) f_hat
}

chisq_independence_test <- function(tbl, pearson=T, likelihood=T) {
  f_hat <- expected_fs(tbl)
  df <- (length(tbl[1,])-1)*(length(tbl[,1])-1)
  chisq <- (tbl-f_hat)^2/f_hat
  like_ratio <- 2*(tbl*log(tbl/f_hat))
  p_value_chi <- pchisq(sum(chisq), df=df, lower.tail=F)
  p_value_g <- pchisq(sum(like_ratio), df=df, lower.tail=F)
  if (pearson) {
    cat(paste0("Pearson chi-square test: X^2 = ",sum(chisq),"\nP(chisq_",df," > X^2) = ",p_value_chi,"\n\n"));
  }
  if (likelihood) {
    cat(paste0("Likelihood ratio chi-square test: G^2 = ",sum(like_ratio),"\nP(chisq_",df," > G^2) = ",p_value_g),"\n\n");
  }
}

cancer <- matrix(c(13,3,2,12), byrow=T, nrow=2)
colnames(cancer) <- c("Survived 10 years+", "Survived less than 10 years")
odds_ratio(cancer[1,1], cancer[1,2], cancer[2,1], cancer[2,2], confidence.level = 0.95)
expected_fs(cancer, check=T)
chisq_independence_test(cancer)

pregnancies <- tibble('alcohol' = c(0,0.055,0.33, 2), '0' = c(105,58,84,47), '8' = c(7,5,37,16), '16' = c(11,13,42,17)) 

pregnancies_long <- pregnancies %>%
  pivot_longer(-alcohol, names_to = 'nicotine', values_to = 'frequency')
pregnancies_long$nicotine <- as.numeric(pregnancies$nicotine)

fishers_exact_test <- function(tbl, sig.level=0.05, lower.tail=F, two.sided=F) {
  row_marginals <- apply(tbl, 1, sum)
  x <- tbl[1,1]
  m <- row_marginals[[1]]
  n <- row_marginals[[2]]
  k <- apply(tbl, 2, sum)[[1]]
  table_p_val <- dhyper(x,m,n,k)
  if (two.sided) { # get all p-values as extreme as the table p value
    p_vals <- dhyper(0:min(m,k), m, n, k)
    p_val <- sum(p_vals[p_vals<=table_p_val])
    mid_p_val <- two_sided_p_val - table_p_val
  } else {
    if (lower.tail) {
      p_val <- phyper(x-1,m,n,k, lower.tail=F)
    } else {
      p_val <- phyper(x,m,n,k)
    }
    mid_p_val <- p_val - 0.5*table_p_val
  }
  cat(paste0("p-value of test: ", p_val, "\nmid p-value of test: ", mid_p_val))
  if (mid_p_val < sig.level) {
    cat("\n> Can reject null hypothesis, sufficient evidence to suggest the two proportions are different.\n\n")
  } else {
    cat("\n> Can't reject null hypothesis, insufficient evidence to suggest the two proportions are different.\n\n")
  }
}

fishers_exact_test(tbl)
fishers_exact_test(tbl, two.sided=T)

conditional_associations <- function(tbls) {
  for (i in seq_along(tbls)) {
    cat(paste0("Z = ", names(tbls)[[i]], ": "))
    m <- tbls[[i]]
    odds_ratio(m[1,1], m[1,2], m[2,1], m[2,2], confidence.level = 0.95)
  }
}

# construct 3 way contingency table
death_penalty <- list(matrix(c(19,132,11,52), byrow=T, nrow=2),
             matrix(c(0,9,6,97), byrow=T, nrow=2)) %>% 
  lapply(function(m) {
    rownames(m) <- c("Cauc defendant", "AA defendant")
    colnames(m) <- c("No death penalty", "Death penalty")
    m
  })
names(death_penalty) <- c("African-American defendant", "Caucasian defendant")

marginal_tbl <- Reduce('+', death_penalty)
odds_ratio(marginal_tbl[1,1], marginal_tbl[1,2], marginal_tbl[2,1], marginal_tbl[2,2], confidence.level = 0.95)

conditional_associations(death_penalty)


graduates <- list(matrix(c(18,12,12,8), byrow=T, nrow=2),
                  matrix(c(2,8,8,32), byrow=T, nrow=2)) %>% 
  lapply(function(m) {
    rownames(m) <- c("Female", "Male")
    colnames(m) <- c("Low income", "High income")
    m
  })
names(graduates) <- c("Arts", "Science")

marginal_tbl <- Reduce('+', graduates)
odds_ratio(marginal_tbl[1,1], marginal_tbl[1,2], marginal_tbl[2,1], marginal_tbl[2,2], confidence.level = 0.95)

conditional_associations(graduates)

smokers <- list(matrix(c(647,2,622,27), byrow=T, nrow=2),
                matrix(c(41,19,28,32), byrow=T, nrow=2)) %>% 
  lapply(function(m) {
    rownames(m) <- c("Lung cancer", "No lung cancer")
    colnames(m) <- c("Smoker", "Non smoker")
    m
  })
names(smokers) <- c("Male", "Female")