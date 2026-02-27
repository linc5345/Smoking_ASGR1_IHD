# Load package
library(data.table)

# Ensure your input file has the following from prior MR analyses:
# a_beta: exposure -> mediator     (e.g., LSI -> ASGR1, excluding ASGR1 gene region)
# a_se: SE of exposure -> mediator 
# b_beta: mediator -> outcome      (e.g., LSI -> MI)
# b_se: SE of mediator -> outcome
# c_beta: exposure -> outcome      (e.g., LSI -> MI)
# c_se: SE of exposure -> outcome

# Optional: inclusion of covariances (cov_ab, cov_ac, cov_bc)

DT  <- fread("path/to/input.csv")

# Define output directory
outfile <- "path/to/output.csv"

# small epsilon to flag unstable denominators
c_eps <- 1e-8

# Ensure optional covariances exist, else set to 0
if (!"cov_ab" %in% names(DT)) DT[, cov_ab := 0]
if (!"cov_ac" %in% names(DT)) DT[, cov_ac := 0]
if (!"cov_bc" %in% names(DT)) DT[, cov_bc := 0]

n <- nrow(DT)
# Prepare output columns
DT[, `:=`(
  prop_med = as.numeric(NA),
  se_prop_med = as.numeric(NA),
  prop_lci = as.numeric(NA),
  prop_uci = as.numeric(NA),
  z = as.numeric(NA),
  pval = as.numeric(NA),
  unstable = logical(n)
)]

# Delta method for each row
for (i in seq_len(n)) {
  a  <- DT$a_beta[i]; sa <- DT$a_se[i]
  b  <- DT$b_beta[i]; sb <- DT$b_se[i]
  c  <- DT$c_beta[i]; sc <- DT$c_se[i]
  cab <- DT$cov_ab[i]; cac <- DT$cov_ac[i]; cbc <- DT$cov_bc[i]
  # skip if any NA in required
  if (any(is.na(c(a,sa,b,sb,c,sc)))) {
    DT$unstable[i] <- TRUE
    next
  }
  # if c is zero/near-zero flag unstable and skip
  if (is.na(c) || abs(c) < c_eps) {
    DT$unstable[i] <- TRUE
    next
  }
  # point estimate: (a * b) / c
  est <- (a * b) / c
  # partial derivatives
  d_a <- b / c
  d_b <- a / c
  d_c <- - (a * b) / (c^2)
  # variance using delta method, including covariances
  var_est <- (d_a^2) * (sa^2) + (d_b^2) * (sb^2) + (d_c^2) * (sc^2) +
    2 * d_a * d_b * cab + 2 * d_a * d_c * cac + 2 * d_b * d_c * cbc
  # numerical stability
  if (is.na(var_est) || var_est < 0) {
    # negative by tiny numeric noise -> set to 0 if very small negative
    if (!is.na(var_est) && abs(var_est) < 1e-12) var_est <- 0 else {
      DT$unstable[i] <- TRUE
      next
    }
  }
  se_est <- sqrt(var_est)
  lci <- est - 1.96 * se_est
  uci <- est + 1.96 * se_est
  zval <- ifelse(se_est == 0, NA_real_, est / se_est)
  pval <- ifelse(is.na(zval), NA_real_, 2 * (1 - pnorm(abs(zval))))
  # store
  DT$prop_med[i]      <- est
  DT$se_prop_med[i]   <- se_est
  DT$prop_lci[i]      <- lci
  DT$prop_uci[i]      <- uci
  DT$z[i]             <- zval
  DT$pval[i]          <- pval
  DT$unstable[i]      <- FALSE
}

# Reorder output columns to be informative
out_cols <- c(names(DT)[!names(DT) %in% c("prop_med","se_prop_med","prop_lci","prop_uci","z","pval","unstable")],
              "prop_med","se_prop_med","prop_lci","prop_uci","z","pval","unstable")
outDT <- DT[, ..out_cols]

# Write output
fwrite(outDT, outfile)