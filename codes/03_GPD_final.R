# Computational Statistics
# Extreme Value Theory on river water levels

# Members of the team:
# Kristóf Andrási
# Gellért Banai
# Hunor Kuti
# Ákos Virág

# analysis

# ------------------------------------------------------------------------------
# importing neccesary libraries
# vector of required packages
pkgs <- c(
  "extRemes"          # main package we use for extreme value analysis (GEV/GPD, return levels, CIs, etc.)
)

# Install any that are missing
to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_install) > 0) {
  install.packages(to_install)   # install only those packages that are not already installed
}

# loading them
invisible(lapply(pkgs, library, character.only = TRUE))  # load each package in 'pkgs' so we can use their functions
rm(to_install)                                           # clean up helper variable, we don't need it later
# ------------------------------------------------------------------------------
# loading the dataframes
river_data <- read.csv("../data/river_data.csv")                 # full combined daily river data (both rivers, all years)
river_data_tisza_02_13 <- read.csv("../data/river_data_tisza_02_13.csv")  # Tisza daily data, 2002–2013 only
river_data_tisza_14_24 <- read.csv("../data/river_data_tisza_14_24.csv")  # Tisza daily data, 2014–2024 only
river_data_duna_02_13  <- read.csv("../data/river_data_duna_02_13.csv")   # Duna daily data, 2002–2013 only
river_data_duna_14_24  <- read.csv("../data/river_data_duna_14_24.csv")   # Duna daily data, 2014–2024 only

sum(is.na(river_data_duna_14_24))  # check how many missing values (NAs) we have in the Duna 2014–2024 dataset (any column)

# find the NA positions in the value column
na_idx <- which(is.na(river_data_duna_14_24$value))  # get the row indices where the 'value' column (water level) is NA

for (i in na_idx) {
  # index range for up to 9 neighbours (4 before, 4 after)
  # here we choose a small window around the missing value:
  # from i-4 to i+4, but cut off at the boundaries of the data frame
  start_i <- max(1, i - 4)                                   # left boundary of the window (at least row 1)
  end_i   <- min(nrow(river_data_duna_14_24), i + 4)         # right boundary of the window (at most last row)
  
  # take neighbour values in that window
  neighbours <- river_data_duna_14_24$value[start_i:end_i]   # extract the 'value' entries in this local window
  
  # drop any NAs inside the window
  neighbours <- neighbours[!is.na(neighbours)]               # keep only those neighbours that are not NA
  
  # replace NA with mean of neighbours
  # so each missing water level is replaced by the average of nearby days (up to 9 neighbours)
  river_data_duna_14_24$value[i] <- mean(neighbours)
}

# check that NAs are gone
sum(is.na(river_data_duna_14_24$value))  # confirm that there are now 0 missing values in the 'value' column

# helper function
# This function fits a Generalized Pareto Distribution (GPD) using the Peaks-Over-Threshold (POT) method
# for either the upper tail (maxima) or the lower tail (minima) of a river level time series.
fit_gpd_pot <- function(dat, name = "", tail = c("max", "min"),
                        thresh_prob = 0.95,
                        time_units  = "365.25/year") {
  tail <- match.arg(tail)  # ensure that 'tail' is either "max" or "min" and pick that value
  
  if (tail == "max") {
    # high water levels
    # for maxima we use the raw water level values as they are
    x <- dat$value
    units_lbl <- "cm (raw values)"         # label for the units to show up in model summaries/plots
  } else {
    # low water levels -> flip sign so low becomes high
    # for minima we multiply by -1, so very low levels become large positive numbers
    # this lets us use the same "exceedance above a high threshold" framework
    x <- -dat$value
    units_lbl <- "cm (flipped, minima)"    # indicate that we are modeling -value instead of value
  }
  
  # threshold = upper quantile of the series we fit (x or -x)
  # we choose the threshold as the 'thresh_prob' quantile (e.g. 0.95 = 95% quantile),
  # i.e. only the top 5% of values (extremes) will be used to fit the GPD
  u <- quantile(x, thresh_prob, na.rm = TRUE)
  
  # fit a GPD (type="GP") to exceedances over u using extRemes::fevd
  fit <- fevd(x,
              threshold  = u,           # high threshold for POT
              type       = "GP",        # Generalized Pareto distribution
              units      = units_lbl,   # just a label, used in output
              time.units = time_units)  # how many observations per "year" (for return periods)
  
  # print some info about the fitted model to the console
  cat("\n=============================\n")
  cat("GPD fit for", name, "tail:", tail, "\n")
  cat("Threshold quantile prob:", thresh_prob, "\n")
  cat("Numeric threshold u:", u, "\n")
  print(fit)   # shows parameter estimates, log-likelihood, AIC, etc.
  
  # return a small list containing:
  # - the quantile level used for threshold (thresh_prob),
  # - the numeric threshold value (threshold),
  # - the actual fevd fit object (fit) that we will use later
  invisible(list(thresh_prob = thresh_prob,
                 threshold   = u,
                 fit         = fit))
}


## --------------------- DUNA 2002–2013 -----------------------------------
# Fit GPD POT models to Duna data for the first period (2002–2013)
# First for maxima (floods)
gpd_duna_02_13_max <- fit_gpd_pot(
  river_data_duna_02_13,
  name        = "Duna 2002–2013",
  tail        = "max",
  thresh_prob = 0.95   # 95% quantile as threshold (top 5% extremes)
)

# Then for minima (low-water events), using flipped values -value
gpd_duna_02_13_min <- fit_gpd_pot(
  river_data_duna_02_13,
  name        = "Duna 2002–2013",
  tail        = "min",
  thresh_prob = 0.95
)


## --------------------- DUNA 2014–2024 -----------------------------------
# Same as above but for the second period of Duna (2014–2024)
gpd_duna_14_24_max <- fit_gpd_pot(
  river_data_duna_14_24,
  name        = "Duna 2014–2024",
  tail        = "max",
  thresh_prob = 0.95
)

gpd_duna_14_24_min <- fit_gpd_pot(
  river_data_duna_14_24,
  name        = "Duna 2014–2024",
  tail        = "min",
  thresh_prob = 0.95
)


## --------------------- TISZA 2002–2013 ----------------------------------
# Now do exactly the same GPD POT fitting for the Tisza river, first period
gpd_tisza_02_13_max <- fit_gpd_pot(
  river_data_tisza_02_13,
  name        = "Tisza 2002–2013",
  tail        = "max",
  thresh_prob = 0.95
)

gpd_tisza_02_13_min <- fit_gpd_pot(
  river_data_tisza_02_13,
  name        = "Tisza 2002–2013",
  tail        = "min",
  thresh_prob = 0.95
)


## --------------------- TISZA 2014–2024 ----------------------------------
# And Tisza, second period
gpd_tisza_14_24_max <- fit_gpd_pot(
  river_data_tisza_14_24,
  name        = "Tisza 2014–2024",
  tail        = "max",
  thresh_prob = 0.95
)

gpd_tisza_14_24_min <- fit_gpd_pot(
  river_data_tisza_14_24,
  name        = "Tisza 2014–2024",
  tail        = "min",
  thresh_prob = 0.95
)

## ---- GPD-based 5000-day return levels (dam heights) ----
# 5000-day return period in years (consistent with time.units = "365.25/year")
# here we convert "5000 days" into "years", because fevd/return.level uses years
rp_5000_days <- 5000 / 365.25
rp_5000_days
# ~ 13.69 years  (this is the return period: we want events that happen on average once every ~13.7 years)

# Duna - Esztergom
# For the Duna maxima models, compute the 5000-day return level (how high the water gets on average every 5000 days)
RL_duna_02_13_GPD  <- return.level(gpd_duna_02_13_max$fit,
                                   return.period = rp_5000_days)
RL_duna_14_24_GPD  <- return.level(gpd_duna_14_24_max$fit,
                                   return.period = rp_5000_days)

# Tisza - Szeged
# Same for Tisza maxima models
RL_tisza_02_13_GPD <- return.level(gpd_tisza_02_13_max$fit,
                                   return.period = rp_5000_days)
RL_tisza_14_24_GPD <- return.level(gpd_tisza_14_24_max$fit,
                                   return.period = rp_5000_days)

# 95% CIs for the 5000-day RLs
# For each of the four maxima fits, we also calculate 95% confidence intervals
# for the 5000-day return level to quantify uncertainty.
CI_duna_02_13_GPD  <- ci(gpd_duna_02_13_max$fit,
                         type = "return.level",
                         return.period = rp_5000_days)
CI_duna_14_24_GPD  <- ci(gpd_duna_14_24_max$fit,
                         type = "return.level",
                         return.period = rp_5000_days)
CI_tisza_02_13_GPD <- ci(gpd_tisza_02_13_max$fit,
                         type = "return.level",
                         return.period = rp_5000_days)
CI_tisza_14_24_GPD <- ci(gpd_tisza_14_24_max$fit,
                         type = "return.level",
                         return.period = rp_5000_days)

# Summarise the 5000-day return levels in a single small data frame for easier comparison
dam_summary <- data.frame(
  river   = c("Duna", "Duna", "Tisza", "Tisza"),          # which river
  period  = c("2002–2013", "2014–2024", "2002–2013", "2014–2024"),  # which time period
  method  = "GPD (POT maxima)",                           # method description (POT GPD on maxima)
  RL_5000 = c(
    as.numeric(RL_duna_02_13_GPD),                        # 5000-day RL for Duna 2002–2013
    as.numeric(RL_duna_14_24_GPD),                        # 5000-day RL for Duna 2014–2024
    as.numeric(RL_tisza_02_13_GPD),                       # 5000-day RL for Tisza 2002–2013
    as.numeric(RL_tisza_14_24_GPD)                        # 5000-day RL for Tisza 2014–2024
  )
  # you can add CI columns if you like, extracting from CI_* objects
)

dam_summary   # look at the table: one row per river-period, with the suggested dam height for 5000-day risk


# comparison:
## -----------------------------------------------------------------------------
## 1) Parameter comparison between periods (scale & shape, with 5% tests)
## -----------------------------------------------------------------------------
# gpd_duna_02_13_max$fit$results$par       # example how to manually inspect parameters
# parcov.fevd(gpd_duna_02_13_max$fit)      # example how to get covariance matrix

# This function compares the GPD parameters (scale and shape) between two fits:
# typically fit1 for 2002–2013, fit2 for 2014–2024.
# It calculates Z-scores and p-values to see whether the change is significant at 5%.
compare_params <- function(fit1, fit2, label1, label2) {
  p1 <- fit1$results$par           # parameter vector (scale, shape) for first fit
  p2 <- fit2$results$par           # parameter vector (scale, shape) for second fit
  
  V1 <- parcov.fevd(fit1)          # covariance matrix of parameters for first fit
  V2 <- parcov.fevd(fit2)          # covariance matrix of parameters for second fit
  
  # basic sanity checks
  if (length(p1) == 0 || length(p2) == 0) {
    cat("\n=====================================\n")
    cat("Parameter comparison:", label1, "vs", label2, "\n")
    cat("ERROR: one of the fits has no parameters.\n")
    print(p1); print(p2)
    return(invisible(NULL))
  }
  
  # keep only parameters that exist in both fits (should be 'scale' and 'shape')
  common_par <- intersect(names(p1), names(p2))
  p1 <- p1[common_par]
  p2 <- p2[common_par]
  V1 <- V1[common_par, common_par, drop = FALSE]
  V2 <- V2[common_par, common_par, drop = FALSE]
  
  # Build a comparison table:
  # est_1: parameter estimate in first period,
  # est_2: parameter estimate in second period,
  # diff:  change (est_2 - est_1),
  # se_1, se_2: standard errors for each period,
  # se_diff: standard error of the difference assuming fits are independent.
  out <- data.frame(
    param   = common_par,
    est_1   = as.numeric(p1),
    est_2   = as.numeric(p2),
    diff    = as.numeric(p2 - p1),
    se_1    = sqrt(diag(V1)),
    se_2    = sqrt(diag(V2)),
    se_diff = sqrt(diag(V1) + diag(V2))  # assume independence of periods
  )
  
  # Wald Z-statistic for "no change" (diff = 0)
  out$z <- out$diff / out$se_diff
  # Two-sided p-value
  out$p <- 2 * pnorm(-abs(out$z))
  
  cat("\n=====================================\n")
  cat("Parameter comparison:", label1, "vs", label2, "\n")
  print(out)
  cat("p-values < 0.05 => significant change at 5% level.\n")
  
  invisible(out)  # return the table invisibly so we can save it if we want
}

## ---- RUN FOR EACH RIVER & TAIL ----

# Duna maxima: compare GPD parameters between 2002–2013 and 2014–2024
cmp_duna_max <- compare_params(
  gpd_duna_02_13_max$fit, gpd_duna_14_24_max$fit,
  "Duna max 2002–2013", "Duna max 2014–2024"
)

# Duna minima (fitted on -value): compare how low-water tail changed
cmp_duna_min <- compare_params(
  gpd_duna_02_13_min$fit, gpd_duna_14_24_min$fit,
  "Duna min 2002–2013", "Duna min 2014–2024"
)

# Tisza maxima
cmp_tisza_max <- compare_params(
  gpd_tisza_02_13_max$fit, gpd_tisza_14_24_max$fit,
  "Tisza max 2002–2013", "Tisza max 2014–2024"
)

# Tisza minima
cmp_tisza_min <- compare_params(
  gpd_tisza_02_13_min$fit, gpd_tisza_14_24_min$fit,
  "Tisza min 2002–2013", "Tisza min 2014–2024"
)



## -----------------------------------------------------------------------------
## 2) MOMENTS of exceedances for each river & period (GPD POT)
##    (mean and sd of exceedances above the threshold)
## -----------------------------------------------------------------------------

# This function calculates "moments" of the GPD exceedances:
# - sigma (scale)
# - xi (shape)
# - mean_exc: mean exceedance above threshold u (E[X - u | X > u])
# - sd_exc: standard deviation of exceedances
gpd_moments <- function(fit, label) {
  # GPD parameters from fevd object
  par_hat <- fit$results$par
  sigma   <- as.numeric(par_hat["scale"])   # scale parameter σ
  xi      <- as.numeric(par_hat["shape"])   # shape parameter ξ
  
  # mean exists if xi < 1, variance if xi < 0.5
  # formula for mean exceedance of GPD: σ / (1 - ξ)
  mean_exc <- if (xi < 1) {
    sigma / (1 - xi)
  } else NA_real_
  
  # formula for variance of exceedance: σ² / ((1 - ξ)² (1 - 2ξ))
  var_exc <- if (xi < 0.5) {
    sigma^2 / ((1 - xi)^2 * (1 - 2 * xi))
  } else NA_real_
  
  # standard deviation is sqrt of variance (if it exists)
  sd_exc <- if (!is.na(var_exc)) sqrt(var_exc) else NA_real_
  
  out <- data.frame(
    label    = label,    # which river & period & tail
    sigma    = sigma,    # GPD scale
    xi       = xi,       # GPD shape
    mean_exc = mean_exc, # mean excess over threshold (in cm)
    sd_exc   = sd_exc    # standard deviation of excess (in cm)
  )
  
  print(out)            # print so we see the numbers immediately
  invisible(out)        # also return them so we can store if needed
}

# Moments for maxima in all river-period combinations
mom_duna_02_13_max  <- gpd_moments(gpd_duna_02_13_max$fit,  "Duna max 2002–2013")
mom_duna_14_24_max  <- gpd_moments(gpd_duna_14_24_max$fit,  "Duna max 2014–2024")
mom_tisza_02_13_max <- gpd_moments(gpd_tisza_02_13_max$fit, "Tisza max 2002–2013")
mom_tisza_14_24_max <- gpd_moments(gpd_tisza_14_24_max$fit, "Tisza max 2014–2024")

# you *can* also do it for minima (on the flipped scale, i.e. for -value):
# here the moments describe exceedances of -value above its threshold,
# which corresponds to extremely low original water levels.
mom_duna_02_13_min  <- gpd_moments(gpd_duna_02_13_min$fit,  "Duna min 2002–2013 (on -value)")
mom_duna_14_24_min  <- gpd_moments(gpd_duna_14_24_min$fit,  "Duna min 2014–2024 (on -value)")
mom_tisza_02_13_min <- gpd_moments(gpd_tisza_02_13_min$fit, "Tisza min 2002–2013 (on -value)")
mom_tisza_14_24_min <- gpd_moments(gpd_tisza_14_24_min$fit, "Tisza min 2014–2024 (on -value)")

## -----------------------------------------------------------------------------
## 3) QUANTILES via return levels (2, 5, 10, 20, 50 years + 5000 days)
## -----------------------------------------------------------------------------

# Return periods in years (because time.units = "365.25/year")
# We consider standard design horizons (2, 5, 10, 20, 50 years) plus the ~13.7-year (5000-day) period.
rp_vec <- c(2, 5, 10, 20, 50, rp_5000_days)

# For maxima: directly call return.level on the GPD fit
get_RL_table_max <- function(fit, label) {
  data.frame(
    label = label,          # which river & period
    RP    = rp_vec,         # return period in years
    RL    = sapply(rp_vec, function(T)
      as.numeric(return.level(fit, return.period = T)))  # return level (cm) for each RP
  )
}

# For minima, we fitted GPD to -value, so the return levels
# we get are on the -value scale. We multiply by -1 to get back to original cm.
get_RL_table_min <- function(fit, label) {
  RL_flip <- sapply(rp_vec, function(T)
    as.numeric(return.level(fit, return.period = T)))
  data.frame(
    label = label,
    RP    = rp_vec,
    RL    = -RL_flip   # more negative RL (after flipping back) => more extreme low water
  )
}

## DUNA
# Put together return levels for Duna maxima:
RL_duna_max <- rbind(
  get_RL_table_max(gpd_duna_02_13_max$fit, "Duna max 2002–2013"),
  get_RL_table_max(gpd_duna_14_24_max$fit, "Duna max 2014–2024")
)

# And for Duna minima:
RL_duna_min <- rbind(
  get_RL_table_min(gpd_duna_02_13_min$fit, "Duna min 2002–2013"),
  get_RL_table_min(gpd_duna_14_24_min$fit, "Duna min 2014–2024")
)

## TISZA
# The same for Tisza maxima:
RL_tisza_max <- rbind(
  get_RL_table_max(gpd_tisza_02_13_max$fit, "Tisza max 2002–2013"),
  get_RL_table_max(gpd_tisza_14_24_max$fit, "Tisza max 2014–2024")
)

# And Tisza minima:
RL_tisza_min <- rbind(
  get_RL_table_min(gpd_tisza_02_13_min$fit, "Tisza min 2002–2013"),
  get_RL_table_min(gpd_tisza_14_24_min$fit, "Tisza min 2014–2024")
)

# Print all RL tables so we can visually compare how high/low extremes differ
# between periods for each river and tail.
# RP = Return Period “How often, on average, does such an extreme happen?”
#RL = Return Level “How high is the event that happens with that frequency?”
# RL = 64.33082 → return level → a water level of about 64 cm that the river is
# expected to drop down to (or below) on average once every 2 years in that period.
RL_duna_max
RL_duna_min
RL_tisza_max
RL_tisza_min
