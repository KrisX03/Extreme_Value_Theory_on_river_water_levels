# Computational Statistics
# Extreme Value Theory on river water levels

# Members of the team:
# Kristóf Andrási
# Gellért Banai
# Hunor Kuti
# Ákos Virág

# analysis GPD

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


# first we started to analyse one database to see how we should decide the threshold
# starting with river_data_duna_02_13

nrow(river_data_duna_02_13[river_data_duna_02_13$value>700,])
nrow(river_data_duna_02_13[river_data_duna_02_13$value>650,])
# there is only 26 values exceed 650cm and only 16 exceed 700 out of 4383 times
# Therefore to ensure we have enough data and to more easily interpret
# the mean residual life (mrl) plot, we will restrict the range to 650.

par(mar = c(4, 4, 2, 1))
threshrange.plot(
  river_data_duna_02_13$value,
  r    = c(300, 650),  
  nint = 20
)
# these are the fitted (reparameterized) scale and shape parameters over
# a range of equally spaced thresholds
# from this plot, we should choose the threshold where the estimates become roughly
# stable (they only change a little within their uncertainty), but the
# threshold is still low enough so that we keep as many data points as possible
# in the extreme value fit

dev.off()
mrlplot(river_data_duna_02_13$value,
        xlim = c(300, 650))  
# here the idea is to choose a threshold where the graph looks roughly like a
# straight line (only small random wiggles) and keeps this shape as the threshold gets higher

# based on these 2 plot it's difficult to say 1 good threshold, after long discussion we picked
# one around 410-430, and it turned out it's basically almost the 95th percentile
quantile(river_data_duna_02_13$value, 0.95)

# for simplification we picked the 95th percentile everywhere so we dont have make so 
# many difficult decisions


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

# Interpretation:
# - We model Duna daily *high* water levels in 2002–2013 using a GPD above 404 cm.
# - scale ≈ 111 cm: the typical size of exceedances over 404 cm is around 100 cm (quite large floods).
# - shape ≈ -0.08 (slightly negative): the upper tail is light and slightly bounded,
#   meaning extremely huge floods become less and less likely.
# - This is our baseline model for old-period Duna maxima.


# Then for minima (low-water events), using flipped values -value
gpd_duna_02_13_min <- fit_gpd_pot(
  river_data_duna_02_13,
  name        = "Duna 2002–2013",
  tail        = "min",
  thresh_prob = 0.95
)

# Interpretation:
# - We model *low* water levels by applying GPD to -value above threshold -51,
#   which corresponds to original water levels below about +51 cm.
# - scale ≈ 27: typical size of low-water exceedances (on -value scale) is smaller than for floods.
# - shape ≈ -0.49: strongly negative shape → very bounded tail for low levels:
#   there is a lower limit to how extreme the low water can get (harder to get extremely low).


## --------------------- DUNA 2014–2024 -----------------------------------
# Same as above but for the second period of Duna (2014–2024)
gpd_duna_14_24_max <- fit_gpd_pot(
  river_data_duna_14_24,
  name        = "Duna 2014–2024",
  tail        = "max",
  thresh_prob = 0.95
)

# Interpretation:
# - In the more recent period, the 95% threshold for Duna is *lower* (344 cm vs 404 cm),
#   meaning fewer very high water levels overall.
# - scale dropped from ~111 to ~70 → typical flood exceedances above threshold are much smaller now.
# - shape moved slightly positive (~0.03), but close to 0 → the tail is still not extremely heavy.
# - Overall: recent Duna floods seem less extreme (smaller scale) than in 2002–2013.

gpd_duna_14_24_min <- fit_gpd_pot(
  river_data_duna_14_24,
  name        = "Duna 2014–2024",
  tail        = "min",
  thresh_prob = 0.95
)

# Interpretation:
# - Low-water model for Duna in 2014–2024.
# - scale ≈ 24 vs 27 previously → typical low-water extremes shrank a bit.
# - shape ≈ -0.36 vs -0.49 → tail is still bounded but a bit “less bounded” than before.
# - Overall: extremely low levels are still limited, but perhaps slightly less extreme structure than in 2002–2013.


## --------------------- TISZA 2002–2013 ----------------------------------
# Now do exactly the same GPD POT fitting for the Tisza river, first period
gpd_tisza_02_13_max <- fit_gpd_pot(
  river_data_tisza_02_13,
  name        = "Tisza 2002–2013",
  tail        = "max",
  thresh_prob = 0.95
)

# Interpretation:
# - Very high threshold: Tisza floods in 2002–2013 were really big (662+ cm).
# - scale ≈ 120 cm: very large typical exceedances above that threshold → extremely big floods.
# - shape ≈ -0.19: moderately negative shape → upper tail is bounded; extremely huge levels are rare.
# - This is our baseline model for old-period Tisza floods.

gpd_tisza_02_13_min <- fit_gpd_pot(
  river_data_tisza_02_13,
  name        = "Tisza 2002–2013",
  tail        = "min",
  thresh_prob = 0.95
)

# Interpretation:
# - Low-water model for Tisza 2002–2013.
# - scale ≈ 12.7: typical exceedances over low-water threshold are modest in size.
# - shape ≈ -0.52: strongly negative → very bounded lower tail (hard limit on how low it goes).
# - Tisza low water in this period is extreme but quite bounded.


## --------------------- TISZA 2014–2024 ----------------------------------
# And Tisza, second period
gpd_tisza_14_24_max <- fit_gpd_pot(
  river_data_tisza_14_24,
  name        = "Tisza 2014–2024",
  tail        = "max",
  thresh_prob = 0.95
)

# Interpretation:
# - Recent Tisza floods: threshold dropped from ~663 cm to ~501 cm → overall lower flood levels.
# - scale dropped from ~120 to ~64 → typical flood exceedances are about half as large as before.
# - shape changed from -0.19 to about -0.54: tail became *much more strongly bounded*,
#   meaning the chance of truly massive floods is much lower now according to this model.

gpd_tisza_14_24_min <- fit_gpd_pot(
  river_data_tisza_14_24,
  name        = "Tisza 2014–2024",
  tail        = "min",
  thresh_prob = 0.95
)

# Interpretation:
# - Recent Tisza low-water extremes.
# - scale decreased from ~12.7 to ~8.0: typical low-water deviations are smaller now.
# - shape from -0.52 to -0.37: tail still bounded, but less strongly than before; the pattern of low extremes changed.


## ---- GPD-based 5000-day return levels (dam heights) ----
# 5000-day return period in years (consistent with time.units = "365.25/year")
# here we convert "5000 days" into "years", because fevd/return.level uses years
rp_5000_days <- 5000 / 365.25
rp_5000_days


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
# We now have the “5000-day design flood height” for each river and period.


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

# These CIs tell us how uncertain the design height is. If CIs for the two periods do not
# overlap much, that suggests a clear change in flood risk between 2002–2013 and 2014–2024.
# im not sure about this overlapping interpretation - Kristóf

# Summarise the 5000-day return levels in a single small data frame for easier comparison
dam_summary <- data.frame(
  river   = c("Duna", "Duna", "Tisza", "Tisza"),          # which river
  period  = c("2002–2013", "2014–2024", "2002–2013", "2014–2024"),  # which time period
  #method  = "GPD (POT maxima)",                           # method description (POT GPD on maxima)
  RL_5000 = c(
    as.numeric(RL_duna_02_13_GPD),                        # 5000-day RL for Duna 2002–2013
    as.numeric(RL_duna_14_24_GPD),                        # 5000-day RL for Duna 2014–2024
    as.numeric(RL_tisza_02_13_GPD),                       # 5000-day RL for Tisza 2002–2013
    as.numeric(RL_tisza_14_24_GPD)                        # 5000-day RL for Tisza 2014–2024
  )
)

dam_summary   # look at the table: one row per river-period, with the suggested dam height for 5000-day risk

# with CI:
dam_summary_CI <- data.frame(
  river     = c("Duna", "Duna", "Tisza", "Tisza"),              # which river
  period    = c("2002–2013", "2014–2024", "2002–2013", "2014–2024"),
  RL_5000   = c(
    as.numeric(RL_duna_02_13_GPD),
    as.numeric(RL_duna_14_24_GPD),
    as.numeric(RL_tisza_02_13_GPD),
    as.numeric(RL_tisza_14_24_GPD)
  ),
  RL_5000_L = c(  # lower CI bound
    as.numeric(CI_duna_02_13_GPD[1]),
    as.numeric(CI_duna_14_24_GPD[1]),
    as.numeric(CI_tisza_02_13_GPD[1]),
    as.numeric(CI_tisza_14_24_GPD[1])
  ),
  RL_5000_U = c(  # upper CI bound
    as.numeric(CI_duna_02_13_GPD[2]),
    as.numeric(CI_duna_14_24_GPD[2]),
    as.numeric(CI_tisza_02_13_GPD[2]),
    as.numeric(CI_tisza_14_24_GPD[2])
  )
)

dam_summary_CI

# Interpretation:
# - Duna:
#   * 2002–2013: a 5000-day flood level is ~904 cm.
#   * 2014–2024: a 5000-day flood level is ~758 cm.
#   → Suggested dam height for this risk level would be much higher in the earlier period.
#     Flood risk (in terms of extreme height) seems to have decreased in 2014–2024.
#
# - Tisza:
#   * 2002–2013: RL_5000 ≈ 1068 cm (!! very high extreme floods).
#   * 2014–2024: RL_5000 ≈ 613 cm, which is dramatically lower.
#   → Huge reduction in modeled extreme flood height compared to the earlier period.
#
# Overall: according to the GPD fits, both rivers show smaller extreme flood levels in 2014–2024
# than in 2002–2013, with a particularly large drop for the Tisza.



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
# once again: WALD TEST!
# Duna maxima: compare GPD parameters between 2002–2013 and 2014–2024
cmp_duna_max <- compare_params(
  gpd_duna_02_13_max$fit, gpd_duna_14_24_max$fit,
  "Duna max 2002–2013", "Duna max 2014–2024"
)

# Interpretation:
# - For Duna flood extremes, the *scale* parameter decreased significantly (p < 0.01),
#   meaning that typical size of extreme exceedances is clearly smaller in 2014–2024.
# - The *shape* parameter change is NOT significant (p ≈ 0.35),
#   so the overall tail heaviness/boundedness is not statistically different at 5%.
# - In short: for Duna maxima, the scale shrank significantly, tail shape stayed similar.

# further explanation: se_1: How uncertain is the estimate 111.44 (scale) in the first period?
# Then we combine them: se_diff = sqrt(se_1^2 + se_2^2) assuming the two fits are independent.
#	z = diff / se_diff

# Duna minima (fitted on -value): compare how low-water tail changed
cmp_duna_min <- compare_params(
  gpd_duna_02_13_min$fit, gpd_duna_14_24_min$fit,
  "Duna min 2002–2013", "Duna min 2014–2024"
)

# Interpretation:
# - For Duna minima, *neither* scale nor shape changes are significant at 5% (p > 0.05).
# - So we cannot statistically claim that low-water tail behaviour changed between periods.
# - There is a hint of shape change (p≈0.07), but it’s not below 0.05.


# Tisza maxima
cmp_tisza_max <- compare_params(
  gpd_tisza_02_13_max$fit, gpd_tisza_14_24_max$fit,
  "Tisza max 2002–2013", "Tisza max 2014–2024"
)

# Interpretation:
# - For Tisza flood extremes, *both* scale and shape changed significantly (p ≪ 0.01).
# - Scale dropped strongly: typical exceedances are much smaller now.
# - Shape became much more negative → tail became much more bounded,
#   i.e. the risk of extremely huge floods decreased a lot in 2014–2024.
# - Statistically, the flood distribution for Tisza changed a lot between the two periods.


# Tisza minima
cmp_tisza_min <- compare_params(
  gpd_tisza_02_13_min$fit, gpd_tisza_14_24_min$fit,
  "Tisza min 2002–2013", "Tisza min 2014–2024"
)

# Interpretation:
# - For Tisza low-water extremes, both scale and shape changes are significant at 5%.
# - scale decreased → typical low-water deviations (on -value scale) are smaller now.
# - shape moved closer to 0 (less negative) → the lower tail is still bounded, but less strongly.
# - So the pattern of extreme low levels of Tisza also changed meaningfully over time.



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
  # originally its mu + σ / (1 - ξ), but we model exceedances over a threshold u:
  # Y = X - u | X > u and for Y the mean is σ / (1 - ξ).
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

# Interpretation:
# - For Duna maxima 2002–2013, the average exceedance above the 404 cm threshold is ~103 cm,
#   with standard deviation ~96 cm. So extreme floods are large and quite variable.

mom_duna_14_24_max  <- gpd_moments(gpd_duna_14_24_max$fit,  "Duna max 2014–2024")

# Interpretation:
# - In 2014–2024, average exceedance is ~72 cm (vs 103 cm before) and variability is smaller (~74 vs 96).
# - So “typical” Duna extreme floods became smaller and slightly less variable.

mom_tisza_02_13_max <- gpd_moments(gpd_tisza_02_13_max$fit, "Tisza max 2002–2013")

# Interpretation:
# - For Tisza maxima 2002–2013, average exceedance above ~663 cm is ~100 cm with SD ~85 cm.
# - So big and quite variable floods in the early period.

mom_tisza_14_24_max <- gpd_moments(gpd_tisza_14_24_max$fit, "Tisza max 2014–2024")

# Interpretation:
# - Recent Tisza floods: mean exceedance dropped to ~41 cm (from 100 cm),
#   and SD shrank to ~29 cm (from 85 cm).
# - So not only are extremes smaller, but they are also much less variable.


# we *can* also do it for minima (on the flipped scale, i.e. for -value):
# here the moments describe exceedances of -value above its threshold,
# which corresponds to extremely low original water levels.
mom_duna_02_13_min  <- gpd_moments(gpd_duna_02_13_min$fit,  "Duna min 2002–2013 (on -value)")

# Interpretation (on -value scale):
# - For Duna low water 2002–2013, the average “depth” of low-water extremes below the threshold
#   corresponds to about 18 units with SD ~13 units. Larger numbers here mean more extreme lows.

mom_duna_14_24_min  <- gpd_moments(gpd_duna_14_24_min$fit,  "Duna min 2014–2024 (on -value)")

# Interpretation:
# - For Duna low water, mean and SD on the -value scale are similar (~18, ~14) between periods.
# - This fits our earlier parameter comparison: no strong change in low-water tail behaviour.

mom_tisza_02_13_min <- gpd_moments(gpd_tisza_02_13_min$fit, "Tisza min 2002–2013 (on -value)")

# Interpretation:
# - Old Tisza low water: mean exceedance ~8.35, SD ~5.86 (on -value scale),
#   so low-water extremes have moderate size and variation.

mom_tisza_14_24_min <- gpd_moments(gpd_tisza_14_24_min$fit, "Tisza min 2014–2024 (on -value)")

# Interpretation:
# - Recent Tisza low water: mean exceedance dropped to ~5.8 and SD to ~4.4.
# - So low-water extremes for Tisza also became milder and less variable.

moments_all <- do.call(
  rbind,
  list(
    mom_duna_02_13_max,
    mom_duna_14_24_max,
    mom_tisza_02_13_max,
    mom_tisza_14_24_max,
    mom_duna_02_13_min,
    mom_duna_14_24_min,
    mom_tisza_02_13_min,
    mom_tisza_14_24_min
  )
)

# optional: clean rownames
rownames(moments_all) <- NULL

moments_all

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

# quick note:
# RP = Return Period “How often, on average, does such an extreme happen?”
# RL = Return Level “How high is the event that happens with that frequency?”

RL_duna_max

# Interpretation:
# - For every RP (2, 5, 10, 20, 50 years and 5000 days), the return level is higher in 2002–2013 than in 2014–2024.
# - Example: 10-year flood:
#   * old period: ~881 cm,
#   * new period: ~733 cm.
# - So the Duna’s extreme flood quantiles systematically decreased in the later period.


RL_duna_min

# Interpretation (!!! smaller / more negative RL = more extreme low water):
# - In 2002–2013, RLs cluster around 0 cm, i.e. low-water extremes are not extremely negative.
# - In 2014–2024, RLs are much more negative (e.g. 10-year low ≈ -15.2 cm).
# - This suggests that Duna low-water events might actually have become more severe (lower levels)
#   in the more recent period, while floods got smaller.


RL_tisza_max

# Interpretation:
# - Early-period Tisza floods are huge: e.g. 10-year RL ≈ 1055 cm.
# - In 2014–2024, all RLs are around 602–616 cm, almost flat with RP.
# - This matches the very negative shape in the recent fit: flood extremes are much less “explosive”,
#   and the difference between 2-year and 50-year floods is small.


RL_tisza_min

# Interpretation:
# - Here, RLs are *higher* in 2014–2024 (e.g. 2-year low water: 64.3 → 67.3 cm).
# - Higher RL for minima means the river does *not* drop as low as before.
# - So for Tisza, low-water extremes became *less* severe in the recent period,
#   while floods also became much smaller and more bounded.

RL_duna_both <- rbind(
  cbind(RL_duna_max, tail = "max"),
  cbind(RL_duna_min, tail = "min")
)

# optional: reset rownames
rownames(RL_duna_both) <- NULL

RL_duna_both


RL_tisza_both <- rbind(
  cbind(RL_tisza_max, tail = "max"),
  cbind(RL_tisza_min, tail = "min")
)

rownames(RL_tisza_both) <- NULL

RL_tisza_both



# Summary for the whole script:
# - Duna:
#   * Flood extremes (maxima): smaller scale, lower RLs in 2014–2024 → less extreme floods.
#   * Low-water extremes (minima): RLs became more negative → possible worsening of drought/low levels.
#   * Significance: flood scale change is significant, low-water parameter changes are not (at 5%).
#
# - Tisza:
#   * Flood extremes: huge drop in scale and strong change in shape → much milder and more bounded floods now.
#   * Low-water extremes: both parameters and moments decreased → low-water events also became milder.
#   * All these changes are statistically significant at the 5% level for both maxima and minima.


