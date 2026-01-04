# Computational Statistics
# Extreme Value Theory on river water levels

# Members of the team:
# Kristóf Andrási
# Gellért Banai
# Hunor Kuti
# Ákos Virág

# analysis GEV

# ------------------------------------------------------------------------------
# importing neccesary libraries
# vector of required packages
pkgs <- c(
  "readr",   # for reading csv files
  "dplyr",   # for data manipulation (group_by, summarise, mutate, etc.)
  "extRemes" # for extreme value modelling (fevd, return.level, lr.test, ci)
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
# quick note: In theory, there should be no need to set the working directory manually, because when you open the file
# R sets it automatically. If this doesn’t work, use setwd to set the “codes” folder as the working directory.

#setwd("C:/Users/User/Desktop/Computational Statistics/Comp_stat_projekt")

# loading the dataframes
river_data <- read.csv("../data/river_data.csv")                 # full combined daily river data (both rivers, all years)
river_data_tisza_02_13 <- read.csv("../data/river_data_tisza_02_13.csv")  # Tisza daily data, 2002–2013 only
river_data_tisza_14_24 <- read.csv("../data/river_data_tisza_14_24.csv")  # Tisza daily data, 2014–2024 only
river_data_duna_02_13  <- read.csv("../data/river_data_duna_02_13.csv")   # Duna daily data, 2002–2013 only
river_data_duna_14_24  <- read.csv("../data/river_data_duna_14_24.csv")   # Duna daily data, 2014–2024 only

# Return period of interest:
# 5000 days ≈ 5000/365 years ≈ 13.7 years
# This will be used in all return.level() calls below.
T_years <- 5000 / 365  


# -----------------------------------------------------------------------------

############################################################
## TISZA 2014–2024
############################################################

# 1. Compute annual minima and maxima -----------------------------------------

minmax_tisza_14_24 <- river_data_tisza_14_24 %>% 
  group_by(year) %>% 
  summarise(
    min_value = min(value, na.rm = TRUE),  # annual minimum water level
    max_value = max(value, na.rm = TRUE)   # annual maximum water level
  )

############################################################
## 2. Fit GEV to annual maxima (Block Maxima)
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "tisza_gev_max_2014_2024_diagnostics.png")

fit_gev_max_tisza_14_24 <- fevd(
  x          = minmax_tisza_14_24$max_value,  # annual maxima
  type       = "GEV",                         # Generalized Extreme Value
  method     = "MLE",                         # maximum likelihood estimation
  time.units = "years"                        # 1 observation = 1 year
)

# Parameter estimates, standard errors, AIC/BIC, etc.
summary(fit_gev_max_tisza_14_24)

# Export clean 2x2 diagnostics plot
png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_max_tisza_14_24, type="qq",      main="", sub="")
plot(fit_gev_max_tisza_14_24, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_max_tisza_14_24, type="density", main="", sub="")
plot(fit_gev_max_tisza_14_24, type="rl",      main="", sub="")

mtext("Tisza River – annual maxima (2014–2024), GEV",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

# Optional: message where it was saved
message("Saved: ", out_file)


############################################################
## 3. Fit GEV to annual minima (via -min trick)
############################################################

############################################################
## TISZA – annual minima (2014–2024): GEV fit on -min + clean diagnostics export (EN)
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "tisza_gev_min_2014_2024_diagnostics.png")

# Fit: for minima model Y = -min_value as maxima (since min X = -max(-X))
fit_gev_min_tisza_14_24 <- fevd(
  x          = -minmax_tisza_14_24$min_value,  # negative annual minima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

# Console output
print(summary(fit_gev_min_tisza_14_24))

# Save the 2x2 diagnostic panel
png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow = c(2, 2),
    mar   = c(4.2, 4.2, 2.0, 1.0),
    oma   = c(0, 0, 3.0, 0),
    cex.lab  = 0.90,
    cex.axis = 0.85)

plot(fit_gev_min_tisza_14_24, type = "qq",      main = "",   sub = "")
plot(fit_gev_min_tisza_14_24, type = "qq2",     main = NULL,
     xlab = "Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_min_tisza_14_24, type = "density", main = "",   sub = "")
plot(fit_gev_min_tisza_14_24, type = "rl",      main = "",   sub = "")

mtext("Tisza River – annual minima (2014–2024), GEV (fit to -min)",
      outer = TRUE, line = 1.2, cex = 1.0)

dev.off()
par(op)

message("Saved: ", out_file)


############################################################
## 4. 5000-day (~13.7-year) return levels
############################################################

# 4a) Return level for annual maxima (e.g. design high water level)
rl_max_5000d_tisza_14_24 <- return.level(
  fit_gev_max_tisza_14_24,
  return.period = T_years
)
rl_max_5000d_tisza_14_24

# 4b) Return level for annual minima
# First we get the return level on the Y = -min_value scale,
# then we negate it back to obtain the level on the original scale.
rl_min_5000d_neg_tisza_14_24 <- return.level(
  fit_gev_min_tisza_14_24,
  return.period = T_years
)

rl_min_5000d_tisza_14_24 <- -rl_min_5000d_neg_tisza_14_24
rl_min_5000d_tisza_14_24


# -----------------------------------------------------------------------------


############################################################
## TISZA 2002–2013
############################################################

# 1. Annual min–max -----------------------------------------------------------

minmax_tisza_02_13 <- river_data_tisza_02_13 %>% 
  group_by(year) %>% 
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  )

############################################################
## 2. GEV fit for annual maxima
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "tisza_gev_max_2002_2013_diagnostics.png")

fit_gev_max_tisza_02_13 <- fevd(
  x          = minmax_tisza_02_13$max_value,  # annual maxima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

summary(fit_gev_max_tisza_02_13)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_max_tisza_02_13, type="qq",      main="", sub="")
plot(fit_gev_max_tisza_02_13, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_max_tisza_02_13, type="density", main="", sub="")
plot(fit_gev_max_tisza_02_13, type="rl",      main="", sub="")

mtext("Tisza River – annual maxima (2002–2013), GEV",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)


############################################################
## 3. Gumbel fit for annual minima (instead of full GEV)
############################################################

# NOTE: Full 3-parameter GEV for minima was numerically unstable here,
# so we switch to the simpler Gumbel model (shape ξ = 0).
# We still work on Y = -min_value to treat minima as maxima.
out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "tisza_gumbel_min_2002_2013_diagnostics.png")

fit_gumbel_min_tisza_02_13 <- fevd(
  x          = -minmax_tisza_02_13$min_value,  # negative annual minima
  type       = "Gumbel",                       # GEV with xi = 0
  method     = "MLE",
  time.units = "years"
)

summary(fit_gumbel_min_tisza_02_13)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gumbel_min_tisza_02_13, type="qq",      main="", sub="")
plot(fit_gumbel_min_tisza_02_13, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gumbel_min_tisza_02_13, type="density", main="", sub="")
plot(fit_gumbel_min_tisza_02_13, type="rl",      main="", sub="")

mtext("Tisza River – annual minima (2002–2013), Gumbel (fit to -min)",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)




############################################################
## 4. 5000-day (~13.7-year) return levels
############################################################

# 4a) Return level for annual maxima
rl_max_5000d_tisza_02_13 <- return.level(
  fit_gev_max_tisza_02_13,
  return.period = T_years
)
rl_max_5000d_tisza_02_13

# 4b) Return level for annual minima (Gumbel fit, on -min scale)
rl_min_5000d_neg_tisza_02_13 <- return.level(
  fit_gumbel_min_tisza_02_13,
  return.period = T_years
)

# Back-transform to original scale
rl_min_5000d_tisza_02_13 <- -rl_min_5000d_neg_tisza_02_13
rl_min_5000d_tisza_02_13


# -----------------------------------------------------------------------------


############################################################
## DUNA 2002–2013
############################################################

# 1. Annual min–max -----------------------------------------------------------

minmax_duna_02_13 <- river_data_duna_02_13 %>% 
  group_by(year) %>% 
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  )

############################################################
## 2. GEV fit for annual maxima
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "danube_gev_max_2002_2013_diagnostics.png")

fit_gev_max_duna_02_13 <- fevd(
  x          = minmax_duna_02_13$max_value,  # annual maxima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

summary(fit_gev_max_duna_02_13)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_max_duna_02_13, type="qq",      main="", sub="")
plot(fit_gev_max_duna_02_13, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_max_duna_02_13, type="density", main="", sub="")
plot(fit_gev_max_duna_02_13, type="rl",      main="", sub="")

mtext("Danube River – annual maxima (2002–2013), GEV",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)


############################################################
## 3. GEV fit for annual minima (via -min trick)
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "danube_gev_min_2002_2013_diagnostics.png")

fit_gev_min_duna_02_13 <- fevd(
  x          = -minmax_duna_02_13$min_value,  # negative annual minima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

summary(fit_gev_min_duna_02_13)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_min_duna_02_13, type="qq",      main="", sub="")
plot(fit_gev_min_duna_02_13, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_min_duna_02_13, type="density", main="", sub="")
plot(fit_gev_min_duna_02_13, type="rl",      main="", sub="")

mtext("Danube River – annual minima (2002–2013), GEV (fit to -min)",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)

############################################################
## 4. 5000-day return levels
############################################################

# 4a) Maxima
rl_max_5000d_duna_02_13 <- return.level(
  fit_gev_max_duna_02_13,
  return.period = T_years
)
rl_max_5000d_duna_02_13

# 4b) Minima (back-transform from -min scale)
rl_min_5000d_neg_duna_02_13 <- return.level(
  fit_gev_min_duna_02_13,
  return.period = T_years
)

rl_min_5000d_duna_02_13 <- -rl_min_5000d_neg_duna_02_13
rl_min_5000d_duna_02_13


# -----------------------------------------------------------------------------


############################################################
## DUNA 2014–2024
############################################################

# 1. Annual min–max -----------------------------------------------------------

minmax_duna_14_24 <- river_data_duna_14_24 %>% 
  group_by(year) %>% 
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  )

############################################################
## 2. GEV fit for annual maxima
############################################################
out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "danube_gev_max_2014_2024_diagnostics.png")

fit_gev_max_duna_14_24 <- fevd(
  x          = minmax_duna_14_24$max_value,  # annual maxima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

summary(fit_gev_max_duna_14_24)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_max_duna_14_24, type="qq",      main="", sub="")
plot(fit_gev_max_duna_14_24, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_max_duna_14_24, type="density", main="", sub="")
plot(fit_gev_max_duna_14_24, type="rl",      main="", sub="")

mtext("Danube River – annual maxima (2014–2024), GEV",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)


############################################################
## 3. GEV fit for annual minima (via -min trick)
############################################################

out_dir <- normalizePath(file.path(getwd(), "..", "figures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "danube_gev_min_2014_2024_diagnostics.png")

fit_gev_min_duna_14_24 <- fevd(
  x          = -minmax_duna_14_24$min_value,  # negative annual minima
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

summary(fit_gev_min_duna_14_24)

png(filename = out_file, width = 2200, height = 1800, res = 300)

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

par(mfrow=c(2,2),
    mar=c(4.2,4.2,2.0,1.0),
    oma=c(0,0,3.0,0),
    cex.lab=0.90, cex.axis=0.85)

plot(fit_gev_min_duna_14_24, type="qq",      main="", sub="")
plot(fit_gev_min_duna_14_24, type="qq2",     main=NULL,
     xlab="Empirical Quantiles")  # keep only "Empirical Quantiles" under panel 2
plot(fit_gev_min_duna_14_24, type="density", main="", sub="")
plot(fit_gev_min_duna_14_24, type="rl",      main="", sub="")

mtext("Danube River – annual minima (2014–2024), GEV (fit to -min)",
      outer=TRUE, line=1.2, cex=1.0)

dev.off()
par(op)

message("Saved: ", out_file)


############################################################
## 4. 5000-day return levels
############################################################

# 4a) Maxima
rl_max_5000d_duna_14_24 <- return.level(
  fit_gev_max_duna_14_24,
  return.period = T_years
)
rl_max_5000d_duna_14_24

# 4b) Minima (back-transform)
rl_min_5000d_neg_duna_14_24 <- return.level(
  fit_gev_min_duna_14_24,
  return.period = T_years
)

rl_min_5000d_duna_14_24 <- -rl_min_5000d_neg_duna_14_24
rl_min_5000d_duna_14_24


# -----------------------------------------------------------------------------


# Add period labels to each min–max data frame for later comparison ----------
minmax_tisza_02_13 <- minmax_tisza_02_13 %>%
  mutate(period = "early")

minmax_tisza_14_24 <- minmax_tisza_14_24 %>%
  mutate(period = "late")

minmax_duna_02_13 <- minmax_duna_02_13 %>%
  mutate(period = "early")

minmax_duna_14_24 <- minmax_duna_14_24 %>%
  mutate(period = "late")


# -----------------------------------------------------------------------------


############################################################
## TISZA – MAXIMA: 2002–2013 vs 2014–2024
############################################################

# Combine both periods into one data frame with a period indicator
tisza_max_all <- dplyr::bind_rows(
  minmax_tisza_02_13 %>% dplyr::select(year, max_value, period),
  minmax_tisza_14_24 %>% dplyr::select(year, max_value, period)
) %>%
  dplyr::mutate(period = factor(period, levels = c("early", "late")))


# H0: stationary GEV (μ, σ, ξ are the same for both periods)
fit_const_tisza_max <- fevd(
  x          = tisza_max_all$max_value,
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

# H1a: only location depends on period (μ_early, μ_late, common σ, ξ)
fit_loc_tisza_max <- fevd(
  x            = tisza_max_all$max_value,
  data         = tisza_max_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,  # μ = β0 + β1 * I(period == "late")
  scale.fun    = ~ 1,       # common scale
  shape.fun    = ~ 1        # common shape
)

# H1b: location and scale both depend on period
#      (μ_early, μ_late, σ_early, σ_late, common ξ)
fit_ns_tisza_max <- fevd(
  x            = tisza_max_all$max_value,
  data         = tisza_max_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ period,
  shape.fun    = ~ 1
)

# Likelihood-ratio tests: compare nested models
# Is a period effect on location significant?
lr_loc_tisza_max <- lr.test(fit_const_tisza_max, fit_loc_tisza_max)
lr_loc_tisza_max

# Is a period effect on location AND scale jointly significant?
lr_ns_tisza_max <- lr.test(fit_const_tisza_max, fit_ns_tisza_max)
lr_ns_tisza_max


# -----------------------------------------------------------------------------


############################################################
## TISZA – MINIMA: 2002–2013 vs 2014–2024
############################################################

tisza_min_all <- dplyr::bind_rows(
  minmax_tisza_02_13 %>% dplyr::select(year, min_value, period),
  minmax_tisza_14_24 %>% dplyr::select(year, min_value, period)
) %>%
  dplyr::mutate(
    period  = factor(period, levels = c("early", "late")),
    neg_min = -min_value   # transform minima to maxima for GEV modelling
  )


# H0: stationary GEV for Y = -min_value
fit_const_tisza_min <- fevd(
  x          = tisza_min_all$neg_min,
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

# H1a: only location depends on period
fit_loc_tisza_min <- fevd(
  x            = tisza_min_all$neg_min,
  data         = tisza_min_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ 1,
  shape.fun    = ~ 1
)

# H1b: location and scale depend on period
fit_ns_tisza_min <- fevd(
  x            = tisza_min_all$neg_min,
  data         = tisza_min_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ period,
  shape.fun    = ~ 1
)

# LR tests for changes in the minimum distribution
lr_loc_tisza_min <- lr.test(fit_const_tisza_min, fit_loc_tisza_min)
lr_loc_tisza_min

lr_ns_tisza_min <- lr.test(fit_const_tisza_min, fit_ns_tisza_min)
lr_ns_tisza_min


# -----------------------------------------------------------------------------


############################################################
## DUNA – MAXIMA: 2002–2013 vs 2014–2024
############################################################

duna_max_all <- dplyr::bind_rows(
  minmax_duna_02_13 %>% dplyr::select(year, max_value, period),
  minmax_duna_14_24 %>% dplyr::select(year, max_value, period)
) %>%
  dplyr::mutate(
    period = factor(period, levels = c("early", "late"))
  )


# H0: stationary GEV
fit_const_duna_max <- fevd(
  x          = duna_max_all$max_value,
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

# H1a: only location depends on period
fit_loc_duna_max <- fevd(
  x            = duna_max_all$max_value,
  data         = duna_max_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ 1,
  shape.fun    = ~ 1
)

# H1b: location and scale depend on period
fit_ns_duna_max <- fevd(
  x            = duna_max_all$max_value,
  data         = duna_max_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ period,
  shape.fun    = ~ 1
)

# LR tests
lr_loc_duna_max <- lr.test(fit_const_duna_max, fit_loc_duna_max)
lr_loc_duna_max

lr_ns_duna_max <- lr.test(fit_const_duna_max, fit_ns_duna_max)
lr_ns_duna_max


# -----------------------------------------------------------------------------


############################################################
## DUNA – MINIMA: 2002–2013 vs 2014–2024
############################################################

duna_min_all <- dplyr::bind_rows(
  minmax_duna_02_13 %>% dplyr::select(year, min_value, period),
  minmax_duna_14_24 %>% dplyr::select(year, min_value, period)
) %>%
  dplyr::mutate(
    period  = factor(period, levels = c("early", "late")),
    neg_min = -min_value  # transform minima to maxima
  )


# H0: stationary GEV for Y = -min_value
fit_const_duna_min <- fevd(
  x          = duna_min_all$neg_min,
  type       = "GEV",
  method     = "MLE",
  time.units = "years"
)

# H1a: only location depends on period
fit_loc_duna_min <- fevd(
  x            = duna_min_all$neg_min,
  data         = duna_min_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ 1,
  shape.fun    = ~ 1
)

# H1b: location and scale depend on period
fit_ns_duna_min <- fevd(
  x            = duna_min_all$neg_min,
  data         = duna_min_all,
  type         = "GEV",
  method       = "MLE",
  time.units   = "years",
  location.fun = ~ period,
  scale.fun    = ~ period,
  shape.fun    = ~ 1
)

# LR tests
lr_loc_duna_min <- lr.test(fit_const_duna_min, fit_loc_duna_min)
lr_loc_duna_min

lr_ns_duna_min <- lr.test(fit_const_duna_min, fit_ns_duna_min)
lr_ns_duna_min


# -----------------------------------------------------------------------------


############################################################
## CONFIDENCE INTERVALS FOR MODEL PARAMETERS
############################################################

## --- TISZA MAXIMA ---

ci_gev_max_tisza_02_13_param <- ci(
  fit_gev_max_tisza_02_13,
  type = "parameter"      # CIs for (μ, σ, ξ)
)
ci_gev_max_tisza_14_24_param <- ci(
  fit_gev_max_tisza_14_24,
  type = "parameter"
)

ci_gev_max_tisza_02_13_param
ci_gev_max_tisza_14_24_param


## --- TISZA MINIMA ---
## NOTE: These refer to the models for Y = -min_value.

ci_gev_min_tisza_02_13_param <- ci(
  fit_gumbel_min_tisza_02_13,
  type = "parameter"
)
ci_gev_min_tisza_14_24_param <- ci(
  fit_gev_min_tisza_14_24,
  type = "parameter"
)

ci_gev_min_tisza_02_13_param
ci_gev_min_tisza_14_24_param


## --- DUNA MAXIMA ---

ci_gev_max_duna_02_13_param <- ci(
  fit_gev_max_duna_02_13,
  type = "parameter"
)
ci_gev_max_duna_14_24_param <- ci(
  fit_gev_max_duna_14_24,
  type = "parameter"
)

ci_gev_max_duna_02_13_param
ci_gev_max_duna_14_24_param


## --- DUNA MINIMA ---

ci_gev_min_duna_02_13_param <- ci(
  fit_gev_min_duna_02_13,
  type = "parameter"
)
ci_gev_min_duna_14_24_param <- ci(
  fit_gev_min_duna_14_24,
  type = "parameter"
)

ci_gev_min_duna_02_13_param
ci_gev_min_duna_14_24_param


