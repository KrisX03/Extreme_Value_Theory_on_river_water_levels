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
  "extRemes"
)

# Install any that are missing
to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_install) > 0) {
  install.packages(to_install)
}

# loading them
invisible(lapply(pkgs, library, character.only = TRUE))
rm(to_install)
# ------------------------------------------------------------------------------
# loading the dataframes
river_data <- read.csv("../data/river_data.csv")
river_data_tisza_02_13 <- read.csv("../data/river_data_tisza_02_13.csv")
river_data_tisza_14_24 <- read.csv("../data/river_data_tisza_14_24.csv")
river_data_duna_02_13  <- read.csv("../data/river_data_duna_02_13.csv")
river_data_duna_14_24  <- read.csv("../data/river_data_duna_14_24.csv")
sum(is.na(river_data_duna_02_13))

# from the example
# Duna 2002–2013 maxima --------------------------------------------------------
# add simple time index

river_data_duna_02_13$t <- 1:nrow(river_data_duna_02_13)
while (!is.null(dev.list())) dev.off()

par(mfrow = c(2, 2))

# 1) time-index vs value: see extremes + dependence over time
plot(river_data_duna_02_13$t, river_data_duna_02_13$value,
     xlab = "Time index (days)",
     ylab = "Duna water level (cm)",
     cex = 1.25, cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)

# 2) year vs value: see long-run changes / trends
plot(river_data_duna_02_13$year, river_data_duna_02_13$value,
     xlab = "Year",
     ylab = "Duna water level (cm)",
     cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)

# 3) Q-Q plot of raw water levels: shows skewness / heavy tails
qqnorm(river_data_duna_02_13$value,
       ylab = "Normal Q-Q of water level",
       cex.lab = 1.25)
qqline(river_data_duna_02_13$value)

# for selecting threshold


nrow(river_data_duna_02_13[river_data_duna_02_13$value>650,])
# We note that only 26 values exceed 650cm and only 16 exceed 700 out of 4383 times.
# Therefore, in order to ensure we have enough data and to more easily interpret
# the mean residual life plot, we will restrict it to the range to 650.

dev.new(width = 8, height = 6) 
threshrange.plot(
  river_data_duna_02_13$value,
  r    = c(300, 650),  
  nint = 20
)

dev.new(width = 8, height = 6) 
mrlplot(river_data_duna_02_13$value,
        xlim = c(300, 650))        

#	In the threshrange.plot, both the reparameterized scale and the shape are fairly
# flat and stable between ~350 and ~430 cm.
#	Around 450 cm there is a clear kink: scale suddenly jumps up and shape drops
# to about −0.6 and keeps drifting, which is a sign the threshold is pushed too high.
# In the mrlplot, the mean excess curve is roughly linear from ~300–430 cm,
# then shows a noticeable change around 450 cm and more curvature afterwards.
# “pick the lowest threshold where parameters look stable and the MRL is roughly linear”
# so we choose 400 as upper boundary
# by the way it's almost the 95th percentile
quantile(river_data_duna_02_13$value, 0.95)

# dimension of the dataset:
nrow(river_data_duna_02_13)/12 #months
# 365.25

fitD <- fevd(value, river_data_duna_02_13, threshold = 400, type = "GP",
             time.units = "365.25/year")
dev.new(width = 8, height = 6) 
plot(fitD)
fitD

pextRemes(fitD, c(800, 900, 1000, 1100), lower.tail = FALSE)
?fevd

# give another t which counts how many days passed in that year (to control for seasonality)
river_data_duna_02_13$tobs <- ave(
  river_data_duna_02_13$year,
  river_data_duna_02_13$year,
  FUN = seq_along
)

fitD2 <- fevd(value, river_data_duna_02_13, threshold = 400, scale.fun =
                    ~ cos(2 * pi * tobs / 365.25) + sin(2 * pi * tobs / 365.25),
                  type = "GP", use.phi = TRUE)
dev.new(width = 8, height = 6) 
plot(fitD2)
lr.test(fitD, fitD2)
# p is 0.003516
# small p value -> the addition of the seasonality terms are statistically significant
dev.new(width = 8, height = 10) 
ci(fitD2, type = "parameter", which.par = 2, method = "proflik",
   xrange = c(-0.5, 0.01), verbose = TRUE)
dev.new(width = 8, height = 10) 
ci(fitD2, type = "parameter", which.par = 3, method = "proflik",
   xrange = c(-0.1, 0.2), verbose = TRUE)
ci(fitD2, type = "parameter")


v <- matrix(1, 730, 5)
v[, 2] <- cos(2 * pi * rep(1:365 / 365.25, 2))
v[, 3] <- sin(2 * pi * rep(1:365 / 365.25, 2))
v[, 5] <- 400
v <- make.qcov(fitD2, vals = v, nr = 730)
FCprobs <- pextRemes(fitD2, q = c(rep(600, 365), rep(700, 365)),
                        qcov = v, lower.tail = FALSE)
# FCprobs[1:365] = daily probabilities that the Duna water level exceeds 600 cm,
# for days 1–365 of the year, under the covariate pattern in v.
# FCprobs[366:730] = daily probabilities that the Duna water level exceeds 700 cm,
# for days 1–365 again.


# if i want the prob. there will be above 600 each day:
n <- nrow(river_data_duna_02_13)

# Build covariate matrix for each observation (1 row per day)
v <- matrix(1, n, 5)

# Seasonal terms based on tobs (day within year)
v[, 2] <- cos(2 * pi * river_data_duna_02_13$tobs / 365.25)
v[, 3] <- sin(2 * pi * river_data_duna_02_13$tobs / 365.25)

# Example: constant 5th covariate (whatever 400 means in your model)
v[, 5] <- 400

# Turn it into qcov in the format that matches fitD2
v_qcov <- make.qcov(fitD2, vals = v, nr = n)

# Probability of exceeding 600 cm on each actual day
probs_600 <- pextRemes(
  fitD2,
  q    = rep(600, n),   # same threshold for each day
  qcov = v_qcov,
  lower.tail = FALSE
)

# 5000-day return period in years
rp_5000_days <- 5000 / 365  # or 365.25

# Duna - Esztergom
RL_duna_02_13  <- return.level(fitD2,  return.period = rp_5000_days)
ci(RL_duna_02_13  <- fitD2,
   type = "return.level",
  return.period = rp_5000_days)

# same pattern for the other fits
