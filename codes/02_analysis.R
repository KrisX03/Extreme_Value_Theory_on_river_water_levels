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
  "dplyr"
)

# Install any that are missing
to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_install) > 0) {
  install.packages(to_install)
}

# loading them
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# loading the dataframes
river_data <- read.csv("../data/river_data.csv")
river_data_tisza_02_13 <- read.csv("../data/river_data_tisza_02_13.csv")
river_data_tisza_14_24 <- read.csv("../data/river_data_tisza_14_24.csv")
river_data_duna_02_13  <- read.csv("../data/river_data_duna_02_13.csv")
river_data_duna_14_24  <- read.csv("../data/river_data_duna_14_24.csv")







