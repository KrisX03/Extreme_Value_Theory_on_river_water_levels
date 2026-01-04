# Computational Statistics
# Extreme Value Theory on river water levels

# Members of the team:
# Kristóf Andrási
# Gellért Banai
# Hunor Kuti
# Ákos Virág

# merging and cleaning

# ------------------------------------------------------------------------------
# importing neccesary libraries
# vector of required packages
pkgs <- c(
  "readxl",
  "dplyr",
  "purrr",
  "stringr",
  "tidyr",
  "lubridate"
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
# quick note: In theory, there should be no need to set the working directory manually, because when you open the file
# R sets it automatically. If this doesn’t work, use setwd to set the “codes” folder as the working directory.

# define rivers and years
rivers <- c("duna", "tisza")
years  <- 2002:2024

# defining all sheet names: "duna_2002", "duna_2003", ..., "tisza_2024"
sheets <- as.vector(outer(rivers, years, paste, sep = "_"))

# ?map_dfr
river_data <- map_dfr(sheets, function(sh) {
  read_excel(
    "../data/river_water_levels.xlsx",
    sheet     = sh,
    range     = "A5:M39", # importing only the important part of the sheet
    col_names = TRUE
  ) %>%
    slice(-c(1, 12, 23)) %>% # drop the empty rows
    # add river and year from sheet name
    mutate(
      river = sub("_.*", "", sh),
      year  = as.integer(str_extract(sh, "\\d{4}"))
    )
})

# now duna_2002–2024 and tisza_2002–2024 stacked together

# renaming the coloumns for data transition
colnames(river_data) <- c("day",rep(1:12), "river","year")

# transforming the data into long format for easier analysis
river_data <- river_data %>% 
  # turn month-columns 1:12 into rows
  pivot_longer(
    cols = `1`:`12`,          # the month columns
    names_to  = "month",
    values_to = "value"       # water level
  ) %>% 
  mutate(
    month = as.integer(month),
    # create ID like 2002_01_23
    date_id = sprintf("%d_%02d_%02d", year, month, day)
  ) %>% 
  select(date_id, river, year, month, day, value)


# check for NAs
sum(is.na(river_data)) # we have to filter for the months that are less than 31 days
river_data <- river_data %>% filter(!is.na(make_date(year, month, day)))

sum(is.na(river_data)) # 2 missing values -> dealing them at the end of the script
river_data[!complete.cases(river_data), ]


# reset the row numbering
# sort rows by date_id
river_data <- river_data[order(river_data$date_id), ]

# reset row numbering to 1:n
rownames(river_data) <- NULL


# for curiosity checking the leap year
nrow(river_data[river_data$year==2002 & river_data$river=="duna",])
nrow(river_data[river_data$year==2003 & river_data$river=="duna",])
nrow(river_data[river_data$year==2004 & river_data$river=="duna",]) # leap year
nrow(river_data[river_data$year==2005 & river_data$river=="duna",])

# now creating the necessary sub-dataframes
# Tisza 2002–2013
river_data_tisza_02_13 <- river_data[
  river_data$river == "tisza" &
    river_data$year >= 2002 &
    river_data$year <= 2013,
]

# Tisza 2014–2024
river_data_tisza_14_24 <- river_data[
  river_data$river == "tisza" &
    river_data$year >= 2014 &
    river_data$year <= 2024,
]

# Duna 2002–2013
river_data_duna_02_13 <- river_data[
  river_data$river == "duna" &
    river_data$year >= 2002 &
    river_data$year <= 2013,
]

# Duna 2014–2024
river_data_duna_14_24 <- river_data[
  river_data$river == "duna" &
    river_data$year >= 2014 &
    river_data$year <= 2024,
]


# dealing with NAs
na_idx <- which(is.na(river_data_duna_14_24$value))  # get the row indices where the 'value' column (water level) is NA
river_data_duna_14_24[na_idx, ]

# --- our initial idea was to replace the NAs with the mean of the 4 neighbouring days (i-4…i+4),
# code:
# for (i in na_idx) {
  # here we choose a small window around the missing value:
  # from i-4 to i+4, but cut off at the boundaries of the data frame
  #start_i <- max(1, i - 4)                                   # left boundary of the window (at least row 1)
  #end_i   <- min(nrow(river_data_duna_14_24), i + 4)         # right boundary of the window (at most last row)
  
  # take neighbour values in that window
  #neighbours <- river_data_duna_14_24$value[start_i:end_i]   # extract the 'value' entries in this local window
  
  # drop any NAs inside the window
  #neighbours <- neighbours[!is.na(neighbours)]               # keep only those neighbours that are not NA
  
  # replace NA with mean of neighbours
  # so each missing water level is replaced by the average of nearby days (up to 9 neighbours)
  #river_data_duna_14_24$value[i] <- mean(neighbours)
#}

# --- but this would not add any real information and could even distort the extremes
# --- since there were only 2 missing values out of 16802 observations, we decided to simply omit
# --- those rows instead of imputing them

river_data_duna_14_24 <- na.omit(river_data_duna_14_24)


# saving the dataframe
write.csv(river_data, file = "../data/river_data.csv", row.names = F)
write.csv(river_data_tisza_02_13, "../data/river_data_tisza_02_13.csv", row.names = F)
write.csv(river_data_tisza_14_24, "../data/river_data_tisza_14_24.csv", row.names = F)
write.csv(river_data_duna_02_13,  "../data/river_data_duna_02_13.csv",  row.names = F)
write.csv(river_data_duna_14_24,  "../data/river_data_duna_14_24.csv",  row.names = F)


