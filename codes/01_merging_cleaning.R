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

# ------------------------------------------------------------------------------
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

sum(is.na(river_data)) # 2 missing values
river_data[!complete.cases(river_data), ]


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


# saving the dataframe
write.csv(river_data, file = "../data/river_data.csv", row.names = F)
write.csv(river_data_tisza_02_13, "../data/river_data_tisza_02_13.csv", row.names = F)
write.csv(river_data_tisza_14_24, "../data/river_data_tisza_14_24.csv", row.names = F)
write.csv(river_data_duna_02_13,  "../data/river_data_duna_02_13.csv",  row.names = F)
write.csv(river_data_duna_14_24,  "../data/river_data_duna_14_24.csv",  row.names = F)


