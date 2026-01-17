# Extreme Value Theory on River Water Levels  
## Analytical comparison of Hungarian rivers
###   Computational Statistics – Group Project



### Team members
- Kristóf Andrási  
- Gellért Banai  
- Hunor Kuti  
- Ákos Virág  

---

## Project overview

This repository contains the code, data and figures for our course project on **flood and drought risk** on Hungary’s two largest rivers:

- **Danube (Duna) at Esztergom**
- **Tisza at Szeged**

We use **Extreme Value Theory (EVT)** on **daily water-level data from 2002–2024** to answer two main questions:

1. Have extreme high and low water levels changed between **2002–2013** and **2014–2024**?
2. What dam height would be needed if we want overflow only about **once every 5000 days** (~13.7 years)?

The results are reported in two formats (see `report.pdf`):

- a **1-page policy brief** for the Hungarian Home Office (plain language, focused on resource allocation), and  
- a **technical 3-page summary** for statisticians (methods, assumptions, tests, diagnostics).

---

## Data

Daily water levels for:

- **Duna – Esztergom, 2002–2024**  
- **Tisza – Szeged, 2002–2024**

Source: *Országos Vízjelző Szolgálat / Hydroinfo* online archive.

Main preprocessing steps:

- convert raw tables into a **tidy daily panel** (river, date, year, value),
- split the series into two periods: **2002–2013** (early) and **2014–2024** (late),
- handle a very small number of missing daily values and remove obvious artefacts,
- create derived data sets of **annual maxima/minima** and **threshold exceedances**.

---

## Methods

We combine the two standard EVT approaches used in the course.

### 1. Block Maxima (GEV)

- For each river and period we compute **annual maxima** (floods) and, after sign-flipping, **annual minima** (droughts).  
- We fit **Generalized Extreme Value (GEV)** models using maximum likelihood (`extRemes::fevd`).  
- To test for changes between 2002–2013 and 2014–2024, we fit  
  - a **stationary GEV** (no time effect), and  
  - a **non-stationary GEV** where **location and/or scale** depend on the period.  
- We use **Likelihood Ratio tests** to check whether adding period-specific parameters improves fit at the 5% level and compute **5000-day return levels** (design dam heights) with confidence intervals.

### 2. Peaks-Over-Threshold (POT, GPD)

- To better use the short time series, we also model **all exceedances over a high threshold**:  
  - high flows → raw levels above a high quantile,  
  - low flows → exceedances of **–water-level** above a high threshold.  
- Thresholds are guided by **Mean Residual Life plots** and **parameter stability plots**; for comparability we mainly use the **95th percentile** in each series.  
- Exceedances are modelled with a **Generalized Pareto Distribution (GPD)** (`type = "GP"` in `fevd`).  
- We compare **scale and shape parameters** between periods using **Wald-type tests** based on the parameter covariance matrices and compute **5000-day return levels** from the GPD fits.

---

## Main findings (very short)

(Full details and plots are in `report.pdf`.)

- **Tisza (Szeged)**  
  - Flood extremes are **clearly lower** in 2014–2024 than in 2002–2013; both GEV and GPD show a marked drop in the 5000-day return level and parameter-change tests strongly reject stationarity.  
  - Low-water extremes show **no strong worsening**; changes in drought risk are small.

- **Danube (Esztergom)**  
  - Flood and low-water extremes are **slightly lower** in the later period, but evidence is weaker: GEV tests do **not** reject stationarity at 5%, and confidence intervals for return levels overlap.

Because each period covers only ~11 years, the **uncertainty is substantial**, especially for 5000-day events. The analysis is therefore best seen as **characterising recent decades**, not as a long-term climate projection.

---

## Repository structure

```text
.
├── codes/
│   ├── 01_merging_cleaning.R      # import, cleaning and reshaping of raw data
│   ├── 02_analysis_GEV.R          # EVT block maxima analysis: GEV fitting + diagnostics plots
│   └── 03_analysis_GPD.R          # EVT POT analysis: GPD fitting + diagnostics/threshold plots
├── data/
│   ├── river_water_levels.xlsx    # raw downloaded river water levels
│   ├── river_data.csv             # full cleaned daily dataset (both rivers, all years)
│   ├── river_data_duna_02_13.csv  # Duna daily data, 2002–2013
│   ├── river_data_duna_14_24.csv  # Duna daily data, 2014–2024
│   ├── river_data_tisza_02_13.csv # Tisza daily data, 2002–2013
│   └── river_data_tisza_14_24.csv # Tisza daily data, 2014–2024
├── figures/
│   ├── danube_gev_max_2002_2013_diagnostics.png
│   ├── danube_gev_max_2014_2024_diagnostics.png
│   ├── danube_gev_min_2002_2013_diagnostics.png
│   ├── danube_gev_min_2014_2024_diagnostics.png
│   ├── mrl_duna_2002_2013.png
│   ├── threshrange_duna_2002_2013.png
│   ├── tisza_gev_max_2002_2013_diagnostics.png
│   ├── tisza_gev_max_2014_2024_diagnostics.png
│   ├── tisza_gev_min_2014_2024_diagnostics.png
│   └── tisza_gumbel_min_2002_2013_diagnostics.png
│
├── report.pdf                     # final written report (main output)
└── README.md
```

### Licence
MIT License (MIT): see the [License File](https://github.com/sensiolabs/GotenbergBundle/blob/1.x/LICENSE) for more details.

