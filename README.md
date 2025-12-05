# Extreme Value Theory on River Water Levels  
## Analytical comparison of Hungarian rivers
###   Computational Statistics – Group Project

### Team members
- Kristóf Andrási  
- Gellért Banai  
- Hunor Kuti  
- Ákos Virág  

---
# egyelőre csak gpt-vel írattam egy ilyet de felülvizsgálásra vár, ha megvagyunk majd átírjuk
# nekem tetszik, de ez innetől kezdve itt is marad a historyban! egyébként szerintem jó, hogy a GPT képes összefoglalni formailag, amit csinálunk és GITHUB szalonképes formában legenerálni.
## Project overview

In this project we study **extreme water levels** of the two largest Hungarian rivers:

- **Duna (Danube) at Esztergom**
- **Tisza at Szeged**

for the period **2002–2024**.  
Using tools from **Extreme Value Theory (EVT)**, our aims are:

- to model **extreme maxima and minima** of daily water levels,
- to fit **Generalized Extreme Value (GEV)** distributions using the **Block Maxima** method,
- to fit **Generalized Pareto (GPD)** distributions using the **Peaks-Over-Threshold (POT)** method,
- to estimate how high a **dam** should be in both cities (for different time periods) so that overflow occurs only about **once every 5000 days**,
- to compare the statistical behaviour of extremes **before and after 2013** and discuss possible links to **climate change**.

We split the data into two subperiods:

- **2002–2013**
- **2014–2024**

and analyse **maxima and minima** separately for each river and period. The final goal is to provide (i) a short, policy-oriented summary for the Hungarian Home Office and (ii) a more technical statistical report.  [oai_citation:0‡Topic04 - Extreme Value Theory.pdf](sediment://file_00000000d998722fa976a337243b5bce)  

---

## Methods (short summary)


- **Block Maxima (BM)**  
  - Aggregate daily data into blocks (e.g. yearly or seasonal blocks).  
  - Extract block **maxima** (for floods) and **minima** (for low water).  
  - Fit a **GEV distribution** to the block extremes using the `extRemes` package.  

- **Peaks-Over-Threshold (POT)**  
  - Choose a high (or low) threshold for water levels.  
  - Model exceedances over the threshold using the **GPD**.  
  - Compare POT-based return levels with BM-based ones.

For both rivers and both periods, we:

- estimate GEV and GPD parameters for maxima and minima,
- compute relevant **return levels** (especially the level corresponding to a **5000-day return period**),
- compare fitted distributions across time (moments, quantiles, density shapes),
- assess whether changes are **statistically significant at the 5% level**.

Details of the modelling choices (block size, threshold selection, diagnostics) are documented in the report.

---

## Data

We use daily water level time series for:

- **Duna – Esztergom**, 2002–2024  
- **Tisza – Szeged**, 2002–2024  

The raw data contain daily water levels; after cleaning, we obtain:

- a tidy panel of daily observations for each river, year, month and day,
- derived datasets of block maxima/minima and threshold exceedances for EVT fitting.

> **Data source:** official Hungarian hydrological / water-level records (https://www.hydroinfo.hu/Html/archivum/archiv_tabla.html).

---

## Repository structure

Planned structure of the project repository:

```text
.
├── codes/
│   ├── 01_merging_cleaning.R      # import, cleaning and reshaping of raw data
│   └── 02_analysis.R              # EVT analysis: GEV/GPD fitting, plots, results
├── data/
│   ├── river_water_levels.xlsx    # raw downloaded river water levels
│   ├── river_data.csv             # full cleaned daily dataset (both rivers, all years)
│   ├── river_data_duna_02_13.csv  # Duna daily data, 2002–2013
│   ├── river_data_duna_14_24.csv  # Duna daily data, 2014–2024
│   ├── river_data_tisza_02_13.csv # Tisza daily data, 2002–2013
│   └── river_data_tisza_14_24.csv # Tisza daily data, 2014–2024
├── figures/
│   └── (plots from 02_analysis.R will be saved here)
├── documents/
│   └── (report, slides, notes, etc.)
└── README.md
