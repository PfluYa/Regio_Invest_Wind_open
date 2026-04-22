# The impact of investor choices on the future spatial development of onshore wind energy

## Overview

**RegioInvest Wind** is a MATLAB-based model for simulating regional onshore wind power deployment in Germany.  
The model projects future capacity expansions across NUTS-3 regions and turbine types under different allocation strategies and regulatory settings.

It supports scenario-based energy system analysis and policy evaluation by integrating:
- geospatial constraints,
- investment economics, and
- historical turbine data.

Key use cases include:
- Calculation of wind energy allocation under spatial constraints  
- Comparison of policy mechanisms (e.g. reference yield model, compensation design)  
- Analysis of regional differences in expansion patterns  

Find methodological details in the following paper: **DOI tba**

---

## Main entry point

### `main_regio_invest.m`

This script is the central control file of the model.  
All simulations are configured and executed from here.

---

## How to use

1. Clone this repository to your local machine  
2. Download and prepare the required input data (see below)  
3. Adjust settings in `main_regio_invest.m`  
4. Run the script in MATLAB  
5. Visualize results using the provided plotting scripts or export shapefiles for external tools (e.g. QGIS)

---

## Input data requirements

Due to licensing restrictions, input files are not completely included.  
The following datasets are required:

### Turbine Data (MaStR)

This repository includes data derived from the **Marktstammdatenregister (MaStR)** of the German Federal Network Agency (Bundesnetzagentur).

- Source: https://www.marktstammdatenregister.de (accessed on 24.03.2025)  
- Legal basis: §111e EnWG  
- License: Data licence Germany – attribution – Version 2.0 (DL-DE-BY-2.0)  
- License text: https://www.govdata.de/dl-de/by-2-0  

Modifications:
- Filtered for onshore wind turbines in Germany  
- Commissioning date before 01.01.2025  
- Cleaned, filtered, and aggregated by the authors  

A reduced example dataset is included for testing.

---

### ERA5 Weather Data (not redistributed)

- Source: Copernicus Climate Data Store (CDS)  
- DOI: https://doi.org/10.24381/cds.adbb2d47  

Required variables:
- `u10`, `v10`, `u100`, `v100`

Specifications:
- Format: NetCDF (`.nc`)  
- Temporal resolution: hourly  
- Spatial resolution: 0.25° grid  

Download:
ERA5 data must be obtained directly from the CDS portal:  
https://cds.climate.copernicus.eu  

Users should select:
- variables: `u10`, `v10`, `u100`, `v100`  
- full-year hourly data  
- spatial domain covering Germany  

---

### Power Curve Data (included)

- CSV or MATLAB file with turbine-specific power curves  
- Based on:  
  Pöstges & Weber (2023), *Energy Economics*  
  https://doi.org/10.1016/j.eneco.2023.106534  

---

### Geodata (included)

- NUTS-3 region polygons and centroids  
- Source: European Union, GISCO (2021)  
- License: CC BY 4.0  

Includes:
- region geometries  
- centroids  
- land-use exclusion zones  

---

## Scenario configuration

The model supports four expansion cases, controlled via  
`paraRegioInvest.expansionCase`:

- **Case 1**: Discrete choice model (nested logit) – baseline methodology  
- **Case 2**: Allocation proportional to existing capacity  
- **Case 3**: Allocation proportional to available land  
- **Case 4**: Merit-order optimization based on NPV  

Key parameters:

- `baseYear` – reference year for turbine data  
- `simYear` – target year for projected capacity  
- `RefYieldInInvestment` – application of EEG §36h compensation logic  

---

## Methodological highlights

- **Nested logit model (Case 1)**  
  Estimates turbine investment choices based on relative profitability  

- **Wind yield estimation**  
  Based on ERA5 reanalysis data  

- **Bias correction of wind yields**  
  Following Pflugfelder et al. (2025):  
  https://doi.org/10.1016/j.apenergy.2024.122890  

- **Spatial constraints**  
  Expansion limited by available land and existing installations  

- **Net Present Value (NPV)**  
  Includes CAPEX, OPEX, compensation, and discounting  

- **NPV transformation**  
  Uses an `asinh` transformation (see `transformNPV.m`)  
  to include both positive and negative values symmetrically  

- **Reference yield adjustment**  
  Enables analysis of spatial allocation and economic efficiency  
  under alternative policy designs  

---

## Outputs

The model generates `.mat` files in the `Results/` folder.

Key output variables:

- `capacityTotal` – Installed capacity per region (MW)  
- `capPerKm2` – Capacity density (total area)  
- `capPerKm2avail` – Capacity density (available area)  
- `capacityToAdd` – Net capacity expansion relative to base year  
- `exhaustionProb` – Degree of land-use saturation  

---

## Reproducibility

A full replication dataset (including processed input data and model results) is available via Zenodo:

**DOI: [Zenodo](https://doi.org/10.5281/zenodo.17062160)**

This dataset includes:
- processed MaStR data  
- model outputs (`.mat` and `.xlsx`)  
- additional input data where redistribution is permitted  

---

## License

- Code: MIT License  
- Data: subject to the original data providers and their respective licenses  

---

## Contact

**Yannik Pflugfelder**  
University of Duisburg-Essen  
House of Energy, Climate and Finance  

📧 yannik.pflugfelder@uni-due.de
