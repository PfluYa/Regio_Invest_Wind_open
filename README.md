# The future spatial distribution of onshore wind energy capacity based on a probabilistic investment calculus

**RegioInvest Wind** is a MATLAB-based model for simulating regional wind power deployment in Germany. The model projects future capacity expansions across NUTS-3 regions and turbine types under different allocation strategies and regulatory setting.
The model supports scenario-based energy system analysis and policy evaluation. It integrates geospatial constraints, investment economics, and historical turbine data. Key use cases include
* Calculate wind energy allocation under spatial constraints
* Comparing policy mechanisms (e.g. reference yield model, compensation design)
* Investigating regional differences in expansion patterns

Find methodological details in the following paper: DOI tba

### _main_regio_invest.m_: Main control script

## How to use
1. Clone this repo to your local machine.
2. Download and preprocess required input data (see below).
3. Adjust settings in main_regio_invest.m.
4. Run the script in MATLAB.
5. Visualize the results using the provided plotting scripts or use the code for extracting shapefiles for better visualizations via QGis etc.

## Input data requirements

Due to licensing restrictions, input files are not completely included.  
The following datasets are required:

* **Turbine Data (MaStR)**
  * This repository includes data derived from the Marktstammdatenregister of the German Federal Network Agency (Bundesnetzagentur).  
Source: https://www.marktstammdatenregister.de (accessed on 24.03.2025).  
© Bundesnetzagentur. The data are publicly available and may be reused in accordance with §111e EnWG.  
Modifications: [Filtered for Wind Turbines in Germany with commission date before 01.01.2025].

* **ERA5 Weather Data (not redistributed)**
  - Source: Copernicus Climate Data Store (CDS), DOI: https://doi.org/10.24381/cds.adbb2d47
  - Required variables: `u10`, `v10`, `u100`, `v100`
  - Format: NetCDF (`.nc`)
  - Temporal resolution: hourly  
    Spatial resolution: 0.25° grid
  - Download: ERA5 data can be obtained directly from the Copernicus Climate Data Store via the
    [download portal](https://cds.climate.copernicus.eu) as NetCDF files. Users should select the
    required variables (`u10`, `v10`, `u100`, `v100`), the temporal coverage (typically full years),
    and the spatial domain covering the study region.

* **Power Curve Data (included)**
  * CSV or MATLAB file with power curves per turbine type

* **Geodata (included)**
  * NUTS-3 region polygons and centroids © European Union, 2021 — GISCO, NUTS 2021 boundaries and points. Licensed under CC BY 4.0.  
  * Land-use exclusion zones

## Scenario configuration
The model supports four main expansion cases. These are controlled via _paraRegioInvest.expansionCase_:
* Case 1: Discrete choice model (nested logit) - standard case and new methodology
* Case 2: Allocation proportional to existing capacity
* Case 3: Allocation proportional to available land
* Case 4: Merit-order optimization based on NPV
Other key parameters
* _baseYear_: Year to use as reference for turbine data
* _simYear_: Target year for projected capacity
* _RefYieldInvestment_: Whether to apply EEG §36h compensation logic

## Methodological highlights
* Nested logit (Case 1): Estimates turbine choice probabilities based on NPV per turbine type and region
* Bias-corrected wind yield: ERA5 wind timeseries are corrected using a method proposed by Pflugfelder et al. (2025) to improve local accuracy. https://doi.org/10.1016/j.apenergy.2024.122890
* Spatial constraints: Expansion is bounded by remaining suitable land per region, considering already installed capacity.
* NPV Calculation: Accounts for CAPEX, OPEX, compensation, and optional correction factors.
To include negative NPVs in regression, the model apllies an _asinh_ transformation, see function _transformed_NPV_. This allows both positive and negative values to enter the model symetrically while avoiding the discontinuity of a log transformation.

## Outputs
The model generates a .mat file containing regional results, e.g.:
* _Capacity total_: Installed capacity in MW
* _capPerKm2_: Capacity density (total area)
* _capPerKm2avail_: Capacity density (available area)
* _newCapToAdd_: Net increase from base year
* _exhaustion_prob_: Utilization of available area

## Contact
Yannik Pflugfelder

University of Duisburg-Essen

House of Energy, Climate and Finance

yannik.pflugfelder@uni-due.de

Code is licensed under MIT. Data files are subject to their respective original licenses.
