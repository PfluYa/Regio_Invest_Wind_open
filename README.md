# RegioInvest Wind

MATLAB replication code for:

> **Spatially differentiated renewable subsidies: Good for ratepayers and bad for investors? Empirical evidence from Onshore Wind in Germany**  
> Yannik Pflugfelder and Christoph Weber

RegioInvest Wind is a spatial investment model for German onshore wind. It
combines NUTS-3 geospatial inputs, wind resources, land availability,
investment economics, turbine choice, and an estimated nested-logit model.
The replication evaluates future deployment and counterfactual changes in the
German reference-yield remuneration mechanism under capacity-based and
energy-based policy targets.

The replication archive is available at:

https://doi.org/10.5281/zenodo.17062160

---

## 1. Quick start

For the complete model-based paper replication, run from the repository root:

```matlab
run_paper_replication
```

The script:

1. runs all model scenarios required for the paper;
2. uses a fixed random-number seed for reproducibility;
3. estimates the nested-logit model;
4. evaluates the capacity- and energy-based reference-yield counterfactuals;
5. creates the model-based paper tables and figures; and
6. runs `validate_paper_replication.m` at the end.

For a single model scenario, edit and run:

```matlab
main_regio_invest
```

No machine-specific repository path needs to be entered manually.

---

## 2. Repository structure

```text
Regio_Invest_Wind_open/
│
├── main_regio_invest.m
├── run_regio_invest_scenario.m
├── run_paper_replication.m
├── validate_paper_replication.m
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── Model/
│   ├── adaptPercentageWind.m
│   ├── calcNetPresentValue.m
│   ├── calcRegioInvestWindOnshore.m
│   ├── energyYieldCorrectionYP.m
│   ├── estimateInvestmentDecision.m
│   ├── estimateLogitParameter.m
│   ├── getGeoShapeData.m
│   ├── getReferenceYieldSchedule.m
│   ├── getSpacePotential.m
│   ├── initLogitParameter.m
│   ├── transformNPV.m
│   │
│   ├── write_table2_estimation_results.m
│   ├── plot_spatial_deployment.m
│   ├── plot_turbine_types.m
│   ├── plot_results.m
│   ├── plot_allocation_comparison.m
│   ├── plot_reference_yield_schedule.m
│   ├── plot_ref_yield_diff.m
│   ├── plot_capacity_vs_energy_based_refYield.m
│   └── plot_probability_vs_npv_producer_surplus.m
│
├── Input_Data/
│   └── README.md
│
└── Results Paper/
    └── README.md
```

Temporary debugging scripts, MATLAB autosave files, old dated code copies, and
legacy shapefile-export scripts are intentionally excluded from the clean
replication repository.

---

## 3. Main entry points

### `main_regio_invest.m`

Convenience script for an individual scenario. Scenario settings are collected
in a `config` struct and passed to `run_regio_invest_scenario.m`.

Example:

```matlab
config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'capacity';
config.RefYieldInInvestment = true;
config.runReferenceYieldSensitivity = false;
config.randomSeed = 1;
```

For an ordinary regional model run,
`runReferenceYieldSensitivity = false` evaluates only the current-policy
equilibrium.

### `run_regio_invest_scenario.m`

Runs one reproducible scenario. It determines the repository path
automatically, sets the random-number seed, initializes the model inputs,
calls `calcRegioInvestWindOnshore.m`, and saves regional results plus a compact
scenario summary.

### `run_paper_replication.m`

Runs the full set of model scenarios required for the paper. The computationally
expensive full reference-yield sensitivity grid is evaluated only for the
2040 nested-logit counterfactuals for which it is required.

### `validate_paper_replication.m`

Runs automatically at the end of the paper replication. It checks that the
required output files are present, verifies the reference-yield grids and
current-policy consistency across target regimes, and reports the main
synthetic-net-benefit maxima and Table-2 fit statistic.

---

## 4. Expansion cases

The model supports four allocation logics through `expansionCase`.

| Case | Model logic | Role in the paper |
|---|---|---|
| 1 | Nested-logit / choice-based investment | Main model |
| 2 | Capacity proportional to existing capacity | Available in the model; not used in the final paper comparison |
| 3 | Capacity proportional to remaining available land | Figure 4 benchmark |
| 4 | NPV merit-order allocation | Figure 4 benchmark |

For Figure 4:

- **Case A** is Case 1, the estimated nested-logit model.
- **Case B** is Case 3, allocation proportional to available land.
- **Case C** is Case 4, NPV merit-order allocation subject to technically eligible land.

Case C abstracts from the reference-yield remuneration mechanism and WindBG
land-quota scaling while retaining underlying technical land eligibility.

---

## 5. Reference-yield policy experiments

The current reference-yield correction schedule is stored centrally in:

```matlab
Model/getReferenceYieldSchedule.m
```

Counterfactual correction factors are constructed as

```text
corr(s) = 1 + s * (corr(1) - 1)
```

where:

- `s = 0` gives uniform remuneration;
- `s = 1` reproduces the current mechanism; and
- `s > 1` strengthens spatial differentiation.


Two target logics are evaluated:

- **capacity-based:** remuneration is calibrated to the installed-capacity target;
- **energy-based:** remuneration is calibrated to annual generation under the
  current-policy (`s = 1`) capacity equilibrium.

The energy target is derived internally from the current-policy equilibrium
rather than being hard-coded.

---

## 6. Probability-based policy evaluation

For each reference-yield counterfactual, the model traces nested-logit
investment probabilities along a remuneration path from zero to the
scenario-specific equilibrium remuneration.

Alternative-specific probability responses are normalized relative to their
zero- and equilibrium-remuneration values and applied to the simulated
equilibrium fleet. This produces a synthetic generation response
`Q(rho,s)`.

At every counterfactual remuneration level, the synthetic fleet is mapped into
its own production-weighted effective remuneration. Probability-based producer
surplus is then calculated numerically as the area under synthetic generation
in this effective-remuneration space.

The associated probability-based synthetic net benefit combines:

- probability-based producer surplus;
- annual support payments; and
- avoided fuel and CO2 costs of displaced gas-fired generation.

These quantities are behaviorally based synthetic policy-evaluation measures,
not exact money-metric compensating-variation or welfare measures.

The model also reports:

- an NPV-based producer-surplus benchmark; and
- a resource-cost-based welfare benchmark.

The formal construction is documented in Appendix VIII.E of the paper.

---

## 7. Parameter estimation and standard errors

The nested-logit parameters are estimated by maximum likelihood using
`fminunc`.

---

## 8. Paper replication scenarios

`run_paper_replication.m` runs:

| Simulation | Expansion case | Policy target | Reference-yield sensitivity |
|---|---:|---|---|
| 2030 | 1 | Capacity | `s = 1` only |
| 2035 | 1 | Capacity | `s = 1` only |
| 2040 | 3 | Capacity | No |
| 2040 | 4 | Capacity | No |
| 2040 | 1 | Capacity | No reference-yield differentiation |
| 2040 | 1 | Capacity | Full `s` grid |
| 2040 | 1 | Energy | Full `s` grid |

The 2040 nested-logit capacity-sensitivity run also supplies the current-policy
regional baseline used by several paper outputs.

---

## 9. Mapping from code to paper outputs

| Paper output | Replication source |
|---|---|
| Table 2: nested-logit estimates | `estimateLogitParameter.m` + `write_table2_estimation_results.m` |
| Figure 2: spatial deployment 2024/2030/2035/2040 | `plot_spatial_deployment.m` |
| Figure 3: installed capacity by turbine type | `plot_turbine_types.m` + `plot_results.m` |
| Figure 4: allocation strategies by federal state | `plot_allocation_comparison.m` |
| Figure 5: reference-yield correction schedules | `plot_reference_yield_schedule.m` |
| Figure 6: regional effect of removing the reference-yield mechanism | `plot_ref_yield_diff.m` |
| Table 3: selected reference-yield counterfactuals | `estimateInvestmentDecision.m` |
| Figure 7: continuous reference-yield comparative statics | `plot_capacity_vs_energy_based_refYield.m` |
| Figure 9: probability-based vs. NPV-based producer surplus | `plot_probability_vs_npv_producer_surplus.m` |

Figure 1 is a conceptual model diagram rather than a direct simulation output.

Figure 8 is an illustrative GIS-input figure and is not generated by the core
paper-replication script.

For Figure 6, relative capacity changes are clipped at +/-50 percent. NUTS-3
regions with less than 5 MW of installed capacity in both the with- and
without-reference-yield scenarios are displayed in grey because a percentage
change is economically uninformative for such small absolute quantities.
Undefined relative differences are also shown in grey.

---

## 10. Result files

Regional simulation outputs are saved in `Results Paper/` with descriptive
filenames such as:

```text
results_nuts3_base2024_target2030_expCase1_capacity.mat
results_nuts3_base2024_target2035_expCase1_capacity.mat
results_nuts3_base2024_target2040_expCase1_capacity.mat
results_nuts3_base2024_target2040_expCase3_capacity.mat
results_nuts3_base2024_target2040_expCase4_capacity.mat
results_nuts3_base2024_target2040_expCase1_capacity_no_reference_yield.mat
```

Each regional result contains:

- `selected_data`: NUTS-3 result table;
- `resultsWindOnshoreTurbinesUnstack`: installed capacity by turbine type and region;
- `scenarioSummary`: scenario-level summary statistics; and
- `scenarioConfig`: configuration used for the run.

The full policy-sensitivity runs additionally create:

```text
Results_scaling_capacity_based.mat
Results_scaling_energy_based.mat
```

These contain the continuous counterfactual results used in Table 3 and
Figure 7. The capacity-based result also contains diagnostics used for Figure 9.

Generated paper figures are exported as vector PDFs. Table 2 is additionally
written as CSV and MAT output.

---

## 11. Important model outputs

The reference-yield sensitivity output contains the variables used in the
paper, including:

- installed capacity;
- annual generation;
- base remuneration;
- production-weighted effective remuneration;
- annual support payments;
- avoided fuel and CO2 costs;
- probability-based producer surplus;
- probability-based synthetic net benefit;
- NPV-based producer surplus benchmark;
- resource-cost-based welfare benchmark; and
- North, Neutral, and South capacity totals.

Important NUTS-3 regional result fields include:

- `capacityTotal`;
- `capPerKm2`;
- `capPerKm2avail`;
- `capacityToAdd`;
- `capacity_baseYear`; and
- `exhaustionProb`.

---

## 12. Input data

Some source datasets cannot be redistributed in full through the public
repository because of licensing restrictions or file size. See
`Input_Data/README.md` and the Zenodo replication archive.

### Marktstammdatenregister

The model uses German onshore-wind installation data derived from the Federal
Network Agency's Marktstammdatenregister (MaStR).

- Source: https://www.marktstammdatenregister.de
- License: Data licence Germany – attribution – Version 2.0 (DL-DE-BY-2.0)
- License text: https://www.govdata.de/dl-de/by-2-0

### ERA5 weather data

ERA5 weather files are not redistributed with this replication package due
to their size. They can be downloaded directly from the Copernicus Climate
Data Store (ERA5 hourly data on single levels).

Source DOI:
https://doi.org/10.24381/cds.adbb2d47

The default paper replication uses the years 2015–2024.

Required variables:
- 10 m u-component of wind (`u10`)
- 10 m v-component of wind (`v10`)
- 100 m u-component of wind (`u100`)
- 100 m v-component of wind (`v100`)

Temporal resolution:
- hourly

Spatial resolution:
- 0.25 degrees
- the downloaded domain must cover Germany and all NUTS-3 centroids used by
  the model

Files must be stored in `Input_Data/` using the following naming convention:

ERA5_data_windspeed_2015.nc
ERA5_data_windspeed_2016.nc
...
ERA5_data_windspeed_2024.nc

The MATLAB code reads the NetCDF variables `u10`, `v10`, `u100`, `v100`,
`latitude`, and `longitude` directly from these files.

### Power curves

Turbine-specific power-curve inputs are based on Pöstges and Weber (2023),
*Energy Economics*:

https://doi.org/10.1016/j.eneco.2023.106534

### Geodata

The model uses NUTS-3 geodata and the land-availability inputs required by the
spatial model. GISCO material is subject to its original licensing terms.

---

## 13. MATLAB requirements

The replication is written in MATLAB and uses functionality including:

- `fminunc` for nested-logit estimation;
- geospatial/shapefile functions such as `shaperead`, `areaint`, and
  `referenceEllipsoid`;
- tables and grouped summaries;
- `polyshape`; and
- `exportgraphics`.

The corresponding MATLAB toolboxes must be installed.

---

## 14. Reproducibility

The scenario runner uses

```matlab
rng(1, 'twister');
```

by default.

A complete paper replication is computationally more demanding than a single
scenario because the capacity- and energy-based policy experiments evaluate
the full reference-yield grid.

For ordinary model use, keep:

```matlab
config.runReferenceYieldSensitivity = false;
```

The full replication run concludes with `validate_paper_replication.m`. A
successful run reports completion of all required output checks and the
current-policy consistency of the capacity- and energy-target simulations.

---

## 15. License

The source code is licensed under the MIT License. External input data remain
subject to their original licenses and terms of use.

See `LICENSE`.

---

## 16. Citation

The replication archive is available at:

https://doi.org/10.5281/zenodo.17062160

See `CITATION.cff` for software citation metadata. Please also cite the
associated paper when using the model. The publication DOI can be added after
publication.

---

## 17. Contact

**Yannik Pflugfelder**  
House of Energy, Climate and Finance  
University of Duisburg-Essen  
Germany

Email: yannik.pflugfelder@uni-due.de
