# RegioInvest Wind

Replication code for:

> **Spatially differentiated renewable subsidies: Good for ratepayers and bad for investors? Empirical evidence from Onshore Wind in Germany**  
> Yannik Pflugfelder and Christoph Weber

**RegioInvest Wind** is a MATLAB model of regional onshore-wind investment in Germany. It combines NUTS-3 geospatial data, wind yields, land availability, investment economics, turbine choice, and a nested-logit investment model. The replication code also evaluates counterfactual versions of the German reference-yield mechanism under capacity-based and energy-based policy targets.

The paper DOI will be added after publication.

---

## 1. Quick start

For a single model run, edit the configuration in:

```matlab
main_regio_invest.m
```

and run the script from MATLAB.

For the paper replication, run:

```matlab
run_paper_replication.m
```

If the model results already exist and only tables/figures need to be
regenerated, run:

```matlab
regenerate_paper_outputs.m
```


The paper-replication script sets the random-number seed to a fixed value, executes the required scenarios, saves the model results, and generates the model-based paper figures.

No machine-specific repository path needs to be entered manually.

---

## 2. Repository structure

The intended repository structure is:

```text
RegioInvestWind/
│
├── main_regio_invest.m
├── run_regio_invest_scenario.m
├── run_paper_replication.m
├── validate_paper_replication.m
├── README.md
├── DEPENDENCY_AUDIT.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── Model/
│   ├── calcRegioInvestWindOnshore.m
│   ├── estimateInvestmentDecision.m
│   ├── estimateLogitParameter.m
│   ├── initLogitParameter.m
│   ├── calcNetPresentValue.m
│   ├── transformNPV.m
│   ├── adaptPercentageWind.m
│   ├── energyYieldCorrectionYP.m
│   ├── getGeoShapeData.m
│   ├── getSpacePotential.m
│   ├── getReferenceYieldSchedule.m
│   │
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
│   └── ...
│
└── Results Paper/
    └── model results and generated paper figures
```

`Results Paper/` is created automatically if it does not yet exist.

MATLAB autosave files (`*.asv`), old dated plotting variants, and obsolete versions of model files are not part of the replication workflow.

---

## 3. Main entry points

### `main_regio_invest.m`

Convenience script for an individual model run. The main scenario settings are collected in a `config` struct and passed to `run_regio_invest_scenario.m`.

Example:

```matlab
config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'capacity';
config.RefYieldInInvestment = true;
config.runReferenceYieldSensitivity = false;
config.randomSeed = 1;
```

For an ordinary regional run, `runReferenceYieldSensitivity = false` evaluates only the current-policy case `s = 1`.

### `run_regio_invest_scenario.m`

Runs one reproducible model scenario. It:

- determines the repository path automatically;
- adds the repository to the MATLAB path;
- sets the random-number seed;
- initializes model parameters and geospatial inputs;
- calls `calcRegioInvestWindOnshore.m`;
- saves regional outputs and a compact scenario summary.

### `run_paper_replication.m`

Runs the model scenarios used for the paper and then calls the dedicated plotting scripts.

The expensive reference-yield sensitivity grid is evaluated only for the 2040 nested-logit policy experiments for which it is required.


At the end of the run, `validate_paper_replication.m` checks that all required
result files and figures are present, verifies the reference-yield grids and
the current-policy consistency of the two target regimes, and reports the
synthetic-net-benefit maxima and Table-2 fit statistic.

---

## 4. Expansion cases

The model supports four allocation logics through `expansionCase`:

| Case | Model logic | Role in the paper |
|---|---|---|
| 1 | Nested-logit / choice-based investment | Main model |
| 2 | Capacity proportional to existing capacity | Available in the model, not used in the final paper comparison |
| 3 | Capacity proportional to remaining available land | Figure 4 benchmark |
| 4 | NPV merit-order allocation using technically eligible land but no WindBG land-quota scaling and no reference-yield remuneration differentiation | Figure 4 benchmark |

The final paper comparison in Figure 4 labels these as:

- **Case A**: Case 1, nested-logit investment;
- **Case B**: Case 3, proportional to available land;
- **Case C**: Case 4, NPV merit order.

---

## 5. Reference-yield policy experiments

The reference-yield correction schedule is stored centrally in:

```matlab
Model/getReferenceYieldSchedule.m
```

Counterfactual correction factors are constructed as

```text
corr(s) = 1 + s * (corr(1) - 1)
```

where:

- `s = 0` gives uniform remuneration;
- `s = 1` reproduces the current reference-yield mechanism;
- `s > 1` strengthens spatial differentiation.

The full sensitivity analysis uses:

```matlab
scaleVec = 0:0.005:1.7;
```

Two policy-target logics are implemented:

- **capacity-based**: remuneration is calibrated to the installed-capacity target;
- **energy-based**: remuneration is calibrated to annual generation under the current-policy (`s = 1`) capacity equilibrium.

The energy target is therefore derived internally rather than hard-coded.

---

## 6. Probability-based policy evaluation

For the reference-yield counterfactuals, the model constructs a synthetic investment response from the estimated nested-logit probabilities.

At every counterfactual base-remuneration level, the synthetic fleet is mapped into its production-weighted effective remuneration. Probability-based producer surplus is calculated as the area under synthetic generation in this effective-remuneration space.

The associated probability-based synthetic net benefit combines:

- probability-based producer surplus;
- support payments; and
- avoided fuel and CO2 costs of displaced gas-fired generation.

These probability-based quantities are behaviorally based synthetic policy-evaluation measures. They are not exact money-metric compensating-variation measures.

The code additionally reports:

- an NPV-based producer-surplus benchmark; and
- a resource-cost-based welfare benchmark.

The formal construction is documented in the paper appendix.

---

## 7. Paper replication: scenarios and outputs

`run_paper_replication.m` currently executes the following model scenarios:

| Simulation | Expansion case | Policy target | Reference-yield sensitivity |
|---|---:|---|---|
| 2030 | 1 | Capacity | No, `s = 1` only |
| 2035 | 1 | Capacity | No, `s = 1` only |
| 2040 | 3 | Capacity | No |
| 2040 | 4 | Capacity | No |
| 2040 | 1 | Capacity | No reference-yield differentiation |
| 2040 | 1 | Capacity | Full `s` grid |
| 2040 | 1 | Energy | Full `s` grid |

The 2040 nested-logit capacity-sensitivity run also supplies the current-policy (`s = 1`) regional baseline used by several paper outputs.

### Mapping from code to paper outputs

| Paper output | Replication source |
|---|---|
| Table 2: nested-logit estimates | `estimateLogitParameter.m` + `write_table2_estimation_results.m` |
| Figure 2: spatial deployment 2024/2030/2035/2040 | `plot_spatial_deployment.m` |
| Figure 3: installed capacity by turbine type | `plot_turbine_types.m` and `plot_results.m` |
| Figure 4: allocation strategies by federal state | `plot_allocation_comparison.m` |
| Figure 5: reference-yield correction schedules | `plot_reference_yield_schedule.m` |
| Figure 6: regional effect of removing the reference-yield mechanism | `plot_ref_yield_diff.m` |
| Table 3: alternative reference-yield scaling | `estimateInvestmentDecision.m` |
| Figure 7: continuous reference-yield comparative statics | `plot_capacity_vs_energy_based_refYield.m` |
| Figure 9: probability-based vs. NPV-based producer surplus | `plot_probability_vs_npv_producer_surplus.m` |

Figure 1 is a conceptual model diagram rather than a direct simulation output.

Figure 8 is an illustrative GIS-input figure showing exclusion and suitable wind areas and is not generated by the core model-replication script.

---

## 8. Result files

Regional simulation results are saved in `Results Paper/` using descriptive filenames such as:

```text
results_nuts3_base2024_target2030_expCase1_capacity.mat
results_nuts3_base2024_target2035_expCase1_capacity.mat
results_nuts3_base2024_target2040_expCase1_capacity.mat
results_nuts3_base2024_target2040_expCase3_capacity.mat
results_nuts3_base2024_target2040_expCase4_capacity.mat
results_nuts3_base2024_target2040_expCase1_capacity_no_reference_yield.mat
```

Each regional result file contains:

- `selected_data`: NUTS-3 result table;
- `resultsWindOnshoreTurbinesUnstack`: installed capacity by turbine type and region;
- `scenarioSummary`: compact scenario-level statistics;
- `scenarioConfig`: configuration used for the run.

The two reference-yield sensitivity runs additionally produce:

```text
Results_scaling_capacity_based.mat
Results_scaling_energy_based.mat
```

These files contain the continuous counterfactual results used in Table 3 and Figure 7. The capacity-based file also contains the diagnostics used for Figure 9.

Generated figures are written to `Results Paper/` as vector PDFs.

---

## 9. Model output variables

Important regional output fields include:

- `capacityTotal`: installed capacity by NUTS-3 region;
- `capPerKm2`: installed capacity per square kilometre of total area;
- `capPerKm2avail`: installed capacity per available wind area;
- `capacityToAdd`: capacity expansion relative to the simulation-year starting fleet;
- `capacity_baseYear`: base-year capacity;
- `exhaustionProb`: land-use saturation indicator.

The reference-yield sensitivity output contains the variables required for the final paper, including:

- installed capacity;
- annual generation;
- base remuneration;
- production-weighted effective remuneration;
- annual support payments;
- avoided fuel and CO2 costs;
- probability-based producer surplus;
- probability-based synthetic net benefit;
- NPV-based producer surplus benchmark;
- resource-cost-based welfare benchmark;
- North, Neutral, and South capacity totals.

---

## 10. Input data requirements

Some source datasets cannot be redistributed in full because of licensing restrictions.

### Turbine data: Marktstammdatenregister (MaStR)

The model uses turbine data derived from the German Federal Network Agency's Marktstammdatenregister.

- Source: https://www.marktstammdatenregister.de
- Access date used for the replication data: 24 March 2025
- Legal basis: §111e EnWG
- License: Data licence Germany – attribution – Version 2.0 (DL-DE-BY-2.0)
- License text: https://www.govdata.de/dl-de/by-2-0

The authors filter the data for German onshore wind turbines, restrict commissioning dates to the relevant historical sample, and clean and aggregate the records for the model.

### ERA5 weather data

ERA5 data are not redistributed with the repository.

- Source: Copernicus Climate Data Store
- DOI: https://doi.org/10.24381/cds.adbb2d47
- Variables used: `u10`, `v10`, `u100`, `v100`
- Temporal resolution: hourly
- Spatial resolution: 0.25 degrees

Users must download the required ERA5 data directly from the Copernicus Climate Data Store:

https://cds.climate.copernicus.eu

### Power curves

The turbine-specific power-curve data are based on:

Pöstges & Weber (2023), *Energy Economics*  
https://doi.org/10.1016/j.eneco.2023.106534

### Geodata

The model uses NUTS-3 polygons and centroids from European Union GISCO (2021), together with the land-availability inputs required by the spatial model.

GISCO data are licensed under CC BY 4.0.

---

## 11. MATLAB requirements

The replication is written in MATLAB.

The code uses functionality including:

- `fminunc` for nested-logit estimation;
- shapefile and geospatial functions such as `shaperead`, `areaint`, and `referenceEllipsoid`;
- tables and grouped summaries;
- `polyshape` for generated maps;
- `exportgraphics` for vector figure output.

Accordingly, the workflow requires MATLAB together with the toolboxes providing the optimization and mapping functionality used by these routines.

---


### Parameter-estimation inference

Parameter estimates are obtained by maximum likelihood using `fminunc`.
For reproducible inference, standard errors are based on a central
finite-difference Hessian recomputed at the final optimum. The
parameter-specific perturbation is

```text
h_i = 3e-5 * max(1, abs(theta_i)).
```

The covariance matrix is the inverse of this positive-definite inference
Hessian. This avoids relying on the optimizer's internally generated
numerical Hessian, whose curvature estimates were found to be sensitive to
numerical settings for this model.

---

## 12. Reproducibility

The scenario runner uses:

```matlab
rng(1, 'twister');
```

by default, making random initialization of the parameter estimation reproducible.

A complete paper replication is computationally more demanding than a standard regional run because the capacity- and energy-based policy experiments evaluate the full reference-yield grid.

For ordinary model use, keep:

```matlab
config.runReferenceYieldSensitivity = false;
```

to evaluate only the current-policy equilibrium.



A replication dataset with processed inputs and model results is available through Zenodo:

https://doi.org/10.5281/zenodo.17062160

Users should consult the archive and the individual source-data licenses for redistribution restrictions.

---

## 13. Dependency audit

`DEPENDENCY_AUDIT.md` documents the model files retained in the final replication workflow and identifies legacy files that are not required.

In particular, the clean replication does **not** require:

- MATLAB autosave files (`*.asv`);
- old dated plotting scripts;
- old versions of `estimateInvestmentDecision.m`;
- `resultScriptShapefile.m`.

The latter was a legacy shapefile-export utility and is not needed to reproduce a final paper output.

---

## 14. License

The source code is licensed under the MIT License.

Input data remain subject to their original licenses, including:

- MaStR: Data licence Germany – attribution – Version 2.0;
- ERA5: Copernicus Climate Data Store terms;
- GISCO: CC BY 4.0.

Users are responsible for complying with the respective licenses for external datasets.

---

## 15. Citation

Please cite the paper when using the model. The final publication DOI will be added to `CITATION.cff` and this README after publication.

For the replication archive, see:

https://doi.org/10.5281/zenodo.17062160

---

## 16. Contact

**Yannik Pflugfelder**  
House of Energy, Climate and Finance  
University of Duisburg-Essen  
Germany

Email: yannik.pflugfelder@uni-due.de
