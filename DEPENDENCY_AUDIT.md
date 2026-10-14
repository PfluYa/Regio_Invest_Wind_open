# Dependency audit for the open-source replication

This file documents which MATLAB source files are required by the final
paper-replication workflow.

## Entry points

- `main_regio_invest.m`: convenience script for one scenario.
- `run_regio_invest_scenario.m`: reproducible scenario runner.
- `run_paper_replication.m`: runs all model scenarios and paper figures.

## Core model files retained

- `Model/calcRegioInvestWindOnshore.m`
  - main data preparation and regional wind-investment model
  - contains local functions for turbine data, weather data, turbine
    production, reference-yield correction, cost reduction, and market
    preparation
- `Model/estimateInvestmentDecision.m`
  - nested-logit investment simulation
  - capacity- and energy-target calibration
  - reference-yield sensitivity analysis
  - probability-based producer surplus and synthetic net benefit
  - NPV-based and resource-cost benchmarks
- `Model/calcNetPresentValue.m`
  - NPV, discounted costs, and discounted revenues
- `Model/transformNPV.m`
  - asinh transformation used in estimation and simulation
- `Model/estimateLogitParameter.m`
  - nested-logit parameter estimation
- `Model/initLogitParameter.m`
  - starting values for parameter estimation
- `Model/adaptPercentageWind.m`
  - regional wind-land share adjustment
- `Model/energyYieldCorrectionYP.m`
  - wind-yield correction
- `Model/getGeoShapeData.m`
  - NUTS geodata and regional polygons
- `Model/getSpacePotential.m`
  - wind-space potential
- `Model/getReferenceYieldSchedule.m`
  - single source for the current reference-yield correction schedule

## Paper-output scripts retained

- `Model/write_table2_estimation_results.m` -> Table 2

- `Model/plot_spatial_deployment.m` -> Figure 2
- `Model/plot_turbine_types.m` / `Model/plot_results.m` -> Figure 3
- `Model/plot_allocation_comparison.m` -> Figure 4
- `Model/plot_reference_yield_schedule.m` -> Figure 5
- `Model/plot_ref_yield_diff.m` -> Figure 6
- `Model/plot_capacity_vs_energy_based_refYield.m` -> Figure 7
- `Model/plot_probability_vs_npv_producer_surplus.m` -> Figure 9

## Files not required in the final replication workflow

The following legacy files should not be part of the clean replication
workflow:

- MATLAB autosave files (`*.asv`)
- old dated plotting variants
- old versions of `estimateInvestmentDecision.m`
- old versions of `plot_capacity_vs_energy_based_refYield.m`
- `resultScriptShapefile.m` (legacy GIS export; not used for a paper result)
- dated intermediate result MAT files, except where intentionally supplied as
  precomputed replication outputs

## Notes

`run_regio_invest_scenario.m` sets a fixed random-number seed before model
estimation. The paper replication therefore does not rely on the user's
machine-specific MATLAB state.

The model still depends on the repository's `Input_Data` directory and the
MATLAB toolboxes/functions already required by the original code (including
mapping/geospatial functionality for shapefiles and maps).


## Inference convention

Standard errors reported by the final open-source replication are computed
from a central finite-difference Hessian evaluated at the final maximum-
likelihood optimum. The relative perturbation is `3e-5`, scaled by
`max(1, abs(theta_i))` for each parameter. The optimizer-provided Hessian is
retained in the model output for transparency but is not used for the final
Table 2 standard errors.
