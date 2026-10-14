# Final repository structure

The clean replication repository should contain the files below.

```text
RegioInvestWind/
├── main_regio_invest.m
├── run_regio_invest_scenario.m
├── run_paper_replication.m
├── validate_paper_replication.m
├── README.md
├── DEPENDENCY_AUDIT.md
├── LICENSE
├── CITATION.cff
├── .gitignore
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
│   ├── write_table2_estimation_results.m
│   ├── plot_spatial_deployment.m
│   ├── plot_turbine_types.m
│   ├── plot_results.m
│   ├── plot_allocation_comparison.m
│   ├── plot_reference_yield_schedule.m
│   ├── plot_ref_yield_diff.m
│   ├── plot_capacity_vs_energy_based_refYield.m
│   └── plot_probability_vs_npv_producer_surplus.m
├── Input_Data/
│   └── [replication inputs]
└── Results Paper/
    └── [generated results]
```

Keep the existing `LICENSE`, `CITATION.cff`, `.gitignore`, and `Input_Data/`
from the repository unless they are intentionally updated separately.

Do not add `.asv` files, dated development copies, or legacy result scripts to
the replication release.
