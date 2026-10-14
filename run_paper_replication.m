%% REGIOINVEST WIND - PAPER REPLICATION
% Runs the model scenarios required for the paper.
%
% The expensive reference-yield sensitivity grid is evaluated only for the
% two 2040 Case-1 policy-target scenarios used in Table 3 and Figure 7.
%
% Dedicated figure scripts are run after the corresponding result files
% have been created.

clc;
clear;

repoPath = ...
    fileparts(mfilename('fullpath'));

addpath(genpath(repoPath));

%% ------------------------------------------------------------------------
% 1. COMMON SETTINGS
% -------------------------------------------------------------------------

baseConfig = struct();

baseConfig.randomSeed = 1;
baseConfig.baseYear = 2024;
baseConfig.iterationLogicInvestment = 'capacity';
baseConfig.RefYieldInInvestment = true;
baseConfig.runReferenceYieldSensitivity = false;
baseConfig.plotResults = false;
baseConfig.saveRegionalResult = true;

%% ------------------------------------------------------------------------
% 2. REGIONAL BASELINE RUNS USED IN THE PAPER
% -------------------------------------------------------------------------

% 2030: nested-logit baseline
config = baseConfig;
config.simYear = 2030;
config.expansionCase = 1;
run_regio_invest_scenario(config);

% 2035: nested-logit baseline
config = baseConfig;
config.simYear = 2035;
config.expansionCase = 1;
run_regio_invest_scenario(config);

% 2040: land-proportional benchmark
config = baseConfig;
config.simYear = 2040;
config.expansionCase = 3;
run_regio_invest_scenario(config);

% 2040: NPV merit-order benchmark
config = baseConfig;
config.simYear = 2040;
config.expansionCase = 4;
run_regio_invest_scenario(config);

% 2040: counterfactual without reference-yield differentiation
config = baseConfig;
config.simYear = 2040;
config.expansionCase = 1;
config.RefYieldInInvestment = false;
run_regio_invest_scenario(config);

%% ------------------------------------------------------------------------
% 3. CAPACITY-BASED REFERENCE-YIELD SENSITIVITY
% -------------------------------------------------------------------------

config = baseConfig;
config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'capacity';
config.RefYieldInInvestment = true;
config.runReferenceYieldSensitivity = true;

run_regio_invest_scenario(config);

%% ------------------------------------------------------------------------
% 4. ENERGY-BASED REFERENCE-YIELD SENSITIVITY
% -------------------------------------------------------------------------

config = baseConfig;
config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'energy';
config.RefYieldInInvestment = true;
config.runReferenceYieldSensitivity = true;

run_regio_invest_scenario(config);

%% ------------------------------------------------------------------------
% 5. PAPER TABLE OUTPUT
% -------------------------------------------------------------------------

run(fullfile(repoPath, 'Model', 'write_table2_estimation_results.m'));


%% ------------------------------------------------------------------------
% 6. PAPER FIGURES GENERATED FROM THE MODEL OUTPUTS
% -------------------------------------------------------------------------

run(fullfile(repoPath, 'Model', 'plot_spatial_deployment.m'));
run(fullfile(repoPath, 'Model', 'plot_turbine_types.m'));
run(fullfile(repoPath, 'Model', 'plot_allocation_comparison.m'));
run(fullfile(repoPath, 'Model', 'plot_reference_yield_schedule.m'));
run(fullfile(repoPath, 'Model', 'plot_ref_yield_diff.m'));
run(fullfile(repoPath, 'Model', 'plot_capacity_vs_energy_based_refYield.m'));
run(fullfile(repoPath, 'Model', 'plot_probability_vs_npv_producer_surplus.m'));


%% ------------------------------------------------------------------------
% 7. VALIDATE REPLICATION OUTPUTS
% -------------------------------------------------------------------------

run(fullfile(repoPath, 'validate_paper_replication.m'));

fprintf('\n');
fprintf('============================================================\n');
fprintf('PAPER REPLICATION MODEL RUNS COMPLETED\n');
fprintf('Figures 2, 3, 4, 5, 6, 7, and 9 have been generated.\n');
fprintf('Replication validation completed.\n');
fprintf('============================================================\n');
