%% REGENERATE PAPER OUTPUTS FROM EXISTING MODEL RESULTS
% Recreates tables/figures and runs validation without rerunning the model.
%
% Use this after all required MAT result files have already been generated.

clc;

repoPath = ...
    fileparts(mfilename('fullpath'));

addpath(genpath(repoPath));

run(fullfile(repoPath, 'Model', 'write_table2_estimation_results.m'));
run(fullfile(repoPath, 'Model', 'plot_spatial_deployment.m'));
run(fullfile(repoPath, 'Model', 'plot_turbine_types.m'));
run(fullfile(repoPath, 'Model', 'plot_allocation_comparison.m'));
run(fullfile(repoPath, 'Model', 'plot_reference_yield_schedule.m'));
run(fullfile(repoPath, 'Model', 'plot_ref_yield_diff.m'));
run(fullfile(repoPath, 'Model', 'plot_capacity_vs_energy_based_refYield.m'));
run(fullfile(repoPath, 'Model', 'plot_probability_vs_npv_producer_surplus.m'));
run(fullfile(repoPath, 'validate_paper_replication.m'));

fprintf('\nPaper outputs regenerated from existing results.\n');
