%% RERUN FIGURE-4 CASE C AND REGENERATE PAPER OUTPUTS
% Recomputes only the corrected 2040 NPV-merit-order benchmark (Case C) and then
% regenerates the paper outputs from the existing result files.
%
% The expensive capacity- and energy-based sensitivity runs are NOT repeated.

clc;
clear;

repoPath = ...
    fileparts(mfilename('fullpath'));

addpath(genpath(repoPath));

config = struct();

config.randomSeed = 1;
config.baseYear = 2024;
config.simYear = 2040;
config.expansionCase = 4;
config.iterationLogicInvestment = 'capacity';
config.RefYieldInInvestment = true;
config.runReferenceYieldSensitivity = false;
config.plotResults = false;
config.saveRegionalResult = true;

run_regio_invest_scenario(config);

run(fullfile(repoPath, 'regenerate_paper_outputs.m'));
