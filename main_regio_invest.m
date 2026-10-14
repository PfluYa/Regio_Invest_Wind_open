%% REGIOINVEST WIND - MAIN SCRIPT
% Convenience entry point for a single reproducible model run.
%
% For the complete paper replication, use:
%   run_paper_replication.m
%
% No machine-specific repository path is required.

clc;
clear;

config = struct();

% Core scenario
config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'capacity';

% Current reference-yield mechanism
config.RefYieldInInvestment = true;

% Ordinary regional runs evaluate only s = 1.
% Set true only when reproducing Table 3 / Figures 7 and 9.
config.runReferenceYieldSensitivity = false;

% Optional generic result plots from plot_results.m
config.plotResults = false;

% Reproducible parameter initialization
config.randomSeed = 1;

[resultsWindOnshoreTurbinesUnstack, ...
 regioDataWind, ...
 resultFile] = ...
    run_regio_invest_scenario(config);
