%% STANDARD-ERROR DIAGNOSTIC RUN
% Runs one 2040 Case-1 scenario with additional Hessian diagnostics.
%
% This does NOT evaluate the full reference-yield sensitivity grid.
% After the estimation, inspect the console table:
%   "STANDARD-ERROR DIAGNOSTICS"
%
% The diagnostic results are also stored in:
%   scenarioSummary.discreteChoiceEstimation.seDiagnostics

clc;
clear;

config = struct();

config.simYear = 2040;
config.expansionCase = 1;
config.iterationLogicInvestment = 'capacity';
config.RefYieldInInvestment = true;

config.runReferenceYieldSensitivity = false;
config.runStandardErrorDiagnostics = true;

config.plotResults = false;
config.saveRegionalResult = true;
config.randomSeed = 1;

run_regio_invest_scenario(config);
