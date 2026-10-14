function [resultsWindOnshoreTurbinesUnstack, regioDataWind, resultFile] = ...
    run_regio_invest_scenario(config)
% RUN_REGIO_INVEST_SCENARIO
%
% Runs one RegioInvest Wind scenario with reproducible settings.
%
% The function is used both by main_regio_invest.m for an individual model
% run and by run_paper_replication.m for the paper replication.
%
% INPUT
%   config : optional struct. Supported fields are documented below.
%
% OUTPUTS
%   resultsWindOnshoreTurbinesUnstack : turbine-level simulation results
%   regioDataWind                     : regional NUTS-3 result table
%   resultFile                        : saved regional result file
%
% The repository root is determined automatically from this file. No
% machine-specific path needs to be edited.

if nargin < 1
    config = struct();
end

%% ------------------------------------------------------------------------
% 0. REPOSITORY AND REPRODUCIBILITY
% -------------------------------------------------------------------------

repoPath = ...
    fileparts(mfilename('fullpath'));

originalFolder = ...
    pwd;

cleanupFolder = ...
    onCleanup(@() cd(originalFolder)); %#ok<NASGU>

cd(repoPath);
addpath(genpath(repoPath));

randomSeed = ...
    getConfig(config, 'randomSeed', 1);

rng(randomSeed, 'twister');

clear global paraRegioInvest optsRegioInvest baseTableRegioInvest
global paraRegioInvest optsRegioInvest baseTableRegioInvest

%% ------------------------------------------------------------------------
% 1. MODEL YEAR AND SCENARIO CONFIGURATION
% -------------------------------------------------------------------------

paraRegioInvest.weatherYear = ...
    getConfig(config, 'weatherYear', 2015:2024);

paraRegioInvest.baseYear = ...
    getConfig(config, 'baseYear', 2024);

paraRegioInvest.simYear = ...
    getConfig(config, 'simYear', 2040);

paraRegioInvest.expansionCase = ...
    getConfig(config, 'expansionCase', 1);

paraRegioInvest.geoScope = ...
    getConfig(config, 'geoScope', {'DE'});

optsRegioInvest.geoScope = ...
    table( ...
        "geoScope.xlsx", ...
        "DE", ...
        'VariableNames', {'filename', 'tablename'});

%% ------------------------------------------------------------------------
% 2. PREPARE BASE GEOGRAPHIC TABLE
% -------------------------------------------------------------------------

baseTableRegioInvest = ...
    getGeoShapeData();

baseTableRegioInvest = ...
    getSpacePotential(baseTableRegioInvest);

%% ------------------------------------------------------------------------
% 3. WIND ONSHORE MODEL OPTIONS
% -------------------------------------------------------------------------

optsRegioInvest.windOnshoreReferenceYieldModel = ...
    getConfig(config, 'windOnshoreReferenceYieldModel', true);

optsRegioInvest.uniqueRefYieldperRegion = ...
    getConfig(config, 'uniqueRefYieldperRegion', false);

optsRegioInvest.RefYieldInInvestment = ...
    getConfig(config, 'RefYieldInInvestment', true);

optsRegioInvest.RefYieldLinearized = ...
    getConfig(config, 'RefYieldLinearized', false);

optsRegioInvest.windOnshoreYieldCorrectionYP = ...
    getConfig(config, 'windOnshoreYieldCorrectionYP', true);

optsRegioInvest.windOnshoreMethodology = ...
    getConfig(config, 'windOnshoreMethodology', 'binomialNestedLogit');

optsRegioInvest.iterationLogicInvestment = ...
    getConfig(config, 'iterationLogicInvestment', 'capacity');

optsRegioInvest.windOnshoreExplanatoryVariable = ...
    getConfig(config, 'windOnshoreExplanatoryVariable', {'netPresentValue'});

optsRegioInvest.windOnshoreNpvAdjustment = ...
    getConfig(config, 'windOnshoreNpvAdjustment', 'NPVperArea');

optsRegioInvest.windOnshoreExplainedVariable = ...
    getConfig(config, 'windOnshoreExplainedVariable', {'relativeUsedWindSpaceDcm'});

% Full reference-yield sensitivity is expensive and is only required for
% Table 3 and Figures 7 and 9. Ordinary regional runs use s = 1 only.
optsRegioInvest.runReferenceYieldSensitivity = ...
    getConfig(config, 'runReferenceYieldSensitivity', false);

%% ------------------------------------------------------------------------
% 4. LAND AND POWER-POTENTIAL ASSUMPTIONS
% -------------------------------------------------------------------------

optsRegioInvest.scalingAreaParameter = ...
    getConfig(config, 'scalingAreaParameter', 'UB21');

optsRegioInvest.scalingAreaInvest = ...
    getConfig(config, 'scalingAreaInvest', 'WindBG');

optsRegioInvest.powerPotential_increase = ...
    getConfig(config, 'powerPotential_increase', true);

optsRegioInvest.powerPotential_individual = ...
    getConfig(config, 'powerPotential_individual', false);

%% ------------------------------------------------------------------------
% 5. ECONOMIC PARAMETERS
% -------------------------------------------------------------------------

paraRegioInvest.windOnshoreTimeHorizonEstimation = ...
    getConfig(config, 'windOnshoreTimeHorizonEstimation', 10);

paraRegioInvest.windOnshoreCompensation = ...
    getConfig(config, 'windOnshoreCompensation', 0.080);

paraRegioInvest.windOnshoreInterestRate = ...
    getConfig(config, 'windOnshoreInterestRate', 0.05);

paraRegioInvest.windOnshorePowerPotential = ...
    getConfig(config, 'windOnshorePowerPotential', 22500);

paraRegioInvest.windOnshorePowerPotential_invest = ...
    getConfig(config, 'windOnshorePowerPotential_invest', 25000);

paraRegioInvest.windOnshoreTurbineLifetime = ...
    getConfig(config, 'windOnshoreTurbineLifetime', 22);

paraRegioInvest.opexPerkW_fix = ...
    getConfig(config, 'opexPerkW_fix', 15);

paraRegioInvest.opexPerkWh_var = ...
    getConfig(config, 'opexPerkWh_var', 0.0075);

paraRegioInvest.costReductionRate = ...
    getConfig(config, 'costReductionRate', 0.0238);

%% ------------------------------------------------------------------------
% 6. RUN MODEL
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf('REGIOINVEST WIND SCENARIO\n');
fprintf('Simulation year: %d\n', paraRegioInvest.simYear);
fprintf('Expansion case:  %d\n', paraRegioInvest.expansionCase);
fprintf('Target logic:    %s\n', optsRegioInvest.iterationLogicInvestment);
fprintf('Ref-yield in investment: %d\n', optsRegioInvest.RefYieldInInvestment);
fprintf('Sensitivity grid:        %d\n', ...
    optsRegioInvest.runReferenceYieldSensitivity);
fprintf('============================================================\n\n');

[resultsWindOnshoreTurbinesUnstack, ...
 regioDataWind, ...
 scenarioSummary] = ...
    calcRegioInvestWindOnshore();

%% ------------------------------------------------------------------------
% 7. OPTIONAL REGIONAL PLOTS
% -------------------------------------------------------------------------

plotResults = ...
    getConfig(config, 'plotResults', false);

if plotResults

    plot_results( ...
        resultsWindOnshoreTurbinesUnstack, ...
        regioDataWind);
end

%% ------------------------------------------------------------------------
% 8. SAVE REGIONAL RESULT
% -------------------------------------------------------------------------

saveRegionalResult = ...
    getConfig(config, 'saveRegionalResult', true);

resultsFolder = ...
    fullfile(repoPath, 'Results Paper');

if ~isfolder(resultsFolder)
    mkdir(resultsFolder);
end

resultFile = ...
    "";

if saveRegionalResult

    if optsRegioInvest.RefYieldInInvestment
        refYieldSuffix = '';
    else
        refYieldSuffix = '_no_reference_yield';
    end

    resultFilename = ...
        sprintf( ...
            'results_nuts3_base%d_target%d_expCase%d_%s%s.mat', ...
            paraRegioInvest.baseYear, ...
            paraRegioInvest.simYear, ...
            paraRegioInvest.expansionCase, ...
            optsRegioInvest.iterationLogicInvestment, ...
            refYieldSuffix);

    resultFile = ...
        fullfile(resultsFolder, resultFilename);

    selected_data = ...
        regioDataWind; %#ok<NASGU>

    scenarioConfig = ...
        config; %#ok<NASGU>

    save( ...
        resultFile, ...
        'selected_data', ...
        'resultsWindOnshoreTurbinesUnstack', ...
        'scenarioSummary', ...
        'scenarioConfig');

    fprintf('Regional results saved to:\n%s\n', resultFile);
end

end


%% ========================================================================
% LOCAL HELPER: CONFIG VALUE WITH DEFAULT
% ========================================================================

function value = getConfig(config, fieldName, defaultValue)

if isfield(config, fieldName)

    value = ...
        config.(fieldName);

else

    value = ...
        defaultValue;
end

end
