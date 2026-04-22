%% REGIOINVEST WIND - MAIN SCRIPT
% RegioInvest allocates onshore wind capacity across regions for a given
% simulation year, accounting for alternative expansion cases and
% land-use restrictions.
%
% INPUT DATA
% Required input data must be placed in the repository's "Input_Data"
% folder. If provided as ZIP archives, unpack them before running this
% script.
%
% OUTPUT
% - Simulation results in MATLAB format (*.mat)
% - Result plots via plot_results(...)
%
% NOTES
% - This script is the main entry point of the model.
% - The model currently supports Germany as geographic scope.
% - Several downstream functions rely on the global variables declared
%   below. These globals should remain the only global variables used.

%% ------------------------------------------------------------------------
% 0. HOUSEKEEPING
% -------------------------------------------------------------------------

clc;
clear;
clear global;

%% ------------------------------------------------------------------------
% 1. PATH SETUP
% -------------------------------------------------------------------------
% Default (current) local repository parent path:
% 'C:\Users\Yannik.Pflugfelder\Documents\Github'
%
% Adjust this only if your local setup differs.

userpath = 'C:\Users\Yannik.Pflugfelder\Documents\Github';

repoName = 'Regio_Invest_Wind_open';
repoPath = fullfile(userpath, repoName);

if ~isfolder(repoPath)
    error(['Repository folder not found: ', repoPath, newline, ...
           'Please update "userpath" in main_regio_invest.m.']);
end

cd(repoPath);
addpath(genpath(repoPath));

%% ------------------------------------------------------------------------
% 2. GLOBAL MODEL CONTAINERS
% -------------------------------------------------------------------------
% These are the only intended global variables in the model.
% Do not introduce additional global variables in other scripts/functions.

global paraRegioInvest optsRegioInvest baseTableRegioInvest

%% ------------------------------------------------------------------------
% 3. MODEL YEAR CONFIGURATION
% -------------------------------------------------------------------------

% Weather years used to compute mean wind speed time series at each
% location.
% Default (current): [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024]
paraRegioInvest.weatherYear = 2015:2024;

% Base year for initializing the input data. Installations commissioned
% after this year are removed from the initial stock.
% Default (current): 2024
paraRegioInvest.baseYear = 2024;

% Target simulation year for the deployment simulation.
% Allowed values in current application: 2030, 2035, 2040
% Default (current): 2040
paraRegioInvest.simYear = 2040;

% Expansion scenario selector:
% 1 = Novel methodology (nested logit)
% 2 = Proportional to existing capacity (no area restriction)
% 3 = Proportional to remaining area (ignores wind yield)
% 4 = Merit-order principle (linear optimization, ignores logit)
% Default (current): 1
paraRegioInvest.expansionCase = 2;

% Geographic scope of the simulation.
% Default (current): {'DE'}
paraRegioInvest.geoScope = {'DE'};

% Mapping of geographic scope to input definition file.
% Default (current):
%   filename  = "geoScope.xlsx"
%   tablename = "DE"
optsRegioInvest.geoScope = table( ...
    "geoScope.xlsx", ...
    "DE", ...
    'VariableNames', {'filename', 'tablename'} ...
);

%% ------------------------------------------------------------------------
% 4. PREPARE BASE GEOGRAPHIC TABLE
% -------------------------------------------------------------------------
% Initializes the base regional table with geographic information such as
% NUTS identifiers, coordinates, polygons, area, and related attributes.
% Then adds the available space potential for wind deployment.

baseTableRegioInvest = getGeoShapeData();
baseTableRegioInvest = getSpacePotential(baseTableRegioInvest);

%% ------------------------------------------------------------------------
% 5. WIND ONSHORE MODEL OPTIONS
% -------------------------------------------------------------------------

% Apply the reference-yield model in remuneration calculations.
% Default (current): true
optsRegioInvest.windOnshoreReferenceYieldModel = true;

% Use one common reference yield per region for all turbine types.
% If false, reference yield may differ by turbine type.
% Default (current): false
optsRegioInvest.uniqueRefYieldperRegion = false;

% Include reference-yield correction already in the investment decision.
% Default (current): true
optsRegioInvest.RefYieldInInvestment = true;

% Use a linearized version of the reference-yield mechanism.
% Default (current): false
optsRegioInvest.RefYieldLinearized = false;

% Apply additional yield correction after Pflugfelder et al. style logic.
% Marked as alternative / experimental in original script.
% Default (current): true
optsRegioInvest.windOnshoreYieldCorrectionYP = true;

% Methodology used for the onshore wind investment model.
% Default (current): 'binomialNestedLogit'
optsRegioInvest.windOnshoreMethodology = 'binomialNestedLogit';

% Iteration logic used to match policy targets:
% 'capacity' = target installed capacity
% 'energy'   = target electricity generation
% Default (current): 'capacity'
optsRegioInvest.iterationLogicInvestment = 'capacity';

% Explanatory variable(s) for the discrete choice model.
% Default (current): {'netPresentValue'}
optsRegioInvest.windOnshoreExplanatoryVariable = {'netPresentValue'};

% Transformation / normalization of NPV in the investment model.
% Options used in code base: 'NPVperArea', 'NPVperInvestCosts', 'none'
% Default (current): 'NPVperArea'
optsRegioInvest.windOnshoreNpvAdjustment = 'NPVperArea';

% Explained variable in the discrete choice model.
% Default (current): {'relativeUsedWindSpaceDcm'}
optsRegioInvest.windOnshoreExplainedVariable = {'relativeUsedWindSpaceDcm'};

%% ------------------------------------------------------------------------
% 6. LAND-SCALING AND AREA ASSUMPTIONS
% -------------------------------------------------------------------------

% Calibration of available area in the base year.
% 'UB21' refers to calibration to designated areas in 2021
% (Bons et al.).
% Default (current): 'UB21'
optsRegioInvest.scalingAreaParameter = 'UB21';

% Area scaling used in future investment simulations.
% 'WindBG' applies the WindBG target logic.
% Alternative mentioned in original script: 'none'
% Default (current): 'WindBG'
optsRegioInvest.scalingAreaInvest = 'WindBG';

%% ------------------------------------------------------------------------
% 7. FUTURE POWER-POTENTIAL ASSUMPTIONS
% -------------------------------------------------------------------------

% Allow higher power density in future years.
% Default (current): true
optsRegioInvest.powerPotential_increase = true;

% Use turbine-specific power potential values.
% Default (current): false
optsRegioInvest.powerPotential_individual = false;

%% ------------------------------------------------------------------------
% 8. ECONOMIC PARAMETERS
% -------------------------------------------------------------------------

% Number of historical years used for parameter estimation.
% Default (current): 10 years
paraRegioInvest.windOnshoreTimeHorizonEstimation = 10;

% Initial remuneration level in EUR/kWh.
% Default (current): 0.080 EUR/kWh
paraRegioInvest.windOnshoreCompensation = 0.080;

% Discount / interest rate.
% Default (current): 0.05
paraRegioInvest.windOnshoreInterestRate = 0.05;

% Current power potential in W/km^2-equivalent model units as used in the
% original code base.
% Default (current): 22500
paraRegioInvest.windOnshorePowerPotential = 22500;

% Future power potential used in investment simulations.
% Default (current): 25000
paraRegioInvest.windOnshorePowerPotential_invest = 25000;

% Technical/economic turbine lifetime in years.
% Default (current): 22
paraRegioInvest.windOnshoreTurbineLifetime = 22;

% Fixed operating expenditures in EUR/kW.
% Default (current): 15
paraRegioInvest.opexPerkW_fix = 15;

% Variable operating expenditures in EUR/kWh.
% Default (current): 0.0075
paraRegioInvest.opexPerkWh_var = 0.0075;

% Annual cost reduction rate for future turbine investment costs.
% Default (current): 0.0238, i.e. 2.38% per year
paraRegioInvest.costReductionRate = 0.0238;

%% ------------------------------------------------------------------------
% 9. RUN WIND ONSHORE MODEL
% -------------------------------------------------------------------------

disp('Calculating onshore wind investments ...');

[resultsWindOnshoreTurbinesUnstack, regioDataWind] = ...
    calcRegioInvestWindOnshore();

%% ------------------------------------------------------------------------
% 10. PLOT RESULTS
% -------------------------------------------------------------------------

plot_results(resultsWindOnshoreTurbinesUnstack, regioDataWind);

%% ------------------------------------------------------------------------
% 11. SAVE RESULTS
% -------------------------------------------------------------------------

currentDate = datestr(now, 'yyyy-mm-dd');

if optsRegioInvest.RefYieldInInvestment
    resultSuffix = '';
else
    resultSuffix = '_wo_refYield';
end

resultFilename = sprintf( ...
    '%s_results_nuts3_base%d_target%d_expCase%d%s.mat', ...
    currentDate, ...
    paraRegioInvest.baseYear, ...
    paraRegioInvest.simYear, ...
    paraRegioInvest.expansionCase, ...
    resultSuffix ...
);

resultsFolder = fullfile(repoPath, 'Results Paper');

if ~isfolder(resultsFolder)
    mkdir(resultsFolder);
end

fullResultPath = fullfile(resultsFolder, resultFilename);

% Save the full regional simulation output.
selected_data = regioDataWind; 
save(fullResultPath, 'selected_data');

fprintf('Results saved to:\n%s\n', fullResultPath);

%% ------------------------------------------------------------------------
% 12. OPTIONAL: CREATE SHAPEFILES FOR EXTERNAL MAPPING
% -------------------------------------------------------------------------
% Uncomment if shapefile export is needed.
%
% resultScriptShapefile()