%% *RegioInvestTool Wind*
% _RegioInvest is a tool to allocate onshore wind capacity to regions
% in a given simulation year, accounting for multiple expansion cases
% and land-use restrictions._
%
% Note on input data:
% Locate the INPUT DATA at XXXXX. Unpack ZIP File in the Input_Data Folder

%% --- User-defined local path ---
userpath = 'C:\Users\Yannik.Pflugfelder\Documents\Github'; % Change this to your local path if necessary
cd([userpath, '\Regio_Invest_Wind\']);
addpath(genpath([userpath, '\Regio_Invest_Wind\']))

%% --- Clean workspace ---
clc
clear
clear global

%% --- General Settings & Global Parameters ---
global paraRegioInvest optsRegioInvest baseTableRegioInvest
% NOTE: These are the *only* global variables. Do not define others globally.

%% --- Model Year Configuration ---
% Years used to compute mean wind speed time series at each location
paraRegioInvest.weatherYear = [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024]; % [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024]

% Base year for initializing input data (removes newer installations)
paraRegioInvest.baseYear = 2024;

% Target simulation year for capacity deployment
% Options: 2030, 2035, 2040
paraRegioInvest.simYear = 2040;

% Expansion scenario (choose 1–4):
% 1: Novel methodology (nested logit)
% 2: Proportional to existing capacity (no area restriction)
% 3: Proportional to remaining area (ignores wind yield)
% 4: Merit-order principle (linear optimization, ignores logit)
paraRegioInvest.expansionCase = 1;

% Geographic scope of simulation (currently only Germany)
paraRegioInvest.geoScope = {'DE'};
optsRegioInvest.geoScope = table("geoScope.xlsx", "DE", 'VariableNames', {'filename', 'tablename'});


%% Prepare Base Table RegioInvest - Get Geo Structure and Information from Shapefiles of GeoScope and Geo Resolution, Add Socio Economic Data and Free SpaceArea for each Technology

baseTableRegioInvest = getGeoShapeData(); % This Function Inits the BaseTable with GeoData -> Nuts ID's, Lat Lons, LatPoly, LanPoly, Area, CityName etc.
baseTableRegioInvest = getSpacePotential(baseTableRegioInvest);

%% --- Wind Onshore Model Settings ---

optsRegioInvest.windOnshoreReferenceYieldModel = true;    % Apply reference yield model
optsRegioInvest.uniqueRefYieldperRegion = false;           % Take same ref yield (from type 4) for all types within a region
optsRegioInvest.RefYieldInInvestment = true;              % Include yield correction in investment decision
optsRegioInvest.windOnshoreYieldCorrectionYP = true;     % Alternative correction approach (experimental)
optsRegioInvest.windOnshoreMethodology = 'binomialNestedLogit';

% Explanatory variable for the discrete choice model
optsRegioInvest.windOnshoreExplanatoryVariable = {'netPresentValue'};
optsRegioInvest.windOnshoreNpvAdjustment = 'NPVperArea';        % Options: 'NPVperArea', 'NPVperInvestCosts', 'none'
optsRegioInvest.windOnshoreExplainedVariable = {'relativeUsedWindSpaceDcm'};

% Scaling settings for regional area allocation
optsRegioInvest.scalingAreaParameter = 'UB21';            % Calibrate to 2021 designated areas (UB21, Bons et al.)
optsRegioInvest.scalingAreaInvest = 'WindBG';             % WindBG law target (2% land area by 2032)

% Future assumptions about power potential
optsRegioInvest.powerPotential_increase = true;           % Allow higher density in future years
optsRegioInvest.powerPotential_individual = false;        % Individual PP per turbine type (if needed)

% Economic settings
paraRegioInvest.windOnshoreTimeHorizonEstimation = 10;    % Years of past data for parameter estimation
paraRegioInvest.windOnshoreCompensation = 0.080;          % Initial remuneration level (€/kWh)
paraRegioInvest.windOnshoreInterestRate = 0.05;          % Discount rate
paraRegioInvest.windOnshorePowerPotential = 22500;        % MW/km² (current)
paraRegioInvest.windOnshorePowerPotential_invest = 25000; % MW/km² (future)
paraRegioInvest.windOnshoreTurbineLifetime = 22;          % Years
paraRegioInvest.opexPerkW_fix = 15;                       % fix OPEX per kW
paraRegioInvest.opexPerkWh_var = 00.0075;                   % var OPEX per kWh
paraRegioInvest.costReductionRate = 0.0238;               % % Apply cost reductions in line with IRENA (2024): approx. -2.38% annually
%% Run the calculations for Wind Onshore

disp('Calculating Investments Wind Onshore ...')
[resultsWindOnshoreTurbinesUnstack, regioDataWind] = calcRegioInvestWindOnshore();

%% Plot results
%plot_results(resultsWindOnshoreTurbinesUnstack, regioDataWind);

%% Save Results
current_date = datestr(now, 'yyyy-mm-dd');

if optsRegioInvest.RefYieldInInvestment
    suffix = '';
else
    suffix = '_wo_refYield';
end

filename = sprintf('%s_results_nuts3_base%d_target%d_expCase%d%s.mat', ...
    current_date, ...
    paraRegioInvest.baseYear, ...
    paraRegioInvest.simYear, ...
    paraRegioInvest.expansionCase, suffix);

% Define the full file path
filepath = fullfile(cd, 'Results Paper');
full_filepath = fullfile(filepath, filename);

% Extract the desired columns
selected_data = regioDataWind;

% Save the selected data in a MAT file
save(full_filepath, 'selected_data');

%% Create shapefiles for futher maps via GitHub
% resultScriptShapefile()

