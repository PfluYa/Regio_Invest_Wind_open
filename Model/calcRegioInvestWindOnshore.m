function [resultsWindOnshoreTurbinesUnstack, regioDataWind] = calcRegioInvestWindOnshore()
% CALCREGIOINVESTWINDONSHORE
% Main model routine for regional onshore wind investment simulations.
%
% DESCRIPTION
% This function prepares all technology, market, land, and weather inputs,
% computes region-technology-specific profitability, and simulates future
% deployment under one of four expansion logics:
%
%   1 = Nested-logit / discrete-choice model (baseline methodology)
%   2 = Capacity scaled proportional to existing capacity
%   3 = Capacity allocated proportional to remaining land
%   4 = Merit-order allocation based on NPV
%
% OUTPUTS
%   resultsWindOnshoreTurbinesUnstack : wide table with estimated capacity
%                                       by NUTS3 region and turbine type
%   regioDataWind                     : regional summary table for plots/maps
%
% GLOBAL INPUTS
%   paraRegioInvest     : parameter struct from main_regio_invest.m
%   optsRegioInvest     : option struct from main_regio_invest.m
%   baseTableRegioInvest: regional base table with spatial information
%
% NOTES
% - This function relies on the global model containers set in the main
%   script.
% - Current application is tailored to Germany.
% - Default values referenced below correspond to the current settings in
%   main_regio_invest.m.

global paraRegioInvest baseTableRegioInvest optsRegioInvest

%% ------------------------------------------------------------------------
% 1. FILE DEFINITIONS
% -------------------------------------------------------------------------

filenameTurbineData        = 'windOnshoreTurbineData.csv';
filenameTurbinePowerCurves = 'windOnshoreTurbinePowerCurves.csv';
filenameDistSea            = 'Nuts3-Abstand.csv';

weatherYears = paraRegioInvest.weatherYear;
fileListWeatherYears = arrayfun( ...
    @(yearValue) sprintf('ERA5_data_windspeed_%d.nc', yearValue), ...
    weatherYears, ...
    'UniformOutput', false);

%% ------------------------------------------------------------------------
% 2. READ OPTIONS AND PARAMETERS
% -------------------------------------------------------------------------

% Reference-yield and investment options
optsReferenceYieldModel       = optsRegioInvest.windOnshoreReferenceYieldModel;   % default: true
optsUniqueRefYieldPerRegion   = optsRegioInvest.uniqueRefYieldperRegion;          % default: false
optsRefYieldInInvestment      = optsRegioInvest.RefYieldInInvestment;             % default: true
optsYieldCorrectionYP         = optsRegioInvest.windOnshoreYieldCorrectionYP;     % default: true
optsExplanatoryVariable       = optsRegioInvest.windOnshoreExplanatoryVariable;   % default: {'netPresentValue'}
optsNpvAdjust                 = optsRegioInvest.windOnshoreNpvAdjustment;         % default: 'NPVperArea'
optsExplainedVariable         = optsRegioInvest.windOnshoreExplainedVariable;     % default: {'relativeUsedWindSpaceDcm'}
optsScaleAreaParameter        = optsRegioInvest.scalingAreaParameter;             % default: 'UB21'
optsScaleAreaInvest           = optsRegioInvest.scalingAreaInvest;                % default: 'WindBG'
optsPowerPotentialIncrease    = optsRegioInvest.powerPotential_increase;          % default: true
optsPowerPotentialIndividual  = optsRegioInvest.powerPotential_individual;        % default: false

% Core parameters
paraTurbineAgeEstimation   = paraRegioInvest.windOnshoreTimeHorizonEstimation;    % default: 10
paraCompensation           = paraRegioInvest.windOnshoreCompensation;             % default: 0.080 EUR/kWh
paraInterestRate           = paraRegioInvest.windOnshoreInterestRate;             % default: 0.05
paraTurbineLifetime        = paraRegioInvest.windOnshoreTurbineLifetime;          % default: 22 years
paraOPEXperkWfix           = paraRegioInvest.opexPerkW_fix;                       % default: 15 EUR/kW
paraOPEXperkWhvar          = paraRegioInvest.opexPerkWh_var;                      % default: 0.0075 EUR/kWh
paraSimYear                = paraRegioInvest.simYear;                             % default: 2040
paraBaseYear               = paraRegioInvest.baseYear;                            % default: 2024
paraExpansionCase          = paraRegioInvest.expansionCase;                       % default: 1
paraInvestFactor           = 1;                                                   % default: 1

%% ------------------------------------------------------------------------
% 3. OPTIONAL TURBINE-SPECIFIC POWER POTENTIAL
% -------------------------------------------------------------------------
% If activated, overwrite the scalar power-potential assumption with
% technology-specific values. In the current default setup, this block is
% inactive because optsPowerPotentialIndividual = false.

if optsPowerPotentialIndividual
    % Current pattern used in the original code
    individualPattern = [23700 16600 26700 26000 26900 21800 31200 36300];

    paraRegioInvest.windOnshorePowerPotential = repmat(individualPattern.', 401, 1);
    paraRegioInvest.windOnshorePowerPotential_invest = ...
        paraRegioInvest.windOnshorePowerPotential;

    if optsPowerPotentialIncrease
        paraRegioInvest.windOnshorePowerPotential_invest = ...
            paraRegioInvest.windOnshorePowerPotential * 1.10;
    end
end

paraPowerPotential        = paraRegioInvest.windOnshorePowerPotential;         % default: 22500
paraPowerPotentialInvest  = paraRegioInvest.windOnshorePowerPotential_invest;  % default: 25000

%% ------------------------------------------------------------------------
% 4. DEFINE ESTIMATION AND SIMULATION HORIZONS
% -------------------------------------------------------------------------

paraCommissionedInSimYear = (paraSimYear - paraTurbineLifetime):paraSimYear;
paraTimeHorizonEstimation = (paraBaseYear - paraTurbineAgeEstimation):paraBaseYear;

%% ------------------------------------------------------------------------
% 5. LOAD MARKET AND TURBINE INPUTS
% -------------------------------------------------------------------------

marketData = getMarketData(paraBaseYear);
[turbineParameter, ~, turbinePowerCurves] = ...
    getTurbineData(filenameTurbineData, filenameTurbinePowerCurves);

marketData = clusterMarketData(marketData, turbineParameter, baseTableRegioInvest);

%% ------------------------------------------------------------------------
% 6. BUILD REGION-TECHNOLOGY PANEL
% -------------------------------------------------------------------------

regioWindOnshore = generateRegioWindonshore( ...
    baseTableRegioInvest, ...
    marketData, ...
    paraPowerPotential, ...
    paraTimeHorizonEstimation, ...
    paraCommissionedInSimYear, ...
    turbineParameter, ...
    optsScaleAreaParameter, ...
    optsScaleAreaInvest, ...
    paraExpansionCase);

paraOPEXfix = paraOPEXperkWfix .* regioWindOnshore.powerTurbine;

%% ------------------------------------------------------------------------
% 7. NATIONAL EXPANSION TARGET
% -------------------------------------------------------------------------

switch paraSimYear
    case 2030
        targetCapacityGermanyMW = 115000;
    case 2035
        targetCapacityGermanyMW = 157000;
    case 2040
        targetCapacityGermanyMW = 160000;
    otherwise
        error('Invalid simulation year. Supported values are 2030, 2035, and 2040.');
end

%% ------------------------------------------------------------------------
% 8. LOAD WEATHER DATA
% -------------------------------------------------------------------------

scenarioWeatherData = getWeatherData(regioWindOnshore, fileListWeatherYears);

%% ------------------------------------------------------------------------
% 9. RUN SELECTED EXPANSION CASE
% -------------------------------------------------------------------------

switch paraExpansionCase

    case 1
        %% =================================================================
        % CASE 1: NESTED-LOGIT / DISCRETE-CHOICE MODEL
        % ==================================================================

        % Step 1: Compute production and wind statistics for each
        % region-technology combination.
        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction( ...
                regioWindOnshore, ...
                turbineParameter, ...
                turbinePowerCurves, ...
                scenarioWeatherData);

        % Optional bias correction following Pflugfelder, Kramer, Weber (2024)
        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP( ...
                filenameDistSea, regioWindOnshore, turbineParameter);
        end

        % Step 2: Calculate reference-yield remuneration adjustment.
        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = calcTurbineCompensationFactor( ...
                regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor = ones(height(regioWindOnshore), 1);
        end

        % Step 3: Calculate NPV.
        [regioWindOnshore.netPresentValue, ~, ~] = calcNetPresentValue( ...
            regioWindOnshore.investCost, ...
            regioWindOnshore.production, ...
            paraTurbineLifetime, ...
            paraInterestRate, ...
            paraCompensation, ...
            regioWindOnshore.compensationFactor, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar);

        % Optional NPV normalization used in estimation/prediction.
        switch optsNpvAdjust
            case 'NPVperInvestCosts'
                regioWindOnshore.netPresentValue = ...
                    regioWindOnshore.netPresentValue ./ regioWindOnshore.investCost;

            case 'NPVperArea'
                regioWindOnshore.netPresentValue = ...
                    regioWindOnshore.netPresentValue ./ ...
                    (regioWindOnshore.powerTurbine ./ paraPowerPotential);

            case 'none'
                % Keep raw NPV.

            otherwise
                error('Invalid value for optsRegioInvest.windOnshoreNpvAdjustment.');
        end

        % Step 4: Prepare discrete-choice model inputs.
        choiceModelInput.explanatoryVariables = cellfun( ...
            @(varName) [varName, 'Dcm'], ...
            optsExplanatoryVariable, ...
            'UniformOutput', false);

        choiceModelInput.explainedVariable = optsExplainedVariable;
        choiceModelInput.nrAlternatives = size(turbineParameter, 1);

        % Transform explanatory variable for DCM estimation.
        regioWindOnshore.netPresentValueDcm = transformNPV(regioWindOnshore.netPresentValue);

        % Step 5: Estimate DCM for Germany.
        idxRelevantCountry = contains(regioWindOnshore.nutsID, 'DE');

        maxEstimationAttempts = 10;
        estimationSuccess = false;
        lastError = [];

        for tryCount = 1:maxEstimationAttempts
            try
                discreteChoiceModel = estimateDiscreteChoiceModel( ...
                    regioWindOnshore(idxRelevantCountry, :), ...
                    choiceModelInput);

                estimationSuccess = true;
                break;

            catch ME
                lastError = ME;
                fprintf('DCM estimation failed. Attempt %d of %d.\n', ...
                    tryCount, maxEstimationAttempts);
            end
        end

        if ~estimationSuccess
            rethrow(lastError);
        end

        % Step 6: Prepare technology/economic parameters for prediction.
        paraTech.paraCompensation          = paraCompensation;
        paraTech.paraInterestRate          = paraInterestRate;
        paraTech.paraTurbineLifetime       = paraTurbineLifetime;
        paraTech.paraInvestFactor          = paraInvestFactor;
        paraTech.paraPowerPotentialInvest  = paraPowerPotentialInvest;
        paraTech.paraPowerPotential        = paraPowerPotential;

        % Optional linearized reference-yield schedule.
        if optsRegioInvest.RefYieldLinearized
            regioWindOnshore.compensationFactor = ...
                calcTurbineCompensationFactor2( ...
                    regioWindOnshore, optsUniqueRefYieldPerRegion);
        end

        % Remove remuneration differentiation from the investment decision
        % if requested.
        if ~optsRefYieldInInvestment
            regioWindOnshore.compensationFactor = ones(height(regioWindOnshore), 1);
        end

        % Step 7: Predict future investment.
        [regioResults, ~, ~] = estimateInvestmentDecision( ...
            discreteChoiceModel, ...
            regioWindOnshore, ...
            targetCapacityGermanyMW, ...
            paraTech, ...
            optsScaleAreaInvest, ...
            optsNpvAdjust, ...
            optsRefYieldInInvestment, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar, ...
            paraRegioInvest, ...
            optsRegioInvest);

    case 2
        %% =================================================================
        % CASE 2: PROPORTIONAL TO EXISTING CAPACITY
        % ==================================================================

        currentCapacity = sum(regioWindOnshore.capacitiesBaseYear);
        scalingFactor = targetCapacityGermanyMW * 1000 / currentCapacity;

        regioWindOnshore.estimatedTotalCapacity = ...
            regioWindOnshore.capacitiesBaseYear .* scalingFactor;

        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction( ...
                regioWindOnshore, ...
                turbineParameter, ...
                turbinePowerCurves, ...
                scenarioWeatherData);

        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP( ...
                filenameDistSea, regioWindOnshore, turbineParameter);
        end

        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = calcTurbineCompensationFactor( ...
                regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor = ones(height(regioWindOnshore), 1);
        end

        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraBaseYear, ...
            paraSimYear);

        regioWindOnshore.netPresentValue = calcNetPresentValue( ...
            regioWindOnshore.investCost, ...
            regioWindOnshore.production, ...
            paraTurbineLifetime, ...
            paraInterestRate, ...
            paraCompensation, ...
            regioWindOnshore.compensationFactor, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar);

        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ regioWindOnshore.powerTurbine) .* ...
            regioWindOnshore.netPresentValue;

        fprintf('Total NPV of wind fleet: %.3f bn EUR\n', ...
            sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9);

        regioResults = regioWindOnshore;

    case 3
        %% =================================================================
        % CASE 3: PROPORTIONAL TO REMAINING LAND AVAILABILITY
        % ==================================================================

        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction( ...
                regioWindOnshore, ...
                turbineParameter, ...
                turbinePowerCurves, ...
                scenarioWeatherData);

        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP( ...
                filenameDistSea, regioWindOnshore, turbineParameter);
        end

        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = calcTurbineCompensationFactor( ...
                regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor = ones(height(regioWindOnshore), 1);
        end

        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraBaseYear, ...
            paraSimYear);

        regioWindOnshore.netPresentValue = calcNetPresentValue( ...
            regioWindOnshore.investCost, ...
            regioWindOnshore.production, ...
            paraTurbineLifetime, ...
            paraInterestRate, ...
            paraCompensation, ...
            regioWindOnshore.compensationFactor, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar);

        totalTargetCapacitykW = targetCapacityGermanyMW * 1000;
        nuts3Regions = unique(regioWindOnshore.nutsID);

        regionSummary = groupsummary(regioWindOnshore, 'nutsID', 'sum', 'capacitiesBaseYear');

        tmpAvailable = groupsummary(regioWindOnshore, 'nutsID', 'mean', 'availableWindSpaceTotal');
        regionSummary.mean_availableWindSpaceTotal = tmpAvailable.mean_availableWindSpaceTotal;

        regionSummary.usedArea = regionSummary.sum_capacitiesBaseYear ./ paraPowerPotential;

        isNearlyFull = ...
            regionSummary.usedArea > 0.9 .* regionSummary.mean_availableWindSpaceTotal;

        regionSummary.mean_availableWindSpaceTotal(isNearlyFull) = ...
            1.1 .* max( ...
                regionSummary.mean_availableWindSpaceTotal(isNearlyFull), ...
                regionSummary.usedArea(isNearlyFull));

        regionSummary.remainingArea = ...
            regionSummary.mean_availableWindSpaceTotal - regionSummary.usedArea;

        regionSummary.regionShare = ...
            regionSummary.remainingArea ./ sum(regionSummary.remainingArea);

        regionSummary.newCapPerRegion = ...
            regionSummary.regionShare .* ...
            (totalTargetCapacitykW - sum(regionSummary.sum_capacitiesBaseYear));

        newCapacity = zeros(height(regioWindOnshore), 1);

        for regionIdx = 1:numel(nuts3Regions)
            idxRegion = strcmp(regioWindOnshore.nutsID, nuts3Regions{regionIdx});
            regionNPV = regioWindOnshore.netPresentValue(idxRegion);

            regionNPVShare = regionNPV ./ sum(regionNPV);
            newCapacity(idxRegion) = ...
                regionNPVShare .* regionSummary.newCapPerRegion(regionIdx);
        end

        regioWindOnshore.estimatedTotalCapacity = ...
            newCapacity + regioWindOnshore.capacitiesBaseYear;

        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ regioWindOnshore.powerTurbine) .* ...
            regioWindOnshore.netPresentValue;

        fprintf('Total NPV of wind fleet: %.3f bn EUR\n', ...
            sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9);

        regioResults = regioWindOnshore;

    case 4
        %% =================================================================
        % CASE 4: MERIT-ORDER ALLOCATION BASED ON NPV
        % ==================================================================

        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction( ...
                regioWindOnshore, ...
                turbineParameter, ...
                turbinePowerCurves, ...
                scenarioWeatherData);

        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP( ...
                filenameDistSea, regioWindOnshore, turbineParameter);
        end

        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = calcTurbineCompensationFactor( ...
                regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor = ones(height(regioWindOnshore), 1);
        end

        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraBaseYear, ...
            paraSimYear);

        regioWindOnshore.usedArea = ...
            regioWindOnshore.capacitiesInSimYear ./ paraPowerPotential;

        isOverused = ...
            regioWindOnshore.usedArea > 0.9 .* regioWindOnshore.availableWindSpaceTotal;

        regioWindOnshore.availableWindSpaceTotal(isOverused) = ...
            1.1 .* max( ...
                regioWindOnshore.availableWindSpaceTotal(isOverused), ...
                regioWindOnshore.usedArea(isOverused));

        regioWindOnshore.remainingAreaFixed = ...
            regioWindOnshore.availableWindSpaceTotal - regioWindOnshore.usedArea;

        % Merit-order logic iterates on remuneration until the marginal
        % chosen project has approximately zero NPV.
        paraCompensation = 0.07;   % current hard-coded initial guess
        stepSize = 0.0005;         % current decrement
        toleranceNPV = 1e3;        % current tolerance
        maxIter = 1000;

        iter = 0;
        npvLast = 1e6;

        disp('Starting iterative optimization of remuneration level ...');

        while (npvLast > toleranceNPV) && (iter < maxIter)
            iter = iter + 1;

            regioWindOnshore.netPresentValue = calcNetPresentValue( ...
                regioWindOnshore.investCost, ...
                regioWindOnshore.production, ...
                paraTurbineLifetime, ...
                paraInterestRate, ...
                paraCompensation, ...
                regioWindOnshore.compensationFactor, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar);

            regioWindOnshore.remainingArea = regioWindOnshore.remainingAreaFixed;

            [~, sortedIndices] = sort(regioWindOnshore.netPresentValue, 'descend');
            regioWindOnshore = regioWindOnshore(sortedIndices, :);

            remainingCapacity = ...
                targetCapacityGermanyMW * 1000 - sum(regioWindOnshore.capacitiesInSimYear);

            regioWindOnshore.newCapPerRegion = zeros(height(regioWindOnshore), 1);

            if isscalar(paraPowerPotentialInvest)
                paraPowerPotentialInvest = repmat( ...
                    paraPowerPotentialInvest, height(regioWindOnshore), 1);
            end

            for rowIdx = 1:height(regioWindOnshore)
                if remainingCapacity <= 0
                    break;
                end

                availableCapacity = min( ...
                    remainingCapacity, ...
                    regioWindOnshore.remainingArea(rowIdx) .* paraPowerPotentialInvest(rowIdx));

                regioWindOnshore.newCapPerRegion(rowIdx) = availableCapacity;
                remainingCapacity = remainingCapacity - availableCapacity;

                % Prevent multiple allocations to the same region after the
                % first selected technology-region entry.
                sameRegion = strcmp( ...
                    regioWindOnshore.nutsID, regioWindOnshore.nutsID(rowIdx));
                regioWindOnshore.remainingArea(sameRegion) = 0;
            end

            idxLast = find(regioWindOnshore.newCapPerRegion > 0, 1, 'last');
            npvLast = -1e6;

            if ~isempty(idxLast)
                npvLast = regioWindOnshore.netPresentValue(idxLast);
            end

            fprintf(['Iteration %d: Remuneration = %.5f EUR/kWh, ', ...
                     'NPV of marginal project = %.0f EUR\n'], ...
                     iter, paraCompensation, npvLast);

            if npvLast > toleranceNPV
                paraCompensation = paraCompensation - stepSize;
            end
        end

        fprintf(['Converged: remuneration = %.5f EUR/kWh ', ...
                 '(marginal NPV approx. %.0f EUR)\n'], ...
                 paraCompensation, npvLast);

        regioWindOnshore.estimatedTotalCapacity = ...
            regioWindOnshore.newCapPerRegion + regioWindOnshore.capacitiesInSimYear;

        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ regioWindOnshore.powerTurbine) .* ...
            regioWindOnshore.netPresentValue;

        fprintf('Total NPV of wind fleet: %.3f bn EUR\n', ...
            sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9);

        regioWindOnshore.share_of_total_generation = ...
            (regioWindOnshore.estimatedTotalCapacity .* regioWindOnshore.production) ./ ...
            sum(regioWindOnshore.estimatedTotalCapacity .* regioWindOnshore.production);

        regioWindOnshore.weighted_comp_share = ...
            regioWindOnshore.share_of_total_generation .* ...
            regioWindOnshore.compensationFactor;

        realCompensation = sum(regioWindOnshore.weighted_comp_share) .* paraCompensation;

        fprintf('Generation-weighted average remuneration: %.3f ct/kWh\n', ...
            realCompensation * 100);

        regioResults = regioWindOnshore;

    otherwise
        error('Unsupported expansion case. Choose 1, 2, 3, or 4.');
end

%% ------------------------------------------------------------------------
% 10. FINALIZE OUTPUT TABLES
% -------------------------------------------------------------------------

resultsWindOnshoreTurbines = regioResults;

resultsWindOnshoreTurbinesUnstack = resultsWindOnshoreTurbines(:, ...
    {'nutsID', 'estimatedTotalCapacity', 'turbineType'});

resultsWindOnshoreTurbinesUnstack = unstack( ...
    resultsWindOnshoreTurbinesUnstack, ...
    'estimatedTotalCapacity', ...
    'turbineType');

resultsWindOnshoreTurbinesUnstack = sortrows( ...
    resultsWindOnshoreTurbinesUnstack, 'nutsID');

%% ------------------------------------------------------------------------
% 11. CREATE REGIONAL SUMMARY TABLE
% -------------------------------------------------------------------------

capPerNuts = sum(table2array(resultsWindOnshoreTurbinesUnstack(:, 2:end)), 2);

regioDataWind = baseTableRegioInvest;
regioDataWind.capacityTotal = capPerNuts;

regioDataWind.capPerKm2 = ...
    regioDataWind.capacityTotal ./ regioDataWind.totalArea;

regioDataWind.capPerKm2avail = ...
    regioDataWind.capacityTotal ./ ...
    (regioDataWind.relativeAvailableWindSpace .* regioDataWind.totalArea);

invalidEntries = isnan(regioDataWind.capPerKm2avail) | ...
                 isinf(regioDataWind.capPerKm2avail);

regioDataWind.capPerKm2avail(invalidEntries) = 0;

%% ------------------------------------------------------------------------
% 12. ADD BASE-YEAR AND EXPANSION METRICS
% -------------------------------------------------------------------------

capInSimYear = groupsummary(regioWindOnshore, "nutsID", "sum", "capacitiesInSimYear");
regioDataWind.capacityToAdd = ...
    regioDataWind.capacityTotal - capInSimYear.sum_capacitiesInSimYear;

capBaseYear = groupsummary(regioWindOnshore, "nutsID", "sum", "capacitiesBaseYear");
regioDataWind.capacity_baseYear = capBaseYear.sum_capacitiesBaseYear;

%% ------------------------------------------------------------------------
% 13. COMPUTE LAND EXHAUSTION METRIC
% -------------------------------------------------------------------------
% Baseline interpretation in the original code uses 22.5 MW/km^2.

regioDataWind.exhaustionProb = ...
    ((regioDataWind.capacityTotal ./ 1000) ./ 22.5) ./ ...
    (regioDataWind.relativeAvailableWindSpace .* regioDataWind.totalArea);

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function marketData = getMarketData(baseYear)
% GETMARKETDATA Load market-register wind turbine data for the selected
% base year.
%
% Supported base years in current code:
%   2020 -> 'MaStR_Wind20.xlsx'
%   2022 -> 'MaStR_Wind22.xlsx'
%   2024 -> 'MaStR_Wind24.xlsx'   (default current setting)

    switch baseYear
        case 2020
            filename = 'MaStR_Wind20.xlsx';
            turbines = readtable(filename);
            turbines = turbines(:, [1 21 8 11 12 13 15 16 17 18 19]);
            turbines = renamevars( ...
                turbines, ...
                ["Nuts3", "Inbetriebnahmedatum", "netPowerRating", "unitGroup"], ...
                ["nutsID", "commissionYear", "netPower", "turbineType"]);
            turbines.commissionYear = year(turbines.commissionYear);
            marketData = turbines;

        case 2022
            filename = 'MaStR_Wind22.xlsx';
            marketData = readtable(filename);

        case 2024
            % Public MaStR extract filtered to German onshore wind units.
            filename = 'MaStR_Wind24.xlsx';
            marketData = readtable(filename);

        otherwise
            error('Unsupported base year in getMarketData: %d', baseYear);
    end
end

function [outputTurbineData, outputTurbineCost, outputPowercurves] = ...
    getTurbineData(inpFilenameData, inpFilenamePowerCurves)
% GETTURBINEDATA Load turbine characteristics and power curves.

    turbineDataLocal = readtable(inpFilenameData);
    powerCurvesLocal = readtable(inpFilenamePowerCurves);

    outputTurbineData = turbineDataLocal;
    outputTurbineCost = turbineDataLocal(:, {'turbineType', 'costs'});
    outputPowercurves = powerCurvesLocal;
end

function outputWindMarket = clusterMarketData(inputWindMarket, inputTurbineData, inputBaseTable)
% CLUSTERMARKETDATA Assign each observed turbine to the closest turbine
% cluster based on hub height, rotor diameter, and net power.

    windMarket = inputWindMarket;
    clusterWind = inputTurbineData;

    [clusterWind{:, {'zHubHeight', 'zDiameter', 'zNetPowerRating'}}, mu, sigma] = ...
        zscore(clusterWind{:, {'hubHeight', 'diameter', 'power'}});

    centroids = clusterWind{:, {'zHubHeight', 'zDiameter', 'zNetPowerRating'}};

    windMarket(isnan(windMarket.hubHeight), :) = [];
    windMarket(isnan(windMarket.diameter), :) = [];

    if isdatetime(windMarket.commissionYear)
        windMarket.commissionYear = year(windMarket.commissionYear);
    end

    if ~ismember('nutsID', windMarket.Properties.VariableNames)
        windMarket.nutsID = strings(height(windMarket), 1);

        for regionIdx = 1:height(inputBaseTable)
            polyLat = inputBaseTable.latPoly{regionIdx};
            polyLon = inputBaseTable.lonPoly{regionIdx};

            inside = inpolygon(windMarket.lon, windMarket.lat, polyLon, polyLat);
            windMarket.nutsID(inside) = inputBaseTable.nutsID(regionIdx);
        end
    end

    windMarket(windMarket.nutsID == "", :) = [];

    windMarket.specPower = ...
        windMarket.netPower ./ (pi .* (windMarket.diameter ./ 2).^2);

    windMarketCentroids = windMarket;
    windMarketCentroids{:, {'zHubHeight', 'zDiameter', 'zNetPowerRating'}} = ...
        (windMarketCentroids{:, {'hubHeight', 'diameter', 'netPower'}} - mu) ./ sigma;

    clusterDistance = pdist2( ...
        centroids, ...
        windMarketCentroids{:, {'zHubHeight', 'zDiameter', 'zNetPowerRating'}});

    [~, idxClusterToMarket] = min(clusterDistance);

    windMarket.turbineType = clusterWind.turbineType(idxClusterToMarket, :);
    windMarket.cHubHeight  = clusterWind.hubHeight(idxClusterToMarket, :);
    windMarket.cDiameter   = clusterWind.diameter(idxClusterToMarket, :);
    windMarket.cPower      = clusterWind.power(idxClusterToMarket, :);
    windMarket.cInvest     = clusterWind.costs(idxClusterToMarket, :);

    outputWindMarket = windMarket;
end

function updatedInvestCost = applyCostReduction(investCost, annualReductionRate, baseYear, targetYear)
% APPLYCOSTREDUCTION Apply compound annual CAPEX reduction from base year
% to target year.
%
% Example current default:
% annualReductionRate = 0.0238 (i.e. 2.38% per year)

    yearsToTarget = targetYear - baseYear;
    reductionFactor = (1 - annualReductionRate) .^ yearsToTarget;
    updatedInvestCost = investCost .* reductionFactor;
end

function outputPastInvestments = generateRegioWindonshore( ...
        inputBaseTable, inputMarketData, inpPowerPotential, inpTimehorizon, ...
        inpCommissionedInSimYear, inpTurbineParameter, ...
        optsScaleAreaParameter, optsScaleAreaInvest, paraExpansionCase)
% GENERATEREGIOWINDONSHORE Construct the region-technology panel used for
% estimation and simulation.
%
% OUTPUT
% Long-format table with one row per region-turbine combination,
% including installed capacities, land availability, and DCM variables.

    baseTableLocal        = inputBaseTable;
    marketDataLocal       = inputMarketData;
    timeHorizon           = inpTimehorizon;
    commissionedSimYear   = inpCommissionedInSimYear;
    powerPotential        = inpPowerPotential;
    turbineParameterLocal = inpTurbineParameter;

    relevantCols = { ...
        'nutsID', 'countryCode', 'cityName', 'latPoint', 'lonPoint', ...
        'totalArea', 'relativeAvailableWindSpace'};
    baseTableLocal = baseTableLocal(:, relevantCols);

    [WindBG, BWE, UB21] = adaptPercentageWind(baseTableLocal);
    baseTableLocal.relativeAvailableWindSpace_WindBG = WindBG;
    baseTableLocal.relativeAvailableWindSpace_BWE    = BWE;
    baseTableLocal.relativeAvailableWindSpace_UB21   = UB21;

    nTurbines = size(turbineParameterLocal, 1);
    nRegions  = size(baseTableLocal, 1);

    localRegioWindonshore = repelem(baseTableLocal, nTurbines, 1);
    localRegioWindonshore.turbineType = ...
        repmat(turbineParameterLocal.turbineType, [nRegions, 1]);
    localRegioWindonshore.referenceProduction = ...
        repmat(turbineParameterLocal.referenceFiveYears ./ 5, [nRegions, 1]);
    localRegioWindonshore.investCost = ...
        repmat(turbineParameterLocal.costs .* turbineParameterLocal.power, [nRegions, 1]);
    localRegioWindonshore.powerTurbine = ...
        repmat(turbineParameterLocal.power, [nRegions, 1]);

    isInHorizon = ismember(marketDataLocal.commissionYear, timeHorizon);
    pastMarket = marketDataLocal(isInHorizon, :);

    isInSimYear = ismember(marketDataLocal.commissionYear, commissionedSimYear);
    futureMarket = marketDataLocal(isInSimYear, :);

    capBaseYear = groupsummary(marketDataLocal, {'nutsID', 'turbineType'}, 'sum', 'netPower');
    capInvestments = groupsummary(pastMarket, {'nutsID', 'turbineType'}, 'sum', 'netPower');
    capFuture = groupsummary(futureMarket, {'nutsID', 'turbineType'}, 'sum', 'netPower');

    capBaseYear = renamevars(capBaseYear, 'sum_netPower', 'capacitiesBaseYear');
    capInvestments = renamevars(capInvestments, 'sum_netPower', 'capacitiesDiscreteChoice');
    capFuture = renamevars(capFuture, 'sum_netPower', 'capacitiesInSimYear');

    localRegioWindonshore = outerjoin( ...
        localRegioWindonshore, capBaseYear, ...
        'Type', 'left', ...
        'Keys', {'nutsID', 'turbineType'}, ...
        'RightVariables', 'capacitiesBaseYear');

    localRegioWindonshore = outerjoin( ...
        localRegioWindonshore, capInvestments, ...
        'Type', 'left', ...
        'Keys', {'nutsID', 'turbineType'}, ...
        'RightVariables', 'capacitiesDiscreteChoice');

    localRegioWindonshore = outerjoin( ...
        localRegioWindonshore, capFuture, ...
        'Type', 'left', ...
        'Keys', {'nutsID', 'turbineType'}, ...
        'RightVariables', 'capacitiesInSimYear');

    localRegioWindonshore = fillmissing( ...
        localRegioWindonshore, ...
        'constant', 0, ...
        'DataVariables', { ...
            'capacitiesBaseYear', 'capacitiesDiscreteChoice', 'capacitiesInSimYear'});

    if paraExpansionCase == 4
        optsScaleAreaParameter = optsScaleAreaInvest;
    end

    switch optsScaleAreaParameter
        case 'none'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* ...
                localRegioWindonshore.relativeAvailableWindSpace;

        case 'BWE'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* ...
                localRegioWindonshore.relativeAvailableWindSpace_BWE;

        case 'WindBG'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* ...
                localRegioWindonshore.relativeAvailableWindSpace_WindBG;

        case 'UB21'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* ...
                localRegioWindonshore.relativeAvailableWindSpace_UB21;

        otherwise
            error('Invalid value for optsScaleAreaParameter: %s', optsScaleAreaParameter);
    end

    localRegioWindonshore.usedWindSpaceDcm = ...
        localRegioWindonshore.capacitiesDiscreteChoice ./ powerPotential;

    localRegioWindonshore.usedWindSpaceSim = ...
        localRegioWindonshore.capacitiesInSimYear ./ powerPotential;

    localRegioWindonshore.availableWindSpaceDcm = ...
        localRegioWindonshore.availableWindSpaceTotal - ...
        localRegioWindonshore.capacitiesBaseYear ./ powerPotential + ...
        localRegioWindonshore.capacitiesDiscreteChoice ./ powerPotential;

    localRegioWindonshore.availableWindSpaceSim = ...
        localRegioWindonshore.availableWindSpaceTotal - ...
        localRegioWindonshore.capacitiesInSimYear ./ powerPotential;

    localRegioWindonshore.relativeUsedWindSpaceDcm = ...
        localRegioWindonshore.usedWindSpaceDcm ./ ...
        localRegioWindonshore.availableWindSpaceDcm;

    invalidRelUse = isnan(localRegioWindonshore.relativeUsedWindSpaceDcm) | ...
                    isinf(localRegioWindonshore.relativeUsedWindSpaceDcm);

    localRegioWindonshore.relativeUsedWindSpaceDcm(invalidRelUse) = 0;

    overIdx = sum(localRegioWindonshore.relativeUsedWindSpaceDcm, 2) > 1;
    localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx, :) = ...
        localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx, :) ./ ...
        sum(localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx, :), 2);

    underIdx = sum(localRegioWindonshore.relativeUsedWindSpaceDcm, 2) < 0;
    localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx, :) = ...
        localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx, :) ./ ...
        sum(localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx, :), 2);

    localRegioWindonshore.relativeUsedWindSpaceSim = ...
        min(localRegioWindonshore.usedWindSpaceSim ./ ...
            localRegioWindonshore.availableWindSpaceTotal, 1);

    localRegioWindonshore = fillmissing( ...
        localRegioWindonshore, ...
        'constant', 0, ...
        'DataVariables', {'relativeUsedWindSpaceDcm', 'relativeUsedWindSpaceSim'});

    localRegioWindonshore.success = ceil( ...
        localRegioWindonshore.capacitiesDiscreteChoice ./ ...
        localRegioWindonshore.powerTurbine);

    localRegioWindonshore.fails = max( ...
        ceil(localRegioWindonshore.availableWindSpaceDcm .* powerPotential ./ ...
             localRegioWindonshore.powerTurbine) - ...
        localRegioWindonshore.success, ...
        0);

    localRegioWindonshore.trials = ...
        localRegioWindonshore.success + localRegioWindonshore.fails;

    outputPastInvestments = localRegioWindonshore;
end

function outputWeatherData = getWeatherData(regionalWindOnshore, fileListWeatherYears)
% GETWEATHERDATA Load ERA5 wind vectors and construct average annual
% duration curves for each location.
%
% ASSUMPTION
% Current implementation assumes 8760 hours per year and therefore ignores
% leap days.

    numYears = numel(fileListWeatherYears);
    numLocations = size(regionalWindOnshore, 1);

    totalWindSpeed100 = zeros(numLocations, 8760);
    totalWindSpeed10  = zeros(numLocations, 8760);

    for yearIdx = 1:numYears
        filename = fileListWeatherYears{yearIdx};

        windU100 = ncread(filename, 'u100');
        windV100 = ncread(filename, 'v100');
        windU10  = ncread(filename, 'u10');
        windV10  = ncread(filename, 'v10');
        latGrid  = ncread(filename, 'latitude');
        lonGrid  = ncread(filename, 'longitude');

        timeIdx = 1:8760;

        for locationIdx = 1:numLocations
            regLat = regionalWindOnshore.latPoint(locationIdx);
            regLon = regionalWindOnshore.lonPoint(locationIdx);

            latRounded = round(regLat / 0.25) * 0.25;
            lonRounded = round(regLon / 0.25) * 0.25;

            [~, latIdx] = min(abs(latGrid - latRounded));
            [~, lonIdx] = min(abs(lonGrid - lonRounded));

            windU100_series = squeeze(windU100(lonIdx, latIdx, timeIdx));
            windV100_series = squeeze(windV100(lonIdx, latIdx, timeIdx));
            windU10_series  = squeeze(windU10(lonIdx, latIdx, timeIdx));
            windV10_series  = squeeze(windV10(lonIdx, latIdx, timeIdx));

            windSpeed100 = sqrt(windU100_series.^2 + windV100_series.^2);
            windSpeed10  = sqrt(windU10_series.^2 + windV10_series.^2);

            if numYears > 1
                windSpeed100 = sort(windSpeed100, 'descend');
                windSpeed10  = sort(windSpeed10, 'descend');
            end

            totalWindSpeed100(locationIdx, :) = ...
                totalWindSpeed100(locationIdx, :) + windSpeed100';

            totalWindSpeed10(locationIdx, :) = ...
                totalWindSpeed10(locationIdx, :) + windSpeed10';
        end
    end

    averageWindSpeed100 = totalWindSpeed100 ./ numYears;
    averageWindSpeed10  = totalWindSpeed10 ./ numYears;

    weatherDataWindOnshore = regionalWindOnshore;
    weatherDataWindOnshore.Wind100 = num2cell(averageWindSpeed100, 2);
    weatherDataWindOnshore.Wind10  = num2cell(averageWindSpeed10, 2);

    for rowIdx = 1:height(weatherDataWindOnshore)
        if isrow(weatherDataWindOnshore.Wind100{rowIdx})
            weatherDataWindOnshore.Wind100{rowIdx} = ...
                weatherDataWindOnshore.Wind100{rowIdx}';
        end

        if isrow(weatherDataWindOnshore.Wind10{rowIdx})
            weatherDataWindOnshore.Wind10{rowIdx} = ...
                weatherDataWindOnshore.Wind10{rowIdx}';
        end
    end

    outputWeatherData = weatherDataWindOnshore;
end

function [outProduction, outFullLoadHours, outMeanWindspeed, outMeanWindspeedOnTurbine] = ...
    calcTurbineProduction(inpPastInvestments, inpTurbineParameter, inpTurbinePowerCurves, inpWeatherData)
% CALCTURBINEPRODUCTION Calculate expected annual production, full-load
% hours, and wind statistics for each region-turbine combination.

    windparkLocal    = inpPastInvestments;
    turbineData      = inpTurbineParameter;
    powerCurvesLocal = inpTurbinePowerCurves;
    weatherData      = inpWeatherData;

    uniqueRegions = unique(windparkLocal.nutsID);
    idxRelevantWeatherData = ismember(weatherData.nutsID, uniqueRegions);
    weatherData = weatherData(idxRelevantWeatherData, :);

    hubHeight = turbineData.hubHeight;
    nrTurbines = size(turbineData, 1);
    nrRegions = numel(uniqueRegions);
    nrTimestamps = size(weatherData.Wind10{1}, 1);

    capacity = repmat(turbineData.power, [1, nrRegions]);

    techAvail = 0.97;  % current default availability factor
    capfactor = 1;     % current default additional scaling factor

    hubHeightMatrix = repmat(hubHeight, nrRegions, nrTimestamps);
    hubHeightMatrix = reshape(hubHeightMatrix', [], 1);

    windSpeed10 = cell2mat(weatherData.Wind10);
    windSpeed100 = cell2mat(weatherData.Wind100);

    alpha = log(windSpeed100 ./ windSpeed10) ./ log(100 / 10);

    capacity = reshape(capacity, 1, nrTurbines, nrRegions);

    windSpeedOnHubheight = windSpeed10 .* ((hubHeightMatrix ./ 10) .^ alpha);
    windSpeedOnHubheight = reshape(windSpeedOnHubheight, nrTimestamps, nrTurbines, nrRegions);

    windSpeed = reshape(windSpeed10, nrTimestamps, nrTurbines, nrRegions);

    xsample = [0 3 4.5 7 11 15 25 35]';
    values  = [0 0 0.035 0.035 0.025 0.015 0 0]';

    scalefactor = interp1(xsample, values, windSpeedOnHubheight);
    windSpeedOnHubheightScaled = (1 - scalefactor) .* windSpeedOnHubheight;

    powerCurve = reshape(powerCurvesLocal.production, [], nrTurbines);
    powerCurveWindspeed = unique(powerCurvesLocal.windspeed);

    windProduction = zeros(size(windSpeedOnHubheight));

    for turbineIdx = 1:nrTurbines
        windProduction(:, turbineIdx, :) = ...
            interp1( ...
                powerCurveWindspeed, ...
                powerCurve(:, turbineIdx), ...
                windSpeedOnHubheightScaled(:, turbineIdx, :)) ./ ...
            max(powerCurve(:, turbineIdx)) .* ...
            capacity(:, turbineIdx, :) .* ...
            techAvail .* capfactor;
    end

    regioProduction = squeeze(sum(windProduction, 1));
    regioFullLoadHours = regioProduction ./ squeeze(capacity);

    outProduction = reshape(regioProduction, [], 1);
    outFullLoadHours = reshape(regioFullLoadHours, [], 1);

    outMeanWindspeed = reshape(squeeze(mean(windSpeed, 1)), [], 1);
    outMeanWindspeedOnTurbine = ...
        reshape(squeeze(mean(windSpeedOnHubheight, 1)), [], 1);
end

function outputCompensationFactor = calcTurbineCompensationFactor(inpPastInvestments, optsUniqueRefYieldPerRegion)
% CALCTURBINECOMPENSATIONFACTOR Calculate EEG reference-yield remuneration
% factors based on site productivity.
%
% INPUT
% optsUniqueRefYieldPerRegion:
%   true  -> use turbine type 4 reference yield for all turbine types in
%            a region
%   false -> use technology-specific values
%
% Default in current main script: false

    windparkLocal = inpPastInvestments;

    referenceFactor = ...
        windparkLocal.production ./ windparkLocal.referenceProduction;

    idxNonDE = ~contains(windparkLocal.nutsID, 'DE');
    referenceFactor(idxNonDE) = 1;

    NUTS3_Suedregion = { ...
        'DE145','DE121','DE146','DE112','DE147','DE132','DE12A','DE133','DE12B','DE113', ...
        'DE131','DE12C','DE114','DE125','DE11C','DE117','DE118','DE119','DE122','DE123', ...
        'DE138','DE139','DE115','DE11B','DE126','DE127','DE11D','DE129','DE124','DE148', ...
        'DE116','DE141','DE128','DE135','DE11A','DE136','DE149','DE111','DE142','DE137', ...
        'DE144','DE13A','DE143','DE275','DE214','DE231','DE234','DE251','DE256','DE261', ...
        'DE264','DE271','DE276','DE216','DE241','DE245','DE242','DE246','DE215','DE235', ...
        'DE217','DE224','DE277','DE22C','DE27D','DE218','DE219','DE21A','DE252','DE257', ...
        'DE248','DE21B','DE225','DE21C','DE253','DE258','DE21D','DE278','DE267','DE211', ...
        'DE272','DE226','DE273','DE268','DE21E','DE221','DE227','DE27A','DE26A','DE274', ...
        'DE21F','DE269','DE21G','DE212','DE21H','DE21I','DE236','DE25A','DE237','DE279', ...
        'DE254','DE259','DE27E','DE27B','DE222','DE228','DE21J','DE229','DE232','DE238', ...
        'DE213','DE21K','DE25B','DE22A','DE255','DE239','DE262','DE26B','DE21L','DE223', ...
        'DE22B','DE23A','DE21M','DE27C','DE233','DE21N','DE25C','DE263','DE26C','DE715', ...
        'DE711','DE716','DE717','DE71B','DE71C','DEB3B','DEB3C','DEB14','DEB22','DEB15', ...
        'DEB3D','DEB23','DEB31','DEB3E','DEB32','DEB3F','DEB3G','DEB33','DEB34','DEB35', ...
        'DEB3J','DEB36','DEB37','DEB1D','DEB3I','DEB38','DEB3H','DEB3K','DEB21','DEB25', ...
        'DEB39','DEB3A','DEC02','DEC03','DEC01','DEC04','DEC05','DEC06'};

    inSouth = ismember(windparkLocal.nutsID, NUTS3_Suedregion);

    x = [0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.50];
    y = [1.55 1.42 1.29 1.16 1.07 1.00 0.94 0.89 0.85 0.81 0.79];

    referenceFactor(inSouth) = max(referenceFactor(inSouth), 0.5);
    referenceFactor(~inSouth) = max(referenceFactor(~inSouth), 0.6);
    referenceFactor = min(referenceFactor, max(x));

    localCompensationFactor = interp1(x, y, referenceFactor);

    if optsUniqueRefYieldPerRegion
        isType4 = strcmp(string(windparkLocal.turbineType), "turbineType4");
        allRegions = unique(windparkLocal.nutsID);

        for regionIdx = 1:numel(allRegions)
            regionId = allRegions{regionIdx};
            idxRegion = strcmp(windparkLocal.nutsID, regionId);
            idxRegionType4 = idxRegion & isType4;

            if any(idxRegionType4)
                factorType4 = localCompensationFactor(find(idxRegionType4, 1, 'first'));
                localCompensationFactor(idxRegion) = factorType4;
            end
        end
    end

    outputCompensationFactor = localCompensationFactor;
end

function outputCompensationFactor = calcTurbineCompensationFactor2( ...
    inpPastInvestments, optsUniqueRefYieldPerRegion)
% CALCTURBINECOMPENSATIONFACTOR2 Linearized alternative to the default EEG
% reference-yield schedule.
%
% NOTE
% optsUniqueRefYieldPerRegion is accepted for interface compatibility, but
% currently not applied in this linearized version.

    if nargin < 2
        optsUniqueRefYieldPerRegion = false; %#ok<NASGU>
    end

    windparkLocal = inpPastInvestments;

    referenceFactor = ...
        windparkLocal.production ./ windparkLocal.referenceProduction;

    idxNonDE = ~contains(windparkLocal.nutsID, 'DE');
    referenceFactor(idxNonDE) = 1;

    x = [0.05 0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.50 3 5];

    xLeft = 0.50;
    yLeft = 1.55;
    xRef  = 1.00;
    yRef  = 1.00;

    slope = (yRef - yLeft) / (xRef - xLeft);
    y = yLeft + slope .* (x - xLeft);

    outputCompensationFactor = interp1(x, y, referenceFactor);
end

function choiceParameter = estimateDiscreteChoiceModel(inpInvestments, inpChoiceModelInput)
% ESTIMATEDISCRETECHOICEMODEL Estimate the discrete-choice model for wind
% investment decisions.

    investments = inpInvestments;
    choiceModelInput = inpChoiceModelInput;

    isGermany = contains(investments.nutsID, 'DE');
    investments = investments(isGermany, :);

    choiceModelInput.nrDecisionMaker = numel(unique(investments.nutsID));
    choiceModelInput.nrExplanatoryVariables = ...
        width(choiceModelInput.explanatoryVariables);

    choiceModelInput = initLogitParameter(choiceModelInput);

    disp('______________________________________')
    disp('Start parameter estimation')
    disp('--------------------------------------')

    choiceParameter = estimateLogitParameter(investments, choiceModelInput);
    choiceParameter.explanatoryVariables = choiceModelInput.explanatoryVariables;
end