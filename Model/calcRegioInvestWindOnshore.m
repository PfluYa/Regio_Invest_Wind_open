function [resultsWindOnshoreTurbinesUnstack, regioDataWind] = calcRegioInvestWindOnshore()
% CALCREGIOINVESTWINDONSHORE
% Main function for regional wind onshore investment modeling

global paraRegioInvest baseTableRegioInvest optsRegioInvest

%% 1. Initialize Options and Parameters

% === Filenames ===
filenameTurbineData         = 'windOnshoreTurbineData.csv';
filenameTurbinePowerCurves  = 'windOnshoreTurbinePowerCurves.csv';

% Weather years as NetCDF files
fileListWeatherYears = cell(1, length(paraRegioInvest.weatherYear));
for i = 1:length(paraRegioInvest.weatherYear)
    year1 = paraRegioInvest.weatherYear(i);
    fileListWeatherYears{i} = sprintf('ERA5_data_windspeed_%d.nc', year1);
end

filenameDistSea = 'Nuts3-Abstand.csv';

% === Options/Settings ===
optsReferenceYieldModel      = optsRegioInvest.windOnshoreReferenceYieldModel;
optsUniqueRefYieldPerRegion  = optsRegioInvest.uniqueRefYieldperRegion;
optsRefYieldInInvestment     = optsRegioInvest.RefYieldInInvestment;
optsYieldCorrectionYP        = optsRegioInvest.windOnshoreYieldCorrectionYP;
optsExplanatoryVariable      = optsRegioInvest.windOnshoreExplanatoryVariable;
optsNpvAdjust                = optsRegioInvest.windOnshoreNpvAdjustment;
optsExplainedVariable        = optsRegioInvest.windOnshoreExplainedVariable;
optsScaleAreaParameter       = optsRegioInvest.scalingAreaParameter;
optsScaleAreaInvest          = optsRegioInvest.scalingAreaInvest;
optsPowerPotential_increase  = optsRegioInvest.powerPotential_increase;
optsPowerPotential_individual= optsRegioInvest.powerPotential_individual;

% === Key Parameters ===
paraTurbineAgeEstimation     = paraRegioInvest.windOnshoreTimeHorizonEstimation;
paraCompensation             = paraRegioInvest.windOnshoreCompensation;
paraInterestRate             = paraRegioInvest.windOnshoreInterestRate;

% Optional: override power potential with type-specific values
if optsPowerPotential_individual
    pattern = [23700 16600 26700 26000 26900 21800 31200 36300];
    paraRegioInvest.windOnshorePowerPotential = repmat(pattern.', 401, 1);
    paraRegioInvest.windOnshorePowerPotential_invest = paraRegioInvest.windOnshorePowerPotential;
    if optsPowerPotential_increase
        paraRegioInvest.windOnshorePowerPotential_invest = paraRegioInvest.windOnshorePowerPotential * 1.10;
    end

end

paraPowerPotential           = paraRegioInvest.windOnshorePowerPotential;
paraPowerPotential_invest    = paraRegioInvest.windOnshorePowerPotential_invest;

paraTurbineLifetime          = paraRegioInvest.windOnshoreTurbineLifetime;
paraOPEXperkWfix             = paraRegioInvest.opexPerkW_fix;
paraOPEXperkWhvar            = paraRegioInvest.opexPerkWh_var;
paraSimYear                  = paraRegioInvest.simYear;
paraBaseYear                 = paraRegioInvest.baseYear;

% Investment scaling factor (typically 1)
paraInvestFactor             = 1;

% Commissioning years for the simulation
paraCommissionedInSimYear    = (paraSimYear - paraTurbineLifetime):paraSimYear;
paraTimeHorizonEstimation    = paraBaseYear - paraTurbineAgeEstimation : paraBaseYear;
paraExpansionCase            = paraRegioInvest.expansionCase;

%% 2. Load and cluster turbine data
marketData = getMarketData(paraBaseYear);
[turbineParameter, ~, turbinePowerCurves] = getTurbineData(filenameTurbineData, filenameTurbinePowerCurves);
marketData = clusterMarketData(marketData, turbineParameter, baseTableRegioInvest);

%% 3. Generate regional investment dataset
regioWindOnshore = generateRegioWindonshore( ...
    baseTableRegioInvest, ...
    marketData, ...
    paraPowerPotential, ...
    paraTimeHorizonEstimation, ...
    paraCommissionedInSimYear, ...
    turbineParameter, ...
    optsScaleAreaParameter, ...
    optsScaleAreaInvest, ...
    paraExpansionCase ...
);

% Fixed OPEX per turbine
paraOPEXfix = paraOPEXperkWfix .* regioWindOnshore.powerTurbine;

%% 4. Set national wind expansion target for Germany (EEG scenario)
switch paraRegioInvest.simYear
    case 2030
        targetCapacity_Germany_MW = 115000;
    case 2035
        targetCapacity_Germany_MW = 157000;
    case 2040
        targetCapacity_Germany_MW = 160000;
    otherwise
        error('Invalid simulation year! Only 2030, 2035 or 2040 are supported.');
end

scenarioWeatherData = getWeatherData(regioWindOnshore, fileListWeatherYears);

%% 5. Choose calculation logic depending on expansion scenario
switch paraExpansionCase
    case 1 % Nested-logit model (Discrete Choice Model / DCM)
        %% ================================
        %  Part 1: Parameter Estimation
        % ================================
        
        % 1. Calculate the expected energy production for all region/turbine combinations
        %    This will be the input for the choice model
        [regioWindOnshore.production, regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction(regioWindOnshore, turbineParameter, turbinePowerCurves, scenarioWeatherData);
        % Apply bias correction from Pflugfelder, Kramer, Weber (2024)
        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter);
        end

        % 2. Compensation Factor Adjustment
        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = calcTurbineCompensationFactor(regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor(:) = 1;
        end
        
        % 3. Calculate Net Present Value (NPV) of each turbine-region combination
        [regioWindOnshore.netPresentValue] = calcNetPresentValue( ...
            regioWindOnshore.investCost, ...
            regioWindOnshore.production, ...
            paraTurbineLifetime, ...
            paraInterestRate, ...
            paraCompensation, ...
            regioWindOnshore.compensationFactor, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar);
        
        % Optionally adjust NPV per investment cost or area
        switch optsNpvAdjust
            case 'NPVperInvestCosts'
                regioWindOnshore.netPresentValue = regioWindOnshore.netPresentValue ./ regioWindOnshore.investCost;
            case 'NPVperArea'
                regioWindOnshore.netPresentValue = regioWindOnshore.netPresentValue ./ ...
                    (regioWindOnshore.powerTurbine ./ paraPowerPotential);
            case 'none'
                % Do nothing
            otherwise
                disp('Invalid entry for optsNpvAdjust');
        end
        
        % 4. Prepare inputs for the discrete choice model estimation
        choiceModelInput.explanatoryVariables = ...
            cellfun(@(x) [x, 'Dcm'], optsExplanatoryVariable, 'UniformOutput', false);
        choiceModelInput.explainedVariable = optsExplainedVariable;
        choiceModelInput.nrAlternatives = size(turbineParameter, 1);
        
        % Scale explanatory variables (log scaling)
        regioWindOnshore.netPresentValueDcm = transformNPV(regioWindOnshore.netPresentValue);
        
        %% Estimate the Discrete Choice Model
        dcmRelevantCountry = 'DE';
        idxDCMRelevantCountry = contains(regioWindOnshore.nutsID, dcmRelevantCountry);
        tryCount = 1;
        estimationSuccess = false;

        while tryCount < 10 && ~estimationSuccess
            try
                while ~estimationSuccess
                    discreteChoiceModel = estimateDiscreteChoiceModel( ...
                        regioWindOnshore(idxDCMRelevantCountry,:), choiceModelInput);

                    userInput = input('Type "r" to repeat DCM or "c" to continue: ', 's');
                    if strcmp(userInput, 'r')
                        continue;
                    elseif strcmp(userInput, 'c')
                        estimationSuccess = true;
                    else
                        disp('Invalid input, continuing anyway...');
                        estimationSuccess = true;
                    end
                end
            catch
                disp(['Could not estimate DCM. Attempt Nr: ', num2str(tryCount)]);
                tryCount = tryCount + 1;
            end
        end

        %% ================================
        %  Part 2: Predict Investment Decision
        % ================================
        % Setup technical parameters for prediction
        paraTech.paraCompensation = paraCompensation;
        paraTech.paraInterestRate = paraInterestRate;
        paraTech.paraTurbineLifetime = paraTurbineLifetime;
        paraTech.paraInvestFactor = paraInvestFactor;
        
        % Choose power potential parameter
        paraTech.paraPowerPotential_invest = paraPowerPotential_invest;
        paraTech.paraPowerPotential = paraPowerPotential;
        
        % Reset compensation factor if reference yield model not used in investment
        if ~optsRefYieldInInvestment
            regioWindOnshore.compensationFactor(:) = 1;
        end
        
        [regioResults, newInvestments, regionalCapacities] = ...
            estimateInvestmentDecision( ...
                discreteChoiceModel, ...
                regioWindOnshore, ...
                targetCapacity_Germany_MW, ...
                paraTech, ...
                optsScaleAreaInvest, ...
                optsNpvAdjust, ...
                optsRefYieldInInvestment, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar, paraRegioInvest);
    case 2  % Case 2: Allocation proportional to existing capacity
    
        %% === Step 1: Scale current capacity to match scenario target ===
        % Calculate the total installed capacity of the base year
        currentCap = sum(regioWindOnshore.capacitiesBaseYear);
    
        % Determine the scaling factor to reach the target scenario capacity (e.g., 160 GW)
        scalingFac = targetCapacity_Germany_MW * 1000 / currentCap;
    
        % Estimate future capacity for each turbine-region pair by scaling base year capacity
        regioWindOnshore.estimatedTotalCapacity = ...
            regioWindOnshore.capacitiesBaseYear * scalingFac;
    
        %% === Step 2: Calculate energy production and wind metrics ===
        % Estimate energy production, full load hours, and wind speed statistics
        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction(regioWindOnshore, ...
                                  turbineParameter, ...
                                  turbinePowerCurves, ...
                                  scenarioWeatherData);
        % Apply bias correction from Pflugfelder, Kramer, Weber (2024)
        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter);
        end
    
        %% === Step 3: Determine compensation factor based on policy ===
        if optsReferenceYieldModel
            % Use EEG reference yield model correction factor
            regioWindOnshore.compensationFactor = ...
                calcTurbineCompensationFactor(regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            % No correction (uniform remuneration)
            regioWindOnshore.compensationFactor(:) = 1;
        end
    
        %% === Step 4: Apply CAPEX reductions for future years ===
        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraRegioInvest.baseYear, ...
            paraRegioInvest.simYear);
    
        %% === Step 5: Compute Net Present Value (NPV) for each unit ===
        regioWindOnshore.netPresentValue = ...
            calcNetPresentValue( ...
                regioWindOnshore.investCost, ...
                regioWindOnshore.production, ...
                paraTurbineLifetime, ...
                paraInterestRate, ...
                paraCompensation, ...
                regioWindOnshore.compensationFactor, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar);
    
        %% === Step 6: Aggregate total NPV across regions and types ===
        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ ...
             regioWindOnshore.powerTurbine) .* ...
             regioWindOnshore.netPresentValue;
    
        % Display total investor returns (NPV in bn EUR)
        disp(['The total NPV of the wind power fleet is: ', ...
              num2str(sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9), ...
              ' bn EUR'])
    
        %% === Step 7: Output result structure ===
        regioResults = regioWindOnshore;
    case 3  % Case 3: Allocation proportional to remaining land availability
  
        %% === Step 1: Estimate energy production and wind characteristics ===
        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction(regioWindOnshore, ...
                                  turbineParameter, ...
                                  turbinePowerCurves, ...
                                  scenarioWeatherData);
        
        % Apply bias correction from Pflugfelder, Kramer, Weber (2024)
        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter);
        end
    
        %% === Step 2: Apply policy-dependent compensation correction ===
        if optsReferenceYieldModel
            % Apply EEG reference yield model
            regioWindOnshore.compensationFactor = ...
                calcTurbineCompensationFactor(regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            % Uniform remuneration across regions
            regioWindOnshore.compensationFactor(:) = 1;
        end
    
        %% === Step 3: Adjust investment cost based on cost learning ===
        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraRegioInvest.baseYear, ...
            paraRegioInvest.simYear);
    
        %% === Step 4: Calculate Net Present Value (NPV) ===
        regioWindOnshore.netPresentValue = ...
            calcNetPresentValue( ...
                regioWindOnshore.investCost, ...
                regioWindOnshore.production, ...
                paraTurbineLifetime, ...
                paraInterestRate, ...
                paraCompensation, ...
                regioWindOnshore.compensationFactor, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar);
    
        %% === Step 5: Distribute new capacity based on remaining land ===
        totalNewCapacity = targetCapacity_Germany_MW * 1000;  % Convert MW to kW
        NUTS3_Regions = unique(regioWindOnshore.nutsID);
    
        % Summarize base year capacity and average available area by region
        regionSummary = groupsummary(regioWindOnshore, 'nutsID', 'sum', 'capacitiesBaseYear');
        regionSummary.mean_availableWindSpaceTotal = ...
            groupsummary(regioWindOnshore, 'nutsID', 'mean', 'availableWindSpaceTotal').mean_availableWindSpaceTotal;
    
        % Calculate area already used for installations
        regionSummary.usedArea = regionSummary.sum_capacitiesBaseYear / paraPowerPotential;
    
        % If region is almost full (>90%), increase usable space slightly to allow expansion
        bufferIndex = regionSummary.usedArea > 0.9 * regionSummary.mean_availableWindSpaceTotal;
        regionSummary.mean_availableWindSpaceTotal(bufferIndex) = ...
            1.1 * max(regionSummary.mean_availableWindSpaceTotal(bufferIndex), ...
                      regionSummary.usedArea(bufferIndex));
    
        % Calculate remaining area available per region
        regionSummary.remainingArea = ...
            regionSummary.mean_availableWindSpaceTotal - regionSummary.usedArea;
    
        % Calculate share of each region in national expansion based on available land
        regionSummary.regionShare = ...
            regionSummary.remainingArea ./ sum(regionSummary.remainingArea);
    
        regionSummary.newCapPerRegion = ...
            regionSummary.regionShare * (totalNewCapacity - sum(regionSummary.sum_capacitiesBaseYear));
    
        %% === Step 6: Distribute regional capacity by turbine NPV ===
        newCapacity = zeros(height(regioWindOnshore), 1);
    
        for i = 1:length(NUTS3_Regions)
            regionIndices = find(strcmp(regioWindOnshore.nutsID, NUTS3_Regions{i}));
            regionNPV = regioWindOnshore.netPresentValue(regionIndices);
            regionNPVShare = regionNPV / sum(regionNPV);
            newCapacity(regionIndices) = regionNPVShare * regionSummary.newCapPerRegion(i);
        end
    
        % Final capacity = base year + new allocation
        regioWindOnshore.estimatedTotalCapacity = ...
            newCapacity + regioWindOnshore.capacitiesBaseYear;
    
        %% === Step 7: Calculate total system NPV ===
        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ regioWindOnshore.powerTurbine) .* ...
             regioWindOnshore.netPresentValue;
    
        disp(['The total NPV of the wind power fleet is: ', ...
              num2str(sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9), ...
              ' bn EUR'])
    
        %% === Step 8: Return structured result ===
        regioResults = regioWindOnshore;
    case 4  % Case 4: Merit-order-based expansion by regional Net Present Value (NPV)
    
        %% === Step 1: Estimate production and wind characteristics ===
        [regioWindOnshore.production, ...
         regioWindOnshore.fullLoadHours, ...
         regioWindOnshore.meanWindspeed, ...
         regioWindOnshore.meanWindspeedOnTurbine] = ...
            calcTurbineProduction(regioWindOnshore, ...
                                  turbineParameter, ...
                                  turbinePowerCurves, ...
                                  scenarioWeatherData);
        
        % Apply bias correction from Pflugfelder, Kramer, Weber (2024)
        if optsYieldCorrectionYP
            regioWindOnshore.production = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter);
        end
    
        %% === Step 2: Apply policy-based compensation correction ===
        if optsReferenceYieldModel
            regioWindOnshore.compensationFactor = ...
                calcTurbineCompensationFactor(regioWindOnshore, optsUniqueRefYieldPerRegion);
        else
            regioWindOnshore.compensationFactor(:) = 1;
        end
    
        %% === Step 3: Apply CAPEX reduction over time ===
        regioWindOnshore.investCost = applyCostReduction( ...
            regioWindOnshore.investCost, ...
            paraRegioInvest.costReductionRate, ...
            paraRegioInvest.baseYear, ...
            paraRegioInvest.simYear);
    
        %% === Step 4: Estimate land usage and remaining area ===
        regioWindOnshore.usedArea = ...
            regioWindOnshore.capacitiesInSimYear ./ paraPowerPotential;
    
        % Buffer adjustment for nearly full regions
        overused = regioWindOnshore.usedArea > ...
                   0.9 * regioWindOnshore.availableWindSpaceTotal;
    
        regioWindOnshore.availableWindSpaceTotal(overused) = ...
            1.1 * max(regioWindOnshore.availableWindSpaceTotal(overused), ...
                      regioWindOnshore.usedArea(overused));
    
        % Calculate fixed remaining area
        regioWindOnshore.remainingAreaFixed = ...
            regioWindOnshore.availableWindSpaceTotal - regioWindOnshore.usedArea;
    
        %% === Step 5: Iteratively adjust compensation level ===
        paraCompensation = 0.07;       % Initial guess [€/kWh]
        stepSize = 0.0005;             % Compensation decrement per iteration
        toleranceNPV = 1e3;            % NPV target threshold
        maxIter = 1000;
        iter = 0;
        npvLast = 1e6;
    
        disp('Starting iterative optimization of compensation level...');
    
        while (npvLast > toleranceNPV) && (iter < maxIter)
            iter = iter + 1;
    
            % Step 5.1: Recalculate NPV with current compensation
            regioWindOnshore.netPresentValue = ...
                calcNetPresentValue(regioWindOnshore.investCost, ...
                                    regioWindOnshore.production, ...
                                    paraTurbineLifetime, ...
                                    paraInterestRate, ...
                                    paraCompensation, ...
                                    regioWindOnshore.compensationFactor, ...
                                    paraOPEXfix, paraOPEXperkWhvar);
    
            % Step 5.2: Sort sites by NPV (Merit Order)
            regioWindOnshore.remainingArea = regioWindOnshore.remainingAreaFixed;
            [~, sortedIndices] = sort(regioWindOnshore.netPresentValue, 'descend');
            regioWindOnshore = regioWindOnshore(sortedIndices, :);
    
            % Step 5.3: Allocate capacity starting from highest NPV
            remainingCapacity = targetCapacity_Germany_MW * 1000 - ...
                                sum(regioWindOnshore.capacitiesInSimYear);
    
            regioWindOnshore.newCapPerRegion = zeros(height(regioWindOnshore), 1);
    
            if isscalar(paraPowerPotential_invest)
                paraPowerPotential_invest = repmat(paraPowerPotential_invest, height(regioWindOnshore), 1);
            end

            for i = 1:height(regioWindOnshore)
                if remainingCapacity <= 0
                    break;
                end
                availableCap = min(remainingCapacity, ...
                                   regioWindOnshore.remainingArea(i) .* paraPowerPotential_invest(i));
    
                regioWindOnshore.newCapPerRegion(i) = availableCap;
                remainingCapacity = remainingCapacity - availableCap;
    
                % Prevent reuse of this region
                regioWindOnshore.remainingArea(strcmp(regioWindOnshore.nutsID, ...
                    regioWindOnshore.nutsID(i))) = 0;
            end
    
            % Step 5.4: Check NPV of last-added site
            idxLast = find(regioWindOnshore.newCapPerRegion > 0, 1, 'last');
            npvLast = -1e6;
            if ~isempty(idxLast)
                npvLast = regioWindOnshore.netPresentValue(idxLast);
            end
    
            fprintf('Iteration %d: Compensation = %.5f €/kWh, Last NPV = %.0f EUR\n', ...
                     iter, paraCompensation, npvLast);
    
            % Step 5.5: Adjust compensation
            if npvLast > toleranceNPV
                paraCompensation = paraCompensation - stepSize;
            end
        end
    
        fprintf('Converged: Optimal compensation = %.5f €/kWh (Last NPV ≈ %.0f EUR)\n', ...
                paraCompensation, npvLast);
    
        %% === Step 6: Final output and indicators ===
        regioWindOnshore.estimatedTotalCapacity = ...
            regioWindOnshore.newCapPerRegion + regioWindOnshore.capacitiesInSimYear;
    
        regioWindOnshore.NPVofTypeAndRegion = ...
            (regioWindOnshore.estimatedTotalCapacity ./ regioWindOnshore.powerTurbine) .* ...
            regioWindOnshore.netPresentValue;
    
        disp(['Total NPV of wind fleet: ', ...
              num2str(sum(regioWindOnshore.NPVofTypeAndRegion) / 1e9), ' bn EUR']);
    
        % Calculate generation-weighted average compensation
        regioWindOnshore.share_of_total_generation = ...
            (regioWindOnshore.estimatedTotalCapacity .* regioWindOnshore.production) / ...
             sum(regioWindOnshore.estimatedTotalCapacity .* regioWindOnshore.production);
    
        regioWindOnshore.weighted_comp_share = ...
            regioWindOnshore.share_of_total_generation .* ...
            regioWindOnshore.compensationFactor;
    
        real_compensation = sum(regioWindOnshore.weighted_comp_share) * paraCompensation;
    
        disp(['Real average remuneration (weighted): ', ...
              num2str(real_compensation * 100), ' ct/kWh']);
    
        %% Return
        regioResults = regioWindOnshore;
end
%% =======================
% Finalize Simulation Outputs
% =======================

% Store full results for output
resultsWindOnshoreTurbines = regioResults;

% Unstack estimated capacities by turbine type
resultsWindOnshoreTurbinesUnstack = resultsWindOnshoreTurbines(:, ...
    {'nutsID','estimatedTotalCapacity','turbineType'});
resultsWindOnshoreTurbinesUnstack = ...
    unstack(resultsWindOnshoreTurbinesUnstack, ...
            'estimatedTotalCapacity', 'turbineType');
resultsWindOnshoreTurbinesUnstack = ...
    sortrows(resultsWindOnshoreTurbinesUnstack, 'nutsID');

%% =======================
% Create Regional Summary Table for Mapping and Plots
% =======================

% Sum all turbine-type capacities per region (in same order as base table)
capPerNuts = sum(table2array(resultsWindOnshoreTurbinesUnstack(:, 2:end)), 2);

% Copy base geographic data
regioDataWind = baseTableRegioInvest;

% Assign total capacity per region
regioDataWind.capacityTotal = capPerNuts;

% Compute density metrics
regioDataWind.capPerKm2 = regioDataWind.capacityTotal ./ regioDataWind.totalArea;
regioDataWind.capPerKm2avail = ...
    regioDataWind.capacityTotal ./ ...
    (regioDataWind.relativeAvailableWindSpace .* regioDataWind.totalArea);

% Clean up NaNs or Infs caused by regions with zero available space
invalidEntries = isnan(regioDataWind.capPerKm2avail) | ...
                 isinf(regioDataWind.capPerKm2avail);
regioDataWind.capPerKm2avail(invalidEntries) = 0;

%% =======================
% Calculate Expansion and Base Year Metrics
% =======================

% Determine installed capacity to be added in each region
capInSimYear = groupsummary(regioWindOnshore, "nutsID", "sum", "capacitiesInSimYear");
regioDataWind.capacityToAdd = ...
    regioDataWind.capacityTotal - capInSimYear.sum_capacitiesInSimYear;

% Store base year capacities
capBaseYear = groupsummary(regioWindOnshore, "nutsID", "sum", "capacitiesBaseYear");
regioDataWind.capacity_baseYear = capBaseYear.sum_capacitiesBaseYear;

%% =======================
% Compute Land Exhaustion Probability
% =======================

% Exhaustion = (installed capacity / potential) / (available area fraction × total area)
% Power potential assumed to be 22.5 MW/km² as baseline (can be adapted)
regioDataWind.exhaustionProb = ...
    ((regioDataWind.capacityTotal ./ 1000) ./ 22.5) ./ ...
    (regioDataWind.relativeAvailableWindSpace .* regioDataWind.totalArea);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%      Local Functions      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION GET MARKETDATA
function marketData = getMarketData(baseYear)
    if baseYear == 2020
        filename = 'MaStR_Wind20.xlsx';                          
        turbines = readtable(filename);
        turbines = turbines(:,[1 21 8 11 12 13 15 16 17 18 19]);
        turbines = renamevars(turbines, ["Nuts3", "Inbetriebnahmedatum", "netPowerRating", "unitGroup"], ["nutsID", "commissionYear", "netPower", "turbineType"]);
        turbines.commissionYear = year(turbines.commissionYear);
        marketData = turbines;
    elseif baseYear == 2022
        filename = 'MaStR_Wind22.xlsx';
        turbines = readtable(filename);
        marketData = turbines;
    elseif baseYear == 2024
        % https://www.marktstammdatenregister.de/MaStR/Einheit/Einheiten/ErweiterteOeffentlicheEinheitenuebersicht?filter=Installierte%20Leistung%20der%20EEG-Anlage~gt~%27100%27~and~Betriebs-Status~eq~%2735%27~and~Inbetriebnahmedatum%20der%20Einheit~lt~%2701.01.2025%27~and~Lage%20der%20Einheit~eq~%27888%27~and~Energietr%C3%A4ger~eq~%272497%27
        filename = 'MaStR_Wind24.xlsx';
        turbines = readtable(filename);
        marketData = turbines;
    end
end
%% FUNCTION GETTURBINEDATA
function [outputTurbineData,outputTurbineCost,outputPowercurves] = getTurbineData(inpFilenameData,inpFilenamePowerCurves)
    % init Filenames and number of Turbines
    filenameData = inpFilenameData;
    filenamePowerCurves = inpFilenamePowerCurves;
    turbineDataLocal = readtable(filenameData);
    powerCurvesLocal = readtable(filenamePowerCurves);
    
    % extract cost from Turbine Data
    turbineCost = turbineDataLocal(:,{'turbineType','costs'});
    
    % write OutputData
    outputTurbineData = turbineDataLocal;
    outputTurbineCost = turbineCost;
    outputPowercurves = powerCurvesLocal;
end
%% FUNCTION CLUSTERMARKETDATA
function [outputWindMarket] = clusterMarketData(inputWindMarket, inputTurbineData, inputBaseTable)
    %CLUSTERMARKETDATA Assigns each wind turbine in the market data to a representative cluster.
    %
    %   INPUTS:
    %     - inputWindMarket   : table of existing turbines (location, specs, net power, etc.)
    %     - inputTurbineData  : predefined turbine clusters with hub height, diameter, power, costs
    %     - inputBaseTable    : spatial geo reference with polygon geometry and NUTS3 IDs
    %
    %   OUTPUT:
    %     - outputWindMarket  : market table enriched with turbineType, cluster match, and NUTS ID
    
    %% Step 1: Initialize inputs
    windMarket = inputWindMarket;
    clusterWind = inputTurbineData;
    
    %% Step 2: Normalize turbine cluster features (Z-score scaling)
    [clusterWind{:,{'zHubHeight','zDiameter','zNetPowerRating'}}, mu, sigma] = ...
        zscore(clusterWind{:,{'hubHeight','diameter','power'}});
    
    centroids = clusterWind{:,{'zHubHeight','zDiameter','zNetPowerRating'}};
    
    %% Step 3: Clean up incomplete market entries
    windMarket(isnan(windMarket.hubHeight), :) = [];
    windMarket(isnan(windMarket.diameter), :) = [];
    
    % Convert datetime to numeric year if necessary
    if isdatetime(windMarket.commissionYear)
        windMarket.commissionYear = year(windMarket.commissionYear);   
    end
    
    %% Step 4: Assign NUTS3 ID based on geographic polygon
    if ~ismember('nutsID', windMarket.Properties.VariableNames)
        windMarket.nutsID = strings(height(windMarket), 1);
    
        for j = 1:height(inputBaseTable)
            polyLat = inputBaseTable.latPoly{j};  % Latitude vertices
            polyLon = inputBaseTable.lonPoly{j};  % Longitude vertices
    
            % Find turbines within polygon
            inside = inpolygon(windMarket.lon, windMarket.lat, polyLon, polyLat);
    
            % Assign NUTS3 ID to matched turbines
            windMarket.nutsID(inside) = inputBaseTable.nutsID(j);
        end
    end
    
    % Remove turbines without valid NUTS3 assignment
    windMarket(windMarket.nutsID == "", :) = [];
    
    %% Step 5: Calculate specific power
    windMarket.specPower = windMarket.netPower ./ (pi * (windMarket.diameter / 2).^2);
    
    %% Step 6: Normalize market turbine data with same mean & std
    windMarketCentroids = windMarket;
    windMarketCentroids{:,{'zHubHeight','zDiameter','zNetPowerRating'}} = ...
        (windMarketCentroids{:,{'hubHeight','diameter','netPower'}} - mu) ./ sigma;
    
    %% Step 7: Match each market turbine to the nearest cluster (Euclidean distance)
    clusterDistance = pdist2(centroids, windMarketCentroids{:,{'zHubHeight','zDiameter','zNetPowerRating'}});
    [~, idxClusterToMarket] = min(clusterDistance);  % Closest cluster index
    
    %% Step 8: Assign cluster values back to market table
    windMarket.turbineType   = clusterWind.turbineType(idxClusterToMarket, :);
    windMarket.cHubHeight    = clusterWind.hubHeight(idxClusterToMarket, :);
    windMarket.cDiameter     = clusterWind.diameter(idxClusterToMarket, :);
    windMarket.cPower        = clusterWind.power(idxClusterToMarket, :);
    windMarket.cInvest       = clusterWind.costs(idxClusterToMarket, :);
    
    %% Output result
    outputWindMarket = windMarket;

end
%% FUNCTION UPDATE (REDUCE) INVEST COSTS FOR FUTURE YEARS
function updatedInvestCost = applyCostReduction(investCost, annualReductionRate, baseYear, targetYear)
    yearsToTarget = targetYear - baseYear;
    reductionFactor = (1 - annualReductionRate) .^ yearsToTarget;
    updatedInvestCost = investCost .* reductionFactor;
end


%% FUNCTION GENERATEREGIOWINDONSHORE
function outputPastInvestments = generateRegioWindonshore( ...
        inputBaseTable, inputMarketData, inpPowerPotential, inpTimehorizon, ...
        inpCommissionedInSimYear, inpTurbineParameter, ...
        optsScaleAreaParameter, optsScaleAreaInvest, paraExpansionCase)
    %GENERATEREGIOWINDONSHORE Generates regional wind investment base table
    %
    % Returns a long-format table of wind technology-region combinations 
    % with added capacity information, land potential, and DCM input values.
    %
    % INPUTS:
    %   - inputBaseTable:       base region information incl. NUTS3 and area
    %   - inputMarketData:      historic market dataset of wind turbines
    %   - inpPowerPotential:    max installable capacity per km² (in kW/km²)
    %   - inpTimehorizon:       years used for past investment estimation
    %   - inpCommissionedInSimYear:  year to evaluate future capacity
    %   - inpTurbineParameter:  turbine cluster specifications
    %   - optsScaleAreaParameter:  source to scale land availability (e.g., 'none', 'UB21')
    %   - optsScaleAreaInvest:      option for scaling (used in expansion case 4)
    %   - paraExpansionCase:    scenario number (1–4)
    
    % === Initialize local variables ===
    baseTableLocal        = inputBaseTable;
    marketDataLocal       = inputMarketData;
    timeHorizon           = inpTimehorizon;
    commissionedSimYear   = inpCommissionedInSimYear;
    powerPotential        = inpPowerPotential; % in kW/km²
    turbineParameterLocal = inpTurbineParameter;
    
    % === Filter base table to essential columns ===
    relevantCols = {'nutsID','countryCode','cityName','latPoint','lonPoint','totalArea','relativeAvailableWindSpace'};
    baseTableLocal = baseTableLocal(:, relevantCols);
    
    % === Compute multiple land scaling options ===
    [WindBG, BWE, UB21] = adaptPercentageWind(baseTableLocal);
    baseTableLocal.relativeAvailableWindSpace_WindBG = WindBG;
    baseTableLocal.relativeAvailableWindSpace_BWE    = BWE;
    baseTableLocal.relativeAvailableWindSpace_UB21   = UB21;
    
    % === Expand baseTable to all turbine-region combinations ===
    nTurbines = size(turbineParameterLocal, 1);
    nRegions  = size(baseTableLocal, 1);
    
    localRegioWindonshore                     = repelem(baseTableLocal, nTurbines, 1);
    localRegioWindonshore.turbineType         = repmat(turbineParameterLocal.turbineType, [nRegions, 1]);
    localRegioWindonshore.referenceProduction = repmat(turbineParameterLocal.referenceFiveYears / 5, [nRegions, 1]);
    localRegioWindonshore.investCost          = repmat(turbineParameterLocal.costs .* turbineParameterLocal.power, [nRegions, 1]);
    localRegioWindonshore.powerTurbine        = repmat(turbineParameterLocal.power, [nRegions, 1]);
    
    % === Filter market data for past and future ===
    isInHorizon      = ismember(marketDataLocal.commissionYear, timeHorizon);
    pastMarket       = marketDataLocal(isInHorizon, :);
    
    isInSimYear      = ismember(marketDataLocal.commissionYear, commissionedSimYear);
    futureMarket     = marketDataLocal(isInSimYear, :);
    
    % === Aggregate installed capacity ===
    capBaseYear      = groupsummary(marketDataLocal, {'nutsID','turbineType'}, 'sum', 'netPower');
    capInvestments   = groupsummary(pastMarket, {'nutsID','turbineType'}, 'sum', 'netPower');
    capFuture        = groupsummary(futureMarket, {'nutsID','turbineType'}, 'sum', 'netPower');
    
    capBaseYear      = renamevars(capBaseYear, 'sum_netPower', 'capacitiesBaseYear');
    capInvestments   = renamevars(capInvestments, 'sum_netPower', 'capacitiesDiscreteChoice');
    capFuture        = renamevars(capFuture, 'sum_netPower', 'capacitiesInSimYear');
    
    % === Join capacity data to base table ===
    localRegioWindonshore = outerjoin(localRegioWindonshore, capBaseYear, 'Type','left', 'Keys',{'nutsID','turbineType'}, 'RightVariables','capacitiesBaseYear');
    localRegioWindonshore = outerjoin(localRegioWindonshore, capInvestments, 'Type','left', 'Keys',{'nutsID','turbineType'}, 'RightVariables','capacitiesDiscreteChoice');
    localRegioWindonshore = outerjoin(localRegioWindonshore, capFuture, 'Type','left', 'Keys',{'nutsID','turbineType'}, 'RightVariables','capacitiesInSimYear');
    
    localRegioWindonshore = fillmissing(localRegioWindonshore, 'constant', 0, ...
        'DataVariables', {'capacitiesBaseYear', 'capacitiesDiscreteChoice', 'capacitiesInSimYear'});
    
    % === Determine available area based on chosen expansion scenario ===
    if paraExpansionCase == 4
        optsScaleAreaParameter = optsScaleAreaInvest;
    end
    
    switch optsScaleAreaParameter
        case 'none'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* localRegioWindonshore.relativeAvailableWindSpace;
        case 'BWE'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* localRegioWindonshore.relativeAvailableWindSpace_BWE;
        case 'WindBG'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* localRegioWindonshore.relativeAvailableWindSpace_WindBG;
        case 'UB21'
            localRegioWindonshore.availableWindSpaceTotal = ...
                localRegioWindonshore.totalArea .* localRegioWindonshore.relativeAvailableWindSpace_UB21;
        otherwise
            warning('Invalid value for optsScaleAreaParameter');
    end
    
    % === Compute area usage and availability ===
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
    
    % === Calculate relative land usage ===
    localRegioWindonshore.relativeUsedWindSpaceDcm = ...
        localRegioWindonshore.usedWindSpaceDcm ./ localRegioWindonshore.availableWindSpaceDcm;
    
    localRegioWindonshore.relativeUsedWindSpaceDcm( ...
        isnan(localRegioWindonshore.relativeUsedWindSpaceDcm) | ...
        isinf(localRegioWindonshore.relativeUsedWindSpaceDcm)) = 0;
    
    % Correct cases with values > 1 or < 0
    overIdx = sum(localRegioWindonshore.relativeUsedWindSpaceDcm,2) > 1;
    localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx,:) = ...
        localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx,:) ./ ...
        sum(localRegioWindonshore.relativeUsedWindSpaceDcm(overIdx,:), 2);
    
    underIdx = sum(localRegioWindonshore.relativeUsedWindSpaceDcm,2) < 0;
    localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx,:) = ...
        localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx,:) ./ ...
        sum(localRegioWindonshore.relativeUsedWindSpaceDcm(underIdx,:), 2);
    
    % === Future year land usage ===
    localRegioWindonshore.relativeUsedWindSpaceSim = ...
        min(localRegioWindonshore.usedWindSpaceSim ./ localRegioWindonshore.availableWindSpaceTotal, 1);
    
    localRegioWindonshore = fillmissing(localRegioWindonshore, 'constant', 0, ...
        'DataVariables', {'relativeUsedWindSpaceDcm', 'relativeUsedWindSpaceSim'});
    
    % === Calculate trial statistics for discrete choice model ===
    localRegioWindonshore.success = ...
        ceil(localRegioWindonshore.capacitiesDiscreteChoice ./ localRegioWindonshore.powerTurbine);
    
    localRegioWindonshore.fails = max( ...
        ceil(localRegioWindonshore.availableWindSpaceDcm .* powerPotential ./ localRegioWindonshore.powerTurbine) - ...
        localRegioWindonshore.success, 0);
    
    localRegioWindonshore.trials = ...
        localRegioWindonshore.success + localRegioWindonshore.fails;
    
    % === Output final table ===
    outputPastInvestments = localRegioWindonshore;

end

%% FUNCTION GETWEATHERDATA
%% FUNCTION: getWeatherData
function outputWeatherData = getWeatherData(regionalWindOnshore, fileListWeatherYears)
    % This function loads hourly wind speed data from ERA5 NetCDF files for 
    % given geolocations and years, then computes the average annual 
    % wind profile (duration curve) at 100m and 10m hub heights.

    % Number of input years and geolocations
    numYears = length(fileListWeatherYears); 
    numLocations = size(regionalWindOnshore, 1);

    % Preallocate wind speed matrices for 100m and 10m heights (non-leap years assumed)
    totalWindSpeed100 = zeros(numLocations, 8760); 
    totalWindSpeed10 = zeros(numLocations, 8760); 

    % Loop over weather years
    for yearIdx = 1:numYears
        filename = fileListWeatherYears{yearIdx};

        % Load wind vector components from NetCDF
        WIND_u100 = ncread(filename, 'u100');
        WIND_v100 = ncread(filename, 'v100');
        WIND_u10 = ncread(filename, 'u10');
        WIND_v10 = ncread(filename, 'v10');
        LAT = ncread(filename, 'latitude');
        LON = ncread(filename, 'longitude');

        % Assume standard year (non-leap) with 8760 time steps
        timeIdx = 1:8760;

        % Loop through all NUTS3 regions or coordinates
        for ct_regio = 1:numLocations
            reglat = regionalWindOnshore.latPoint(ct_regio);
            reglon = regionalWindOnshore.lonPoint(ct_regio);

            % Snap to nearest ERA5 grid
            latRounded = round(reglat / 0.25) * 0.25;
            lonRounded = round(reglon / 0.25) * 0.25;

            % Find nearest indices
            [~, latIdx] = min(abs(LAT - latRounded));
            [~, lonIdx] = min(abs(LON - lonRounded));

            % Extract time series for this point
            WindU_100 = squeeze(WIND_u100(lonIdx, latIdx, timeIdx));
            WindV_100 = squeeze(WIND_v100(lonIdx, latIdx, timeIdx));
            WindU_10 = squeeze(WIND_u10(lonIdx, latIdx, timeIdx));
            WindV_10 = squeeze(WIND_v10(lonIdx, latIdx, timeIdx));

            % Calculate scalar wind speed
            windSpeed100 = sqrt(WindU_100.^2 + WindV_100.^2);
            windSpeed10 = sqrt(WindU_10.^2 + WindV_10.^2);

            % Sort to get duration curves (if averaging over multiple years)
            if numYears > 1
                windSpeed100 = sort(windSpeed100, 'descend');
                windSpeed10 = sort(windSpeed10, 'descend');
            end

            % Accumulate
            totalWindSpeed100(ct_regio, :) = totalWindSpeed100(ct_regio, :) + windSpeed100';
            totalWindSpeed10(ct_regio, :) = totalWindSpeed10(ct_regio, :) + windSpeed10';
        end
    end

    % Compute average wind profiles
    averageWindSpeed100 = totalWindSpeed100 / numYears;
    averageWindSpeed10 = totalWindSpeed10 / numYears;

    % Store in output table
    weatherDataWindOnshore = regionalWindOnshore;
    weatherDataWindOnshore.Wind100 = num2cell(averageWindSpeed100, 2);
    weatherDataWindOnshore.Wind10 = num2cell(averageWindSpeed10, 2);

    % Transpose vectors to column format (8760x1)
    for i = 1:height(weatherDataWindOnshore)
        if isrow(weatherDataWindOnshore.Wind100{i})
            weatherDataWindOnshore.Wind100{i} = weatherDataWindOnshore.Wind100{i}';
        end
        if isrow(weatherDataWindOnshore.Wind10{i})
            weatherDataWindOnshore.Wind10{i} = weatherDataWindOnshore.Wind10{i}';
        end
    end

    outputWeatherData = weatherDataWindOnshore;
end

%% FUNCTION CALCTURBINEPRODUCTION
function [outProduction, outFullLoadHours, outMeanWindspeed, outMeanWindspeedOnTurbine] = calcTurbineProduction(inpPastInvestments, inpTurbineParameter, inpTurbinePowerCurves, inpWeatherData)
    %CALCTURBINEPRODUCTION Calculates turbine production, full-load hours and mean wind speed
    %   Computes hourly production time series based on turbine specs, weather input,
    %   and site-specific wind speed adjustment using power curves.
    
    % === Initialization ===
    windparkLocal       = inpPastInvestments;
    turbineData         = inpTurbineParameter;
    powerCurvesLocal    = inpTurbinePowerCurves;
    weatherData         = inpWeatherData;
    
    % Identify all unique regions to model
    uniqueRegions = unique(windparkLocal.nutsID);
    
    % Filter relevant weather data for the selected regions
    idxRelevantWeatherData  = ismember(weatherData.nutsID, uniqueRegions);
    weatherData             = weatherData(idxRelevantWeatherData, :);
    
    % === Define dimensions ===
    hubHeight       = turbineData.hubHeight;
    nrTurbines      = size(turbineData, 1);
    nrRegions       = length(uniqueRegions);
    nrTimestamps    = size(weatherData.Wind10{1}, 1);
    capacity        = repmat(turbineData.power, [1, nrRegions]);
    
    % Define availability and capacity adjustment factors
    techAvail       = 0.97;  % Availability factor
    capfactor       = 1;     % Additional capacity scaling (default: 1)
    
    % === Wind speed extrapolation to hub height ===
    hubHeightMatrix = repmat(hubHeight, nrRegions, nrTimestamps);
    hubHeightMatrix = reshape(hubHeightMatrix', [], 1);
    
    windSpeed10 = cell2mat(weatherData.Wind10);
    windSpeed100 = cell2mat(weatherData.Wind100);
    alpha = log(windSpeed100 ./ windSpeed10) ./ log(100 / 10);
    
    % Reshape capacity matrix
    capacity = reshape(capacity, 1, nrTurbines, nrRegions);
    
    % Calculate hub height wind speed
    windSpeedOnHubheight = windSpeed10 .* ((hubHeightMatrix / 10) .^ alpha);
    windSpeedOnHubheight = reshape(windSpeedOnHubheight, nrTimestamps, nrTurbines, nrRegions);
    
    % Store base 10m wind speed for later averaging
    windSpeed = reshape(windSpeed10, nrTimestamps, nrTurbines, nrRegions);
    
    % === Apply turbulence scaling ===
    xsample = [0 3 4.5 7 11 15 25 35]';
    values = [0 0 0.035 0.035 0.025 0.015 0 0]';
    scalefactor = interp1(xsample, values, windSpeedOnHubheight);
    windSpeedOnHubheightScaled = (1 - scalefactor) .* windSpeedOnHubheight;
    
    % === Power curve interpolation ===
    powerCurve = reshape(powerCurvesLocal.production, [], nrTurbines);
    powerCurveWindspeed = unique(powerCurvesLocal.windspeed);
    
    windProduction = zeros(size(windSpeedOnHubheight));
    for i = 1:nrTurbines
        windProduction(:, i, :) = interp1(powerCurveWindspeed, powerCurve(:, i), windSpeedOnHubheightScaled(:, i, :)) ./ ...
            max(powerCurve(:, i)) .* capacity(:, i, :) .* techAvail .* capfactor;
    end
    
    % === Aggregate results ===
    regioProduction = squeeze(sum(windProduction, 1));
    regioFullLoadHours = regioProduction ./ squeeze(capacity);
    
    % Reshape outputs to column vectors
    outProduction = reshape(regioProduction, [], 1);
    outFullLoadHours = reshape(regioFullLoadHours, [], 1);
    
    % Mean wind speeds (10m and hub height)
    meanWindspeed = reshape(squeeze(mean(windSpeed, 1)), [], 1);
    meanWindspeedOnTurbine = reshape(squeeze(mean(windSpeedOnHubheight, 1)), [], 1);
    
    outMeanWindspeed = meanWindspeed;
    outMeanWindspeedOnTurbine = meanWindspeedOnTurbine;

end

%% FUNCTION CALCTURBINECOMPENSATIONFACTOR
function outputCompensationFactor = calcTurbineCompensationFactor(inpPastInvestments, optsUniqueRefYieldPerRegion)
    %CALCTURBINECOMPENSATIONFACTOR Calculates the compensation factor for wind turbines
    % according to the German Renewable Energy Sources Act (EEG), depending on
    % their location and relative production.
    
    % Copy input to local variable
    windparkLocal = inpPastInvestments;
    
    % Step 1: Calculate the relative production compared to reference
    referenceFactor = windparkLocal.production ./ windparkLocal.referenceProduction;
    
    % Step 2: Set compensation factor to 1 for non-German turbines
    idxNonDE = ~contains(windparkLocal.nutsID,'DE');
    referenceFactor(idxNonDE) = 1;
    
    % Step 3: Define the list of German southern NUTS3 regions eligible for special compensation
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
        'DEB39','DEB3A','DEC02','DEC03','DEC01','DEC04','DEC05','DEC06' ...
    };
    
    % Step 4: Identify which regions are in the southern support zone
    in_sued = ismember(windparkLocal.nutsID, NUTS3_Suedregion);
    
    % Step 5: Define support levels according to EEG (Stützpunkte)
    x = [0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.50];
    y = [1.55 1.42 1.29 1.16 1.07 1.00 0.94 0.89 0.85 0.81 0.79];
    
    % Step 6: Clamp the reference factors based on region
    referenceFactor(in_sued) = max(referenceFactor(in_sued), 0.5);
    referenceFactor(~in_sued) = max(referenceFactor(~in_sued), 0.6);
    referenceFactor = min(referenceFactor, max(x));
    
    % Step 7: Interpolate to find the compensation factor
    localCompensationFactor = interp1(x, y, referenceFactor);

    if optsUniqueRefYieldPerRegion
        isType4 = strcmp(string(windparkLocal.turbineType), "turbineType4");
    
        allRegions = unique(windparkLocal.nutsID);
        for r = 1:numel(allRegions)
            rid  = allRegions{r};
            idxR = strcmp(windparkLocal.nutsID, rid);
    
            idxR_type4 = idxR & isType4;
            if any(idxR_type4)
                % Nimm den (einzigen) Faktor von Typ 4 in dieser Region
                f4 = localCompensationFactor(find(idxR_type4, 1, 'first'));
                localCompensationFactor(idxR) = f4;
            end
        end
    end
    
    % Step 8: Return result
    outputCompensationFactor = localCompensationFactor;

end

%% FUNCTION ESTIMATEDISCRETECHOICEMODEL
function choiceParameter = estimateDiscreteChoiceModel(inpInvestments, inpChoiceModelInput)
% ESTIMATEDISCRETECHOICEMODEL estimates the parameters of a discrete choice
% model (logit) for wind turbine investments in Germany.
%
% INPUTS:
%   - inpInvestments: Table of past regional turbine investments
%   - inpChoiceModelInput: Struct containing settings and variable names
%
% OUTPUT:
%   - choiceParameter: Struct with estimated parameters and metadata

    % Assign input arguments
    investments = inpInvestments;
    choiceModelInput = inpChoiceModelInput;

    % Filter for Germany only (NUTS ID starts with 'DE')
    isGermany = contains(investments.nutsID, 'DE');
    investments = investments(isGermany, :);

    % Set number of decision makers (regions) and explanatory variables
    choiceModelInput.nrDecisionMaker = numel(unique(investments.nutsID));
    choiceModelInput.nrExplanatoryVariables = width(choiceModelInput.explanatoryVariables);

    % Initialize parameters (e.g., zero start values, matrix setup, etc.)
    choiceModelInput = initLogitParameter(choiceModelInput);

    % Logging
    disp('______________________________________')
    disp('Start Parameter Estimation')
    disp('--------------------------------------')

    % Estimate model parameters using log-likelihood maximization
    choiceParameter = estimateLogitParameter(investments, choiceModelInput);

    % Store the used explanatory variables
    choiceParameter.explanatoryVariables = choiceModelInput.explanatoryVariables;
end

