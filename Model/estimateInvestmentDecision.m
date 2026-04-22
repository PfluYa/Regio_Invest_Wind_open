function [marketData, newInvestments, regionalCapacities, resultsScaling] = ...
    estimateInvestmentDecision( ...
        inpDiscreteChoiceModel, ...
        inpMarket, ...
        inpTargetInvestments, ...
        inpParaTech, ...
        optsScaleAreaInvest, ...
        optsNpvAdjust, ...
        optsRefYieldInInvestment, ...
        paraOPEXfix, ...
        paraOPEXperkWhvar, ...
        paraRegioInvest, ...
        optsRegioInvest)
% ESTIMATEINVESTMENTDECISION
% Estimate regional wind investment outcomes using the nested-logit model.
%
% DESCRIPTION
% This function:
%   1. prepares the market data and land constraints,
%   2. performs a sensitivity analysis over reference-yield scaling factors,
%   3. stores the scaling results to disk,
%   4. runs the final baseline simulation for scaleRefYield = 1.
%
% OUTPUTS
%   marketData          : full result table including estimated capacities
%   newInvestments      : total newly installed capacity
%   regionalCapacities  : total installed capacity after investment
%   resultsScaling      : sensitivity results over ref-yield scaling factors
%
% IMPORTANT
% resultsScaling is saved automatically to:
%   Results Paper/Results_scaling_capacity_based.mat
%   or
%   Results Paper/Results_scaling_energy_based.mat
% depending on optsRegioInvest.iterationLogicInvestment.
%
% CURRENT DEFAULTS (from main script)
%   optsRegioInvest.iterationLogicInvestment = 'capacity'
%   optsNpvAdjust                            = 'NPVperArea'
%   optsRefYieldInInvestment                 = true

%% ------------------------------------------------------------------------
% 0. GLOBALS USED BY THE NESTED-LOGIT PREDICTION
% -------------------------------------------------------------------------

global beta alpha gamma lambda nrAlt nrDecisionMaker

%% ------------------------------------------------------------------------
% 1. LOAD INPUTS
% -------------------------------------------------------------------------

discreteChoiceModel  = inpDiscreteChoiceModel;
dcmParam             = discreteChoiceModel.discreteChoiceParam;
nrAlt                = discreteChoiceModel.nrAlternatives;
explanatoryVariables = discreteChoiceModel.explanatoryVariables;

marketData       = inpMarket;
nrDecisionMaker  = numel(unique(marketData.nutsID));

% Technical parameters
compensation          = inpParaTech.paraCompensation;
interestRate          = inpParaTech.paraInterestRate;
powerPotential        = inpParaTech.paraPowerPotential;
powerPotentialInvest  = inpParaTech.paraPowerPotentialInvest;
turbineLifetime       = inpParaTech.paraTurbineLifetime;
investFactor          = inpParaTech.paraInvestFactor;

production    = marketData.production;
powerTurbine  = marketData.powerTurbine;

origCompensation = compensation;

%% ------------------------------------------------------------------------
% 2. CHECK / ADAPT AREA CONSTRAINTS
% -------------------------------------------------------------------------

[groups, nutsID] = findgroups(marketData.nutsID);

summedCapacitiesBaseYear = splitapply(@sum, marketData.capacitiesBaseYear, groups);
summedCapacitiesDiscreteChoice = splitapply(@sum, marketData.capacitiesDiscreteChoice, groups);
summedCapacitiesInSimYear = splitapply(@sum, marketData.capacitiesInSimYear, groups);

totalArea = splitapply(@(x) x(1), marketData.totalArea, groups);
relativeAvailableWindSpace = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace, groups);
relativeAvailableWindSpaceWindBG = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_WindBG, groups);
relativeAvailableWindSpaceBWE = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_BWE, groups);
relativeAvailableWindSpaceUB21 = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_UB21, groups);

tableAdaptWindBG = table( ...
    nutsID, ...
    totalArea, ...
    relativeAvailableWindSpace, ...
    relativeAvailableWindSpaceWindBG, ...
    relativeAvailableWindSpaceBWE, ...
    relativeAvailableWindSpaceUB21, ...
    summedCapacitiesBaseYear, ...
    summedCapacitiesDiscreteChoice, ...
    summedCapacitiesInSimYear);

tableAdaptWindBG.usedAreaBase = ...
    splitapply(@sum, marketData.capacitiesBaseYear ./ powerPotential, groups);

tableAdaptWindBG.usedAreaBaseRelative = ...
    tableAdaptWindBG.usedAreaBase ./ tableAdaptWindBG.totalArea;

tableAdaptWindBG.index = ...
    tableAdaptWindBG.usedAreaBaseRelative > ...
    0.9 .* tableAdaptWindBG.relativeAvailableWindSpaceWindBG;

tableAdaptWindBG.relativeAvailableWindSpaceWindBG(tableAdaptWindBG.index == 1) = ...
    1.1 .* max( ...
        tableAdaptWindBG.relativeAvailableWindSpaceWindBG(tableAdaptWindBG.index == 1), ...
        tableAdaptWindBG.usedAreaBaseRelative(tableAdaptWindBG.index == 1));

% NOTE:
% Current model assumes 8 turbine clusters when expanding from regional to
% region-technology level.
expandedWindBG = repelem(tableAdaptWindBG.relativeAvailableWindSpaceWindBG, 8);
marketData.relativeAvailableWindSpace_WindBG = expandedWindBG;

%% ------------------------------------------------------------------------
% 3. SET AVAILABLE LAND UNDER CHOSEN SCALING RULE
% -------------------------------------------------------------------------

switch optsScaleAreaInvest
    case 'none'
        marketData.availableWindSpaceSim = ...
            marketData.relativeAvailableWindSpace .* marketData.totalArea;

    case 'BWE'
        marketData.availableWindSpaceSim = ...
            marketData.relativeAvailableWindSpace_BWE .* marketData.totalArea;

    case 'WindBG'
        marketData.availableWindSpaceSim = ...
            marketData.relativeAvailableWindSpace_WindBG .* marketData.totalArea;

    case 'UB21'
        marketData.availableWindSpaceSim = ...
            marketData.relativeAvailableWindSpace_UB21 .* marketData.totalArea;

    otherwise
        error('Invalid value for optsScaleAreaInvest.');
end

availableSpace = max(marketData.availableWindSpaceSim, 0);
compensationFactor = marketData.compensationFactor;

%% ------------------------------------------------------------------------
% 4. DEFINE REGIONAL CLUSTERS: NORTH / SOUTH / NEUTRAL
% -------------------------------------------------------------------------

nuts3 = string(marketData.nutsID);

southList = [ ...
"DE145","DE121","DE146","DE112","DE147","DE132","DE12A","DE133","DE12B","DE113", ...
"DE131","DE12C","DE114","DE125","DE11C","DE117","DE118","DE119","DE122","DE123", ...
"DE138","DE139","DE115","DE11B","DE126","DE127","DE11D","DE129","DE124","DE148", ...
"DE116","DE141","DE128","DE135","DE11A","DE136","DE149","DE111","DE142","DE137", ...
"DE144","DE13A","DE143","DE275","DE214","DE231","DE234","DE251","DE256","DE261", ...
"DE264","DE271","DE276","DE216","DE241","DE245","DE242","DE246","DE215","DE235", ...
"DE217","DE224","DE277","DE22C","DE27D","DE218","DE219","DE21A","DE252","DE257", ...
"DE248","DE21B","DE225","DE21C","DE253","DE258","DE21D","DE278","DE267","DE211", ...
"DE272","DE226","DE273","DE268","DE21E","DE221","DE227","DE27A","DE26A","DE274", ...
"DE21F","DE269","DE21G","DE212","DE21H","DE21I","DE236","DE25A","DE237","DE279", ...
"DE254","DE259","DE27E","DE27B","DE222","DE228","DE21J","DE229","DE232","DE238", ...
"DE213","DE21K","DE25B","DE22A","DE255","DE239","DE262","DE26B","DE21L","DE223", ...
"DE22B","DE23A","DE21M","DE27C","DE233","DE21N","DE25C","DE263","DE26C","DE715", ...
"DE711","DE716","DE717","DE71B","DE71C","DEB3B","DEB3C","DEB14","DEB22","DEB15", ...
"DEB3D","DEB23","DEB31","DEB3E","DEB32","DEB3F","DEB3G","DEB33","DEB34","DEB35", ...
"DEB3J","DEB36","DEB37","DEB1D","DEB3I","DEB38","DEB3H","DEB3K","DEB21","DEB25", ...
"DEB39","DEB3A","DEC02","DEC03","DEC01","DEC04","DEC05","DEC06" ];

northList = [ ...
"DE80M","DE804","DE80K","DE803","DE80L","DE80N", ...
"DE942","DE94C","DE94H","DE945","DE94A","DE946","DE94G","DE943","DE94D","DE941", ...
"DE501","DE936","DE502","DE932","DE939","DE600", ...
"DEF09","DEF0D","DEF03","DEF08","DEF0E","DEF05","DEF07","DEF0C", ...
"DEF01","DEF0B","DEF02","DEF0A","DEF04" ];

cluster = repmat("Neutral", height(marketData), 1);
cluster(ismember(nuts3, southList)) = "South";
cluster(ismember(nuts3, northList)) = "North";

overlap = ismember(nuts3, southList) & ismember(nuts3, northList);
if any(overlap)
    warning('Some NUTS3 regions are in both South and North. Assigned to South.');
    cluster(overlap) = "South";
end

marketData.Cluster = categorical(cluster);

origMarketData = marketData;
origCompensationFactor = compensationFactor;

%% ------------------------------------------------------------------------
% 5. DEFINE TARGETS
% -------------------------------------------------------------------------

targetCapacityWind = inpTargetInvestments * 1000;  % convert MW to kW

% Current hard-coded energy target used in the original code.
% Corresponds to the 2040 reference case.
targetEnergyProduction = 353.6;

% Current tolerance used for both capacity and energy logic.
tolEnergy = 0.002;

%% ------------------------------------------------------------------------
% 6. PREPARE LOGIT PARAMETERS
% -------------------------------------------------------------------------

alpha = dcmParam{:, startsWith(dcmParam.Properties.VariableNames, 'alpha')};
beta  = dcmParam{:, startsWith(dcmParam.Properties.VariableNames, 'beta')};
gamma = dcmParam{:, startsWith(dcmParam.Properties.VariableNames, 'gamma')};
lambda = dcmParam{:, startsWith(dcmParam.Properties.VariableNames, 'lambda')};

% Reinsert base-category alpha = 0 for turbine type 4
alpha = [alpha(1:3) 0 alpha(4:end)];

alpha  = repmat(alpha', [nrDecisionMaker, 1]);
beta   = repmat(beta,   [nrDecisionMaker * nrAlt, 1]);
gamma  = repmat(gamma,  [nrDecisionMaker * nrAlt, 1]);
lambda = repmat(lambda, [nrDecisionMaker * nrAlt, 1]);

%% ------------------------------------------------------------------------
% 7. SENSITIVITY ANALYSIS OVER REFERENCE-YIELD SCALING
% -------------------------------------------------------------------------

refValue = 1.00;
scaleVec = 0:0.005:1.7;   % current grid of scaling factors
nScale   = numel(scaleVec);

resultsScaling = table( ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    'VariableNames', { ...
        'scaleRefYield', ...
        'installedCapacityGW', ...
        'targetCapacityGW', ...
        'finalComp_ctperkWh', ...
        'realComp_ctperkWh', ...
        'producerSurplus_bnEUR', ...
        'subsidyPV_bnEUR', ...
        'welfareNoGrid_bnEUR', ...
        'totalGeneration_annual', ...
        'capNorth_GW', ...
        'capSouth_GW', ...
        'capNeutral_GW', ...
        'rent_per_kWh', ...
        'avoidedVarCostGasCO2_bnEUR', ...
        'welfarePlusGasCO2_bnEUR'});

for iScale = 1:nScale

    % Reset to baseline state
    marketData = origMarketData;
    compensationFactor = origCompensationFactor;
    compensation = origCompensation;

    scaleRefYield = scaleVec(iScale);

    % Scale compensation factor around 1
    compensationFactor = refValue + scaleRefYield .* (compensationFactor - refValue);

    if ~optsRefYieldInInvestment
        compensationFactor = ones(size(compensationFactor));
    end

    % Apply future CAPEX reduction
    marketData.investCost = applyCostReduction( ...
        marketData.investCost, ...
        paraRegioInvest.costReductionRate, ...
        paraRegioInvest.baseYear, ...
        paraRegioInvest.simYear);

    [estimatedNewCap, compensation, installedCapacity, totalGenerationAnnual, ...
        ~, npv, npvCostsTotalLifetime, npvRevenuesTotalLifetime] = ...
        runInvestmentLoop( ...
            marketData, ...
            production, ...
            turbineLifetime, ...
            interestRate, ...
            compensation, ...
            compensationFactor, ...
            paraOPEXfix, ...
            paraOPEXperkWhvar, ...
            optsNpvAdjust, ...
            explanatoryVariables, ...
            powerPotential, ...
            powerPotentialInvest, ...
            powerTurbine, ...
            investFactor, ...
            availableSpace, ...
            optsRegioInvest.iterationLogicInvestment, ...
            targetCapacityWind, ...
            targetEnergyProduction, ...
            tolEnergy);

    % Store capacities
    marketData.estimatedNewCapacity = estimatedNewCap;
    marketData.estimatedTotalCapacity = estimatedNewCap + marketData.capacitiesInSimYear;

    % Production-weighted effective remuneration
    marketData.share_of_total_generation = ...
        ((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* marketData.production) ./ ...
        sum((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* marketData.production);

    marketData.weighted_comp_share = ...
        marketData.share_of_total_generation .* compensationFactor;

    realCompensation = sum(marketData.weighted_comp_share) .* compensation;

    % Annualized producer surplus
    annualizationFactor = ...
        (interestRate .* (1 + interestRate)^turbineLifetime) ./ ...
        ((1 + interestRate)^turbineLifetime - 1);

    producerSurplusOneYear = sum( ...
        (marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* ...
        (npvRevenuesTotalLifetime .* annualizationFactor - ...
         npvCostsTotalLifetime .* annualizationFactor));

    % Annual generation
    nTurbines = marketData.estimatedTotalCapacity ./ marketData.powerTurbine;
    annualGeneration = nTurbines .* marketData.production;
    totalGenerationAnnual = sum(annualGeneration);

    % Subsidy payments
    subsidyOneYear = totalGenerationAnnual .* realCompensation;

    % Welfare excluding grid/system costs
    totalWelfare = producerSurplusOneYear - subsidyOneYear;

    % Avoided variable gas + CO2 generation costs
    avoidVarCost_EURperMWh = 76.38;  % current appendix value
    gen_MWh_el = totalGenerationAnnual ./ 1000;
    avoidedVarCostOneYear = gen_MWh_el .* avoidVarCost_EURperMWh;

    totalWelfarePlusGasCO2 = totalWelfare + avoidedVarCostOneYear;

    rentPerkWh = producerSurplusOneYear ./ totalGenerationAnnual;

    % Aggregate capacity by cluster
    sumCap = groupsummary(marketData, "Cluster", "sum", "estimatedTotalCapacity");

    capNorth = 0;
    capSouth = 0;
    capNeutral = 0;

    idxNorth = sumCap.Cluster == 'North';
    if any(idxNorth)
        capNorth = sumCap.sum_estimatedTotalCapacity(idxNorth);
    end

    idxSouth = sumCap.Cluster == 'South';
    if any(idxSouth)
        capSouth = sumCap.sum_estimatedTotalCapacity(idxSouth);
    end

    idxNeutral = sumCap.Cluster == 'Neutral';
    if any(idxNeutral)
        capNeutral = sumCap.sum_estimatedTotalCapacity(idxNeutral);
    end

    capNorth_GW   = capNorth   / 1e6;
    capSouth_GW   = capSouth   / 1e6;
    capNeutral_GW = capNeutral / 1e6;

    % Write output row
    resultsScaling.scaleRefYield(iScale)               = scaleRefYield;
    resultsScaling.installedCapacityGW(iScale)         = installedCapacity / 1e6;
    resultsScaling.targetCapacityGW(iScale)            = targetCapacityWind / 1e6;
    resultsScaling.finalComp_ctperkWh(iScale)          = compensation * 100;
    resultsScaling.realComp_ctperkWh(iScale)           = realCompensation * 100;
    resultsScaling.producerSurplus_bnEUR(iScale)       = producerSurplusOneYear / 1e9;
    resultsScaling.subsidyPV_bnEUR(iScale)             = subsidyOneYear / 1e9;
    resultsScaling.welfareNoGrid_bnEUR(iScale)         = totalWelfare / 1e9;
    resultsScaling.avoidedVarCostGasCO2_bnEUR(iScale)  = avoidedVarCostOneYear / 1e9;
    resultsScaling.welfarePlusGasCO2_bnEUR(iScale)     = totalWelfarePlusGasCO2 / 1e9;
    resultsScaling.totalGeneration_annual(iScale)      = totalGenerationAnnual / 1e9;
    resultsScaling.capNorth_GW(iScale)                 = capNorth_GW;
    resultsScaling.capSouth_GW(iScale)                 = capSouth_GW;
    resultsScaling.capNeutral_GW(iScale)               = capNeutral_GW;
    resultsScaling.rent_per_kWh(iScale)                = rentPerkWh;
end

%% ------------------------------------------------------------------------
% 8. SAVE SCALING RESULTS
% -------------------------------------------------------------------------

resultsFolder = fullfile(pwd, 'Results Paper');
if ~isfolder(resultsFolder)
    mkdir(resultsFolder);
end

switch optsRegioInvest.iterationLogicInvestment
    case 'energy'
        scalingFilename = 'Results_scaling_energy_based.mat';
    case 'capacity'
        scalingFilename = 'Results_scaling_capacity_based.mat';
    otherwise
        error('Invalid optsRegioInvest.iterationLogicInvestment. Use ''capacity'' or ''energy''.');
end

save(fullfile(resultsFolder, scalingFilename), 'resultsScaling');

%% ------------------------------------------------------------------------
% 9. PLOT SCALING RESULTS
% -------------------------------------------------------------------------

x = resultsScaling.scaleRefYield;

figure('Position', [100 100 1000 800]);

subplot(3,2,1);
plot(x, resultsScaling.welfarePlusGasCO2_bnEUR, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('Welfare (bn EUR / year)');
title('Welfare');

subplot(3,2,2);
plot(x, resultsScaling.finalComp_ctperkWh, 'LineWidth', 2); hold on;
plot(x, resultsScaling.realComp_ctperkWh, '--', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('ct/kWh');
title('Remuneration');
legend('Required (auction price)', 'Real (production-weighted)', 'Location', 'best');

subplot(3,2,3);
plot(x, resultsScaling.producerSurplus_bnEUR, 'LineWidth', 2); hold on;
plot(x, resultsScaling.subsidyPV_bnEUR, '--', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('bn EUR / year');
title('Rents and fiscal costs');
legend('Producer surplus', 'Subsidy', 'Location', 'best');

subplot(3,2,4);
plot(x, resultsScaling.totalGeneration_annual, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('TWh / year');
title('Total annual generation');

subplot(3,2,5);
plot(x, resultsScaling.capNorth_GW, 'LineWidth', 2); hold on;
plot(x, resultsScaling.capNeutral_GW, '--', 'LineWidth', 2);
plot(x, resultsScaling.capSouth_GW, ':', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('GW');
title('Regional capacities');
legend('North', 'Neutral', 'South', 'Location', 'best');

subplot(3,2,6);
plot(x, resultsScaling.rent_per_kWh, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('EUR / kWh');
title('Annualized rent per kWh');

sgtitle('Sensitivity to reference-yield scaling', 'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% 10. FINAL BASELINE RUN WITH scaleRefYield = 1
% -------------------------------------------------------------------------

marketData = origMarketData;
compensationFactor = origCompensationFactor;
compensation = origCompensation;

scaleRefYieldFinal = 1.0;
refValue = 1.0;

compensationFactor = refValue + scaleRefYieldFinal .* (compensationFactor - refValue);

if ~optsRefYieldInInvestment
    compensationFactor = ones(size(compensationFactor));
end

marketData.investCost = applyCostReduction( ...
    marketData.investCost, ...
    paraRegioInvest.costReductionRate, ...
    paraRegioInvest.baseYear, ...
    paraRegioInvest.simYear);

initialNPV = calcNetPresentValue( ...
    marketData.investCost, ...
    production, ...
    turbineLifetime, ...
    interestRate, ...
    compensation, ...
    compensationFactor, ...
    paraOPEXfix, ...
    paraOPEXperkWhvar);

meanNPV = mean(initialNPV);
disp(['Average net present value at ', num2str(compensation * 100), ...
      ' ct/kWh: ', num2str(floor(meanNPV)), ' EUR'])

[estimatedNewCap, compensation, installedCapacity, totalGenerationAnnual, ...
    ~, npv] = runInvestmentLoop( ...
        marketData, ...
        production, ...
        turbineLifetime, ...
        interestRate, ...
        compensation, ...
        compensationFactor, ...
        paraOPEXfix, ...
        paraOPEXperkWhvar, ...
        optsNpvAdjust, ...
        explanatoryVariables, ...
        powerPotential, ...
        powerPotentialInvest, ...
        powerTurbine, ...
        investFactor, ...
        availableSpace, ...
        optsRegioInvest.iterationLogicInvestment, ...
        targetCapacityWind, ...
        targetEnergyProduction, ...
        tolEnergy);

marketData.estimatedNewCapacity = estimatedNewCap;
marketData.estimatedTotalCapacity = estimatedNewCap + marketData.capacitiesInSimYear;

marketData.share_of_total_generation = ...
    ((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* marketData.production) ./ ...
    sum((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* marketData.production);

marketData.weighted_comp_share = ...
    marketData.share_of_total_generation .* compensationFactor;

realCompensation = sum(marketData.weighted_comp_share) .* compensation;

disp(['Estimated capacity: ', num2str(installedCapacity / 1e6), ' GW'])
disp(['Target capacity: ', num2str(targetCapacityWind / 1e6), ' GW'])
disp(['Final compensation: ', num2str(compensation * 100), ' ct/kWh'])
disp(['Real avg. compensation: ', num2str(realCompensation * 100), ' ct/kWh'])
disp(['Total NPV of the power plants: ', ...
      num2str(sum((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* npv) / 1e9), ...
      ' bn EUR'])

newInvestments = sum(marketData.estimatedNewCapacity);
regionalCapacities = sum(marketData.estimatedTotalCapacity);

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function y = binomNestLog(x)
% BINOMNESTLOG
% Compute predicted market shares from the nested-logit model.

    global beta alpha gamma lambda nrAlt nrDecisionMaker

    vni = exp((sum(x .* beta, 2) + alpha) ./ lambda);
    vniSum = repelem(sum(reshape(vni, [nrAlt, nrDecisionMaker]), 1), nrAlt)';
    inkVal = log(vniSum);

    propNest = exp(gamma + lambda .* inkVal) ./ ...
               (1 + exp(gamma + lambda .* inkVal));

    propAlternative = vni ./ vniSum;
    y = propNest .* propAlternative;
end

function updatedInvestCost = applyCostReduction(investCost, annualReductionRate, baseYear, targetYear)
% APPLYCOSTREDUCTION
% Apply compound annual CAPEX reduction between base year and target year.
%
% Current default annualReductionRate in main script: 0.0238

    yearsToTarget = targetYear - baseYear;
    reductionFactor = (1 - annualReductionRate) .^ yearsToTarget;
    updatedInvestCost = investCost .* reductionFactor;
end

function [estimatedNewCap, compensation, installedCapacity, totalGenerationAnnual, ...
          countWhile, npv, npvCostsTotalLifetime, npvRevenuesTotalLifetime] = ...
    runInvestmentLoop( ...
        marketData, production, turbineLifetime, interestRate, compensation, compensationFactor, ...
        paraOPEXfix, paraOPEXperkWhvar, optsNpvAdjust, explanatoryVariables, ...
        powerPotential, powerPotentialInvest, powerTurbine, ...
        investFactor, availableSpace, ...
        iterationLogicInvestment, ...
        targetCapacityWind, targetEnergyProduction, tolEnergy)
% RUNINVESTMENTLOOP
% Iteratively adjust remuneration until either the capacity target or the
% energy target is met, depending on iterationLogicInvestment.
%
% INPUT
% iterationLogicInvestment:
%   'capacity' -> match installed capacity target
%   'energy'   -> match annual generation target
%
% CURRENT DEFAULT
%   'capacity'

    maxIter = 100000;
    countWhile = 0;

    estimatedNewCap = zeros(height(marketData), 1);

    function gen = computeAnnualGeneration(totalCap)
        nTurb = totalCap ./ powerTurbine;
        gen = sum(nTurb .* marketData.production);
    end

    installedCapacity = sum(marketData.capacitiesInSimYear);
    totalGenerationAnnual = computeAnnualGeneration(marketData.capacitiesInSimYear);

    if strcmp(iterationLogicInvestment, 'capacity')

        while ((installedCapacity <= targetCapacityWind || ...
                installedCapacity / targetCapacityWind > (1 + tolEnergy)) && ...
                countWhile < maxIter)

            [npv, npvCostsTotalLifetime, npvRevenuesTotalLifetime] = calcNetPresentValue( ...
                marketData.investCost, ...
                production, ...
                turbineLifetime, ...
                interestRate, ...
                compensation, ...
                compensationFactor, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar);

            switch optsNpvAdjust
                case 'NPVperInvestCosts'
                    npvAdjusted = npv ./ marketData.investCost;
                case 'NPVperArea'
                    npvAdjusted = npv ./ (marketData.powerTurbine ./ powerPotential);
                case 'none'
                    npvAdjusted = npv;
                otherwise
                    error('Invalid optsNpvAdjust value.');
            end

            npvTransformed = transformNPV(npvAdjusted);

            explVar = zeros(height(marketData), numel(explanatoryVariables));
            explVar(:,1) = npvTransformed;

            if numel(explanatoryVariables) > 1
                explVar(:,2:end) = marketData{:, explanatoryVariables(2:end)};
            end

            predictedMarketShare = binomNestLog(explVar);

            estimatedNewCap = ...
                floor(investFactor .* availableSpace .* powerPotentialInvest ./ powerTurbine .* ...
                predictedMarketShare) .* powerTurbine;

            installedCapacity = sum(estimatedNewCap + marketData.capacitiesInSimYear);

            if installedCapacity / targetCapacityWind > (1 + tolEnergy)
                compensation = compensation - 0.00005;
            else
                compensation = compensation + 0.001;
            end

            countWhile = countWhile + 1;
        end

        totalGenerationAnnual = ...
            computeAnnualGeneration(estimatedNewCap + marketData.capacitiesInSimYear);

    elseif strcmp(iterationLogicInvestment, 'energy')

        totalGenerationAnnual = computeAnnualGeneration(marketData.capacitiesInSimYear);

        while ((totalGenerationAnnual <= targetEnergyProduction || ...
                totalGenerationAnnual / targetEnergyProduction > (1 + tolEnergy)) && ...
                countWhile < maxIter)

            [npv, npvCostsTotalLifetime, npvRevenuesTotalLifetime] = calcNetPresentValue( ...
                marketData.investCost, ...
                production, ...
                turbineLifetime, ...
                interestRate, ...
                compensation, ...
                compensationFactor, ...
                paraOPEXfix, ...
                paraOPEXperkWhvar);

            switch optsNpvAdjust
                case 'NPVperInvestCosts'
                    npvAdjusted = npv ./ marketData.investCost;
                case 'NPVperArea'
                    npvAdjusted = npv ./ (marketData.powerTurbine ./ powerPotential);
                case 'none'
                    npvAdjusted = npv;
                otherwise
                    error('Invalid optsNpvAdjust value.');
            end

            npvTransformed = transformNPV(npvAdjusted);

            explVar = zeros(height(marketData), numel(explanatoryVariables));
            explVar(:,1) = npvTransformed;

            if numel(explanatoryVariables) > 1
                explVar(:,2:end) = marketData{:, explanatoryVariables(2:end)};
            end

            predictedMarketShare = binomNestLog(explVar);

            estimatedNewCap = ...
                floor(investFactor .* availableSpace .* powerPotentialInvest ./ powerTurbine .* ...
                predictedMarketShare) .* powerTurbine;

            totalCap = estimatedNewCap + marketData.capacitiesInSimYear;
            totalGenerationAnnual = computeAnnualGeneration(totalCap);

            if totalGenerationAnnual / targetEnergyProduction > (1 + tolEnergy)
                compensation = compensation - 0.00005;
            else
                compensation = compensation + 0.001;
            end

            countWhile = countWhile + 1;
        end

        installedCapacity = sum(estimatedNewCap + marketData.capacitiesInSimYear);

    else
        error('Unknown iterationLogicInvestment. Use ''capacity'' or ''energy''.');
    end
end