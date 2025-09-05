function [marketData, newInvestments, regionalCapacities] = estimateInvestmentDecision(inpDiscreteChoiceModel, inpMarket, inpTargetInvestments, inpParaTech, optsScaleAreaInvest, optsNpvAdjust, optsRefYieldInInvestment, paraOPEXfix, paraOPEXperkWhvar, paraRegioInvest)
% ESTIMATEINVESTMENTDECISION estimates regional wind power investments based on a discrete choice model.
% 
% INPUTS:
%   inpDiscreteChoiceModel     - Struct containing estimated logit parameters and variable names
%   inpMarket                  - Table with market characteristics
%   inpTargetInvestments       - Struct with investment targets (in MW)
%   inpParaTech                - Struct of technical parameters
%   optsScaleAreaInvest        - Area scaling method ('none', 'BWE', 'WindBG', 'UB21')
%   optsNpvAdjust              - NPV normalization ('NPVperInvestCosts', 'NPVperArea', 'none')
%   optsRefYieldInInvestment   - Boolean, whether to apply compensation correction factor
%   paraOPEXfix                - Fixed OPEX per turbine (EUR/year)
%   paraOPEXperkWhvar          - Variable OPEX per kWh (EUR/kWh)
% 
% OUTPUTS:
%   marketData                 - Updated market table with investment predictions
%   newInvestments             - Total new installed capacity (MW)
%   regionalCapacities         - Estimated total capacity (MW)

% === Load and prepare input data ===
global beta alpha gamma lambda nrAlt nrDecisionMaker

discreteChoiceModel  = inpDiscreteChoiceModel;
dcmParam             = discreteChoiceModel.discreteChoiceParam;
nrAlt                = discreteChoiceModel.nrAlternatives;
explanatoryVariables = discreteChoiceModel.explanatoryVariables;

marketData           = inpMarket;
nrDecisionMaker      = numel(unique(marketData.nutsID));

% === Extract technical parameters ===
compensation     = inpParaTech.paraCompensation;
interestRate     = inpParaTech.paraInterestRate;
powerPotential   = inpParaTech.paraPowerPotential;
powerPotential_invest = inpParaTech.paraPowerPotential_invest;
turbineLifetime  = inpParaTech.paraTurbineLifetime;
investFactor     = inpParaTech.paraInvestFactor;

investCost       = marketData.investCost;
production       = marketData.production;
powerTurbine     = marketData.powerTurbine;

% === Check for area constraints violations ===
% If yes, replace WindBG Limit. Interpret WindBG as a lower limit.

% Group rows by nutsID
[groups, nutsID] = findgroups(marketData.nutsID);

% Aggregate required columns
summed_capacitiesBaseYear = splitapply(@sum, marketData.capacitiesBaseYear, groups);
summed_capacitiesDiscreteChoice = splitapply(@sum, marketData.capacitiesDiscreteChoice, groups);
summed_capacitiesInSimYear = splitapply(@sum, marketData.capacitiesInSimYear, groups);

% Extract representative values for specific columns
totalArea = splitapply(@(x) x(1), marketData.totalArea, groups);
relativeAvailableWindSpace = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace, groups);
relativeAvailableWindSpace_WindBG = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_WindBG, groups);
relativeAvailableWindSpace_BWE = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_BWE, groups);
relativeAvailableWindSpace_UB21 = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_UB21, groups);

% Create tableAdaptWindBG
tableAdaptWindBG = table(nutsID, totalArea, relativeAvailableWindSpace, ...
    relativeAvailableWindSpace_WindBG, relativeAvailableWindSpace_BWE, ...
    relativeAvailableWindSpace_UB21, summed_capacitiesBaseYear, ...
    summed_capacitiesDiscreteChoice, summed_capacitiesInSimYear);

% Compute UsedAreaBase and usedAreaBaseRelative
tableAdaptWindBG.UsedAreaBase = splitapply(@sum, marketData.capacitiesBaseYear./powerPotential, groups);
tableAdaptWindBG.usedAreaBaseRelative = tableAdaptWindBG.UsedAreaBase ./ tableAdaptWindBG.totalArea;

% Add index column where usedAreaBaseRelative > 0.9 * relativeAvailableWindSpace_WindBG
tableAdaptWindBG.index = tableAdaptWindBG.usedAreaBaseRelative > ...
    0.9 * tableAdaptWindBG.relativeAvailableWindSpace_WindBG;

% Update relativeAvailableWindSpace_WindBG where index is 1
tableAdaptWindBG.relativeAvailableWindSpace_WindBG(tableAdaptWindBG.index == 1) = ...
    1.1 * max(tableAdaptWindBG.relativeAvailableWindSpace_WindBG(tableAdaptWindBG.index == 1), tableAdaptWindBG.usedAreaBaseRelative(tableAdaptWindBG.index == 1));
% Expand updated values to match marketData
expanded_WindBG = repelem(tableAdaptWindBG.relativeAvailableWindSpace_WindBG, 8);

% Update the marketData table
marketData.relativeAvailableWindSpace_WindBG = expanded_WindBG;

%%

switch optsScaleAreaInvest
    case 'none'
        marketData.availableWindSpaceSim = marketData.relativeAvailableWindSpace .* marketData.totalArea;
    case 'BWE'
        marketData.availableWindSpaceSim = marketData.relativeAvailableWindSpace_BWE .* marketData.totalArea;
    case 'WindBG'
        marketData.availableWindSpaceSim = marketData.relativeAvailableWindSpace_WindBG .* marketData.totalArea;
    case 'UB21'
        marketData.availableWindSpaceSim = marketData.relativeAvailableWindSpace_UB21 .* marketData.totalArea;
    otherwise
        % Code, der ausgeführt wird, wenn optsScaleAreaParameter keinen der oben genannten Werte hat
        disp('Invalid entry for optsScaleAreaInvest.');
end

availableSpace          = max(marketData.availableWindSpaceSim, 0);
availablePowerPotential = availableSpace .* powerPotential;
compensationFactor      = marketData.compensationFactor;

if ~optsRefYieldInInvestment
    compensationFactor = 1;
end

targetCapacityWind = inpTargetInvestments * 1000;

% === Set logit parameters ===
alpha   = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'alpha')};
beta    = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'beta')};
gamma   = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'gamma')};
lambda  = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'lambda')};

alpha = [alpha(1:3) 0 alpha(4:end)];
alpha = repmat(alpha', [nrDecisionMaker,1]);
beta = repmat(beta, [nrDecisionMaker*nrAlt,1]);
gamma = repmat(gamma, [nrDecisionMaker*nrAlt,1]);
lambda = repmat(lambda, [nrDecisionMaker*nrAlt,1]);

% Add reduction of CAPEX according to https://www.irena.org/Publications/2024/Sep/Renewable-Power-Generation-Costs-in-2023
% In Germany: yearly decrease of 2.4 %
%% === Apply CAPEX reduction over time ===
marketData.investCost = applyCostReduction( ...
    marketData.investCost, ...
    paraRegioInvest.costReductionRate, ...
    paraRegioInvest.baseYear, ...
    paraRegioInvest.simYear);

explVar(:,1) = calcNetPresentValue(marketData.investCost,production,turbineLifetime,interestRate,compensation,compensationFactor, paraOPEXfix, paraOPEXperkWhvar);
meanNPV = mean(explVar(:,1));
% explVar(:,1) = (explVar(:,1) - rescaleParaNPV(1)) ./(rescaleParaNPV(2)- rescaleParaNPV(1))+1;
disp(['Average Net Present Value considering a Compensation of ',num2str(compensation.*100),'ct/kWh : ',num2str(floor(meanNPV)) , '€'])
if numel(explanatoryVariables)>1
    explVar(:,2:numel(explanatoryVariables)) = marketData{:,explanatoryVariables(2:end)};
end

% === Investment decision loop ===
countWhile = 0;
installedCapacity = sum(marketData.capacitiesInSimYear);

while (installedCapacity <= targetCapacityWind || installedCapacity / targetCapacityWind > 1.005) && countWhile < 100000

    npv = calcNetPresentValue(marketData.investCost, production, turbineLifetime, interestRate, compensation, compensationFactor, paraOPEXfix, paraOPEXperkWhvar);

    % Optional NPV adjustment
    switch optsNpvAdjust
        case 'NPVperInvestCosts'
            npv2 = npv ./ marketData.investCost;
        case 'NPVperArea'
            npv2 = npv ./ (marketData.powerTurbine ./ powerPotential);
        case 'none'
        otherwise
            error('Invalid optsNpvAdjust value.');
    end
    npv_transformed = transformNPV(npv2);

    explVar = zeros(height(marketData), numel(explanatoryVariables));
    explVar(:,1) = npv_transformed;
    if numel(explanatoryVariables) > 1
        explVar(:,2:end) = marketData{:,explanatoryVariables(2:end)};
    end

    % Calculate market shares
    predictedMarketShare = binomNestLog(explVar);

    % Update capacity and compensation
    estimatedNewCap = floor(investFactor .* availableSpace .* powerPotential_invest ./ powerTurbine .* predictedMarketShare) .* powerTurbine;
    installedCapacity = sum(estimatedNewCap + marketData.capacitiesInSimYear);

    if installedCapacity / targetCapacityWind > 1.005
        compensation = compensation - 0.00005;
    else
        compensation = compensation + 0.001;
    end
    countWhile = countWhile + 1;
end

% === Final outputs ===
marketData.estimatedNewCapacity = estimatedNewCap;
marketData.estimatedTotalCapacity = estimatedNewCap + marketData.capacitiesInSimYear;

% Weighted compensation
marketData.share_of_total_generation = (marketData.estimatedTotalCapacity .* marketData.production) / sum(marketData.estimatedTotalCapacity .* marketData.production);
marketData.weighted_comp_share = marketData.share_of_total_generation .* marketData.compensationFactor;
real_compensation = sum(marketData.weighted_comp_share) * compensation;

disp(['Estimated capacity: ', num2str(installedCapacity/1e6), ' GW'])
disp(['Target capacity: ', num2str(targetCapacityWind/1e6), ' GW'])
disp(['Final compensation: ', num2str(compensation*100), ' ct/kWh'])
disp(['Real avg. compensation: ', num2str(real_compensation*100), ' ct/kWh'])
disp(['Total NPV of the power plants is: ' num2str(sum((marketData.estimatedTotalCapacity./marketData.powerTurbine).*npv)/10^9) ' bn. EUR'])

newInvestments = sum(marketData.estimatedNewCapacity);
regionalCapacities = sum(marketData.estimatedTotalCapacity);

end

%% FUNCTION BINOMNESTLOG
function y = binomNestLog(x)

    global beta alpha gamma lambda nrAlt nrDecisionMaker
    vni = exp((sum(x.*beta,2)+alpha)./lambda);
    
    vniSum = repelem(sum(reshape(vni,[nrAlt,nrDecisionMaker]),1),nrAlt)';
    
    InkVal = log(vniSum);
    
    propNest = exp(gamma+lambda.*InkVal)./(1+exp(gamma+lambda.*InkVal));
    propAlternative = vni./vniSum;
    
    y = propNest.*propAlternative;
end

%% FUNCTION UPDATE (REDUCE) INVEST COSTS FOR FUTURE YEARS
function updatedInvestCost = applyCostReduction(investCost, annualReductionRate, baseYear, targetYear)
    yearsToTarget = targetYear - baseYear;
    reductionFactor = (1 - annualReductionRate) .^ yearsToTarget;
    updatedInvestCost = investCost .* reductionFactor;
end

