function [marketData, newInvestments, regionalCapacities, resultsScaling] = estimateInvestmentDecision(inpDiscreteChoiceModel, inpMarket, inpTargetInvestments, inpParaTech, optsScaleAreaInvest, optsNpvAdjust, optsRefYieldInInvestment, paraOPEXfix, paraOPEXperkWhvar, paraRegioInvest, optsRegioInvest)
% ESTIMATEINVESTMENTDECISION estimates regional wind power investments based on a discrete choice model.

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

investCost       = marketData.investCost; %#ok<NASGU>
production       = marketData.production;
powerTurbine     = marketData.powerTurbine;

origCompensation = compensation;

% === Check for area constraints violations ===
[groups, nutsID] = findgroups(marketData.nutsID);

summed_capacitiesBaseYear       = splitapply(@sum, marketData.capacitiesBaseYear, groups);
summed_capacitiesDiscreteChoice = splitapply(@sum, marketData.capacitiesDiscreteChoice, groups);
summed_capacitiesInSimYear      = splitapply(@sum, marketData.capacitiesInSimYear, groups);

totalArea                        = splitapply(@(x) x(1), marketData.totalArea, groups);
relativeAvailableWindSpace       = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace, groups);
relativeAvailableWindSpace_WindBG= splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_WindBG, groups);
relativeAvailableWindSpace_BWE   = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_BWE, groups);
relativeAvailableWindSpace_UB21  = splitapply(@(x) x(1), marketData.relativeAvailableWindSpace_UB21, groups);

tableAdaptWindBG = table(nutsID, totalArea, relativeAvailableWindSpace, ...
    relativeAvailableWindSpace_WindBG, relativeAvailableWindSpace_BWE, ...
    relativeAvailableWindSpace_UB21, summed_capacitiesBaseYear, ...
    summed_capacitiesDiscreteChoice, summed_capacitiesInSimYear);

tableAdaptWindBG.UsedAreaBase = splitapply(@sum, marketData.capacitiesBaseYear./powerPotential, groups);
tableAdaptWindBG.usedAreaBaseRelative = tableAdaptWindBG.UsedAreaBase ./ tableAdaptWindBG.totalArea;

tableAdaptWindBG.index = tableAdaptWindBG.usedAreaBaseRelative > ...
    0.9 * tableAdaptWindBG.relativeAvailableWindSpace_WindBG;

tableAdaptWindBG.relativeAvailableWindSpace_WindBG(tableAdaptWindBG.index == 1) = ...
    1.1 * max(tableAdaptWindBG.relativeAvailableWindSpace_WindBG(tableAdaptWindBG.index == 1), tableAdaptWindBG.usedAreaBaseRelative(tableAdaptWindBG.index == 1));

expanded_WindBG = repelem(tableAdaptWindBG.relativeAvailableWindSpace_WindBG, 8);
marketData.relativeAvailableWindSpace_WindBG = expanded_WindBG;

%% === Available area depending on chosen scaling ===
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
        disp('Invalid entry for optsScaleAreaInvest.');
end

availableSpace          = max(marketData.availableWindSpaceSim, 0);
availablePowerPotential = availableSpace .* powerPotential;
compensationFactor      = marketData.compensationFactor;

%% === Cluster: North / South / Neutral auf Basis NUTS1 ===
nuts1 = cellfun(@(x) x(1:3), marketData.nutsID, 'UniformOutput', false);
nuts1 = string(nuts1);

southStates   = ["DE1", "DE2"];
neutralStates = ["DEB", "DEC", "DE7", "DEG"];

cluster = repmat("North", height(marketData), 1);
cluster(ismember(nuts1, southStates))   = "South";
cluster(ismember(nuts1, neutralStates)) = "Neutral";

marketData.Cluster = categorical(cluster);

% Originalzustände sichern
origMarketData         = marketData;
origCompensationFactor = compensationFactor;

%% === Zielgrößen ===
targetCapacityWind = inpTargetInvestments * 1000;
targetEnergyProduction = 384e9; % for 2040, corresponds to acutal reference yield model
tolEnergy = 0.002;

% === Set logit parameters ===
alpha   = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'alpha')};
beta    = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'beta')};
gamma   = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'gamma')};
lambda  = dcmParam{:, startsWith(dcmParam.Properties.VariableNames,'lambda')};

alpha = [alpha(1:3) 0 alpha(4:end)];
alpha = repmat(alpha', [nrDecisionMaker,1]);
beta  = repmat(beta,  [nrDecisionMaker*nrAlt,1]);
gamma = repmat(gamma,[nrDecisionMaker*nrAlt,1]);
lambda= repmat(lambda,[nrDecisionMaker*nrAlt,1]);

% ====================================================================
%  SENSITIVITY ANALYSIS: Reference-Yield Scaling
% ====================================================================

refValue  = 1.00;
scaleVec  = 0:0.005:1.7;   % Vector of scaling factors to test (can be adapted)
nScale    = numel(scaleVec);

% Pre-allocate results table including regional cluster capacities
resultsScaling = table( ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), zeros(nScale,1), ...
    zeros(nScale,1),  ...
    'VariableNames', {'scaleRefYield', ...
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
                      'rent_per_kWh'});              % annualized producer surplus per kWh

for iScale = 1:nScale

    % ---------------------------------------------------------------
    % Reset the state to baseline for each scaling factor
    % ---------------------------------------------------------------
    marketData         = origMarketData;
    compensationFactor = origCompensationFactor;
    compensation       = origCompensation;

    scaleRefYield = scaleVec(iScale);

    % ---------------------------------------------------------------
    % Apply refYield scaling: shifts compensationFactor around 1
    % ---------------------------------------------------------------
    compensationFactor = refValue + scaleRefYield * (compensationFactor - refValue);
    if ~optsRefYieldInInvestment
        compensationFactor = 1;
    end

    % ---------------------------------------------------------------
    % Apply CAPEX cost reduction for future years
    % ---------------------------------------------------------------
    marketData.investCost = applyCostReduction( ...
        marketData.investCost, ...
        paraRegioInvest.costReductionRate, ...
        paraRegioInvest.baseYear, ...
        paraRegioInvest.simYear);

    [estimatedNewCap, compensation, installedCapacity, totalGeneration_annual, countWhile, npv, npv_costs_total_lifetime, npv_revenues_total_lifetime] = ...
        runInvestmentLoop( ...
            marketData, production, turbineLifetime, interestRate, compensation, compensationFactor, ...
            paraOPEXfix, paraOPEXperkWhvar, optsNpvAdjust, explanatoryVariables, ...
            powerPotential, powerPotential_invest, powerTurbine, ...
            investFactor, availableSpace, ...
            optsRegioInvest.iterationLogicInvestment, ...
            targetCapacityWind, targetEnergyProduction, tolEnergy);

    % ---------------------------------------------------------------
    % Outputs for this scale factor
    % ---------------------------------------------------------------
    marketData.estimatedNewCapacity   = estimatedNewCap;
    marketData.estimatedTotalCapacity = estimatedNewCap + marketData.capacitiesInSimYear;

    % Weighted effective remuneration (real compensation)
    marketData.share_of_total_generation = ...
        ((marketData.estimatedTotalCapacity./marketData.powerTurbine) .* marketData.production) / ...
        sum((marketData.estimatedTotalCapacity ./marketData.powerTurbine) .* marketData.production);

    marketData.weighted_comp_share = marketData.share_of_total_generation .* compensationFactor;
    real_compensation = sum(marketData.weighted_comp_share) * compensation;  % EUR/kWh

    %%%%%%%%%%%%%%%%% calculate in annual values
    annualizationFactor = (interestRate * (1 + interestRate)^turbineLifetime) / ...
                          ((1 + interestRate)^turbineLifetime - 1);
    producer_surplus_one_year = sum((marketData.estimatedTotalCapacity ./ marketData.powerTurbine) .* (npv_revenues_total_lifetime .* annualizationFactor - npv_costs_total_lifetime .* annualizationFactor));
    % "fiscal costs"
    % Annual physical generation (kWh or MWh per year, depending on units)
    nTurbines               = marketData.estimatedTotalCapacity ./ marketData.powerTurbine;
    annualGeneration        = nTurbines .* marketData.production;
    totalGeneration_annual  = sum(annualGeneration);

    subsidy_one_year = totalGeneration_annual * real_compensation;

    totalWelfare = producer_surplus_one_year - subsidy_one_year;

    rent_per_kWh = producer_surplus_one_year ./ totalGeneration_annual;
    %%%%%%%%%%%%%%%%%
    
    % ---------------------------------------------------------------
    % Cluster capacities (North / South / Neutral) – total capacity
    % ---------------------------------------------------------------
    sumCap = groupsummary(marketData, "Cluster", "sum", "estimatedTotalCapacity");

    capNorth   = 0;
    capSouth   = 0;
    capNeutral = 0;

    idxN = sumCap.Cluster == 'North';
    if any(idxN), capNorth = sumCap.sum_estimatedTotalCapacity(idxN); end
    idxS = sumCap.Cluster == 'South';
    if any(idxS), capSouth = sumCap.sum_estimatedTotalCapacity(idxS); end
    idxZ = sumCap.Cluster == 'Neutral';
    if any(idxZ), capNeutral = sumCap.sum_estimatedTotalCapacity(idxZ); end

    capNorth_GW   = capNorth   / 1e6;
    capSouth_GW   = capSouth   / 1e6;
    capNeutral_GW = capNeutral / 1e6;

    % ---------------------------------------------------------------
    % Write all outputs into results table
    % ---------------------------------------------------------------
    resultsScaling.scaleRefYield(iScale)             = scaleRefYield;
    resultsScaling.installedCapacityGW(iScale)       = installedCapacity/1e6;
    resultsScaling.targetCapacityGW(iScale)          = targetCapacityWind/1e6;
    resultsScaling.finalComp_ctperkWh(iScale)        = compensation*100;
    resultsScaling.realComp_ctperkWh(iScale)         = real_compensation*100;
    resultsScaling.producerSurplus_bnEUR(iScale)     = producer_surplus_one_year/1e9;
    resultsScaling.subsidyPV_bnEUR(iScale)           = subsidy_one_year/1e9;
    resultsScaling.welfareNoGrid_bnEUR(iScale)       = totalWelfare/1e9;
    resultsScaling.totalGeneration_annual(iScale)    = totalGeneration_annual/1e9;
    %resultsScaling.totalGenerationPV(iScale)         = totalGeneration_PV;
    resultsScaling.capNorth_GW(iScale)               = capNorth_GW;
    resultsScaling.capSouth_GW(iScale)               = capSouth_GW;
    resultsScaling.capNeutral_GW(iScale)             = capNeutral_GW;
    resultsScaling.rent_per_kWh(iScale)              = rent_per_kWh;
end

%% ================================================================
%  Plot results over reference-yield scaling
% ================================================================

x = resultsScaling.scaleRefYield;

figure('Position',[100 100 1000 800]);

% ------------------------------------------------
% 1. Welfare
% ------------------------------------------------
subplot(3,2,1);
plot(x, resultsScaling.welfareNoGrid_bnEUR, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('Welfare (bn EUR / year)');
title('Welfare');

% ------------------------------------------------
% 2. Remuneration levels
% ------------------------------------------------
subplot(3,2,2);
plot(x, resultsScaling.finalComp_ctperkWh, 'LineWidth', 2); hold on;
plot(x, resultsScaling.realComp_ctperkWh, '--', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('ct/kWh');
title('Remuneration');
legend('Required (auction price)', 'Real (production-weighted)', 'Location','best');

% ------------------------------------------------
% 3. Producer surplus & subsidy
% ------------------------------------------------
subplot(3,2,3);
plot(x, resultsScaling.producerSurplus_bnEUR, 'LineWidth', 2); hold on;
plot(x, resultsScaling.subsidyPV_bnEUR, '--', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('bn EUR / year');
title('Rents and fiscal costs');
legend('Producer surplus', 'Subsidy', 'Location','best');

% ------------------------------------------------
% 4. Total generation
% ------------------------------------------------
subplot(3,2,4);
plot(x, resultsScaling.totalGeneration_annual, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('TWh / year');
title('Total annual generation');

% ------------------------------------------------
% 5. Regional capacities
% ------------------------------------------------
subplot(3,2,5);
plot(x, resultsScaling.capNorth_GW, 'LineWidth', 2); hold on;
plot(x, resultsScaling.capNeutral_GW, '--', 'LineWidth', 2);
plot(x, resultsScaling.capSouth_GW, ':', 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('GW');
title('Regional capacities');
legend('North','Neutral','South','Location','best');

% ------------------------------------------------
% 6. Rent per kWh
% ------------------------------------------------
subplot(3,2,6);
plot(x, resultsScaling.rent_per_kWh, 'LineWidth', 2);
grid on;
xlabel('Scaling factor s');
ylabel('EUR / kWh');
title('Annualized rent per kWh');

sgtitle('Sensitivity to reference-yield scaling', 'FontWeight','bold');


%% ===================================================================
%  FINAL RUN for standard result: scaleRefYield = 1
% ====================================================================
marketData         = origMarketData;
compensationFactor = origCompensationFactor;
compensation       = origCompensation;

scaleRefYield_final = 1.0;
refValue            = 1.0;

compensationFactor = refValue + scaleRefYield_final * (compensationFactor - refValue);
if ~optsRefYieldInInvestment
    compensationFactor = 1;
end

marketData.investCost = applyCostReduction( ...
    marketData.investCost, ...
    paraRegioInvest.costReductionRate, ...
    paraRegioInvest.baseYear, ...
    paraRegioInvest.simYear);

explVar_initNPV(:,1) = calcNetPresentValue(marketData.investCost,production,turbineLifetime,interestRate,compensation,compensationFactor, paraOPEXfix, paraOPEXperkWhvar);
meanNPV = mean(explVar_initNPV(:,1));
disp(['Average Net Present Value considering a Compensation of ',num2str(compensation.*100),'ct/kWh : ',num2str(floor(meanNPV)) , '€'])
if numel(explanatoryVariables)>1
    explVar_initNPV(:,2:numel(explanatoryVariables)) = marketData{:,explanatoryVariables(2:end)};
end

[estimatedNewCap, compensation, installedCapacity, totalGeneration_annual, countWhile, npv] = ...
    runInvestmentLoop( ...
        marketData, production, turbineLifetime, interestRate, compensation, compensationFactor, ...
        paraOPEXfix, paraOPEXperkWhvar, optsNpvAdjust, explanatoryVariables, ...
        powerPotential, powerPotential_invest, powerTurbine, ...
        investFactor, availableSpace, ...
        optsRegioInvest.iterationLogicInvestment, ...
        targetCapacityWind, targetEnergyProduction, tolEnergy);

marketData.estimatedNewCapacity   = estimatedNewCap;
marketData.estimatedTotalCapacity = estimatedNewCap + marketData.capacitiesInSimYear;

marketData.share_of_total_generation = ((marketData.estimatedTotalCapacity./marketData.powerTurbine) .* marketData.production) / ...
    sum((marketData.estimatedTotalCapacity ./marketData.powerTurbine) .* marketData.production);
marketData.weighted_comp_share = marketData.share_of_total_generation .* compensationFactor;
real_compensation = sum(marketData.weighted_comp_share) * compensation;

disp(['Estimated capacity: ', num2str(installedCapacity/1e6), ' GW'])
disp(['Target capacity: ', num2str(targetCapacityWind/1e6), ' GW'])
disp(['Final compensation: ', num2str(compensation*100), ' ct/kWh'])
disp(['Real avg. compensation: ', num2str(real_compensation*100), ' ct/kWh'])
disp(['Total NPV of the power plants is: ' num2str(sum((marketData.estimatedTotalCapacity./marketData.powerTurbine).*npv)/10^9) ' bn. EUR'])

newInvestments     = sum(marketData.estimatedNewCapacity);
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
    yearsToTarget   = targetYear - baseYear;
    reductionFactor = (1 - annualReductionRate) .^ yearsToTarget;
    updatedInvestCost = investCost .* reductionFactor;
end

%% FUNCTION
function [estimatedNewCap, compensation, installedCapacity, totalGeneration_annual, countWhile, npv, npv_costs_total_lifetime, npv_revenues_total_lifetime] = ...
    runInvestmentLoop( ...
        marketData, production, turbineLifetime, interestRate, compensation, compensationFactor, ...
        paraOPEXfix, paraOPEXperkWhvar, optsNpvAdjust, explanatoryVariables, ...
        powerPotential, powerPotential_invest, powerTurbine, ...
        investFactor, availableSpace, ...
        iterationLogicInvestment, ...
        targetCapacityWind, targetEnergyProduction, tolEnergy)

    %RUNINVESTMENTLOOP Runs the endogenous-compensation loop to hit either a
    %capacity target or an annual energy target, depending on iterationLogicInvestment.
    %
    % Returns:
    %   estimatedNewCap         - vector (MW) new capacity additions by row
    %   compensation            - final compensation level (EUR/kWh)
    %   installedCapacity       - total installed capacity after additions (MW)
    %   totalGeneration_annual  - annual generation after additions (same unit as production)
    %   countWhile              - iteration counter

    maxIter = 100000;
    countWhile = 0;

    % Initialize
    estimatedNewCap = zeros(height(marketData), 1);

    % Helper to compute annual generation from total capacity
    function gen = computeAnnualGeneration(totalCap)
        nTurb = totalCap ./ powerTurbine;
        gen   = sum(nTurb .* marketData.production);
    end

    % Initialize reporting variables
    installedCapacity      = sum(marketData.capacitiesInSimYear);
    totalGeneration_annual = computeAnnualGeneration(marketData.capacitiesInSimYear);

    if strcmp(iterationLogicInvestment, 'capacity')

        while ( (installedCapacity <= targetCapacityWind || ...
                 installedCapacity / targetCapacityWind > (1 + tolEnergy)) ...
                && countWhile < maxIter )

            % NPV
            [npv, npv_costs_total_lifetime, npv_revenues_total_lifetime] = calcNetPresentValue( ...
                marketData.investCost, production, turbineLifetime, ...
                interestRate, compensation, compensationFactor, ...
                paraOPEXfix, paraOPEXperkWhvar);

            % Optional NPV adjustment
            switch optsNpvAdjust
                case 'NPVperInvestCosts'
                    npv2 = npv ./ marketData.investCost;
                case 'NPVperArea'
                    npv2 = npv ./ (marketData.powerTurbine ./ powerPotential);
                case 'none'
                    npv2 = npv;
                otherwise
                    error('Invalid optsNpvAdjust value.');
            end
            npv_transformed = transformNPV(npv2);

            % Expl. vars
            explVar = zeros(height(marketData), numel(explanatoryVariables));
            explVar(:,1) = npv_transformed;
            if numel(explanatoryVariables) > 1
                explVar(:,2:end) = marketData{:,explanatoryVariables(2:end)};
            end

            % Market shares
            predictedMarketShare = binomNestLog(explVar);

            % Capacity update (MW)
            estimatedNewCap = floor(investFactor .* availableSpace .* ...
                                    powerPotential_invest ./ powerTurbine .* ...
                                    predictedMarketShare) .* powerTurbine;

            installedCapacity = sum(estimatedNewCap + marketData.capacitiesInSimYear);

            % Adjust compensation
            if installedCapacity / targetCapacityWind > (1 + tolEnergy)
                compensation = compensation - 0.00005;
            else
                compensation = compensation + 0.001;
            end

            countWhile = countWhile + 1;
        end

        % compute energy for reporting
        totalGeneration_annual = computeAnnualGeneration(estimatedNewCap + marketData.capacitiesInSimYear);

    elseif strcmp(iterationLogicInvestment, 'energy')

        % initialize based on baseline
        totalGeneration_annual = computeAnnualGeneration(marketData.capacitiesInSimYear);

        while ( (totalGeneration_annual <= targetEnergyProduction || ...
                 totalGeneration_annual / targetEnergyProduction > (1 + tolEnergy)) ...
                && countWhile < maxIter )

            % NPV
            [npv, npv_costs_total_lifetime, npv_revenues_total_lifetime] = calcNetPresentValue( ...
                marketData.investCost, production, turbineLifetime, ...
                interestRate, compensation, compensationFactor, ...
                paraOPEXfix, paraOPEXperkWhvar);

            switch optsNpvAdjust
                case 'NPVperInvestCosts'
                    npv2 = npv ./ marketData.investCost;
                case 'NPVperArea'
                    npv2 = npv ./ (marketData.powerTurbine ./ powerPotential);
                case 'none'
                    npv2 = npv;
                otherwise
                    error('Invalid optsNpvAdjust value.');
            end
            npv_transformed = transformNPV(npv2);

            explVar = zeros(height(marketData), numel(explanatoryVariables));
            explVar(:,1) = npv_transformed;
            if numel(explanatoryVariables) > 1
                explVar(:,2:end) = marketData{:,explanatoryVariables(2:end)};
            end

            predictedMarketShare = binomNestLog(explVar);

            estimatedNewCap = floor(investFactor .* availableSpace .* ...
                                    powerPotential_invest ./ powerTurbine .* ...
                                    predictedMarketShare) .* powerTurbine;

            % Update energy implied by total capacity
            totalCap = estimatedNewCap + marketData.capacitiesInSimYear;
            totalGeneration_annual = computeAnnualGeneration(totalCap);

            % Adjust compensation
            if totalGeneration_annual / targetEnergyProduction > (1 + tolEnergy)
                compensation = compensation - 0.00005;
            else
                compensation = compensation + 0.001;
            end

            countWhile = countWhile + 1;
        end

        installedCapacity = sum(estimatedNewCap + marketData.capacitiesInSimYear);

    else
        error('Unknown optsRegioInvest.iterationLogicInvestment. Use ''capacity'' or ''energy''.');
    end
end


