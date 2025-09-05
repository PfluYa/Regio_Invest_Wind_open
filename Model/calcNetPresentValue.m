function [netPresentValue] = calcNetPresentValue( ...
    inpInvestCost, inpProduction, inpLifetime, inpInterestRate, ...
    inpCompensation, inpCompensationFactor, inpOPEXfix, inpOPEXperkWhvar)

    % Rename inputs for readability
    investCost           = inpInvestCost;
    production           = inpProduction;
    interestRate         = inpInterestRate;
    lifetime             = inpLifetime;
    compensation         = inpCompensation;
    compCorrectionFactor = inpCompensationFactor;
    opexFix              = inpOPEXfix;
    opexVar              = inpOPEXperkWhvar;

    % Annual profit = revenue - variable OPEX - fixed OPEX
    profitInfeed = ...
        compensation .* compCorrectionFactor .* production ...
        - opexVar .* production ...
        - opexFix;

    % Annuity factor for present value calculation
    annuityFactor = ((1 + interestRate) ^ lifetime - 1) / ...
                    ((1 + interestRate) ^ lifetime * interestRate);

    % NPV formula: NPV = -CAPEX + Annual Profit * Annuity Factor
    netPresentValue = -investCost + profitInfeed * annuityFactor;
end