function [netPresentValue, npv_costs_total_lifetime, npv_revenues_total_lifetime] = calcNetPresentValue( ...
    inpInvestCost, inpProduction, inpLifetime, inpInterestRate, ...
    inpCompensation, inpCompensationFactor, inpOPEXfix, inpOPEXperkWhvar)

    % Rename inputs for readability
    investCost           = inpInvestCost;           % EUR (one-time)
    production           = inpProduction;           % kWh/year (or MWh/year)
    interestRate         = inpInterestRate;
    lifetime             = inpLifetime;
    compensation         = inpCompensation;         % EUR/kWh
    compCorrectionFactor = inpCompensationFactor;   % unitless
    opexFix              = inpOPEXfix;              % EUR/year
    opexVar              = inpOPEXperkWhvar;        % EUR/kWh

    % ------------------------------------------------------------
    % Annual components
    % ------------------------------------------------------------

    % Annual revenues (EUR/year)
    annualRevenues = compensation .* compCorrectionFactor .* production;

    % Annual operating costs (EUR/year)
    annualOpex = opexVar .* production + opexFix;

    % ------------------------------------------------------------
    % Discounting over lifetime
    % ------------------------------------------------------------

    % Annuity factor for present value calculation
    annuityFactor = ((1 + interestRate)^lifetime - 1) / ...
                    ((1 + interestRate)^lifetime * interestRate);

    % Present value of revenues over lifetime (EUR)
    npv_revenues_total_lifetime = annualRevenues * annuityFactor;

    % Present value of operating costs over lifetime (EUR)
    npv_opex_total_lifetime = annualOpex * annuityFactor;

    % Present value of total costs = CAPEX + discounted OPEX (EUR)
    npv_costs_total_lifetime = investCost + npv_opex_total_lifetime;

    % ------------------------------------------------------------
    % Net present value
    % ------------------------------------------------------------

    netPresentValue = npv_revenues_total_lifetime - npv_costs_total_lifetime;
end
