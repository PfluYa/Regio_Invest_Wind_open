function [netPresentValue, npvCostsTotalLifetime, npvRevenuesTotalLifetime] = ...
    calcNetPresentValue( ...
        inpInvestCost, ...
        inpProduction, ...
        inpLifetime, ...
        inpInterestRate, ...
        inpCompensation, ...
        inpCompensationFactor, ...
        inpOPEXfix, ...
        inpOPEXperkWhvar)
% CALCNETPRESENTVALUE
% Computes the net present value (NPV) of a wind turbine investment.
%
% DESCRIPTION
% The function calculates the discounted lifetime revenues and costs of a
% turbine-location combination and returns the resulting net present value.
%
% INPUTS (current defaults from main script)
%   inpInvestCost           : investment cost [EUR]
%   inpProduction           : annual electricity production [kWh/year]
%   inpLifetime             : turbine lifetime [years] (default: 22)
%   inpInterestRate         : discount rate [-] (default: 0.05)
%   inpCompensation         : remuneration level [EUR/kWh] (default: 0.080)
%   inpCompensationFactor   : location-specific correction factor [-]
%   inpOPEXfix              : fixed O&M costs [EUR/year] (default: 15 EUR/kW * capacity)
%   inpOPEXperkWhvar        : variable O&M costs [EUR/kWh] (default: 0.0075)
%
% OUTPUTS
%   netPresentValue             : NPV [EUR]
%   npvCostsTotalLifetime       : discounted total costs (CAPEX + OPEX) [EUR]
%   npvRevenuesTotalLifetime    : discounted total revenues [EUR]

%% ------------------------------------------------------------------------
% 1. ASSIGN INPUTS (for readability)
% -------------------------------------------------------------------------

investCost           = inpInvestCost;           % EUR
production           = inpProduction;           % kWh/year
interestRate         = inpInterestRate;         % -
lifetime             = inpLifetime;             % years
compensation         = inpCompensation;         % EUR/kWh
compensationFactor   = inpCompensationFactor;   % -
opexFix              = inpOPEXfix;              % EUR/year
opexVar              = inpOPEXperkWhvar;        % EUR/kWh

%% ------------------------------------------------------------------------
% 2. ANNUAL CASH FLOWS
% -------------------------------------------------------------------------

% Annual revenues (EUR/year)
annualRevenues = compensation .* compensationFactor .* production;

% Annual operating costs (EUR/year)
annualOpex = opexVar .* production + opexFix;

%% ------------------------------------------------------------------------
% 3. DISCOUNTING OVER LIFETIME
% -------------------------------------------------------------------------

% Present value factor for an annuity (sum of discounted annual flows)
if interestRate == 0
    % Edge case: no discounting
    annuityFactor = lifetime;
else
    annuityFactor = ((1 + interestRate)^lifetime - 1) / ...
                    ((1 + interestRate)^lifetime * interestRate);
end

% Discounted lifetime revenues (EUR)
npvRevenuesTotalLifetime = annualRevenues .* annuityFactor;

% Discounted lifetime OPEX (EUR)
npvOpexTotalLifetime = annualOpex .* annuityFactor;

% Total discounted costs (CAPEX + OPEX)
npvCostsTotalLifetime = investCost + npvOpexTotalLifetime;

%% ------------------------------------------------------------------------
% 4. NET PRESENT VALUE
% -------------------------------------------------------------------------

netPresentValue = npvRevenuesTotalLifetime - npvCostsTotalLifetime;

end