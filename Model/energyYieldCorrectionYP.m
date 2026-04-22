function correctedYield = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter)
% ENERGYYIELDCORRECTIONYP
% Apply the empirical yield correction proposed in Pflugfelder, Kramer,
% and Weber (2024).
%
% DESCRIPTION
% This function adjusts modeled annual production values using an empirical
% correction based on full-load-hour deviations. The correction depends on:
%   - distance to sea,
%   - hub height,
%   - rotor diameter,
%   - turbine capacity.
%
% The correction is implemented as:
%   1. predict the deviation in full-load hours from the regression model,
%   2. adjust observed/model-based full-load hours,
%   3. scale annual production accordingly.
%
% INPUTS
%   filenameDistSea   : file containing distance-to-sea information
%   regioWindOnshore  : table with regional turbine data; must contain
%                       at least:
%                           - turbineType
%                           - powerTurbine
%                           - fullLoadHours
%                           - production
%   turbineParameter  : table with turbine specifications; must contain
%                           - turbineType
%                           - diameter
%                           - hubHeight
%
% OUTPUT
%   correctedYield    : corrected annual production [same unit as input
%                       production]
%
% NOTE
% The current baseline in the main script sets:
%   optsRegioInvest.windOnshoreYieldCorrectionYP = true
%
% REFERENCE
% Pflugfelder, Kramer, Weber (2024)

%% ------------------------------------------------------------------------
% 1. MERGE TURBINE CHARACTERISTICS
% -------------------------------------------------------------------------

regioWindOnshore = join( ...
    regioWindOnshore, ...
    turbineParameter, ...
    'Keys', 'turbineType', ...
    'RightVariables', {'diameter', 'hubHeight'});

%% ------------------------------------------------------------------------
% 2. LOAD DISTANCE-TO-SEA DATA
% -------------------------------------------------------------------------

distSeaTable = readtable(filenameDistSea);

if ~ismember('HubDist', distSeaTable.Properties.VariableNames)
    error('Column "HubDist" not found in distance-to-sea file: %s', filenameDistSea);
end

distSea = distSeaTable.HubDist;

if numel(distSea) ~= height(regioWindOnshore)
    error(['Distance-to-sea vector length does not match the number of ', ...
           'rows in regioWindOnshore.']);
end

%% ------------------------------------------------------------------------
% 3. EXTRACT REQUIRED VARIABLES
% -------------------------------------------------------------------------

hubHeight = regioWindOnshore.hubHeight;
diameter  = regioWindOnshore.diameter;
capacity  = regioWindOnshore.powerTurbine;
flhActual = regioWindOnshore.fullLoadHours;

%% ------------------------------------------------------------------------
% 4. PREDICT FULL-LOAD-HOUR DEVIATION
% -------------------------------------------------------------------------
% Regression-based correction from Pflugfelder et al. (2024)

flhDeviation = ...
    844.67225 ...
    - 1.5175  .* distSea ...
    - 1.87725 .* hubHeight ...
    - 2.2045  .* diameter ...
    + 0.09875 .* capacity;

%% ------------------------------------------------------------------------
% 5. APPLY CORRECTION
% -------------------------------------------------------------------------

flhCorrected = flhActual - flhDeviation;
correctionFactor = flhCorrected ./ flhActual;

% Prevent invalid values if flhActual is zero or missing
invalidIdx = isnan(correctionFactor) | isinf(correctionFactor);
correctionFactor(invalidIdx) = 0;

correctedYield = regioWindOnshore.production .* correctionFactor;

end