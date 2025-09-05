function yield_new = energyYieldCorrectionYP(filenameDistSea, regioWindOnshore, turbineParameter)
% ENERGYYIELDCORRECTIONYP computes energy yield correction factors based on
% turbine characteristics and distance to sea, according to
% Pflugfelder, Kramer, Weber (2024).
%
% INPUTS:
%   filenameDistSea      - Path to file containing sea distance data
%   regioWindOnshore     - Table of regional turbine data (must include turbineType)
%   turbineParameter     - Table of turbine specs (diameter, hubHeight)
%
% OUTPUT:
%   correctionFactor     - Vector of correction factors per turbine

    % Merge turbine specs into regional dataset
    regioWindOnshore = join(regioWindOnshore, turbineParameter, ...
                            'Keys', 'turbineType', ...
                            'RightVariables', {'diameter', 'hubHeight'});

    % Load sea distance values
    distSea = readtable(filenameDistSea);
    distSea = distSea.HubDist;

    % Extract required parameters
    hubHeight   = regioWindOnshore.hubHeight;
    diameter    = regioWindOnshore.diameter;
    capacity    = regioWindOnshore.powerTurbine;
    flhActual   = regioWindOnshore.fullLoadHours;

    % === Calculate predicted FLH deviation based on regression model ===
    % Source: Pflugfelder et al. (2024)
    flhDeviation = ...
        844.67225 ...
        - 1.5175  * distSea ...
        - 1.87725 * hubHeight ...
        - 2.2045  * diameter ...
        + 0.09875 * capacity;

    % === Apply correction ===
    flhCorrected = flhActual - flhDeviation;
    correctionFactor = flhCorrected ./ flhActual;

    yield_new = regioWindOnshore.production .* correctionFactor;

end
