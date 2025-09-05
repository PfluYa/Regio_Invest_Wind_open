function [relativeWind_BG, relativeWind_BWE, relativeWind_UB21] = adaptPercentageWind(baseTable)
% adaptPercentageWind adjusts available wind space per NUTS1 region 
% according to different policy studies (WindBG, BWE, UB21).
%
% Input:
%   - baseTable: Table with NUTS3 regions, total area, and relative available wind space
%
% Output:
%   - relativeWind_BG: Adjusted wind space based on WindBG law
%   - relativeWind_BWE: Adjusted wind space based on BWE study
%   - relativeWind_UB21: Adjusted wind space based on UBA (UB21) study

% Extract NUTS1 region code from NUTS3
baseTable.nuts1 = extractBefore(baseTable.nutsID, 4);

% --- WindBG scenario ---
relativeWind_BG = applyCorrection(baseTable, [
    0.018, 0.018, 0.005, 0.022, 0.005, 0.005, 0.022, 0.021, ...
    0.022, 0.018, 0.022, 0.018, 0.020, 0.022, 0.020, 0.022]);

% --- BWE scenario ---
relativeWind_BWE = applyCorrection(baseTable, [
    0.036, 0.046, 0.007, 0.083, 0.010, 0.006, 0.052, 0.065, ...
    0.078, 0.027, 0.046, 0.027, 0.054, 0.111, 0.052, 0.096]);

% --- UB21 scenario ---
relativeWind_UB21 = applyCorrection(baseTable, [
    0.005, 0.007, 0.000, 0.008, 0.008, 0.002, 0.018, 0.002, ...
    0.010, 0.010, 0.015, 0.018, 0.002, 0.008, 0.020, 0.030]);

end

function corrected = applyCorrection(table, windAreaTargets)
% Applies correction factors to relative available wind space
% to match target area shares per NUTS1 region

    nuts1Codes = unique(table.nuts1, 'stable');
    corrected = zeros(height(table),1);

    for j = 1:length(nuts1Codes)
        idx = strcmp(table.nuts1, nuts1Codes{j});

        totalArea = sum(table.totalArea(idx));
        windArea = sum(table.totalArea(idx) .* table.relativeAvailableWindSpace(idx));
        currentShare = windArea / totalArea;

        correctionFactor = windAreaTargets(j) / currentShare;
        corrected(idx) = table.relativeAvailableWindSpace(idx) * correctionFactor;
    end
end
