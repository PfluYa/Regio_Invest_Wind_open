%% Plot relative capacity difference: without RYM vs. with RYM
% This script compares two simulation results and visualizes the relative
% difference in installed capacity at NUTS3 level.
%
% Color logic:
%   - red   -> lower capacity without RYM
%   - green -> higher capacity without RYM
%   - grey  -> missing values
%
% Visualization is clipped to +/- 50%.

%% ------------------------------------------------------------------------
% 1. LOAD RESULTS
% -------------------------------------------------------------------------

withRYM = load('2025-08-12_results_nuts3_base2024_target2040_expCase1.mat');
withoutRYM = load('2025-08-12_results_nuts3_base2024_target2040_expCase1_wo_refYield.mat');

withRYM = withRYM.selected_data;
withoutRYM = withoutRYM.selected_data;

%% ------------------------------------------------------------------------
% 2. ALIGN DATA BY NUTS3 ID
% -------------------------------------------------------------------------

[isMatched, idxWith] = ismember(withoutRYM.nutsID, withRYM.nutsID);

if any(~isMatched)
    error('Mismatch in NUTS3 IDs between the two result files.');
end

withoutRYM.capacityTotal_withRYM = withRYM.capacityTotal(idxWith);

%% ------------------------------------------------------------------------
% 3. CALCULATE RELATIVE DIFFERENCE
% -------------------------------------------------------------------------

withoutRYM.perc_diff_to_with = ...
    (withoutRYM.capacityTotal - withoutRYM.capacityTotal_withRYM) ./ ...
     withoutRYM.capacityTotal_withRYM;

values = withoutRYM.perc_diff_to_with;

%% ------------------------------------------------------------------------
% 4. PREPARE VISUALIZATION
% -------------------------------------------------------------------------

% Clip values to +/- 50% for readability
clipLimit = 0.5;
valuesClipped = min(max(values, -clipLimit), clipLimit);

% Original color mapping:
% - negative values -> more red
% - positive values -> more green
% - blue channel fixed to keep colors bright
colorMapFun = @(v) [ ...
    1 - max(0, v / clipLimit), ...
    1 - max(0, -v / clipLimit), ...
    1];

%% ------------------------------------------------------------------------
% 5. PLOT MAP
% -------------------------------------------------------------------------

figure();
hold on;
title('');

for i = 1:height(withoutRYM)
    regionShape = polyshape(withoutRYM.lonPoly{i}, withoutRYM.latPoly{i});
    valueHere = valuesClipped(i);

    if isnan(values(i))
        faceCol = [0.5 0.5 0.5];
    else
        faceCol = colorMapFun(valueHere);
    end

    plot(regionShape, ...
        'FaceColor', faceCol, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', [0.7 0.7 0.7]);
end

%% ------------------------------------------------------------------------
% 6. CUSTOM COLORBAR
% -------------------------------------------------------------------------

nSteps = 100;
vals = linspace(-clipLimit, clipLimit, nSteps);
colors = zeros(nSteps, 3);

for i = 1:nSteps
    colors(i, :) = colorMapFun(vals(i));
end

colormap(colors);

cb = colorbar;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'-50%', '0%', '+50%'};
cb.Label.String = '';

%% ------------------------------------------------------------------------
% 7. AXIS SETTINGS
% -------------------------------------------------------------------------

axis off;
axis equal;
xlim([min(withoutRYM.lonPoint) - 1, max(withoutRYM.lonPoint) + 1]);
ylim([min(withoutRYM.latPoint) - 1, max(withoutRYM.latPoint) + 1]);
hold off;