withRYM    = load('2025-08-12_results_nuts3_base2024_target2040_expCase1').selected_data;
withoutRYM = load('2025-08-12_results_nuts3_base2024_target2040_expCase1_wo_refYield').selected_data;

withoutRYM.capacityTotal_withRYM = withRYM.capacityTotal;
withoutRYM.perc_diff_to_with = (withoutRYM.capacityTotal - withRYM.capacityTotal) ./ withRYM.capacityTotal;

plot_perc_diff_map(withoutRYM);


function plot_perc_diff_map(regioData)
% PLOT_PERC_DIFF_MAP Visualizes percent difference between scenarios (red/green)
%   Colors:
%   - Red for < 0 (less capacity w/o RYM)
%   - Green for > 0 (more capacity w/o RYM)
%   - Grey for NaN

values = regioData.perc_diff_to_with;

% Limit to ±0.5 for visualization
valuesClipped = min(max(values, -0.5), 0.5);  % clip to [-0.5, 0.5]

% Color mapping function
colorMap = @(v) [ ...
    1 - max(0, v / 0.5), ...  % Red fades from 1 to 0 as v goes from -0.5 to 0.5
    1 - max(0, -v / 0.5), ... % Green fades from 0 to 1 as v goes from -0.5 to 0.5
    1];                       % Blue stays high to keep background bright


% Begin plot
figure(); hold on;
title('');%'Relative capacity difference (w/o RYM vs. with RYM)');

for i = 1:height(regioData)
    shape = polyshape(regioData.lonPoly{i}, regioData.latPoly{i});
    val = valuesClipped(i);
    if isnan(values(i))
        faceCol = [0.5 0.5 0.5];  % gray for NaN
    else
        faceCol = colorMap(val);
    end
    plot(shape, 'FaceColor', faceCol, 'FaceAlpha', 0.8, 'EdgeColor', [.7 .7 .7]);
end

% Add custom legend with colorbar
% Build colormap matching colorMap logic
nSteps = 100;
vals = linspace(-0.5, 0.5, nSteps);
colors = zeros(nSteps, 3);
for i = 1:nSteps
    colors(i, :) = colorMap(vals(i));
end
colormap(colors);

cb = colorbar;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'-50%', '0%', '+50%'};
cb.Label.String = '';%'Capacity difference (w/o RYM vs. with RYM)';

axis off
xlim([min(regioData.lonPoint)-1, max(regioData.lonPoint)+1]);
ylim([min(regioData.latPoint)-1, max(regioData.latPoint)+1]);
end
