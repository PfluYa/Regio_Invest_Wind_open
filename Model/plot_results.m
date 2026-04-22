function plot_results(resultsWindOnshoreTurbinesUnstack, regioDataWind)
% PLOT_RESULTS
% Visualize simulated onshore wind deployment across turbine types and
% regions.
%
% INPUTS
%   resultsWindOnshoreTurbinesUnstack : table with installed capacities by
%                                       NUTS3 region and turbine type
%   regioDataWind                     : regional table with polygons and
%                                       capacity indicators
%
% OUTPUT
%   This function creates:
%     1. a bar chart of installed capacity by turbine type,
%     2. a map of total installed capacity,
%     3. a map of installed capacity per km^2,
%     4. a map of installed capacity per available wind area.
%
% NOTES
% - Capacities are internally assumed to be stored in kW and converted to
%   GW for the bar chart.
% - Map plots use alpha transparency to represent intensity.

%% ------------------------------------------------------------------------
% 1. INSTALLED CAPACITY BY TURBINE TYPE
% -------------------------------------------------------------------------

turbineColumns = resultsWindOnshoreTurbinesUnstack.Properties.VariableNames(2:end);
turbineCapacity = table2array(resultsWindOnshoreTurbinesUnstack(:, turbineColumns));

figure();
bar(sum(turbineCapacity, 1) ./ 1e6);

% Create readable x-axis labels from column names if possible
defaultLabels = {'Type 1', 'Type 2', 'Type 3', 'Type 4', ...
                 'Type 5', 'Type 6', 'Type 7', 'Type 8'};

if numel(turbineColumns) == numel(defaultLabels)
    set(gca, 'XTickLabel', defaultLabels);
else
    set(gca, 'XTickLabel', turbineColumns);
end

ylabel('Installed capacity [GW]');
title('Installed wind capacity by turbine type');
grid on;

%% ------------------------------------------------------------------------
% 2. MAP: TOTAL INSTALLED CAPACITY
% -------------------------------------------------------------------------

valuesTotalCapacity = regioDataWind.capacityTotal;
valuesTotalCapacity(valuesTotalCapacity <= 500) = NaN;

plotCapacityMap( ...
    regioDataWind, ...
    valuesTotalCapacity, ...
    'Installed capacity by NUTS3 region in target year', ...
    'MW');

%% ------------------------------------------------------------------------
% 3. MAP: INSTALLED CAPACITY PER TOTAL AREA
% -------------------------------------------------------------------------

valuesCapacityPerKm2 = regioDataWind.capPerKm2;
valuesCapacityPerKm2(valuesCapacityPerKm2 <= 0.0001) = NaN;

plotCapacityMap( ...
    regioDataWind, ...
    valuesCapacityPerKm2, ...
    'Installed capacity per km^2 by NUTS3 region in target year', ...
    'MW/km^2');

%% ------------------------------------------------------------------------
% 4. MAP: INSTALLED CAPACITY PER AVAILABLE WIND AREA
% -------------------------------------------------------------------------

valuesCapacityPerKm2Available = regioDataWind.capPerKm2avail;
valuesCapacityPerKm2Available(valuesCapacityPerKm2Available <= 0.0001) = NaN;

plotCapacityMap( ...
    regioDataWind, ...
    valuesCapacityPerKm2Available, ...
    'Installed capacity per available wind area by NUTS3 region in target year', ...
    'MW/km^2');

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function plotCapacityMap(regioData, values, titleStr, unitStr)
% PLOTCAPACITYMAP
% Plot a choropleth-style regional map using polygon transparency.
%
% INPUTS
%   regioData : table containing
%       - lonPoly, latPoly : polygon coordinates
%       - lonPoint, latPoint : point coordinates used for axis limits
%   values    : numeric vector of values to visualize
%   titleStr  : plot title
%   unitStr   : unit label shown in legend text
%
% NOTE
% Regions with NaN values are drawn in black.

    maxVal = max(values, [], 'omitnan');
    minVal = min(values, [], 'omitnan');

    if all(isnan(values))
        warning('All map values are NaN. Plot will show only missing regions.');
        alphaValues = values;
    else
        alphaValues = rescale(values, 0.1, 1.0);
    end

    figure();
    hold on;

    warning('off', 'all');

    for regionIdx = 1:height(regioData)
        regionShape = polyshape( ...
            regioData.lonPoly{regionIdx}, ...
            regioData.latPoly{regionIdx});

        if ~isnan(alphaValues(regionIdx))
            plot(regionShape, ...
                'FaceColor', 'blue', ...
                'FaceAlpha', alphaValues(regionIdx), ...
                'EdgeColor', 'none');
        else
            plot(regionShape, ...
                'FaceColor', 'black', ...
                'FaceAlpha', 1, ...
                'EdgeColor', 'none');
        end
    end

    warning('on', 'all');

    title(titleStr);

    % Minimal legend via dummy objects
    hMax = plot(nan, nan, 's', ...
        'MarkerFaceColor', 'blue', ...
        'MarkerEdgeColor', 'blue', ...
        'MarkerSize', 8);
    hMin = plot(nan, nan, 's', ...
        'MarkerFaceColor', [0.7 0.7 1.0], ...
        'MarkerEdgeColor', [0.7 0.7 1.0], ...
        'MarkerSize', 8);
    hNaN = plot(nan, nan, 's', ...
        'MarkerFaceColor', 'black', ...
        'MarkerEdgeColor', 'black', ...
        'MarkerSize', 8);

    legend( ...
        [hMax, hMin, hNaN], ...
        {['max = ', num2str(round(maxVal, 2)), ' ', unitStr], ...
         ['min = ', num2str(round(minVal, 2)), ' ', unitStr], ...
         'no / negligible value'}, ...
        'Location', 'best');

    xlim([min(regioData.lonPoint) - 1, max(regioData.lonPoint) + 1]);
    ylim([min(regioData.latPoint) - 1, max(regioData.latPoint) + 1]);

    axis off;
    axis equal;
    hold off;
end