function plot_results(resultsWindOnshoreTurbinesUnstack, regioDataWind)
% PLOT_RESULTS Visualizes onshore wind capacity results across turbine types and regions.
%
% INPUTS:
%   resultsWindOnshoreTurbinesUnstack - Table with installed capacities for each turbine type
%   regioDataWind                     - Table with geographic and capacity data per NUTS3 region

%% 1. Bar plot: Installed capacity per turbine type
turbineData = table2array(resultsWindOnshoreTurbinesUnstack(:, 2:9));
bar(sum(turbineData) / 1e6);
set(gca, "XTickLabel", {'Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5', 'Type 6', 'Type 7', 'Type 8'});
ylabel('Installed capacity [GW]');
title('Installed Wind Capacity by Turbine Type');

%% 2. Map plot: Total installed capacity per NUTS3
resultsToPlot = regioDataWind.capacityTotal;
totalCapacityGW = sum(resultsToPlot) / 1e6;
resultsToPlot(resultsToPlot <= 500) = NaN;

% disp(['Total Capacity: ', num2str(totalCapacityGW), ' GW']);

plotCapacityMap(regioDataWind, resultsToPlot, 'Installed capacity per NUTS 3 in target year', 'MW');

%% 3. Map plot: Capacity density per km² (absolute)
resultsToPlot = regioDataWind.capPerKm2;
resultsToPlot(resultsToPlot <= 0.0001) = NaN;

plotCapacityMap(regioDataWind, resultsToPlot, ...
    'Installed capacity per area and NUTS 3 in target year', 'MW/km²');

%% 4. Map plot: Capacity density per available area
resultsToPlot = regioDataWind.capPerKm2avail;
totalCapacityGW = sum(resultsToPlot) / 1e6;
resultsToPlot(resultsToPlot <= 0.0001) = NaN;

% disp(['Total Capacity: ', num2str(totalCapacityGW), ' GW']);

plotCapacityMap(regioDataWind, resultsToPlot, ...
    'Installed capacity per area available and NUTS 3 in target year', 'MW/km²');

end

%% Helper Function: Plot regional map
function plotCapacityMap(regioData, values, titleStr, unitStr)
    maxVal = max(values, [], 'omitnan');
    minVal = min(values, [], 'omitnan');
    rescaled = rescale(values);  % Normalize for transparency

    figure(); hold on;
    plot(polyshape(), 'FaceColor', 'blue', 'FaceAlpha', 1);  % Dummy legend
    plot(polyshape(), 'FaceColor', 'white');
    plot(polyshape(), 'FaceColor', 'black');

    warning('off', 'all');
    for i = 1:height(regioData)
        shape = polyshape(regioData.lonPoly{i}, regioData.latPoly{i});
        if ~isnan(rescaled(i))
            plot(shape, 'FaceColor', 'blue', 'FaceAlpha', rescaled(i));
        else
            plot(shape, 'FaceColor', 'black');
        end
    end
    warning('on', 'all');

    title(titleStr);
    legend(...
        ['max = ', num2str(round(maxVal)), ' ', unitStr], ...
        ['min = ', num2str(round(minVal)), ' ', unitStr], ...
        'Location', 'best');

    xlim([min(regioData.lonPoint)-1, max(regioData.lonPoint)+1]);
    ylim([min(regioData.latPoint)-1, max(regioData.latPoint)+1]);
    axis off;
    hold off;
end
