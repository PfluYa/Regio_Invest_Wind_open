function fig = plot_results(resultsWindOnshoreTurbinesUnstack, ~, outputFile)
% PLOT_RESULTS
% Plot installed onshore wind capacity by turbine type.
%
% The open-source version is intentionally limited to the output used in
% Figure 3 of the paper. Generic regional maps that are not used in the
% final paper have been removed.
%
% INPUTS
%   resultsWindOnshoreTurbinesUnstack : wide table with NUTS3 regions in
%                                       rows and turbine types in columns
%   second input                       : retained for backward compatibility
%   outputFile                         : optional PDF output path
%
% OUTPUT
%   fig                                : figure handle

if nargin < 3
    outputFile = "";
end

turbineColumns = ...
    resultsWindOnshoreTurbinesUnstack.Properties.VariableNames(2:end);

turbineCapacityGW = ...
    sum( ...
        table2array( ...
            resultsWindOnshoreTurbinesUnstack(:, turbineColumns)), ...
        1) ./ ...
    1e6;

fig = ...
    figure( ...
        'Position', [200 200 850 520], ...
        'Color', 'w');

bar(turbineCapacityGW);

nTurbines = ...
    numel(turbineColumns);

labels = ...
    compose('Type %d', 1:nTurbines);

set(gca, ...
    'XTick', 1:nTurbines, ...
    'XTickLabel', labels, ...
    'FontSize', 11, ...
    'LineWidth', 0.8);

ylabel('Installed capacity (GW)');
xlabel('Turbine type');

grid on;
box on;

if strlength(string(outputFile)) > 0

    exportgraphics( ...
        fig, ...
        outputFile, ...
        'ContentType', 'vector');
end

end
