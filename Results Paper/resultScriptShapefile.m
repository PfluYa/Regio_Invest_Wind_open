% Define the names of the result files
resultsFile_2030 = '2025-08-12_results_nuts3_base2024_target2030_expCase1.mat';
resultsFile_2035 = '2025-08-12_results_nuts3_base2024_target2035_expCase1.mat';
resultsFile_2040 = '2025-08-12_results_nuts3_base2024_target2040_expCase4.mat';

%% Load results for capacities per km^2
% Load data for 2030
load(resultsFile_2030);
results_2030 = selected_data;
results_2030.Properties.VariableNames{'capPerKm2'} = 'capPerKm2_2030';

% Load data for 2035
load(resultsFile_2035);
results_2035 = selected_data;
results_2035.Properties.VariableNames{'capPerKm2'} = 'capPerKm2_2035';

% Load data for 2040
load(resultsFile_2040);
results_2040 = selected_data;
results_2040.Properties.VariableNames{'capPerKm2'} = 'capPerKm2_2040';

% Combine data from different years
results_capPerKm2 = results_2030(:, 1:13);
% in case of expanding cases 2-4, only 2040 is relevant
% results_capPerKm2 = results_2040;
results_capPerKm2.capPerKm2_2035 = results_2035.capPerKm2_2035;
results_capPerKm2.capPerKm2_2040 = results_2040.capPerKm2_2040;
results_capPerKm2.capPerKm2_baseYear = results_2040.capacity_baseYear ./ results_2040.totalArea;

%% Load results for exhaustion probability
% Load data for 2030
load(resultsFile_2030);
results_2030 = selected_data;
results_2030.Properties.VariableNames{'exhaustionProb'} = 'exhaustionProb_2030';

% Load data for 2035
load(resultsFile_2035);
results_2035 = selected_data;
results_2035.Properties.VariableNames{'exhaustionProb'} = 'exhaustionProb_2035';

% Load data for 2040
load(resultsFile_2040);
results_2040.exhaustionProb_base = ((results_2040.capacity_baseYear ./ 1000) ./ 22.5) ./ (results_2040.relativeAvailableWindSpace .* results_2040.totalArea);
results_2040.Properties.VariableNames{'exhaustionProb'} = 'exhaustionProb_2040';

% Combine data for exhaustion probability
results_exh_prob = results_2030(:, [1:10, 17]);
results_exh_prob.exh_prob_2035 = results_2035.exhaustionProb_2035;
results_exh_prob.exh_prob_2040 = results_2040.exhaustionProb_2040;
results_exh_prob.exh_prob_base = results_2040.exhaustionProb_base;
results_exh_prob.Properties.VariableNames{'exhaustionProb_2030'} = 'exh_prob_2030';

% Replace NaN and Inf with 1
fieldsToCheck = {'exh_prob_2030', 'exh_prob_2035', 'exh_prob_2040', 'exh_prob_base'};
for field = fieldsToCheck
    results_exh_prob.(field{:})(isnan(results_exh_prob.(field{:}))) = 1;
    results_exh_prob.(field{:})(isinf(results_exh_prob.(field{:}))) = 1;
    results_exh_prob.(field{:})(results_exh_prob.(field{:}) > 1) = 1;
end

%% Create a shapefile from the data
% Assuming your table is named 'results_capPerKm2'
geoData = geoshape([]); % Create an empty geospatial object

% Loop over the rows of the table and add polygons
for i = 1:size(results_capPerKm2, 1)
    % Create the polygon
    disp(i);
    polyShape = geoshape(results_capPerKm2.latPoly(i), results_capPerKm2.lonPoly(i), 'Geometry', 'polygon');
    
    % Assign desired attributes
    polyShape.nutsID = results_capPerKm2.nutsID{i};
    polyShape.countryCode = results_capPerKm2.countryCode{i};
    polyShape.cityName = results_capPerKm2.cityName{i};
    
    % Add the capacity per km² as attributes
    % polyShape.km2cap_2030 = results_capPerKm2.capPerKm2_2030(i);
    % polyShape.km2cap_2035 = results_capPerKm2.capPerKm2_2035(i);
    polyShape.km2cap_2040 = results_capPerKm2.capPerKm2_2040(i);
    % polyShape.km2cap_base = results_capPerKm2.capPerKm2_baseYear(i);
    
    % Append the polygon to the geoData structure
    geoData = [geoData; polyShape];
end

% Save the shapefile
shapeFile = fullfile(cd, 'Results Paper\2025-08-12_shapeResults_cap_per_km2_base24_Case4');
shapewrite(geoData, shapeFile);


%% Create a shapefile for exhaustion probability
geoData_exh_prob = geoshape([]); % Create an empty geospatial object

% Loop over the rows of the table and add polygons
for i = 1:size(results_exh_prob, 1)
    % Create the polygon
    disp(i);
    polyShape = geoshape(results_exh_prob.latPoly(i), results_exh_prob.lonPoly(i), 'Geometry', 'polygon');
    
    % Assign desired attributes
    polyShape.nutsID = results_exh_prob.nutsID{i};
    polyShape.countryCode = results_exh_prob.countryCode{i};
    polyShape.cityName = results_exh_prob.cityName{i};
    
    % Add the exhaustion probabilities as attributes
    polyShape.exh_prob_2030 = results_exh_prob.exh_prob_2030(i);
    polyShape.exh_prob_2035 = results_exh_prob.exh_prob_2035(i);
    polyShape.exh_prob_2040 = results_exh_prob.exh_prob_2040(i);
    polyShape.exh_prob_base = results_exh_prob.exh_prob_base(i);
    
    % Append the polygon to the geoData structure
    geoData_exh_prob = [geoData_exh_prob; polyShape];
end

% Save the shapefile for exhaustion probability
shapeFile_exh_prob = fullfile(cd, 'Results Paper\shapeResults_exh_prob_2025_01_15');
shapewrite(geoData_exh_prob, shapeFile_exh_prob);
