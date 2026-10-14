%% TABLE 2: NESTED-LOGIT PARAMETER ESTIMATES
% Writes estimation results from the 2040 Case-1 run.
%
% Outputs:
%   Results Paper/Table_2_nested_logit_estimates.csv
%   Results Paper/Table_2_nested_logit_estimates.mat

scriptPath = mfilename('fullpath');
scriptFolder = fileparts(scriptPath);
repoFolder = fileparts(scriptFolder);
resultsFolder = fullfile(repoFolder, 'Results Paper');

inputFile = fullfile( ...
    resultsFolder, ...
    'results_nuts3_base2024_target2040_expCase1_capacity.mat');

if ~isfile(inputFile)
    error('2040 Case-1 result file not found: %s', inputFile);
end

S = load(inputFile, 'scenarioSummary');

if ~isfield(S, 'scenarioSummary') || ...
        ~isfield(S.scenarioSummary, 'discreteChoiceEstimation')
    error([ ...
        'The result file does not contain discrete-choice estimation ', ...
        'results. Re-run the 2040 Case-1 scenario first.']);
end

E = S.scenarioSummary.discreteChoiceEstimation;

parameterNames = string(E.parameterNames(:));
estimates = E.estimates(:);
standardErrors = E.standardErrors(:);

if numel(parameterNames) ~= numel(estimates) || ...
        numel(estimates) ~= numel(standardErrors)
    error('Inconsistent dimensions in stored estimation results.');
end

table2Results = table( ...
    parameterNames, ...
    estimates, ...
    standardErrors, ...
    'VariableNames', {'Parameter', 'Estimate', 'StdError'});

table2Results = [ ...
    table2Results; ...
    table( ...
        "Pseudo-R2", ...
        E.pseudoR2, ...
        NaN, ...
        'VariableNames', table2Results.Properties.VariableNames)];

nrDecisionMaker = E.nrDecisionMaker;

fprintf('\n');
fprintf('============================================================\n');
fprintf('TABLE 2: NESTED-LOGIT PARAMETER ESTIMATES\n');
fprintf('============================================================\n\n');
disp(table2Results);
fprintf('Number of NUTS3 decision makers: %d\n', nrDecisionMaker);

csvFile = fullfile( ...
    resultsFolder, ...
    'Table_2_nested_logit_estimates.csv');

matFile = fullfile( ...
    resultsFolder, ...
    'Table_2_nested_logit_estimates.mat');

writetable(table2Results, csvFile);

save( ...
    matFile, ...
    'table2Results', ...
    'nrDecisionMaker');

fprintf('Table 2 CSV saved to:\n%s\n', csvFile);
fprintf('Table 2 MAT file saved to:\n%s\n', matFile);
