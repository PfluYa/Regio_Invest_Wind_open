function initParameter = initLogitParameter(logitSettings)
% INITLOGITPARAMETER
% Create starting values for nested-logit parameter estimation.
%
% DESCRIPTION
% This function initializes the parameter vector used in the discrete-choice
% estimation of the wind investment model. The parameter vector consists of:
%
%   - beta   : coefficients on explanatory variables
%   - gamma  : upper-level nest parameter
%   - lambda : dissimilarity parameter
%   - alpha  : alternative-specific constants for non-reference alternatives
%
% INPUT
%   logitSettings : struct with fields
%       .nrAlternatives       : number of alternatives (current default: 8)
%       .explanatoryVariables : cell array of explanatory variable names
%
% OUTPUT
%   initParameter : struct extending logitSettings with fields
%       .parameter       : initial parameter vector
%       .parameterNames  : parameter names corresponding to .parameter
%
% MODEL CONVENTION
% The current specification fixes one alternative-specific constant as the
% reference category. In the present application, turbine type 4 serves as
% the reference alternative, so alpha4 is not estimated.

%% ------------------------------------------------------------------------
% 1. SETTINGS
% -------------------------------------------------------------------------

% Current default:
%   'none' -> do not estimate the reference-category alpha
%   'estimateAlpha0' -> estimate alpha for all alternatives
alphaEstimation = 'none';

%% ------------------------------------------------------------------------
% 2. READ MODEL DIMENSIONS
% -------------------------------------------------------------------------

nAlt = logitSettings.nrAlternatives;
nExpl = numel(logitSettings.explanatoryVariables);

%% ------------------------------------------------------------------------
% 3. INITIALIZE START VALUES
% -------------------------------------------------------------------------
% Current implementation uses random starting values.
% For reproducible estimation across runs, consider setting rng(...) in the
% main estimation workflow.

if strcmp(alphaEstimation, 'estimateAlpha0')
    alpha = rand(1, nAlt);
else
    alpha = rand(1, nAlt - 1);
end

beta = rand(1, nExpl);

% Current convention:
% scale the first beta to allow a wider initial range
if ~isempty(beta)
    beta(1) = beta(1) * 10;
end

% Current initialization:
% gamma starts negative, lambda starts between 0 and 1
gamma = -rand(1, 1);
lambda = rand(1, 1);

%% ------------------------------------------------------------------------
% 4. STORE INITIAL PARAMETER VECTOR
% -------------------------------------------------------------------------

initParameter = logitSettings;
initParameter.parameter = [beta, gamma, lambda, alpha];

%% ------------------------------------------------------------------------
% 5. ASSIGN PARAMETER NAMES
% -------------------------------------------------------------------------

paramNames = cell(1, numel(initParameter.parameter));

for i = 1:nExpl
    paramNames{i} = ['beta', num2str(i)];
end

paramNames{nExpl + 1} = 'gamma';
paramNames{nExpl + 2} = 'lambda';

if strcmp(alphaEstimation, 'estimateAlpha0')
    for i = 1:nAlt
        paramNames{nExpl + 2 + i} = ['alpha', num2str(i)];
    end
else
    % Current model convention:
    % turbine type 4 is the omitted reference category
    alphaIDs = [1, 2, 3, 5, 6, 7, 8];

    for j = 1:numel(alpha)
        paramNames{nExpl + 2 + j} = ['alpha', num2str(alphaIDs(j))];
    end
end

initParameter.parameterNames = paramNames(:);

%% ------------------------------------------------------------------------
% 6. DISPLAY INITIAL VALUES
% -------------------------------------------------------------------------

disp(' ');
disp('INIT PARAMETER');
disp('-------------------------');

for i = 1:numel(paramNames)
    fprintf('%-10s %10.4f\n', paramNames{i}, initParameter.parameter(i));
end

end