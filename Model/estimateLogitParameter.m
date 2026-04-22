function discreteChoiceData = estimateLogitParameter(inputData, discreteChoiceSettings)
% ESTIMATELOGITPARAMETER
% Estimate nested-logit parameters by maximum likelihood.
%
% DESCRIPTION
% This function estimates the parameters of the discrete-choice model used
% in the wind investment framework. It prepares the estimation inputs,
% calls MATLAB's unconstrained optimizer, and returns parameter estimates,
% fit measures, and approximate standard errors based on the Hessian.
%
% INPUTS
%   inputData : table containing observed choices and explanatory variables
%
%   discreteChoiceSettings : struct with fields
%       .explainedVariable
%       .explanatoryVariables
%       .parameter
%       .parameterNames
%       .nrDecisionMaker
%       .nrAlternatives
%       .nrExplanatoryVariables
%
% OUTPUT
%   discreteChoiceData : struct containing
%       .discreteChoiceParam
%       .nll_beta
%       .nll_0
%       .RSquare
%       .exitflag
%       .output
%       .grad
%       .hessian
%       .stderr
%       .nrDecisionMaker
%       .nrExplanatoryVariables
%       .nrAlternatives
%
% MODEL CONVENTIONS
% - The current specification fixes turbine type 4 as the reference
%   alternative, i.e. alpha4 = 0.
% - The current implementation uses lambdaDivision = true.
% - No turbine-specific alpha covariate is used in the baseline
%   specification, i.e. z_alpha = 1.

%% ------------------------------------------------------------------------
% 0. GLOBALS USED INSIDE THE LOG-LIKELIHOOD
% -------------------------------------------------------------------------

global y x z_alpha NDecisionMaker NAlt NExpl alphaEstimation lambdaDivision

%% ------------------------------------------------------------------------
% 1. PREPARE ESTIMATION INPUTS
% -------------------------------------------------------------------------

y = inputData{:, discreteChoiceSettings.explainedVariable};
x = inputData{:, discreteChoiceSettings.explanatoryVariables};

% Current baseline: no additional turbine-specific explanatory variable
explanatoryVariableTurbine = 'none';

if strcmp(explanatoryVariableTurbine, 'none')
    z_alpha = ones(size(inputData, 1), 1);
else
    z_alpha = inputData{:, explanatoryVariableTurbine};
end

% Current baseline settings
alphaEstimation = 'none';
lambdaDivision  = true;

NDecisionMaker = discreteChoiceSettings.nrDecisionMaker;
NAlt           = discreteChoiceSettings.nrAlternatives;
NExpl          = discreteChoiceSettings.nrExplanatoryVariables;

%% ------------------------------------------------------------------------
% 2. MAXIMUM-LIKELIHOOD ESTIMATION
% -------------------------------------------------------------------------
% The current implementation uses fminunc with numerical derivatives.

options = optimset( ...
    'LargeScale', 'on', ...
    'Display', 'on', ...
    'GradObj', 'off', ...
    'MaxFunEvals', 1e5, ...
    'MaxIter', [], ...
    'DerivativeCheck', 'on');

[paramhat, fval, exitflag, output, grad, hessian] = ...
    fminunc(@negloglike, discreteChoiceSettings.parameter, options);

%% ------------------------------------------------------------------------
% 3. ASSEMBLE ESTIMATION OUTPUT
% -------------------------------------------------------------------------

discreteChoiceData.discreteChoiceParam = array2table(paramhat);
discreteChoiceData.discreteChoiceParam.Properties.VariableNames = ...
    discreteChoiceSettings.parameterNames';

discreteChoiceData.nll_beta = fval;
discreteChoiceData.nll_0 = negloglike(zeros(size(paramhat)));
discreteChoiceData.RSquare = ...
    1 - discreteChoiceData.nll_beta / discreteChoiceData.nll_0;

discreteChoiceData.exitflag = exitflag;
discreteChoiceData.output   = output;
discreteChoiceData.grad     = grad;
discreteChoiceData.hessian  = hessian;

discreteChoiceData.nrDecisionMaker         = NDecisionMaker;
discreteChoiceData.nrExplanatoryVariables  = NExpl;
discreteChoiceData.nrAlternatives          = NAlt;

%% ------------------------------------------------------------------------
% 4. INFERENCE BASED ON REGULARIZED HESSIAN
% -------------------------------------------------------------------------

disp(['McFadden R^2: ', num2str(discreteChoiceData.RSquare)])
disp(['Negative log-likelihood (estimated): ', num2str(fval)])
disp(['Negative log-likelihood (null model): ', num2str(discreteChoiceData.nll_0)])

hessianRegularized = ensure_positive_definite(hessian);
invHessian = inv(hessianRegularized);
stderr = sqrt(diag(invHessian));

discreteChoiceData.stderr = stderr;

disp(['Gradient'' * inv(Hessian) * Gradient = ', ...
      num2str(grad' * invHessian * grad)])
disp(' ')
disp('ESTIMATION RESULTS')
disp('----------------------------')
disp('Parameter        Est         SE        t-stat')

for i = 1:numel(paramhat)
    fprintf('%-10s %10.4f %10.4f %10.4f\n', ...
        discreteChoiceSettings.parameterNames{i}, ...
        paramhat(i), ...
        stderr(i), ...
        paramhat(i) / stderr(i));
end

disp(' ')

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function hessianRegularized = ensure_positive_definite(hessianMatrix)
% ENSURE_POSITIVE_DEFINITE
% Regularize the Hessian until it is positive definite.
%
% INPUT
%   hessianMatrix : Hessian from the optimizer
%
% OUTPUT
%   hessianRegularized : positive-definite Hessian used for inference
%
% NOTE
% The current implementation applies diagonal ridge regularization with a
% sequence of increasing epsilon values.

    epsilons = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2];

    if all(eig(hessianMatrix) > 0)
        disp('Hessian is already positive definite.');
        hessianRegularized = hessianMatrix;
        return;
    end

    for epsValue = epsilons
        candidate = hessianMatrix + epsValue .* eye(size(hessianMatrix));

        if all(eig(candidate) > 0)
            disp(['Regularized Hessian with epsilon = ', num2str(epsValue)]);
            hessianRegularized = candidate;
            return;
        end
    end

    error('Failed to regularize Hessian to positive definiteness.');
end

function nll = negloglike(paramLogLike)
% NEGLOGLIKE
% Evaluate the negative log-likelihood of the nested-logit model.
%
% PARAMETER VECTOR ORDER
%   [beta, gamma, lambda, alpha]
%
% CURRENT MODEL CONVENTION
%   alpha4 = 0 (turbine type 4 is the reference alternative)

    global y x z_alpha NDecisionMaker NAlt NExpl alphaEstimation lambdaDivision

    %% --------------------------------------------------------------------
    % 1. PREPARE DATA
    % ---------------------------------------------------------------------

    explVar = x;
    explVarAlpha = z_alpha;

    ani = y(:, 1);
    ani = reshape(ani, [NAlt, NDecisionMaker]);

    an0 = 1 - sum(ani, 1);

    %% --------------------------------------------------------------------
    % 2. UNPACK PARAMETERS
    % ---------------------------------------------------------------------

    beta = paramLogLike(1:NExpl);
    gamma = paramLogLike(NExpl + 1);
    lambda = paramLogLike(NExpl + 2);

    if strcmp(alphaEstimation, 'estimateAlpha0')
        alpha = paramLogLike(end - (NAlt - 1):end);
    else
        % Current convention: alpha4 = 0
        alpha = [ ...
            paramLogLike(end - (NAlt - 2):end - (NAlt - 2) + 2), ...
            0, ...
            paramLogLike(end - (NAlt - 2) + 3:end)];
    end

    alpha = repmat(alpha, [1, NDecisionMaker]);

    %% --------------------------------------------------------------------
    % 3. COMPUTE CHOICE PROBABILITIES
    % ---------------------------------------------------------------------

    systematicUtility = explVar * beta' + alpha' .* explVarAlpha;

    if lambda ~= 0 && lambdaDivision
        vni = exp(systematicUtility ./ lambda);
    else
        vni = exp(systematicUtility);
    end

    vni = reshape(vni, [NAlt, NDecisionMaker]);
    vniSum = sum(vni, 1);

    inclusiveValue = log(vniSum);
    vnk = exp(-(gamma + lambda .* inclusiveValue));

    conditionalProb = vni ./ vniSum;
    nestProb = 1 ./ (1 + vnk);
    noInvestmentProb = 1 ./ (1 + exp(gamma + lambda .* inclusiveValue));

    %% --------------------------------------------------------------------
    % 4. LOG-LIKELIHOOD
    % ---------------------------------------------------------------------

    p1 = sum(ani .* log(conditionalProb .* nestProb), 1);
    p2 = an0 .* log(noInvestmentProb);

    nll = -sum(p1 + p2, 2);
end