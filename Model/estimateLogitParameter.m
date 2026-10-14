function discreteChoiceData = estimateLogitParameter(inputData, discreteChoiceSettings)
% ESTIMATELOGITPARAMETER
% Estimate nested-logit parameters by maximum likelihood.
%
% DESCRIPTION
% This function estimates the parameters of the discrete-choice model used
% in the wind investment framework. It prepares the estimation inputs,
% calls MATLAB's unconstrained optimizer, and returns parameter estimates,
% fit measures, and standard errors based on a central finite-difference Hessian evaluated at the final optimum.
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
%       .inferenceHessian
%       .inferenceHessianRelativeStep
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
discreteChoiceData.optimizerHessian = hessian;

discreteChoiceData.nrDecisionMaker         = NDecisionMaker;
discreteChoiceData.nrExplanatoryVariables  = NExpl;
discreteChoiceData.nrAlternatives          = NAlt;

%% ------------------------------------------------------------------------
% 4. INFERENCE FROM A CENTRAL FINITE-DIFFERENCE HESSIAN
% -------------------------------------------------------------------------
% The optimizer's numerical Hessian can be sensitive to internal finite-
% difference choices in this model. For reproducible inference, the
% Hessian is therefore recomputed explicitly at the final optimum using a
% central finite-difference scheme and a documented relative step size.
%
% Parameter-specific steps are:
%   h_i = relativeStep * max(1, abs(theta_i))
%
% A sensitivity check showed stable standard errors around this step size.

relativeStep = ...
    3e-5;

inferenceHessian = ...
    finiteDifferenceHessianCentral( ...
        @negloglike, ...
        paramhat, ...
        relativeStep);

inferenceHessian = ...
    (inferenceHessian + inferenceHessian') ./ 2;

inferenceEigenvalues = ...
    eig(inferenceHessian);

minimumEigenvalue = ...
    min(inferenceEigenvalues);

if minimumEigenvalue <= 0

    error([ ...
        'Finite-difference inference Hessian is not positive definite. ', ...
        'Minimum eigenvalue: %.6e.'], ...
        minimumEigenvalue);
end

covarianceMatrix = ...
    inferenceHessian \ ...
    eye(size(inferenceHessian));

covarianceMatrix = ...
    (covarianceMatrix + covarianceMatrix') ./ 2;

varianceDiagonal = ...
    diag(covarianceMatrix);

if any(~isfinite(varianceDiagonal)) || ...
        any(varianceDiagonal <= 0)

    error( ...
        'Invalid variance estimate from finite-difference Hessian.');
end

stderr = ...
    sqrt(varianceDiagonal);

discreteChoiceData.stderr = ...
    stderr;

discreteChoiceData.inferenceHessian = ...
    inferenceHessian;

discreteChoiceData.inferenceHessianRelativeStep = ...
    relativeStep;

discreteChoiceData.inferenceHessianMinEigenvalue = ...
    minimumEigenvalue;

discreteChoiceData.inferenceHessianConditionNumber = ...
    cond(inferenceHessian);


%% ------------------------------------------------------------------------
% 5. REPORT ESTIMATION RESULTS
% -------------------------------------------------------------------------

disp(['McFadden R^2: ', num2str(discreteChoiceData.RSquare)])
disp(['Negative log-likelihood (estimated): ', num2str(fval)])
disp(['Negative log-likelihood (null model): ', num2str(discreteChoiceData.nll_0)])
disp(['Inference Hessian relative step: ', num2str(relativeStep)])
disp(['Inference Hessian minimum eigenvalue: ', num2str(minimumEigenvalue)])
disp(['Inference Hessian condition number: ', ...
      num2str(discreteChoiceData.inferenceHessianConditionNumber)])
disp(['Gradient'' * inv(Hessian) * Gradient: ', ...
      num2str(grad' * covarianceMatrix * grad)])

disp(' ')
disp('ESTIMATION RESULTS')
disp('----------------------------')
disp('Parameter        Est         SE        t-stat')

for i = 1:numel(paramhat)

    fprintf( ...
        '%-10s %10.4f %10.4f %10.4f\n', ...
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

function hessianMatrix = finiteDifferenceHessianCentral( ...
    objectiveFunction, ...
    parameterVector, ...
    relativeStep)
% FINITEDIFFERENCEHESSIANCENTRAL
% Central finite-difference Hessian evaluated at parameterVector.
%
% Parameter-specific steps are scaled by max(1, abs(theta_i)).

    theta = ...
        parameterVector(:);

    nParameter = ...
        numel(theta);

    step = ...
        relativeStep .* ...
        max(1, abs(theta));

    hessianMatrix = ...
        zeros(nParameter);

    f0 = ...
        objectiveFunction(theta');

    %% Diagonal terms

    for i = 1:nParameter

        thetaPlus = ...
            theta;

        thetaMinus = ...
            theta;

        thetaPlus(i) = ...
            thetaPlus(i) + step(i);

        thetaMinus(i) = ...
            thetaMinus(i) - step(i);

        fPlus = ...
            objectiveFunction(thetaPlus');

        fMinus = ...
            objectiveFunction(thetaMinus');

        hessianMatrix(i, i) = ...
            ( ...
                fPlus - ...
                2 .* f0 + ...
                fMinus ...
            ) ./ ...
            (step(i) .^ 2);
    end

    %% Off-diagonal terms

    for i = 1:nParameter

        for j = (i + 1):nParameter

            thetaPP = theta;
            thetaPM = theta;
            thetaMP = theta;
            thetaMM = theta;

            thetaPP([i, j]) = ...
                thetaPP([i, j]) + ...
                [step(i); step(j)];

            thetaPM(i) = ...
                thetaPM(i) + step(i);

            thetaPM(j) = ...
                thetaPM(j) - step(j);

            thetaMP(i) = ...
                thetaMP(i) - step(i);

            thetaMP(j) = ...
                thetaMP(j) + step(j);

            thetaMM([i, j]) = ...
                thetaMM([i, j]) - ...
                [step(i); step(j)];

            mixedDerivative = ...
                ( ...
                    objectiveFunction(thetaPP') - ...
                    objectiveFunction(thetaPM') - ...
                    objectiveFunction(thetaMP') + ...
                    objectiveFunction(thetaMM') ...
                ) ./ ...
                (4 .* step(i) .* step(j));

            hessianMatrix(i, j) = ...
                mixedDerivative;

            hessianMatrix(j, i) = ...
                mixedDerivative;
        end
    end

    hessianMatrix = ...
        (hessianMatrix + hessianMatrix') ./ 2;
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