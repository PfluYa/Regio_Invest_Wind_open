function [discreteChoiceData] = estimateLogitParameter(inputData, discreteChoiceSettings)
% ESTIMATELOGITPARAMETER Estimates parameters of a discrete choice model via MLE.
%
% INPUTS:
%   inputData              - Table with observed investment decisions and explanatory variables
%   discreteChoiceSettings - Struct with settings:
%       .explainedVariable
%       .explanatoryVariables
%       .parameter
%       .parameterNames
%       .nrDecisionMaker
%       .nrAlternatives
%       .nrExplanatoryVariables
%
% OUTPUT:
%   discreteChoiceData     - Struct containing estimation results (parameters, fit metrics, etc.)

    % Global variables needed for log-likelihood calculation
    global y x z_alpha NDecisionMaker NAlt NExpl alphaEstimation lambdaDivision

    % === Prepare inputs ===
    y = inputData{:, discreteChoiceSettings.explainedVariable};
    x = inputData{:, discreteChoiceSettings.explanatoryVariables};

    % Optional: turbine-specific explanatory variable
    explanatoryVariable_turbine = 'none';
    if strcmp(explanatoryVariable_turbine, 'none')
        z_alpha = ones(size(inputData, 1), 1);
    else
        z_alpha = inputData{:, explanatoryVariable_turbine};
    end

    % Model settings
    alphaEstimation = 'none';
    lambdaDivision  = true;

    NDecisionMaker  = discreteChoiceSettings.nrDecisionMaker;
    NAlt            = discreteChoiceSettings.nrAlternatives;
    NExpl           = discreteChoiceSettings.nrExplanatoryVariables;

    % === Estimation via Maximum Likelihood ===
    options = optimset('LargeScale','on','Display','on','GradObj','off',...
                       'MaxFunEvals',1e5,'MaxIter',[],'DerivativeCheck','on');

    [paramhat,fval,exitflag,output,grad,hessian] = ...
        fminunc(@negloglike, discreteChoiceSettings.parameter, options);

    % === Assemble output ===
    discreteChoiceData.discreteChoiceParam = array2table(paramhat);
    discreteChoiceData.discreteChoiceParam.Properties.VariableNames = ...
        discreteChoiceSettings.parameterNames';

    discreteChoiceData.nll_beta = fval;
    discreteChoiceData.nll_0 = negloglike(zeros(size(paramhat)));
    discreteChoiceData.RSquare = 1 - discreteChoiceData.nll_beta / discreteChoiceData.nll_0;

    discreteChoiceData.exitflag = exitflag;
    discreteChoiceData.output = output;
    discreteChoiceData.grad = grad;
    discreteChoiceData.hessian = hessian;

    discreteChoiceData.nrDecisionMaker = NDecisionMaker;
    discreteChoiceData.nrExplanatoryVariables = NExpl;
    discreteChoiceData.nrAlternatives = NAlt;

    % === Diagnostics and uncertainty ===
    disp(['The McFadden R² is: ', num2str(discreteChoiceData.RSquare)])
    disp(['Negative log-likelihood (estimated): ', num2str(fval)])
    disp(['Negative log-likelihood (null model): ', num2str(discreteChoiceData.nll_0)])

    hessian_reg = ensure_positive_definite(hessian);
    ihess = inv(hessian_reg);
    stderr = sqrt(diag(ihess));

    disp(['Gradient * inv(Hessian) * Gradient = ', num2str(grad' * ihess * grad)])
    disp(' ')
    disp('ESTIMATION RESULTS')
    disp('----------------------------')
    disp('Parameter        Est       SE        t-stat')
    for i = 1:length(paramhat)
        fprintf('%-10s %10.4f %10.4f %10.4f\n', ...
            discreteChoiceSettings.parameterNames{i}, ...
            paramhat(i), stderr(i), paramhat(i) / stderr(i));
    end
    disp(' ');
end


function hessian_regularized = ensure_positive_definite(Hessian)
%ENSURE_POSITIVE_DEFINITE Ensures Hessian matrix is positive definite by regularizing.
%
% INPUT:
%   Hessian - Original Hessian matrix
%
% OUTPUT:
%   hessian_regularized - Regularized Hessian matrix

    epsilons = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2];
    
    for eps = epsilons
        if all(eig(Hessian) > 0)
            disp('Hessian is already positive definite.');
            hessian_regularized = Hessian;
            return;
        else
            Hessian = Hessian + eps * eye(size(Hessian));
            disp(['Regularizing Hessian with epsilon = ', num2str(eps)]);
        end
    end

    if any(eig(Hessian) <= 0)
        error('Failed to regularize Hessian. Not positive definite.');
    end

    hessian_regularized = Hessian;
    disp('Hessian successfully regularized.');
end




function [nll] = negloglike(paramLogLike)
    %% This Script constructs a negative log likelihood function for a given problem
    %  Input parameters are the parameters to estimate in param, the explanatory
    %  variables x and the obsvered dicision y
    %  param is a table with beta, gamma and lambda – gamma and lambda are 1x1
    %  values, the size of beta depends on the amount of explanatory variables,
    %  therefore beta is 1xNumberX vector
    %  x is a matrix n*i x NrExplanatory ( i are the alternatives and n is
    %  the usergroup oder user n -> in our Case NrExplanatory = 1  -> size(x) = n*i x 1)
    %  y is a matrix with the size i*n
    
    global y x z_alpha NDecisionMaker NAlt NExpl alphaEstimation lambdaDivision
    
    % init explanatory variables and the obsveration
    explVar = x;
    explVar_alpha = z_alpha;
    ani = y(:,1);
    ani = reshape(ani,[NAlt,NDecisionMaker]);
    an0 = 1-sum(ani,1);
    %an0 = reshape(an0,[NAlt, NDecisionMaker]);
    %an0 = floor(sum(reshape(y(:,2),[NAlt,NDecisionMaker]),1)/NAlt)-sum(ani,1);
    % init the parameter to estimate
    
    beta = paramLogLike(1:NExpl); %Nx x 1 Vector
    gamma = paramLogLike(NExpl+1);
    lambda = paramLogLike(NExpl+2);
    if strcmp(alphaEstimation,'estimateAlpha0')
        alpha = paramLogLike(end-(NAlt-1):end); %Alpha0 ~=0
    else
        alpha = [paramLogLike(end-(NAlt-2):end-(NAlt-2)+2) 0 paramLogLike(end-(NAlt-2)+3:end)]; %Alpha4 = 0 - Turbine 4 now reference turbine
    end
    alpha = repmat(alpha,[1,NDecisionMaker]);
    %
    %% Equation of the log likelihood in General
    
    if lambda ~= 0
        if lambdaDivision
            vni = exp((explVar*beta' + alpha'.*explVar_alpha)./lambda);
        else
            vni = exp((explVar*beta' + alpha'.*explVar_alpha));
        end
    else
        vni = exp(explVar*beta' + alpha'.*explVar_alpha);
    end
    vni = reshape(vni,[NAlt,NDecisionMaker]);
    vni_sum = sum(vni,1);
    
    wnk = gamma;
    Ink = log(vni_sum);
    vnk = exp(-(wnk+lambda.*Ink));
    
    p1 = sum(ani.*log((vni./vni_sum).*(1./(1+vnk))),1); 
    p2 = an0.*log(1./(1+exp(wnk+lambda.*Ink)));
    p = p1+p2;
    nll = -(sum(p,2));

end