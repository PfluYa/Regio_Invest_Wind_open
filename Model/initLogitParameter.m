function initParameter = initLogitParameter(logitSettings)
    % INITLOGITPARAMETER Initializes starting parameters for logit estimation
    %
    % This function sets up random initial values for a discrete choice model
    % (nested logit) commonly used in renewable energy modeling, particularly
    % onshore wind investment with 8 turbine alternatives.
    %
    % INPUT:
    %   logitSettings - Struct containing:
    %       .nrAlternatives        - Number of choice alternatives (e.g., turbine types)
    %       .explanatoryVariables  - Cell array of explanatory variable names
    %
    % OUTPUT:
    %   initParameter  - Struct with fields:
    %       .parameter        - Initial values for beta, gamma, lambda, alpha
    %       .parameterNames   - Corresponding parameter names
    %
    % Notes:
    % - Alpha0 (the reference alternative) is fixed (e.g. turbine 4)
    % - Only alpha for non-reference turbines is estimated
    
    % --- Settings ---
    alphaEstimation = 'none';  % Options: 'none' or 'estimateAlpha0'
    
    % --- Read settings ---
    NAlt = logitSettings.nrAlternatives;
    NExpl = numel(logitSettings.explanatoryVariables);
    
    % --- Initialize parameters ---
    if strcmp(alphaEstimation, 'estimateAlpha0')
        alpha = rand(1, NAlt);  % Estimate alpha for all
    else
        alpha = rand(1, NAlt-1);  % Fix one alternative (e.g., turbine 4)
    end
    
    beta = rand(1, NExpl);
    beta(1) = beta(1) * 10;  % Scale first beta parameter
    
    gamma = -1 * rand(1,1);   % Negative nesting parameter
    lambda = rand(1,1);       % Nest elasticity
    
    % --- Store in struct ---
    initParameter = logitSettings;
    initParameter.parameter = [beta, gamma, lambda, alpha];
    
    % --- Assign parameter names ---
    paramNames = cell(1, numel(initParameter.parameter));
    
    for i = 1:NExpl
        paramNames{i} = ['beta', num2str(i)];
    end
    paramNames{NExpl + 1} = 'gamma';
    paramNames{NExpl + 2} = 'lambda';
    
    if strcmp(alphaEstimation, 'estimateAlpha0')
        for i = 1:NAlt
            paramNames{NExpl + 2 + i} = ['alpha', num2str(i)];
        end
    else
        % Reference: Turbine 4 is base, not estimated
        alphaIDs = [1, 2, 3, 5, 6, 7, 8];
        for j = 1:numel(alpha)
            paramNames{NExpl + 2 + j} = ['alpha', num2str(alphaIDs(j))];
        end
    end
    
    initParameter.parameterNames = paramNames(:);
    
    % --- Display output ---
    disp(' ');
    disp('INIT PARAMETER');
    disp('-------------------------');
    for i = 1:numel(paramNames)
        fprintf('%-10s %10.4f\n', paramNames{i}, initParameter.parameter(i));
    end

end
