function transformedNPV = transformNPV(NPV)
    %TRANSFORMNPV Applies conditional log transformation to NPV values.
    %   For abs(NPV) >= 1: sgn(NPV) * ln(abs(NPV))
    %   For -1 < NPV < 1: 0
    % 
    % transformedNPV = zeros(size(NPV));  % Initialize
    % 
    % idx = abs(NPV) >= 1;
    % transformedNPV(idx) = sign(NPV(idx)) .* log(abs(NPV(idx)));

    transformedNPV = asinh(NPV);
    % transformedNPV = log(max(NPV,0.00000000000000000000000001));
end
