function transformedNPV = transformNPV(npv)
% TRANSFORMNPV
% Transform net present value (NPV) values for use in the discrete-choice
% model.
%
% DESCRIPTION
% The current baseline specification applies the inverse hyperbolic sine
% transformation:
%
%   transformedNPV = asinh(npv)
%
% This transformation is similar to a logarithm for large positive values,
% but remains well-defined for zero and negative values. It is therefore
% suitable for profitability measures such as NPV that may take both
% positive and negative values.
%
% INPUT
%   npv : net present value [EUR]
%
% OUTPUT
%   transformedNPV : transformed NPV used in estimation and simulation [-]
%
% NOTE
% The current model baseline uses the asinh transformation throughout.
% Earlier alternative transformations (e.g. sign-log or log of positive
% values only) are intentionally not used here.

    transformedNPV = asinh(npv);

end