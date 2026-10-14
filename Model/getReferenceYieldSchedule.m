function [relativeYield, correctionFactor] = getReferenceYieldSchedule()
% GETREFERENCEYIELDSCHEDULE
% Return the reference-yield schedule used throughout the model.
%
% OUTPUTS
%   relativeYield     : yield relative to the reference yield [-]
%   correctionFactor  : remuneration correction factor at s = 1 [-]
%
% Counterfactual schedules are constructed as
%
%   corr(s) = 1 + s * (corr(1) - 1)
%
% so that s = 0 implies uniform remuneration, s = 1 reproduces the
% current schedule, and s > 1 strengthens spatial differentiation.

relativeYield = ...
    [0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.50];

correctionFactor = ...
    [1.55 1.42 1.29 1.16 1.07 1.00 0.94 0.89 0.85 0.81 0.79];

end
