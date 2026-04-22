%% Compare reference-yield scaling under energy-based vs capacity-based targets
% This script compares sensitivity results from two policy logics:
%   - energy-based target
%   - capacity-based target
%
% Optional:
%   plotLinearized = true  -> also plot linearized reference-yield variants
%   plotLinearized = false -> plot only baseline variants
%
% Styling:
%   - Energy-based:   blue tones
%   - Capacity-based: red tones
%   - Baseline:       solid line
%   - Linearized:     dashed line

clear;
clc;

%% ------------------------------------------------------------------------
% 1. USER SETTINGS
% -------------------------------------------------------------------------

plotLinearized = false;   % default: false

% File names
fileEnergyBased   = "Results_scaling_energy_based.mat";
fileCapacityBased = "Results_scaling_capacity_based.mat";

fileEnergyLinear   = "Results scaling energy based linear.mat";
fileCapacityLinear = "Results scaling capacity based linear.mat";

%% ------------------------------------------------------------------------
% 2. LOAD DATA
% -------------------------------------------------------------------------

S_EN = load(fileEnergyBased);
S_CN = load(fileCapacityBased);

assert(isfield(S_EN, "resultsScaling"), ...
    'File "%s" does not contain variable "resultsScaling".', fileEnergyBased);
assert(isfield(S_CN, "resultsScaling"), ...
    'File "%s" does not contain variable "resultsScaling".', fileCapacityBased);

T_EN = S_EN.resultsScaling;
T_CN = S_CN.resultsScaling;

if plotLinearized
    S_EL = load(fileEnergyLinear);
    S_CL = load(fileCapacityLinear);

    assert(isfield(S_EL, "resultsScaling"), ...
        'File "%s" does not contain variable "resultsScaling".', fileEnergyLinear);
    assert(isfield(S_CL, "resultsScaling"), ...
        'File "%s" does not contain variable "resultsScaling".', fileCapacityLinear);

    T_EL = S_EL.resultsScaling;
    T_CL = S_CL.resultsScaling;
end

%% ------------------------------------------------------------------------
% 3. EXTRACT SERIES
% -------------------------------------------------------------------------

% X-axis
xEN = T_EN.scaleRefYield;
xCN = T_CN.scaleRefYield;

% Welfare
wEN = T_EN.welfarePlusGasCO2_bnEUR;
wCN = T_CN.welfarePlusGasCO2_bnEUR;

% Subsidy expenditure
sEN = T_EN.subsidyPV_bnEUR;
sCN = T_CN.subsidyPV_bnEUR;

% Inframarginal rent per kWh
rEN = T_EN.rent_per_kWh;
rCN = T_CN.rent_per_kWh;

% Total installed capacity
cEN = getTotalCapacityGW(T_EN);
cCN = getTotalCapacityGW(T_CN);

if plotLinearized
    xEL = T_EL.scaleRefYield;
    xCL = T_CL.scaleRefYield;

    wEL = T_EL.welfarePlusGasCO2_bnEUR;
    wCL = T_CL.welfarePlusGasCO2_bnEUR;

    sEL = T_EL.subsidyPV_bnEUR;
    sCL = T_CL.subsidyPV_bnEUR;

    rEL = T_EL.rent_per_kWh;
    rCL = T_CL.rent_per_kWh;

    cEL = getTotalCapacityGW(T_EL);
    cCL = getTotalCapacityGW(T_CL);
end

%% ------------------------------------------------------------------------
% 4. COLORS
% -------------------------------------------------------------------------

blueDark  = [0.00 0.45 0.74];
blueLight = [0.40 0.65 0.90];

redDark   = [0.85 0.20 0.20];
redLight  = [0.95 0.50 0.50];

%% ------------------------------------------------------------------------
% 5. CREATE FIGURE
% -------------------------------------------------------------------------

figure('Position', [100 100 1100 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ------------------------------------------------------------------------
% 5.1 Welfare
% -------------------------------------------------------------------------

nexttile(1);
plot(xEN, wEN, '-', 'Color', blueDark, 'LineWidth', 2); hold on;
plot(xCN, wCN, '-', 'Color', redDark,  'LineWidth', 2);

if plotLinearized
    plot(xEL, wEL, '--', 'Color', blueLight, 'LineWidth', 2);
    plot(xCL, wCL, '--', 'Color', redLight,  'LineWidth', 2);
end

grid on;
xlabel('Scaling factor s');
ylabel('Welfare (bn EUR / year)');
title('Welfare');

if plotLinearized
    legend( ...
        'Energy-based (baseline)', ...
        'Capacity-based (baseline)', ...
        'Energy-based (linearized)', ...
        'Capacity-based (linearized)', ...
        'Location', 'best');
else
    legend('Energy-based', 'Capacity-based', 'Location', 'best');
end

%% ------------------------------------------------------------------------
% 5.2 Installed capacity
% -------------------------------------------------------------------------

nexttile(2);
plot(xEN, cEN, '-', 'Color', blueDark, 'LineWidth', 2); hold on;
plot(xCN, cCN, '-', 'Color', redDark,  'LineWidth', 2);

if plotLinearized
    plot(xEL, cEL, '--', 'Color', blueLight, 'LineWidth', 2);
    plot(xCL, cCL, '--', 'Color', redLight,  'LineWidth', 2);
end

grid on;
xlabel('Scaling factor s');
ylabel('Installed capacity (GW)');
title('Installed capacity');

%% ------------------------------------------------------------------------
% 5.3 Subsidy expenditure
% -------------------------------------------------------------------------

nexttile(3);
plot(xEN, sEN, '-', 'Color', blueDark, 'LineWidth', 2); hold on;
plot(xCN, sCN, '-', 'Color', redDark,  'LineWidth', 2);

if plotLinearized
    plot(xEL, sEL, '--', 'Color', blueLight, 'LineWidth', 2);
    plot(xCL, sCL, '--', 'Color', redLight,  'LineWidth', 2);
end

grid on;
xlabel('Scaling factor s');
ylabel('bn EUR / year');
title('Annual subsidy payments');

%% ------------------------------------------------------------------------
% 5.4 Inframarginal rent per kWh
% -------------------------------------------------------------------------

nexttile(4);
plot(xEN, rEN, '-', 'Color', blueDark, 'LineWidth', 2); hold on;
plot(xCN, rCN, '-', 'Color', redDark,  'LineWidth', 2);

if plotLinearized
    plot(xEL, rEL, '--', 'Color', blueLight, 'LineWidth', 2);
    plot(xCL, rCL, '--', 'Color', redLight,  'LineWidth', 2);
end

grid on;
xlabel('Scaling factor s');
ylabel('EUR / kWh');
title('Inframarginal rent per kWh (annualized)');

%% ------------------------------------------------------------------------
% 6. HARMONIZE AXES
% -------------------------------------------------------------------------

allAxes = findall(gcf, 'type', 'axes');
for iAx = 1:numel(allAxes)
    xlim(allAxes(iAx), [0 1.6]);
end

if plotLinearized
    sgtitle( ...
        'Reference-yield scaling: energy-based vs capacity-based targets (baseline and linearized)', ...
        'FontWeight', 'bold');
else
    sgtitle( ...
        'Reference-yield scaling: energy-based vs capacity-based targets', ...
        'FontWeight', 'bold');
end

%% ========================================================================
% LOCAL HELPER
% ========================================================================

function cap = getTotalCapacityGW(T)
% Determine total installed capacity in GW from the resultsScaling table.

    vars = string(T.Properties.VariableNames);

    if ismember("installedCapacityGW", vars)
        cap = T.installedCapacityGW;
        return;
    end

    requiredVars = ["capNorth_GW", "capNeutral_GW", "capSouth_GW"];
    if all(ismember(requiredVars, vars))
        cap = T.capNorth_GW + T.capNeutral_GW + T.capSouth_GW;
        return;
    end

    error('Could not determine total installed capacity from resultsScaling.');
end