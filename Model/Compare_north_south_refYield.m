% === Load simulation results (with and without RYM) ====================
withRYM    = load('2025-08-12_results_nuts3_base2024_target2040_expCase1').selected_data(:, [1 12]);
withoutRYM = load('2025-08-12_results_nuts3_base2024_target2040_expCase1_wo_refYield').selected_data(:, [1 12]);

% === 1. Extract NUTS1 codes from NUTS3 =================================
nuts1_with    = cellfun(@(x) x(1:3), withRYM.nutsID,    'UniformOutput', false);
nuts1_without = cellfun(@(x) x(1:3), withoutRYM.nutsID, 'UniformOutput', false);
nuts1_with    = string(nuts1_with);
nuts1_without = string(nuts1_without);

% === 2. Define cluster regions =========================================
southStates   = ["DE1", "DE2"];
neutralStates = ["DEB", "DEC", "DE7", "DEG"];

% === 3. Assign clusters to NUTS1 regions ===============================
cluster_with    = repmat("North", height(withRYM), 1);
cluster_without = repmat("North", height(withoutRYM), 1);

cluster_with(ismember(nuts1_with, southStates))   = "South";
cluster_with(ismember(nuts1_with, neutralStates)) = "Neutral";

cluster_without(ismember(nuts1_without, southStates))   = "South";
cluster_without(ismember(nuts1_without, neutralStates)) = "Neutral";

withRYM.Cluster    = categorical(cluster_with);
withoutRYM.Cluster = categorical(cluster_without);

% === 4. Aggregate total capacity per cluster ===========================
sum_with    = groupsummary(withRYM,    "Cluster", "sum", "capacityTotal");
sum_without = groupsummary(withoutRYM, "Cluster", "sum", "capacityTotal");

sum_with.Properties.VariableNames(end)    = "WithRYM";
sum_without.Properties.VariableNames(end) = "WithoutRYM";

% === 5. Merge results for comparison ===================================
result = outerjoin(sum_with(:, ["Cluster", "WithRYM"]), ...
                   sum_without(:, ["Cluster", "WithoutRYM"]), ...
                   'Keys', 'Cluster', 'MergeKeys', true);

% === 6. Display result =================================================
disp(result)

