%%  step08_subsystem_sweep_mass_balanced

%  Repeats the subsystem sweep analysis from step07 but with mass-balanced
%  perturbations (null-space optimization so S*J = 0 is maintained). This
%  provides more biologically realistic subsystem sensitivity curves. Also
%  generates inter-subsystem correlograms.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (predict, fitlm, corr)
%    - Optimization Toolbox         (optimvar, optimproblem, solve)
%
%  Inputs:
%    - out/iCardio_optimized.mat                     — ModelOpt with S
%                                                      matrix for null-space
%    - out/deltaFlux_stats.mat                       — RxnStatsMedMed
%    - out/ROR_results.mat                           — regresTable
%    - data/drugMap.xlsx                             — Drug mapping table
%    - out/DICT_L0X/c1_<date>/Rxn_models_<a>.mat    — 400 classifiers
%
%  Outputs:
%    - out/driver8_results.mat       — Mass-balanced subsystem correlations,
%                                      sweep scores
%    - figures/Fig8_SubsystemPositiveCorr.png/.fig
%    - figures/Fig8_SubsystemNegativeCorr.png/.fig
%    - figures/Fig8_TopPosRegression.png/.fig
%    - figures/Fig8_TopNegRegression.png/.fig
%    - figures/Fig8_FavRegression.png/.fig
%    - figures/Fig5_FavSweepScores.png/.fig
%    - figures/AllSubsysFullSweeps.png/.fig
%    - figures/Correlogram_Sweep_<subsys>.png/.fig

clear; clc; close all;

%% Parameters
expName = "DICT_L0X";
currentDateString = "10_6_30hr";
iRange = 1:100;
cMax = 4;
saveTF = false;
sweepRange = [-0.1, 0, 0.1];
topN = 10; % Number of top subsystems to display

%% Load Model and Prepare Null Space
fprintf('Loading optimized model and computing null space...\n');
load(fullfile('out','iCardio_optimized.mat'), 'ModelOpt');
Model = ModelOpt;
clear ModelOpt;

% Compute null space if not already present
if ~isfield(Model, 'N')
    fprintf('Computing null space of S...\n');
    Model.N = null(full(Model.S));
    Model.NNp = Model.N * Model.N';
    fprintf('Null space computed. Rank(S) = %d, Rank(N) = %d\n', ...
        rank(full(Model.S)), size(Model.N, 2));
else
    fprintf('Using pre-computed null space\n');
end

fprintf('Model loaded: %d reactions, %d metabolites\n', ...
    length(Model.rxns), length(Model.mets));

%% Load Feature Names and Trained Models
fprintf('Loading feature matrix and trained classifiers...\n');
load(fullfile('out','deltaFlux_stats.mat'), 'RxnStatsMedMed');
load(fullfile('out','ROR_results.mat'), 'regresTable');
drugMap = readtable(fullfile("data", "drugMap.xlsx"));

predictorNames = "x"+string(int2str([1:length(Model.rxns)]'));
uniqueSubSystems = unique(string(Model.subSystems));
uniqueSubSystems(uniqueSubSystems == "") = [];

%% Load All Models for Correlation Analysis
fprintf('\nLoading all trained models...\n');
basePath = fullfile("out", expName, "c1_"+ currentDateString);
numModels = 400;
allModels = cell(numModels, 1);
modelCounter = 1;

for idx = 1:length(iRange)
    a = iRange(idx);
    modelFile = fullfile(basePath, sprintf('Rxn_models_%d.mat', a));
    if ~exist(modelFile, 'file'), continue; end

    load(modelFile, 'tRxn_models');
    for c = 1:cMax
        if modelCounter <= numModels
            allModels{modelCounter} = tRxn_models{1,c};
            modelCounter = modelCounter + 1;
        end
    end
end

allModels = allModels(~cellfun(@isempty, allModels));
fprintf('Loaded %d models for analysis\n', length(allModels));

%% MASS-BALANCED SUBSYSTEM CORRELATION SWEEP
fprintf('\n=== PERFORMING MASS-BALANCED SUBSYSTEM SWEEP ===\n');
fprintf('This may take a while due to optimization at each sweep point...\n');

subsystemCorr = table();
subsystemCorr.SubSystem = uniqueSubSystems;
subsystemCorr.MeanCorrelation = zeros(length(uniqueSubSystems), 1);
subsystemCorr.SEMCorrelation = zeros(length(uniqueSubSystems), 1);
subsystemCorr.NumReactions = zeros(length(uniqueSubSystems), 1);

% Initialize storage for all sweep scores
numSubsys = length(uniqueSubSystems);
allSweepScores = nan(numSubsys, length(allModels), length(sweepRange));

% Pre-compute baseline flux
J0 = Model.NNp * ones(length(Model.rxns), 1);
fprintf('Baseline flux J0 computed (min=%.2e, max=%.2e)\n', min(J0), max(J0));

% Precompute subsystem reaction indices for averaging
subsIdxList = cell(numSubsys, 1);
for j = 1:numSubsys
    subsysJ = uniqueSubSystems(j);
    inSubsysJ = string(Model.subSystems) == subsysJ;
    subsIdxList{j} = find(inSubsysJ);
end

   allSubsysAverages = cell(numSubsys, 1);

% Loop over subsystems
for i = 1:length(uniqueSubSystems)
    subsys = uniqueSubSystems(i);

    % Find reactions in this subsystem
    inSubsys = string(Model.subSystems) == subsys;
    subsysIdx = find(inSubsys);
    subsystemCorr.NumReactions(i) = length(subsysIdx);

    if isempty(subsysIdx)
        fprintf(' [%d/%d] %s: EMPTY - skipping\n', i, numSubsys, subsys);
        continue;
    end

    fprintf(' [%d/%d] %s (%d rxns): ', i, numSubsys, subsys, length(subsysIdx));

    % Create mass-balanced sweep inputs for this subsystem
    inputMatrix_MB = create_mass_balanced_sweep_local(Model, subsysIdx, sweepRange, J0);

    % Compute average % changes per subsystem during this sweep (for inter-subsystem analysis)
    subsys_changes = zeros(length(sweepRange), numSubsys);
    for s = 1:length(sweepRange)
        for j = 1:numSubsys
            rxnIdx = subsIdxList{j};
            if ~isempty(rxnIdx)
                avg_change = mean(inputMatrix_MB(s, rxnIdx));
                % Handle NaN/Inf (replace with 0, as in create_mass_balanced_sweep_local)
                if isnan(avg_change) || isinf(avg_change)
                    avg_change = 0;
                end
                subsys_changes(s, j) = avg_change;
            end
        end
    end
    allSubsysAverages{i} = subsys_changes;



    fprintf('sweep done, predicting...');

    % Store correlations for each model
    corrPerModel = nan(length(allModels), 1);

    for m = 1:length(allModels)
        mdl = allModels{m};

        % Predict cardiotoxicity scores for all sweep points
        [~, scores] = predict(mdl, inputMatrix_MB);
        sweepScores = scores(:, 2); % Class 1 (toxic) probability

        % Store full sweep scores
        allSweepScores(i, m, :) = sweepScores;

        % Calculate Spearman correlation between sweep % and cardiotoxicity
        validScores = ~isnan(sweepScores);
        if sum(validScores) > 2
            corrPerModel(m) = corr(sweepRange(validScores)', ...
                sweepScores(validScores), 'Type', 'Spearman');
        end
    end

    % Aggregate across models
    subsystemCorr.MeanCorrelation(i) = mean(corrPerModel, 'omitnan');
    subsystemCorr.SEMCorrelation(i) = std(corrPerModel, 'omitnan') / ...
        sqrt(sum(~isnan(corrPerModel)));

    fprintf(' r=%.3f ± %.3f\n', subsystemCorr.MeanCorrelation(i), ...
        subsystemCorr.SEMCorrelation(i));
end

fprintf('\nMass-balanced correlation sweep complete!\n');

%% Sort Subsystems by Correlation
posCorr = subsystemCorr(subsystemCorr.MeanCorrelation > 0, :);
posCorr = sortrows(posCorr, 'MeanCorrelation', 'descend');

negCorr = subsystemCorr(subsystemCorr.MeanCorrelation < 0, :);
negCorr = sortrows(negCorr, 'MeanCorrelation', 'ascend');

%% Create Figures Directory
if ~exist('figures', 'dir'), mkdir('figures'); end

%% FIGURE 1: Subsystems Positively Correlated with Cardiotoxicity
fprintf('\nGenerating Figure 1: Positive correlations...\n');
topPos = posCorr(1:min(topN, height(posCorr)), :);

fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
hold on;

yPos = 1:height(topPos);
barh(yPos, topPos.MeanCorrelation, 'FaceColor', [0.8 0.3 0.3], ...
    'EdgeColor', 'k', 'LineWidth', 1);

% Add error bars
errorbar(topPos.MeanCorrelation, yPos, topPos.SEMCorrelation, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 6);

% Format axes
set(gca, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames_local(topPos.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('Mean Correlation with Cardiotoxicity Score', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystems Positively Correlated with Cardiotoxicity (Mass-Balanced)', ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(topPos.MeanCorrelation + topPos.SEMCorrelation) * 1.1]);
ylim([0.5 height(topPos) + 0.5]);
grid on; box on;

% Add reaction counts
for i = 1:height(topPos)
    numRxns = topPos.NumReactions(i);
    text(0.02 * max(xlim), i, sprintf('n=%d', numRxns), ...
        'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

if saveTF
    saveas(fig1, fullfile('figures', 'Fig8_SubsystemPositiveCorr.png'));
    savefig(fig1, fullfile('figures', 'Fig8_SubsystemPositiveCorr.fig'));
    fprintf('Saved Figure 1\n');
end

%% FIGURE 2: Subsystems Negatively Correlated with Cardiotoxicity
fprintf('\nGenerating Figure 2: Negative correlations...\n');
topNeg = negCorr(1:min(topN, height(negCorr)), :);

fig2 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
hold on;

yPos = 1:height(topNeg);
barh(yPos, abs(topNeg.MeanCorrelation), 'FaceColor', [0.3 0.7 0.3], ...
    'EdgeColor', 'k', 'LineWidth', 1);

errorbar(abs(topNeg.MeanCorrelation), yPos, topNeg.SEMCorrelation, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 6);

set(gca, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames_local(topNeg.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('Mean |Correlation| with Cardiotoxicity Score', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystems Negatively Correlated with Cardiotoxicity (Mass-Balanced)', ...
    'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(abs(topNeg.MeanCorrelation) + topNeg.SEMCorrelation) * 1.1]);
ylim([0.5 height(topNeg) + 0.5]);
grid on; box on;

for i = 1:height(topNeg)
    numRxns = topNeg.NumReactions(i);
    text(0.02 * max(xlim), i, sprintf('n=%d', numRxns), ...
        'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

if saveTF
    saveas(fig2, fullfile('figures', 'Fig8_SubsystemNegativeCorr.png'));
    savefig(fig2, fullfile('figures', 'Fig8_SubsystemNegativeCorr.fig'));
    fprintf('Saved Figure 2\n');
end

%% FIGURE 3: Top Positive Subsystems - Regression Plots
fprintf('\nGenerating Figure 3: Top positive subsystem regressions...\n');

fig3 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 14], 'Color', 'w');
tiledlayout('flow', 'TileSpacing', 'compact');

% Focus on narrower range for visualization
minRangeIdx = find(abs(sweepRange + 0.1) < 1e-6, 1);
maxRangeIdx = find(abs(sweepRange - 0.1) < 1e-6, 1);
if isempty(minRangeIdx), minRangeIdx = 1; end
if isempty(maxRangeIdx), maxRangeIdx = length(sweepRange); end

cutRange = sweepRange(minRangeIdx:maxRangeIdx) .* 100; % Convert to percentage
cutSweepScores = allSweepScores(:, :, minRangeIdx:maxRangeIdx);

for a = 1:min(topN, height(topPos))
    nexttile;
    topPosSubsys = topPos.SubSystem(a);
    topPosIdx = find(uniqueSubSystems == topPosSubsys);

    % Extract and process data
    topPosData = squeeze(cutSweepScores(topPosIdx, :, :)); % models × sweeps
    processedData = process_sweep_data_local(topPosData, cutRange);

    % Compute statistics
    meanSweep = mean(processedData, 1, 'omitnan')';
    validPerSweep = sum(~isnan(processedData), 1)';
    semSweep = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

    % Linear fit
    p_mdl = polyfit(cutRange, meanSweep, 1);
    yFit = polyval(p_mdl, cutRange);

    % Plot
    plot(cutRange, meanSweep, 'b-', 'LineWidth', 2);
    hold on;
    fill([cutRange fliplr(cutRange)], ...
        [meanSweep + semSweep; flipud(meanSweep - semSweep)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(cutRange, yFit, 'r--', 'LineWidth', 2);

    xlabel('Sweep % Change', 'FontSize', 10);
    ylabel('Norm. ΔScore', 'FontSize', 10);
    title(sprintf('%s\nr=%.3f', cleanSubsystemNames_local(topPosSubsys), ...
        topPos.MeanCorrelation(a)), 'FontSize', 9, 'Interpreter', 'none');
    grid on; box on;
    ylim([-0.3 1]);
end

if saveTF
    saveas(fig3, fullfile('figures', 'Fig8_TopPosRegression.png'));
    savefig(fig3, fullfile('figures', 'Fig8_TopPosRegression.fig'));
    fprintf('Saved Figure 3\n');
end

%% FIGURE 4: Top Negative Subsystems - Regression Plots
if height(topNeg) > 0
    fprintf('\nGenerating Figure 4: Top negative subsystem regressions...\n');

    fig4 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 14], 'Color', 'w');
    tiledlayout('flow', 'TileSpacing', 'compact');

    for b = 1:min(topN, height(topNeg))
        nexttile;
        topNegSubsys = topNeg.SubSystem(b);
        topNegIdx = find(uniqueSubSystems == topNegSubsys);

        % Extract and process data
        topNegData = squeeze(cutSweepScores(topNegIdx, :, :));
        processedData = process_sweep_data_local(topNegData, cutRange);

        % Compute statistics
        meanSweep = mean(processedData, 1, 'omitnan')';
        validPerSweep = sum(~isnan(processedData), 1)';
        semSweep = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

        % Linear fit
        p_mdl = polyfit(cutRange, meanSweep, 1);
        yFit = polyval(p_mdl, cutRange);

        % Plot
        plot(cutRange, meanSweep, 'g-', 'LineWidth', 2);
        hold on;
        fill([cutRange fliplr(cutRange)], ...
            [meanSweep + semSweep; flipud(meanSweep - semSweep)], ...
            'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(cutRange, yFit, 'r--', 'LineWidth', 2);

        xlabel('Sweep % Change', 'FontSize', 10);
        ylabel('Norm. ΔScore', 'FontSize', 10);
        title(sprintf('%s\nr=%.3f', cleanSubsystemNames_local(topNegSubsys), ...
            topNeg.MeanCorrelation(b)), 'FontSize', 9, 'Interpreter', 'none');
        grid on; box on;
        ylim([-0.3 1]);
    end

    if saveTF
        saveas(fig4, fullfile('figures', 'Fig8_TopNegRegression.png'));
        savefig(fig4, fullfile('figures', 'Fig8_TopNegRegression.fig'));
        fprintf('Saved Figure 4\n');
    end
end

%% FIGURE 5: User-Specified Favorite Subsystems
FavSubSystems = ["Omega-3 fatty acid metabolism";
    "Omega-6 fatty acid metabolism";
    "Inositol phosphate metabolism";
    "Carbohydrate metabolism";
    "Beta oxidation of fatty acids";
    "Nucleotide metabolism";
    "Glycine serine and threonine metabolism";
    "Exchange";
    "Purine metabolism";
    "Central carbon metabolism"];

subsystemCorr.Properties.RowNames = subsystemCorr.SubSystem;
FavCor = subsystemCorr(FavSubSystems, :);

fprintf('\nGenerating Figure 5: Favorite subsystems regression...\n');
fig5 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 14], 'Color', 'w');
tiledlayout(2, 5, 'TileSpacing', 'compact');

for b = 1:height(FavCor)
    nexttile;
    favSubsys = FavCor.SubSystem(b);
    favIdx = find(uniqueSubSystems == favSubsys);

    if isempty(favIdx)
        continue;
    end

    % Extract and process data
    favData = squeeze(cutSweepScores(favIdx, :, :));
    processedData = process_sweep_data_local(favData, cutRange);

    % Compute statistics
    meanSweep = mean(processedData, 1, 'omitnan')';
    validPerSweep = sum(~isnan(processedData), 1)';
    semSweep = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

    % Linear fit with statistics
    p_mdl = fitlm(cutRange, meanSweep, 'linear');
    r_pearson = corr(cutRange', meanSweep);
    yFit = predict(p_mdl, cutRange(:));

    % Plot
    plot(cutRange, meanSweep, 'k-', 'LineWidth', 2);
    hold on;
    fill([cutRange fliplr(cutRange)], ...
        [meanSweep + semSweep; flipud(meanSweep - semSweep)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(cutRange, yFit, 'r--', 'LineWidth', 2);

    xlabel('Sweep %', 'FontSize', 10);
    ylabel('Norm. ΔScore', 'FontSize', 10);
    title(sprintf('%s\nr=%.2f', cleanSubsystemNames_local(favSubsys), r_pearson), ...
        'FontSize', 9, 'Interpreter', 'none');
    grid on; box on;
    ylim([-0.3 0.6]);
end

if saveTF
    saveas(fig5, fullfile('figures', 'Fig8_FavRegression.png'));
    savefig(fig5, fullfile('figures', 'Fig8_FavRegression.fig'));
    fprintf('Saved Figure 5\n');
end

%% Figure 5 better

%% FIGURE 5: User-Specified Favorite Subsystems (Modified for Single Plot)
FavSubSystems = ["Omega-3 fatty acid metabolism";
    "Omega-6 fatty acid metabolism";
    "Inositol phosphate metabolism";
    "Carbohydrate metabolism";
    "Beta oxidation of fatty acids";
    "Nucleotide metabolism";
    "Glycine serine and threonine metabolism";
    "Exchange";
    "Purine metabolism";
    "Central carbon metabolism"];
subsystemCorr.Properties.RowNames = subsystemCorr.SubSystem;
FavCor = subsystemCorr(FavSubSystems, :);
fprintf('\nGenerating Figure 5: Favorite subsystems sweep scores...\n');

% Compute statistics for all favorite subsystems
numSubs = height(FavCor);
meanSweeps = cell(numSubs, 1);
semSweeps = cell(numSubs, 1);
subIdxs = zeros(numSubs, 1);
for b = 1:numSubs
    favSubsys = FavCor.SubSystem{b};
    favIdx = find(strcmp(uniqueSubSystems, favSubsys));
    if isempty(favIdx)
        continue;
    end
    subIdxs(b) = favIdx;
    favData = squeeze(cutSweepScores(favIdx, :, :));
    processedData = process_sweep_data_local(favData, cutRange);
    meanSweep = mean(processedData, 1, 'omitnan')';
    validPerSweep = sum(~isnan(processedData), 1)';
    semSweep = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));
    meanSweeps{b} = meanSweep;
    semSweeps{b} = semSweep;
end

% Find index for +10%
pos10_idx = find(cutRange == 10, 1);
if isempty(pos10_idx)
    error('No +10%% point found in cutRange');
end

% Compute mean at +10% for sorting
mean_at10 = zeros(numSubs, 1);
for b = 1:numSubs
    mean_at10(b) = meanSweeps{b}(pos10_idx);
end

% Sort indices descending by mean_at10 (highest first)
[~, sort_order] = sort(mean_at10, 'descend');

% Define low to high order for fills and second pass
low_to_high_order = sort_order(end:-1:1);

% Create figure
fig5 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 14], 'Color', 'w');
hold on;

% Define colors (using MATLAB's lines colormap for 10 distinct colors)
colors = [
    0.1216 0.4667 0.7059;  % Blue
    0.8392 0.1529 0.1569;  % Red
    0.1725 0.6275 0.1725;  % Green
    0.9804 0.5020 0.4471;  % Orange
    0.5804 0.4039 0.7412;  % Purple
    0.3010 0.7450 0.9330;  % Teal
    1.0000 0.8392 0.0000;  % Yellow
    0.65098 0.4627 0.1137; % Brown
    0.4940 0.1840 0.5560;  % Deep Purple
    0.4660 0.6740 0.1880   % Forest Green
];


% Plot error fills first, from lowest to highest (highest fill on top)
for bb = 1:10
    orig_idx = low_to_high_order(bb);
    meanSweep = meanSweeps{orig_idx};
    semSweep = semSweeps{orig_idx};
    col = colors(orig_idx, :);
    x_fill = [cutRange fliplr(cutRange)];
    y_fill = [meanSweep + semSweep; flipud(meanSweep - semSweep)];
    fill(x_fill, y_fill, col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% First pass: plot lines from highest to lowest (for legend handles, highest first)
handles = gobjects(10, 1);
legend_names = cell(10, 1);
for bb = 1:10
    orig_idx = sort_order(bb);
    favSubsys = FavCor.SubSystem{orig_idx};
    meanSweep = meanSweeps{orig_idx};
    col = colors(orig_idx, :);
    legend_names{bb} = cleanSubsystemNames_local(favSubsys);
    handles(bb) = plot(cutRange, meanSweep, '-', 'LineWidth', 2, 'Color', col);
end

% Second pass: plot lines from lowest to highest (puts highest line on top)
for bb = 1:10
    orig_idx = low_to_high_order(bb);
    meanSweep = meanSweeps{orig_idx};
    col = colors(orig_idx, :);
    plot(cutRange, meanSweep, '-', 'LineWidth', 2, 'Color', col);
end

% Formatting
title('Metabolic Perturbation vs Predicted Toxicity', 'FontSize', 12);
xlabel('Subsytem Delta Flux %', 'FontSize', 11);
ylabel('Norm. ΔScore', 'FontSize', 11);
legend(handles, legend_names, 'Location', 'bestoutside', 'FontSize', 9, 'Interpreter', 'none');
grid on;
box on;
%ylim([-0.3 0.6]);

if saveTF
    saveas(fig5, fullfile('figures', 'Fig5_FavSweepScores.png'));
    savefig(fig5, fullfile('figures', 'Fig5_FavSweepScores.fig'));
    fprintf('Saved Figure 5\n');
end


%% VISUALIZATION: All Subsystem Sweeps with Full Range
fprintf('\nGenerating tiled visualization for all subsystems (full range)...\n');

% Create figure for tiled layout
fig_all = figure();

t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

% Process full sweep data for all subsystems (no cutRange)
numSubsys = length(uniqueSubSystems);
processedFull = nan(numSubsys, length(allModels), length(sweepRange));

zeroIdx = find(abs(sweepRange) < 1e-6, 1);  % Index for 0% sweep

for i = 1:numSubsys
    rawData = squeeze(allSweepScores(i, :, :));  % models x sweeps

    for m = 1:size(rawData, 1)
        modelScores = rawData(m, :);

        % Handle inf values (similar to process_sweep_data_local)
        isInfPos = modelScores == inf;
        isInfNeg = modelScores == -inf;
        if any(isInfPos) || any(isInfNeg)
            finite = modelScores(~(isInfPos | isInfNeg));
            if ~isempty(finite)
                maxF = max(finite);
                minF = min(finite);
                rngF = max(1, maxF - minF);
                modelScores(isInfPos) = maxF + rngF;
                modelScores(isInfNeg) = minF - rngF;
            else
                modelScores(:) = 0;
            end
        end

        % Compute delta relative to 0% change
        baseScore = modelScores(zeroIdx);
        if isfinite(baseScore)
            modelScores = modelScores - baseScore;
        else
            finiteBase = modelScores(isfinite(modelScores));
            if ~isempty(finiteBase)
                baseScore = mean(finiteBase);
                modelScores = modelScores - baseScore;
            end
        end

        % Normalize by max absolute deviation
        maxDev = max(abs(modelScores));
        if maxDev > 0
            modelScores = modelScores / maxDev;
        end

        processedFull(i, m, :) = modelScores;
    end
end

% Generate plots for all subsystems
for i = 1:numSubsys
    nexttile;

    subsys = uniqueSubSystems(i);
    data = squeeze(processedFull(i, :, :));  % models x sweeps

    % Compute mean and SEM across models
    meanSweep = mean(data, 1, 'omitnan')';
    validPerSweep = sum(~isnan(data), 1)';
    semSweep = std(data, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

    % Convert sweepRange to percentage for x-axis
    xPercent = sweepRange * 100;

    % Plot mean line
    plot(xPercent, meanSweep, 'k-', 'LineWidth', 1.5);
    hold on;

    % Shaded SEM
    fill([xPercent fliplr(xPercent)], ...
        [meanSweep + semSweep; flipud(meanSweep - semSweep)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Title with cleaned subsystem name
    title(cleanSubsystemNames_local(subsys), 'FontSize', 9, 'Interpreter', 'none');

    % No xlabel, ylabel, or legend
    grid on; box on;

    % Set consistent y-limits to visualize plateaus
    %ylim([-1 1]);
    %xlim([-50 50]);
end

% Save the figure
if saveTF
    saveas(fig_all, fullfile('figures', 'AllSubsysFullSweeps.png'));
    savefig(fig_all, fullfile('figures', 'AllSubsysFullSweeps.fig'));
end

fprintf('Saved all-subsystem visualization (full range)\n');

%% INTER-SUBSYSTEM CORRELOGRAMS
% Compute and plot correlograms showing correlations between subsystem responses
% during sweeps of each target subsystem (sweeps as "samples")
fprintf('\n=== GENERATING INTER-SUBSYSTEM CORRELOGRAMS ===\n');
fprintf('This analyzes how other subsystems co-vary during each target sweep.\n');

numValidSubsys = sum(subsystemCorr.NumReactions > 0);  % Skip empty
for i = 1:numSubsys
    if subsystemCorr.NumReactions(i) == 0
        fprintf('Skipping empty subsystem %d/%d: %s\n', i, numSubsys, uniqueSubSystems(i));
        continue;
    end
    
    targetSubsys = uniqueSubSystems(i);
    subsys_changes = allSubsysAverages{i};  % sweeps x numSubsys (fractional changes)
    
    % Remove columns for empty subsystems (all NaN or 0) to avoid singular corr
    validCols = true(1, numSubsys);
    for j = 1:numSubsys
        if all(isnan(subsys_changes(:, j))) || all(subsys_changes(:, j) == 0)
            validCols(j) = false;
        end
    end
    subsys_changes_valid = subsys_changes(:, validCols);
    validSubsysNames = uniqueSubSystems(validCols);
    
    if sum(validCols) < 2
        fprintf(' [%d/%d] %s: Insufficient valid subsystems - skipping\n', i, numValidSubsys, targetSubsys);
        continue;
    end
    
    % Compute Pearson correlation matrix (rows=sweep steps as observations/samples)
    % corr() uses pairwise complete obs, handles NaNs
    corr_matrix = corr(subsys_changes_valid);
    
    % Plot correlogram as heatmap
    fig_corr = figure('Units', 'centimeters', 'Color', 'w');
    imagesc(corr_matrix);
    colormap(bwr);  %TODO define colormap. Red-blue for positive/negative correlations (-1 to +1)
    caxis([-1 1]);
    colorbar('Ticks', [-1 -0.5 0 0.5 1], 'FontSize', 10);
    
    % Labels (use cleaned names, rotate x for readability)
    cleanedNames = cleanSubsystemNames_local(validSubsysNames);
    numValid = length(validSubsysNames);
    set(gca, 'XTick', 1:numValid, 'XTickLabel', cleanedNames, ...
             'YTick', 1:numValid, 'YTickLabel', cleanedNames, ...
             'FontSize', 8, 'TickLabelInterpreter', 'none');
    xtickangle(45);
    xlabel('Affected Subsystem', 'FontSize', 11, 'FontName', 'Arial');
    ylabel('Affected Subsystem', 'FontSize', 11, 'FontName', 'Arial');
    title(sprintf('Correlogram: Subsystem Co-Variations During %s Sweep (-10%% to +10%%)', ...
                  cleanSubsystemNames_local(targetSubsys)), ...
          'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
    axis square; grid on;
    
    if saveTF
        safeName = strrep(string(targetSubsys), {' ', '-', '(', ')', '/', '\\'}, {'_', '_', '', '', '_', '_'});
        saveas(fig_corr, fullfile('figures', sprintf('Correlogram_Sweep_%s.png', safeName)));
        savefig(fig_corr, fullfile('figures', sprintf('Correlogram_Sweep_%s.fig', safeName)));
        close(fig_corr);
        fprintf(' [%d/%d] Saved correlogram for %s\n', i, numValidSubsys, targetSubsys);
    else
        fprintf(' [%d/%d] Generated correlogram for %s (close to proceed)\n', i, numValidSubsys, targetSubsys);
    end
end
fprintf('Inter-subsystem correlograms complete! (%d generated)\n');


%% Save Results
fprintf('\nSaving analysis results...\n');
if saveTF
    save(fullfile('out', 'driver8_results.mat'), 'subsystemCorr', ...
        'predictorNames', 'allSweepScores', 'sweepRange', 'topPos', 'topNeg', ...
        'FavCor', 'uniqueSubSystems','allSubsysAverages', '-v7.3');
end

fprintf('\n=== DRIVER 8 COMPLETED SUCCESSFULLY ===\n');
fprintf('Generated 5 figures with mass-balanced subsystem sweeps:\n');
fprintf(' - figures/Fig8_SubsystemPositiveCorr.png\n');
fprintf(' - figures/Fig8_SubsystemNegativeCorr.png\n');
fprintf(' - figures/Fig8_TopPosRegression.png\n');
fprintf(' - figures/Fig8_TopNegRegression.png\n');
fprintf(' - figures/Fig8_FavRegression.png\n');

%close all;

%% LOCAL HELPER FUNCTIONS -------------------------------
%% create_mass_balanced_sweep_local
% Creates mass-balanced sweep input matrix using rigorous optimization
function inputMatrix_MB = create_mass_balanced_sweep_local(Model, subsysIdx, sweepRange, J0)
n = length(Model.rxns);
numSweeps = length(sweepRange);
inputMatrix_MB = zeros(numSweeps, n);

% For each sweep value, solve for mass-balanced flux distribution
for i = 1:numSweeps
    sweepVal = sweepRange(i);

    % Create perturbation vector V_g
    V_g = ones(n, 1);
    V_g(subsysIdx) = sweepVal+1;

    % Solve for V_d to restore mass balance
    V_d = solve_V_d_for_sweep_local(Model, V_g, J0);

    % Calculate new flux
    J_new = J0 .* (V_d + V_g);

    % Convert to percent change
    deltaJ = (J_new - J0) ./ J0;
    deltaJ(isnan(deltaJ) | isinf(deltaJ)) = 0;

    inputMatrix_MB(i, :) = deltaJ';
end
%delete(gcp('nocreate'))
end

%% solve_V_d_for_sweep_local
% Solves optimization: minimize ||V_d||^2 subject to S*(J0.*(V_d + V_g)) = 0
function V_d = solve_V_d_for_sweep_local(Model, V_g, J0)
% Set up optimization problem
xOpty = optimvar('x', Model.rxnNames);
quadprobN = optimproblem();

% Objective: minimize ||V_d||^2
quadprobN.Objective = sum(xOpty .* xOpty);

% Constraint: S * (J0 .* (V_d + V_g)) = 0 (steady-state mass balance)
quadprobN.Constraints.sj0 = full(Model.S) * (J0 .* (xOpty + V_g)) == 0;

linsol = solve(quadprobN, 'Options', struct('Display', 'none'));

V_d = linsol.x;
end

%% process_sweep_data_local
% Processes sweep score data: handles inf, computes delta, normalizes
function processedData = process_sweep_data_local(rawData, sweepRange)
% rawData: models × sweeps
zeroIdx = find(abs(sweepRange) < 1e-6, 1);
if isempty(zeroIdx)
    zeroIdx = round(size(rawData, 2) / 2);
end

numModels = size(rawData, 1);
numSweeps = size(rawData, 2);
processedData = nan(numModels, numSweeps);

for m = 1:numModels
    modelScores = rawData(m, :);

    % Handle inf values
    isInfPos = modelScores == inf;
    isInfNeg = modelScores == -inf;
    if any(isInfPos) || any(isInfNeg)
        finite = modelScores(~(isInfPos | isInfNeg));
        if ~isempty(finite)
            maxF = max(finite);
            minF = min(finite);
            rngF = max(1, maxF - minF);
            modelScores(isInfPos) = maxF + rngF;
            modelScores(isInfNeg) = minF - rngF;
        else
            modelScores(:) = 0;
        end
    end

    % Compute delta relative to 0% change
    baseScore = modelScores(zeroIdx);
    if isfinite(baseScore)
        modelScores = modelScores - baseScore;
    else
        finiteBase = modelScores(isfinite(modelScores));
        if ~isempty(finiteBase)
            baseScore = mean(finiteBase);
            modelScores = modelScores - baseScore;
        end
    end

    % Normalize by max absolute deviation
    maxDev = max(abs(modelScores));
    if maxDev > 0
        modelScores = modelScores / maxDev;
    end

    processedData(m, :) = modelScores;
end
end

%% cleanSubsystemNames_local
% Cleans subsystem names for better visualization
function cleaned = cleanSubsystemNames_local(subsystemNames)
cleaned = string(subsystemNames);

% Remove/shorten common words
cleaned = regexprep(cleaned, ' metabolism', '');
cleaned = regexprep(cleaned, ' biosynthesis', ' synth.');
cleaned = regexprep(cleaned, 'Oxidative phosphorylation', 'Ox. phos.');
cleaned = regexprep(cleaned, 'Tricarboxylic acid cycle', 'TCA cycle');
cleaned = regexprep(cleaned, 'Fatty acid', 'FA');

% Capitalize first letter
for i = 1:length(cleaned)
    if strlength(cleaned(i)) > 0
        cleaned(i) = upper(extractBefore(cleaned(i), 2)) + ...
            extractAfter(cleaned(i), 1);
    end
end
end
