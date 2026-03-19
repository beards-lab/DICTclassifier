%%  step07_subsystem_importance

%  Analyzes predictor importance by metabolic subsystem using
%  hypergeometric enrichment. Performs subsystem correlation sweeps (+/-50%)
%  to identify subsystems positively/negatively correlated with predicted
%  cardiotoxicity. Generates bubble plots and regression tiled layouts.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (predictorImportance, predict, hygecdf,
%                                    fitlm, bubblechart)
%
%  Inputs:
%    - out/iCardio_optimized.mat                     — ModelOpt with
%                                                      subsystem annotations
%    - out/deltaFlux_stats.mat                       — RxnStatsMedMed
%    - out/ROR_results.mat                           — regresTable
%    - data/drugMap.xlsx                             — Drug mapping table
%    - out/DICT_L0X/c1_<date>/Rxn_models_<a>.mat    — 400 trained Rxn
%                                                      classifiers
%
%  Outputs:
%    - out/driver7_results.mat       — subsystemStats, subsystemCorr,
%                                      allSweepScores, allImportance
%    - figures/Fig7_SubsystemImportance.png/.fig
%    - figures/Fig1b_SubsystemImportance_Bubble_Clean.png/.fig
%    - figures/Fig7_SubsystemPositiveCorr.png/.fig
%    - figures/Fig7_SubsystemNegativeCorr.png/.fig
%    - figures/Fig7_TopPosRegression.png/.fig
%    - figures/Fig7_TopNegRegression.png/.fig
%    - figures/Fig7_FavRegression.png/.fig

clear; clc; close all;

%% Parameters
expName = "DICT_L0X";
currentDateString = "10_6_30hr";
iRange = 1:100;
cMax = 4;
sweepRange = -0.5:0.001:0.5; % Sweep from -0.5 to +0.5 in 0.01 steps
topN = 10; % Number of top subsystems to display

%% Load Model and Trained Classifiers
fprintf('Loading optimized model and trained classifiers...\n');
load(fullfile('out','iCardio_optimized.mat'), 'ModelOpt');
Model = ModelOpt;
clear ModelOpt;
fprintf('Model loaded: %d reactions, %d metabolites\n', length(Model.rxns), length(Model.mets));

%% Load Feature Matrix to Get Reaction Names
fprintf('Loading feature matrix for reaction name mapping...\n');
load(fullfile('out','deltaFlux_stats.mat'), 'RxnStatsMedMed');
load(fullfile('out','ROR_results.mat'), 'regresTable');
drugMap = readtable(fullfile("data", "drugMap.xlsx"));
predictorNames = "x"+string(int2str([1:length(Model.rxns)]'));

%% Extract Predictor Importance from All 400 Models
fprintf('Extracting predictor importance from 400 trained models...\n');
basePath = fullfile("out", expName, "c1_"+ currentDateString);
numModels = length(iRange) * cMax; % 100 iterations × 4 folds = 400 models
% Initialize storage
allImportance = nan(length(predictorNames), numModels);
modelCounter = 1;
for idx = 1:length(iRange)
    a = iRange(idx);
    % Load models for this iteration
    modelFile = fullfile(basePath, sprintf('Rxn_models_%d.mat', a));
    if ~exist(modelFile, 'file')
        warning('Model file not found: %s', modelFile);
        continue;
    end
    load(modelFile, 'tRxn_models');
    for c = 1:cMax
        try
            mdl = tRxn_models{1,c};
            % Get predictor importance
            imp = predictorImportance(mdl);
            % Store importance (pad with NaN if sizes don't match)
            if length(imp) <= size(allImportance, 1)
                allImportance(1:length(imp), modelCounter) = imp;
            else
                warning('Importance vector larger than expected for model %d', modelCounter);
            end
            modelCounter = modelCounter + 1;
        catch ME
            warning('Error processing model %d-%d: %s', a, c, ME.message);
        end
    end
    if mod(idx, 10) == 0
        fprintf(' Processed %d/%d iterations\n', idx, length(iRange));
    end
end
fprintf('Successfully extracted importance from %d models\n', modelCounter-1);

%% Aggregate Importance by SubSystem
fprintf('Aggregating importance by subsystem...\n');
uniqueSubSystems = unique(string(Model.subSystems));
uniqueSubSystems(uniqueSubSystems == "") = [];
subsystemStats = table();
subsystemStats.SubSystem = uniqueSubSystems;
subsystemStats.NumReactions = zeros(length(uniqueSubSystems), 1);
subsystemStats.MeanImportance = zeros(length(uniqueSubSystems), 1);
subsystemStats.SEMImportance = zeros(length(uniqueSubSystems), 1);
subsystemStats.EnrichmentPval = ones(length(uniqueSubSystems), 1);
% For each subsystem, calculate statistics
totalPredictors = height(predictorNames);
allImportance(ismissing(allImportance)) = 0;
topPredictorThreshold = prctile(allImportance(:), 95, 'all'); % Top 5% important predictors
for i = 1:length(uniqueSubSystems)
    subsys = uniqueSubSystems(i);
    % Find predictors in this subsystem
    inSubsys = string(Model.subSystems) == subsys;
    subsystemStats.NumReactions(i) = sum(inSubsys);
    % Get importance values for predictors in this subsystem
    subsysIdx = find(inSubsys);
    if ~isempty(subsysIdx) && subsysIdx(end) <= size(allImportance, 1)
        subsysImportance = allImportance(subsysIdx, :);
        subsystemStats.MeanImportance(i) = mean(subsysImportance(:), 'omitnan');
        subsystemStats.SEMImportance(i) = std(subsysImportance(:), 'omitnan') / sqrt(sum(~isnan(subsysImportance(:))));
        % Hypergeometric enrichment test
        % Test if top-important predictors are enriched in this subsystem
        numTopInSubsys = sum(subsysImportance > topPredictorThreshold, 'all');
        subsystemStats.NumTopTot(i) = numTopInSubsys;
        numInSubsys = numel(subsysImportance);
        subsystemStats.NumTot(i) = numInSubsys;
        numTopTotal = sum(allImportance > topPredictorThreshold, 'all');
        numTotal = numel(allImportance);
        % Hypergeometric test: P(X >= k) where k = numTopInSubsys
        try
            pval = hygecdf(numTopInSubsys-1, numTotal, numTopTotal, numInSubsys, 'upper');
            subsystemStats.EnrichmentPval(i) = min(pval.*44,1);
        catch
            subsystemStats.EnrichmentPval(i) = 1;
        end
    end
end
% Sort by enrichment p-value
subsystemStats = sortrows(subsystemStats, 'EnrichmentPval', 'ascend');
subsystemStats.NormImportance = normalize(subsystemStats.MeanImportance);
fprintf('Calculated statistics for %d subsystems\n', height(subsystemStats));

%% FIGURE 1: Overall Subsystem Importance with Hypergeometric Enrichment
fprintf('\nGenerating Figure 1: Overall subsystem importance...\n');
% Select top N subsystems by enrichment
topSubsystems = subsystemStats(1:min(topN, height(subsystemStats)), :);
fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
hold on;
% Create horizontal bar plot
yPos = 1:height(topSubsystems);
barh(yPos, topSubsystems.MeanImportance, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1);
% Add error bars (SEM)
errorbar(topSubsystems.MeanImportance, yPos, topSubsystems.SEMImportance, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 6);
% Format axes
set(gca, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames(topSubsystems.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('Mean Predictor Importance', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystem Importance for Cardiotoxicity Prediction', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(topSubsystems.MeanImportance + topSubsystems.SEMImportance) * 1.1]);
ylim([0.5 height(topSubsystems) + 0.5]);
grid on;
box on;
% Add significance stars
for i = 1:height(topSubsystems)
    pval = topSubsystems.EnrichmentPval(i);
    if pval < 0.001
        sigText = '***';
    elseif pval < 0.01
        sigText = '**';
    elseif pval < 0.05
        sigText = '*';
    else
        sigText = '';
    end
    if ~isempty(sigText)
        xPos = topSubsystems.MeanImportance(i) + topSubsystems.SEMImportance(i) * 1.2;
        text(xPos, i, sigText, 'FontSize', 14, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    % Add reaction count
    numRxns = topSubsystems.NumReactions(i);
    text(0.02 * max(xlim), i, sprintf('n=%d (%.2g)', numRxns, topSubsystems.NumTopTot(i)./topSubsystems.NumTot(i).*100), ...
        'FontSize', 8, 'Color', 'k', 'FontWeight', 'normal', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
% Save figure
if ~exist('figures', 'dir'), mkdir('figures'); end
saveas(fig1, fullfile('figures', 'Fig7_SubsystemImportance.png'));
savefig(fig1, fullfile('figures', 'Fig7_SubsystemImportance.fig'));
fprintf('Saved Figure 1\n');

%% FIGURE 1b: Subsystem Importance Bubble Plot (Cleaned with bubblechart)
fprintf('\nGenerating Figure 1b: Subsystem importance bubble plot...\n');
% Prepare data: Sort subsystems by mean importance descending for better ordering
topSubsystems = sortrows(topSubsystems, 'MeanImportance', 'descend');
topSubsystems.PropTop = topSubsystems.NumTopTot ./ topSubsystems.NumTot; % Proportion of top reactions
topSubsystems.NegLogP = -log10(topSubsystems.EnrichmentPval); % -log10(p-value) for x-axis (significance)
topSubsystems.NegLogP(isinf(topSubsystems.NegLogP)) = max(topSubsystems.NegLogP(~isinf(topSubsystems.NegLogP)));
% Create figure
fig1b = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
% Y positions for subsystems
yPos = 1:height(topSubsystems);
% X values: -log10(p-value) to emphasize significance
xVals = topSubsystems.NegLogP;
% Bubble sizes: Use NumReactions directly for relative sizing
sz = topSubsystems.NumReactions;
% Colors: Map shifted importance to upper parula indices (51-100) for medium-to-high tones
minImp = min(subsystemStats.NormImportance);
maxImp = max(subsystemStats.NormImportance);
propColors = normcdf(topSubsystems.NormImportance)*100;
c_indices = propColors;
% Create bubble chart
bc = bubblechart(xVals, yPos, sz, c_indices, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
line([-log10(0.05) -log10(0.05)], [0 11])
% Set colormap to full parula, but colors will use upper range
colormap(parula(100));
% Adjust bubble size range for readability (diameters in points)
bubblesize([10 50]); % Small to large diameters, adjustable as needed
% Add bubble size legend
bl = bubblelegend('Number of Reactions', 'Location', 'northwest', 'FontName', 'Arial', 'FontSize', 10, 'Style','telescopic', 'NumBubbles',2);
% Colorbar with adjusted labels to reflect high importance
cbar = colorbar;
minC = 15.9;
maxC = 84.1;
midC = round((minC + maxC) / 2);
cbar.Ticks = [minC, midC, maxC];
cbar.TickLabels = {'Low', 'Average', 'High'};
cbar.Label.String = 'Predictor Importance';
cbar.FontName = 'Arial';
cbar.FontSize = 11;
clim([0 100])
% Format axes
ax = gca;
set(ax, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames(topSubsystems.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('-log_{10}(P_{adj})', 'FontName', 'Arial', 'FontSize', 12, ...
    'Interpreter', 'tex');
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystem Enrichment and Importance (Top 10)', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, max(xVals)*1.1]);
ylim([0.5, height(topSubsystems) + 0.5]);
grid on;
box on;
% Save figure
saveas(fig1b, fullfile('figures', 'Fig1b_SubsystemImportance_Bubble_Clean.png'));
savefig(fig1b, fullfile('figures', 'Fig1b_SubsystemImportance_Bubble_Clean.fig'));
fprintf('Saved cleaned Figure 1b\n');

%% Load All Models for Correlation Analysis
fprintf('\nLoading all models for correlation analysis...\n');
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
fprintf('Loaded %d models for correlation analysis\n', length(allModels));

%% SUBSYSTEM CORRELATION SWEEP (MODIFIED: Save Full Sweep Scores)
fprintf('\nPerforming subsystem correlation sweep...\n');
subsystemCorr = table();
subsystemCorr.SubSystem = uniqueSubSystems;
subsystemCorr.MeanCorrelation = zeros(length(uniqueSubSystems), 1);
subsystemCorr.SEMCorrelation = zeros(length(uniqueSubSystems), 1);
subsystemCorr.NumReactions = zeros(length(uniqueSubSystems), 1);
% MODIFIED: Initialize 3D array for all sweep scores (subsys x model x sweep)
numSubsys = length(uniqueSubSystems);
allSweepScores = nan(numSubsys, length(allModels), length(sweepRange));
% For each subsystem
for i = 1:length(uniqueSubSystems)
    disp(i);
    subsys = uniqueSubSystems(i);
    % Find predictors in this subsystem
    inSubsys = string(Model.subSystems) == subsys;
    subsysIdx = find(inSubsys);
    subsystemCorr.NumReactions(i) = length(subsysIdx);
    if isempty(subsysIdx) || subsysIdx(end) > size(allImportance, 1)
        continue;
    end
    % Store correlations for each model
    corrPerModel = nan(length(allModels), 1);
    for m = 1:length(allModels)
        %disp(i + "_" + m)
        mdl = allModels{m};
        numSweeps = length(sweepRange);
        inputMatrix = zeros(numSweeps, length(predictorNames)); % Rows: sweeps; Columns: predictors
        inputMatrix(:, subsysIdx) = repmat(sweepRange', 1, length(subsysIdx));
        [~, scores] = predict(mdl, inputMatrix);
        sweepScores = scores(:, 2); % Class 1 (toxic) probability for all sweeps
        % MODIFIED: Store full sweepScores
        allSweepScores(i, m, :) = sweepScores;
        if any(isnan(sweepScores))
            disp("Error")
        end
        validScores = ~isnan(sweepScores);
        corrPerModel(m) = corr(sweepRange(validScores)', sweepScores(validScores), 'Type', 'Spearman');
    end
    % Aggregate across models
    subsystemCorr.MeanCorrelation(i) = mean(corrPerModel, 'omitnan');
    subsystemCorr.SEMCorrelation(i) = std(corrPerModel, 'omitnan') / sqrt(sum(~isnan(corrPerModel)));
    if mod(i, 10) == 0
        fprintf(' Processed %d/%d subsystems\n', i, length(uniqueSubSystems));
    end
end
fprintf('Correlation sweep complete\n');

%% Identify Top Subsystems for Regression Plots (MODIFIED: New)
% Sort for positive and negative
posCorr = subsystemCorr(subsystemCorr.MeanCorrelation > 0, :);
posCorr = sortrows(posCorr, 'MeanCorrelation', 'descend');
negCorr = subsystemCorr(subsystemCorr.MeanCorrelation < 0, :);
negCorr = sortrows(negCorr, 'MeanCorrelation', 'ascend');

%% FIGURE 2: Subsystems Positively Correlated with Cardiotoxicity
fprintf('\nGenerating Figure 2: Positive correlations...\n');
topPos = posCorr(1:min(topN, height(posCorr)), :);
fig2 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
hold on;
yPos = 1:height(topPos);
barh(yPos, topPos.MeanCorrelation, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1);
% Add error bars
errorbar(topPos.MeanCorrelation, yPos, topPos.SEMCorrelation, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 6);
% Format axes
set(gca, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames(topPos.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('Mean Correlation with Cardiotoxicity Score', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystems Positively Correlated with Cardiotoxicity', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(topPos.MeanCorrelation + topPos.SEMCorrelation) * 1.1]);
ylim([0.5 height(topPos) + 0.5]);
grid on;
box on;
% Add reaction counts
for i = 1:height(topPos)
    numRxns = topPos.NumReactions(i);
    text(0.02 * max(xlim), i, sprintf('n=%d', numRxns), ...
        'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
% Save figure
saveas(fig2, fullfile('figures', 'Fig7_SubsystemPositiveCorr.png'));
savefig(fig2, fullfile('figures', 'Fig7_SubsystemPositiveCorr.fig'));
fprintf('Saved Figure 2\n');

%% FIGURE 3: Subsystems Negatively Correlated with Cardiotoxicity
fprintf('\nGenerating Figure 3: Negative correlations...\n');
topNeg = negCorr(1:min(topN, height(negCorr)), :);
fig3 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color', 'w');
hold on;
yPos = 1:height(topNeg);
barh(yPos, abs(topNeg.MeanCorrelation), 'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'k', 'LineWidth', 1);
% Add error bars
errorbar(abs(topNeg.MeanCorrelation), yPos, topNeg.SEMCorrelation, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 6);
% Format axes
set(gca, 'YTick', yPos, 'YTickLabel', cleanSubsystemNames(topNeg.SubSystem), ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel('Mean |Correlation| with Cardiotoxicity Score', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Metabolic Subsystem', 'FontName', 'Arial', 'FontSize', 12);
title('Subsystems Negatively Correlated with Cardiotoxicity', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(abs(topNeg.MeanCorrelation) + topNeg.SEMCorrelation) * 1.1]);
ylim([0.5 height(topNeg) + 0.5]);
grid on;
box on;
% Add reaction counts
for i = 1:height(topNeg)
    numRxns = topNeg.NumReactions(i);
    text(0.02 * max(xlim), i, sprintf('n=%d', numRxns), ...
        'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
% Save figure
saveas(fig3, fullfile('figures', 'Fig7_SubsystemNegativeCorr.png'));
savefig(fig3, fullfile('figures', 'Fig7_SubsystemNegativeCorr.fig'));
fprintf('Saved Figure 3\n');

%% MODIFIED: New Figures for Top Subsystem Regression Plots
fig4 = figure('Units', 'centimeters', 'Position', [2, 2, 12, 10], 'Color', 'w');
tiledlayout('flow');
for a = 1:height(topPos)
    nexttile
    topPosSubsys = topPos.SubSystem(a);
    topPosIdx = find(uniqueSubSystems == topPosSubsys);
    fprintf('\nGenerating Figure 4: Top positive subsystem regression...\n');
    % Extract mean and SEM across models for top positive
    topPosData = squeeze(allSweepScores(topPosIdx, :, :)); % models x sweeps
    topPosData(topPosData == inf) = max(allSweepScores(~isinf(allSweepScores)),[], 'all');
    topPosData(topPosData == -inf) = min(allSweepScores(~isinf(allSweepScores)),[], 'all');
    meanSweepPos = mean(topPosData, 1, "omitmissing")'; % sweeps x 1
    semSweepPos = std(topPosData, 0, 1, "omitmissing")' / sqrt(sum(~isnan(topPosData), "all")); % sweeps x 1
    % Linear fit (polyfit degree 1)
    p_pos = polyfit(sweepRange, meanSweepPos, 1);
    yFitPos = polyval(p_pos, sweepRange);
    % Plot
    plot(sweepRange, meanSweepPos, 'b-', 'LineWidth', 2); % Mean
    hold on;
    % ±SEM shading
    fill([sweepRange fliplr(sweepRange)], [meanSweepPos + semSweepPos; flipud(meanSweepPos - semSweepPos)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(sweepRange, yFitPos, 'r--', 'LineWidth', 2); % Linear fit
    xlabel('Sweep % Change', 'FontName', 'Arial'); %, 'FontSize', 12);
    ylabel('Predicted Cardiotoxicity Score', 'FontName', 'Arial', 'FontSize', 12);
    title(sprintf('Top Positive Subsystem: %s (r = %.3f)', topPosSubsys, posCorr.MeanCorrelation(1)), ...
        'FontName', 'Arial'); %, 'FontSize', 14, 'FontWeight', 'bold');
    legend('Mean Score', '±SEM', 'Linear Fit', 'Location', 'best');
    grid on; box on; xlim([-0.5 0.5]);
    % Save
    saveas(fig4, fullfile('figures', 'Fig7_TopPosRegression.png'));
    savefig(fig4, fullfile('figures', 'Fig7_TopPosRegression.fig'));
    fprintf('Saved Figure 4\n');
end

%%
% MODIFIED: New Figures for Top Subsystem Regression Plots
% Process each top positive and negative subsystem with delta scores and normalization

fig4 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 14], 'Color', 'w');
tiledlayout('flow', 'TileSpacing', 'compact');
minRangeIdx = find(abs(sweepRange + 0.1) < 1e-12, 1);
maxRangeIdx = find(abs(sweepRange - 0.1) < 1e-12, 1);
cutRange = sweepRange(minRangeIdx:maxRangeIdx).*100;
cutSweepScores = allSweepScores(:,:, minRangeIdx:maxRangeIdx);
for a = 1:height(topPos)
    nexttile
    topPosSubsys = topPos.SubSystem(a);
    topPosIdx = find(uniqueSubSystems == topPosSubsys);
    fprintf('Generating Figure 4: Top positive subsystem %d/%d regression...\n', a, height(topPos));

    % Extract data across models for this subsystem
    topPosData = squeeze(cutSweepScores(topPosIdx, :, :)); % models x sweeps

    % Find zero index
    zeroIdx = find(cutRange == 0, 1);
    if isempty(zeroIdx)
        warning('Zero sweep not found for subsystem %s', topPosSubsys);
        continue;
    end

    % Handle inf and compute delta and normalize per model
    numModels = size(topPosData, 1);
    numSweeps = length(cutRange);
    processedData = nan(numModels, numSweeps);
    for m = 1:numModels
        modelScores = topPosData(m, :);

        % Handle inf per model
        isInfPos = modelScores == inf;
        isInfNeg = modelScores == -inf;
        if any(isInfPos) || any(isInfNeg)
            finite = modelScores(~(isInfPos | isInfNeg));
            if ~isempty(finite)
                maxF = max(finite);
                minF = min(finite);
                rngF = max(1, maxF - minF); % Ensure positive range
                modelScores(isInfPos) = maxF + rngF;
                modelScores(isInfNeg) = minF - rngF;
            else
                % All inf, set to 0 (unlikely)
                modelScores(:) = 0;
            end
        end

        % Compute delta score relative to 0% change
        baseScore = modelScores(zeroIdx);
        if isfinite(baseScore)
            modelScores = modelScores - baseScore;
        else
            % If base is non-finite (after handling, unlikely), use mean of finite
            finiteBase = modelScores(isfinite(modelScores));
            if ~isempty(finiteBase)
                baseScore = mean(finiteBase);
                modelScores = modelScores - baseScore;
            end
        end

        % Standardize by max absolute deviation to preserve zero and scale range
        maxDev = max(abs(modelScores));
        if maxDev > 0
            modelScores = modelScores / maxDev;
        end % else remains zero

        processedData(m, :) = modelScores;
    end

    % Compute mean and SEM across models (fix SEM to be per sweep)
    meanSweepPos = mean(processedData, 1, 'omitnan')';
    validPerSweep = sum(~isnan(processedData), 1)';
    semSweepPos = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep)); % Avoid div by zero

    % Linear fit on mean
    p_pos = polyfit(cutRange, meanSweepPos, 1);
    yFitPos = polyval(p_pos, cutRange);

    % Plot
    plot(cutRange, meanSweepPos, 'b-', 'LineWidth', 2); % Mean
    hold on;
    % ±SEM shading
    upper = meanSweepPos + semSweepPos;
    lower = meanSweepPos - semSweepPos;
    fill([cutRange fliplr(cutRange)], [upper; flipud(lower)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(cutRange, yFitPos, 'r--', 'LineWidth', 2); % Linear fit
    xlabel('Sweep % Change', 'FontName', 'Arial', 'FontSize', 12);
    ylabel('Normalized Delta Score', 'FontName', 'Arial', 'FontSize', 12);
    title(sprintf('Top Positive: %s (r = %.3f)', topPosSubsys, topPos.MeanCorrelation(a)), ...
        'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
    %legend('Mean Delta', '±SEM', 'Linear Fit', 'Location', 'best');
    grid on; box on; %xlim([-0.5 0.5]);
    ylim([-0.3 1])
end
% Save after all positive plots
saveas(fig4, fullfile('figures', 'Fig7_TopPosRegression.png'));
savefig(fig4, fullfile('figures', 'Fig7_TopPosRegression.fig'));
fprintf('Saved Figure 4 (all top positive regressions)\n');

% Now for top negative, similar loop
if height(topNeg) > 0
    fprintf('\nGenerating Figure 5: Top negative subsystems regression...\n');
    fig5 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 14], 'Color', 'w');
    tiledlayout('flow', 'TileSpacing', 'compact');
    for b = 1:height(topNeg)
        nexttile
        topNegSubsys = topNeg.SubSystem(b);
        topNegIdx = find(uniqueSubSystems == topNegSubsys);
        fprintf('Top negative subsystem %d/%d regression...\n', b, height(topNeg));

        % Extract data across models for this subsystem
        topNegData = squeeze(cutSweepScores(topNegIdx, :, :)); % models x sweeps

        % Find zero index (same as above)
        zeroIdx = find(cutRange == 0, 1);
        if isempty(zeroIdx)
            warning('Zero sweep not found for subsystem %s', topNegSubsys);
            continue;
        end

        % Handle inf and compute delta and normalize per model (same as positive)
        numModels = size(topNegData, 1);
        numSweeps = length(cutRange);
        processedData = nan(numModels, numSweeps);
        for m = 1:numModels
            modelScores = topNegData(m, :);

            % Handle inf per model
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

            % Compute delta score relative to 0% change
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

            % Standardize by max absolute deviation
            maxDev = max(abs(modelScores));
            if maxDev > 0
                modelScores = modelScores / maxDev;
            end

            processedData(m, :) = modelScores;
        end

        % Compute mean and SEM across models (per sweep)
        meanSweepNeg = mean(processedData, 1, 'omitnan')';
        validPerSweep = sum(~isnan(processedData), 1)';
        semSweepNeg = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

        % Linear fit on mean
        p_neg = polyfit(cutRange, meanSweepNeg, 1);
        yFitNeg = polyval(p_neg, cutRange);

        % Plot
        plot(cutRange, meanSweepNeg, 'g-', 'LineWidth', 2); % Mean
        hold on;
        % ±SEM shading
        upper = meanSweepNeg + semSweepNeg;
        lower = meanSweepNeg - semSweepNeg;
        fill([cutRange fliplr(cutRange)], [upper; flipud(lower)], ...
            'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(cutRange, yFitNeg, 'r--', 'LineWidth', 2); % Linear fit
        xlabel('Sweep % Change', 'FontName', 'Arial', 'FontSize', 12);
        ylabel('Normalized Delta Score', 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('Top Negative: %s (r = %.3f)', topNegSubsys, topNeg.MeanCorrelation(b)), ...
            'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
        %legend('Mean Delta', '±SEM', 'Linear Fit', 'Location', 'best');
        grid on; box on; %xlim([-0.5 0.5]);
        ylim([-0.3 1])
    end
    % Save after all negative plots
    saveas(fig5, fullfile('figures', 'Fig7_TopNegRegression.png'));
    savefig(fig5, fullfile('figures', 'Fig7_TopNegRegression.fig'));
    fprintf('Saved Figure 5 (all top negative regressions)\n');
end
%%
FavSubSystems = ["Omega-3 fatty acid metabolism";
    "Omega-6 fatty acid metabolism";
    "Inositol phosphate metabolism"
    "Carbohydrate metabolism";
    "Beta oxidation of fatty acids";
    "Nucleotide metabolism";
    "Glycine serine and threonine metabolism";
    "Exchange";
    "Purine metabolism";
    "Central carbon metabolism"
    ];
subsystemCorr.Properties.RowNames = subsystemCorr.SubSystem;
FavCor = subsystemCorr(FavSubSystems, :);
fprintf('\nGenerating Figure 6: Favorite subsystems regression...\n');
fig6 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 14], 'Color', 'w');
tiledlayout(2,5, 'TileSpacing', 'compact');
for b = 1:height(FavCor)
    nexttile
    topFavSubsys = FavCor.SubSystem(b);
    topFavIdx = find(uniqueSubSystems == topFavSubsys);
    fprintf('Fav subsystem %d/%d regression...\n', b, height(FavCor));

    % Extract data across models for this subsystem
    topFavData = squeeze(cutSweepScores(topFavIdx, :, :)); % models x sweeps

    % Find zero index (same as above)
    zeroIdx = find(cutRange == 0, 1);
    if isempty(zeroIdx)
        warning('Zero sweep not found for subsystem %s', topFavSubsys);
        continue;
    end

    % Handle inf and compute delta and normalize per model (same as positive)
    numModels = size(topFavData, 1);
    numSweeps = length(cutRange);
    processedData = nan(numModels, numSweeps);
    for m = 1:numModels
        modelScores = topFavData(m, :);

        % Handle inf per model
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

        % Compute delta score relative to 0% change
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

        % Standardize by max absolute deviation
        maxDev = max(abs(modelScores));
        if maxDev > 0
            modelScores = modelScores / maxDev;
        end

        processedData(m, :) = modelScores;
    end

    % Compute mean and SEM across models (per sweep)
    meanSweepNeg = mean(processedData, 1, 'omitnan')';
    validPerSweep = sum(~isnan(processedData), 1)';
    semSweepNeg = std(processedData, 0, 1, 'omitnan')' ./ sqrt(max(1, validPerSweep));

    % Linear fit on mean
    p_fav_mdl = fitlm(cutRange, meanSweepNeg, 'linear');
    r2_adj = p_fav_mdl.Rsquared.Adjusted;
    p_val = p_fav_mdl.Coefficients.pValue(2).*10;  % Slope p-value
    r_pearson = corr(cutRange', meanSweepNeg);  % Keep for title if preferred, or use sqrt(mdl.Rsquared.Ordinary)

    % Prediction for full range (replaces polyval)
    yFitFav = predict(p_fav_mdl, cutRange(:));

    % Plot
    plot(cutRange, meanSweepNeg, 'k-', 'LineWidth', 2); % Mean
    hold on;
    % ±SEM shading
    upper = meanSweepNeg + semSweepNeg;
    lower = meanSweepNeg - semSweepNeg;
    fill([cutRange fliplr(cutRange)], [upper; flipud(lower)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(cutRange, yFitFav, 'r--', 'LineWidth', 2); % Linear fit
    xlabel('Sweep % Change', 'FontName', 'Arial', 'FontSize', 12);
    ylabel('Normalized Delta Score', 'FontName', 'Arial', 'FontSize', 12);
    title(sprintf('%s (r = %.2f)', topFavSubsys, r_pearson), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
    %legend('Mean Delta', '±SEM', 'Linear Fit', 'Location', 'best');
    grid on; box on; %xlim([-0.5 0.5]);
    ylim([-0.3 0.6])
end
% Save after all negative plots
saveas(fig6, fullfile('figures', 'Fig7_FavRegression.png'));
savefig(fig6, fullfile('figures', 'Fig7_FavRegression.fig'));
fprintf('Saved Figure 6 (fav regressions)\n');

%% Save Results (MODIFIED: Include allSweepScores)
fprintf('\nSaving analysis results...\n');
save(fullfile('out', 'driver7_results.mat'), 'subsystemStats', 'subsystemCorr', ...
    'predictorNames', 'allImportance', 'allSweepScores', 'sweepRange', ...
    'topPosSubsys', 'topFavSubsys', '-v7.3');
fprintf('\nDriver 7 completed successfully!\n');
fprintf('Generated 5 figures:\n');
fprintf(' - figures/Fig7_SubsystemImportance.png\n');
fprintf(' - figures/Fig7_SubsystemPositiveCorr.png\n');
fprintf(' - figures/Fig7_SubsystemNegativeCorr.png\n');
fprintf(' - figures/Fig7_TopPosRegression.png\n');
fprintf(' - figures/Fig7_TopNegRegression.png\n');

%% House-keeping
close all;

%% LOCAL HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makeXYtable - Convert statistics to ML-ready feature matrix and labels
function [X, y] = makeXYtable(UseStats, drugMap, regresTable)
% Initialize data table from median difference statistics
dataTable = table();
pathTableT_rowNames = string(UseStats.MedianDiff.Properties.VariableNames);
pathTableT_varNames = string(UseStats.MedianDiff.Properties.RowNames);
% Transpose and assign proper variable names
dataTable{:,:} = UseStats.MedianDiff{:,:}';
dataTable.Properties.VariableNames = pathTableT_varNames;
dataTable.Properties.RowNames = pathTableT_rowNames;
% Prepare target labels
regresTable.Drug = upper(string(regresTable.Drug));
dataTable.y = zeros(height(dataTable), 1);
% Map three-letter drug codes to full drug names
allThreeLet = extractBefore(pathTableT_rowNames, 4)';
for a = 1:height(regresTable)
    threeLet = drugMap.ThreeLet(drugMap.Drug == regresTable.Drug(a));
    drugSpots = contains(allThreeLet, string(threeLet));
    % Ensure unique mapping
    assert(isscalar(unique(allThreeLet(drugSpots))) || ...
        isempty(unique(allThreeLet(drugSpots))), ...
        'Multiple mappings found for drug %s', regresTable.Drug(a));
    dataTable.y(drugSpots) = regresTable.sigTF(a);
end
% Remove unlabeled samples
dataTable(dataTable.y == 0, :) = [];
% Separate features and labels
X = dataTable(:, 1:end-1);
Xtemp = X{:,:};
X{:,:} = fillmissing(Xtemp, 'constant', 0); % Fill missing values with zero
y = dataTable(:, end);
y.y(y.y == -1) = 0; % Convert -1 labels to 0 for binary classification
end

%% cleanSubsystemNames - Remove common suffixes and clean up subsystem names
function cleaned = cleanSubsystemNames(subsystemNames)
% Clean up subsystem names for better display
cleaned = string(subsystemNames);
% Common words to remove or shorten
cleaned = regexprep(cleaned, ' metabolism', '');
cleaned = regexprep(cleaned, ' biosynthesis', ' synth.');
cleaned = regexprep(cleaned, 'Oxidative phosphorylation', 'Ox. phos.');
cleaned = regexprep(cleaned, 'Tricarboxylic acid cycle', 'TCA cycle');
cleaned = regexprep(cleaned, 'Fatty acid', 'FA');
% Capitalize first letter
for i = 1:length(cleaned)
    if strlength(cleaned(i)) > 0
        cleaned(i) = upper(extractBefore(cleaned(i), 2)) + extractAfter(cleaned(i), 1);
    end
end
end
