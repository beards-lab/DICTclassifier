%%  step12_dose_response_heatmaps

%  Loads pre-computed grid flux data from step11, runs ML inference across
%  all 400 classifiers for each grid point, normalizes scores, and
%  generates dose-response heatmaps, dot plots, and bubble plots with
%  FDR-corrected significance. Also creates an average anthracycline
%  delta-toxicity summary plot.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (predict, signtest)
%    - Bioinformatics Toolbox       (mafdr)
%    - Parallel Computing Toolbox   (parfor over models)
%
%  Inputs:
%    - outPf_GL/<Drug>_<Drug2>_combo/<combo>_d<x>_dd<y>.mat
%                                    — Grid flux data from step11
%    - out/iCardio_optimized.mat     — ModelOpt
%    - out/DICT_L0X/c1_<date>/Rxn_models_<a>.mat
%                                    — 400 trained Rxn classifiers
%
%  Outputs:
%    - outPf_GL/driver12_phase1_flux_and_scores.mat
%                                    — Cached flux + inference scores
%    - Figures displayed:            — Per-combo heatmaps, dot plots,
%                                      bubble plots, and average
%                                      anthracycline delta-toxicity plot

clear; clc; close all;

%% ========================== 2. CONFIGURATION ========================

% Drug combinations
DrugCombList = { ...
    {"DOX","OLM"}, ...
    {"IDA","OLM"}, ...
    {"DAU","OLM"}, ...
    {"EPI","OLM"} ...
    };

% Cell lines (must exist in grid files)
cellLines = {'MSN01','MSN02','MSN03','MSN04','MSN05','MSN06'};

numCombos   = numel(DrugCombList);
numCells    = numel(cellLines);

% Grid parameters
drug1_concs = 0:0.1:2;
drug2_concs = 0:0.1:2;
numD1      = numel(drug1_concs);
numD2      = numel(drug2_concs);
numGridPts = numD1 * numD2;

% ML configuration
expName           = "DICT_L0X";
currentDateString = "10_6_30hr";
iRange            = 1:100;  % Rxn_models_#.mat indices
cMax              = 4;      % 4 classifiers per file

% File/figure options
SAVE_FIG        = true;
SAVE_PHASE1_MAT = true;     % save flux + scores to .mat
phase1_mat_file = fullfile('outPf_GL','driver12_phase1_flux_and_scores.mat');

% Paths
if ~exist('outPf_GL','dir'),     mkdir('outPf_GL');     end
if ~exist('figures','dir'), mkdir('figures'); end

fprintf('\n========================================\n');
fprintf('DRIVER 12 v2: Grid Heatmaps with ML Risk\n');
fprintf('========================================\n\n');

%% ======================= 3. LOAD ML MODELS =========================

if exist('allModels','var') && ~isempty(allModels)
    % Already in memory from previous run in this session
    numModels = numel(allModels);
    fprintf('Using ML models already in memory (%d models)\n', numModels);
else
    fprintf('Loading ML models from disk...\n');
    basePath  = fullfile("out", expName, "c1_" + currentDateString);
    numModels = 400;
    allModels = cell(numModels,1);
    modelCount = 1;

    for idx = 1:numel(iRange)
        a = iRange(idx);
        modelFile = fullfile(basePath, sprintf('Rxn_models_%d.mat', a));
        if ~exist(modelFile,'file')
            continue;
        end
        S = load(modelFile, 'tRxn_models');
        tRxn_models = S.tRxn_models;
        for c = 1:cMax
            if modelCount <= numModels
                allModels{modelCount} = tRxn_models{1,c};
                modelCount = modelCount + 1;
            end
        end
    end

    allModels = allModels(~cellfun(@isempty, allModels));
    numModels = numel(allModels);
    fprintf('Loaded %d ML models\n', numModels);
end

%% =================== 4. LOAD MODEL STRUCTURE =======================

if exist('Model','var') && isstruct(Model) && isfield(Model,'rxns')
    numReactions = numel(Model.rxns);
    fprintf('Using iCardio model already in memory (%d reactions)\n', numReactions);
else
    fprintf('Loading iCardio model...\n');
    load(fullfile('out','iCardio_optimized.mat'), 'ModelOpt');
    Model        = ModelOpt;
    clear ModelOpt;
    numReactions = numel(Model.rxns);
    fprintf('Model loaded: %d reactions\n', numReactions);
end

%% ============== 5. INFERENCE LOOPS (PHASE 1) =======================

numCombos   = numel(DrugCombList);
numCells    = numel(cellLines);

% Precompute mapping grid index → linear index (1..121)
gridLinearIdx = reshape(1:numGridPts, [numD1, numD2]);

% Check if we already have flux + scores in memory or on disk
havePhase1InMemory = exist('fluxMat_allCombos','var') && ...
    exist('scores_allModels','var')   && ...
    numel(fluxMat_allCombos)  == numCombos && ...
    numel(scores_allModels)   == numCombos;

if havePhase1InMemory
    fprintf('Phase 1 data already in memory; skipping inference.\n');
elseif exist(phase1_mat_file,'file')
    fprintf('Phase 1 .mat found, loading from %s\n', phase1_mat_file);
    load(phase1_mat_file, 'fluxMat_allCombos','scores_allModels', ...
        'drug1_concs','drug2_concs','gridLinearIdx','DrugCombList');
    fprintf('Phase 1 data loaded from disk.\n');
else
    fprintf('Running Phase 1 inference (building fluxMat + scores)...\n');

    fluxMat_allCombos = cell(numCombos,1);   % each: [nRxn × 121]
    scores_allModels  = cell(numCombos,1);   % each: [numModels × 121]

    for comboIdx = 1:numCombos
        drug1     = char(DrugCombList{comboIdx}{1});
        drug2     = char(DrugCombList{comboIdx}{2});
        comboPair = sprintf('%s_vs_%s%s', drug1, drug1, drug2);
        comboLabel= sprintf('%s+%s', drug1, drug2);
        comboDir = sprintf('%s_%s_combo', drug1, drug2);

        fprintf('\n=== COMBO %d/%d: %s – building flux matrix ===\n', ...
            comboIdx, numCombos, comboLabel);

        % fluxMat: [numReactions × numGridPts]
        fluxMat = NaN(numReactions, numGridPts);

        for d1_idx = 1:numD1
            for d2_idx = 1:numD2
                d1_conc = drug1_concs(d1_idx);
                d2_conc = drug2_concs(d2_idx);
                linIdx  = gridLinearIdx(d1_idx, d2_idx);

                gridFile = fullfile('outPf_GL', comboDir, sprintf('%s_d%.1f_dd%.1f.mat', comboDir, d1_conc, d2_conc));
                if exist(gridFile, 'file')
                    S = load(gridFile);
                    flux = S.Conv_MedianDiff;
                    flux(ismissing(flux)) = 0;
                    flux(isnan(flux)) = 0;
                    col = zeros(numReactions, 1);
                    maxIdx = min(numel(flux), numReactions);
                    col(1:maxIdx) = flux(1:maxIdx);
                    fluxMat(:, linIdx) = col;
                end

            end
        end

        fluxMat_allCombos{comboIdx} = fluxMat;

        % ---------- inference over all models, vectorized over 121 points ----------
        fprintf('Running inference for %s ...\n', comboLabel);

        % X is [121 × numReactions]; fluxMat is [numReactions × 121]
        X = fluxMat';
        scoresMat = NaN(numModels, numGridPts);  % each row: model, col: grid

        parfor (m = 1:numModels)
            mdl = allModels{m};
            mdl.ScoreTransform = "doublelogit";
            [~, score] = predict(mdl, X);    % score: [121 × 2]
            scoresMat(m,:) = score(:,2)';    % toxic-class scores
        end

        scores_allModels{comboIdx} = scoresMat;

        fprintf('  Stored fluxMat (%d×%d) and scoresMat (%d×%d)\n', ...
            size(fluxMat,1), size(fluxMat,2), size(scoresMat,1), size(scoresMat,2));
    end

    if SAVE_PHASE1_MAT
        save(phase1_mat_file, ...
            'fluxMat_allCombos','scores_allModels', ...
            'drug1_concs','drug2_concs','gridLinearIdx','DrugCombList', ...
            '-v7.3');
        fprintf('\nPhase 1 data saved to %s\n', phase1_mat_file);
    end
end
%%
comboScoreMat = [];
for a = 1:numCombos
    comboScoreMat = [comboScoreMat, scores_allModels{a}];
end

comboScoreMat_zerod = comboScoreMat - comboScoreMat(:,1);
denom = max(comboScoreMat_zerod,[],2);
denom(denom == 0| isnan(denom)) = 1;
comboScoreMat_norm = comboScoreMat_zerod./denom;
%% ========= 6. NORMALIZATION + VISUALIZATION LOOPS (PHASE 2) ========

fprintf('\nRunning Phase 2 normalization + visualization...\n');

for comboIdx = 1:4
    drug1     = char(DrugCombList{comboIdx}{1});
    drug2     = char(DrugCombList{comboIdx}{2});
    comboPair = sprintf('%s_vs_%s%s', drug1, drug1, drug2);
    comboLabel= sprintf('%s+%s', drug1, drug2);

    fprintf('\n=== COMBO %d/%d: %s – normalizing & plotting ===\n', ...
        comboIdx, numCombos, comboLabel);

    scoresMat = scores_allModels{comboIdx};   % [numModels × 121]
    [numModels_local, numGridPts_local] = size(scoresMat);
    if numGridPts_local ~= numGridPts
        error('scoresMat size mismatch for combo %s', comboLabel);
    end

    splitDrugGrid = comboScoreMat_norm(:,((comboIdx-1)*numGridPts_local + 1): (comboIdx*numGridPts_local));
    splitDrugGrid_delta = splitDrugGrid.*0;
    pvalMat = [];
    for d1_idx = 1:numD1
        for d2_idx = 1:numD2
            linIdx = gridLinearIdx(d1_idx, d2_idx);
            zeroIdx = gridLinearIdx(d1_idx, 1);
            splitDrugGrid_delta(:, linIdx) = splitDrugGrid(:,linIdx) - splitDrugGrid(:, zeroIdx);
            pvalMat(d1_idx, d2_idx) = signtest(splitDrugGrid_delta(:,linIdx));
        end
    end

    pvals_vec = pvalMat(:);  % Flatten all tests
    pvalMat_fdr = mafdr(pvals_vec);
    pvalMat_fdr = reshape(pvalMat_fdr, size(pvalMat));  % Back to matrix
    pvalMat_fdr(pvalMat_fdr > 1) = 1;  % Clamp FDR <=1


    medianRiskPerGrid_delta = mean(splitDrugGrid_delta, 1);
    medianRiskPerGrid = mean(splitDrugGrid, 1);

    % reshape back to [numD1 × numD2] using same linear ordering
    heatmap_data = NaN(numD1, numD2);
    heatmap_data_delta = NaN(numD1, numD2);
    for d1_idx = 1:numD1
        for d2_idx = 1:numD2
            linIdx = gridLinearIdx(d1_idx, d2_idx);
            heatmap_data_delta(d1_idx, d2_idx) = medianRiskPerGrid_delta(linIdx);
            heatmap_data(d1_idx, d2_idx) = medianRiskPerGrid(linIdx);
        end
    end

    figure('Position', [10, 10, 1200, 1000]);
    tiledlayout('flow')
    nexttile
    [X, Y] = meshgrid(drug2_concs(1:11), flipud(drug1_concs(1:11)));
    pval_subset = pvalMat(1:11,1:11);  % Validate first
    pval_subset(pval_subset < 0 | pval_subset > 1 | ~isfinite(pval_subset)) = 1;  % Safe pvals
    pval_subset(1,:) = 1;
    sz = max(min(-log10(pval_subset), 20), 0);  % [0,20]; handles all cases
    subset_mask = (sz > 0) & (sz < Inf);  % Only significant/valid points
    subset_mask(:,1) = false;  % Exclude anchor
    clr = max(min(heatmap_data_delta(1:11,1:11), 0.01), -0.01);
    scatter(X(subset_mask), Y(subset_mask), 300 * (sz(subset_mask) / 20), ...
        clr(subset_mask), 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
    hold on
    colormap(redblue); clim([-0.05 0.05]);
    cb = colorbar; cb.Label.String = 'Mean ΔToxicity';
    xlabel(sprintf('Drug2 (%s) Conc. ', drug2)); ylabel(sprintf('Drug1 (%s) Conc. ', drug1));
    title(sprintf('%s + %s: OLM ΔToxicity Dot Plot', drug1, drug2));
    xlim([0 1.1]); ylim([0 1.1]); grid on; box on;
    hold off

    nexttile;
    [X, Y] = meshgrid(drug2_concs(1:11), flipud(drug1_concs(1:11)));
    pval_subset = pvalMat_fdr(1:11,1:11);
    pval_subset(pval_subset <= 0 | pval_subset > 0.05 | ~isfinite(pval_subset)) = 1;
    pval_subset(1,:) = 1;
    sz_raw = max(min(-log10(pval_subset), 20), 0);
    subset_mask = sz_raw > 0; subset_mask(:,1) = false;
    clr_raw = sign(heatmap_data_delta(1:11,1:11));  % -1/0/+1 discrete
    % Increased (red)
    pos_mask = subset_mask & (clr_raw > 0);
    bubblechart(X(pos_mask), Y(pos_mask), sz_raw(pos_mask), 'r', 'MarkerEdgeColor','k');
    hold on;
    % Decreased (blue)
    neg_mask = subset_mask & (clr_raw < 0);
    bubblechart(X(neg_mask), Y(neg_mask), sz_raw(neg_mask), 'b', 'MarkerEdgeColor','k');
    % Neutral (white, optional/small)
    neut_mask = subset_mask & (clr_raw == 0);
    bubblechart(X(neut_mask), Y(neut_mask), sz_raw(neut_mask)*0.5, 'w', 'MarkerEdgeColor','k');  % Smaller
    bubblesize(gca, [1 20]);
    bubblelegend('-log_{10}(P_{adj})', 'Location','northeastoutside', 'FontSize',10);
    legend('Increased Risk', 'Decreased Risk', 'Neutral Effect', ...
        'Location','southeastoutside', 'FontSize',10);
    xlabel(sprintf('Drug2 (%s) Conc. ', drug2)); ylabel(sprintf('Drug1 (%s) Conc. ', drug1));
    title(sprintf('%s + %s: OLM ΔToxicity Bubbles', drug1, drug2));
    xlim([0 1.1]); ylim([0 1.1]); grid on; box on; hold off;



    %  ------------- Visualization -------------
    fprintf('  Creating heatmap figure...\n');
    % figure('Units','centimeters','Position',[2,2,14,12],'Color','w');
    nexttile
    heatmap_data_delta(1,:) = 0;
    imagesc(drug2_concs(1:11), drug1_concs(1:11), heatmap_data_delta(1:11,1:11));
    hold on
    set(gca,'YDir','normal');
    colormap(redblue);   % 0=blue, 0.5=white, 1=red
    clim([-0.05, 0.05]);
    %
    cb = colorbar;
    cb.Label.String     = 'Mean \Delta Toxicity Risk Score';
    cb.Label.FontSize   = 11;
    cb.Label.FontWeight = 'bold';

    xlabel(sprintf('Drug2 (%s) Concentration ', drug2), ...
        'FontSize',11,'FontWeight','bold');
    ylabel(sprintf('Drug1 (%s) Concentration (normalized)', drug1), ...
        'FontSize',11,'FontWeight','bold');
    title(sprintf('%s + %s: ML-Predicted Toxicity Risk', drug1, drug2), ...
        'FontSize',12,'FontWeight','bold');

    xticks(drug2_concs);
    yticks(drug1_concs);
    set(gca,'FontSize',9);

    hold on;
    %contour(drug2_concs(1:11), drug1_concs(1:11), heatmap_data_delta(1:11, 1:11), [0.5,0.5], 'k--', 'LineWidth',1.5);
    hold off;

    nexttile
    imagesc(drug2_concs(1:11), drug1_concs(1:11), heatmap_data(1:11,1:11));
    hold on
    set(gca,'YDir','normal');
    colormap(redblue);   % 0=blue, 0.5=white, 1=red
    clim([0, 1]);
    %
    cb = colorbar;
    cb.Label.String     = 'Toxicity Risk Score';
    cb.Label.FontSize   = 11;
    cb.Label.FontWeight = 'bold';

    xlabel(sprintf('Drug2 (%s) Concentration ', drug2), ...
        'FontSize',11,'FontWeight','bold');
    ylabel(sprintf('Drug1 (%s) Concentration (normalized)', drug1), ...
        'FontSize',11,'FontWeight','bold');
    title(sprintf('%s + %s: ML-Predicted Toxicity Risk', drug1, drug2), ...
        'FontSize',12,'FontWeight','bold');

    xticks(drug2_concs);
    yticks(drug1_concs);
    set(gca,'FontSize',9);

end
%% ============== 6.5 AVERAGE ANTHRACYCLINE DELTA PLOT ===============
fprintf('\n=== Creating Average Anthracycline ΔToxicity Plot ===\n');

% Collect ALL delta values from all 4 anthracyclines
% Store as [numModels × numGridPts × numCombos]
all_delta_values = NaN(numModels, numGridPts, numCombos);

for comboIdx = 1:numCombos
    drug1 = char(DrugCombList{comboIdx}{1});
    drug2 = char(DrugCombList{comboIdx}{2});
    
    scoresMat = scores_allModels{comboIdx};
    splitDrugGrid = comboScoreMat_norm(:, ((comboIdx-1)*numGridPts + 1):(comboIdx*numGridPts));
    splitDrugGrid_delta = splitDrugGrid * 0;
    
    % Compute delta for this combo
    for d1_idx = 1:numD1
        for d2_idx = 1:numD2
            linIdx = gridLinearIdx(d1_idx, d2_idx);
            zeroIdx = gridLinearIdx(d1_idx, 1);
            splitDrugGrid_delta(:, linIdx) = splitDrugGrid(:, linIdx) - splitDrugGrid(:, zeroIdx);
        end
    end
    
    % Store all delta values for this combo
    all_delta_values(:, :, comboIdx) = splitDrugGrid_delta;
end

% Now compute p-values using ALL 400 delta values per grid location
pvalMat_combined = NaN(numD1, numD2);
effect_direction = NaN(numD1, numD2);

for d1_idx = 1:numD1
    for d2_idx = 1:numD2
        linIdx = gridLinearIdx(d1_idx, d2_idx);
        
        % Get all 400 delta values (400 models × 4 combos) for this grid point
        all_deltas_this_point = all_delta_values(:, linIdx, :); % [400 × 1 × 4]
        all_deltas_vec = all_deltas_this_point(:); % Flatten to [1600 × 1]
        
        % Signtest on all 400 values
        pvalMat_combined(d1_idx, d2_idx) = signtest(all_deltas_vec);
        
        % Mean of signs for direction
        effect_direction(d1_idx, d2_idx) = mean(sign(all_deltas_vec), 'omitnan');
    end
end
pvalMat_combined = pvalMat_combined(1:11, 1:11);
% FDR correction on the combined p-values
pvals_vec = pvalMat_combined(:);
pvalMat_fdr = mafdr(pvals_vec);
pvalMat_fdr = reshape(pvalMat_fdr, size(pvalMat_combined));
pvalMat_fdr(pvalMat_fdr > 1) = 1;

% Threshold effect to binary (blue/red) and mask non-significant
effect_sign = sign(effect_direction); % -1, 0, or +1
significant_mask = (pvalMat_fdr <= 0.05) & isfinite(pvalMat_fdr);
significant_mask(1, :) = false; % Exclude first row (reference)

% Create bubble plot
figure('Position', [100, 100, 800, 700]);
[X, Y] = meshgrid(drug2_concs(1:11), flipud(drug1_concs(1:11)));

% Prepare bubble size from combined p-value
sz_combined = max(min(-log10(pvalMat_fdr(1:11, 1:11)), 20), 0);
sz_combined(~significant_mask(1:11, 1:11)) = 0; % Zero out non-significant

% Plot increased risk (red)
pos_mask = significant_mask(1:11, 1:11) & (effect_sign(1:11, 1:11) > 0);
bubblechart(X(pos_mask), Y(pos_mask), sz_combined(pos_mask), 'r', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');
hold on;

% Plot decreased risk (blue)
neg_mask = significant_mask(1:11, 1:11) & (effect_sign(1:11, 1:11) < 0);
bubblechart(X(neg_mask), Y(neg_mask), sz_combined(neg_mask), 'b', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');

% Plot neutral (white) - optional, typically very small effect
neut_mask = significant_mask(1:11, 1:11) & (effect_sign(1:11, 1:11) == 0);
if any(neut_mask(:))
    bubblechart(X(neut_mask), Y(neut_mask), sz_combined(neut_mask)*0.5, 'w', 'MarkerEdgeColor', 'k');
end

bubblesize(gca, [1 20]);
bubblelegend('-log_{10}(P_{FDR})', 'Location', 'northeastoutside', 'FontSize', 10);
legend('Increased Risk', 'Decreased Risk', 'Location', 'southeastoutside', 'FontSize', 10);

xlabel('Drug2 (OLM) Conc.', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Anthracycline Conc.', 'FontSize', 11, 'FontWeight', 'bold');
title('Average Anthracycline + OLM: ΔToxicity', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 1.1]); ylim([0 1.1]); grid on; box on;
hold off;

fprintf('Average anthracycline ΔToxicity plot created.\n');

%% ======================= 7. PRINT DONE SECTION =====================

fprintf('\n========================================\n');
fprintf('DRIVER 12 v2 completed successfully!\n');
fprintf('Generated 4 heatmaps in figures/\n');
fprintf('========================================\n\n');

%% =========================== 8. HELPER FUNCTIONS ===================

function cmap = redblue
% Red-Blue diverging colormap: Blue (0) -> White (0.5) -> Red (1)
n = 256;
redValues   = [linspace(0,1,n/2),   linspace(1,1,n/2)];
greenValues = [linspace(0,1,n/2),   linspace(1,0,n/2)];
blueValues  = [linspace(1,1,n/2),   linspace(1,0,n/2)];
cmap = [redValues(:), greenValues(:), blueValues(:)];
end
