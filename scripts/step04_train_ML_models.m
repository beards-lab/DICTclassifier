%%  step04_train_ML_models

%  Trains RusBoost ensemble classifiers for cardiotoxicity prediction using
%  both reaction-flux and RNA-seq features. Runs 100 iterations x 4-fold
%  cross-validation (400 total models). Also evaluates comprehensive
%  metrics (accuracy, F1, AUC, AUPRC) and saves ROC/PRC curve objects.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (fitcensemble, templateTree, crossval,
%                                    perfcurve, predict, rocmetrics)
%
%  Inputs:
%    - out/ROR_results.mat           — regresTable (cardiotoxicity labels)
%    - out/deltaFlux_stats.mat       — RxnStatsMedMed, RnaStatsMedMed
%    - data/drugMap.xlsx             — Drug name mapping table
%
%  Outputs:
%    - out/DICT_L0X/c1_<date>/Rxn_models_<a>.mat   — Per-iteration Rxn
%                                                      models (a = 1..100)
%    - out/DICT_L0X/c1_<date>/Rna_models_<a>.mat   — Per-iteration Rna
%                                                      models
%    - out/DICT_L0X/c1_<date>/CVpartition_<a>.mat   — CV partition per iter
%    - out/DICT_L0X/c1_<date>/Acca_<a>.mat, F1a_<a>.mat, AUCa_<a>.mat,
%                                    AUPRCa_<a>.mat — Per-iteration metrics
%    - out/trainedModels_comprehensive.mat           — Aggregated metrics,
%                                                      ROC/PRC objects
%    - out/DICT_L0X/c1_<date>/Plot_<a>.png          — Swarmer diagnostic
%                                                      plots per iteration

clear; clc; close all;

%% Load previous stage outputs
fprintf('Loading ROR results and flux statistics...\n');
load(fullfile('out','ROR_results.mat'), 'regresTable');
load(fullfile('out','deltaFlux_stats.mat'), 'RxnStatsMedMed', 'RnaStatsMedMed');
drugMap = readtable(fullfile("data", "drugMap.xlsx"));

%% Parameters
expName = "DICT_L0X";
currentDateString = "10_6_30hr";
saveVec = ["Acca", "F1a", "AUCa", "AUPRCa", "models", "CVpartition", "iRange"];
duoIdx = 1:5; % Indices for metrics to save separately for Rxn and Rna
iRange = 1:100; % Number of iterations for robust evaluation
cMax = 4; % Number of cross-validation folds

%% Prepare drug mapping and data integration
fprintf('Preparing drug mapping and feature matrices...\n');
drugMap.Drug = string(drugMap.Drug);
drugMap.ThreeLet = string(drugMap.ThreeLet);
regresTable.Drug = string(upper(regresTable.Drug));

%% Create feature matrices and target labels
[Rxn_X, Rxn_y] = makeXYtable(RxnStatsMedMed, drugMap, regresTable);
[Rna_X, Rna_y] = makeXYtable(RnaStatsMedMed, drugMap, regresTable);

% Verify target consistency between feature sets
assert(all(Rxn_y{:,:} == Rna_y{:,:}), 'Target labels must match between Rxn and Rna datasets');
y = Rxn_y;
clear Rna_y Rxn_y

%% Create output directories
if ~exist('out', 'dir'), mkdir('out'); end
if ~exist(fullfile("out", expName, "c1_" + currentDateString), 'dir')
    mkdir(fullfile("out", expName, "c1_" + currentDateString));
end

%% Main training loop with cross-validation
fprintf('Starting ML model training with %d iterations...\n', length(iRange));
close all force;

parfor a = iRange
    rng(a); % Set seed for reproducibility

    fprintf('Iteration %d/%d\n', a, max(iRange));

    % Initialize model storage
    Rxn_models = cell(1, cMax);
    Rna_models = cell(1, cMax);

    % Create stratified CV partition
    CVpartition = cvpartition(y{:,:}, 'KFold', cMax, "Stratify", true);

    % Initialize metric arrays for this iteration
    tRxn_Acca = zeros(1, cMax);
    tRxn_F1a = zeros(1, cMax);
    tRxn_AUCa = zeros(1, cMax);
    tRxn_AUPRCa = zeros(1, cMax);

    tRna_Acca = zeros(1, cMax);
    tRna_F1a = zeros(1, cMax);
    tRna_AUCa = zeros(1, cMax);
    tRna_AUPRCa = zeros(1, cMax);

    % K-fold cross-validation
    for c = 1:cMax
        % Split training data
        Rxn_X_train = Rxn_X(training(CVpartition, c), :);
        Rna_X_train = Rna_X(training(CVpartition, c), :);
        y_train = y(training(CVpartition, c), :);

        % Train ensemble models with hyperparameter optimization
        Rxn_models{1,c} = fitcOPTYS(Rxn_X_train, y_train);
        Rna_models{1,c} = fitcOPTYS(Rna_X_train, y_train);

        % Split testing data
        Rxn_X_test = Rxn_X(test(CVpartition, c), :);
        Rna_X_test = Rna_X(test(CVpartition, c), :);
        y_test = y(test(CVpartition, c), :);

        % Evaluate models
        Rxn_outStruct = testcOPTYS(Rxn_models{1,c}, Rxn_X_test, y_test);
        Rna_outStruct = testcOPTYS(Rna_models{1,c}, Rna_X_test, y_test);

        % Store metrics
        tRxn_Acca(1,c) = Rxn_outStruct.Acca;
        tRxn_F1a(1,c) = Rxn_outStruct.F1a;
        tRxn_AUCa(1,c) = Rxn_outStruct.AUCa;
        tRxn_AUPRCa(1,c) = Rxn_outStruct.AUPRCa;

        tRna_Acca(1,c) = Rna_outStruct.Acca;
        tRna_F1a(1,c) = Rna_outStruct.F1a;
        tRna_AUCa(1,c) = Rna_outStruct.AUCa;
        tRna_AUPRCa(1,c) = Rna_outStruct.AUPRCa;
    end

    % Save iteration results and generate plots
    saveTempMats(a, expName, currentDateString, saveVec, duoIdx, ...
        tRxn_Acca, tRxn_F1a, tRxn_AUCa, tRxn_AUPRCa, ...
        tRna_Acca, tRna_F1a, tRna_AUCa, tRna_AUPRCa, ...
        CVpartition, iRange, Rxn_models, Rna_models);

    plotSwarmerLoader(a, expName, currentDateString, saveVec);
end
%% Define comprehensive metrics and initialize storage arrays
% This defines the set of 22 metrics to compute for all CV iterations/folds
metrics = {'Accuracy', 'macroRecall', 'macroPrecision', 'macroSensitivity', ...
    'macroSpecificity', 'macroF1', 'macroAUC', 'macroAUPRC', ...
    'class1Recall', 'class1Precision', 'class1Sensitivity', 'class1Specificity', ...
    'class1F1', 'class1AUC', 'class1AUPRC', ...
    'class2Recall', 'class2Precision', 'class2Sensitivity', 'class2Specificity', ...
    'class2F1', 'class2AUC', 'class2AUPRC'};

Rxn = struct();
Rna = struct();
for i = 1:length(metrics)
    Rxn.(metrics{i}) = zeros(length(iRange), cMax);
    Rna.(metrics{i}) = zeros(length(iRange), cMax);
end
RxnROCObj = cell(length(iRange), cMax);
RnaROCObj = cell(length(iRange), cMax);
RxnPRCObj = cell(length(iRange), cMax);
RnaPRCObj = cell(length(iRange), cMax);
RxnROCObjStruct = cell(length(iRange), cMax);
RnaROCObjStruct = cell(length(iRange), cMax);
RxnPRCObjStruct = cell(length(iRange), cMax);
RnaPRCObjStruct = cell(length(iRange), cMax);

%% Loop over iRange and folds to calculate metrics for all saved models
basePath = fullfile('out', expName, ['c1' currentDateString]);
fprintf('Calculating comprehensive metrics for all models...\n');
for idx = 1:length(iRange)
    a = iRange(idx);
    rng(a);
    fprintf('Processing iteration %d/%d\n', a, max(iRange));
    % Load partition and models for this iteration
    load(fullfile(basePath, sprintf('CVpartition_%d.mat', a)), 'CVpartition');
    load(fullfile(basePath, sprintf('Rxn_models_%d.mat', a)), 'tRxnmodels');
    load(fullfile(basePath, sprintf('Rna_models_%d.mat', a)), 'tRnamodels');
    for c = 1:cMax
        % Prepare test data for fold
        RxnXtest = RxnX(test(CVpartition, c), :);
        RnaXtest = RnaX(test(CVpartition, c), :);
        ytest = y(test(CVpartition, c), :);
        % Compute metrics for Rxn classifier
        [allMetrics, RxnROCObj{idx, c}, RxnPRCObj{idx, c}, RxnROCObjStruct{idx, c}, RxnPRCObjStruct{idx, c}] = ...
            calculateComprehensiveMetrics(tRxnmodels{1,c}, RxnXtest, ytest);
        for i = 1:length(metrics)
            Rxn.(metrics{i})(idx, c) = allMetrics.(metrics{i});
        end
        % Compute metrics for Rna classifier
        [allMetrics, RnaROCObj{idx, c}, RnaPRCObj{idx, c}, RnaROCObjStruct{idx, c}, RnaPRCObjStruct{idx, c}] = ...
            calculateComprehensiveMetrics(tRnamodels{1,c}, RnaXtest, ytest);
        for i = 1:length(metrics)
            Rna.(metrics{i})(idx, c) = allMetrics.(metrics{i});
        end
    end
end

%% Save comprehensive metrics, ROC/PRC objects
save(fullfile('out', 'trainedModels_comprehensive.mat'), ...
    'metrics','Rxn','Rna','RxnROCObj','RnaROCObj','RxnPRCObj','RnaPRCObj','RxnROCObjStruct','RnaROCObjStruct','RxnPRCObjStruct','RnaPRCObjStruct','-v7.3');
fprintf('Saved comprehensive metrics to trainedModels_comprehensive.mat\n');

%% House-keeping
clear; clc;

%% LOCAL HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fitcOPTYS - Train ensemble classifier with Bayesian optimization
function [Mdl] = fitcOPTYS(X_train, y_train)

X_train = X_train{:,:};
y_train = y_train{:,:};

% Configure hyperparameter optimization
hpoOptions = hyperparameterOptimizationOptions(...
    'UseParallel', false, ...
    'ShowPlots', false, ...
    'Verbose', 0, ...
    'MaxObjectiveEvaluations', 100, ...
    'MaxTime', Inf, ...
    'Repartition', true, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Optimizer', 'bayesopt');

% Train ensemble with optimized hyperparameters
[Mdl] = fitcensemble(X_train, y_train, 'Learners', 'tree', ...
    'OptimizeHyperparameters', {'Method','NumLearningCycles','LearnRate', ...
    'MaxNumSplits', 'MinLeafSize', 'NumVariablesToSample', 'SplitCriterion'}, ...
    'HyperparameterOptimizationOptions', hpoOptions);

end

%% testcOPTYS - Evaluate model performance with comprehensive metrics
function outStruct = testcOPTYS(Mdl, X_test, y_test)

X_test = X_test{:,:};
y_test = y_test{:,:};

% Make predictions
[yPred_test, yScore_test] = predict(Mdl, X_test);

% Calculate ROC and PR metrics
Testing_rocObj = rocmetrics(y_test, yScore_test(:,2), Mdl.ClassNames(2), ...
    "AdditionalMetrics", "all");

% Confusion matrix analysis
testConfMat = confusionmat(y_test, yPred_test);
TotPos = sum(y_test);
TP = testConfMat(2,2);
TN = testConfMat(1,1);
FN = testConfMat(2,1);
FP = testConfMat(1,2);
numSamp = sum(sum(testConfMat));
TotNeg = numSamp - TotPos;

% Handle edge cases
if TotPos == 0
    fprintf('Warning: No positive samples in test set\n');
elseif TotNeg == 0
    fprintf('Warning: No negative samples in test set\n');
end

% Calculate F1 score with error handling
try
    F1a = (TP)/(TP + 1/2*(FP + FN));
    assert(~isnan(F1a));
catch
    F1a = 0;
    fprintf('Warning: F1 score calculation failed\n');
end

% Calculate performance metrics
Acca = (TP + TN)/numSamp;
AUCa = auc(Testing_rocObj, "roc");
AUPRCa = auc(Testing_rocObj, "pr");

% Package results
outStruct.F1a = F1a;
outStruct.Acca = Acca;
outStruct.AUCa = AUCa;
outStruct.AUPRCa = AUPRCa;

end

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

%% saveTempMats - Save iteration results to individual files
function saveTempMats(a, expName, currentDateString, saveVec, duoIdx, ...
    tRxn_Acca, tRxn_F1a, tRxn_AUCa, tRxn_AUPRCa, ...
    tRna_Acca, tRna_F1a, tRna_AUCa, tRna_AUPRCa, ...
    CVpartition, iRange, tRxn_models, tRna_models)

% Save metrics and models for each save variable
for sN = 1:length(saveVec)
    if any(sN == duoIdx)
        % Save separately for Rxn and Rna features
        save(fullfile("out", expName, "c1_" + currentDateString, ...
            "Rna_" + saveVec(sN) + "_" + a + ".mat"), ...
            "tRna_" + saveVec(sN));
        save(fullfile("out", expName, "c1_" + currentDateString, ...
            "Rxn_" + saveVec(sN) + "_" + a + ".mat"), ...
            "tRxn_" + saveVec(sN));
    else
        % Save common variables
        save(fullfile("out", expName, "c1_" + currentDateString, ...
            saveVec(sN) + "_" + a + ".mat"), saveVec(sN));
    end
end

end

%% plotSwarmerLoader - Generate swarm plots of performance metrics
function plotSwarmerLoader(a, expName, currentDateString, saveVec)

basePath = fullfile("out", expName, "c1_" + currentDateString);

% Initialize metric storage structures
metrics = {'Acca', 'F1a', 'AUCa', 'AUPRCa'};
Rxn = struct();
Rna = struct();

for m = metrics
    Rxn.(m{1}) = [];
    Rna.(m{1}) = [];
end

% Load and aggregate performance data
for metric = metrics
    % Load Rxn metric files
    rxFiles = dir(fullfile(basePath, sprintf('Rxn_%s_*.mat', metric{1})));
    for f = rxFiles'
        try
            loaded = load(fullfile(f.folder, f.name));
            Rxn.(metric{1}) = [Rxn.(metric{1}); loaded.(['tRxn_' metric{1}])];
        catch
            warning('Error loading %s', f.name);
        end
    end

    % Load Rna metric files
    rnFiles = dir(fullfile(basePath, sprintf('Rna_%s_*.mat', metric{1})));
    for f = rnFiles'
        try
            loaded = load(fullfile(f.folder, f.name));
            Rna.(metric{1}) = [Rna.(metric{1}); loaded.(['tRna_' metric{1}])];
        catch
            warning('Error loading %s', f.name);
        end
    end
end

% Generate performance comparison plots
close all force;
figure('Visible', 'on');
tiledlayout(2, 1, "TileSpacing", "tight", "Padding", "tight");

% Reaction-based features plot
nexttile();
RxnAcca = Rxn.Acca(:);
RxnF1a = Rxn.F1a(:);
RxnAUCa = Rxn.AUCa(:);
RxnAUPRCa = Rxn.AUPRCa(:);
numRxnAcca = length(RxnAcca);
numRxnF1a = length(RxnF1a);
numRxnAUCa = length(RxnAUCa);
numRxnAUPRCa = length(RxnAUPRCa);
numReadFiles = min([numRxnAcca, numRxnF1a, numRxnAUCa, numRxnAUPRCa]);
swarmData = [RxnAcca(1:numReadFiles), RxnF1a(1:numReadFiles), RxnAUCa(1:numReadFiles), RxnAUPRCa(1:numReadFiles)];
swarmchart(reshape(repmat(1:4, size(swarmData,1), 1), [], 1), swarmData(:));
ylim([0 1.01]);
xticks(1:4);
title('Reaction-based Features');
xticklabels([string("Accuracy " + sprintf('%.3f', mean(Rxn.Acca, 'all'))), ...
    string("F1 " + sprintf('%.3f', mean(Rxn.F1a, 'all'))), ...
    string("AUROC " + sprintf('%.3f', mean(Rxn.AUCa, 'all'))), ...
    string("AUPRC " + sprintf('%.3f', mean(Rxn.AUPRCa, 'all')))]);
box off;

% RNA-based features plot
nexttile();
RnaAcca = Rna.Acca(:);
RnaF1a = Rna.F1a(:);
RnaAUCa = Rna.AUCa(:);
RnaAUPRCa = Rna.AUPRCa(:);
numRnaAcca = length(RnaAcca);
numRnaF1a = length(RnaF1a);
numRnaAUCa = length(RnaAUCa);
numRnaAUPRCa = length(RnaAUPRCa);
numReadFiles = min([numRnaAcca, numRnaF1a, numRnaAUCa, numRnaAUPRCa]);
swarmData = [RnaAcca(1:numReadFiles), RnaF1a(1:numReadFiles), RnaAUCa(1:numReadFiles), RnaAUPRCa(1:numReadFiles)];
swarmchart(reshape(repmat(1:4, size(swarmData,1), 1), [], 1), swarmData(:));
ylim([0 1.01]);
xticks(1:4);
title('RNA-based Features');
xticklabels([string("Accuracy " + sprintf('%.3f', mean(Rna.Acca, 'all'))), ...
    string("F1 " + sprintf('%.3f', mean(Rna.F1a, 'all'))), ...
    string("AUROC " + sprintf('%.3f', mean(Rna.AUCa, 'all'))), ...
    string("AUPRC " + sprintf('%.3f', mean(Rna.AUPRCa, 'all')))]);
box off;

drawnow;
saveas(gcf, fullfile(basePath, sprintf('Plot_%d.png', a)));

end
%% calculateComprehensiveMetrics
function [allMetrics, rocObj, prcObj, rocObjStruct, prcObjStruct] = calculateComprehensiveMetrics(model, Xtest, ytest)
% Calculate all 22 metrics: macro, class-1, and class-2
Xtest = Xtest{:,:};
ytest = ytest{:,:};
[yPredtest, yScoretest] = predict(model, Xtest);
allMetrics = struct();
numClasses = 2;
numSamp = length(ytest);
classRecalls = zeros(numClasses, 1);
classPrecisions = zeros(numClasses, 1);
classSpecificities = zeros(numClasses, 1);
classF1s = zeros(numClasses, 1);
for c = 1:numClasses
    if c == 1
        TNtemp = sum((ytest == 0) & (yPredtest == 0));
        FPtemp = sum((ytest == 0) & (yPredtest == 1));
        FNtemp = sum((ytest == 1) & (yPredtest== 0));
        TPtemp = sum((ytest == 1) & (yPredtest == 1));
        TP = TNtemp; FP = FNtemp; FN = FPtemp; TN = TPtemp;
    elseif c == 2
        TN = sum((ytest == 0) & (yPredtest == 0));
        FP = sum((ytest == 0) & (yPredtest == 1));
        FN = sum((ytest == 1) & (yPredtest== 0));
        TP = sum((ytest == 1) & (yPredtest == 1));
    end
    try classRecalls(c) = TP / (TP + FN); catch, classRecalls(c) = 0; end
    if isnan(classRecalls(c)), classRecalls(c) = 0; end
    try classPrecisions(c) = TP / (TP + FP); catch, classPrecisions(c) = 0; end
    if isnan(classPrecisions(c)), classPrecisions(c) = 0; end
    try classSpecificities(c) = TN / (TN + FP); catch, classSpecificities(c) = 0; end
    if isnan(classSpecificities(c)), classSpecificities(c) = 0; end
    try classF1s(c) = TP./(TP + 1./2.*(FP + FN)); catch, classF1s(c) = 0; end
    if isnan(classF1s(c)), classF1s(c) = 0; end
    TestingrocObj = rocmetrics(ytest, yScoretest(:,c), model.ClassNames(c), "AdditionalMetrics","all");
    classAUCs = auc(TestingrocObj, "roc");
    classAUPRC = auc(TestingrocObj, "pr");
    if c == 1
        allMetrics.class1AUC = classAUCs;
        allMetrics.class1AUPRC = classAUPRC;
    elseif c == 2
        allMetrics.class2AUC = classAUCs;
        allMetrics.class2AUPRC = classAUPRC;
    end
end
allMetrics.Accuracy = (TP + TN)./numSamp;
allMetrics.macroRecall = mean(classRecalls);
allMetrics.macroPrecision = mean(classPrecisions);
allMetrics.macroSensitivity = mean(classRecalls);
allMetrics.macroSpecificity = mean(classSpecificities);
allMetrics.macroF1 = mean(classF1s);
assert(all(model.ClassNames == [0;1]), "Class name switched");
try
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames);
    [~, ~, ~, allMetrics.macroAUC] = average(rocObj, "macro");
catch
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames);
    [~, ~, ~, allMetrics.macroAUC] = average(rocObj, "macro");
end
try
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames);
    [~, ~, ~, allMetrics.macroAUPRC] = average(rocObj, "macro", "recall", "precision");
catch
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames, 'NaNFlag', 'includenan');
    [~, ~, ~, allMetrics.macroAUPRC] = average(rocObj, "macro", "recall", "precision");
end
allMetrics.class1Recall = classRecalls(1);
allMetrics.class1Precision = classPrecisions(1);
allMetrics.class1Sensitivity = classRecalls(1);
allMetrics.class1Specificity = classSpecificities(1);
allMetrics.class1F1 = classF1s(1);
allMetrics.class2Recall = classRecalls(2);
allMetrics.class2Precision = classPrecisions(2);
allMetrics.class2Sensitivity = classRecalls(2);
allMetrics.class2Specificity = classSpecificities(2);
allMetrics.class2F1 = classF1s(2);
rocObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all");
try
    prcObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all");
catch
    prcObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all", 'NaNFlag', 'includenan');
end
[rocObjStruct, prcObjStruct] = calculateCustomROCPRC(ytest, yScoretest);
end

%% calculateCustomROCPRC
function [rocObj, prcObj] = calculateCustomROCPRC(ytest, yScoretest)
% Calculate custom ROC/PRC for binary classification (returns struct curves)
minScore = min(yScoretest(:,2)); maxScore = max(yScoretest(:,2));
thresholds = linspace(minScore-0.01*(maxScore-minScore), maxScore+0.01*(maxScore-minScore), 1000);
numThresholds = length(thresholds);
fpr_class0 = zeros(numThresholds, 1); tpr_class0 = zeros(numThresholds, 1);
fpr_class1 = zeros(numThresholds, 1); tpr_class1 = zeros(numThresholds, 1);
precision_class0 = zeros(numThresholds, 1); recall_class0 = zeros(numThresholds, 1);
precision_class1 = zeros(numThresholds, 1); recall_class1 = zeros(numThresholds, 1);
for i = 1:numThresholds
    thresh = thresholds(i);
    yPred = (yScoretest(:, 2) >= thresh);
    TN = sum((ytest == 0) & (yPred == 0)); FP = sum((ytest == 0) & (yPred == 1));
    FN = sum((ytest == 1) & (yPred == 0)); TP = sum((ytest == 1) & (yPred == 1));
    TP_0 = TN; FP_0 = FN; FN_0 = FP; TN_0 = TP;
    fpr_class0(i) = FP_0 / max(FP_0 + TN_0, 1);
    tpr_class0(i) = TP_0 / max(TP_0 + FN_0, 1);
    precision_class0(i) = TP_0 / max(TP_0 + FP_0, 1);
    recall_class0(i) = tpr_class0(i);
    fpr_class1(i) = FP / max(FP + TN, 1);
    tpr_class1(i) = TP / max(TP + FN, 1);
    precision_class1(i)= TP / max(TP + FP, 1);
    recall_class1(i) = tpr_class1(i);
end
fpr_macro = (fpr_class0 + fpr_class1) / 2;
tpr_macro = (tpr_class0 + tpr_class1) / 2;
precision_macro = (precision_class0 + precision_class1) / 2;
recall_macro = (recall_class0 + recall_class1) / 2;
[fpr_class0, idx0] = sort(fpr_class0); tpr_class0 = tpr_class0(idx0);
[fpr_class1, idx1] = sort(fpr_class1); tpr_class1 = tpr_class1(idx1);
[fpr_macro, idxm] = sort(fpr_macro); tpr_macro = tpr_macro(idxm);
[recall_class0, idx0p] = sort(recall_class0, 'descend'); precision_class0 = precision_class0(idx0p);
[recall_class1, idx1p] = sort(recall_class1, 'descend'); precision_class1 = precision_class1(idx1p);
[recall_macro, idxmp] = sort(recall_macro, 'descend'); precision_macro = precision_macro(idxmp);
rocObj = struct();
rocObj.FPR_class0 = fpr_class0; rocObj.TPR_class0 = tpr_class0;
rocObj.FPR_class1 = fpr_class1; rocObj.TPR_class1 = tpr_class1;
rocObj.FPR_macro = fpr_macro; rocObj.TPR_macro = tpr_macro;
rocObj.thresholds = thresholds;
prcObj = struct();
prcObj.Precision_class0 = precision_class0; prcObj.Recall_class0 = recall_class0;
prcObj.Precision_class1 = precision_class1; prcObj.Recall_class1 = recall_class1;
prcObj.Precision_macro = precision_macro; prcObj.Recall_macro = recall_macro;
prcObj.thresholds = thresholds;
end