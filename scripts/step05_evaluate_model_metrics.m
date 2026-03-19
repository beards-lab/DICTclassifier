%%  step05_evaluate_model_metrics

%  Recomputes comprehensive performance metrics (22 metrics including
%  accuracy, F1, AUC, AUPRC per class) across all 400 trained models.
%  Generates summary statistics and ROC/PRC curve objects.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (perfcurve, predict, rocmetrics)
%
%  Inputs:
%    - out/ROR_results.mat                           — regresTable (labels)
%    - out/deltaFlux_stats.mat                       — RxnStatsMedMed,
%                                                      RnaStatsMedMed
%    - data/drugMap.xlsx                             — Drug mapping table
%    - out/DICT_L0X/c1_<date>/CVpartition_<a>.mat   — CV partitions
%    - out/DICT_L0X/c1_<date>/Rxn_models_<a>.mat    — Trained Rxn models
%    - out/DICT_L0X/c1_<date>/Rna_models_<a>.mat    — Trained Rna models
%
%  Outputs:
%    - out/trainedModels_comprehensive.mat           — Rxn, Rna metric
%                                                      structs plus
%                                                      RxnROCObj, RnaPRCObj

clear; clc; close all;

%% Load previous stage outputs
fprintf('Loading ROR results and flux statistics...\n');
load(fullfile('out','ROR_results.mat'), 'regresTable');
load(fullfile('out','deltaFlux_stats.mat'), 'RxnStatsMedMed', 'RnaStatsMedMed');
drugMap = readtable(fullfile("data", "drugMap.xlsx"));

%% Parameters (must match driver4)
expName = "DICT_L0X";
currentDateString = "10_6_30hr";
iRange = 1:100;
cMax = 4;

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
basePath = fullfile("out", expName, "c1_"+ currentDateString);
fprintf('Calculating comprehensive metrics for all models...\n');

for idx = 1:length(iRange)
    a = iRange(idx);
    rng(a);
    fprintf('Processing iteration %d/%d\n', a, max(iRange));
    
    % Check if required files exist
    cvFile = fullfile(basePath, sprintf('CVpartition_%d.mat', a));
    rxnFile = fullfile(basePath, sprintf('Rxn_models_%d.mat', a));
    rnaFile = fullfile(basePath, sprintf('Rna_models_%d.mat', a));
    
    if ~exist(cvFile, 'file')
        fprintf('Warning: CVpartition file missing for iteration %d\n', a);
        continue;
    end
    if ~exist(rxnFile, 'file')
        fprintf('Warning: Rxn models file missing for iteration %d\n', a);
        continue;
    end
    if ~exist(rnaFile, 'file')
        fprintf('Warning: Rna models file missing for iteration %d\n', a);
        continue;
    end
    
    % Load partition and models for this iteration
    load(cvFile, 'CVpartition');
    rxnData = load(rxnFile);
    rnaData = load(rnaFile);
    
    % Handle variable name variations
    if isfield(rxnData, 'tRxnmodels')
        tRxnmodels = rxnData.tRxnmodels;
    elseif isfield(rxnData, 'tRxn_models')
        tRxnmodels = rxnData.tRxn_models;
    else
        fprintf('Warning: Cannot find Rxn models in file for iteration %d\n', a);
        continue;
    end
    
    if isfield(rnaData, 'tRnamodels')
        tRnamodels = rnaData.tRnamodels;
    elseif isfield(rnaData, 'tRna_models')
        tRnamodels = rnaData.tRna_models;
    else
        fprintf('Warning: Cannot find Rna models in file for iteration %d\n', a);
        continue;
    end
    
    for c = 1:cMax
        try
            % Prepare test data for fold
            RxnXtest = Rxn_X(test(CVpartition, c), :);
            RnaXtest = Rna_X(test(CVpartition, c), :);
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
            
        catch ME
            fprintf('Error processing iteration %d, fold %d: %s\n', a, c, ME.message);
            continue;
        end
    end
end

%% Save comprehensive metrics, ROC/PRC objects
save(fullfile('out', 'trainedModels_comprehensive.mat'), ...
    'metrics','Rxn','Rna','RxnROCObj','RnaROCObj','RxnPRCObj','RnaPRCObj',...
    'RxnROCObjStruct','RnaROCObjStruct','RxnPRCObjStruct','RnaPRCObjStruct',...
    'iRange','cMax','-v7.3');
fprintf('Saved comprehensive metrics to trainedModels_comprehensive.mat\n');

%% House-keeping
clear; clc;

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

%% calculateComprehensiveMetrics - Calculate all 22 metrics plus ROC/PRC objects
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

% Calculate class-specific metrics
for c = 1:numClasses
    if c == 1
        % For class 0: swap TP/TN and FP/FN
        TNtemp = sum((ytest == 0) & (yPredtest == 0));
        FPtemp = sum((ytest == 0) & (yPredtest == 1));
        FNtemp = sum((ytest == 1) & (yPredtest== 0));
        TPtemp = sum((ytest == 1) & (yPredtest == 1));
        TP = TNtemp; FP = FNtemp; FN = FPtemp; TN = TPtemp;
    elseif c == 2
        % For class 1: standard confusion matrix
        TN = sum((ytest == 0) & (yPredtest == 0));
        FP = sum((ytest == 0) & (yPredtest == 1));
        FN = sum((ytest == 1) & (yPredtest== 0));
        TP = sum((ytest == 1) & (yPredtest == 1));
    end
    
    % Calculate metrics with error handling
    try classRecalls(c) = TP / (TP + FN); catch, classRecalls(c) = 0; end
    if isnan(classRecalls(c)), classRecalls(c) = 0; end
    
    try classPrecisions(c) = TP / (TP + FP); catch, classPrecisions(c) = 0; end
    if isnan(classPrecisions(c)), classPrecisions(c) = 0; end
    
    try classSpecificities(c) = TN / (TN + FP); catch, classSpecificities(c) = 0; end
    if isnan(classSpecificities(c)), classSpecificities(c) = 0; end
    
    try classF1s(c) = TP./(TP + 1./2.*(FP + FN)); catch, classF1s(c) = 0; end
    if isnan(classF1s(c)), classF1s(c) = 0; end
    
    % Calculate AUC metrics for each class
    try
        TestingrocObj = rocmetrics(ytest, yScoretest(:,c), model.ClassNames(c), "AdditionalMetrics","all");
        classAUCs = auc(TestingrocObj, "roc");
        classAUPRC = auc(TestingrocObj, "pr");
    catch
        classAUCs = 0;
        classAUPRC = 0;
    end
    
    if c == 1
        allMetrics.class1AUC = classAUCs;
        allMetrics.class1AUPRC = classAUPRC;
    elseif c == 2
        allMetrics.class2AUC = classAUCs;
        allMetrics.class2AUPRC = classAUPRC;
    end
end

% Calculate overall metrics
TN = sum((ytest == 0) & (yPredtest == 0));
FP = sum((ytest == 0) & (yPredtest == 1));
FN = sum((ytest == 1) & (yPredtest== 0));
TP = sum((ytest == 1) & (yPredtest == 1));

allMetrics.Accuracy = (TP + TN)./numSamp;
allMetrics.macroRecall = mean(classRecalls);
allMetrics.macroPrecision = mean(classPrecisions);
allMetrics.macroSensitivity = mean(classRecalls);
allMetrics.macroSpecificity = mean(classSpecificities);
allMetrics.macroF1 = mean(classF1s);

% Calculate macro-averaged AUC metrics
assert(all(model.ClassNames == [0;1]), "Class names must be [0;1]");

try
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames);
    [~, ~, ~, allMetrics.macroAUC] = average(rocObj, "macro");
catch
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames);
    allMetrics.macroAUC = 0.5; % Default for failed calculations
end

try
    rocObjTemp = rocmetrics(ytest, yScoretest, model.ClassNames);
    [~, ~, ~, allMetrics.macroAUPRC] = average(rocObjTemp, "macro", "recall", "precision");
catch
    try
        rocObjTemp = rocmetrics(ytest, yScoretest, model.ClassNames, 'NaNFlag', 'includenan');
        [~, ~, ~, allMetrics.macroAUPRC] = average(rocObjTemp, "macro", "recall", "precision");
    catch
        allMetrics.macroAUPRC = 0.5; % Default for failed calculations
    end
end

% Assign class-specific metrics to output structure
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

% Calculate ROC and PRC objects
try
    rocObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all");
catch
    rocObj = [];
end

try
    prcObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all");
catch
    try
        prcObj = rocmetrics(ytest, yScoretest, model.ClassNames,"FixedMetric","tpr", "AdditionalMetrics","all", 'NaNFlag', 'includenan');
    catch
        prcObj = [];
    end
end

% Calculate custom ROC/PRC structures
try
    [rocObjStruct, prcObjStruct] = calculateCustomROCPRC(ytest, yScoretest);
catch
    rocObjStruct = [];
    prcObjStruct = [];
end
end

%% calculateCustomROCPRC - Calculate custom ROC/PRC curves
function [rocObj, prcObj] = calculateCustomROCPRC(ytest, yScoretest)
% Calculate custom ROC/PRC for binary classification (returns struct curves)
minScore = min(yScoretest(:,2)); 
maxScore = max(yScoretest(:,2));
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
    
    % Class 0 metrics (inverted)
    TP_0 = TN; FP_0 = FN; FN_0 = FP; TN_0 = TP;
    fpr_class0(i) = FP_0 / max(FP_0 + TN_0, 1);
    tpr_class0(i) = TP_0 / max(TP_0 + FN_0, 1);
    precision_class0(i) = TP_0 / max(TP_0 + FP_0, 1);
    recall_class0(i) = tpr_class0(i);
    
    % Class 1 metrics (standard)
    fpr_class1(i) = FP / max(FP + TN, 1);
    tpr_class1(i) = TP / max(TP + FN, 1);
    precision_class1(i)= TP / max(TP + FP, 1);
    recall_class1(i) = tpr_class1(i);
end

% Calculate macro averages
fpr_macro = (fpr_class0 + fpr_class1) / 2;
tpr_macro = (tpr_class0 + tpr_class1) / 2;
precision_macro = (precision_class0 + precision_class1) / 2;
recall_macro = (recall_class0 + recall_class1) / 2;

% Sort for proper curve construction
[fpr_class0, idx0] = sort(fpr_class0); tpr_class0 = tpr_class0(idx0);
[fpr_class1, idx1] = sort(fpr_class1); tpr_class1 = tpr_class1(idx1);
[fpr_macro, idxm] = sort(fpr_macro); tpr_macro = tpr_macro(idxm);

[recall_class0, idx0p] = sort(recall_class0, 'descend'); precision_class0 = precision_class0(idx0p);
[recall_class1, idx1p] = sort(recall_class1, 'descend'); precision_class1 = precision_class1(idx1p);
[recall_macro, idxmp] = sort(recall_macro, 'descend'); precision_macro = precision_macro(idxmp);

% Package ROC results
rocObj = struct();
rocObj.FPR_class0 = fpr_class0; rocObj.TPR_class0 = tpr_class0;
rocObj.FPR_class1 = fpr_class1; rocObj.TPR_class1 = tpr_class1;
rocObj.FPR_macro = fpr_macro; rocObj.TPR_macro = tpr_macro;
rocObj.thresholds = thresholds;

% Package PRC results
prcObj = struct();
prcObj.Precision_class0 = precision_class0; prcObj.Recall_class0 = recall_class0;
prcObj.Precision_class1 = precision_class1; prcObj.Recall_class1 = recall_class1;
prcObj.Precision_macro = precision_macro; prcObj.Recall_macro = recall_macro;
prcObj.thresholds = thresholds;
end