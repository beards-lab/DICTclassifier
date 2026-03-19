%%  step01_calculate_ROR

%  Calculates Reporting Odds Ratios (ROR) from FDA Adverse Event Reporting
%  System (FAERS) data to quantify cardiotoxicity signal for each drug.
%
%  Requirements:
%    - MATLAB (no special toolboxes)
%
%  Inputs:
%    - data/meddra_annotate.csv                        — MedDRA annotations
%    - data/aersMineExploreDataSet_176_all54.tsv.xlsx  — FAERS adverse event
%                                                        reports (176 PTs,
%                                                        54 drugs)
%
%  Outputs:
%    - out/ROR_results.mat           — regresTable with per-drug ROR values
%                                      and significance labels

clear; clc; close all;

%% Parameters
meddraFile = fullfile('data', 'meddra_annotate.csv');
aersFile   = fullfile('data', 'aersMineExploreDataSet_176_all54.tsv.xlsx');

%% Load MEDDRA Cardiac Terms (Level 4)
meddra       = readtable(meddraFile);
cardiacTerms = meddra.level4(ismember(meddra.level2, ["heart failures", ...
    "cardiac disorder signs and symptoms", "myocardial disorders"]));

%% Process FAERS Data
aersRaw = readtable(aersFile);

aersRaw(:, "amiodarone") = [];
aersRaw(:, "dimethylSulfoxide") = [];
aersRaw(:, "delavirdine") = [];
aersRaw(:, "diclofenac") = [];
aersRaw(:, "dobutamine") = [];
aersRaw(:, "estradiol") = [];
aersRaw(:, "flecainide") = [];
aersRaw(:, "isoprenaline") = [];
aersRaw(:, "mecasermin") = [];
aersRaw(:, "milrinone") = [];
aersRaw(:, "olmesartanMedoxomil") = [];
aersRaw(:, "phenylephrine") = [];
aersRaw(:, "prednisolone") = [];
aersRaw(:, "pioglitazone") = [];
aersRaw(:, "rosiglitazone") = [];
aersRaw(:, "saxagliptin") = [];
aersRaw(:, "tumorNecrosisFactorAlpha_tnf_alpha_Inhibitors") = [];
aersRaw(:, "verapamil") = [];

drugs = string(aersRaw.Properties.VariableNames(3:end));
aersRaw.Properties.VariableNames(3:end) = drugs;

% Reshape to long format
aersLong     = aersRaw;
aersLongNum  = aersLong{:, 2:end};
aersLongNum(isnan(aersLongNum)) = 0;
aersLong{:, 2:end} = aersLongNum;
aersLong(1,:) = [];
aersLong(:,2) = [];

% Filter cardiac events
isCardiac      = ismember(aersLong.AdverseEvents, cardiacTerms);
aersCardiac    = aersLong(isCardiac,:);
aersNonCardiac = aersLong(~isCardiac,:);
aersRaw(1,:)   = [];

%% Calculate ROR for Each Drug
results = table();
drugs   = string(aersLong.Properties.VariableNames(2:end));
for i = 1:numel(drugs)
    drug = drugs{i};
    fprintf('Processing %s (%d/%d)...\n', drug, i, numel(drugs));

    % Target drug counts
    drugCardiac = aersCardiac(:, drug);
    dycy        = sum(drugCardiac{:,:});
    drugNonCardiac = aersNonCardiac(:, drug);
    dycn        = sum(drugNonCardiac{:,:});

    % Other drugs counts
    otherCardiac    = sum(sum(aersRaw{isCardiac,"Total_AdverseEvents_Reports"})) - dycy;
    otherNonCardiac = sum(sum(aersRaw{~isCardiac,"Total_AdverseEvents_Reports"})) - dycn;

    % Calculate ROR components
    orr = (dycy / dycn) / (otherCardiac / otherNonCardiac);
    se  = sqrt(1/dycy + 1/dycn + 1/otherCardiac + 1/otherNonCardiac);
    ci  = exp(log(orr) + [1 -1]*1.96*se);

    % Store results
    results = [results; table({drug}, dycy, dycn, otherCardiac, otherNonCardiac, ...
        orr, se, ci(1), ci(2), 'VariableNames', ...
        {'Drug','dycy','dycn','dncy','dncn','ROR','SElogROR','CIhigh', ...
        'CIlow'})];
end

%% Post-Processing (debug figures)
figure();
tiledlayout(1,1, "TileSpacing","tight", "Padding","tight");

% Sort and filter TKIs
resultsPLOT = sortrows(results, 'ROR', 'ascend');
x = log(resultsPLOT.ROR);
y = height(resultsPLOT) - [1:height(resultsPLOT)] + 1;
ci = [log(resultsPLOT.CIlow) log(resultsPLOT.CIhigh)];
forestPlot(x, y, ci, resultsPLOT, "log ROR");

regresTable = results;
regresTable = regresTable(:, [1, 6,8,9]);
regresTable.sigTF = zeros(height(regresTable), 1);
regresTable.sigTF = double(results.ROR > median(results.ROR) & results.CIlow > median(results.ROR));
regresTable.sigTF(results.ROR < median(results.ROR) & results.CIhigh < median(results.ROR)) = -1;
regresTable(regresTable.sigTF == 0, :) = [];

%% Persist results for downstream drivers
if ~exist('out','dir'), mkdir out; end
save(fullfile('out','ROR_results.mat'),'regresTable','-v7.3');

%% House-keeping
clear;

%% LOCAL HELPER: forestPlot
function forestPlot(xA, yA, ciA, resultsPLOT, type)
% Visualization
hold on;
for i = 1:height(resultsPLOT)
    x = xA(i);
    y = yA(i);
    ci = ciA(i,:);
    line([ci(1) ci(2)], [y y], 'Color',[0.5 0.5 0.5], 'LineWidth',2);
    plot(x, y, 'o','MarkerFaceColor',[0.7 0 0], 'MarkerEdgeColor','k');
end
set(gca,'YDir','reverse','YTick',1:height(resultsPLOT),...
    'YTickLabel',resultsPLOT.Drug(end:-1:1),...
    'XGrid','on');
xlabel(type);
ylim([0 height(resultsPLOT)])
title('Cardiotoxicity Risk');
end
