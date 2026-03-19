%%  step02_optimize_GEM

%  Performs thermodynamic optimization of the iCardio genome-scale metabolic
%  model (GEM) using irreversibility constraints. Identifies and enforces
%  reaction directionality.
%
%  Requirements:
%    - MATLAB
%    - Optimization Toolbox         (linprog, optimvar, optimproblem)
%
%  Inputs:
%    - models/objRxns.mat            — Objective reaction list
%    - models/HeartModel.mat         — Base iCardio heart model
%
%  Outputs:
%    - out/iCardio_optimized.mat     — ModelOpt (optimized model with S, N,
%                                      NNp, bounds, reversibility flags)
%                                      and boundTable

clear; clc; close all;

%% Parameters
allRxnFile    = fullfile('models', 'objRxns.mat');
modelFile     = fullfile('models', 'HeartModel.mat');
resetBounds   = 2;      % 1 = all [-1000 1000]; 2 = as loaded bounds
dispReactions = false;  % Set to true for debugging
minValue      = eps * 10^9;  % Minimum flux constraint

%% Load model and reaction table
fprintf('Loading model and reaction data...\n');
loaded_allRxnFile = load(allRxnFile);
load(modelFile);
loadedModel = heart_model_curation;
loadedModel.allRxn = loaded_allRxnFile.allRxn;

%% Add purine degradation reactions
fprintf('Adding purine degradation reactions...\n');
Model = addNewReaction(loadedModel, "Aden2Ino", "Adenosine2Ino", ...
    {"adenosine", "H2O"}, {"c", "c"}, {"inosine", "NH3"}, {"c", "c"}, ...
    [-1 -1 1 1], -1000, 1000, 0, "Purine Deg");

Model = addNewReaction(Model, "AMP2IMP", "AMP Deaminase", ...
    {"AMP", "H2O"}, {"c", "c"}, {"IMP", "NH3"}, {"c", "c"}, ...
    [-1 -1 1 1], 0, 1000, 0, "Purine Deg");

Model.allRxn = allRxn_reIndexer(Model);
%% Reset bounds and track initial reaction directionality
if resetBounds == 1
    Model.ub(:) = 1000;
    Model.lb(:) = -1000;
end

Model = cleanSfunc(Model);
Model.allRxn = allRxn_reIndexer(Model);
S = full(Model.S);

% Track which reactions get reversed (for pathway analysis later)
isReversed = Model.FlippedIrrev;

%% Set metabolic constraints for optimization
constraintRxns = ["Glucose Importation", "Oxygen Import", "Hypoxanthine Export", "H2O Import", "CO2 removal", ...
    "ATP Hydrolosis", "ATPsynthase"];

constrainIdxs = [
    Model.allRxn{constraintRxns{1}, "IDX"}, ...
    Model.allRxn{constraintRxns{2}, "IDX"}, ...
    Model.allRxn{constraintRxns{3}, "IDX"}, ...
    Model.allRxn{constraintRxns{4}, "IDX"}, ...
    Model.allRxn{constraintRxns{5}, "IDX"}, ...
    Model.allRxn{constraintRxns{6}, "IDX"}, ...
    Model.allRxn{constraintRxns{7}, "IDX"}];

fprintf('Starting epsilon optimization to determine thermodynamic feasibility...\n');
numIrrevOrig = sum((sign(Model.ub) == sign(Model.lb)) | ...
    (sign(Model.ub) == 0) | (sign(Model.lb) == 0));

%% Epsilon optimization loop
boundTable = table();
boundTable.Name = string(Model.rxnNames);
boundTable.Properties.RowNames = boundTable.Name;
boundTable.irrev = false(height(boundTable), 1);

noChange = false;
cycles = 0;

while ~noChange
    cycles = cycles + 1;
    fprintf('Optimization cycle %d...\n', cycles);

    boundTableOLD = boundTable;
    tempMax = zeros(1, width(S));
    tempMin = zeros(1, width(S));
    tempNames = strings(1, width(S));

    %%% Test each reaction for thermodynamic feasibility
    parfor a = 1:width(S)
        if mod(a, 100) == 0
            fprintf('  Testing reaction %d/%d\n', a, width(S));
        end

        %%% Apply metabolic constraints during optimization
        temp_Model = Model;
        temp_Model.lb(Model.allRxn{constraintRxns{1}, "IDX"}) = minValue;  % Glucose uptake
        temp_Model.lb(Model.allRxn{constraintRxns{2}, "IDX"}) = minValue;  % Oxygen uptake
        temp_Model.ub(Model.allRxn{constraintRxns{3}, "IDX"}) = -minValue; % Hypoxanthine Export
        temp_Model.ub(Model.allRxn{constraintRxns{5}, "IDX"}) = -minValue; % CO2 export
        temp_Model.lb(Model.allRxn{constraintRxns{6}, "IDX"}) = minValue;  % ATP Hydrolosis
        temp_Model.lb(Model.allRxn{constraintRxns{7}, "IDX"}) = minValue;  % ATP synthesis

        rxnName = temp_Model.rxns(a);
        if dispReactions
            fprintf('  %s\n', rxnName);
        end
        tempNames(a) = rxnName;

        % Calculate maximum feasible flux
        try
            if temp_Model.ub(a) > 0 && temp_Model.lb(a) >= 0
                tempMax(a) = temp_Model.ub(a);
            else
                max_tempFlux = GEMoptyFUNC_noName(temp_Model, S, a, "max", dispReactions);
                tempMax(a) = max_tempFlux(a);
            end
        catch
            fprintf('  Warning: Could not calculate max flux for rxn %d\n', a);
            tempMax(a) = NaN;
        end

        % Calculate minimum feasible flux
        try
            if temp_Model.lb(a) < 0 && temp_Model.ub(a) <= 0
                tempMin(a) = temp_Model.lb(a);
            else
                min_tempFlux = GEMoptyFUNC_noName(temp_Model, S, a, "min", dispReactions);
                tempMin(a) = min_tempFlux(a);
            end
        catch
            fprintf('  Warning: Could not calculate min flux for rxn %d\n', a);
            tempMin(a) = NaN;
        end
    end

    %%% Update bound table and determine irreversibility
    boundTable.Name = tempNames';
    boundTable.Max = tempMax';
    boundTable.Min = tempMin';

    % Reaction is irreversible if min and max have same sign
    boundTable.irrev = (sign(boundTable.Max) == sign(boundTable.Min)) | ...
        (sign(boundTable.Max) == 0) | (sign(boundTable.Min) == 0);

    %%% Check for convergence
    if all(boundTableOLD.irrev == boundTable.irrev)
        noChange = true;
        fprintf('Convergence achieved!\n');
    else
        % Update bounds and check for reactions to reverse
        Model = updateBounds(Model, boundTable);
    end
end

%% Final cleanup - remove constraints and set irreversible bounds
Model = updateBounds(Model, boundTable);
Model.lb(Model.allRxn{constraintRxns{1}, "IDX"}) = minValue;  % Glucose uptake
Model.lb(Model.allRxn{constraintRxns{2}, "IDX"}) = minValue;  % Oxygen uptake
Model.ub(Model.allRxn{constraintRxns{3}, "IDX"}) = -minValue; % Hypoxanthine Export
Model.ub(Model.allRxn{constraintRxns{5}, "IDX"}) = -minValue; % CO2 export
Model.lb(Model.allRxn{constraintRxns{6}, "IDX"}) = minValue;  % ATP Hydrolosis
Model.lb(Model.allRxn{constraintRxns{7}, "IDX"}) = minValue;  % ATP synthesis

%% Report optimization results
new_numIrrev = sum((sign(Model.ub) == sign(Model.lb)) | ...
    (sign(Model.ub) == 0) | (sign(Model.lb) == 0));

fprintf('\n=== OPTIMIZATION RESULTS ===\n');
fprintf('Cycles completed: %d\n', cycles);
fprintf('Irreversible reactions: %d → %d (net: +%d)\n', ...
    numIrrevOrig, new_numIrrev, new_numIrrev - numIrrevOrig);
fprintf('Reactions reversed: %d\n', sum(isReversed));

%% Store optimized model with reversal tracking
ModelOpt = Model;
ModelOpt.description = sprintf("Epsilon-optimized iCardio (%d cycles)", cycles);
ModelOpt.isReversed = isReversed;  % Track which rxns were reversed for pathway analysis

%% Persist results for downstream drivers
if ~exist('out', 'dir'), mkdir('out'); end
save(fullfile('out', 'iCardio_optimized.mat'), 'ModelOpt', 'boundTable', '-v7.3');

%% House-keeping
clear;

%% LOCAL HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cleanSfunc
function [ModelOut] = cleanSfunc(Model)

killNucEmptyAnd00 = true;
killAllImportersBut = true;

mets = convertStringsToChars(string(Model.mets));
metNames = string(Model.metNames);
killMet = false(height(mets),1);

rxns = string(Model.rxns);
rxnNames = string(Model.rxnNames);
killRxns = false(width(Model.S),1);
S = full(Model.S);

onlyBackwards = (Model.lb < 0 & Model.ub == 0) |  all(S ~= 1, 1)';

S_new = full(Model.S);
S_new(:,onlyBackwards) = S_new(:,onlyBackwards).*-1;


lb_new = Model.lb;
lb_new(onlyBackwards) = Model.ub(onlyBackwards).*-1;

ub_new = Model.ub;
ub_new(onlyBackwards) = Model.lb(onlyBackwards).*-1;

% Finds index of reactiosn that contain nuclear metabolites, "Bad" metabolites,
% metabolites not involved in reactions, reactions with no reactants or
% products, and reactions with lb = ub = 0

if killNucEmptyAnd00
    tbRemovedMets = [
        "DNA"
        "misc biomass metabolites"];%

    for a = 1:height(mets)
        tempEND = mets{a}(end-2:end);
        if all(tempEND == '[n]')
            killMet(a) = true;
        end
        if sum(abs(S_new(a,:))) == 0
            killMet(a) = true;
        end
    end
    for a = 1:height(tbRemovedMets)
        killMet(metNames == tbRemovedMets(a)) = 1;
    end
    for a = 1:width(S_new)
        tempS = S_new(:,a);
        tempS(tempS ~= 0) = 1;
        for b = 1:height(tempS)
            if tempS(b) && killMet(b)
                killRxns(a) = true;
            end
        end
        if sum(abs(S_new(:,a))) == 0
            killRxns(a) = true;
        end
    end
    killRxns(ub_new == lb_new) = 1;
end

% Removes all but __ importers
if killAllImportersBut
    keepPorterList = [
        "H2O"
        "O2"
        "CO2"
        "glucose"
        "palmitate"
        "L-lactate"
        "glycine"
        "glutamine"
        "aspartate"
        ];

    isImporter = zeros(width(S_new),1);
    for a = 1:width(S_new)
        tempS = S_new(:,a);
        if (sum(abs(tempS)) == 1) && (sum(tempS) == 1) && (ub_new(a) > 0)
            isImporter(a) = 1;
            tempM = metNames(logical(abs(tempS)));
        elseif (sum(abs(tempS)) == 1) && (sum(tempS) == -1) && (lb_new(a) < 0)
            isImporter(a) = -1;
            tempM = metNames(logical(abs(tempS)));
        else
            tempM = "er";
        end
        if ~ismember(tempM, keepPorterList)
            if isImporter(a) == 1
                ub_new(a) = 0;
            elseif isImporter(a) == -1
                lb_new(a) = 0;
            end
        end
    end
end

% Remove mets and rxns
rxns(killRxns) = [];
S_new(:,killRxns') = [];
S_new(killMet,:) = [];
lb_new(killRxns) = [];
ub_new(killRxns) = [];
Model.c(killRxns) = [];
mets = string(mets);
mets(killMet) = [];
Model.b(killMet) = [];
Model.rules(killRxns) = [];
Model.rxnGeneMat(killRxns,:) = [];
metNames(killMet) = [];
Model.subSystems(killRxns) = [];
Model.grRules(killRxns) = [];
rxnNames(killRxns) = [];
Model.rxnConfidenceScores(killRxns) = [];
Model.rxnECNumbers(killRxns) = [];
Model.rxnNotes(killRxns) = [];
Model.rxnReferences(killRxns) = [];
Model.metFormulas(killMet) = [];
Model.metCharges(killMet) = [];
Model.metSmiles(killMet) = [];
Model.metKEGGID(killMet) = [];
Model.metInChIString(killMet) = [];
Model.metPubChemID(killMet) = [];
Model.metChEBIID(killMet) = [];
Model.description = "Post 1st cleanup: " + string(Model.description);

ub_new(ub_new > 0) = 10000;
lb_new(lb_new < 0) = -10000;

% Save updated mets and rxns and bounds
ModelOut = Model;
ModelOut.S = sparse(S_new);
ModelOut.mets = mets;
ModelOut.metNames = metNames;
ModelOut.rxns = rxns;
ModelOut.rxnNames = rxnNames;
ModelOut.ub = ub_new;
ModelOut.lb = lb_new;
onlyBackwards(killRxns) = [];
try
    ModelOut.FlippedIrrev(onlyBackwards == 1) = 1;
catch
    ModelOut.FlippedIrrev = onlyBackwards;
end
end

%% updateBounds - Set reaction bounds based on feasible flux ranges
function Model = updateBounds(Model, boundTable)
% Create masks for bounds to update
Model.lb(boundTable.Min >= 0) = 0;
Model.ub(boundTable.Max <= 0) = 0;
end

%% addNewReaction - Add metabolic reaction to model
function updatedModel = addNewReaction(model, newRxnId, newRxnName, ...
    reactants, reactantComps, products, productComps, ...
    stoichiometry, lb, ub, c, subsystem)
updatedModel = model;

% Add reaction identifiers
updatedModel.rxns = [model.rxns; newRxnId];
updatedModel.rxnNames = [model.rxnNames; newRxnName];

% Build stoichiometric column
newColumn = zeros(length(updatedModel.mets), 1);

% Add reactants (negative stoichiometry)
for i = 1:length(reactants)
    metIndex = find(ismember(model.metNames, reactants{i}) & ...
        secondToLastCharMatch(model.mets, reactantComps{i}));
    if isempty(metIndex)
        error('Reactant %s[%s] not found', reactants{i}, reactantComps{i});
    end
    newColumn(metIndex) = stoichiometry(i);
end

% Add products (positive stoichiometry)
for i = 1:length(products)
    metIndex = find(ismember(model.metNames, products{i}) & ...
        secondToLastCharMatch(model.mets, productComps{i}));
    if isempty(metIndex)
        error('Product %s[%s] not found', products{i}, productComps{i});
    end
    newColumn(metIndex) = stoichiometry(length(reactants) + i);
end

% Update model matrices
updatedModel.S = [model.S, newColumn];
updatedModel.lb = [model.lb; lb];
updatedModel.ub = [model.ub; ub];
updatedModel.c = [model.c; c];
updatedModel.subSystems = [model.subSystems; {subsystem}];
updatedModel.rxnGeneMat = [model.rxnGeneMat; zeros(1, size(model.rxnGeneMat, 2))];

% Initialize other fields
updatedModel.rules = [model.rules; {''}];
updatedModel.grRules = [model.grRules; {''}];
updatedModel.rxnConfidenceScores = [model.rxnConfidenceScores; NaN];
updatedModel.rxnECNumbers = [model.rxnECNumbers; {''}];
updatedModel.rxnNotes = [model.rxnNotes; {''}];
updatedModel.rxnReferences = [model.rxnReferences; {''}];
end

%% secondToLastCharMatch - Check metabolite compartment
function match = secondToLastCharMatch(mets, comp)
match = false(size(mets));
for j = 1:length(mets)
    if length(mets{j}) >= 2 && mets{j}(end-1) == comp
        match(j) = true;
    end
end
end

%% GEMoptyFUNC_noName - Linear programming optimization for flux bounds
function J_out = GEMoptyFUNC_noName(Model, S, RxnIDX, optyDir, showDebug)
lb0 = Model.lb;
ub0 = Model.ub;
n = length(lb0);

% Set up optimization problem
do_not_optimize_aux = false(n, 1);
do_not_optimize_aux(RxnIDX) = true;

% Objective function
f = zeros(2 * n, 1);
f(n + 1:end) = 1;  % Minimize auxiliary variables

if optyDir == "max"
    f(RxnIDX) = -10000;  % Maximize target reaction
elseif optyDir == "min"
    f(RxnIDX) = 10000;   % Minimize target reaction
end

f(n + find(do_not_optimize_aux)) = 0;

% Constraints
Aeq = [S, zeros(size(S))];  % Mass balance
beq = zeros(size(S, 1), 1);

% Bounds
lb_ext = [lb0; zeros(n, 1)];
ub_ext = [ub0; inf(n, 1)];

% Auxiliary variable constraints: t >= |J|
A = [eye(n), -eye(n); -eye(n), -eye(n)];
b = zeros(2 * n, 1);

% Solve
if showDebug
    optionStruct = struct("Display", "final");
else
    optionStruct = struct("Display", "off");
end
[J_t, ~, ~] = linprog(f, A, b, Aeq, beq, lb_ext, ub_ext, optionStruct);

% Extract flux vector
J_out = J_t(1:n);
end

%% allRxn_reIndexer
function [allRxn_new] = allRxn_reIndexer(Model_old)
allRxn_new = Model_old.allRxn;
removeKeyRxn = false(height(allRxn_new),1);
for a = 1:height(allRxn_new)
    if allRxn_new.RATCON(a) ~= "Blank"
        try
            allRxn_new.IDX(a) = find(Model_old.rxnNames == allRxn_new.RATCON(a));
        catch
            removeKeyRxn(a) = true;
        end
    end
end
allRxn_new(removeKeyRxn,:) = [];
end