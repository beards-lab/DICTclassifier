%%  step03_calculate_delta_flux

%  Calculates delta-flux for every drug vs. DMSO comparison across all 54
%  drugs and 6 cell lines. Uses the Noah linprog method with null-space
%  projection for mass-balanced flux prediction.
%
%  Requirements:
%    - MATLAB
%    - Optimization Toolbox         (linprog, optimvar, optimproblem)
%    - Bioinformatics Toolbox       (mattest, mafdr)
%    - Parallel Computing Toolbox   (parpool, parfor)
%    - HPC cluster recommended      (266 RNA-seq files x flux calculations)
%
%  Inputs:
%    - out/iCardio_optimized.mat     — Optimized GEM (ModelOpt)
%    - <pathRNA>/*_mapped.csv        — RNA-seq expression CSVs
%                                      (pathRNA set by user; %TODO)
%
%  Outputs:
%    - out/deltaFlux_stats.mat       — PathStatsMedMed, RxnStatsMedMed,
%                                      RnaStatsMedMed (+ ByCell, ByIndv)
%                                      with median diffs, p-values, and FDR

clear; clc; close all;

%% Configuration and Settings
doAllinParTF = true;            % Enable parallel processing
M = 34;                         % Parallel pool size
extension = '_mapped.csv';       % File extension for RNA-seq data
pathRNA = ; %TODO
addpath(pathRNA)
Method = "Noah";                % Use Noah linprog method only
numExpectedFiles = 266;         % Expected number of files
numExpectedDrugs = 54;          % Expected number of unique drugs
numExpectedCells = 6;           % Expected number of unique cell lines

%% Load Optimized GEM
disp('Loading optimized iCardio GEM...');
load(fullfile('out','iCardio_optimized.mat'), 'ModelOpt');
Model = ModelOpt;
clear ModelOpt;

%% Prepare Model Matrices
disp('Preparing model matrices...');
Model.N = null(full(Model.S));
Model.NNp = Model.N * Model.N';
disp('Model matrices prepared.');

%% Initialize Parallel Pool
if doAllinParTF
    delete(gcp('nocreate'));
    parpool(M);
    disp('Parallel pool initiated.');
else
    disp('No parallel pool - running sequentially.');
end

%% Discover and Validate RNA-seq Files
disp('Discovering RNA-seq files...');
fileList = findBaseNames(pathRNA, extension);
assert(length(fileList) == numExpectedFiles, 'Missing RNA-seq files');

% Extract cell lines and drugs from filenames
cellList = string(cellfun(@(x) x(7:11), cellstr(fileList), 'UniformOutput', false));
assert(length(unique(cellList)) == numExpectedCells, 'Missing cell lines');

drugList = string(cellfun(@(x) x(end-2:end), cellstr(fileList), 'UniformOutput', false));
assert(length(unique(drugList)) == numExpectedDrugs, 'Missing drugs');

disp(['Found ' num2str(length(fileList)) ' RNA-seq files']);
disp(['Cell lines: ' num2str(length(unique(cellList)))]);
disp(['Drugs: ' num2str(length(unique(drugList)))]);

%% Initialize Output Tables
% Create template tables for results
tempTable = zeros(numExpectedCells, numExpectedDrugs);
varNames = unique(drugList);
rowNames = unique(cellList);
cellName = "cell";
tableType = repelem(cellName, length(varNames));
tempTable = array2table(tempTable);
tempTable.Properties.VariableNames = varNames;
tempTable.Properties.RowNames = rowNames;
tempTable.Properties.VariableTypes = tableType;

% Initialize result tables
drugPathFluxSJ0 = tempTable;
drugRxnFluxSJ0 = tempTable;
CTRLPathFluxSJ0 = tempTable;
CTRLRxnFluxSJ0 = tempTable;
drugRNASJ0 = tempTable;
CTRLRNASJ0 = tempTable;

% Initialize parallel result arrays
drugPathFluxSJ0_par = cell(length(fileList), 1);
drugRxnFluxSJ0_par = cell(length(fileList), 1);
CTRLPathFluxSJ0_par = cell(length(fileList), 1);
CTRLRxnFluxSJ0_par = cell(length(fileList), 1);
drugRNASJ0_par = cell(length(fileList), 1);
CTRLRNASJ0_par = cell(length(fileList), 1);

%% Main Processing Loop
disp('Starting flux calculations...');

if doAllinParTF
    % Parallel processing
    parfor a = 1:numExpectedFiles
        [drugRNAtable, CTRLRNAtable, CellLine, DrugExp] = ...
            loadDEGfile(fileList(a), a, length(fileList));

        % Store RNA data
        rnaTable = table;
        rnaTable.GeneName = string(CTRLRNAtable.Properties.RowNames);
        rnaTable.Properties.RowNames = CTRLRNAtable.Properties.RowNames;
        CTRLRNASJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', [rnaTable, CTRLRNAtable]);

        rnaTable = table;
        rnaTable.GeneName = string(drugRNAtable.Properties.RowNames);
        rnaTable.Properties.RowNames = drugRNAtable.Properties.RowNames;
        drugRNASJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', [rnaTable, drugRNAtable]);

        % Calculate flux using Noah method
        [tempRes1, tempRes2, tempRes3, tempRes4] = runSJ0_DEG_local(...
            Model, drugRNAtable, CTRLRNAtable, CellLine, DrugExp);

        % Store flux results
        CTRLPathFluxSJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', tempRes1);
        CTRLRxnFluxSJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', tempRes2);
        drugPathFluxSJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', tempRes3);
        drugRxnFluxSJ0_par{a,1} = struct('CellLine', CellLine, 'DrugExp', DrugExp, ...
            'Data', tempRes4);
    end
else
    % Sequential processing
    for a = 1:length(fileList)
        [drugRNAtable, CTRLRNAtable, CellLine, DrugExp] = ...
            loadDEGfile(fileList(a), a, length(fileList));

        % Store RNA data directly in tables
        rnaTable = table;
        rnaTable.GeneName = string(CTRLRNAtable.Properties.RowNames);
        rnaTable.Properties.RowNames = CTRLRNAtable.Properties.RowNames;
        CTRLRNASJ0{CellLine, DrugExp} = {[rnaTable, CTRLRNAtable]};

        rnaTable = table;
        rnaTable.GeneName = string(drugRNAtable.Properties.RowNames);
        rnaTable.Properties.RowNames = drugRNAtable.Properties.RowNames;
        drugRNASJ0{CellLine, DrugExp} = {[rnaTable, drugRNAtable]};

        % Calculate flux using Noah method
        [tempRes1, tempRes2, tempRes3, tempRes4] = runSJ0_DEG_local(...
            Model, drugRNAtable, CTRLRNAtable, CellLine, DrugExp);

        % Store flux results directly in tables
        CTRLPathFluxSJ0{CellLine, DrugExp} = {tempRes1};
        CTRLRxnFluxSJ0{CellLine, DrugExp} = {tempRes2};
        drugPathFluxSJ0{CellLine, DrugExp} = {tempRes3};
        drugRxnFluxSJ0{CellLine, DrugExp} = {tempRes4};
    end
end

%% Process Parallel Results into Tables
if doAllinParTF
    delete(gcp('nocreate'));

    % Convert parallel results to tables
    [drugPathFluxSJ0, drugRxnFluxSJ0, CTRLPathFluxSJ0, CTRLRxnFluxSJ0] = ...
        fluxStruct2table(fileList, drugPathFluxSJ0_par, drugRxnFluxSJ0_par, ...
        CTRLPathFluxSJ0_par, CTRLRxnFluxSJ0_par, ...
        drugPathFluxSJ0, drugRxnFluxSJ0, CTRLPathFluxSJ0, CTRLRxnFluxSJ0);

    [CTRLRNASJ0, drugRNASJ0] = rnaStruct2table(fileList, CTRLRNASJ0_par, ...
        drugRNASJ0_par, CTRLRNASJ0, drugRNASJ0);
end

disp('Files processed. Starting statistical analysis...');

%% Statistical Analysis

% MedMed Statistics (Drug vs Mean DMSO across all conditions)
PathStatsMedMed = DToXs_statizizer_MedMed(drugPathFluxSJ0, CTRLPathFluxSJ0);
RxnStatsMedMed = DToXs_statizizer_MedMed(drugRxnFluxSJ0, CTRLRxnFluxSJ0);
RnaStatsMedMed = DToXs_statizizer_MedMed(drugRNASJ0, CTRLRNASJ0);
disp('MedMed Stats Done');

% ByCell Statistics (Drug vs DMSO within same cell line)
PathStatsByCell = DToXs_statizizer_ByCell(drugPathFluxSJ0, CTRLPathFluxSJ0);
RxnStatsByCell = DToXs_statizizer_ByCell(drugRxnFluxSJ0, CTRLRxnFluxSJ0);
RnaStatsByCell = DToXs_statizizer_ByCell(drugRNASJ0, CTRLRNASJ0);
disp('ByCell Stats Done');

% ByIndv Statistics (Individual sample comparisons)
PathStatsByIndv = DToXs_statizizer_ByIndv(drugPathFluxSJ0, CTRLPathFluxSJ0);
RxnStatsByIndv = DToXs_statizizer_ByIndv(drugRxnFluxSJ0, CTRLRxnFluxSJ0);
RnaStatsByIndv = DToXs_statizizer_ByIndv(drugRNASJ0, CTRLRNASJ0);
disp('ByIndv Stats Done');

%% Save Results
disp('Saving results...');
if ~exist('out','dir'), mkdir('out'); end

save(fullfile('out','deltaFlux_stats.mat'), ...
    'PathStatsMedMed','PathStatsByCell','PathStatsByIndv', ...
    'RxnStatsMedMed','RxnStatsByCell','RxnStatsByIndv', ...
    'RnaStatsMedMed','RnaStatsByCell','RnaStatsByIndv', '-v7.3');

disp('Driver 3 completed successfully!');

%% Local Helper Functions

function nameList = findBaseNames(rootPath, extension)
% Get all files with specified extension recursively
currentFiles = dir(fullfile(rootPath, ['*' extension]));
isFile = ~[currentFiles.isdir];
fileNames = {currentFiles(isFile).name};

% Extract base names without extension
baseNames = cellfun(@(f) erase(f, extension), fileNames, 'UniformOutput', false);

% Recursively process subdirectories
subDirs = dir(rootPath);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name}, {'.', '..'}));

for i = 1:numel(subDirs)
    subPath = fullfile(rootPath, subDirs(i).name);
    baseNames = [baseNames, findBaseNames(subPath, extension)];
end

nameList = string(baseNames);
end

function [drugRNAtable, CTRLRNAtable, CellLine, DrugExp] = loadDEGfile(file, a, allNum)
% Extract cell line and drug from filename
CellLine = string(unique(cellfun(@(x) x(7:11), cellstr(file), 'UniformOutput', false)));
DrugExp = string(unique(cellfun(@(x) x(end-2:end), cellstr(file), 'UniformOutput', false)));

disp([num2str(a) '/' num2str(allNum) ': ' char(CellLine) ' - ' char(DrugExp)]);

% Load RNA-seq data
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
fullRNAtable = readtable([char(file) '_mapped.csv'], 'VariableNamingRule', 'modify');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

% Remove genes with missing EntrezID
fullRNAtable(ismissing(fullRNAtable.ncbiID), :) = [];

% Handle duplicated gene IDs by keeping highest expression
[C, ~, ic] = unique(fullRNAtable.ncbiID, 'stable');
EntrezFreqs = [C, accumarray(ic, 1)];
EntrezFreqs = EntrezFreqs(EntrezFreqs(:,2) > 1, :);

if ~isempty(EntrezFreqs)
    multTable = fullRNAtable(ismember(fullRNAtable.ncbiID, EntrezFreqs(:,1)), :);
    multTable = sortrows(multTable, 'logCPM', 'descend');
    fullRNAtable(ismember(fullRNAtable.ncbiID, EntrezFreqs(:,1)), :) = [];

    % Keep only the highest expressing duplicate for each gene
    dupdEntrez = unique(multTable.ncbiID);
    for c = 1:length(dupdEntrez)
        dupdTemp = multTable(multTable.ncbiID == dupdEntrez(c), :);
        newRow = dupdTemp(1, :);
        newRow{1, 2:end-7} = sum(dupdTemp{:, 2:end-7}, 1);
        fullRNAtable = [fullRNAtable; newRow];
    end
end

fullRNAtable = sortrows(fullRNAtable, 'ncbiID');
assert(length(unique(fullRNAtable.ncbiID)) == height(fullRNAtable), 'Duplicate genes remain');

fullRNAtable.Properties.RowNames = string(fullRNAtable.ncbiID);

% Split by experimental vs control conditions
drugRNAtable = fullRNAtable(:, contains(string(fullRNAtable.Properties.VariableNames), CellLine) & ...
    contains(string(fullRNAtable.Properties.VariableNames), DrugExp) & ...
    contains(string(fullRNAtable.Properties.VariableNames), "_Norm"));

CTRLRNAtable = fullRNAtable(:, contains(string(fullRNAtable.Properties.VariableNames), CellLine) & ...
    contains(string(fullRNAtable.Properties.VariableNames), "CTRL") & ...
    contains(string(fullRNAtable.Properties.VariableNames), "_Norm"));

% Calculate average control expression for normalization
avgCTRL = mean(CTRLRNAtable, 2);
%assert(~any(ismissing(avgCTRL)), 'Missing control expression');

MINctrlExpres = min(avgCTRL{avgCTRL{:,:} ~= 0,1});
avgCTRL{avgCTRL{:,:} == 0,1} = MINctrlExpres;
%assert(all(avgCTRL > 0), 'Zero control expression found');

% Calculate fold changes

drugRNAtable{:,:} = drugRNAtable{:,:} ./ avgCTRL{:,:};
CTRLRNAtable{:,:} = CTRLRNAtable{:,:} ./ avgCTRL{:,:};

end

function [CTRLPathFlux_SJ0, CTRLRxnFlux_SJ0, drugPathFlux_SJ0, drugRxnFlux_SJ0] = ...
    runSJ0_DEG_local(Model, drugRNA_table, CTRLRNA_table, CellLine, DrugExp)

% Initialize output tables
drugPathFlux_SJ0 = table();
drugPathFlux_SJ0.Pathways = unique(string(Model.subSystems));
drugPathFlux_SJ0.Pathways(1) = "blank";
drugPathFlux_SJ0.Properties.RowNames = drugPathFlux_SJ0.Pathways;

drugRxnFlux_SJ0 = table();
drugRxnFlux_SJ0.RATCON = Model.rxnNames;
drugRxnFlux_SJ0.Properties.RowNames = drugRxnFlux_SJ0.RATCON;

CTRLPathFlux_SJ0 = drugPathFlux_SJ0;
CTRLRxnFlux_SJ0 = drugRxnFlux_SJ0;

% Process gene expression data
GeneNames_exp = double(string(drugRNA_table.Properties.RowNames));
GeneNames_ctrl = double(string(CTRLRNA_table.Properties.RowNames));
assert(all(GeneNames_exp == GeneNames_ctrl), 'Gene lists do not match');

GeneNames = intersect(GeneNames_exp, GeneNames_ctrl);
numDrug = width(drugRNA_table);
numCTRL = width(CTRLRNA_table);

disp("Starting N*N'*l_f loop");

% Map genes to model
EntrezKeeps = intersect(GeneNames, str2double(Model.genes));
C_drugRNA_table = drugRNA_table(string(EntrezKeeps), :);
C_CTRLRNA_table = CTRLRNA_table(string(EntrezKeeps), :);

% Create gene expression tables for model
drugModelGenes = table();
drugModelGenes.genes = Model.genes;
drugModelGenes.Properties.RowNames = drugModelGenes.genes;
CTRLModelGenes = drugModelGenes;

% Map gene expression to model genes
for a = 1:height(drugModelGenes)
    tempEntrez = string(str2double(drugModelGenes.genes(a)));
    if any(string(C_drugRNA_table.Properties.RowNames) == tempEntrez)
        tempDrugFC = C_drugRNA_table(tempEntrez, :);
        tempCTRLFC = C_CTRLRNA_table(tempEntrez, :);

        if width(drugModelGenes) == 1
            for b = 1:width(tempDrugFC)
                drugModelGenes.(string(tempDrugFC.Properties.VariableNames{b})) = ...
                    ones(height(drugModelGenes), 1);

                drugModelGenes(a, end) = tempDrugFC(:, b);
            end
            for b = 1:width(tempCTRLFC)
                CTRLModelGenes.(string(tempCTRLFC.Properties.VariableNames{b})) = ...
                    ones(height(CTRLModelGenes), 1);

                CTRLModelGenes(a, end) = tempCTRLFC(:, b);
            end
        else
            drugModelGenes(a, end-width(tempDrugFC)+1:end) = tempDrugFC;
            CTRLModelGenes(a, end-width(tempCTRLFC)+1:end) = tempCTRLFC;
        end
    end
end

% Process drug samples
for a = 1:numDrug
    disp(newline() + "Exp: Starting " + CellLine + " for " + DrugExp + ". " + a + "/" + numDrug);
    [drugPathFlux_SJ0, drugRxnFlux_SJ0] = SJ0_Flux_DEG_linprog_local(...
        Model, string(drugRNA_table.Properties.VariableNames{a}), ...
        drugPathFlux_SJ0, drugRxnFlux_SJ0, drugModelGenes(:, [1, 1+a]));
end

% Process control samples
for a = 1:numCTRL
    disp(newline() + "CTRL: Starting " + CellLine + " for " + DrugExp + ". " + a + "/" + numCTRL);
    [CTRLPathFlux_SJ0, CTRLRxnFlux_SJ0] = SJ0_Flux_DEG_linprog_local(...
        Model, string(CTRLRNA_table.Properties.VariableNames{a}), ...
        CTRLPathFlux_SJ0, CTRLRxnFlux_SJ0, CTRLModelGenes(:, [1, 1+a]));
end
end

function [PathFlux_SJ0, RxnFlux_SJ0] = SJ0_Flux_DEG_linprog_local(...
    Model, sampleID, PathFlux_SJ0, RxnFlux_SJ0, ModelGenes)

% Set up reaction rules and gene expression matrix
ModelRxns = table();
ModelRxns.rxns = string(Model.rxns);
ModelRxns.rules = string(Model.rules);
ModelRxns.rxnFCmat = Model.rxnGeneMat .* ModelGenes{:,end}';

% Evaluate gene expression rules for each reaction
for a = 1:height(ModelRxns)
    if ModelRxns.rules(a) ~= ""
        tempRule = ModelRxns.rules(a);
        x = ModelRxns.rxnFCmat(a, :);

        if contains(tempRule, '&') || contains(tempRule, '|')
            finalstring = ParseMyRule_max_local(tempRule, x);
            ModelRxns.Rules_simple(a) = strcat('x', finalstring);
            ModelRxns.rxn_Fold(a) = eval(strcat('x', finalstring));
        else
            ModelRxns.Rules_simple(a) = tempRule;
            ModelRxns.rxn_Fold(a) = eval(tempRule);
        end
    end
end

% Calculate null space and delta flux
l_f = full(ModelRxns.rxn_Fold);
l_f(Model.ub > 0 & Model.lb < 0) = 1;
l_fones = ones(height(l_f), 1);

J0 = Model.NNp * l_fones;
Vd = solve_V_d_local(Model, l_f, J0);
newMult = l_f + Vd;
Jnew = J0 .* newMult;
perChange = (Jnew - J0) ./ J0;
deltaJ = perChange;
deltaJ(logical(Model.isReversed)) = deltaJ(logical(Model.isReversed)) .* -1;

disp("Done calculating N*N'*l_f of " + sampleID);

% Group by subsystems
[G, ID] = findgroups(string(Model.subSystems));
Group_Mean = zeros(length(ID), 1);
for a = 1:length(ID)
    Group_Mean(a) = mean(deltaJ(G == a), 'omitmissing');
end

ID(1) = "blank";
Subsystem = table();
Subsystem.(sampleID) = Group_Mean;
Subsystem.Properties.RowNames = ID;
PathFlux_SJ0 = [PathFlux_SJ0, Subsystem];

% Store reaction-level results
RxnChange = table();
RxnChange.(sampleID) = deltaJ;
RxnChange.Properties.RowNames = string(Model.rxnNames);
RxnFlux_SJ0 = [RxnFlux_SJ0, RxnChange];
end

function [V_d_opty] = solve_V_d_local(S, V_g, Jo)
% Solve optimization problem using linprog approach
xOpty = optimvar('x', S.rxnNames);
quadprobN = optimproblem();
quadprobN.Objective = sum(xOpty .* xOpty);
quadprobN.Constraints.sj0 = full(S.S) * (Jo .* (xOpty + V_g)) == 0;
linsol = solve(quadprobN, 'Options', struct('Display', 'none'));
V_d_opty = linsol.x;
end

function finalstring = ParseMyRule_max_local(S, x)
%detect where brackets occur
S=char(S);
[oc]=detect_brackets_local(S);

%Determine the largest brackets that do not contain subbrackets
hold_strings=strings(1); %will contain substrings
leftovers=strings(1); % will contain characters not in substrings (should be opperators)
oc2=oc;
if oc2(end,1)==1
    oc2(end,:)=[];
end
j=1;
i=1;
while i<max(oc2(:,2))
    [~,index]=min(oc2(:,1));
    hold_strings(j)=string(S(oc2(index,1):oc2(index,2)));
    if i>1
        leftovers=strcat(leftovers,string(S(i:oc2(index,1))));
    end
    i=oc2(index,2);
    oc2(oc2(:,1)<oc2(index,2),:)=[];
    j=j+1;
end

for i=1:length(hold_strings)
    if  any(regexp(hold_strings(i) ,'[0-9]')) == 1 && contains(hold_strings(i),'x')==1
        hold_strings(i)=ParseMyRule_max_local(hold_strings(i),x);
    end
end
opp=0;
if contains(leftovers, "|")==1 && contains(leftovers, "&")==1
    disp('ERROR AT BASE LEVEL')
elseif contains(leftovers, "|")==1 && contains(leftovers, "&")==0
    opp=1; %USE OR statement to take MAXIMUM
elseif contains(leftovers, "|")==0 && contains(leftovers, "&")==1
    opp=2; %USE AND statement to take MINIMUM
end
hold_values=zeros(length(hold_strings),1);
for i=1:length(hold_strings)
    hold_values(i)=eval(strcat('x',hold_strings(i)));
end
if opp==1 % Maximum
    [~,index]=max(hold_values);
    finalstring=hold_strings(index);
elseif opp==2  % Minimum
    [~,index]=min(hold_values);
    finalstring=hold_strings(index);
elseif opp==0 & isscalar(hold_strings)
    finalstring=hold_strings;
else
    fprintf("more than two substrings with no opperator \n")
end
end

function [oc] = detect_brackets_local(str)
oc = [];
% find all opening and closing brackets in the string
op=strfind(str,'(');
cl=strfind(str,')');
% search for pairs until all are identified
while ~~any(op | cl, "all") %was ~isempty
    % find opening bracket for first closing bracket
    idx = find(op < cl(1),1,'last');
    % append this pair to function output
    oc = [oc;op(idx) cl(1)];
    % remove found opening bracket from vector
    op(idx) = [];
    % remove found closing bracket from vector
    cl(1) = [];
end
end

function [drugPathFluxSJ0, drugRxnFluxSJ0, CTRLPathFluxSJ0, CTRLRxnFluxSJ0] = ...
    fluxStruct2table(filename, drugPathFluxSJ0_par, drugRxnFluxSJ0_par, ...
    CTRLPathFluxSJ0_par, CTRLRxnFluxSJ0_par, ...
    drugPathFluxSJ0, drugRxnFluxSJ0, CTRLPathFluxSJ0, CTRLRxnFluxSJ0)
% Convert parallel cell arrays to tables
for i = 1:length(filename)
    if ~isempty(drugPathFluxSJ0_par{i})
        CellLine = drugPathFluxSJ0_par{i}.CellLine;
        DrugExp = drugPathFluxSJ0_par{i}.DrugExp;
        drugPathFluxSJ0(CellLine, DrugExp) = {drugPathFluxSJ0_par{i}.Data};
        drugRxnFluxSJ0(CellLine, DrugExp) = {drugRxnFluxSJ0_par{i}.Data};
        CTRLPathFluxSJ0(CellLine, DrugExp) = {CTRLPathFluxSJ0_par{i}.Data};
        CTRLRxnFluxSJ0(CellLine, DrugExp) = {CTRLRxnFluxSJ0_par{i}.Data};
    end
end
end

function [CTRLRNASJ0, drugRNASJ0] = rnaStruct2table(...
    filename, CTRLRNASJ0_par, drugRNASJ0_par, CTRLRNASJ0, drugRNASJ0)
% Convert parallel RNA cell arrays to tables
for i = 1:length(filename)
    if ~isempty(CTRLRNASJ0_par{i})
        CellLine = CTRLRNASJ0_par{i}.CellLine;
        DrugExp = CTRLRNASJ0_par{i}.DrugExp;
        CTRLRNASJ0(CellLine, DrugExp) = {CTRLRNASJ0_par{i}.Data};
        drugRNASJ0(CellLine, DrugExp) = {drugRNASJ0_par{i}.Data};
    end
end
end

function StatStruct = DToXs_statizizer_MedMed(DRUGtableIN, CTRLtableIN)
StatStruct = struct();
StatTable = table();
DrugNames = string(DRUGtableIN.Properties.VariableNames);
CellLines = string(DRUGtableIN.Properties.RowNames);
for a = 1:width(DRUGtableIN)
    for b = 1:height(DRUGtableIN)
        try
            COMPexist(b,a) = double(DRUGtableIN{b,a}{1} ~= 0);
        catch
            COMPexist(b,a) = 1;
            if isempty(StatTable)
                StatTable.name = DRUGtableIN{b,a}{1}{:,1};
                StatName = string(DRUGtableIN{b,a}{1}.Properties.VariableNames{1});
            end
        end
    end
end


StatTable.Properties.RowNames = StatTable.name;
StatTable(:,1) = [];
FDRtable = StatTable;
Meantable_Drug = StatTable;
SDtable_Drug = StatTable;
Meantable_CTRL = StatTable;
SDtable_CTRL = StatTable;
pvaltable = StatTable;
Mediantable_Drug = StatTable;
Mediantable_CTRL = StatTable;
tempTableStructure = StatTable;

for a = 1:width(COMPexist)
    DRUGstatTable_temp = tempTableStructure;
    CTRLstatTable_temp = tempTableStructure;
    drugExistTF = 0;
    for b = 1:height(COMPexist)
        if COMPexist(b,a)
            drugExistTF = 1;
            if StatName == "GeneName"
                DRUGstatTable_temp = outerjoin(DRUGstatTable_temp, DRUGtableIN{b,a}{1}(:,2:end), 'Keys', 'Row', 'MergeKeys', true);
                CTRLstatTable_temp = outerjoin(CTRLstatTable_temp, CTRLtableIN{b,a}{1}(:,2:end), 'Keys', 'Row', 'MergeKeys', true);
            else
                DRUGstatTable_temp = [DRUGstatTable_temp, DRUGtableIN{b,a}{1}(:,2:end)];
                CTRLstatTable_temp = [CTRLstatTable_temp, CTRLtableIN{b,a}{1}(:,2:end)];
            end
        end
    end
    if drugExistTF
        try
            pvals = mattest(DRUGstatTable_temp{:,:},CTRLstatTable_temp{:,:},'VarType','unequal','permute',true);
        catch
            pvals = ones(height(DRUGstatTable_temp),1);
        end
        if isempty(pvals) | ismissing(pvals)
            pvals = ones(height(DRUGstatTable_temp),1);
        end
        try
            fdr = mafdr(pvals,'Method','bootstrap','Showplot',false);
        catch
            fdr = pvals.*height(DRUGstatTable_temp);
        end

        if StatName == "GeneName"
            disp("MAFDR'd " + a + "/" + width(COMPexist))
            StatTable = table();
            StatTable.names = string(DRUGstatTable_temp.Properties.RowNames);
            StatTable.Properties.RowNames = string(DRUGstatTable_temp.Properties.RowNames);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = fdr;
            FDRtable = outerjoin(FDRtable, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = pvals;
            pvaltable = outerjoin(pvaltable, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = mean(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
            Meantable_Drug = outerjoin(Meantable_Drug, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = mean(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
            Meantable_CTRL = outerjoin(Meantable_CTRL, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = median(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
            Mediantable_Drug = outerjoin(Mediantable_Drug, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = median(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
            Mediantable_CTRL = outerjoin(Mediantable_CTRL, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = std(DRUGstatTable_temp{:,:}, [], 2, 'omitmissing');
            SDtable_Drug = outerjoin(SDtable_Drug, StatTable, "Keys","Row", "MergeKeys",true);
            StatTable(:,1) = [];
            StatTable.(DrugNames(a)) = std(CTRLstatTable_temp{:,:}, [], 2, 'omitmissing');
            SDtable_CTRL = outerjoin(SDtable_CTRL, StatTable, "Keys","Row", "MergeKeys",true);
        else
            disp("MAFDR'd " + a + "/" + width(COMPexist))
            FDRtable.(DrugNames(a)) = fdr;
            pvaltable.(DrugNames(a)) = pvals;
            Meantable_Drug.(DrugNames(a)) = mean(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
            Meantable_CTRL.(DrugNames(a)) = mean(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
            Mediantable_Drug.(DrugNames(a)) = median(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
            Mediantable_CTRL.(DrugNames(a)) = median(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
            SDtable_Drug.(DrugNames(a)) = std(DRUGstatTable_temp{:,:}, [], 2, 'omitmissing');
            SDtable_CTRL.(DrugNames(a)) = std(CTRLstatTable_temp{:,:}, [], 2, 'omitmissing');
        end
    end
end
StatStruct.FDR = FDRtable;
StatStruct.pVal = pvaltable;
StatStruct.MeanDrug = Meantable_Drug;
StatStruct.MeanCTRL = Meantable_CTRL;
StatStruct.MeanDiff = Meantable_Drug - Meantable_CTRL;
StatStruct.SDDrug = SDtable_Drug;
StatStruct.SDCTRL = SDtable_CTRL;
StatStruct.MedianDrug = Mediantable_Drug;
StatStruct.MedianCTRL = Mediantable_CTRL;
StatStruct.MedianDiff = Mediantable_Drug - Mediantable_CTRL;
end

function StatStruct = DToXs_statizizer_ByCell(DRUGtableIN, CTRLtableIN)
StatStruct = struct();
StatTable = table();
DrugNames = string(DRUGtableIN.Properties.VariableNames);
CellLines = string(DRUGtableIN.Properties.RowNames);
for a = 1:width(DRUGtableIN)
    for b = 1:height(DRUGtableIN)
        try
            COMPexist(b,a) = double(DRUGtableIN{b,a}{1} ~= 0);
        catch
            COMPexist(b,a) = 1;
            if isempty(StatTable)
                StatTable.name = DRUGtableIN{b,a}{1}{:,1};
                StatName = string(DRUGtableIN{b,a}{1}.Properties.VariableNames{1});
            end
        end
    end
end
StatTable.Properties.RowNames = StatTable.name;
StatTable(:,1) = [];
FDRtable = StatTable;
Meantable_Drug = StatTable;
SDtable_Drug = StatTable;
Meantable_CTRL = StatTable;
SDtable_CTRL = StatTable;
pvaltable = StatTable;
Mediantable_Drug = StatTable;
Mediantable_CTRL = StatTable;
tempTableStructure = StatTable;

for a = 1:width(COMPexist)
    DRUGstatTable_temp = tempTableStructure;
    CTRLstatTable_temp = tempTableStructure;
    drugExistTF = 0;
    for b = 1:height(COMPexist)
        if COMPexist(b,a)
            drugExistTF = 1;
            DRUGstatTable_temp = DRUGtableIN{b,a}{1}(:,2:end);
            CTRLstatTable_temp = CTRLtableIN{b,a}{1}(:,2:end);
            try
                pvals = mattest(DRUGstatTable_temp{:,:},CTRLstatTable_temp{:,:},'VarType','unequal','permute',true);
            catch
                pvals = ones(height(DRUGstatTable_temp),1);
            end
            if isempty(pvals) | ismissing(pvals)
                pvals = ones(height(DRUGstatTable_temp),1);
            end
            try
                fdr = mafdr(pvals,'Method','bootstrap','Showplot',false);
            catch
                fdr = pvals.*height(DRUGstatTable_temp);
            end
            disp("MAFDR'd " + a + "/" + width(COMPexist))
            if StatName == "GeneName"
                StatTable = table();
                StatTable.name = string(DRUGstatTable_temp.Properties.RowNames);
                StatTable.Properties.RowNames = StatTable.name;
                StatTable(:,1) = [];
                %outerjoin(  ,   , 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = fdr;
                FDRtable = outerjoin(FDRtable,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = pvals;
                pvaltable = outerjoin(pvaltable,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = mean(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
                Meantable_Drug = outerjoin(Meantable_Drug,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = mean(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
                Meantable_CTRL = outerjoin(Meantable_CTRL,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = median(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
                Mediantable_Drug = outerjoin(Mediantable_Drug,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = median(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
                Mediantable_CTRL = outerjoin(Mediantable_CTRL,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = std(DRUGstatTable_temp{:,:}, [], 2, 'omitmissing');
                SDtable_Drug = outerjoin(SDtable_Drug,StatTable, 'Keys', 'Row', 'MergeKeys', true);

                StatTable.(string(DrugNames(a)+"_"+CellLines(b))) = std(CTRLstatTable_temp{:,:}, [], 2, 'omitmissing');
                SDtable_CTRL = outerjoin(SDtable_CTRL,StatTable, 'Keys', 'Row', 'MergeKeys', true);

            else
                FDRtable.(string(DrugNames(a)+"_"+CellLines(b))) = fdr;
                pvaltable.(string(DrugNames(a)+"_"+CellLines(b))) = pvals;
                Meantable_Drug.(string(DrugNames(a)+"_"+CellLines(b))) = mean(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
                Meantable_CTRL.(string(DrugNames(a)+"_"+CellLines(b))) = mean(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
                Mediantable_Drug.(string(DrugNames(a)+"_"+CellLines(b))) = median(DRUGstatTable_temp{:,:}, 2, 'omitmissing');
                Mediantable_CTRL.(string(DrugNames(a)+"_"+CellLines(b))) = median(CTRLstatTable_temp{:,:}, 2, 'omitmissing');
                SDtable_Drug.(string(DrugNames(a)+"_"+CellLines(b))) = std(DRUGstatTable_temp{:,:}, [], 2, 'omitmissing');
                SDtable_CTRL.(string(DrugNames(a)+"_"+CellLines(b))) = std(CTRLstatTable_temp{:,:}, [], 2, 'omitmissing');
            end
        end
    end
end
StatStruct.FDR = FDRtable;
StatStruct.pVal = pvaltable;
StatStruct.MeanDrug = Meantable_Drug;
StatStruct.MeanCTRL = Meantable_CTRL;
StatStruct.MeanDiff = Meantable_Drug - Meantable_CTRL;
StatStruct.SDDrug = SDtable_Drug;
StatStruct.SDCTRL = SDtable_CTRL;
StatStruct.MedianDrug = Mediantable_Drug;
StatStruct.MedianCTRL = Mediantable_CTRL;
StatStruct.MedianDiff = Mediantable_Drug - Mediantable_CTRL;
end

function StatStruct = DToXs_statizizer_ByIndv(DRUGtableIN, CTRLtableIN)
StatStruct = struct();
StatTable = table();
DrugNames = string(DRUGtableIN.Properties.VariableNames);
CellLines = string(DRUGtableIN.Properties.RowNames);
for a = 1:width(DRUGtableIN)
    for b = 1:height(DRUGtableIN)
        try
            COMPexist(b,a) = double(DRUGtableIN{b,a}{1} ~= 0);
        catch
            COMPexist(b,a) = 1;
            if isempty(StatTable)
                StatTable.name = DRUGtableIN{b,a}{1}{:,1};
                StatName = string(DRUGtableIN{b,a}{1}.Properties.VariableNames{1});
            end
        end
    end
end
StatTable.Properties.RowNames = StatTable.name;
StatTable(:,1) = [];
FC_Drug = StatTable;
FC_CTRL = StatTable;
FC_Drug_minMed = StatTable;
FC_Drig_minMean = StatTable;
tempTableStructure = StatTable;

for a = 1:width(COMPexist)
    DRUGstatTable_temp = tempTableStructure;
    CTRLstatTable_temp = tempTableStructure;
    for b = 1:height(COMPexist)
        if COMPexist(b,a)
            DRUGstatTable_temp = DRUGtableIN{b,a}{1}(:,2:end);
            CTRLstatTable_temp = CTRLtableIN{b,a}{1}(:,2:end);
            disp("No MAFDR'd for Indv stats")
            if StatName == "GeneName"
                for c = 1:width(DRUGstatTable_temp)
                    StatTable = table();
                    StatTable.name = string(DRUGstatTable_temp.Properties.RowNames);
                    StatTable.Properties.RowNames  = DRUGstatTable_temp.Properties.RowNames;
                    StatTable(:,1) = [];
                    StatTable.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c};
                    FC_Drug = outerjoin(FC_Drug,  StatTable , 'Keys', 'Row', 'MergeKeys', true);
                    StatTable.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c} - median(CTRLstatTable_temp{:,:},2, "omitmissing");
                    FC_Drug_minMed = outerjoin(FC_Drug_minMed,  StatTable , 'Keys', 'Row', 'MergeKeys', true);
                    StatTable.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c} - mean(CTRLstatTable_temp{:,:},2, 'omitmissing');
                    FC_Drig_minMean = outerjoin(FC_Drig_minMean,  StatTable , 'Keys', 'Row', 'MergeKeys', true);
                end
                for c = 1:width(CTRLstatTable_temp)
                    StatTable = table();
                    StatTable.name = string(DRUGstatTable_temp.Properties.RowNames);
                    StatTable.Properties.RowNames  = DRUGstatTable_temp.Properties.RowNames;
                    StatTable(:,1) = [];
                    StatTable.(string(DrugNames(a)+"_"+CellLines(b)+"_CTRL"+c)) = CTRLstatTable_temp{:,c};
                    FC_CTRL = outerjoin(FC_CTRL,  StatTable , 'Keys', 'Row', 'MergeKeys', true);
                end
            else
                for c = 1:width(DRUGstatTable_temp)
                    FC_Drug.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c};
                    FC_Drug_minMed.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c} - median(CTRLstatTable_temp{:,:},2, "omitmissing");
                    FC_Drig_minMean.(string(DrugNames(a)+"_"+CellLines(b)+"_Drug"+c)) = DRUGstatTable_temp{:,c} - mean(CTRLstatTable_temp{:,:},2, 'omitmissing');
                end
                for c = 1:width(CTRLstatTable_temp)
                    FC_CTRL.(string(DrugNames(a)+"_"+CellLines(b)+"_CTRL"+c)) = CTRLstatTable_temp{:,c};
                end
            end
        end
    end
end
StatStruct.FCDrug = FC_Drug;
StatStruct.FCCTRL = FC_CTRL;
StatStruct.MeanDiff = FC_Drig_minMean;
StatStruct.MedianDiff = FC_Drug_minMed;
end