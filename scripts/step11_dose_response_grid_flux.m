%%  step11_dose_response_grid_flux

%  Calculates delta-flux for a dose-response grid of drug combinations
%  (anthracycline x OLM concentrations from 0 to 2 in 0.1 steps =
%  21x21 = 441 grid points per combo). Uses convolution-based expression
%  estimates.
%
%  Requirements:
%    - MATLAB
%    - Optimization Toolbox         (linprog, optimvar, optimproblem)
%    - Bioinformatics Toolbox       (mattest, mafdr)
%    - Parallel Computing Toolbox   (parpool, parfor)
%    - HPC cluster strongly recommended
%                                   (441 pts x 4 combos x 6 cells
%                                    = ~10,584 flux calculations)
%
%  Inputs:
%    - out/iCardio_optimized.mat     — ModelOpt (optimized GEM)
%    - <pathRNA>/*_mapped.csv        — RNA-seq expression CSVs
%                                      (pathRNA set by user; %TODO)
%
%  Outputs:
%    - outPf_GL/<Drug1>_<Drug2>_combo/<combo>_d<x>_dd<y>.mat
%                                    — One .mat per grid point with
%                                      Drug_MedianDiff and Conv_MedianDiff

clear; clc; close all;

%% Configuration and Settings

extension = '_mapped.csv'; % File extension for RNA-seq data

pathRNA = ;%TODO

addpath(pathRNA)
if ~exist('outPf_GL','dir'), mkdir('outPf_GL'); end

numExpectedFiles = 266; % Expected number of files

numExpectedDrugs = 54; % Expected number of unique drugs

numExpectedCells = 6; % Expected number of unique cell lines

% User-defined drug combinations (only these pairs will be tested)
DrugCombList = {{"DOX", "OLM"}, {"IDA", "OLM"}, {"DAU", "OLM"},{"EPI", "OLM"}};
DrugCombList_unique = unique(["DOX", "OLM","IDA", "OLM", "DAU", "OLM","EPI", "OLM"], "stable");

drug1_concs = 0:0.1:2;
drug2_concs = 0:0.1:2;

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

%% Build Drug-to-File Map

% Create a mapping of which drugs exist in which cell lines
% drugFileMap{cell_idx, drug_idx} = full file path (or empty if not found)
uniqueCells = unique(cellList);
drugFileMap = cell(length(uniqueCells), length(DrugCombList_unique));

for f = 1:length(fileList)
    [~, cellIdx] = ismember(cellList(f), uniqueCells);
    [~, drugIdx] = ismember(drugList(f), DrugCombList_unique);
    if cellIdx ~= 0 && drugIdx ~= 0
        drugFileMap{cellIdx, drugIdx} = fileList(f);
    end
end

%% Pre-load all Drug RNA Data

% Load and cache RNA data for all drugs in all cell lines to avoid repeated I/O
disp('Pre-loading all drug RNA data...');
drugRNAcache = cell(length(uniqueCells), length(DrugCombList_unique));
ctrlRNAcache = cell(length(uniqueCells), length(DrugCombList_unique));

for cl = 1:length(uniqueCells)
    for d = 1:length(DrugCombList_unique)
        if ~isempty(drugFileMap{cl, d})
            fileName = drugFileMap{cl, d};
            [drugRNA, ctrlRNA, ~, ~] = loadDEGfile(fileName, (cl-1)*length(DrugCombList_unique)+d, numExpectedFiles);
            drugRNAcache{cl, d} = drugRNA;
            ctrlRNAcache{cl, d} = ctrlRNA;
        end
    end
end
disp('RNA data cached.');

%% Main Processing Loop - Nested: Combos and Cell Lines

% Create one row per combination pair
numCombos = length(DrugCombList);

comboPairStrs = cell(numCombos, 1);
for c = 1:numCombos
    comboPairStrs{c} = [char(DrugCombList{c}{1}) '_vs_' char(DrugCombList{c}{1}) char(DrugCombList{c}{2})];
end

numDrug1_concs = length(drug1_concs);
numDrug2_concs = length(drug2_concs);

disp('Starting combination flux calculations...');
for a = 1:numCombos
    % Get the two drugs in this combination
    drug1 = DrugCombList{a}{1};
    drug2 = DrugCombList{a}{2};

    comboName = drug1 + "_" + drug2 + "_combo";

    if ~exist(fullfile('outPf_GL', comboName),'dir'), mkdir(fullfile('outPf_GL', comboName)); end

    % Find indices in the drug list
    [~, drug1_idx] = ismember(drug1, DrugCombList_unique);
    [~, drug2_idx] = ismember(drug2, DrugCombList_unique);

    disp(['Processing combo ' num2str(a) '/' num2str(numCombos) ': ' char(drug1) ' vs ' char(drug2)]);

    tempTable = zeros(numExpectedCells, 1);
    varNames = [string(drug1 + drug2)];
    rowNames = unique(cellList);
    cellName = "cell";
    tableType = repelem(cellName, length(varNames));
    tempTable = array2table(tempTable);
    tempTable.Properties.VariableNames = varNames;
    tempTable.Properties.RowNames = rowNames;
    tempTable.Properties.VariableTypes = tableType;

    poolobj = gcp('nocreate');
    delete(poolobj);
    parpool('local');
    parfor pair_idx = 1:(numDrug1_concs * numDrug2_concs)
        [conc2_idx, conc1_idx] = ind2sub([numDrug2_concs, numDrug1_concs], pair_idx);
        conc1 = drug1_concs(conc1_idx);
        conc2 = drug2_concs(conc2_idx);

        cellLine_results = struct();
        s = [];
        cellLine_results.ctrl1_RxnFlux = tempTable;
        cellLine_results.drug1_RxnFlux= tempTable;
        cellLine_results.ctrl2_RxnFlux = tempTable;
        cellLine_results.conv2_RxnFlux = tempTable;
        
        for cl = 1:length(uniqueCells)
            cellLine = uniqueCells(cl);

            % Check if both drugs exist for this cell line
            if isempty(drugFileMap{cl, drug1_idx}) || isempty(drugFileMap{cl, drug2_idx})
                disp(['  Skipping ' char(cellLine) ' - drug pair not available']);
                continue;
            end

            % Load Drug1 data
            toxDrug_RNA_table = drugRNAcache{cl, drug1_idx};
            toxDrug_CTRL_table = ctrlRNAcache{cl, drug1_idx};

            % Load Drug2 data
            helpDrug_RNA_table = drugRNAcache{cl, drug2_idx};
            helpDrug_CTRL_table = ctrlRNAcache{cl, drug2_idx};

            % Correct Drug Missing Genes
            tempKeyJoinTable_RNA = outerjoin(toxDrug_RNA_table, helpDrug_RNA_table, 'Keys', 'Row');
            tempKeyJoinTable_RNA = fillmissing(tempKeyJoinTable_RNA, 'constant', 1);

            numToxDrug = width(toxDrug_RNA_table);
            numHelpDrug = width(helpDrug_RNA_table);

            toxDrug_RNA_table = tempKeyJoinTable_RNA(:,1:numToxDrug);
            helpDrug_RNA_table = tempKeyJoinTable_RNA(:,numToxDrug+1:end);
            assert(width(helpDrug_RNA_table) == numHelpDrug, "Bad correct missing genes");

            % Correct CTRL Missing Genes
            tempKeyJoinTable_CTRL = outerjoin(toxDrug_CTRL_table, helpDrug_CTRL_table, 'Keys', 'Row');
            tempKeyJoinTable_CTRL = fillmissing(tempKeyJoinTable_CTRL, 'constant', 1);

            numToxCTRL = width(toxDrug_CTRL_table);
            numHelpCTRL = width(helpDrug_CTRL_table);

            toxDrug_CTRL_table = tempKeyJoinTable_CTRL(:,1:numToxCTRL);
            helpDrug_CTRL_table = tempKeyJoinTable_CTRL(:,numToxCTRL+1:end);
            assert(width(helpDrug_CTRL_table) == numHelpCTRL, "Bad correct missing genes");

            % Use toxDrug as control for flux prediction
            medToxDrug_RNA = median(toxDrug_RNA_table,2,"omitmissing");

            % Compute convolution estimate: FC_combo = FC_drug1 * FC_drug2
            conv_FC = helpDrug_RNA_table.^conc2 .* medToxDrug_RNA{:,:}.^conc1;

            % Calculate flux for Drug1 (single treatment)
            [ctrl1_PathFlux, ctrl1_RxnFlux, drug1_PathFlux, drug1_RxnFlux] = runSJ0_DEG_local(...
                Model, toxDrug_RNA_table.^conc1, toxDrug_CTRL_table, cellLine, sprintf('%s_%.1f', char(drug1), conc1));

            % Calculate flux for Convolution (combo estimate)
            [ctrl2_PathFlux, ctrl2_RxnFlux, conv2_PathFlux, conv2_RxnFlux] = runSJ0_DEG_local(...
                Model, conv_FC, toxDrug_CTRL_table, cellLine, sprintf('%s_%.1f_%s_%.1f', char(drug1), conc1, char(drug2), conc2));

            cellLine_results.ctrl1_RxnFlux{cellLine, string(drug1 + drug2)} = {ctrl1_RxnFlux};
            cellLine_results.drug1_RxnFlux{cellLine, string(drug1 + drug2)} = {drug1_RxnFlux};
            cellLine_results.ctrl2_RxnFlux{cellLine, string(drug1 + drug2)} = {ctrl2_RxnFlux};
            cellLine_results.conv2_RxnFlux{cellLine, string(drug1 + drug2)} = {conv2_RxnFlux};
        end

        % MedMed Statistics (Drug vs Mean DMSO across all conditions)
        RxnMedMed_Drug = DToXs_statizizer_MedMed(cellLine_results.drug1_RxnFlux, cellLine_results.ctrl1_RxnFlux);
        Drug_MedianDiff = RxnMedMed_Drug.MedianDiff{:,:};
        % MedMed Statistics (Drug vs Mean DMSO across all conditions)
        RxnMedMed_Conv = DToXs_statizizer_MedMed(cellLine_results.conv2_RxnFlux, cellLine_results.ctrl2_RxnFlux);
        Conv_MedianDiff = RxnMedMed_Conv.MedianDiff{:,:};
        s = struct("Drug_MedianDiff", Drug_MedianDiff,"Conv_MedianDiff", Conv_MedianDiff);
        save(fullfile('outPf_GL', comboName, sprintf('%s_d%.1f_dd%.1f.mat', comboName, conc1, conc2)),"-fromstruct" ,s);
    end
end



%% Done

disp('Driver 11 completed successfully!');

%% ============= Local Helper Functions ========
%% findBaseNames
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

%%  loadDEGfile
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

%% runSJ0_DEG_local
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

%% SJ0_Flux_DEG_linprog_local
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

%% solve_V_d_local
function [V_d_opty] = solve_V_d_local(S, V_g, Jo)
% Solve optimization problem using linprog approach
xOpty = optimvar('x', S.rxnNames);
quadprobN = optimproblem();
quadprobN.Objective = sum(xOpty .* xOpty);
quadprobN.Constraints.sj0 = full(S.S) * (Jo .* (xOpty + V_g)) == 0;
linsol = solve(quadprobN, 'Options', struct('Display', 'none'));
V_d_opty = linsol.x;
end

%% ParseMyRule_max_local
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

%% detect_brackets_local
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
%% DToXs_statizizer_MedMed
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