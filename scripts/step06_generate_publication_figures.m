%%  step06_generate_publication_figures

%  Generates all main publication figures (metric plots, ROC/PRC curves,
%  t-SNE, correlograms, irreversibility comparison). Uses self-contained
%  plotting functions.
%
%  Requirements:
%    - MATLAB
%    - Statistics and Machine Learning Toolbox
%                                   (tsne, boxchart, swarmchart)
%
%  Inputs:
%    - out/trainedModels_comprehensive.mat   — Metrics and ROC/PRC objects
%    - out/deltaFlux_stats.mat               — RxnStatsMedMed,
%                                              RnaStatsMedMed,
%                                              RxnStatsByIndv, RnaStatsByIndv
%    - out/ROR_results.mat                   — regresTable (ROR data)
%    - data/drugMap.xlsx                     — Drug mapping table
%    - models/HeartModel.mat                 — Original heart model
%                                              (irreversibility panel)
%    - models/objRxns.mat                    — Objective reaction list
%                                              (irreversibility panel)
%    - out/iCardio_optimized.mat             — ModelOpt
%                                              (irreversibility panel)
%
%  Outputs:
%    - figures/Figure_<N>_metrics.png/.svg/.fig    — Per-metric panels
%    - figures/NEW_ROC_PRC_curves.png/.svg/.fig    — ROC and PRC curves
%    - figures/tSNE_Rxn_analysis.png/.svg/.fig     — t-SNE reaction plot
%    - figures/tSNE_RNA_analysis.png/.svg/.fig     — t-SNE RNA plot
%    - figures/Combined_publication_summary.png/.svg/.fig
%    - figures/Panel4_irreversibility_comparison.png/.fig
%    - figures/Correlogram_Reaction_31drugs.png/.svg/.fig
%    - figures/Correlogram_RNA_31drugs.png/.svg/.fig
%    - figures/Correlogram_Rxn_ClassMarginals.png/.svg
%    - figures/Correlogram_RNA_ClassMarginals.png/.svg
%    - figures/generation_complete.mat             — metrics, iRange, cMax

clear; clc; close all;

%% Load data required for plotting

fprintf('Loading comprehensive model results and statistics...\n');

% Load comprehensive metrics from driver5
load(fullfile('out','trainedModels_comprehensive.mat'));

% Load flux statistics from driver3
load(fullfile('out','deltaFlux_stats.mat'), 'RxnStatsMedMed', 'RnaStatsMedMed', 'RxnStatsByIndv', 'RnaStatsByIndv');

% Load ROR results from driver1
load(fullfile('out','ROR_results.mat'), 'regresTable');

% Load drug mapping
drugMap = readtable(fullfile('data', 'drugMap.xlsx'));

%% Create output folder

if ~exist('figures','dir'), mkdir figures; end

%% Panel 1 — Forest plots & violin plots with optimized Pub styling

fprintf('Generating Panel 1: Forest plots and violin plots...\n');

% Generate optimized Pub plots
generateOptimizedPubPlots(Rxn, Rna, metrics);

%% Panel 2 — NEW IMPROVED ROC/PRC curves

fprintf('Generating Panel 2: NEW ROC and PRC curves...\n');

% Plot NEW comprehensive ROC curves
generateNewROCPRCPlots(RnaROCObj, RxnROCObj, RnaPRCObj, RxnPRCObj);

%% Panel 3 — t-SNE reaction space plots

fprintf('Generating Panel 3: t-SNE reaction space analysis...\n');

% Create t-SNE plots
generateTSNEReactionPlots(RxnStatsByIndv, RnaStatsByIndv, drugMap, regresTable);

%% Panel 4 - Ab Initio Result
generateIrreversibilityComparisonPlot()

%% Panel 5 — Drug Similarity Correlograms (Reaction vs RNA features)
fprintf('Generating Panel 5: Correlograms for drug similarity analysis...\n');

% Generate correlograms comparing drug clustering
generateDrugCorrelogramComparison(RxnStatsMedMed, RnaStatsMedMed, drugMap, regresTable);

%% Save completion status

fprintf('All publication figures generated successfully!\n');
save(fullfile('figures','generation_complete.mat'), 'metrics', 'iRange', 'cMax');

%% House-keeping

%clear; clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL HELPER FUNCTIONS - All self-contained, no external dependencies
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generateOptimizedPubPlots - Create publication quality box/violin plots

function generateOptimizedPubPlots(Rxn, Rna, metrics)

% Configuration
metricNames = ["_{ }Accuracy_{ }"	"Recall_{macro}"	"Precision_{macro}"	"Sensitivity_{macro}"	"Specificity_{macro}"	"F1_{macro}"	"_{ }AUROC_{ }"	"AUPRC_{macro}", ...
    "Recall"	"Precision"	"Sensitivity"	"Specificity"	"F1"	"AUROC"	"AUPRC", ...
    "Recall"	"Precision"	"Sensitivity"	"Specificity"	"F1"	"AUROC"	"AUPRC"];
metricNums = {[1,6,3,2,7,8]; [13, 10, 9, 15]; [20, 17, 16, 22]};

% Color palette
blue_color = [0 0.4470 0.7410];
orange_color = [0.8500 0.3250 0.0980];

% Font settings
font_name = 'Arial';
title_size = 11;
label_size = 14;
text_size = 9;

% Generate first 3 optimized figures
for outerA = 1:length(metricNums)

    % Create figure with Pub Communications dimensions
    fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 12], 'Color','white');
    set(fig, 'DefaultTextFontName', font_name);
    set(fig, 'DefaultTextFontSize', text_size);
    set(fig, 'DefaultAxesFontName', font_name);
    set(fig, 'DefaultAxesFontSize', label_size);

    tlpltO = tiledlayout(4, length(metricNums{outerA}), "TileSpacing", "compact", ...
        "Padding", "compact", "TileIndexing", "columnmajor");
    numStats = length(metricNums{outerA});

    for a = metricNums{outerA}

        % Prepare data
        RxnT = mean(Rxn.(metrics{a})(:,:), 2, 'omitnan');
        RnaT = mean(Rna.(metrics{a})(:,:), 2, 'omitnan');

        Rxn_all = Rxn.(metrics{a})(:);
        Rna_all = Rna.(metrics{a})(:);

        % Main comparison plot
        nexttile([3,1])

        % Box plots with optimized styling
        boxchart(ones(length(RnaT),1), RnaT, 'Notch', 'on', ...
            'BoxFaceColorMode', 'manual', 'BoxFaceColor', blue_color, ...
            'BoxEdgeColorMode', 'manual', 'BoxEdgeColor', 'k', ...
            'BoxMedianLineColorMode', 'manual', 'BoxMedianLineColor', 'k', ...
            'MarkerStyle', 'none', 'BoxFaceAlpha', 1, 'LineWidth', 1, ...
            'WhiskerLineStyle', '-')
        hold on

        boxchart(2*ones(length(RxnT),1), RxnT, 'Notch', 'on', ...
            'BoxFaceColorMode', 'manual', 'BoxFaceColor', orange_color, ...
            'BoxEdgeColorMode', 'manual', 'BoxEdgeColor', 'k', ...
            'BoxMedianLineColorMode', 'manual', 'BoxMedianLineColor', 'k', ...
            'MarkerStyle', 'none', 'BoxFaceAlpha', 1, 'LineWidth', 1, ...
            'WhiskerLineStyle', '-')

        % Add jittered points with improved visibility
        n = length(Rxn_all);
        jitter_amount = 0.3;
        density_jitter = (rand(n,1) - 0.5) * jitter_amount;
        x1 = ones(n,1) + density_jitter;
        x2 = 2*ones(n,1) + density_jitter;

        scatter(x1, Rna_all, 25, blue_color, 'filled', 'MarkerFaceAlpha', 0.03);
        scatter(x2, Rxn_all, 25, orange_color, 'filled', 'MarkerFaceAlpha', 0.03);

        % Connect paired points with improved styling
        for i = 1:n
            plot([x1(i), x2(i)], [Rna_all(i), Rxn_all(i)], '-', 'Color', [0 0 0 0.05], 'LineWidth', 0.5);
        end

        % Redraw boxcharts on top
        boxchart(ones(length(RnaT),1), RnaT, 'Notch', 'on', ...
            'BoxFaceColorMode', 'manual', 'BoxFaceColor', blue_color, ...
            'BoxEdgeColorMode', 'manual', 'BoxEdgeColor', 'k', ...
            'BoxMedianLineColorMode', 'manual', 'BoxMedianLineColor', 'k', ...
            'MarkerStyle', 'none', 'BoxFaceAlpha', 1, 'LineWidth', 1, ...
            'WhiskerLineStyle', '-')

        boxchart(2*ones(length(RxnT),1), RxnT, 'Notch', 'on', ...
            'BoxFaceColorMode', 'manual', 'BoxFaceColor', orange_color, ...
            'BoxEdgeColorMode', 'manual', 'BoxEdgeColor', 'k', ...
            'BoxMedianLineColorMode', 'manual', 'BoxMedianLineColor', 'k', ...
            'MarkerStyle', 'none', 'BoxFaceAlpha', 1, 'LineWidth', 1, ...
            'WhiskerLineStyle', '-')

        % Perform statistical analysis
        diff_data = RxnT - RnaT;
        [nbp, ~, obs_diff, ~, ci_boot] = permutation_test_with_bootstrap(RxnT, RnaT);
        nbp = double(nbp);

        % Professional statistical annotation
        if nbp < 0.001/numStats
            sig_note = "***";
        elseif nbp < 0.01/numStats
            sig_note = "**";
        elseif nbp < 0.05/numStats
            sig_note = "*";
        else
            sig_note = "ns";
        end

        text(1.5, 0.95, sig_note, 'HorizontalAlignment', 'center', ...
            'FontName', font_name, 'FontSize', text_size, 'FontWeight', 'normal');

        % Professional effect size reporting
        effect_text = sprintf('ES: %+.2f\n[%.2f, %.2f]', ...
            obs_diff, ci_boot(1), ci_boot(2));
        text(0.1, 0.01, effect_text, 'Units', 'normalized', ...
            'FontName', font_name, 'FontSize', text_size, 'Color', 'k', ...
            'HorizontalAlignment','left', 'VerticalAlignment','bottom');

        % Formatting
        box on
        xlabel(sprintf('%s', metricNames(a)), 'FontName', font_name, 'FontSize', title_size);
        xticks([1 2])
        xticklabels([])
        xtickangle(45)
        ylim([0 1])
        yticks(0:0.1:1)
        ax = gca;
        ax.YAxis.FontSize = label_size-1;
        ax.XAxis.TickLength = [0 0];
        ax.YAxis.TickLength = [0.015 0.025];

        if all(a ~= [1, 13, 20])
            yticklabels(["", "", "", "", "", "", "", "", "", "", ""])
        else
            ylabel("Metric value", 'FontName', font_name, 'FontSize', label_size+1);
        end

        if any(a == [8, 15, 22])
            bar(1,0)
            bar(2,0)
            lgd = legend({"RNA", "Rxn"}, 'FontSize', text_size, 'FontName', font_name, ...
                'Units', 'normalized', 'Position', [0.82 0.40 0.09 0.05]);
            legend('boxoff')
        end

        xlim([0.6 2.4])

        % Difference histogram
        nexttile()

        [counts, edges] = histcounts(diff_data, [-1:0.05:1]);
        bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

        % Color bars based on sign
        colors = zeros(length(bin_centers), 3);
        colors(bin_centers < 0, :) = repmat(blue_color, sum(bin_centers < 0), 1);
        colors(bin_centers > 0, :) = repmat(orange_color, sum(bin_centers > 0), 1);

        b = bar(bin_centers, counts./100, 'FaceColor', 'flat');
        b.BarWidth = 1;
        b.CData = colors;

        xlim([-0.45 0.45])
        xticks([-0.2:0.4:0.2]);
        xtickangle(0)
        ax = gca;
        ax.XAxis.FontSize = label_size-1;
        ylim([0 .45])
        yticks([0:.10:.40])
        ax = gca;
        ax.YAxis.FontSize = label_size-1;

        if all(a ~= [1, 13, 20])
            yticklabels(["", "", "", "", ""])
        else
            ylabel("Proportion", 'FontName', font_name, 'FontSize', label_size+1)
        end

        ax.TickLength = [0.04, 0.025];

        % Sign test
        p = signtest(diff_data);

        if p*numStats < 0.001
            pval_text = "***";
        elseif p*numStats < 0.01
            pval_text = "**";
        elseif p*numStats < 0.05
            pval_text = "*";
        else
            pval_text = "ns";
        end

        text(0.5, 1, pval_text, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontName', font_name, 'FontSize', text_size, 'Color', 'k');

        if any(a == [3, 10, 17])
            xlabel(" Paired difference (sign test)", 'FontName', font_name, 'FontSize', title_size)
        end

    end

    if outerA == 1
        title(tlpltO, "Model Statistics", 'FontName', font_name, 'FontSize', title_size+1)
    elseif outerA == 2
        title(tlpltO, "Class 0 statistics (non-toxic)")
    elseif outerA == 3
        title(tlpltO, "Class 1 statistics (toxic)")
    end

    % Export figure
    saveas(fig, fullfile('figures', sprintf('Figure_%d_metrics.png', outerA)));
    saveas(fig, fullfile('figures','svg_figs', sprintf('Figure_%d_metrics.svg', outerA)), 'svg');
    savefig(fig, fullfile('figures', sprintf('Figure_%d_metrics.fig', outerA)));

end

end

%% generateNewROCPRCPlots - NEW IMPROVED ROC/PRC visualization with your code

function generateNewROCPRCPlots(Rna_ROCObj, Rxn_ROCObj, Rna_PRCObj, Rxn_PRCObj)

% Create figure with performance comparison plots
fig = figure('Units', 'centimeters', 'Position', [-20, 2, 18, 16.5], 'Color','white');
tl1 = tiledlayout(3, 2, "TileSpacing", "compact", "Padding", "compact");

% Font settings
font_name = 'Arial';
title_size = 11;
label_size = 10;
text_size = 9;
set(fig, 'DefaultTextFontName', font_name);
set(fig, 'DefaultAxesFontName', font_name);
set(fig, 'DefaultAxesFontSize', label_size);

% Configuration parameters
FPR_grid = [0 0.25 1/3 0.5 2/3 0.75 1];

% Define plot configurations
plot_configs = {
    struct('type', 'ROC', ...
    'data_sources', {{'Rna_ROCObj', 'Rxn_ROCObj'}}, ...
    'field_names', {{'TruePositiveRate', 'FalsePositiveRate'}}, ...
    'title_prefix', '_{ }ROC curve_{ }', ...
    'real_values', [0.6456, 0.7303], ...
    'class', 1, ...
    'xlabel', 'False Positive Rate', ...
    'ylabel', 'True Positive Rate', ...
    'FPR_grid_short', [0 0.25 0.5 0.75 1], ...
    'FPR_grid_fine', [0 0.25-1e-7 0.25 0.5-1e-7 0.5 0.75-1e-7 0.75 1]);
    struct('type', 'PRC_class0', ...
    'data_sources', {{'Rna_PRCObj', 'Rxn_PRCObj'}}, ...
    'field_names', {{'PositivePredictiveValue', 'TruePositiveRate'}}, ...
    'title_prefix', '_{ }Class 0 PRC_{ }', ...
    'real_values', [0.6469, 0.7333], ...
    'class', 0, ...
    'xlabel', 'Recall', ...
    'ylabel', 'Precision', ...
    'FPR_grid_short', [0 0.25 1/3 0.5 2/3 0.75 1], ...
    'FPR_grid_fine', [0 0.25-1e-7 0.25 1/3-1e-7 1/3 0.5-1e-7 0.5 2/3-1e-7 2/3 0.75-1e-7 0.75 1]);
    struct('type', 'PRC_class1', ...
    'data_sources', {{'Rna_PRCObj', 'Rxn_PRCObj'}}, ...
    'field_names', {{'PositivePredictiveValue', 'TruePositiveRate'}}, ...
    'title_prefix', '_{  }Class 1 PRC_{  }', ...
    'real_values', [0.6878, 0.7516], ...
    'class', 1, ...
    'xlabel', 'Recall', ...
    'ylabel', 'Precision', ...
    'FPR_grid_short', [0 0.25 1/3 0.5 2/3 0.75 1], ...
    'FPR_grid_fine', [0 0.25-1e-7 0.25 1/3-1e-7 1/3 0.5-1e-7 0.5 2/3-1e-7 2/3 0.75-1e-7 0.75 1])
    };

color_idx = lines(2);
plotValues = cell(3,2);
auc_values = [];
err = [];
spots = [1, 3, 4, 2, 5];

for plot_idx = 1:length(plot_configs)
    config = plot_configs{plot_idx};
    nexttile(spots(plot_idx));

    % Generate plots for both RNA and Rxn
    for source_idx = 1:2
        [err(plot_idx, source_idx), auc_values(plot_idx, source_idx), plotValues{plot_idx, source_idx}] = generate_performance_plot(config, source_idx, config.FPR_grid_short);
    end

    % Set plot properties
    ylim([0 1]);
    xlim([0 1]);
    ax = gca;
    ax.YAxis.FontSize = label_size-1;
    ax.XAxis.FontSize = label_size-1;
    xlabel(config.xlabel, "FontName", font_name, "FontSize",  title_size);
    ylabel(config.ylabel, "FontName", font_name, "FontSize",  title_size);
    title(config.title_prefix, "FontName", font_name, "FontSize",  title_size);

    grid on;
    set(ax, 'TickLength', [0 0])

    % Add AUC values as text
    text(0.983, 0.03, sprintf('             RNA: %.2f\n            Rxn:  %.2f', config.real_values(1), config.real_values(2)),'EdgeColor','k','BackgroundColor','white','FontName',font_name, 'FontSize',text_size,'HorizontalAlignment','right', 'VerticalAlignment','bottom');
    rectangle('Position', [0.62, 0.165, 0.14, 0.03], 'FaceColor', color_idx(1,:), 'EdgeColor','none');
    rectangle('Position', [0.62, 0.06, 0.14, 0.03], 'FaceColor', color_idx(2,:), 'EdgeColor','none');
end

% Macro PRC average plot
nexttile(spots(4));
for a = 1:2
    X_fine = plotValues{2,1}.Xval;
    mean_curve = mean([plotValues{2,a}.Yval; plotValues{3,a}.Yval],1);
    std_curve = mean([plotValues{2,a}.er; plotValues{3,a}.er],1);
    plot_performance_curve(X_fine, mean_curve, std_curve, a)
end

ylim([0 1]);
xlim([0 1]);
ax = gca;
ax.YAxis.FontSize = label_size-1;
ax.XAxis.FontSize = label_size-1;
xlabel(config.xlabel, "FontName", font_name, "FontSize",  title_size);
ylabel(config.ylabel, "FontName", font_name, "FontSize",  title_size);
title("PRC_{macro}", "FontName", font_name, "FontSize",  title_size);
grid on;
set(ax, 'TickLength', [0 0])

% Add AUC values as text
text(0.983, 0.03, sprintf('             RNA: %.2f\n            Rxn:  %.2f', 0.6663, 0.7442), 'EdgeColor','k','BackgroundColor','white','FontName',font_name, 'FontSize',text_size,'HorizontalAlignment','right', 'VerticalAlignment','bottom');
rectangle('Position', [0.62, 0.165, 0.14, 0.03], 'FaceColor', color_idx(1,:), 'EdgeColor','none');
rectangle('Position', [0.62, 0.06, 0.14, 0.03], 'FaceColor', color_idx(2,:), 'EdgeColor','none');

% Text summary plot
nexttile(spots(5));
axis off

textPlot1 = sprintf('AUROC: trapz vs (actual)   (mean error = %.1g)\nRNA: %.2f (%.2f)   Rxn:  %.2f (%.2f)\n', mean(abs(err(1,:))), auc_values(1,1), plot_configs{1}.real_values(1), auc_values(1,2), plot_configs{1}.real_values(2));
textPlot2 = sprintf('\n0-AUPRC: trapz vs (actual)   (mean error  = %.1g)\nRNA: %.2f (%.2f)   Rxn:  %.2f (%.2f)\n', mean(abs(err(2,:))), auc_values(2,1), plot_configs{2}.real_values(1), auc_values(2,2), plot_configs{2}.real_values(2));
textPlot3 = sprintf('\n1-AUPRC: trapz vs (actual)   (mean error  = %.1g)\nRNA: %.2f (%.2f)   Rxn:  %.2f (%.2f)\n', mean(abs(err(3,:))), auc_values(3,1), plot_configs{3}.real_values(1), auc_values(3,2), plot_configs{3}.real_values(2));
textPlot4 = sprintf('\nmac-AUPRC: trapz vs (actual)   (mean error  = %.1g)\nRNA: %.2f (%.2f)   Rxn:  %.2f (%.2f)\n', mean(abs(err(2:3,:)),"all"), mean(auc_values(2:3,1)), 0.6719, mean(auc_values(2:3,2)), 0.7362);

text(0, 1, [textPlot1, textPlot4, textPlot2, textPlot3], 'FontName',font_name, 'FontSize',text_size, 'HorizontalAlignment','left', 'VerticalAlignment','baseline');

% Save main ROC figure
saveas(fig, fullfile('figures', 'NEW_ROC_PRC_curves.png'));
saveas(fig, fullfile('figures','svg_figs', 'NEW_ROC_PRC_curves.svg'), 'svg');
savefig(fig, fullfile('figures', 'NEW_ROC_PRC_curves.fig'));

end

%% generateTSNEReactionPlots - Create t-SNE visualization of reaction space

function generateTSNEReactionPlots(RxnStats, RnaStats, drugMap, regresTable)

regresTable(regresTable.sigTF == 0, :) = [];

colors = [
    0.000, 0.000, 0.000;  % Black
    0.412, 0.412, 0.412;  % Dark gray
    0.533, 0.533, 0.533;  % Medium gray
    0.600, 0.600, 0.600;  % Gray
    0.647, 0.165, 0.165;  % Brown
    0.545, 0.000, 0.000;  % Dark red
    0.400, 0.067, 0.000;  % Dark brown
    0.835, 0.369, 0.000;  % Vermillion
    1.000, 0.647, 0.000;  % Dark orange
    0.902, 0.624, 0.000;  % Orange
    0.867, 0.800, 0.467;  % Light yellow
    0.933, 0.867, 0.510;  % Khaki
    0.941, 0.894, 0.259;  % Yellow
    0.600, 0.600, 0.200;  % Olive
    0.000, 0.392, 0.000;  % Forest green
    0.067, 0.467, 0.200;  % Dark green
    0.000, 0.620, 0.451;  % Bluish green
    0.267, 0.667, 0.600;  % Teal
    0.184, 0.310, 0.310;  % Dark slate gray
    0.000, 0.502, 0.502;  % Dark teal
    0.533, 0.800, 0.933;  % Light blue
    0.000, 0.447, 0.698;  % Blue
    0.337, 0.706, 0.914;  % Sky blue
    0.690, 0.769, 0.871;  % Light steel blue
    0.400, 0.400, 0.800;  % Light purple
    0.200, 0.157, 0.533;  % Dark blue
    0.502, 0.000, 0.502;  % Dark purple
    0.667, 0.267, 0.600;  % Purple
    0.800, 0.475, 0.655;  % Reddish purple
    0.533, 0.157, 0.333;  % Maroon
    0.800, 0.400, 0.467;  % Pink
    ];


for a = [1,2]
    if a ==1
        UseStats = RxnStats;
    elseif a ==2
        UseStats = RnaStats;
    end

    [X,y] = makeXYtable(UseStats,drugMap, regresTable, a);
    Xtemp = fillmissing(X{:,:}, 'constant', 0);
    X{:,:} = Xtemp;
    y.Drug = categorical(extractBefore(string(y.Properties.RowNames), 4));
    y.Name = strings(height(y),1);
    for b = 1:height(y)
        tempTHREELET = string(y.Drug(b));
        tempName = convertStringsToChars(lower(string(drugMap.Drug(string(drugMap.ThreeLet) == tempTHREELET))));
        tempName(1) = upper(tempName(1));
        y.Name(b) = string(tempName);
        y.threeLet(b) = string(tempTHREELET);
    end
    %y.ROR = 2.^(y.ROR);
    [v, I] = sort(y.ROR, 'ascend');
    y = y(I, :);
    X = X(I,:);

    group = y.Drug;
    names = y.Name;
    uniqueGroups = unique(group, 'stable');
    uniqueNames = unique(names, 'stable');
    n_groups = length(uniqueGroups);
    numGroups = zeros(height(group),1);
    %sizeSym = (y.ROR - min(y.ROR)+0.5).*10;
    sizeSym = numGroups.*0 + 10;
    clr = colors(1:n_groups, :);
    %clr = turbo(n_groups);
    % symbols = 'o+*xsd^v><ph';  % 13 different symbols available in MATLAB
    symbols = 'o<sd>h^pv';
    sym = repmat(symbols, 1, ceil(n_groups/length(symbols)));
    sym = sym(1:n_groups);  % Take exactly 31 symbols
    numClr = zeros(height(group),3);
    numSym = zeros(height(group),1);
    for b = 1:length(uniqueGroups)
        groupIDX = find(group == uniqueGroups(b));
        numGroups(groupIDX) = b;
        for c = 1:length(groupIDX)
            numClr(groupIDX(c),:) = clr(b,:);
        end
    end




    rng default
    tsneMat = tsne(X{:,:}, 'Standardize', true, "Algorithm","exact", 'Distance', 'euclidean');



    xminLim = min(min(tsneMat(:,1)))-1;
    xmaxLim = max(max(tsneMat(:,1)))+1;
    yminLim = min(min(tsneMat(:,2)))-1;
    ymaxLim = max(max(tsneMat(:,2)))+1;
    figure('color','white');

    % X-axis marginal plot (top)
    subplot(4,4,[5, 6]);
    % Bin tSNE1 values and calculate mean ROR
    n_bins = 10;
    tVal = tinv([0.025  0.975],height(y));
    tVal = tVal(2);
    %[~, edges] = histcounts(tsneMat(:,1), n_bins);
    %bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    quantile_points = linspace(0, 1, n_bins+1);
    edges = quantile(tsneMat(:,1), quantile_points);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    mean_ror_x = zeros(size(bin_centers));
    nf_ror_x = zeros(size(bin_centers));
    bin_data_x = cell(n_bins, 1);  % Store data for each bin

    for i = 1:length(bin_centers)
        if i == 1
            idx = tsneMat(:,1) <= edges(i+1);
        elseif i == length(bin_centers)
            idx = tsneMat(:,1) > edges(i);
        else
            idx = tsneMat(:,1) > edges(i) & tsneMat(:,1) <= edges(i+1);
        end

        bin_data_x{i} = y.ROR(idx);

        if sum(idx) > 0
            mean_ror_x(i) = mean(y.ROR(idx));
            nf_ror_x(i) = tVal.*std(y.ROR(idx))./sqrt(sum(idx));
        elseif i > 1
            mean_ror_x(i) = mean_ror_x(i-1);
            nf_ror_x(i) = nf_ror_x(i-1);
        end
    end

    % LINEAR REGRESSION FOR X-AXIS
    mdl_x = fitlm([1:n_bins]', mean_ror_x', "linear");
    r2_x = mdl_x.Rsquared.Adjusted;
    p_x = mdl_x.Coefficients.pValue(2); % p-value for slope coefficient
    mdl_x = fitlm(bin_centers', mean_ror_x', "linear");
    x_range = linspace(min(bin_centers), max(bin_centers), 100);
    y_pred_x = predict(mdl_x, x_range');
    if p_x < 0.05
        plot(x_range, y_pred_x, 'r--', 'LineWidth', 1);
        hold on;
    end
    errorbar(bin_centers, mean_ror_x,nf_ror_x,'k', 'LineStyle','-', 'Marker','.');
    hold on
    if p_x <= 0.05
        text(0.05, 0.95, sprintf('R² = %.2f\np = %.2g', r2_x, p_x), ...
            'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'none', ...
            'VerticalAlignment', 'top');
    end
    % Perform t-tests for X-axis bins
    for i = 1:n_bins
        current_bin_data = bin_data_x{i};
        if length(current_bin_data) > 1
            % Test vs left bins (bins 1 to i-1)
            if i > 1
                left_data = [];
                for j = 1:(i-1)
                    left_data = [left_data; bin_data_x{j}];
                end
                if length(left_data) > 1
                    [~, p_left] = ttest2(current_bin_data, left_data);
                    if p_left < 0.05./(4.*n_bins)
                        text(bin_centers(i), mean_ror_x(i) + nf_ror_x(i) + 0.4, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
                    end
                end
            end
            % Test vs right bins (bins i+1 to n_bins)
            if i < n_bins
                right_data = [];
                for j = (i+1):n_bins
                    right_data = [right_data; bin_data_x{j}];
                end
                if length(right_data) >= 2
                    [~, p_right] = ttest2(current_bin_data, right_data);
                    if p_right < 0.05./(4.*n_bins)
                        text(bin_centers(i), mean_ror_x(i) + nf_ror_x(i) + 0.2, '#', 'FontSize', 10, 'HorizontalAlignment', 'center');
                    end
                end
            end
        else
        end
    end

    xlim([xminLim xmaxLim]);
    ylim([-0.5 2])
    ylabel('Mean log_{2}(ROR)')
    box off
    set(gca, 'XTickLabel', []);  % Remove x-axis labels to align with main plot
    hold off

    % Y-axis marginal plot (right)
    subplot(4,4,[11, 15]);
    % Bin tSNE2 values and calculate mean ROR
    % [~, edges] = histcounts(tsneMat(:,2), n_bins);
    % bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    quantile_points = linspace(0, 1, n_bins+1);
    edges = quantile(tsneMat(:,2), quantile_points);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    mean_ror_y = zeros(size(bin_centers));
    nf_ror_y = zeros(size(bin_centers));
    bin_data_y = cell(n_bins, 1);  % Store data for each bin
    for i = 1:length(bin_centers)
        if i == 1
            idx = tsneMat(:,2) <= edges(i+1);
        elseif i == length(bin_centers)
            idx = tsneMat(:,2) > edges(i);
        else
            idx = tsneMat(:,2) > edges(i) & tsneMat(:,2) <= edges(i+1);
        end

        bin_data_y{i} = y.ROR(idx);

        if sum(idx) > 0
            mean_ror_y(i) = mean(y.ROR(idx));
            nf_ror_y(i) = tVal.*std(y.ROR(idx))./sqrt(sum(idx));
        elseif i > 1
            mean_ror_y(i) = mean_ror_y(i-1);
            nf_ror_y(i) = nf_ror_y(i-1);
        end
    end


    mdl_y = fitlm([1:n_bins]', mean_ror_y', "linear");
    r2_y = mdl_y.Rsquared.Adjusted;
    p_y = mdl_y.Coefficients.pValue(2); % p-value for slope coefficient
    mdl_y = fitlm(bin_centers', mean_ror_y', "linear");
    y_range = linspace(min(bin_centers), max(bin_centers), 100);
    x_pred_y = predict(mdl_y, y_range');

    if p_y <= 0.05
        % Plot regression line first (behind errorbar)
        plot(x_pred_y, y_range, 'r--', 'LineWidth', 1);
        hold on;
    end
    errorbar(mean_ror_y, bin_centers,nf_ror_y,'horizontal','k','LineStyle','-', 'Marker','.');
    hold on
    if p_y < 0.05
        text(0.65, 0.15, sprintf('R² = %.2f\np = %.2g', r2_y, p_y), ...
            'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white', ...
            'VerticalAlignment', 'top');
    end
    % Perform t-tests for Y-axis bins
    for i = 1:n_bins
        current_bin_data = bin_data_y{i};
        if length(current_bin_data) >1

            % Test vs left bins (bins 1 to i-1) - lower bins in Y direction
            if i > 1
                left_data = [];
                for j = 1:(i-1)
                    left_data = [left_data; bin_data_y{j}];
                end
                if length(left_data) >= 2
                    [~, p_left] = ttest2(current_bin_data, left_data);
                    if p_left < 0.05./(4.*n_bins)
                        text(mean_ror_y(i) + nf_ror_y(i) + 0.2, bin_centers(i), '#', 'FontSize', 10, 'HorizontalAlignment', 'center');
                    end
                end
            end

            % Test vs right bins (bins i+1 to n_bins) - upper bins in Y direction
            if i < n_bins
                right_data = [];
                for j = (i+1):n_bins
                    right_data = [right_data; bin_data_y{j}];
                end
                if length(right_data) >= 2
                    [~, p_right] = ttest2(current_bin_data, right_data);
                    if p_right < 0.05./(4.*n_bins)
                        text(mean_ror_y(i) + nf_ror_y(i) + 0.4, bin_centers(i), '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
                    end
                end
            end
        else
        end
    end

    ylim([yminLim ymaxLim]);
    xlim([-0.5 2])
    xlabel('Mean log_{2}(ROR)')
    box off
    set(gca, 'YTickLabel', []);  % Remove y-axis labels to align with main plot
    hold off

    % LEGEND SUBPLOT - Position 3 (top-right)
    subplot(4,4,7);
    % Create dummy plots for legend
    hold on;
    uniqueROR = unique(y.ROR, 'stable');
    gscatter(tsneMat(:,1).*NaN, tsneMat(:,2).*NaN, names, clr, sym, sizeSym, "off", "filled");
    legend(uniqueNames,'Location', 'eastoutside', 'NumColumns',3)
    legend('boxoff')
    axis off;  % Hide the axes
    hold off;

    % Main tSNE plot (center)
    subplot(4,4,[9, 10, 13, 14]);
    %[~, x_edges] = histcounts(tsneMat(:,1), n_bins);
    % [~, y_edges] = histcounts(tsneMat(:,2), n_bins);
    x_edges = quantile(tsneMat(:,1), quantile_points);
    y_edges = quantile(tsneMat(:,2), quantile_points);

    hold on;
    % Add vertical lines for tSNE1 bin edges (exclude first and last edge)
    for i = 2:(length(x_edges)-1)
        xline(x_edges(i), ':k', 'LineWidth', 1);
    end
    % Add horizontal lines for tSNE2 bin edges (exclude first and last edge)
    for i = 2:(length(y_edges)-1)
        yline(y_edges(i), ':k', 'LineWidth', 1);
    end

    gscatter(tsneMat(:,1), tsneMat(:,2), names, clr, sym, 7, "off", "filled")
    if a == 1
        title("Rxn data")
    elseif a == 2
        title("RNA data")
    end
    % set(gca, 'Color', [0.8 0.8 0.8]);
    xlim([xminLim xmaxLim]);
    ylim([yminLim, ymaxLim]);
    %legend(uniqueGroups, 'Location','bestoutside')
    box off
    ylabel("t-SNE_2")
    xlabel("t-SNE_1")


    subplot(4,4,[1, 2, 3]);
    resultTable = sortrows(regresTable, "ROR", "descend");
    forestPlot(log(resultTable.ROR), 1:height(resultTable), [log(resultTable.CIlow) log(resultTable.CIhigh)], resultTable, "log_{2}(ROR)", clr, sym, uniqueGroups)
    % Adjust subplot spacing
    set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 .6 0.80 ]);

    % Save figure
    if a == 1
        saveas(gcf, fullfile('figures', 'tSNE_Rxn_analysis.png'));
        saveas(gcf, fullfile('figures','svg_figs', 'tSNE_Rxn_analysis.svg'), 'svg');
        savefig(gcf, fullfile('figures', 'tSNE_Rxn_analysis.fig'));
    elseif a == 2
        saveas(gcf, fullfile('figures', 'tSNE_RNA_analysis.png'));
        saveas(gcf, fullfile('figures','svg_figs', 'tSNE_RNA_analysis.svg'), 'svg');
        savefig(gcf, fullfile('figures', 'tSNE_RNA_analysis.fig'));
    end
end
end


%% Supporting statistical functions

function [p_value, perm_diff, obs_diff, ci_perm, ci_boot] = permutation_test_with_bootstrap(group1, group2)

% Permutation test with bootstrap confidence intervals
n_permutations = 10000;
n = length(group1);
differences = group1 - group2;
obs_diff = mean(differences);

% Permutation test
perm_diffs = zeros(n_permutations, 1);

for i = 1:n_permutations
    signs = 2 * (rand(n, 1) > 0.5) - 1;
    perm_differences = differences .* signs;
    perm_diffs(i) = mean(perm_differences);
end

perm_diff = mean(perm_diffs);
p_value = (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_permutations + 1);

% Bootstrap CI (BCa method)
try
    ci_boot = bootci(10000, {@mean, differences}, 'Type', 'bca', 'Alpha', 0.05/14);
catch
    % Fallback to percentile method if BCa fails
    ci_boot = bootci(10000, {@mean, differences}, 'Type', 'percentile', 'Alpha', 0.05/14);
end

ci_perm = []; % Not needed for permutation distribution

end

%% NEW ROC/PRC Supporting Functions

function [err, auc_value, plotValues] = generate_performance_plot(config, source_idx, X_grid)
% Generate performance plot for given configuration and data source

data_source_name = config.data_sources{source_idx};
X_fine = linspace(0, 1, 999);

% Storage for all performance curves
performance_curves = nan(400, length(X_fine));

% Extract data from all iterations
curve_idx = 1;
for iteration = 1:100
    for sub_iteration = 1:4
        try
            % Get data structure
            data_struct = evalin('caller', sprintf('%s{%d,%d}', ...
                data_source_name, iteration, sub_iteration));

            metrics = data_struct.Metrics;

            % Filter by class if specified
            if config.class == 0
                metrics = metrics(metrics.ClassName == 0, :);
            elseif config.class == 1
                metrics = metrics(metrics.ClassName == 1, :);
            end

            % Extract performance metrics
            Y_data = metrics.(config.field_names{1});
            X_data = metrics.(config.field_names{2});

            % Clean and process data
            [Y_clean, X_clean] = clean_performance_data(Y_data, X_data);

            % Sort by X values
            [X_sorted, sort_idx] = sort(X_clean, 'ascend');
            Y_sorted = Y_clean(sort_idx);

            % Apply boundary conditions
            [Y_final, X_final] = apply_boundary_conditions(Y_sorted, X_sorted, config.type);

            % Add small noise to prevent interpolation issues
            Y_final = Y_final + randn(length(Y_final), 1) * 1e-7;
            X_final = X_final + randn(length(X_final), 1) * 1e-7;

            % Interpolate to fine grid
            performance_curves(curve_idx, :) = interp1(X_final, Y_final, X_fine,'nearest', 'extrap');
            curve_idx = curve_idx + 1;
        catch
        end
    end
end

% Calculate statistics and plot
[mean_curve, std_curve] = calculate_binned_statistics(X_fine, performance_curves, X_grid, config.FPR_grid_fine);
plot_performance_curve(X_fine, mean_curve, std_curve, source_idx);

% Calculate AUC
valid_idx = ~isnan(mean_curve);
auc_value = trapz(X_fine(valid_idx), mean_curve(valid_idx));
err = abs(config.real_values(source_idx) - auc_value(end));

plotValues.Xval = X_fine;
plotValues.Yval = mean_curve;
plotValues.er = std_curve;

end

function [Y_clean, X_clean] = clean_performance_data(Y_data, X_data)
% Remove missing values and zeros from performance data

% Handle missing X values
Y_clean = Y_data;
X_clean = X_data;
Y_clean(ismissing(X_clean)) = 0;
X_clean(ismissing(X_clean)) = 0;

% Remove invalid entries
invalid_idx = ismissing(Y_clean) | Y_clean == 0;
Y_clean(invalid_idx) = [];
X_clean(invalid_idx) = [];

end

function [Y_final, X_final] = apply_boundary_conditions(Y_data, X_data, plot_type)
% Apply appropriate boundary conditions based on plot type

switch plot_type
    case 'ROC'
        % ROC curves: ensure (0,0) to (1,1) boundary
        X_final = [X_data];
        Y_final = [Y_data];
    case {'PRC_class0', 'PRC_class1'}
        % PRC curves: start at (0,1)
        X_final = [X_data];
        Y_final = [Y_data];
    otherwise
        error('Unknown plot type: %s', plot_type);
end

end

function [mean_output, std_output] = calculate_binned_statistics(X, Y, X_grid, X_interp_grid)
% Calculate binned statistics for performance curves

X_fine = linspace(0, 1, 999);

% Bin the data
binned_mean = nan(1, length(X_interp_grid));
binned_std = nan(1, length(X_interp_grid));

for bin_idx = 2:length(X_grid)
    if bin_idx == length(X_grid)
        in_bin = X <= X_grid(bin_idx) & X >= X_grid(bin_idx-1);
    else
        in_bin = X < X_grid(bin_idx) & X >= X_grid(bin_idx-1);
    end

    bin_start = (bin_idx-1)*2-1;
    bin_end = (bin_idx-1)*2;
    binned_mean(bin_start:bin_end) = mean(Y(:, in_bin), 'all', 'omitmissing');
    binned_std(bin_start:bin_end) = std(Y(:, in_bin), [], 'all', 'omitmissing');
end

% Interpolate to fine grid
mean_output = interp1(X_interp_grid, binned_mean, X_fine, 'previous', 'extrap');
std_output = interp1(X_interp_grid, binned_std, X_fine, 'previous', 'extrap');

% Remove NaN values
valid_idx = ~isnan(mean_output);
mean_output = mean_output(valid_idx);
std_output = std_output(valid_idx);

end

function plot_performance_curve(X_fine, mean_curve, std_curve, source_idx)
% Plot performance curve with confidence intervals

% Calculate confidence intervals
alpha = 0.05;
df = 400;
t_critical = abs(tinv(1 - alpha/14, df));
margin_error = t_critical * std_curve / 20;
ci_upper = min(mean_curve + margin_error, 1);
ci_lower = max(mean_curve - margin_error, 0);

color_idx = lines(2);

% Plot confidence interval and mean curve
fill([X_fine fliplr(X_fine)], [ci_upper fliplr(ci_lower)], color_idx(source_idx,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(X_fine, mean_curve, 'LineWidth', 2, 'Color', color_idx(source_idx,:));

end

%% Data processing functions

function [X, y] = makeXYtable(UseStats, drugMap, regresTable, type)

% Convert statistics to ML-ready feature matrix and labels
dataTable = table();
pathTableT_rowNames = string(UseStats.MedianDiff.Properties.VariableNames);
pathTableT_varNames = string(UseStats.MedianDiff.Properties.RowNames);

% Transpose and assign proper variable names
if type == 6
    dataTable{:,:} = UseStats.MedianDiff{:,:}';
elseif type == 9
    dataTable{:,:} = UseStats.FCDrug{:,:}';
end
dataTable.Properties.VariableNames = pathTableT_varNames;
dataTable.Properties.RowNames = pathTableT_rowNames;

regresTable.Drug = upper(string(regresTable.Drug));

dataTable.sigTF = zeros(height(dataTable),1).*NaN;

allThreeLet = extractBefore(pathTableT_rowNames, 4)';

for a = 1:height(regresTable)
    threeLet = drugMap.ThreeLet(drugMap.Drug == regresTable.Drug(a));
    drugSpots = contains(allThreeLet, string(threeLet));
    assert(isscalar(unique(allThreeLet(drugSpots))) | isempty(unique(allThreeLet(drugSpots))));
    dataTable.sigTF(drugSpots) = regresTable.sigTF(a);
    dataTable.ROR(drugSpots) = log2(regresTable.ROR(a));
end

dataTable(isnan(dataTable.sigTF), :) = [];

X = dataTable(:,1:end-2);
X{:,:} = normalize(X{:,:});
y = dataTable(:,end-1:end);
y.sigTF(y.sigTF == -1) = 0;

end

function forestPlot(xA, yA, ciA, resultsPLOT, type, clr, sym, uniqueGroups)
% Visualization
hold on;


toxIDX = find(resultsPLOT.sigTF == 1);
nontoxIDX = find(resultsPLOT.sigTF == -1);
maxToxIDX = max(toxIDX);
minnontoxIDX = min(nontoxIDX);
cutIDX = mean([maxToxIDX, minnontoxIDX]);
plot([cutIDX, cutIDX], [-2 2], '--k')
hold on
cutROR =  mean([log(resultsPLOT.ROR(maxToxIDX)), log(resultsPLOT.ROR(minnontoxIDX))]);
plot([0 height(resultsPLOT)+1], [cutROR, cutROR], '--k');

for i = 1:height(resultsPLOT)
    x = xA(i);
    y = yA(i);
    ci = ciA(i,:);

    errorbar(y,x, ci(1)-x, ci(2)-x, 'Color','k', 'LineWidth',1, 'LineStyle', 'none');
    plot(y, x,sym(end-i+1),'MarkerFaceColor',clr(end-i+1,:), 'MarkerEdgeColor','k');
end

set(gca,'XDir','reverse','XTick',1:height(resultsPLOT),...
    'XTickLabel',string(uniqueGroups(end:-1:1)),...
    'YGrid','on');
ylabel(type);
xlim([0 height(resultsPLOT)+1])
title('Cardiotoxicity Risk');


%set(gca, 'Color', [0.8 0.8 0.8]);
end

function generateIrreversibilityComparisonPlot()
% generateIrreversibilityComparisonPlot
% Compare subsystem-level irreversibility between the original and
% optimized cardiac metabolic models and visualize the top 31 changes
% as a horizontal bar chart.
%
% File names, struct field names, and data loading remain unchanged.
%
% 2025-10-08  Cleaned & annotated for clarity.  No '%%' sections used.

% -------------------------------------------------------------------------
% Load models and reaction lists
% -------------------------------------------------------------------------
fprintf('Loading model data...\n');
load(fullfile('models', 'HeartModel.mat'));        % originalModel → heart_model_curation
originalModel = heart_model_curation;              %#ok<NASGU>
load(fullfile('models', 'objRxns.mat'));           % allRxn list
originalModel.allRxn = allRxn;
load(fullfile('out', 'iCardio_optimized.mat'), 'ModelOpt');

% -------------------------------------------------------------------------
% Identify irreversible reactions (either lb or ub is zero OR same sign)
% -------------------------------------------------------------------------
origIrrev = (sign(originalModel.ub) == sign(originalModel.lb)) | ...
    (sign(originalModel.ub) == 0) | (sign(originalModel.lb) == 0);
optIrrev  = (sign(ModelOpt.ub) == sign(ModelOpt.lb)) | ...
    (sign(ModelOpt.ub) == 0)   | (sign(ModelOpt.lb) == 0);

% Gene-associated reaction indices
origGene = ~cellfun(@isempty, originalModel.rules);
optGene  = ~cellfun(@isempty, ModelOpt.rules);
optGene(end-1:end) = [];   % last two entries are memos in iCardio

% Subsystem names (remove placeholders)
optSubs  = string(ModelOpt.subSystems);
optSubs(optSubs == "Purine Deg") = [];
origSubs = string(originalModel.subSystems);

% Unique, non-empty subsystems
uSubs = unique(optSubs);
uSubs = uSubs(~strcmp(uSubs, ""));

% -------------------------------------------------------------------------
% Compute irreversibility percentages and reaction counts
% -------------------------------------------------------------------------
n = numel(uSubs);
origPct     = zeros(n,1);
optPct      = zeros(n,1);
origIrCt    = zeros(n,1);
optIrCt     = zeros(n,1);
origTotCt   = zeros(n,1);
optTotCt    = zeros(n,1);

for k = 1:n
    name = uSubs{k};

    idxO = strcmp(origSubs, name) & origGene;
    idxP = strcmp(optSubs , name) & optGene;

    if any(idxO)
        origIrCt(k)  = sum(origIrrev(idxO));
        origTotCt(k) = sum(idxO);
        origPct(k)   = 100 * origIrCt(k) / origTotCt(k);
    end
    if any(idxP)
        optIrCt(k)   = sum(optIrrev(idxP));
        optTotCt(k)  = sum(idxP);
        optPct(k)    = 100 * optIrCt(k)  / optTotCt(k);
    end
end

% -------------------------------------------------------------------------
% Fisher exact test (one-sided: optimized > original)
% -------------------------------------------------------------------------
p = nan(n,1);
for k = 1:n
    if origTotCt(k) && optTotCt(k)
        tbl = [origIrCt(k), origTotCt(k) - origIrCt(k); ...
            optIrCt(k) , optTotCt(k) - optIrCt(k)];
        [~, p(k)] = fishertest(tbl, 'Tail','left');
    end
end

% -------------------------------------------------------------------------
% Sort by Δ% and keep top 31
% -------------------------------------------------------------------------
delta   = optPct - origPct;
[delta, sIdx] = sort(delta,'descend');
origPct = origPct(sIdx); optPct = optPct(sIdx);
uSubs   = uSubs(sIdx);   p      = p(sIdx);
optTotCt = optTotCt(sIdx);

keep = min(31, numel(uSubs));
origPct = origPct(1:keep); optPct = optPct(1:keep);
uSubs   = uSubs(1:keep);   p      = p(1:keep);
optTotCt = optTotCt(1:keep);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fprintf('Creating plot...\n');
fig = figure('Units','centimeters','Position',[2 2 8 8],'Color','w'); hold on;

fontN = 'Arial'; set(fig,'DefaultAxesFontName',fontN,'DefaultTextFontName',fontN);

clr = lines(2); optC = clr(2,:); origC = clr(1,:);

% Horizontal bars (optimized on left for overlay clarity)
b1 = barh(1:keep, optPct, 'FaceColor', optC , 'EdgeColor','k', 'LineWidth',1);
b2 = barh(1:keep, origPct,'FaceColor', origC,'EdgeColor','k', 'LineWidth',1, ...
    'FaceAlpha',0.7);


yticks(1:keep);
yticklabels(labelCleaner(uSubs, ["metabolism", "biosynthesis"]));
ax = gca;
ax.FontSize = 7;
set(gca,'YDir','reverse'); % top = largest Δ
xlim([0 100]); ylim([0.5 keep+0.5]); grid on; box on;
xlabel('% Reactions Irreversible','FontSize',10);
title('Irreversibility by Subsystem – Optimized vs Original','FontSize',10);

lgd = legend([b1, b2], {'Optimized','Original'}, ...
    'Location','southeast','FontSize',6);

% -------------------------------------------------------------------------
% Significance stars (x ≈ 5%) and pathway size (x ≈ 20%)
% -------------------------------------------------------------------------
for k = 1:keep
    % Stars
    if ~isnan(p(k)) && p(k) < 0.05
        stars = '*'; if p(k) < 0.001, stars = '***';
        elseif p(k) < 0.01,  stars = '**'; end
        text(8, k+0.3, stars, 'FontSize',11, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Color', 'w');
    end

    % Pathway size
    txt = sprintf('(%d)', optTotCt(k));
    text(30, k, txt, 'FontSize',5, 'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', 'EdgeColor','none', 'Color','w');
end

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
saveas(fig, fullfile('figures','Panel4_irreversibility_comparison.png'));
savefig(fig, fullfile('figures','Panel4_irreversibility_comparison.fig'));
fprintf('Done.\n');
end
%%
function cleaned = labelCleaner(strArray, wordsToRemove)
% Pre-convert wordsToRemove to lowercase for efficiency
wordsToRemoveLower = lower(wordsToRemove);

% Initialize output
cleaned = strings(size(strArray));

% Loop over each label
for i = 1:numel(strArray)
    % Split on whitespace, collapsing multiples and trimming
    parts = split(strArray(i));  % Returns 1xN string array

    % Identify parts to keep (case-insensitive match)
    keepMask = ~ismember(lower(parts), wordsToRemoveLower);

    % Join kept parts with single space (original casing preserved)
    cleaned(i) = join(parts(keepMask), " ");
end
end

%% Updated generateDrugCorrelogramComparison function

function generateDrugCorrelogramComparison(RxnStats, RnaStats, drugMap, regresTable)

% Font settings for publication quality
font_name = 'Arial';
title_size = 12;
label_size = 10;

%% ===== GET THE 31 SIGNIFICANT DRUGS (MATCHES t-SNE FILTERING) =====

% Get unique significant drugs from regresTable (same as makeXYtable)
regresTable.Drug = upper(string(regresTable.Drug));  % Match makeXYtable preprocessing
sigDrugs = unique(regresTable.Drug(regresTable.sigTF ~= 0));

% Map to 3-letter codes (exactly like makeXYtable)
sigThreeLet = string(missing);
for i = 1:length(sigDrugs)
    idx = find(strcmp(drugMap.Drug, sigDrugs(i)));
    if ~isempty(idx)
        sigThreeLet(i) = drugMap.ThreeLet{idx};
    end
end
sigThreeLet = sigThreeLet(~ismissing(sigThreeLet));  % Clean missing
nDrugs = length(sigThreeLet);

fprintf('  Found %d significant drugs (matching t-SNE filtering)\n', nDrugs);

%% ===== EXTRACT AND PREPARE DATA (FILTERED TO 31 DRUGS) =====

% Reaction data (RxnStats.MedianDiff - features x samples)
rxnVarNames = string(RxnStats.MedianDiff.Properties.VariableNames);
rxnThreeLet = extractBefore(rxnVarNames, 4);
rxnData = RxnStats.MedianDiff{:,:};

% RNA data (RnaStats.MeanDrug - assuming this is your field; change if FCDrug)
rnaVarNames = string(RnaStats.MedianDiff.Properties.VariableNames);
rnaThreeLet = extractBefore(rnaVarNames, 4);
rnaData = RnaStats.MedianDiff{:,:};

%% ===== COMPUTE MEDIAN PROFILES PER DRUG (ONLY 31 DRUGS) =====

rxnDrugProfiles = zeros(size(rxnData, 1), nDrugs);
rnaDrugProfiles = zeros(size(rnaData, 1), nDrugs);

for i = 1:nDrugs
    threeLet = sigThreeLet(i);

    % Reaction samples for this drug
    rxnSamples = strcmp(rxnThreeLet, threeLet);
    assert(sum(rxnSamples) == 1)
    rxnDrugProfiles(:, i) = rxnData(:, rxnSamples);

    % RNA samples for this drug
    rnaSamples = strcmp(rnaThreeLet, threeLet);
    assert(sum(rnaSamples) == 1)
    rnaDrugProfiles(:, i) = rnaData(:, rnaSamples);
    rnaDrugProfiles(ismissing(rnaDrugProfiles)) = 0;
end

%% ===== COMPUTE CORRELATION MATRICES =====
fprintf('  Computing drug-drug correlation matrices (31 drugs)...\n');

rxnCorr = corr(rxnDrugProfiles);
rnaCorr = corr(rnaDrugProfiles);

%% ===== CREATE DRUG LABELS FOR 31 DRUGS =====

drugLabels = cell(nDrugs, 1);
for i = 1:nDrugs
    threeLet = sigThreeLet(i);
    idx = find(strcmp(string(drugMap.ThreeLet), threeLet), 1);
    if ~isempty(idx)
        drugName = char(drugMap.Drug(idx));
        drugName(1) = upper(drugName(1));
        drugLabels{i} = sprintf('%s (%s)', drugName, threeLet);
    else
        drugLabels{i} = threeLet;
    end
end

%% ===== PLOT CORRELOGRAM 1: REACTION FEATURES =====
fprintf('  Creating Reaction-based correlogram (31 drugs)...\n');

fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 18], 'Color', 'white');
set(fig1, 'DefaultTextFontName', font_name);
set(fig1, 'DefaultAxesFontName', font_name);
set(fig1, 'DefaultAxesFontSize', label_size);

correlogram(rxnCorr, 'AxisLabels', drugLabels, 'Sorting', true);
sgtitle('Drug Similarity: Reaction Features (31 sig. drugs)', ...
    'FontName', font_name, 'FontSize', title_size + 2, 'FontWeight', 'bold');

saveas(fig1, fullfile('figures', 'Correlogram_Reaction_31drugs.png'));
saveas(fig1, fullfile('figures', 'svg_figs', 'Correlogram_Reaction_31drugs.svg'), 'svg');
savefig(fig1, fullfile('figures', 'Correlogram_Reaction_31drugs.fig'));

%% ===== PLOT CORRELOGRAM 2: RNA FEATURES =====
fprintf('  Creating RNA-based correlogram (31 drugs)...\n');

fig2 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 18], 'Color', 'white');
set(fig2, 'DefaultTextFontName', font_name);
set(fig2, 'DefaultAxesFontName', font_name);
set(fig2, 'DefaultAxesFontSize', label_size);

correlogram(rnaCorr, 'AxisLabels', drugLabels, 'Sorting', true);
sgtitle('Drug Similarity: RNA Features (31 sig. drugs)', ...
    'FontName', font_name, 'FontSize', title_size + 2, 'FontWeight', 'bold');

saveas(fig2, fullfile('figures', 'Correlogram_RNA_31drugs.png'));
saveas(fig2, fullfile('figures', 'svg_figs', 'Correlogram_RNA_31drugs.svg'), 'svg');
savefig(fig2, fullfile('figures', 'Correlogram_RNA_31drugs.fig'));

%% ===== SUMMARY STATISTICS =====
fprintf('\n  === Correlogram Summary (31 drugs only) ===\n');
fprintf('  Reaction features:\n');
fprintf('    Shape: %dx%d\n', size(rxnCorr));
fprintf('    Mean corr: %.3f\n', mean(rxnCorr(triu(true(size(rxnCorr)),1)), 'omitnan'));
fprintf('    Median corr: %.3f\n', median(rxnCorr(triu(true(size(rxnCorr)),1)), 'omitnan'));

fprintf('\n  RNA features:\n');
fprintf('    Shape: %dx%d\n', size(rnaCorr));
fprintf('    Mean corr: %.3f\n', mean(rnaCorr(triu(true(size(rnaCorr)),1)), 'omitnan'));
fprintf('    Median corr: %.3f\n', median(rnaCorr(triu(true(size(rnaCorr)),1)), 'omitnan'));

fprintf('\n  Correlograms saved (31 drugs only)!\n\n');
figure(fig1); corr_ax_rxn = gca;  % Rxn correlogram axes

% EXTRACT CLUSTERED ORDER DIRECTLY FROM CORRELOGRAM (100% match)
clustered_labels_rxn = string(corr_ax_rxn.XTickLabel);
threeLet_rxn = extractAfter(clustered_labels_rxn, '(');
threeLet_rxn = extractBefore(threeLet_rxn, ')');
threeLet_rxn = upper(threeLet_rxn);  % Safe upper

% Assign classes & colors
[classID_rxn, classNames_rxn] = assignDrugClass(threeLet_rxn);
classColors = [
    1.0, 0.2, 0.2;   % 1. Anthracyclines: Red
    0.2, 0.2, 1.0;   % 2. BCR-ABL TKIs: Blue
    0.2, 1.0, 0.2;   % 3. EGFR/HER: Green
    1.0, 0.6, 0.2;   % 4. Multi/VEGFR: Orange
    0.8, 0.2, 0.8;   % 5. Other TKIs: Magenta
    1.0, 0.5, 0.8;   % 6. mAbs: Pink
    0.6, 0.6, 0.2;   % 7. Misc: Olive
    0.5, 0.5, 0.5    % Other: Gray
    ];
marginalColors_rxn = classColors(classID_rxn, :);

% LEFT MARGINAL: Vertical class color bar (matches left axis)
ax_left_rxn = axes('Position', [0.035, 0.12, 0.025, 0.76]);  % Adjusted for fit
barh((1:nDrugs), ones(nDrugs,1)*0.6, 'FaceColor', 'flat', 'CData', marginalColors_rxn, ...
    'EdgeColor', 'none', 'BarWidth', 1);
set(ax_left_rxn, 'YDir', 'reverse', 'YTick', [], 'XTick', [], 'XLim', [0 1], ...
    'YLim', [0.5 nDrugs+0.5], 'Color', 'none');

% BOTTOM MARGINAL: Horizontal class color bar (matches bottom axis)
ax_bottom_rxn = axes('Position', [0.13, 0.035, 0.71, 0.025]);
bar((1:nDrugs), ones(nDrugs,1)*0.6, 'FaceColor', 'flat', 'CData', marginalColors_rxn, ...
    'EdgeColor', 'none', 'BarWidth', 1);
set(ax_bottom_rxn, 'XTick', [], 'YTick', [], 'XLim', [0.5 nDrugs+0.5], 'YLim', [0 1], 'Color', 'none');

% LEGEND (top-right)
legend_handles = arrayfun(@(i) patch([0 1 1 0], [0 0 1 1], classColors(i,:), ...
    'EdgeColor', 'none'), 1:7, 'Unif', 0);
lgd_rxn = legend([legend_handles{:}], classNames_rxn, ...
    'Position', [0.78, 0.75, 0.18, 0.22], 'Box', 'off', 'FontSize', 9);

% Re-save enhanced Rxn
saveas(fig1, fullfile('figures', 'Correlogram_Rxn_ClassMarginals.png'));
saveas(fig1, fullfile('figures', 'svg_figs', 'Correlogram_Rxn_ClassMarginals.svg'), 'svg');

%% ===== RNA (identical, exact from its correlogram) =====
figure(fig2); corr_ax_rna = gca;

clustered_labels_rna = string(corr_ax_rna.XTickLabel);
threeLet_rna = extractAfter(clustered_labels_rna, '(');
threeLet_rna = extractBefore(threeLet_rna, ')');
threeLet_rna = upper(threeLet_rna);

[classID_rna, ~] = assignDrugClass(threeLet_rna);
marginalColors_rna = classColors(classID_rna, :);

% LEFT MARGINAL RNA
ax_left_rna = axes('Position', [0.035, 0.12, 0.025, 0.76]);
barh((1:nDrugs), ones(nDrugs,1)*0.6, 'FaceColor', 'flat', 'CData', marginalColors_rna, ...
    'EdgeColor', 'none');
set(ax_left_rna, 'YDir', 'reverse', 'YTick', [], 'XTick', [], 'XLim', [0 1], ...
    'YLim', [0.5 nDrugs+0.5], 'Color', 'none');

% BOTTOM MARGINAL RNA
ax_bottom_rna = axes('Position', [0.13, 0.035, 0.71, 0.025]);
bar((1:nDrugs), ones(nDrugs,1)*0.6, 'FaceColor', 'flat', 'CData', marginalColors_rna, ...
    'EdgeColor', 'none');
set(ax_bottom_rna, 'XTick', [], 'YTick', [], 'XLim', [0.5 nDrugs+0.5], 'YLim', [0 1], 'Color', 'none');

% LEGEND RNA
lgd_rna = legend([legend_handles{:}], classNames_rxn, ...
    'Position', [0.78, 0.75, 0.18, 0.22], 'Box', 'off', 'FontSize', 9);

% Re-save enhanced RNA
saveas(fig2, fullfile('figures', 'Correlogram_RNA_ClassMarginals.png'));
saveas(fig2, fullfile('figures', 'svg_figs', 'Correlogram_RNA_ClassMarginals.svg'), 'svg');

end

function [classID, classNames] = assignDrugClass(threeLet)
% Table-like STRUCT: class names as fields, loop columns (drugs per class)
% Edit/add your EXACT 3-let codes here - super easy!

classDrugs(1).name = 'Anthracyclines';
classDrugs(1).codes = {'DOX', 'DAU', 'IDA', 'EPI'};

classDrugs(2).name = 'BCR-ABL TKI';
classDrugs(2).codes = {'IMA', 'DAS', 'NIL', 'BOS'};

classDrugs(3).name = 'EGFR/HER TKI';
classDrugs(3).codes = {'ERL', 'LAP', 'AFA'};

classDrugs(4).name = 'Multi/VEGFR TKI';
classDrugs(4).codes = {'SUN', 'SOR', 'REG', 'VAN', 'CAB'};

classDrugs(5).name = 'Other TKI';
classDrugs(5).codes = {'DAB', 'CER'};

classDrugs(6).name = 'mAbs';
classDrugs(6).codes = {'RIT', 'TRA', 'BEV', 'CET'};

classDrugs(7).name = 'Misc';
classDrugs(7).codes = {'TOF', 'BOR', 'CYC', 'DAC'};

classNames = {classDrugs.name};

% Convert to strings if needed
threeLet = upper(string(threeLet));
nDrugs = length(threeLet);
classID = 8 * ones(nDrugs, 1);  % 8 = Gray "Other/Unassigned"

% LOOP THROUGH EACH CLASS "COLUMN" (drugs per class)
for c = 1:length(classDrugs)
    drugs_c = upper(string(classDrugs(c).codes));
    match = ismember(threeLet, drugs_c);
    classID(match) = c;
end
end
