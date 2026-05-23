function summary = plot_penalty(optsIn)
%Plot penalty and conditioning data.

assert(exist('optsIn', 'var') == 1, 'plot_penalty requires an options structure.');

clc;
close all;
format short g;

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

exampleDir = fileparts(fileparts(mfilename('fullpath')));
rootDir = fullfile(exampleDir, 'data');
cd(rootDir);

Nc = get_opt(optsIn, 'Nc', 25);
t = get_opt(optsIn, 't', 1);
pdeg = 1 + t;
refineList = get_opt(optsIn, 'refineList', 1:7);
CsigmaList = get_opt(optsIn, 'CsigmaList', [10 20 30]);
sweepTag = sanitize_tag(get_opt(optsIn, 'sweepTag', join_numbers(CsigmaList)));

savePng = get_opt(optsIn, 'save_png', true);
savePdf = get_opt(optsIn, 'save_pdf', true);
pngDPI = get_opt(optsIn, 'pngDPI', 600);
plotOnly = logical(get_opt(optsIn, 'plot_only', false));
tauShift = get_opt(optsIn, 'tau_shift', 0.9 * 4.96999274451623);
epsDiag = get_opt(optsIn, 'eps_diag', 1e-12);
ifaceReg = get_opt(optsIn, 'iface_reg', 1e-12);
condMaxN = get_opt(optsIn, 'prec_cond_max_n', 20000);
condFullMaxN = get_opt(optsIn, 'cond_full_max_n', 1400);
condEigsTol = get_opt(optsIn, 'cond_eigs_tol', 1e-8);
condEigsMaxit = get_opt(optsIn, 'cond_eigs_maxit', 2000);
condOpts = struct( ...
    'cond_full_max_n', condFullMaxN, ...
    'cond_eigs_tol', condEigsTol, ...
    'cond_eigs_maxit', condEigsMaxit);

baseCandidates = get_opt(optsIn, 'baseCsigmaCandidates', [20 CsigmaList(:).']);

resultRoot = fullfile(rootDir, 'result', 'penalty');
if ~exist(resultRoot, 'dir')
    error('Result root not found: %s', resultRoot);
end

if isempty(sweepTag)
    fileSuffix = '';
    plotName = 'figures';
else
    fileSuffix = ['_' sweepTag];
    plotName = ['figures_' sweepTag];
end

summaryMat = fullfile(resultRoot, ['condition' fileSuffix '.mat']);
summaryCsv = fullfile(resultRoot, ['condition' fileSuffix '.csv']);
plotDir = fullfile(resultRoot, plotName);
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

if plotOnly
    if ~exist(summaryMat, 'file')
        error('Summary file not found: %s', summaryMat);
    end
    plot_condition_summary(summaryMat, struct( ...
        'outDir', plotDir, ...
        'save_png', savePng, ...
        'save_pdf', savePdf, ...
        'pngDPI', pngDPI));
    S = load(summaryMat, 'summary');
    summary = S.summary;
    fprintf('\n[PLOT ] redrew C_sigma condition figure from %s\n', summaryMat);
    fprintf('[SAVED] %s\n', plotDir);
    return;
end

rows = [];
for refine = refineList(:).'
    runFile = find_saved_run(resultRoot, refine, baseCandidates);
    S = load(runFile, 'run');
    run = S.run;
    fprintf('[LOAD] refine=%d %s\n', refine, runFile);

    for Csigma = CsigmaList(:).'
        [Ksigma, M, P] = reconstruct_Ksigma(run, Csigma);
        Ktau = Ksigma - tauShift * M;
        Ktau = 0.5 * (Ktau + Ktau');

        Pnone = build_condition_preconditioner_matrix_local( ...
            Ktau, P, 'none', epsDiag);
        condNone = shifted_generalized_abs_condition_local( ...
            Ktau, Pnone, condMaxN, 'none', condOpts);

        Pjacobi = build_condition_preconditioner_matrix_local( ...
            Ktau, P, 'purediag', epsDiag);
        condJacobi = shifted_generalized_abs_condition_local( ...
            Ktau, Pjacobi, condMaxN, 'Jacobi', condOpts);

        Pib = build_condition_preconditioner_matrix_local( ...
            Ktau, P, 'interfaceblock', epsDiag);
        condIB = shifted_generalized_abs_condition_local( ...
            Ktau, Pib, condMaxN, 'PITB-DG', condOpts);

        row = struct();
        row.Csigma = Csigma;
        row.refine = run.meta.Refinement;
        row.h = run.meta.h;
        row.cond_none = condNone;
        row.cond_jacobi = condJacobi;
        row.cond_interfaceblock = condIB;
        rows = [rows; row]; %#ok<AGROW>

        fprintf('[COND] C_sigma=%g refine=%d none=%.6e jacobi=%.6e PITB-DG=%.6e\n', ...
            Csigma, row.refine, condNone, condJacobi, condIB);
    end
end

T = struct2table(rows);
T = sortrows(T, {'Csigma', 'refine'});

summary = struct();
summary.Nc = Nc;
summary.pdeg = pdeg;
summary.t = t;
summary.refineList = refineList;
summary.CsigmaList = CsigmaList;
summary.quantity = ['kappa_tau(P)=sigma_max(P^{-1/2}(K_sigma-tau*M)P^{-1/2})/' ...
    'sigma_min(P^{-1/2}(K_sigma-tau*M)P^{-1/2})'];
summary.tau_shift = tauShift;
summary.eps_diag = epsDiag;
summary.iface_reg = ifaceReg;
summary.prec_cond_max_n = condMaxN;
summary.cond_full_max_n = condFullMaxN;
summary.cond_eigs_tol = condEigsTol;
summary.cond_eigs_maxit = condEigsMaxit;
summary.table = T;

save(summaryMat, 'summary');
writetable(T, summaryCsv);
plot_condition_summary(summaryMat, struct( ...
    'outDir', plotDir, ...
    'save_png', savePng, ...
    'save_pdf', savePdf, ...
    'pngDPI', pngDPI));

fprintf('\n[SAVED] %s\n', summaryMat);
fprintf('[SAVED] %s\n', summaryCsv);
fprintf('[SAVED] %s\n', plotDir);
end

function runFile = find_saved_run(resultRoot, refine, baseCandidates)
%Locate an index or object used by the computation.
baseCandidates = unique(baseCandidates(:).', 'stable');
for Csigma = baseCandidates
    candidate = fullfile(resultRoot, sprintf('base_%g', Csigma), ...
        sprintf('refine_%02d', refine), 'run.mat');
    if exist(candidate, 'file')
        runFile = candidate;
        return;
    end
end

matches = dir(fullfile(resultRoot, 'base_*', ...
    sprintf('refine_%02d', refine), 'run.mat'));
if ~isempty(matches)
    runFile = fullfile(matches(1).folder, matches(1).name);
    return;
end

error('No saved run.mat found for refine=%d under %s', refine, resultRoot);
end

function [Ksigma, M, P] = reconstruct_Ksigma(run, Csigma)
%Reconstruct cutoff-weighted penalty values.
M = 0.5 * (run.M + run.M');
P = run.P;
Sdg = run.S;

H = run.Mat + 0.5 * Sdg + 0.5 * Sdg' - run.meta.sigma * P;
H = 0.5 * (H + H');

sigma = Csigma * (1 / run.meta.h + run.meta.Nc);
Ksigma = H - 0.5 * Sdg - 0.5 * Sdg' + sigma * P;
Ksigma = 0.5 * (Ksigma + Ksigma');
end

function Pmat = build_condition_preconditioner_matrix_local(A_tau, penaltyMat, solveMode, epsD)
%Build condition preconditioner matrix.
n = size(A_tau, 1);
solveMode = lower(string(solveMode));
gamma = find(sum(abs(penaltyMat), 2) ~= 0);
eta = setdiff((1:n).', gamma);

switch solveMode
    case "none"
        Pmat = speye(n);

    case "purediag"
        d = max(abs(diag(A_tau)), epsD);
        Pmat = spdiags(d, 0, n, n);

    case "interfaceblock"
        dEta = max(abs(diag(A_tau(eta, eta))), epsD);
        Ag = A_tau(gamma, gamma);
        Ag = 0.5 * (Ag + Ag');
        [~, AgShifted] = make_block_positive_definite_local(Ag, epsD);

        Pmat = sparse(n, n);
        Pmat(eta, eta) = spdiags(dEta, 0, numel(eta), numel(eta));
        Pmat(gamma, gamma) = AgShifted;

    otherwise
        error('Unsupported solve_mode: %s', solveMode);
end

Pmat = 0.5 * (Pmat + Pmat');
end

function [delta, AgShifted] = make_block_positive_definite_local(Ag, epsD)
%Build block positive definite.
n = size(Ag, 1);
if n == 0
    delta = 0;
    AgShifted = sparse(0, 0);
    return;
end

scale = max(1, norm(Ag, 1));
floorEig = max(epsD, epsD * scale);
if n <= 3000
    lamMin = min(real(eig(full(Ag))));
else
    eigOpts = struct('issym', true, 'isreal', isreal(Ag), 'tol', 1e-8, 'maxit', 2000);
    lamMin = real(eigs(Ag, 1, 'smallestreal', eigOpts));
end

delta = max(0, floorEig - lamMin);
I = speye(n);
AgShifted = Ag + delta * I;

for attempt = 1:8
    [~, pflag] = chol(AgShifted);
    if pflag == 0
        return;
    end
    delta = max(10 * max(delta, floorEig), floorEig * 10^attempt);
    AgShifted = Ag + delta * I;
end
error('PITB-DG gamma block is not numerically positive definite after delta regularization.');
end

function kappa = shifted_generalized_abs_condition_local(A, B, maxN, label, opts)
%Compute generalized abs condition.
n = size(A, 1);
assert(n <= maxN, '%s: matrix size %d exceeds max_n = %g.', label, n, maxN);

fullMaxN = get_opt(opts, 'cond_full_max_n', 1400);
tol = get_opt(opts, 'cond_eigs_tol', 1e-8);
maxit = get_opt(opts, 'cond_eigs_maxit', 2000);

if n <= fullMaxN
    mu = eig(full(A), full(B), 'vector');
else
    eigOpts = struct('tol', tol, 'maxit', maxit);
    muMax = eigs(A, B, 1, 'largestabs', eigOpts);
    muMin = eigs(A, B, 1, 'smallestabs', eigOpts);
    mu = [muMax; muMin];
end

absMu = abs(mu(isfinite(mu)));
absMu = absMu(absMu > 0);
assert(~isempty(absMu), '%s: no nonzero generalized eigenvalues found for condition estimate.', label);
kappa = max(absMu) / min(absMu);
end

function figs = plot_condition_summary(summaryFile, opts)
%Plot condition summary.
S = load(summaryFile, 'summary');
summary = S.summary;
T = sortrows(summary.table, {'Csigma', 'refine'});
Tplot = filter_plot_table(T, opts);

style = default_style();
style = merge_struct(style, get_opt(opts, 'style', struct()));

outDir = get_opt(opts, 'outDir', fullfile(fileparts(summaryFile), 'figures'));
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

savePng = get_opt(opts, 'save_png', true);
savePdf = get_opt(opts, 'save_pdf', true);
pngDPI = get_opt(opts, 'pngDPI', 600);

CsigmaList = summary.CsigmaList(:).';
labels = arrayfun(@(x) sprintf('$C_\\sigma=%g$', x), CsigmaList, 'UniformOutput', false);

figs = struct();
fig = figure('Color', style.fig_color, ...
    'Units', style.fig_unit, ...
    'Position', style.fig_pos, ...
    'Renderer', style.fig_renderer);
ax = axes(fig);
set_axes_layout(ax, style);
plot_condition_comparison_panel(ax, Tplot, CsigmaList, labels, style);
export_figure(fig, fullfile(outDir, 'condition'), savePng, savePdf, pngDPI);
figs.none_vs_interfaceblock = fig;

fprintf('[PLOT ] saved C_sigma condition figure to %s\n', outDir);
end

function Tplot = filter_plot_table(T, opts)
%Compute plot table.
plotRefines = get_opt(opts, 'plot_refine_list', []);
if isempty(plotRefines)
    plotRefines = 2:min(6, max(T.refine));
end
Tplot = T(ismember(T.refine, plotRefines), :);
end

function plot_condition_panel(ax, Tplot, CsigmaList, labels, style, fieldName, yLabelStr)
%Plot condition panel.
hold(ax, 'on');
hLines = gobjects(numel(CsigmaList), 1);
Yall = [];
for k = 1:numel(CsigmaList)
    Tk = sortrows(Tplot(Tplot.Csigma == CsigmaList(k), :), 'h');
    y = Tk.(fieldName);
    y(y <= 0) = eps;
    hLines(k) = plot(ax, Tk.h, y, '-', ...
        'LineWidth', style.lw, ...
        'Color', style.lineColors(k,:), ...
        'Marker', style.markers{k}, ...
        'MarkerSize', style.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.lineColors(k,:), ...
        'DisplayName', labels{k});
    Yall = [Yall; y(:)]; %#ok<AGROW>
end

set_axes_style(ax, style);
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
xlab = xlabel(ax, "$h$", 'Interpreter', 'latex', 'FontSize', style.label_fs);
xlab.Units = 'normalized';
xlab.Position(2) = style.xlabel_norm_y;
ylabel(ax, yLabelStr, ...
    'Interpreter', 'latex', 'FontSize', style.label_fs);

apply_log_padding(ax, Tplot.h(:), Yall, style.padX, style.padYLow, style.padYHigh);
legend(ax, hLines, labels, ...
    'Location', style.leg_loc, ...
    'Interpreter', style.leg_interpreter, ...
    'FontSize', style.leg_fs, ...
    'Box', style.leg_box);
end

function plot_condition_comparison_panel(ax, Tplot, CsigmaList, labels, style)
%Plot condition comparison panel.
hold(ax, 'on');
hNoneLines = gobjects(numel(CsigmaList), 1);
hIBLines = gobjects(numel(CsigmaList), 1);
legendLabelsNone = cell(numel(CsigmaList), 1);
legendLabelsIB = cell(numel(CsigmaList), 1);
Yall = [];
for k = 1:numel(CsigmaList)
    Tk = sortrows(Tplot(Tplot.Csigma == CsigmaList(k), :), 'h');

    yNone = Tk.cond_none;
    yNone(yNone <= 0) = eps;
    hNoneLines(k) = plot(ax, Tk.h, yNone, '--', ...
        'LineWidth', style.lw, ...
        'Color', style.lineColors(k,:), ...
        'Marker', style.markers{k}, ...
        'MarkerSize', style.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.lineColors(k,:), ...
        'DisplayName', labels{k});
    legendLabelsNone{k} = labels{k};

    yIB = Tk.cond_interfaceblock;
    yIB(yIB <= 0) = eps;
    hIBLines(k) = plot(ax, Tk.h, yIB, '-', ...
        'LineWidth', style.lw, ...
        'Color', style.lineColors(k,:), ...
        'Marker', style.markers{k}, ...
        'MarkerSize', style.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.lineColors(k,:), ...
        'DisplayName', sprintf('$C_\\sigma=%g$, TB-DG', CsigmaList(k)));
    legendLabelsIB{k} = sprintf('$C_\\sigma=%g$, TB-DG', CsigmaList(k));

    Yall = [Yall; yNone(:); yIB(:)]; %#ok<AGROW>
end

set_axes_style(ax, style);
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
xlab = xlabel(ax, "$h$", 'Interpreter', 'latex', 'FontSize', style.label_fs);
xlab.Units = 'normalized';
xlab.Position(2) = style.xlabel_norm_y;
ylabel(ax, "Condition number", ...
    'Interpreter', 'latex', 'FontSize', style.label_fs);

apply_log_padding(ax, Tplot.h(:), Yall, style.padX, style.padYLow, style.padYHigh);
Ypos = Yall(isfinite(Yall) & Yall > 0);
yTop = max(ylim(ax));
if ~isempty(Ypos)
    yTop = max(yTop, max(Ypos) * style.comp_ylim_high_factor);
end
ylim(ax, [style.comp_ylim_low yTop]);
yExp = floor(log10(style.comp_ylim_low)):style.comp_ytick_step:ceil(log10(max(ylim(ax))));
if numel(yExp) < 5
    yExp = floor(log10(style.comp_ylim_low)):ceil(log10(max(ylim(ax))));
end
set(ax, 'YTick', 10 .^ yExp);

legendHandles = [hNoneLines; hIBLines];
legendLabels = [legendLabelsNone; legendLabelsIB];
leg = legend(ax, legendHandles, legendLabels, ...
    'Location', style.comp_leg_loc, ...
    'Interpreter', style.leg_interpreter, ...
    'FontSize', style.comp_leg_fs, ...
    'Box', style.leg_box, ...
    'NumColumns', style.comp_leg_num_columns);
leg.ItemTokenSize = style.comp_leg_token_size;
end

function set_axes_layout(ax, style)
%Apply axes layout settings.
ax.Units = 'normalized';
ax.Position = [ ...
    style.layout_left, ...
    style.layout_bottom, ...
    1 - style.layout_left - style.layout_right, ...
    1 - style.layout_bottom - style.layout_top];
end

function set_axes_style(ax, style)
%Apply axes style settings.
set(ax, ...
    'FontName', style.ax_fontname, ...
    'FontSize', style.ax_fontsize, ...
    'LineWidth', style.ax_lw, ...
    'TickDir', style.tickdir, ...
    'Box', style.box, ...
    'XMinorTick', style.xminor, ...
    'YMinorTick', style.yminor, ...
    'Layer', 'top');
grid(ax, 'off');
end

function apply_log_padding(ax, x, y, padX, padYLow, padYHigh)
%Apply log padding.
x = x(isfinite(x) & x > 0);
y = y(isfinite(y) & y > 0);
if isempty(x) || isempty(y)
    return;
end
set(ax, 'XLim', [min(x) / (1 + padX), max(x) * (1 + padX)]);
set(ax, 'YLim', [min(y) / (1 + padYLow), max(y) * (1 + padYHigh)]);
end

function export_figure(fig, baseName, savePng, savePdf, pngDPI)
%Export figure.
drawnow;
if savePdf
    exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
end
if savePng
    exportgraphics(fig, [baseName '.png'], 'Resolution', pngDPI);
end
end

function style = default_style()
%Return plotting style values.
style = struct();
style.fig_pos   = [1 1 4.8 3.0];
style.fig_color = 'w';
style.fig_unit  = 'inches';
style.fig_renderer = 'painters';
style.layout_left = 0.14;
style.layout_right = 0.04;
style.layout_bottom = 0.16;
style.layout_top = 0.08;

style.lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;
style.markers = {'o','s','^','d','x','+'};
style.lw = 1.8;
style.ms = 8;

style.ax_fontname = 'Times New Roman';
style.ax_fontsize = 11;
style.ax_lw       = 1.0;
style.tickdir     = 'out';
style.box         = 'on';
style.xminor      = 'off';
style.yminor      = 'off';
style.label_fs    = 13;
style.xlabel_norm_y = -0.075;
style.padX         = 0.4;
style.padYLow      = 6000.0;
style.padYHigh     = 2.5;
style.comp_ylim_low = 1.0;
style.comp_ylim_high_factor = 1e3;
style.comp_ytick_step = 2;

style.leg_loc         = 'southwest';
style.leg_box         = 'off';
style.leg_fs          = 10;
style.leg_interpreter = 'latex';
style.legend_line_color = [0 0 0];
style.comp_leg_loc = 'northeast';
style.comp_leg_fs = 10;
style.comp_leg_num_columns = 2;
style.comp_leg_token_size = [28 10];
end

function val = get_opt(s, fieldName, defaultVal)
%Return one option value.
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function out = merge_struct(base, override)
%Merge struct.
out = base;
if isempty(override)
    return;
end
names = fieldnames(override);
for i = 1:numel(names)
    out.(names{i}) = override.(names{i});
end
end

function tag = sanitize_tag(tag)
%Convert a value to a filename tag.
if isnumeric(tag)
    tag = join_numbers(tag);
end
if isstring(tag)
    tag = char(tag);
end
tag = strtrim(tag);
tag = regexprep(tag, '[^\w.-]+', '_');
tag = regexprep(tag, '^_+|_+$', '');
end

function tag = join_numbers(values)
%Format numeric values for text output.
parts = arrayfun(@(x) sprintf('%g', x), values(:).', 'UniformOutput', false);
tag = strjoin(parts, '_');
end
