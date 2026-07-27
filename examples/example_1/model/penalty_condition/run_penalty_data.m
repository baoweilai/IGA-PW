function summary = run_penalty_data(optsIn)
% Generate retained penalty and condition-number results.

assert(exist('optsIn', 'var') == 1, 'run_penalty_data requires an options structure.');

clc;
close all;
format short g;

% Set plotting, workflow, and sweep parameters.
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

activate_example_workflow('penalty_condition', ...
    {'nurbs', 'iga', 'assembly', 'operators', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
rootDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(rootDir);

Example = 'Example_1';
Nc = get_opt(optsIn, 'Nc', 25);
t = get_opt(optsIn, 't', 1);
pdeg = 1 + t;
refineList = get_opt(optsIn, 'refineList', 1:7);
CsigmaList = get_opt(optsIn, 'CsigmaList', [10 100 1000]);
baseCsigma = get_opt(optsIn, 'baseCsigma', 20);
nEigenvalues = 2;

lambdaRef = [4.96999274451623, 6.374026300804];
tauShift = get_opt(optsIn, 'tau_shift', 0.9 * lambdaRef(1));
epsDiag = get_opt(optsIn, 'eps_diag', 1e-12);
ifaceReg = get_opt(optsIn, 'iface_reg', 1e-12);

primmeTol = get_opt(optsIn, 'primme_tol', 1e-8);
primmeMaxit = get_opt(optsIn, 'primme_maxit', 1e7);
primmeMethod = get_opt(optsIn, 'primme_method', 'DEFAULT_MIN_MATVECS');
primmeTarget = get_opt(optsIn, 'primme_target', 'SA');
primmeReportLevel = get_opt(optsIn, 'primme_reportLevel', 0);
seed = get_opt(optsIn, 'seed', 20260404);

sweepTag = sanitize_tag(get_opt(optsIn, 'sweepTag', ''));
if isempty(sweepTag)
    fileSuffix = '';
    plotName = 'figures';
else
    fileSuffix = ['_' sweepTag];
    plotName = ['figures_' sweepTag];
end

resultRoot = fullfile(rootDir, 'result', 'penalty');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

summaryMat = fullfile(resultRoot, ['summary' fileSuffix '.mat']);
summaryCsv = fullfile(resultRoot, ['summary' fileSuffix '.csv']);
plotDir = fullfile(resultRoot, plotName);

cacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRoot, 'dir')
    mkdir(cacheRoot);
end

rows = [];

% Evaluate every refinement and penalty coefficient.
for refine = refineList(:).'
    fprintf('\n============================================================\n');
    fprintf('[BASE] refine=%d, Nc=%d, p=%d, base C_sigma=%g\n', refine, Nc, pdeg, baseCsigma);

    baseRun = load_or_build_base_run(refine, t, Nc, nEigenvalues, ...
        baseCsigma, lambdaRef(1), tauShift, epsDiag, ifaceReg, ...
        primmeTol, primmeMaxit, primmeMethod, primmeTarget, ...
        primmeReportLevel, cacheRoot, resultRoot, Example);

    for Csigma = CsigmaList(:).'
        if abs(Csigma - baseCsigma) <= eps(max(1, abs(baseCsigma))) ...
                && isfield(baseRun.run.result, 'None') ...
                && isfield(baseRun.run.result, 'InterfaceBlock')
            row = pack_base_row(baseRun.run, Csigma, lambdaRef);
        else
            row = solve_saved_csigma(baseRun.run, Csigma, lambdaRef, ...
                tauShift, epsDiag, ifaceReg, primmeTol, primmeMaxit, ...
                primmeMethod, primmeTarget, primmeReportLevel, seed);
        end
        rows = [rows; row]; %#ok<AGROW>
    end
end

% Save the sweep data and generate the penalty figure.
T = struct2table(rows);
T = sortrows(T, {'Csigma', 'refine'});

summary = struct();
summary.Example = Example;
summary.Nc = Nc;
summary.pdeg = pdeg;
summary.t = t;
summary.refineList = refineList;
summary.CsigmaList = CsigmaList;
summary.baseCsigma = baseCsigma;
summary.sweepTag = sweepTag;
summary.lambdaRef = lambdaRef;
summary.tau_shift = tauShift;
summary.table = T;

save(summaryMat, 'summary');
writetable(T, summaryCsv);

plot_penalty_sweep(summaryMat, struct('outDir', plotDir));

fprintf('\n[SAVED] %s\n', summaryMat);
fprintf('[SAVED] %s\n', summaryCsv);
fprintf('[SAVED] %s\n', plotDir);
end

function baseRun = load_or_build_base_run(refine, t, Nc, nEigenvalues, ...
    Csigma, lambdaRef1, tauShift, epsDiag, ifaceReg, primmeTol, ...
    primmeMaxit, primmeMethod, primmeTarget, primmeReportLevel, ...
    cacheRoot, resultRoot, Example)
% Load or build the base run.

baseDir = fullfile(resultRoot, sprintf('base_%g', Csigma), ...
    sprintf('refine_%02d', refine));
runFile = fullfile(baseDir, 'run.mat');

if exist(runFile, 'file')
    S = load(runFile, 'run');
    if has_required_saved_matrices(S.run, Csigma)
        fprintf('[LOAD] %s\n', runFile);
        baseRun = S;
        return;
    end
    fprintf('[REBUILD] saved run lacks full None/TB-DG data or matching C_sigma; rebuilding.\n');
end

if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

runOpts = struct();
runOpts.Example = Example;
runOpts.beta = Csigma;
runOpts.n_gp = 10;
runOpts.L = 4;
runOpts.a = 0.2;
runOpts.lambda_ref = lambdaRef1;
runOpts.tau_shift = tauShift;
runOpts.eps_diag = epsDiag;
runOpts.iface_reg = ifaceReg;
runOpts.primme_tol = primmeTol;
runOpts.primme_maxit = primmeMaxit;
runOpts.primme_method = primmeMethod;
runOpts.primme_target = primmeTarget;
runOpts.primme_reportLevel = primmeReportLevel;
runOpts.enabled_variants = {'None', 'InterfaceBlock'};
runOpts.profile_mode = true;
runOpts.use_pw_cache = true;
runOpts.cacheRoot = cacheRoot;
runOpts.inner_cheb_n = 200;
runOpts.pw_fft_grid_n = 256;
runOpts.save_eigenvectors = true;
runOpts.save_matrices = true;
runOpts.save_mat = true;
runOpts.outDir = baseDir;

[~, ~, ~] = solve_iga_pw_dg(refine, t, Nc, nEigenvalues, runOpts);

S = load(runFile, 'run');
baseRun = S;
end

function tf = has_required_saved_matrices(run, Csigma)
% Check whether a saved run contains the matrices required for this penalty.
tf = isfield(run, 'M') && isfield(run, 'Mat') && isfield(run, 'P') ...
    && isfield(run, 'S') && isfield(run, 'meta') ...
    && isfield(run, 'result') && isfield(run.result, 'None') ...
    && isfield(run.result, 'InterfaceBlock') ...
    && isfield(run.meta, 'beta') && abs(run.meta.beta - Csigma) <= eps(max(1, abs(Csigma)));
end

function row = pack_base_row(run, Csigma, lambdaRef)
% Pack the common fields for one penalty-result row.
resNone = run.result.None;
resIB = run.result.InterfaceBlock;
assert(numel(resNone.lambda) >= 2, 'None result must contain two eigenvalues.');
assert(numel(resIB.lambda) >= 2, 'InterfaceBlock result must contain two eigenvalues.');
lambdaNone = resNone.lambda(1:2);
lambdaIB = resIB.lambda(1:2);

row = common_row_fields(run, Csigma, lambdaIB, lambdaRef);
row.lambda1_none = lambdaNone(1);
row.lambda2_none = lambdaNone(2);
row.err1_none = abs(lambdaNone(1) - lambdaRef(1));
row.err2_none = abs(lambdaNone(2) - lambdaRef(2));
row.total_s_none = resNone.time_total;
row.primme_s_none = resNone.time_primme;
row.lambda1_ib = lambdaIB(1);
row.lambda2_ib = lambdaIB(2);
row.err1_ib = abs(lambdaIB(1) - lambdaRef(1));
row.err2_ib = abs(lambdaIB(2) - lambdaRef(2));
row.total_s_ib = resIB.time_total;
row.build_prec_s_ib = resIB.time_build_prec;
row.primme_s_ib = resIB.time_primme;

fprintf('[CSIG] C_sigma=%g refine=%d none_err1=%.3e TB-DG_err1=%.3e\n', ...
    Csigma, row.refine, row.err1_none, row.err1_ib);
end

function row = solve_saved_csigma(run, Csigma, lambdaRef, tauShift, epsDiag, ...
    ifaceReg, primmeTol, primmeMaxit, primmeMethod, primmeTarget, ...
    primmeReportLevel, seed)
% Solve one saved penalty case.

% Form the symmetric operator for the selected penalty value.
M = 0.5 * (run.M + run.M');
P = run.P;
Sdg = run.S;
n = run.n_dofs_total;
nI = run.n_dofs_nurbs;
idxI = 1:nI;
idxP = (nI + 1):n;

H = run.Mat + 0.5 * Sdg + 0.5 * Sdg' - run.meta.sigma * P;
H = 0.5 * (H + H');

sigma = Csigma * (1 / run.meta.h + run.meta.Nc);
Mat = H - 0.5 * Sdg - 0.5 * Sdg' + sigma * P;
Mat = 0.5 * (Mat + Mat');
Atau = Mat - tauShift * M;

Hpw = H(idxP, idxP);
Mpw = M(idxP, idxP);
diagP = full(diag(P));

% Build diagonal and interface-block preconditioner components.
tBuild = tic;
dII = abs(diag(Atau(idxI, idxI)));
dII(dII < epsDiag) = 1;
dPP = abs(diag(Hpw) + sigma * diagP(idxP) - tauShift * diag(Mpw));
dPP(dPP < epsDiag) = 1;

bd = struct();
bd.dinv = [1 ./ dII; 1 ./ dPP];

gamma = find(sum(abs(P), 2) ~= 0);
Ag = 0.5 * (Atau(gamma, gamma) + Atau(gamma, gamma)');
delta = ifaceReg * max(1, norm(Ag, 1));
solveGamma = decomposition(Ag + delta * speye(size(Ag)), 'chol');
timeBuild = toc(tBuild);

precfun = @(x) apply_interface_block_prec(x, bd, gamma, solveGamma);

% Set a reproducible starting vector and PRIMME options.
rng(seed, 'twister');
v0 = randn(n, 1) + 1i * randn(n, 1);
v0 = v0 / max(norm(v0), eps);

ops = struct();
ops.tol = primmeTol;
ops.maxit = primmeMaxit;
ops.reportLevel = primmeReportLevel;
ops.v0 = v0;

% Solve the unpreconditioned and TB-DG eigenproblems.
tSolveNone = tic;
[~, Dnone] = primme_eigs(Mat, M, 2, primmeTarget, ops, primmeMethod);
timeSolveNone = toc(tSolveNone);

lambdaNone = sort(real(diag(Dnone)), 'ascend');
lambdaNone = lambdaNone(:).';
assert(numel(lambdaNone) >= 2, 'No-preconditioner PRIMME result must contain two eigenvalues.');

tSolve = tic;
[~, D] = primme_eigs(Mat, M, 2, primmeTarget, ops, primmeMethod, precfun);
timeSolve = toc(tSolve);

lambdaIB = sort(real(diag(D)), 'ascend');
lambdaIB = lambdaIB(:).';
assert(numel(lambdaIB) >= 2, 'TB-DG PRIMME result must contain two eigenvalues.');

% Collect eigenvalue errors and timings for the summary row.
row = common_row_fields(run, Csigma, lambdaIB, lambdaRef);
row.lambda1_none = lambdaNone(1);
row.lambda2_none = lambdaNone(2);
row.err1_none = abs(lambdaNone(1) - lambdaRef(1));
row.err2_none = abs(lambdaNone(2) - lambdaRef(2));
row.total_s_none = timeSolveNone;
row.primme_s_none = timeSolveNone;
row.lambda1_ib = lambdaIB(1);
row.lambda2_ib = lambdaIB(2);
row.err1_ib = abs(lambdaIB(1) - lambdaRef(1));
row.err2_ib = abs(lambdaIB(2) - lambdaRef(2));
row.total_s_ib = timeBuild + timeSolve;
row.build_prec_s_ib = timeBuild;
row.primme_s_ib = timeSolve;

fprintf('[CSIG] C_sigma=%g refine=%d none_err1=%.3e TB-DG_err1=%.3e solve=%.6fs\n', ...
    Csigma, row.refine, row.err1_none, row.err1_ib, timeSolveNone + timeBuild + timeSolve);
end

function y = apply_interface_block_prec(x, bd, gamma, solveGamma)
% Apply the diagonal and interface-block preconditioner.
y = bsxfun(@times, x, bd.dinv);
y(gamma, :) = solveGamma \ x(gamma, :);
end

function row = common_row_fields(run, Csigma, lambda, lambdaRef)
% Pack shared fields for one penalty-result row.
row = struct();
row.Csigma = Csigma;
row.refine = run.meta.Refinement;
row.h = run.meta.h;
row.dof = run.n_dofs_total;
row.sigma = Csigma * (1 / run.meta.h + run.meta.Nc);
row.lambda1 = lambda(1);
row.lambda2 = lambda(2);
row.err1 = abs(lambda(1) - lambdaRef(1));
row.err2 = abs(lambda(2) - lambdaRef(2));
end

function figs = plot_penalty_sweep(summaryFile, opts)
% Plot csigma sweep.
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

CsigmaList = summary.CsigmaList(:).';
labels = arrayfun(@(x) sprintf('$C_\\sigma=%g$, TB-DG', x), CsigmaList, 'UniformOutput', false);

figs = struct();
figs.err = figure('Color', style.fig_color, ...
    'Units', style.fig_unit, ...
    'Position', style.fig_pos, ...
    'Renderer', style.fig_renderer);
ax = axes(figs.err);
set_axes_layout(ax, style);
plot_lambda1_panel(ax, Tplot, CsigmaList, labels, style);

export_figure(figs.err, fullfile(outDir, 'error_zoom'));
export_figure(figs.err, fullfile(outDir, 'error'));

fprintf('[PLOT ] saved C_sigma figures to %s\n', outDir);
end

function Tplot = filter_plot_table(T, opts)
% Filter summary rows to the requested refinements.
plotRefines = get_opt(opts, 'plot_refine_list', []);
if isempty(plotRefines)
    plotRefines = 2:min(6, max(T.refine));
end
Tplot = T(ismember(T.refine, plotRefines), :);
end

function plot_lambda1_panel(ax, Tplot, CsigmaList, labels, style)
% Plot lambda1 panel.
% Draw the available penalty-error curve families.
hold(ax, 'on');
Yall = [];

hasNone = ismember('err1_none', Tplot.Properties.VariableNames) ...
    && any(isfinite(Tplot.err1_none));
if hasNone
    assert(ismember('err1_ib', Tplot.Properties.VariableNames), ...
        'The comparison table is missing err1_ib.');
    hNone = gobjects(numel(CsigmaList), 1);
    hIB = gobjects(numel(CsigmaList), 1);
    labelsNone = arrayfun(@(x) sprintf('$C_\\sigma=%g$', x), ...
        CsigmaList, 'UniformOutput', false);
    labelsIB = arrayfun(@(x) sprintf('$C_\\sigma=%g$, TB-DG', x), ...
        CsigmaList, 'UniformOutput', false);

    for k = 1:numel(CsigmaList)
        Tk = sortrows(Tplot(Tplot.Csigma == CsigmaList(k), :), 'h');
        yIB = Tk.err1_ib;
        yIB(yIB <= 0) = eps;
        hIB(k) = plot(ax, Tk.h, yIB, '-', ...
            'LineWidth', style.lw, ...
            'Color', style.lineColors(k,:), ...
            'Marker', style.markers{k}, ...
            'MarkerSize', style.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', style.lineColors(k,:), ...
            'DisplayName', labelsIB{k});

        yNone = Tk.err1_none;
        yNone(yNone <= 0) = eps;
        hNone(k) = plot(ax, Tk.h, yNone, '--', ...
            'LineWidth', style.lw, ...
            'Color', style.lineColors(k,:), ...
            'Marker', style.markers{k}, ...
            'MarkerSize', style.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', style.lineColors(k,:), ...
            'DisplayName', labelsNone{k});

        Yall = [Yall; yNone(:); yIB(:)]; %#ok<AGROW>
    end

    legendHandles = [hNone; hIB];
    legendLabels = [labelsNone(:); labelsIB(:)];
else
    hLines = gobjects(numel(CsigmaList), 1);
    for k = 1:numel(CsigmaList)
        Tk = sortrows(Tplot(Tplot.Csigma == CsigmaList(k), :), 'h');
        y = Tk.err1;
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
    legendHandles = hLines;
    legendLabels = labels;
end

% Apply logarithmic axes and data-dependent padding.
set_axes_style(ax, style);
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
xlab = xlabel(ax, "$h$", 'Interpreter', 'latex', 'FontSize', style.label_fs);
xlab.Units = 'normalized';
xlab.Position(2) = style.xlabel_norm_y;
ylabel(ax, "$|\lambda_1-\lambda_{1}^{\mathrm{DG}}|$", ...
    'Interpreter', 'latex', 'FontSize', style.label_fs);

% Add the combined legend for the displayed curves.
apply_log_padding(ax, Tplot.h(:), Yall, style.padX, style.padYLow, style.padYHigh);
leg = legend(ax, legendHandles, legendLabels, ...
    'Location', style.leg_loc, ...
    'Interpreter', style.leg_interpreter, ...
    'FontSize', style.leg_fs, ...
    'Box', style.leg_box);
leg.ItemTokenSize = style.leg_token_size;
end

function set_axes_style(ax, style)
% Apply axes style settings.
set(ax, ...
    'FontName', style.ax_fontname, ...
    'FontSize', style.ax_fontsize, ...
    'LineWidth', style.ax_lw, ...
    'TickDir', style.tickdir, ...
    'Box', style.box, ...
    'XMinorTick', style.xminor, ...
    'YMinorTick', style.yminor);
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
grid(ax, 'off');
end

function set_axes_layout(ax, style)
% Apply axes layout settings.
ax.Units = 'normalized';
ax.Position = [ ...
    style.layout_left, ...
    style.layout_bottom, ...
    1 - style.layout_left - style.layout_right, ...
    1 - style.layout_bottom - style.layout_top];
end

function apply_log_padding(ax, x, y, padX, padYLow, padYHigh)
% Apply log padding.
x = x(isfinite(x) & x > 0);
y = y(isfinite(y) & y > 0);
if isempty(x) || isempty(y)
    return;
end
set(ax, 'XLim', [min(x) / (1 + padX), max(x) * (1 + padX)]);
set(ax, 'YLim', [min(y) / (1 + padYLow), max(y) * (1 + padYHigh)]);
end

function export_figure(fig, baseName)
% Export the figure as a PDF.
drawnow;
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
end

function style = default_style()
% Return plotting style values.
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
style.padYLow      = 2.0;
style.padYHigh     = 1.0;
style.leg_loc         = 'southeast';
style.leg_box         = 'off';
style.leg_fs          = 10;
style.leg_interpreter = 'latex';
style.leg_token_size = [34 10];
end

function val = get_opt(s, fieldName, defaultVal)
% Read one option with a default value.
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function out = merge_struct(base, override)
% Merge override fields into the base structure.
out = base;
if isempty(override)
    return;
end
f = fieldnames(override);
for k = 1:numel(f)
    out.(f{k}) = override.(f{k});
end
end

function tag = sanitize_tag(tag)
% Convert a value to a filename tag.
if isstring(tag)
    tag = char(tag);
end
if isempty(tag)
    return;
end
tag = regexprep(char(tag), '[^A-Za-z0-9_=-]', '_');
end
