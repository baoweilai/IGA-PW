function plot_cutoff_convergence()
% Plot cutoff-convergence data.
clc; close all;
format short g;

figDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(figDir);
addpath(fullfile(exampleDir, 'model', 'cutoff_convergence', 'nurbs'));

%% --------------------- Global interpreter ---------------------
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

% User settings
Example          = 'Example_2';
Nc_list          = [4 6 8 10 12];
p_list           = 2;
Refine_fixed     = 7;
eig_list         = [1 2 3];
n_eigenvalues    = 4;
dx               = 5e-3;

refNc            = 45;
refP             = 3;
refRefine        = 8;

refTag = sprintf('refine_%02d', Refine_fixed);

outSubDirName = 'cache_Nc';

cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

cfg.axes.fontSize   = 11;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.labelSize  = 13;
cfg.axes.box        = 'on';

useCustomColors = true;
lineColors255 = [ ...
    223 122  94;   % lambda_1
    60  64  91;   % lambda_2
    130 178 154;   % lambda_3
    242 204 142;   % lambda_4
    ];
lineStyles = {'-','-','-','-'};
markers    = {'o','s','^','d','x','+'};
lw = 1.8;
ms = 6;

cfg.legend.fontSize = 11;
legendLocation = 'northeast';
legendBox      = 'off';

xTickMode = 'Nc';        % 'Nc' | 'auto'
yTickMode = 'manual';    % 'auto' | 'manual'
yTickExp_p1 = [-8 -6 -4 -2 0];
yTickExp_p2 = [-7 -6 -5 -4 -3 -2 -1 0];

padX_left  = 0.06;
padX_right = 0.06;
padY_low   = 1;
padY_high  = 1;

% Paths
resultRoot = fullfile(exampleDir, 'data', 'result', Example);
if ~exist(resultRoot,'dir')
    error('Cannot find resultRoot: %s', resultRoot);
end

refCsvFile = fullfile(resultRoot, sprintf('Nc_%d', refNc), sprintf('p_%d', refP), 'summary.csv');
assert(exist(refCsvFile, 'file') == 2, 'Missing latest reference summary: %s', refCsvFile);
Tref = readtable(refCsvFile);
Tref = Tref(Tref.refine == refRefine, :);
assert(height(Tref) == 1, 'Expected exactly one latest reference row for refine=%d in %s.', refRefine, refCsvFile);
lambda_ref_latest = [Tref.lambda1, Tref.lambda2, Tref.lambda3, Tref.lambda4];
fprintf('[REF ] latest reference Nc=%d p=%d refine=%d lambda=%s\n', ...
    refNc, refP, refRefine, mat2str(lambda_ref_latest, 12));
refRun = fullfile(resultRoot, sprintf('Nc_%d', refNc), sprintf('p_%d', refP), ...
    sprintf('refine_%02d', refRefine), 'run.mat');
assert(exist(refRun, 'file') == 2, 'Missing reference run: %s', refRun);

outRoot = fullfile(resultRoot, outSubDirName, refTag);
if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

Nc_list = Nc_list(:).';

% Main loop
for pp = p_list
    fprintf('\n==================== [P=%d] fixed refine=%d, sweep Nc ====================\n', ...
        pp, Refine_fixed);

    [Nc_ok, lam_ok, dof_ok] = ...
        read_runs_over_Nc(resultRoot, Nc_list, pp, refTag, n_eigenvalues);

    lambda_ref = lambda_ref_latest;

    err = abs(lam_ok(:, eig_list) - lambda_ref(1, eig_list));
    labels = arrayfun(@(i) sprintf('$i=%d$', i), eig_list, 'UniformOutput', false);
    fig = plot_semilogy_data( ...
        Nc_ok, err, labels, '$|\lambda_i-\lambda_i^{\mathrm{DG}}|$', ...
        useCustomColors, lineColors255, lineStyles, markers, lw, ms, ...
        legendLocation, legendBox, cfg.legend.fontSize, ...
        xTickMode, yTickMode, yTickExp_p1, yTickExp_p2, ...
        padX_left, padX_right, padY_low, padY_high, ...
        cfg);

    pOutRoot = fullfile(outRoot, sprintf('p_%d', pp));
    if ~exist(pOutRoot, 'dir'), mkdir(pOutRoot); end
    save_plot(fig, pOutRoot, 'cutoff');

    files = arrayfun(@(K) fullfile(resultRoot, sprintf('Nc_%d', K), ...
        sprintf('p_%d', pp), refTag, 'run.mat'), Nc_ok, 'UniformOutput', false);
    [eigErr, sigma] = cutoff_eigerr(files, refRun, eig_list, dx);
    figEig = plot_semilogy_data( ...
        Nc_ok, eigErr, labels, '$\|u_i-u_i^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
        useCustomColors, lineColors255, lineStyles, markers, lw, ms, ...
        legendLocation, legendBox, cfg.legend.fontSize, ...
        xTickMode, yTickMode, yTickExp_p1, yTickExp_p2, ...
        padX_left, padX_right, padY_low, padY_high, ...
        cfg);
    save_plot(figEig, pOutRoot, 'cutoff_eigfun');

    Teig = table(Nc_ok(:), sigma, 'VariableNames', {'Nc','sigma'});
    for j = 1:numel(eig_list)
        Teig.(sprintf('u%d_DG', eig_list(j))) = eigErr(:, j);
    end
    eigCsv = fullfile(pOutRoot, 'eigfun_DG.csv');
    writetable(Teig, eigCsv);
    fprintf('[SAVED] %s\n', eigCsv);

    % Save the cutoff summary for the manuscript plot.
    Tsum = table(Nc_ok(:), dof_ok(:), 'VariableNames', {'Nc','dof'});
    for k = 1:n_eigenvalues
        Tsum.(sprintf('lambda%d',k)) = lam_ok(:,k);
    end
    out_csv = fullfile(pOutRoot, 'summary.csv');
    writetable(Tsum, out_csv);
    fprintf('[SAVED] %s\n', out_csv);

    % orders
    assert(all(err(:) > 0), 'Cutoff-convergence errors must be positive.');
    rate = local_exp_rates_Nc(Nc_ok(:), err);

    varNames = {'p','eig','Nc_from','Nc_to','err_to','rate_perNc'};
    if numel(Nc_ok) < 2
        OrderNcTable = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
    else
        Rows = cell((numel(Nc_ok)-1)*numel(eig_list), numel(varNames));
        rc = 0;
        for j = 1:numel(eig_list)
            eigIdx = eig_list(j);
            for i = 1:(numel(Nc_ok)-1)
                rc = rc + 1;
                Rows(rc,:) = {pp, eigIdx, Nc_ok(i), Nc_ok(i+1), ...
                    round_sig(err(i+1,j), 10), round_sig(rate(i,j), 10)};
            end
        end
        OrderNcTable = cell2table(Rows, 'VariableNames', varNames);
    end
    out_order_csv = fullfile(pOutRoot, 'orders.csv');
    writetable(OrderNcTable, out_order_csv);
    fprintf('[SAVED] %s\n', out_order_csv);
end

fprintf('\n[DONE] All outputs in: %s\n', outRoot);
end

function save_plot(fig, outDir, name)
% Save one figure as a PDF.
file = fullfile(outDir, [name '.pdf']);
exportgraphics(fig, file, 'ContentType', 'vector');
fprintf('[SAVED] %s\n', file);
end

function [Nc_ok, lam_ok, dof_ok] = ...
read_runs_over_Nc(resultRoot, Nc_list, pdeg, refTag, n_eigs)

% Load cutoff results across the requested plane-wave cutoffs.
nNc   = numel(Nc_list);
lamNc = zeros(nNc, n_eigs);
dofNc = zeros(nNc, 1);

for i = 1:nNc
    Nc = Nc_list(i);
    runMat = fullfile(resultRoot, sprintf('Nc_%d',Nc), sprintf('p_%d',pdeg), refTag, 'run.mat');
    assert(exist(runMat, 'file') == 2, 'Missing run file: %s', runMat);

    S = load(runMat, 'run');
    assert(isfield(S, 'run'), 'Missing run structure in %s.', runMat);
    run = S.run;
    assert(isfield(run, 'lambda'), 'Missing eigenvalues in %s.', runMat);
    assert(isfield(run, 'n_dofs_total'), 'Missing total DOFs in %s.', runMat);

    lam = run.lambda(:).';
    assert(numel(lam) >= n_eigs, 'The run file has too few eigenvalues: %s', runMat);
    lamNc(i,:) = lam(1:n_eigs);
    dofNc(i) = run.n_dofs_total;
    fprintf('[LOAD] Nc=%d p=%d  lambda(1:%d)=%s\n', Nc, pdeg, n_eigs, mat2str(lamNc(i,:), 12));
end

Nc_ok     = Nc_list(:);
lam_ok    = lamNc;
dof_ok    = dofNc;

[Nc_ok, idx] = sort(Nc_ok, 'ascend');
lam_ok    = lam_ok(idx,:);
dof_ok    = dof_ok(idx);
end

function fig = plot_semilogy_data( ...
Nc_ok, err, legLabels, yLabel, ...
    useCustomColors, lineColors255, lineStyles, markers, lw, ms, ...
    legendLocation, legendBox, legendFontSize, ...
    xTickMode, yTickMode, yTickExp_p1, yTickExp_p2, ...
    padX_left, padX_right, padY_low, padY_high, ...
    cfg)

% Plot and style a semilog convergence figure.
% Validate the error data and create the plotting axes.
assert(all(err(:) > 0), 'Cutoff-convergence errors must be positive.');

fig = figure( ...
    'Color',    cfg.fig.bgColor, ...
    'Units',    'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);

ax = axes(fig);
hold(ax,'on');
box(ax,'on');

ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1-cfg.layout.left-cfg.layout.right, ...
    1-cfg.layout.bottom-cfg.layout.top];

% Select line styles and draw every eigenvalue curve.
if useCustomColors
    colors = double(lineColors255)/255;
else
    colors = lines(max(8, size(err, 2)));
end
nC  = size(colors,1);
nLS = numel(lineStyles);
nMK = numel(markers);

hEig = gobjects(1, size(err, 2));

for j = 1:size(err, 2)
    col = colors(mod(j-1,nC)+1,:);
    ls  = lineStyles{mod(j-1,nLS)+1};
    mk  = markers{mod(j-1,nMK)+1};

    hEig(j) = semilogy(ax, Nc_ok, err(:,j), ls, ...
        'LineWidth', lw, ...
        'Color', col, ...
        'Marker', mk, ...
        'MarkerSize', ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', col);
end

% Apply logarithmic axes, labels, and legend styling.
set(ax, ...
    'XScale',     'linear', ...
    'YScale',     'log', ...
    'FontSize',   cfg.axes.fontSize, ...
    'LineWidth',  cfg.axes.lineWidth, ...
    'TickDir',    cfg.axes.tickDir, ...
    'Box',        cfg.axes.box, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick);

ax.TickLabelInterpreter = 'latex';

grid(ax,'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

ax.XRuler.MinorTick = 'off';
ax.YRuler.MinorTick = 'off';

xlabel(ax, '$K$', ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);

ylabel(ax, yLabel, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);

lgd = legend(ax, hEig, legLabels, ...
    'Location', legendLocation, ...
    'Interpreter', 'latex', ...
    'FontSize', legendFontSize);
lgd.Box = legendBox;

% Derive axis limits and ticks from the plotted data.
xMin = min(Nc_ok);
xMax = max(Nc_ok);
xRng = xMax - xMin;
if xRng == 0
    xRng = 1;
end

axisSpec.XLim = [xMin - padX_left*xRng, xMax + padX_right*xRng];
switch lower(char(xTickMode))
    case 'nc'
        axisSpec.XTick = Nc_ok(:).';
    case 'auto'
        axisSpec.XTick = ax.XTick;
    otherwise
        error('Unknown xTickMode=%s (use ''Nc'' or ''auto'')', xTickMode);
end

yAll = err(:);
yMin = min(yAll);
yMax = max(yAll);

switch lower(char(yTickMode))
    case 'auto'
        emin = floor(log10(yMin));
        emax = ceil(log10(yMax));
        if emin == emax
            exps = emin + (-2:2);
        else
            exps = round(linspace(emin, emax, 5));
            exps = unique(exps,'stable');
        end
        exps = sort(exps(:).');
    case 'manual'
        if pdeg == 1
            exps = sort(yTickExp_p1(:).');
        else
            exps = sort(yTickExp_p2(:).');
        end
    otherwise
        error('Unknown yTickMode=%s (use ''auto'' or ''manual'')', yTickMode);
end

axisSpec.YTick = 10.^exps;
axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), exps, 'UniformOutput', false);
axisSpec.YLim = [yMin/(1+padY_low), yMax*(1+padY_high)];

ax.XLim = axisSpec.XLim;
ax.XTick = axisSpec.XTick;
ax.YTick = axisSpec.YTick;
ax.YTickLabel = axisSpec.YTickLabel;
ax.YLim = axisSpec.YLim;
end

function rate = local_exp_rates_Nc(Nc, err)
% Estimate exponential convergence rates between successive cutoffs.
Nc = Nc(:);
n  = numel(Nc);
ne = size(err,2);
rate = zeros(n-1, ne);
for i = 1:n-1
    dNc = Nc(i+1) - Nc(i);
    assert(dNc ~= 0, 'Cutoff values must be distinct.');
    rate(i,:) = -(log10(err(i+1,:)) - log10(err(i,:))) / dNc;
end
end

function y = round_sig(x, nSig)
% Round a value to significant digits.
y = x;
mask = isfinite(x) & (x ~= 0);
ax = abs(x(mask));
p  = floor(log10(ax));
scale = 10.^(nSig-1-p);
y(mask) = round(x(mask).*scale)./scale;
end
