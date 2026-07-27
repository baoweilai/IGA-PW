function plot_cutoff_convergence()
% Plot cutoff-convergence data.

clc; close all;
format short g;

figDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(figDir);
addpath(fullfile(exampleDir, 'model', 'cutoff_convergence', 'nurbs'));

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

% 1. User parameters

% -------------------- problem setup --------------------
cfg.Example = 'Example_1';
cfg.Nc_list = [4 6 8 10 12];
cfg.p_list  = [1 2];

cfg.Refine_fixed  = 7;
cfg.eig_list      = [1 2 3 4];
cfg.n_eigenvalues = 4;
cfg.ref_Nc        = 45;
cfg.ref_p         = 3;
cfg.ref_refine    = 8;
cfg.dx_in         = 5e-3;
cfg.dx_out        = 5e-3;

% -------------------- reference eigenvalues --------------------
cfg.lambda_ref = [ ...
    4.969971740613, ...
    6.374026300804, ...
    6.847114865135, ...
    6.847114865135 ...
    ];

% -------------------- output settings --------------------
cfg.refTag        = sprintf('refine_%02d', cfg.Refine_fixed);
cfg.outSubDirName = 'cache_Nc';

% -------------------- figure settings --------------------
cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

% -------------------- layout settings --------------------
cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

% -------------------- axis style --------------------
cfg.axes.fontSize   = 11;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.labelSize  = 13;
cfg.axes.box        = 'on';
cfg.axes.xScale     = 'linear';
cfg.axes.yScale     = 'log';

% -------------------- curve style --------------------
cfg.lineColors255 = [ ...
    223 122  94;
    60  64  91;
    130 178 154;
    242 204 142
    ];

cfg.lineStyles   = {'-','-','-','-'};
cfg.markers      = {'o','s','^','d','x','+'};
cfg.lw           = 1.8;
cfg.ms           = 8;

% -------------------- legend style --------------------
cfg.legend.location   = 'northeast';
cfg.legend.fontSize   = 11;
cfg.legend.box        = 'off';
cfg.legend.numColumns = 1;

% -------------------- y-axis manual control --------------------
cfg.yTickExp_p1 = [-1 -2 -3 -4 -5 -6 -7];
cfg.yTickExp_p2 = [-1 -2 -3 -4 -5 -6 -7];

cfg.manualYLim_p1 = [1e-7, 1e-1];
cfg.manualYLim_p2 = [1e-7, 1e-1];

% -------------------- plot padding --------------------
cfg.padX_left  = 0.10;
cfg.padX_right = 0.10;
cfg.padY_low   = 1;
cfg.padY_high  = 1;

% 2. Build paths

resultRoot = fullfile(exampleDir, 'data', 'result', cfg.Example);
if ~exist(resultRoot, 'dir')
    error('Cannot find resultRoot: %s', resultRoot);
end

outRoot = fullfile(resultRoot, cfg.outSubDirName, cfg.refTag);
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end

cfg.Nc_list = cfg.Nc_list(:).';

runFiles = cell(numel(cfg.p_list), numel(cfg.Nc_list));
for ip = 1:numel(cfg.p_list)
    for iNc = 1:numel(cfg.Nc_list)
        runFiles{ip, iNc} = fullfile(resultRoot, ...
            sprintf('Nc_%d', cfg.Nc_list(iNc)), ...
            sprintf('p_%d', cfg.p_list(ip)), cfg.refTag, 'run.mat');
        assert(exist(runFiles{ip, iNc}, 'file') == 2, ...
            'Missing run file: %s', runFiles{ip, iNc});
    end
end

refFile = fullfile(resultRoot, sprintf('Nc_%d', cfg.ref_Nc), ...
    sprintf('p_%d', cfg.ref_p), sprintf('refine_%02d', cfg.ref_refine), 'run.mat');
assert(exist(refFile, 'file') == 2, 'Missing reference run file: %s', refFile);

cacheFile = fullfile(outRoot, 'eigenfunction_dg_errors.mat');
[eigfunDG, sigmaDG] = load_or_compute_eigenfunction_errors( ...
    cacheFile, runFiles, refFile, cfg);

% 3. Main loop: fixed refine, sweep Nc

for ip = 1:numel(cfg.p_list)
    pp = cfg.p_list(ip);
    fprintf('\n==================== [P=%d] fixed refine=%d, sweep Nc ====================\n', ...
        pp, cfg.Refine_fixed);

    [Nc_ok, lam_ok, dof_ok] = ...
        read_runs_over_Nc(resultRoot, cfg.Nc_list, pp, cfg.refTag, cfg.n_eigenvalues);

    lambda_ref = cfg.lambda_ref;

    fig = plot_semilogy_err_vs_Nc(Nc_ok, lam_ok, lambda_ref, cfg.eig_list, pp, cfg);

    pOutRoot = fullfile(outRoot, sprintf('p_%d', pp));
    if ~exist(pOutRoot, 'dir'), mkdir(pOutRoot); end
    save_plot(fig, pOutRoot, 'cutoff');

    eigfunDG_p = reshape(eigfunDG(ip, :, :), numel(cfg.Nc_list), []);
    assert(isequal(Nc_ok(:), cfg.Nc_list(:)), ...
        'The eigenvalue and eigenfunction cutoff grids do not agree.');
    figEigfun = plot_semilogy_eigfun_vs_Nc( ...
        Nc_ok, eigfunDG_p, cfg.eig_list, pp, cfg);
    save_plot(figEigfun, pOutRoot, 'cutoff_eigenfunction');

    Teigfun = table(Nc_ok(:), sigmaDG(ip, :).', 'VariableNames', {'Nc','sigma'});
    for j = 1:numel(cfg.eig_list)
        Teigfun.(sprintf('u%d_DG', cfg.eig_list(j))) = eigfunDG_p(:, j);
    end
    outEigfunCsv = fullfile(pOutRoot, 'eigenfunction_DG.csv');
    writetable(Teigfun, outEigfunCsv);
    fprintf('[SAVED] %s\n', outEigfunCsv);

    Tsum = table(Nc_ok(:), dof_ok(:), 'VariableNames', {'Nc','dof'});

    for k = 1:cfg.n_eigenvalues
        Tsum.(sprintf('lambda%d', k)) = lam_ok(:,k);
    end

    out_csv = fullfile(pOutRoot, 'summary.csv');
    writetable(Tsum, out_csv);

    fprintf('[SAVED] %s\n', out_csv);

    err = abs(lam_ok(:, cfg.eig_list) - lambda_ref(1, cfg.eig_list));
    assert(all(err(:) > 0), 'Cutoff-convergence errors must be positive.');

    rate = local_exp_rates_Nc(Nc_ok(:), err);

    varNames = {'p','eig','Nc_from','Nc_to','err_to','rate_perNc'};

    if numel(Nc_ok) < 2
        OrderNcTable = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
    else
        Rows = cell((numel(Nc_ok)-1) * numel(cfg.eig_list), numel(varNames));
        rc = 0;

        for j = 1:numel(cfg.eig_list)
            eigIdx = cfg.eig_list(j);
            for i = 1:(numel(Nc_ok)-1)
                rc = rc + 1;
                Rows(rc,:) = { ...
                    pp, ...
                    eigIdx, ...
                    Nc_ok(i), ...
                    Nc_ok(i+1), ...
                    round_sig(err(i+1,j), 10), ...
                    round_sig(rate(i,j), 10) ...
                    };
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

function [errDG, sigmaDG] = load_or_compute_eigenfunction_errors( ...
    cacheFile, runFiles, refFile, cfg)
% The cached aligned errors correspond to the current runs and integration code.
helperFile = fullfile(fileparts(mfilename('fullpath')), ...
    'cutoff_eigenfunction_errors.m');
sourceFiles = [runFiles(:); {refFile}; {helperFile}];
newestSource = 0;
for i = 1:numel(sourceFiles)
    info = dir(sourceFiles{i});
    assert(~isempty(info), 'Missing eigenfunction-error source: %s', sourceFiles{i});
    newestSource = max(newestSource, info.datenum);
end

if exist(cacheFile, 'file') == 2
    cacheInfo = dir(cacheFile);
    C = load(cacheFile);
    required = {'errDG','sigmaDG','NcCache','pCache','refineCache','refFileCache'};
    valid = all(isfield(C, required)) && cacheInfo.datenum >= newestSource && ...
        isequal(C.NcCache, cfg.Nc_list) && isequal(C.pCache, cfg.p_list) && ...
        C.refineCache == cfg.Refine_fixed && strcmp(C.refFileCache, refFile);
    if valid
        errDG = C.errDG;
        sigmaDG = C.sigmaDG;
        fprintf('[CACHE] %s\n', cacheFile);
        return;
    end
end

eigBlocks = {1, 2, 3:4};
[errDG, sigmaDG] = cutoff_eigenfunction_errors( ...
    runFiles, refFile, eigBlocks, cfg.dx_in, cfg.dx_out);
NcCache = cfg.Nc_list;
pCache = cfg.p_list;
refineCache = cfg.Refine_fixed;
refFileCache = refFile;
save(cacheFile, 'errDG', 'sigmaDG', 'NcCache', 'pCache', ...
    'refineCache', 'refFileCache');
fprintf('[SAVED] %s\n', cacheFile);
end

function save_plot(fig, outDir, baseName)
% Export one cutoff figure.
outPdf = fullfile(outDir, [baseName '.pdf']);
exportgraphics(fig, outPdf, 'ContentType', 'vector');
fprintf('[SAVED] %s\n', outPdf);
end

function [Nc_ok, lam_ok, dof_ok] = ...
read_runs_over_Nc(resultRoot, Nc_list, pdeg, refTag, n_eigs)

% Read the cutoff-convergence runs for one refinement level.
nNc   = numel(Nc_list);
lamNc = zeros(nNc, n_eigs);
dofNc = zeros(nNc, 1);

for i = 1:nNc
    Nc = Nc_list(i);
    runMat = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg), refTag, 'run.mat');

    assert(exist(runMat, 'file') == 2, 'Missing run file: %s', runMat);

    S = load(runMat, 'run');
    assert(isfield(S, 'run'), 'Missing run structure in %s.', runMat);

    run = S.run;
    assert(isfield(run, 'lambda'), 'Missing eigenvalues in %s.', runMat);

    lam = run.lambda(:).';
    assert(numel(lam) >= n_eigs, 'The run file has too few eigenvalues: %s', runMat);
    assert(isfield(run, 'n_dofs_total'), 'Missing total DOFs in %s.', runMat);

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

function fig = plot_semilogy_err_vs_Nc(Nc_ok, lam_ok, lambda_ref, eig_list, pdeg, cfg)
% Plot eigenvalue errors versus cutoff.
err = abs(lam_ok(:, eig_list) - lambda_ref(1, eig_list));
labels = arrayfun(@(q) sprintf('$i={%d}$', q), eig_list, 'UniformOutput', false);
fig = plot_semilogy_curves(Nc_ok, err, labels, ...
    '$|\lambda_i-\lambda_{i}^{\mathrm{DG}}|$', pdeg, cfg);
end

function fig = plot_semilogy_eigfun_vs_Nc(Nc_ok, errDG, eig_list, pdeg, cfg)
% Plot only aligned DG eigenfunction errors versus cutoff.
labels = arrayfun(@(q) sprintf('$i={%d}$', q), eig_list, 'UniformOutput', false);
cfg.yTickExp_p1 = [];
cfg.yTickExp_p2 = [];
cfg.manualYLim_p1 = [];
cfg.manualYLim_p2 = [];
fig = plot_semilogy_curves(Nc_ok, errDG, labels, ...
    '$\|u_i-u_{i}^{\mathrm{DG}}\|_{\mathrm{DG}}$', pdeg, cfg);
end

function fig = plot_semilogy_curves(Nc_ok, err, legLabels, yLabel, pdeg, cfg)
% Apply the common Example 1 cutoff-convergence plot style.
% Validate the error data and create the plotting axes.
assert(all(err(:) > 0), 'Cutoff-convergence errors must be positive.');
fig = figure( ...
    'Color',    cfg.fig.bgColor, ...
    'Units',    'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);

ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');

ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];

% Draw one styled curve for each eigenvalue error.
colors = double(cfg.lineColors255) / 255;

nC  = size(colors,1);
nLS = numel(cfg.lineStyles);
nMK = numel(cfg.markers);

hEig = gobjects(1, size(err, 2));

for j = 1:size(err, 2)
    col = colors(mod(j-1, nC) + 1, :);
    ls  = cfg.lineStyles{mod(j-1, nLS) + 1};
    mk  = cfg.markers{mod(j-1, nMK) + 1};

    hEig(j) = semilogy(ax, Nc_ok, err(:,j), ls, ...
        'LineWidth', cfg.lw, ...
        'Color', col, ...
        'Marker', mk, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', col);
end

% Apply logarithmic axes, labels, and legend styling.
set(ax, ...
    'XScale',     cfg.axes.xScale, ...
    'YScale',     cfg.axes.yScale, ...
    'FontSize',   cfg.axes.fontSize, ...
    'LineWidth',  cfg.axes.lineWidth, ...
    'TickDir',    cfg.axes.tickDir, ...
    'Box',        cfg.axes.box, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick);

ax.TickLabelInterpreter = 'latex';

grid(ax, 'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.XRuler.MinorTick = 'off';
ax.YRuler.MinorTick = 'off';

xlabel(ax, '$K$', ...
    'Interpreter', 'latex', ...
    'FontSize',    cfg.axes.labelSize);

ylabel(ax, yLabel, ...
    'Interpreter', 'latex', ...
    'FontSize',    cfg.axes.labelSize);

lgd = legend(ax, hEig, legLabels, ...
    'Location',    cfg.legend.location, ...
    'Interpreter', 'latex', ...
    'FontSize',    cfg.legend.fontSize, ...
    'NumColumns',  cfg.legend.numColumns);
lgd.Box = cfg.legend.box;

% Derive axis limits and ticks from the plotted data.
xMin = min(Nc_ok);
xMax = max(Nc_ok);
xRng = xMax - xMin;
if xRng == 0
    xRng = 1;
end

axisSpec.XTick = Nc_ok(:).';
axisSpec.XLim  = [xMin - cfg.padX_left * xRng, xMax + cfg.padX_right * xRng];

yAll = err(:);
yMin = min(yAll);
yMax = max(yAll);

if pdeg == 1
    expsInput = cfg.yTickExp_p1;
    manualYLim = cfg.manualYLim_p1;
else
    expsInput = cfg.yTickExp_p2;
    manualYLim = cfg.manualYLim_p2;
end

if ~isempty(expsInput)
    exps = sort(expsInput(:).');
    axisSpec.YTick = 10.^exps;
    axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), exps, 'UniformOutput', false);
else
    exps = floor(log10(yMin)):ceil(log10(yMax));
    axisSpec.YTick = 10.^exps;
    axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), exps, 'UniformOutput', false);
end

if numel(manualYLim) == 2 && all(manualYLim > 0)
    axisSpec.YLim = manualYLim;
else
    axisSpec.YLim = [yMin / (1 + cfg.padY_low), yMax * (1 + cfg.padY_high)];
end

ax.XTick = axisSpec.XTick;
ax.YTick = axisSpec.YTick;
ax.YTickLabel = axisSpec.YTickLabel;
ax.XLim = axisSpec.XLim;
ax.YLim = axisSpec.YLim;

end

function rate = local_exp_rates_Nc(Nc, err)
% Estimate exponential convergence rates between successive cutoffs.

Nc = Nc(:);
n  = numel(Nc);
ne = size(err, 2);

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
scale = 10.^(nSig - 1 - p);

y(mask) = round(x(mask) .* scale) ./ scale;
end
