function plot_cutoff_convergence()
%Plot cutoff-convergence data.

clc; close all;
format short g;

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

cfg.savePNG = true;
cfg.savePDF = true;
cfg.saveFIG = false;
cfg.pngDPI  = 600;

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

% -------------------- grid style --------------------
cfg.gridOn      = false;
cfg.minorGridOn = false;

% -------------------- curve style --------------------
cfg.lineColors255 = [ ...
    223 122  94;
    60  64  91;
    130 178 154;
    242 204 142
    ];

cfg.lineStyles   = {'-','-','-','-'};
cfg.markers      = {'o','s','^','d','x','+'};
cfg.markerFilled = false;
cfg.lw           = 1.8;
cfg.ms           = 8;

% -------------------- legend style --------------------
cfg.legend.location   = 'northeast';
cfg.legend.fontSize   = 11;
cfg.legend.box        = 'off';
cfg.legend.numColumns = 1;

% -------------------- y-axis manual control --------------------
cfg.useManualYTicks = true;
cfg.yTickExp_p1 = [-1 -2 -3 -4 -5 -6 -7];
cfg.yTickExp_p2 = [-1 -2 -3 -4 -5 -6 -7];

cfg.useManualYLim = true;
cfg.manualYLim_p1 = [1e-7, 1e-1];
cfg.manualYLim_p2 = [1e-7, 1e-1];

% -------------------- plot padding --------------------
cfg.padX_left  = 0.10;
cfg.padX_right = 0.10;
cfg.padY_low   = 1;
cfg.padY_high  = 1;

% 2. Build paths

resultRoot = fullfile(pwd, 'result', cfg.Example);
if ~exist(resultRoot, 'dir')
    error('Cannot find resultRoot: %s', resultRoot);
end

outRoot = fullfile(resultRoot, cfg.outSubDirName, cfg.refTag);
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end

cfg.Nc_list = cfg.Nc_list(:).';

% 3. Main loop: fixed refine, sweep Nc

for pp = cfg.p_list
    fprintf('\n==================== [P=%d] fixed refine=%d, sweep Nc ====================\n', ...
        pp, cfg.Refine_fixed);

    [Nc_ok, lam_ok, dof_ok] = ...
        read_runs_over_Nc(resultRoot, cfg.Nc_list, pp, cfg.refTag, cfg.n_eigenvalues);

    lambda_ref = cfg.lambda_ref;

    fig = plot_semilogy_err_vs_Nc(Nc_ok, lam_ok, lambda_ref, cfg.eig_list, pp, cfg);

    pOutRoot = fullfile(outRoot, sprintf('p_%d', pp));
    if ~exist(pOutRoot, 'dir'), mkdir(pOutRoot); end
    figBase = 'cutoff';

    if cfg.savePNG
        outPng = fullfile(pOutRoot, [figBase '.png']);
        exportgraphics(fig, outPng, 'Resolution', cfg.pngDPI);
        fprintf('[SAVED] %s\n', outPng);
    end

    if cfg.savePDF
        outPdf = fullfile(pOutRoot, [figBase '.pdf']);
        exportgraphics(fig, outPdf, 'ContentType', 'vector');
        fprintf('[SAVED] %s\n', outPdf);
    end

    if cfg.saveFIG
        outFig = fullfile(pOutRoot, [figBase '.fig']);
        saveas(fig, outFig);
        fprintf('[SAVED] %s\n', outFig);
    end

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
%Plot semilogy err vs nc.

err = abs(lam_ok(:, eig_list) - lambda_ref(1, eig_list));
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

colors = double(cfg.lineColors255) / 255;

nC  = size(colors,1);
nLS = numel(cfg.lineStyles);
nMK = numel(cfg.markers);

hEig = gobjects(1, numel(eig_list));

for j = 1:numel(eig_list)
    col = colors(mod(j-1, nC) + 1, :);
    ls  = cfg.lineStyles{mod(j-1, nLS) + 1};
    mk  = cfg.markers{mod(j-1, nMK) + 1};

    if cfg.markerFilled
        mfc = col;
    else
        mfc = 'w';
    end

    hEig(j) = semilogy(ax, Nc_ok, err(:,j), ls, ...
        'LineWidth', cfg.lw, ...
        'Color', col, ...
        'Marker', mk, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', mfc, ...
        'MarkerEdgeColor', col);
end

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

ylabel(ax, '$|\lambda_i-\lambda_{i}^{\mathrm{DG}}|$', ...
    'Interpreter', 'latex', ...
    'FontSize',    cfg.axes.labelSize);

legLabels = arrayfun(@(q) sprintf('$i={%d}$', q), eig_list, 'UniformOutput', false);

lgd = legend(ax, hEig, legLabels, ...
    'Location',    cfg.legend.location, ...
    'Interpreter', 'latex', ...
    'FontSize',    cfg.legend.fontSize, ...
    'NumColumns',  cfg.legend.numColumns);
lgd.Box = cfg.legend.box;

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

if cfg.useManualYTicks
    exps = sort(expsInput(:).');
    axisSpec.YTick = 10.^exps;
    axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), exps, 'UniformOutput', false);
else
    exps = floor(log10(yMin)):ceil(log10(yMax));
    axisSpec.YTick = 10.^exps;
    axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), exps, 'UniformOutput', false);
end

if cfg.useManualYLim && numel(manualYLim) == 2 && all(manualYLim > 0)
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
%Compute exp rates nc.

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
%Round a value to significant digits.

y = x;
mask = isfinite(x) & (x ~= 0);

ax = abs(x(mask));
p  = floor(log10(ax));
scale = 10.^(nSig - 1 - p);

y(mask) = round(x(mask) .* scale) ./ scale;
end

function t = onoff(flag)
%Convert a logical value to on/off text.
if flag
    t = 'on';
else
    t = 'off';
end
end
