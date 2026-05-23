function plot_h_convergence()
%Plot h-convergence data.

clc; close all;
format short g;

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ========================= User parameters =============================
Example  = 'Example_1';
Nc       = 30;
p_list   = [1 2];
eig_list = [1 2 3 4];
refine_list = 2:6;

ref_Nc     = 45;
ref_p      = 3;
ref_refine = 8;
refCsv = fullfile(pwd, 'result', Example, sprintf('Nc_%d', ref_Nc), ...
    sprintf('p_%d', ref_p), 'summary.csv');
assert(exist(refCsv, 'file') == 2, 'Missing reference summary: %s', refCsv);

RefT = readtable(refCsv);
RefT = RefT(RefT.refine == ref_refine, :);
assert(height(RefT) == 1, ...
    'Expected one reference row for Nc=%d, p=%d, refine=%d in %s.', ...
    ref_Nc, ref_p, ref_refine, refCsv);
lambda_ref = [RefT.lambda1(1), RefT.lambda2(1), RefT.lambda3(1), RefT.lambda4(1)];

resultRoot = fullfile(pwd, 'result', Example, sprintf('Nc_%d', Nc));

%% ========================= Unified style params ========================
cfg = struct();

cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

cfg.axes.fontSize   = 11;
cfg.axes.labelSize  = 13;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';

cfg.legend.fontSize = 11;
cfg.legend.box      = 'off';

cfg.lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;

cfg.order2Color = [033 158 188] / 255;
cfg.order4Color = [239 065 067] / 255;

cfg.markers = {'o','s','^','d','x','+'};
cfg.lw      = 1.8;
cfg.ms      = 8;
cfg.orderLW = 1.8;
cfg.orderLS = '--';

cfg.xlabel = '$h$';
cfg.ylabel = '$|\lambda_i-\lambda_{i}^{\mathrm{DG}}|$';

cfg.a = 0.2;

cfg.p2_fac2 = 0.04;
cfg.p2_fac4 = 0.05;

cfg.padX     = 0.4;
cfg.padYLow  = 2;
cfg.padYHigh = 1;

cfg.xL1  = 0.40;
cfg.xL2  = 0.50;
cfg.xLT  = 0.52;
cfg.xR1  = 0.70;
cfg.xR2  = 0.88;
cfg.xRT  = 0.90;
cfg.rowY = [0.4, 0.3, 0.2, 0.1];
cfg.fakeLegendFontSize = 10;

%% ======================= Read data from summary.csv ====================
data = struct([]);

for k = 1:numel(p_list)
    pp = p_list(k);

    pDir = fullfile(resultRoot, sprintf('p_%d', pp));
    csvFile = fullfile(pDir, 'summary.csv');

    assert(exist(csvFile, 'file') == 2, 'Missing h-convergence summary: %s', csvFile);
    T = readtable(csvFile);
    T = sortrows(T, 'refine');
    T = T(ismember(T.refine, refine_list), :);
    assert(height(T) == numel(refine_list), ...
        'Expected refinements %s in %s.', mat2str(refine_list), csvFile);

    refine = T.refine(:);
    h = (2 * cfg.a) ./ (2 .^ refine);
    lam = [T.lambda1, T.lambda2, T.lambda3, T.lambda4];

    data(k).p      = pp;
    data(k).refine = refine;
    data(k).h      = h;
    data(k).lambda = lam;
    data(k).pDir   = pDir;
end

%% ======================= Prepare output table ==========================
Rows = {};
rowCount = 0;

%% ======================= Plot and compute orders =======================
for k = 1:numel(data)
    pp  = data(k).p;
    h   = data(k).h;
    lam = data(k).lambda;
    ref = data(k).refine;

    fig = plot_one_p(h, lam, lambda_ref, pp, eig_list, cfg);

    figBase = 'h';
    outPng  = fullfile(data(k).pDir, [figBase '.png']);
    outPdf  = fullfile(data(k).pDir, [figBase '.pdf']);

    exportgraphics(fig, outPng, 'Resolution', 600);
    exportgraphics(fig, outPdf, 'ContentType', 'vector');
    drawnow;

    err = abs(lam(:, eig_list) - lambda_ref(1, eig_list));
    ord = local_orders_matrix(h, err);

    fprintf('\n================ Local Orders: Nc=%d, p=%d ================\n', Nc, pp);
    for j = 1:numel(eig_list)
        eigIdx = eig_list(j);
        fprintf('lambda_%d: ', eigIdx);
        for ii = 1:size(ord,1)
            fprintf('r%d->r%d: %.4g  ', ...
                ref(ii), ref(ii+1), round_sig(ord(ii,j), 4));
        end
        fprintf('\n');
    end

    for j = 1:numel(eig_list)
        eigIdx = eig_list(j);
        for ii = 1:size(ord,1)
            rowCount = rowCount + 1;

            err_to  = round_sig(err(ii+1, j), 8);
            ord_loc = round_sig(ord(ii, j), 8);

            Rows(rowCount, :) = {Nc, pp, eigIdx, ref(ii+1), err_to, ord_loc};
        end
    end
end

%% ======================= Save order table ==============================
OrderTable = cell2table(Rows, ...
    'VariableNames', {'Nc','p','eig','refine','err','order'});

out_csv = fullfile(resultRoot, 'orders.csv');
out_mat = fullfile(resultRoot, 'orders.mat');

writetable(OrderTable, out_csv);
save(out_mat, 'OrderTable');

fprintf('\n[Saved] %s\n', out_csv);
fprintf('[Saved] %s\n', out_mat);

end

function fig = plot_one_p(h, lam, lambda_ref, pdeg, eig_list, cfg)
%Plot one p.

err = abs(lam(:, eig_list) - lambda_ref(1, eig_list));

[h, idx] = sort(h(:), 'ascend');
err = err(idx, :);

i1 = find(eig_list == 1, 1);
i2 = find(eig_list == 2, 1);
i3 = find(eig_list == 3, 1);

%% --------------------------- Figure setup ------------------------------
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
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

%% --------------------------- Plot error lines --------------------------
hEig = gobjects(1, numel(eig_list));

for j = 1:numel(eig_list)
    c = cfg.lineColors(1 + mod(j-1, size(cfg.lineColors,1)), :);
    mk = cfg.markers{min(j, numel(cfg.markers))};

    hEig(j) = plot(ax, h, err(:,j), '-', ...
        'LineWidth', cfg.lw, ...
        'Color', c, ...
        'Marker', mk, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', c);
end

%% --------------------------- Axis style -------------------------------
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick, ...
    'Box', 'on');
grid(ax, 'off');

xlabel(ax, cfg.xlabel, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);

ylabel(ax, cfg.ylabel, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);

%% ---------------------- Reference slope lines -------------------------
n = numel(h);
anchor = max(2, min(n-1, round(n/2)));
h0 = h(anchor);

make_ref = @(q, y0) (y0 / (h0^q)) * (h.^q);

hOrder2 = gobjects(0);
hOrder4 = gobjects(0);
yref2 = [];
yref4 = [];

if pdeg == 1 && ~isempty(i2) && ~isempty(i3)
    y0 = sqrt(err(anchor, i2) * err(anchor, i3));
    yref2 = make_ref(2, y0);
    hOrder2 = plot(ax, h, yref2, cfg.orderLS, ...
        'Color', cfg.order2Color, ...
        'LineWidth', cfg.orderLW);
end

if pdeg == 2
    if ~isempty(i1) && ~isempty(i2)
        y0_2 = cfg.p2_fac2 * max(err(anchor, i1), err(anchor, i2));
        yref2 = make_ref(2, y0_2);
        hOrder2 = plot(ax, h, yref2, cfg.orderLS, ...
            'Color', cfg.order2Color, ...
            'LineWidth', cfg.orderLW);
    end

    if ~isempty(i2) && ~isempty(i3)
        y0_4 = cfg.p2_fac4 * sqrt(err(anchor, i2) * err(anchor, i3));
        yref4 = make_ref(4, y0_4);
        hOrder4 = plot(ax, h, yref4, cfg.orderLS, ...
            'Color', cfg.order4Color, ...
            'LineWidth', cfg.orderLW);
    end
end

%% --------------------------- Axis limits -------------------------------
xMin = min(h);
xMax = max(h);

yAll = err(:);
if ~isempty(yref2), yAll = [yAll; yref2(:)]; end
if ~isempty(yref4), yAll = [yAll; yref4(:)]; end

yMin = min(yAll);
yMax = max(yAll);

set(ax, ...
    'XLim', [xMin/(1+cfg.padX), xMax*(1+cfg.padX)], ...
    'YLim', [yMin/(1+cfg.padYLow), yMax*(1+cfg.padYHigh)]);

%% ------------------------ Manual y-ticks -------------------------------
if pdeg == 1
    ax.YTick = 10.^[-8 -6 -4 -2];
    ax.YTickLabel = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$'};
    yl = ax.YLim;
    ax.YLim = [min(yl(1), 1e-8), max(yl(2), 1e-2)];
end

if pdeg == 2
    ax.YTick = 10.^[-12 -9 -6 -3 0 ];
    ax.YTickLabel = {'$10^{-12}$','$10^{-9}$', ...
        '$10^{-6}$','$10^{-3}$','$10^{0}$'};

    yl = ax.YLim;
    ax.YLim = [min(yl(1), 1e-12), max(yl(2), 1e-0)];
end

%% ------------------------ Fake legend in same axes ---------------------
draw_fake_legend_in_axes(ax, pdeg, hEig, hOrder2, hOrder4, eig_list, cfg);

end

function draw_fake_legend_in_axes(ax, pdeg, hEig, hOrder2, hOrder4, eig_list, cfg)
%Draw fake legend in axes.

for j = 1:numel(eig_list)
    lbl = sprintf('$i=%d$', eig_list(j));
    draw_one_fake_entry_in_axes(ax, hEig(j), cfg.xR1, cfg.xR2, cfg.xRT, cfg.rowY(j), lbl, cfg);
end

if pdeg == 1
    if isgraphics(hOrder2)
        draw_one_fake_entry_in_axes(ax, hOrder2, cfg.xL1, cfg.xL2, cfg.xLT, cfg.rowY(4), '$\mathrm{Slope} = 2$', cfg);
    end

elseif pdeg == 2
    if isgraphics(hOrder2)
        draw_one_fake_entry_in_axes(ax, hOrder2, cfg.xL1, cfg.xL2, cfg.xLT, cfg.rowY(3), '$\mathrm{Slope} = 2$', cfg);
    end
    if isgraphics(hOrder4)
        draw_one_fake_entry_in_axes(ax, hOrder4, cfg.xL1, cfg.xL2, cfg.xLT, cfg.rowY(4), '$\mathrm{Slope} = 4$', cfg);
    end
end

end

function draw_one_fake_entry_in_axes(ax, hLine, x1n, x2n, xtn, yn, labelStr, cfg)
%Draw one fake entry in axes.

c  = get(hLine, 'Color');
ls = get(hLine, 'LineStyle');
lw = get(hLine, 'LineWidth');

mk  = 'none';
ms  = 6;
mfc = 'none';
mec = c;

if isprop(hLine, 'Marker')
    mk = get(hLine, 'Marker');
end
if isprop(hLine, 'MarkerSize')
    ms = get(hLine, 'MarkerSize');
end
if isprop(hLine, 'MarkerFaceColor')
    mfc = get(hLine, 'MarkerFaceColor');
end
if isprop(hLine, 'MarkerEdgeColor')
    mec = get(hLine, 'MarkerEdgeColor');
end

[x1, y1] = axes_norm_to_data(ax, x1n, yn);
[x2, y2] = axes_norm_to_data(ax, x2n, yn);
[xt, yt] = axes_norm_to_data(ax, xtn, yn);

xm = 0.49 * (x1 + x2);
ym = y1;

line(ax, [x1 x2], [y1 y2], ...
    'LineStyle', ls, ...
    'Color', c, ...
    'LineWidth', lw, ...
    'Marker', 'none', ...
    'Clipping', 'off');

line(ax, xm, ym, ...
    'LineStyle', 'none', ...
    'Color', c, ...
    'Marker', mk, ...
    'MarkerSize', ms, ...
    'MarkerFaceColor', mfc, ...
    'MarkerEdgeColor', mec, ...
    'LineWidth', lw, ...
    'Clipping', 'off');

text(ax, xt, yt, labelStr, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.fakeLegendFontSize, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', ...
    'Margin', 0.5, ...
    'Clipping', 'off');

end

function [x, y] = axes_norm_to_data(ax, xn, yn)
%Compute norm to data.

xlimv = ax.XLim;
ylimv = ax.YLim;

if strcmpi(ax.XScale, 'log')
    lx1 = log10(xlimv(1));
    lx2 = log10(xlimv(2));
    x = 10^(lx1 + xn * (lx2 - lx1));
else
    x = xlimv(1) + xn * (xlimv(2) - xlimv(1));
end

if strcmpi(ax.YScale, 'log')
    ly1 = log10(ylimv(1));
    ly2 = log10(ylimv(2));
    y = 10^(ly1 + yn * (ly2 - ly1));
else
    y = ylimv(1) + yn * (ylimv(2) - ylimv(1));
end
end

function ord = local_orders_matrix(h, err)
%Compute orders matrix.

n  = numel(h);
ne = size(err, 2);

ord = zeros(n-1, ne);

for i = 1:n-1
    ord(i,:) = log(err(i+1,:) ./ err(i,:)) ./ log(h(i+1) / h(i));
end

end

function y = round_sig(x, nSig)
%Round a value to significant digits.

y = x;

mask = isfinite(x) & (x ~= 0);
ax = abs(x(mask));
p = floor(log10(ax));
scale = 10.^(nSig - 1 - p);

y(mask) = round(x(mask) .* scale) ./ scale;

end
