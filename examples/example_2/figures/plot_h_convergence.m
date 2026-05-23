function plot_h_convergence()
%Plot h-convergence data.

clc; close all; format short g;

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ========================= User parameters =============================
Example  = 'Example_2';
Nc       = 30;
p        = 2;
eig_list = [1 2 3];
refines_to_plot = 2:6;

refNc     = 45;
refP      = 3;
refRefine = 8;

resultRoot = fullfile(pwd, 'result', Example, sprintf('Nc_%d', Nc));
pDir       = fullfile(resultRoot, sprintf('p_%d', p));
csvFile    = fullfile(pDir, 'summary.csv');
refCsvFile = fullfile(pwd, 'result', Example, ...
    sprintf('Nc_%d', refNc), sprintf('p_%d', refP), 'summary.csv');

%% ========================= Style parameters ============================
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

cfg.lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154] / 255;

cfg.order2Color = [033 158 188] / 255;
cfg.order4Color = [239 065 067] / 255;

cfg.markers = {'o','s','^'};
cfg.lw      = 1.8;
cfg.ms      = 8;
cfg.orderLW = 1.8;
cfg.orderLS = '--';

cfg.xlabel = '$h$';
cfg.ylabel = '$|\lambda_i-\lambda_{i}^{\mathrm{DG}}|$';

cfg.p2_fac2 = 0.02;
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
if ~exist(csvFile, 'file')
    error('Missing file: %s', csvFile);
end
if ~exist(refCsvFile, 'file')
    error('Missing reference file: %s', refCsvFile);
end

Tref = readtable(refCsvFile);
Tref = Tref(Tref.refine == refRefine, :);
assert(height(Tref) == 1, ...
    'Expected exactly one reference row with refine=%d in %s.', refRefine, refCsvFile);
lambda_ref = [Tref.lambda1, Tref.lambda2, Tref.lambda3, Tref.lambda4];

T = readtable(csvFile);
T = sortrows(T, 'refine');
T = T(ismember(T.refine, refines_to_plot), :);
assert(height(T) == numel(refines_to_plot), ...
    'Expected h-convergence rows refine=%s in %s.', mat2str(refines_to_plot), csvFile);

refine = T.refine(:);
h      = 0.4 ./ (2.^refine);
lam    = [T.lambda1, T.lambda2, T.lambda3, T.lambda4];

%% ======================= Plot and compute orders =======================
fig = plot_one_p2(h, lam, lambda_ref, eig_list, cfg);

figBase = 'h';
exportgraphics(fig, fullfile(pDir, [figBase '.png']), 'Resolution', 600);
exportgraphics(fig, fullfile(pDir, [figBase '.pdf']), 'ContentType', 'vector');

err = abs(lam(:, eig_list) - lambda_ref(1, eig_list));
ord = local_orders_matrix(h, err);

fprintf('\n================ Local Orders: Nc=%d, p=%d ================\n', Nc, p);
for j = 1:numel(eig_list)
    fprintf('lambda_%d: ', eig_list(j));
    for ii = 1:size(ord,1)
        fprintf('r%d->r%d: %.4g  ', refine(ii), refine(ii+1), round_sig(ord(ii,j), 4));
    end
    fprintf('\n');
end

Rows = {};
rowCount = 0;

for j = 1:numel(eig_list)
    for ii = 1:size(ord,1)
        rowCount = rowCount + 1;
        Rows(rowCount, :) = {Nc, p, eig_list(j), refine(ii+1), ...
            round_sig(err(ii+1,j), 8), ...
            round_sig(ord(ii,j), 8)};
    end
end

OrderTable = cell2table(Rows, ...
    'VariableNames', {'Nc','p','eig','refine','err','order'});

out_csv = fullfile(resultRoot, 'orders.csv');
out_mat = fullfile(resultRoot, 'orders.mat');

writetable(OrderTable, out_csv);
save(out_mat, 'OrderTable');

fprintf('\n[Saved] %s\n', out_csv);
fprintf('[Saved] %s\n', out_mat);

end

function fig = plot_one_p2(h, lam, lambda_ref, eig_list, cfg)
%Plot one p2.

err = abs(lam(:, eig_list) - lambda_ref(1, eig_list));
[h, idx] = sort(h(:), 'ascend');
err = err(idx, :);

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

hEig = gobjects(1, numel(eig_list));
for j = 1:numel(eig_list)
    hEig(j) = plot(ax, h, err(:,j), '-', ...
        'LineWidth', cfg.lw, ...
        'Color', cfg.lineColors(j,:), ...
        'Marker', cfg.markers{j}, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cfg.lineColors(j,:));
end

set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'Box', 'on');

grid(ax, 'off');

xlabel(ax, cfg.xlabel, 'FontSize', cfg.axes.labelSize);
ylabel(ax, cfg.ylabel, 'FontSize', cfg.axes.labelSize);

anchor = max(2, min(numel(h)-1, round(numel(h)/2)));
h0 = h(anchor);
make_ref = @(q, y0) (y0 / h0^q) * h.^q;

yref2 = make_ref(2, cfg.p2_fac2 * max(err(anchor,1), err(anchor,2)));
yref4 = make_ref(4, cfg.p2_fac4 * sqrt(err(anchor,2) * err(anchor,3)));

hOrder2 = plot(ax, h, yref2, cfg.orderLS, ...
    'Color', cfg.order2Color, ...
    'LineWidth', cfg.orderLW);

hOrder4 = plot(ax, h, yref4, cfg.orderLS, ...
    'Color', cfg.order4Color, ...
    'LineWidth', cfg.orderLW);

yAll = [err(:); yref2(:); yref4(:)];

set(ax, ...
    'XLim', [min(h)/(1+cfg.padX), max(h)*(1+cfg.padX)], ...
    'YLim', [min(yAll)/(1+cfg.padYLow), max(yAll)*(1+cfg.padYHigh)]);

ax.YTick = 10.^[-12 -9 -6 -3 0];
ax.YTickLabel = {'$10^{-12}$','$10^{-9}$','$10^{-6}$','$10^{-3}$','$10^{0}$'};

yl = ax.YLim;
ax.YLim = [min(yl(1), 1e-12), max(yl(2), 1e0)];

draw_fake_entry(ax, hOrder2, cfg.xL1, cfg.xL2, cfg.xLT, cfg.rowY(3), '$\mathrm{Slope} = 2$', cfg);
draw_fake_entry(ax, hOrder4, cfg.xL1, cfg.xL2, cfg.xLT, cfg.rowY(4), '$\mathrm{Slope} = 4$', cfg);

for j = 1:numel(eig_list)
    draw_fake_entry(ax, hEig(j), cfg.xR1, cfg.xR2, cfg.xRT, cfg.rowY(j+1), ...
        sprintf('$i={%d}$', eig_list(j)), cfg);
end

end

function draw_fake_entry(ax, hLine, x1n, x2n, xtn, yn, labelStr, cfg)
%Draw fake entry.

c  = get(hLine, 'Color');
ls = get(hLine, 'LineStyle');
lw = get(hLine, 'LineWidth');
mk = get(hLine, 'Marker');

[x1, y1] = axes_norm_to_data(ax, x1n, yn);
[x2, ~ ] = axes_norm_to_data(ax, x2n, yn);
[xt, yt] = axes_norm_to_data(ax, xtn, yn);

xm = 0.49 * (x1 + x2);

line(ax, [x1 x2], [y1 y1], ...
    'LineStyle', ls, ...
    'Color', c, ...
    'LineWidth', lw, ...
    'Marker', 'none', ...
    'Clipping', 'off');

line(ax, xm, y1, ...
    'LineStyle', 'none', ...
    'Color', c, ...
    'Marker', mk, ...
    'MarkerSize', cfg.ms, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', c, ...
    'LineWidth', lw, ...
    'Clipping', 'off');

text(ax, xt, yt, labelStr, ...
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

x = 10^(log10(xlimv(1)) + xn * diff(log10(xlimv)));
y = 10^(log10(ylimv(1)) + yn * diff(log10(ylimv)));

end

function ord = local_orders_matrix(h, err)
%Compute orders matrix.

ord = log(err(2:end,:) ./ err(1:end-1,:)) ./ log(h(2:end) ./ h(1:end-1));

end

function y = round_sig(x, nSig)
%Round a value to significant digits.

y = x;
mask = isfinite(x) & (x ~= 0);
p = floor(log10(abs(x(mask))));
scale = 10.^(nSig - 1 - p);
y(mask) = round(x(mask) .* scale) ./ scale;

end
