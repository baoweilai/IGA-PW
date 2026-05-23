function out = plot_h_convergence()
%Plot h-convergence data.

clc;clear;
rootDir = fileparts(fileparts(mfilename('fullpath')));
add_workflow_paths(fullfile(rootDir, 'model', 'h_convergence'), ...
    {'config', 'core', 'operators', 'solver'});
finalDir = fullfile(rootDir, 'data', 'result', 'IGA', 'convergence');
T = readtable(fullfile(finalDir, 'values.csv'));

apply_paper_plot_style();

cfg = plot_cfg_local();
pList = [1 2];

eigFig = make_figure_local(cfg);
eigAx = axes(eigFig); hold(eigAx, 'on'); box(eigAx, 'on');
style_axes_local(eigAx, cfg);
yEigAll = [];
for ip = 1:numel(pList)
    p = pList(ip);
    [h, err] = series_local(T, p, 'eigAbsError');
    yEigAll = [yEigAll; err(:)]; %#ok<AGROW>
    plot(eigAx, h, err, '-', ...
        'LineWidth', cfg.lw, ...
        'Marker', cfg.markers{ip}, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cfg.lineColors(ip, :), ...
        'Color', cfg.lineColors(ip, :), ...
        'DisplayName', sprintf('$p=%d$', p));
end
s1 = add_slope_local(eigAx, T, 1, 'eigAbsError', 2, cfg.order2Color, ...
    '$\mathrm{Slope}=2$', 10, cfg);
s2 = add_slope_local(eigAx, T, 2, 'eigAbsError', 3, cfg.order3Color, ...
    '$\mathrm{Slope}=3$', 0.05, cfg);
yEigAll = [yEigAll; s1(:); s2(:)];
finish_axes_local(eigAx, T.h, yEigAll, cfg);
xlabel(eigAx, cfg.xlabel, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(eigAx, '$|\lambda_1-\lambda_{1}^{\mathrm{DG}}|$', ...
    'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
legend(eigAx, 'show', 'Location', 'southeast', 'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize, 'Box', cfg.legend.box);

dgFig = make_figure_local(cfg);
dgAx = axes(dgFig); hold(dgAx, 'on'); box(dgAx, 'on');
style_axes_local(dgAx, cfg);
yDgAll = [];
for ip = 1:numel(pList)
    p = pList(ip);
    [h, err] = series_local(T, p, 'dgError_g160');
    yDgAll = [yDgAll; err(:)]; %#ok<AGROW>
    plot(dgAx, h, err, '-', ...
        'LineWidth', cfg.lw, ...
        'Marker', cfg.markers{ip}, ...
        'MarkerSize', cfg.ms, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cfg.lineColors(ip, :), ...
        'Color', cfg.lineColors(ip, :), ...
        'DisplayName', sprintf('$p=%d$', p));
end
s3 = add_slope_local(dgAx, T, 1, 'dgError_g160', 1, cfg.order2Color, ...
    '$\mathrm{Slope}=1$', 5, cfg);
s4 = add_slope_local(dgAx, T, 2, 'dgError_g160', 1.5, cfg.order3Color, ...
    '$\mathrm{Slope}=1.5$', 0.3, cfg);
yDgAll = [yDgAll; s3(:); s4(:)];
finish_axes_local(dgAx, T.h, yDgAll, cfg);
xlabel(dgAx, cfg.xlabel, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(dgAx, '$\|u_1-u_{1}^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
    'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
legend(dgAx, 'show', 'Location', 'southeast', 'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize, 'Box', cfg.legend.box);

eigPng = fullfile(finalDir, 'eig.png');
eigPdf = fullfile(finalDir, 'eig.pdf');
dgPng = fullfile(finalDir, 'dg.png');
dgPdf = fullfile(finalDir, 'dg.pdf');
exportgraphics(eigFig, eigPng, 'Resolution', 600);
exportgraphics(eigFig, eigPdf, 'ContentType', 'vector');
exportgraphics(dgFig, dgPng, 'Resolution', 600);
exportgraphics(dgFig, dgPdf, 'ContentType', 'vector');

out = struct('eigFig', eigFig, 'dgFig', dgFig, ...
    'eigPng', eigPng, 'eigPdf', eigPdf, 'dgPng', dgPng, 'dgPdf', dgPdf);
end

function cfg = plot_cfg_local()
%Return plotting constants for this figure.
cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

cfg.axes.fontName   = 'Times New Roman';
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
cfg.order3Color = [239 065 067] / 255;

cfg.markers = {'o','s','^','d','x','+'};
cfg.lw      = 1.8;
cfg.ms      = 8;
cfg.orderLW = 1.8;
cfg.orderLS = '--';

cfg.xlabel = '$h$';
cfg.padX     = 0.4;
cfg.padYLow  = 2;
cfg.padYHigh = 1;
cfg.xTicks = [0.025 0.05 0.10 0.20];
cfg.xTickLabels = {'0.025','0.05','0.10','0.20'};
end

function fig = make_figure_local(cfg)
%Create a figure with the saved size.
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
set(fig, 'ToolBar', 'none', 'MenuBar', 'none');
end

function style_axes_local(ax, cfg)
%Apply axis styling.
ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontName', cfg.axes.fontName, ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick, ...
    'Box', 'on');
grid(ax, 'off');
end

function [h, err] = series_local(T, p, fieldName)
%Return one plotting series.
S = sortrows(T(T.p == p, :), 'h', 'ascend');
h = S.h;
err = S.(fieldName);
end

function y = add_slope_local(ax, T, p, fieldName, order, color, labelText, scale, cfg)
%Add a reference slope marker.
[h, err] = series_local(T, p, fieldName);
anchor = max(2, min(numel(h) - 1, round(numel(h) / 2)));
h0 = h(anchor);
y0 = err(anchor) * scale;
y = (y0 / h0 ^ order) * h .^ order;
plot(ax, h, y, cfg.orderLS, ...
    'Color', color, ...
    'LineWidth', cfg.orderLW, ...
    'DisplayName', labelText);
end

function finish_axes_local(ax, hAll, yAll, cfg)
%Finalize axis labels and limits.
xMin = min(hAll);
xMax = max(hAll);
yMin = min(yAll(yAll > 0));
yMax = max(yAll);
set(ax, ...
    'XLim', [xMin / (1 + cfg.padX), xMax * (1 + cfg.padX)], ...
    'YLim', [yMin / (1 + cfg.padYLow), yMax * (1 + cfg.padYHigh)], ...
    'XTick', cfg.xTicks, ...
    'XTickLabel', cfg.xTickLabels);
end
