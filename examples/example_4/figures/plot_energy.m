function out = plot_energy(userCfg)
%Plot energy-convergence data.
assert(exist('userCfg', 'var') == 1, 'plot_energy requires userCfg.');

clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
add_workflow_paths(fullfile(rootDir, 'model', 'energy'), ...
    {'config', 'core', 'operators', 'solver'});

apply_paper_plot_style();

cfg.fig.width = 4.8; cfg.fig.height = 3.0;
cfg.fig.renderer = 'painters'; cfg.fig.bgColor = 'w';
cfg.layout.left = 0.14; cfg.layout.right = 0.04;
cfg.layout.bottom = 0.16; cfg.layout.top = 0.08;
cfg.axes.fontSize = 10; cfg.axes.lineWidth = 1.0;
cfg.axes.tickDir = 'out'; cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off'; cfg.axes.labelSize = 12;
cfg.legend.fontSize = 11;
cfg.lineColors = [223 122 094; 060 064 091; 130 178 154; 242 204 142] / 255;
cfg.padX = 0.2;
cfg.padY = 0.12;
cfg.labelYOffset = 0.1;

cfgRun = default_config(userCfg);
H = example_helpers(cfgRun);

referenceEnergy = load_reference_energy_local(H.cfg.resultRoot);
pList = [1 2];
KList = [2 4 8 12 16];
refineList = [1 2 3 4 5];
yTickList = [22.357 22.358 22.359 22.360 22.361 22.362];
if isfield(H.cfg, 'energy')
    if isfield(H.cfg.energy, 'p_list')
        pList = reshape(H.cfg.energy.p_list, 1, []);
    end
    if isfield(H.cfg.energy, 'K_list')
        KList = reshape(H.cfg.energy.K_list, 1, []);
    end
    if isfield(H.cfg.energy, 'refine_list')
        refineList = reshape(H.cfg.energy.refine_list, 1, []);
    end
end
if numel(KList) ~= numel(refineList)
    error('energy.K_list and energy.refine_list must have the same length.');
end

nP = numel(pList);
nCase = numel(KList);
energyTotal = zeros(nP, nCase);
hListByP = zeros(nP, nCase);

for ip = 1:nP
    for i = 1:nCase
        runFile = H.refinement_case_run_file('energy', pList(ip), KList(i), refineList(i));
        if ~exist(runFile, 'file')
            error('Energy run file not found: %s', runFile);
        end
        S = load(runFile, 'run');
        if isfield(S.run, 'meta') && isfield(S.run.meta, 'energy_total')
            energyTotal(ip, i) = S.run.meta.energy_total;
        else
            error('energy_total not found in run.meta: %s', runFile);
        end
        if isfield(S.run, 'meta') && isfield(S.run.meta, 'hmin')
            hListByP(ip, i) = S.run.meta.hmin;
        else
            error('hmin not found in run.meta: %s', runFile);
        end
    end
end
hList = hListByP(1, :);

fig = figure('Color', cfg.fig.bgColor, 'Renderer', cfg.fig.renderer, ...
    'Units', 'inches', 'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'PaperPositionMode', 'auto');
ax = axes('Parent', fig); hold(ax, 'on'); box(ax, 'on');
ax.XScale = 'log';
ax.Position = [cfg.layout.left, cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
set(ax, 'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick, ...
    'Box', 'on', ...
    'Layer', 'top');
grid(ax, 'off'); ax.XMinorGrid = 'off'; ax.YMinorGrid = 'off';
ax.XRuler.MinorTick = 'off'; ax.YRuler.MinorTick = 'off';

for ip = 1:nP
    semilogx(ax, hList, energyTotal(ip, :), '-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', 'w', ...
        'Color', cfg.lineColors(ip, :), ...
        'DisplayName', sprintf('$p=%d$', pList(ip)));
end
yline(ax, referenceEnergy, '--', ...
    'LineWidth', 1.5, ...
    'Color', cfg.lineColors(3, :), ...
    'DisplayName', 'Reference');
yTickRange = max(yTickList) - min(yTickList);
for i = 1:nCase
    labelText = sprintf('$(%d,%d)$', KList(i), refineList(i));
    yText = min(energyTotal(:, i)) - cfg.labelYOffset * yTickRange;
    text(ax, hList(i), yText, labelText, ...
        'Interpreter', 'latex', ...
        'FontSize', cfg.axes.fontSize - 1, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Clipping', 'off');
end
xlabel(ax, '$h$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(ax, 'Energy', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ax.XTick = sort(hList);
ax.XLim = [min(hList), max(hList)];
ax.YTick = yTickList;
ax.YLim = [min(yTickList), max(yTickList)];
xl = ax.XLim;
ax.XLim = [xl(1) / (1 + cfg.padX), xl(2) * (1 + cfg.padX)];
yr = max(yTickList) - min(yTickList);
ax.YLim = [min(yTickList) - cfg.padY * yr, max(yTickList) + cfg.padY * yr];
ytickformat(ax, '%.3f');
draw_energy_inset_local(ax, hList, energyTotal, referenceEnergy, ...
    KList, refineList, cfg);
lgd = legend(ax, 'show', 'Location', 'northwest', 'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize);
lgd.Box = 'off';

pngFile = fullfile(H.cfg.plotRoot, 'energy.png');
pdfFile = fullfile(H.cfg.plotRoot, 'energy.pdf');
exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

out = struct('fig', fig, 'png', pngFile, 'pdf', pdfFile, ...
    'pList', pList, 'KList', KList, 'refineList', refineList, ...
    'hList', hList, 'energyTotal', energyTotal, ...
    'referenceEnergy', referenceEnergy);
end

function referenceEnergy = load_reference_energy_local(resultRoot)
%Load the reference energy value.
refFile = fullfile(resultRoot, 'REFERENCE', 'K_45', 'p_2', ...
    'nelem_32', 'run.mat');
referenceEnergy = double(h5read(refFile, '/run/meta/energy_total'));
end

function draw_energy_inset_local(axMain, hList, energyTotal, referenceEnergy, ...
KList, refineList, cfg)
target = find(KList == 16 & refineList == 5, 1, 'first');
assert(~isempty(target), 'Missing the energy inset target case.');
idx = max(1, target - 1):target;
xData = hList(idx);
yData = energyTotal(:, idx);
yWindow = [yData(:); referenceEnergy];
ySpan = max(yWindow) - min(yWindow);
if ySpan <= eps(max(abs(yWindow)))
    ySpan = max(1, max(abs(yWindow))) * 1e-5;
end

xPad = 1.15;
zoomYLim = [min(yWindow) - 0.20 * ySpan, max(yWindow) + 0.20 * ySpan];
mainYLim = ylim(axMain);
boxYLim = [max(zoomYLim(1), mainYLim(1)), min(zoomYLim(2), mainYLim(2))];
if boxYLim(1) >= boxYLim(2)
    boxYLim = zoomYLim;
end
xLim = [min(xData) / xPad, max(xData) * xPad];

fig = ancestor(axMain, 'figure');
drawnow;
[xLeft, yTop] = energy_data_to_fig_norm_local(axMain, xLim(1), boxYLim(2));
[xRight, ~] = energy_data_to_fig_norm_local(axMain, xLim(2), boxYLim(2));
oldUnits = axMain.Units;
axMain.Units = 'normalized';
mainPos = axMain.Position;
axMain.Units = oldUnits;
insetWidth = min(max(xRight - xLeft, 0.28 * mainPos(3)), 0.34 * mainPos(3));
insetHeight = 0.26 * mainPos(4);
insetLeft = 0.5 * (xLeft + xRight) - 0.5 * insetWidth + 0.18 * mainPos(3);
insetLeft = min(max(insetLeft, mainPos(1) + 0.02 * mainPos(3)), ...
    mainPos(1) + mainPos(3) - insetWidth - 0.02 * mainPos(3));
insetBottom = yTop + 0.23 * mainPos(4);
insetBottom = min(insetBottom, mainPos(2) + mainPos(4) - insetHeight - 0.05 * mainPos(4));
insetPos = [insetLeft, insetBottom, insetWidth, insetHeight];
axInset = axes('Parent', fig, 'Position', insetPos, 'Color', 'w');
box(axInset, 'on'); hold(axInset, 'on');

for ip = 1:size(energyTotal, 1)
    semilogx(axInset, hList, energyTotal(ip, :), '-o', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4.5, ...
        'MarkerFaceColor', 'w', ...
        'Color', cfg.lineColors(ip, :));
end
yline(axInset, referenceEnergy, '--', ...
    'LineWidth', 1.0, ...
    'Color', cfg.lineColors(3, :));

set(axInset, ...
    'XScale', 'log', ...
    'XLim', xLim, ...
    'YLim', zoomYLim, ...
    'FontSize', cfg.axes.fontSize - 2, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'Box', 'on');
grid(axInset, 'off');
axInset.XRuler.MinorTick = 'off';
axInset.YRuler.MinorTick = 'off';
ytickformat(axInset, '%.4f');

drawnow;
draw_energy_zoom_box_and_connectors_local(axMain, axInset, xLim, boxYLim);
uistack(axInset, 'top');
end

function draw_energy_zoom_box_and_connectors_local(axMain, axInset, xLim, yLim)
%Draw the energy inset connectors.
oldUnits = axInset.Units;
axInset.Units = 'normalized';
insetPos = axInset.Position;
axInset.Units = oldUnits;

[xInsetLeft, yInsetBottom] = energy_fig_norm_to_data_local( ...
    axMain, insetPos(1), insetPos(2));
[xInsetRight, ~] = energy_fig_norm_to_data_local( ...
    axMain, insetPos(1) + insetPos(3), insetPos(2));

line(axMain, [xLim(1) xLim(2) xLim(2) xLim(1) xLim(1)], ...
    [yLim(1) yLim(1) yLim(2) yLim(2) yLim(1)], ...
    'Color', [0.70 0.70 0.70], ...
    'LineStyle', '-', ...
    'LineWidth', 0.7, ...
    'HandleVisibility', 'off', ...
    'Clipping', 'on');
line(axMain, [xLim(1), xInsetLeft], [yLim(2), yInsetBottom], ...
    'Color', [0.70 0.70 0.70], ...
    'LineStyle', '-', ...
    'LineWidth', 0.7, ...
    'HandleVisibility', 'off', ...
    'Clipping', 'on');
line(axMain, [xLim(2), xInsetRight], [yLim(2), yInsetBottom], ...
    'Color', [0.70 0.70 0.70], ...
    'LineStyle', '-', ...
    'LineWidth', 0.7, ...
    'HandleVisibility', 'off', ...
    'Clipping', 'on');
end

function [xf, yf] = energy_data_to_fig_norm_local(ax, x, y)
%Map energy data coordinates to figure coordinates.
[xn, yn] = energy_data_to_axes_norm_local(ax, x, y);
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
xf = pos(1) + xn * pos(3);
yf = pos(2) + yn * pos(4);
end

function [xn, yn] = energy_data_to_axes_norm_local(ax, x, y)
%Map energy data coordinates to axes coordinates.
xLim = ax.XLim;
yLim = ax.YLim;
if strcmp(ax.XScale, 'log')
    xn = (log10(x) - log10(xLim(1))) / (log10(xLim(2)) - log10(xLim(1)));
else
    xn = (x - xLim(1)) / (xLim(2) - xLim(1));
end
yn = (y - yLim(1)) / (yLim(2) - yLim(1));
xn = min(max(xn, 0), 1);
yn = min(max(yn, 0), 1);
end

function [x, y] = energy_fig_norm_to_data_local(ax, xf, yf)
%Map figure coordinates back to energy data.
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
xn = (xf - pos(1)) / pos(3);
yn = (yf - pos(2)) / pos(4);
xLim = ax.XLim;
yLim = ax.YLim;
if strcmp(ax.XScale, 'log')
    x = 10 .^ (log10(xLim(1)) + xn * (log10(xLim(2)) - log10(xLim(1))));
else
    x = xLim(1) + xn * (xLim(2) - xLim(1));
end
y = yLim(1) + yn * (yLim(2) - yLim(1));
end
