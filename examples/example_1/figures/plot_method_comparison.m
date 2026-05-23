function plot_method_comparison(optsIn)
%Plot method-comparison data.

assert(exist('optsIn', 'var') == 1, 'plot_method_comparison requires an options structure.');

exampleDir = fileparts(fileparts(mfilename('fullpath')));
rootDir = fullfile(exampleDir, 'data');
outRoot = get_opt_local(optsIn, 'outRoot', fullfile(rootDir, ...
    'method_comparison'));
plotPointsFile = get_opt_local(optsIn, 'plotPointsFile', ...
    fullfile(outRoot, 'summary_plot_points.csv'));
figDir = get_opt_local(optsIn, 'figDir', fullfile(outRoot, 'figures'));

assert(exist(plotPointsFile, 'file') == 2, ...
    'Missing plot point table: %s', plotPointsFile);
ensure_dir_local(figDir);

T = read_plot_table_local(plotPointsFile);
requiredVars = {'method', 'err1', 'dof', 'total_time_s', 'figure'};
missingVars = setdiff(requiredVars, T.Properties.VariableNames);
assert(isempty(missingVars), 'Missing columns in %s: %s', ...
    plotPointsFile, strjoin(missingVars, ', '));

T.method = string(T.method);
T.figure = string(T.figure);

Tdof = T(T.figure == "dof", :);
Ttime = T(T.figure == "time", :);
assert(~isempty(Tdof), 'No rows marked figure=dof in %s.', plotPointsFile);
assert(~isempty(Ttime), 'No rows marked figure=time in %s.', plotPointsFile);

plot_method_compare_local(Tdof, Ttime, figDir);

fprintf('[SAVED] %s\n', fullfile(figDir, 'dof.pdf'));
fprintf('[SAVED] %s\n', fullfile(figDir, 'time.pdf'));
end

function T = read_plot_table_local(plotPointsFile)
%Read plot table.
T = readtable(plotPointsFile, 'TextType', 'string');
end

function plot_method_compare_local(Tdof, Ttime, figDir)
%Plot method compare.
cfg = default_style_local();
methods = {'PW', 'IGA (p=1)', 'IGA (p=2)', 'IGA-PW'};
labels = {'Plane wave', 'IGA, $p=1$', 'IGA, $p=2$', 'IGA-PW'};
colors = cfg.lineColors([1 2 3 4], :);
markers = cfg.markers(1:numel(methods));

plot_one_compare_local(Tdof, methods, labels, colors, markers, cfg, ...
    'dof', 'Degrees of Freedom', ...
    fullfile(figDir, 'dof'));
plot_one_compare_local(Ttime, methods, labels, colors, markers, cfg, ...
    'total_time_s', 'Total Time', ...
    fullfile(figDir, 'time'));
end

function plot_one_compare_local(T, methods, labels, colors, markers, cfg, xField, xLabelText, baseName)
%Plot one compare.
errLabel = '$|\lambda_{1}-\lambda_{1}^{\mathrm{DG}}|$';
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
ax = axes(fig);
set_axes_layout_local(ax, cfg);
h = gobjects(numel(methods), 1);
hold(ax, 'on');
for k = 1:numel(methods)
    Tk = T(T.method == methods{k}, :);
    Tk = sortrows(Tk, xField);
    h(k) = plot(ax, Tk.(xField), Tk.err1, '-', ...
        'LineWidth', cfg.lineWidth, ...
        'Color', colors(k, :), ...
        'Marker', markers{k}, ...
        'MarkerSize', cfg.markerSize, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', colors(k, :), ...
        'DisplayName', labels{k});
end
style_axes_local(ax, cfg);
xlabel(ax, xLabelText, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(ax, errLabel, 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
apply_log_limits_local(ax, T.(xField), T.err1, cfg.padX, cfg.padYLow, cfg.padYHigh);
lgd = legend(ax, ...
    'Location', 'southwest', ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize);
lgd.Box = 'off';
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [baseName '.png'], 'Resolution', 600);
close(fig);
end

function style_axes_local(ax, cfg)
%Apply axis styling.
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'Box', 'on', ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick);
grid(ax, 'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.XRuler.MinorTick = 'off';
ax.YRuler.MinorTick = 'off';
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end

function set_axes_layout_local(ax, cfg)
%Apply axes layout settings.
ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
end

function apply_log_limits_local(ax, x, y, padX, padYLow, padYHigh)
%Apply log limits.
x = x(isfinite(x) & x > 0);
y = y(isfinite(y) & y > 0);
set(ax, 'XLim', [min(x) / (1 + padX), max(x) * (1 + padX)]);
set(ax, 'YLim', [min(y) / (1 + padYLow), max(y) * (1 + padYHigh)]);
end

function cfg = default_style_local()
%Return plotting style values.
cfg = struct();
cfg.fig.width = 4.8;
cfg.fig.height = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor = 'w';
cfg.layout.left = 0.14;
cfg.layout.right = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top = 0.08;
cfg.axes.fontSize = 10;
cfg.axes.lineWidth = 1.0;
cfg.axes.tickDir = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.labelSize = 12;
cfg.legend.fontSize = 11;
cfg.lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;
cfg.markers = {'o', 's', '^', 'd', 'x', '+'};
cfg.lineWidth = 1.8;
cfg.markerSize = 8;
cfg.padX = 0.4;
cfg.padYLow = 2;
cfg.padYHigh = 1;
end

function val = get_opt_local(s, fieldName, defaultVal)
%Return one option value.
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function ensure_dir_local(p)
%Create a directory when it is missing.
if ~exist(p, 'dir')
    mkdir(p);
end
end
