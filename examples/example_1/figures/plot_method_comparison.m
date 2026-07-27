function plot_method_comparison(optsIn)
% Plot method-comparison data.

assert(exist('optsIn', 'var') == 1, 'plot_method_comparison requires an options structure.');

exampleDir = fileparts(fileparts(mfilename('fullpath')));
rootDir = fullfile(exampleDir, 'data');
outRoot = get_opt_local(optsIn, 'outRoot', fullfile(rootDir, ...
    'method_comparison'));
plotPointsFile = get_opt_local(optsIn, 'plotPointsFile', ...
    fullfile(outRoot, 'summary_plot_points.csv'));
summaryFile = get_opt_local(optsIn, 'summaryFile', ...
    fullfile(outRoot, 'summary_method_compare.mat'));
figDir = get_opt_local(optsIn, 'figDir', fullfile(outRoot, 'figures'));

ensure_dir_local(figDir);

T = read_plot_table_local(plotPointsFile, summaryFile);
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
[Tdof, Ttime] = use_time_igapw_points_local(Tdof, Ttime);

plot_method_compare_local(Tdof, Ttime, figDir);

fprintf('[SAVED] %s\n', fullfile(figDir, 'dof.pdf'));
fprintf('[SAVED] %s\n', fullfile(figDir, 'time.pdf'));
end

function T = read_plot_table_local(plotPointsFile, summaryFile)
% Read plot table.
if exist(plotPointsFile, 'file') == 2
    T = readtable(plotPointsFile, 'TextType', 'string');
    return;
end

assert(exist(summaryFile, 'file') == 2, ...
    'Missing plot point table and summary file: %s; %s', ...
    plotPointsFile, summaryFile);
S = load(summaryFile, 'summary');
assert(isfield(S, 'summary') && isfield(S.summary, 'plot_points'), ...
    'Missing summary.plot_points in %s.', summaryFile);
T = S.summary.plot_points;
end

function [Tdof, Ttime] = use_time_igapw_points_local(Tdof, Ttime)
% Both figures use the time-selected IGA-PW configurations.
isIgapwDof = Tdof.method == "IGA-PW";
isIgapwTime = Ttime.method == "IGA-PW";
TigapwTime = Ttime(Ttime.method == "IGA-PW", :);
assert(height(TigapwTime) == 6, ...
    'Expected six time-selected IGA-PW points, found %d.', height(TigapwTime));
removePoint = TigapwTime.Nc == 10 ...
    & TigapwTime.refine == 5 ...
    & TigapwTime.dof == 1473;
assert(nnz(removePoint) == 1, ...
    'Expected one IGA-PW point with split (317,1156).');
TigapwTime = TigapwTime(~removePoint, :);

Ttime = [Ttime(~isIgapwTime, :); TigapwTime];
TigapwDof = TigapwTime;
TigapwDof.figure(:) = "dof";
Tdof = [Tdof(~isIgapwDof, :); TigapwDof];
end

function plot_method_compare_local(Tdof, Ttime, figDir)
% Plot method compare.
cfg = default_style_local();
methods = {'PW', 'IGA (p=1)', 'IGA (p=2)', 'IGA-PW'};
labels = {'Plane wave', 'IGA, $p=1$', 'IGA, $p=2$', 'IGA-PW'};
colors = cfg.lineColors([1 2 3 4], :);
markers = cfg.markers(1:numel(methods));

plot_one_compare_local(Tdof, methods, labels, colors, markers, cfg, ...
    'dof', 'Degrees of Freedom', ...
    fullfile(figDir, 'dof'));
plot_one_compare_local(Ttime, methods, labels, colors, markers, cfg, ...
    'total_time_s', 'Total Time (s)', ...
    fullfile(figDir, 'time'));
end

function plot_one_compare_local(T, methods, labels, colors, markers, cfg, xField, xLabelText, baseName)
% Plot one method-comparison panel against the selected x variable.
errLabel = '$|\lambda_{1}-\lambda_{1}^{\mathrm{DG}}|$';
Tigapw = sortrows(T(T.method == "IGA-PW", :), xField);
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
apply_log_limits_local(ax, T.(xField), T.err1, ...
    cfg.padXLow, cfg.padXHigh, cfg.padYLow, cfg.padYHigh);
annotate_igapw_dofs_local(ax, Tigapw, xField, cfg);
lgd = legend(ax, ...
    'Location', 'southwest', ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize);
lgd.Box = 'off';
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
close(fig);
end

function annotate_igapw_dofs_local(ax, T, xField, cfg)
% Annotate each IGA-PW point as (PW DOFs, IGA DOFs).
assert(height(T) == 5, 'Expected five IGA-PW points for annotation.');
[dx, dy, hAlign, vAlign] = annotation_layout_local();
xLim = ax.XLim;
yLim = ax.YLim;
xDen = log10(xLim(2)) - log10(xLim(1));
yDen = log10(yLim(2)) - log10(yLim(1));

for i = 1:height(T)
    nPw = count_pw_dofs_local(T.Nc(i));
    nIga = round(T.dof(i) - nPw);
    assert(nIga > 0, 'Invalid IGA/PW DOF split at row %d.', i);

    xNorm = (log10(T.(xField)(i)) - log10(xLim(1))) / xDen;
    yNorm = (log10(T.err1(i)) - log10(yLim(1))) / yDen;
    labelText = sprintf('$(%d,%d)$', nPw, nIga);
    text(ax, xNorm + dx(i), yNorm + dy(i), labelText, ...
        'Units', 'normalized', ...
        'Interpreter', 'latex', ...
        'FontSize', cfg.annotation.fontSize, ...
        'FontWeight', 'normal', ...
        'Color', cfg.annotation.color, ...
        'HorizontalAlignment', hAlign{i}, ...
        'VerticalAlignment', vAlign{i}, ...
        'Clipping', 'on');
end
end

function [dx, dy, hAlign, vAlign] = annotation_layout_local()
% Return tuned offsets that keep labels away from the IGA-PW line.
dx = [-0.002, -0.028, -0.028, -0.028, -0.028];
dy = [0.020, -0.008, -0.008, -0.008, -0.008];
hAlign = {'right', 'right', 'right', 'right', 'right'};
vAlign = {'bottom', 'middle', 'middle', 'middle', 'middle'};
end

function nPw = count_pw_dofs_local(Nc)
% Count integer wave vectors inside the radius-Nc disk.
N = floor(Nc);
nPw = 0;
for k1 = -N:N
    m = floor(sqrt(N ^ 2 - k1 ^ 2));
    nPw = nPw + 2 * m + 1;
end
end

function style_axes_local(ax, cfg)
% Apply axis styling.
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
% Apply axes layout settings.
ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
end

function apply_log_limits_local(ax, x, y, padXLow, padXHigh, padYLow, padYHigh)
% Apply log limits.
x = x(isfinite(x) & x > 0);
y = y(isfinite(y) & y > 0);
set(ax, 'XLim', [min(x) / (1 + padXLow), max(x) * (1 + padXHigh)]);
set(ax, 'YLim', [min(y) / (1 + padYLow), max(y) * (1 + padYHigh)]);
end

function cfg = default_style_local()
% Return plotting style values.
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
cfg.annotation.fontSize = 10;
cfg.annotation.color = [0.16 0.16 0.16];
cfg.padXLow = 2.0;
cfg.padXHigh = 0.4;
cfg.padYLow = 2;
cfg.padYHigh = 1;
end

function val = get_opt_local(s, fieldName, defaultVal)
% Read one option with a default value.
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function ensure_dir_local(p)
% Create a directory when it is missing.
if ~exist(p, 'dir')
    mkdir(p);
end
end
