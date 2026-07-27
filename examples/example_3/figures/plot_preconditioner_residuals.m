function plot_preconditioner_residuals()
% Plot Example 3 residual histories.

% Resolve data paths and initialize the figure style.
close all; clc;

exampleDir = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(exampleDir, 'data', 'preconditioner', 'residuals');
outDir = dataRoot;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cfg = default_style();
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

% Load residual histories for all displayed methods and degrees.
none = readtable(fullfile(dataRoot, 'p1_refine05', 'none', ...
    'solver_history_last_scf_refine_05_none.csv'), 'TextType', 'string');
jacobi = readtable(fullfile(dataRoot, 'p1_refine05', 'purediag', ...
    'solver_history_last_scf_refine_05_jacobi.csv'), 'TextType', 'string');
tb1 = readtable(fullfile(dataRoot, 'p1_refine05', 'interfaceblock', ...
    'solver_history_last_scf_refine_05_pitbdg.csv'), 'TextType', 'string');

[x2, r2] = load_tb_history(fullfile(dataRoot, 'p_02', 'interfaceblock', 'run.mat'));
[x3, r3] = load_tb_history(fullfile(dataRoot, 'p_03', 'interfaceblock', 'run.mat'));

% Create the axes and draw every residual curve.
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
ax = axes(fig, 'Position', [0.15, 0.18, 0.80, 0.75]);
hold(ax, 'on');

colors = cfg.line.colors;
extraColors = [244 201 125; 032 161 185] / 255;

plot(ax, none.solver_iter, none.solver_error, '-', ...
    'Color', colors(1, :), 'LineWidth', cfg.line.width);
plot(ax, jacobi.solver_iter, jacobi.solver_error, '-', ...
    'Color', colors(2, :), 'LineWidth', cfg.line.width);
plot(ax, tb1.solver_iter, tb1.solver_error, '-', ...
    'Color', colors(3, :), 'LineWidth', cfg.line.width);
plot(ax, x2, r2, '-', ...
    'Color', extraColors(1, :), 'LineWidth', cfg.line.width);
plot(ax, x3, r3, '-', ...
    'Color', extraColors(2, :), 'LineWidth', cfg.line.width);

% Apply logarithmic residual scaling and export the PDF.
set(ax, ...
    'XScale', 'linear', ...
    'YScale', 'log', ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'Box', 'on', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
grid(ax, 'off');
xlim(ax, [0, 1000]);
ylim(ax, [1e-5, 1e2]);
xlabel(ax, 'Iteration number', 'FontSize', cfg.axes.labelSize);
ylabel(ax, 'Euclidean norm of residuals', 'FontSize', cfg.axes.labelSize);
legend(ax, ...
    {'Unpreconditioned, $p=1$', 'Jacobi, $p=1$', ...
    'TB-DG, $p=1$', 'TB-DG, $p=2$', 'TB-DG, $p=3$'}, ...
    'Location', 'south', ...
    'Box', 'off', ...
    'FontSize', cfg.legend.fontSize, ...
    'Interpreter', 'latex');

export_figure(fig, fullfile(outDir, 'Residuals_Nc=20_refine=5'), cfg);
close(fig);
end

function [x, r] = load_tb_history(runMat)
% Read the last SCF solver history.

S = load(runMat, 'run');
assert(isfield(S, 'run'), 'Missing run structure in %s.', runMat);
hist = S.run.result.state_results(1).scf_last_solver_history;
x = (1:size(hist, 1)).';
r = hist(:, 6);
end

function export_figure(fig, baseName, cfg)
% Save the final figure.

set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 cfg.fig.width cfg.fig.height], ...
    'PaperSize', [cfg.fig.width cfg.fig.height], ...
    'InvertHardcopy', 'off');
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
end

function cfg = default_style()
% Return the fixed residual-figure style.

cfg = struct();
cfg.fig.width = 4.8;
cfg.fig.height = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor = 'w';
cfg.axes.fontSize = 10;
cfg.axes.lineWidth = 1.0;
cfg.axes.tickDir = 'out';
cfg.axes.labelSize = 12;
cfg.legend.fontSize = 11;
cfg.line.colors = [ ...
    223 122 094;
    060 064 091;
    130 178 154] / 255;
cfg.line.width = 1.8;
end
