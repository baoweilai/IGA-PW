function out = plot_cutoff_convergence(userCfg)
% Plot cutoff-convergence data.
assert(exist('userCfg', 'var') == 1, 'plot_cutoff_convergence requires userCfg.');

% Resolve workflow paths and read the cutoff configuration.
clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
add_workflow_paths(fullfile(rootDir, 'model', 'cutoff_convergence'), ...
    {'config', 'core', 'operators', 'solver'});

apply_paper_plot_style();

cfgRun = default_config(userCfg);
H = example_helpers(cfgRun);
pw = H.cfg.pw;

finalDir = fullfile(rootDir, 'data', 'result', 'PW', 'convergence');
if ~exist(finalDir, 'dir'), mkdir(finalDir); end

% Load the reference eigenvalue and allocate case errors.
refRunFile = H.case_run_file('reference', H.cfg.reference.p, ...
    H.cfg.reference.Nc, H.cfg.reference.Nelement);
lambdaRef = read_lambda_local(refRunFile);

KList = reshape(pw.Nc_list, 1, []);
nCase = numel(KList);
lambda = zeros(nCase, 1);
eigErr = zeros(nCase, 1);
dgErr = zeros(nCase, 1);

dgOpt = struct('innerGridN', H.cfg.state_error_grid_n, ...
    'outerGridN', H.cfg.state_error_grid_n, ...
    'faceGridN', H.cfg.state_error_grid_n, ...
    'chunkSize', 1200, 'Csigma', 10);

% Evaluate eigenvalue and DG errors for every cutoff.
for i = 1:nCase
    runFile = H.case_run_file('pw', pw.fixed_p, KList(i), pw.fixed_Nelement);
    lambda(i) = read_lambda_local(runFile);
    eigErr(i) = abs(lambda(i) - lambdaRef);
    assert(eigErr(i) > 0, 'Eigenvalue error must be positive for %s.', runFile);

    cacheFile = fullfile(finalDir, 'dg_cache', ...
        sprintf('ref_K_%d_p_%d_nelem_%02d_case_K_%d_p_%d_nelem_%02d.mat', ...
        H.cfg.reference.Nc, H.cfg.reference.p, H.cfg.reference.Nelement, ...
        KList(i), pw.fixed_p, pw.fixed_Nelement));
    E = cached_dg_error(cacheFile, refRunFile, runFile, dgOpt);
    dgErr(i) = E.errDG;
    assert(dgErr(i) > 0, 'DG error must be positive for %s.', runFile);
end

[eigDecay, dgDecay, eigReduction, dgReduction] = adjacent_rates_local(KList, eigErr, dgErr);

% Build and save value and adjacent-rate tables.
values = table(KList(:), repmat(pw.fixed_p, nCase, 1), ...
    repmat(pw.fixed_Nelement, nCase, 1), lambda, ...
    repmat(lambdaRef, nCase, 1), eigErr, dgErr, ...
    'VariableNames', {'K','p','nElem','lambda', ...
    'referenceLambda','eigAbsError','dgError'});

orders = table(KList(1:end-1).', KList(2:end).', ...
    compose('%d->%d', KList(1:end-1).', KList(2:end).'), ...
    eigDecay, dgDecay, eigReduction, dgReduction, ...
    'VariableNames', {'K_from','K_to','interval','eigDecayRate', ...
    'dgDecayRate','eigReductionFactor','dgReductionFactor'});

writetable(values, fullfile(finalDir, 'values.csv'));
writetable(orders, fullfile(finalDir, 'orders.csv'));

% Draw the combined convergence plot and export its data.
fig = plot_combined_local(KList, eigErr, dgErr);

pdfFile = fullfile(finalDir, 'cutoff.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

save(fullfile(finalDir, 'data.mat'), 'KList', 'lambda', 'lambdaRef', ...
    'eigErr', 'dgErr', 'eigDecay', 'dgDecay', 'eigReduction', ...
    'dgReduction', '-v7.3');

out = struct('fig', fig, 'pdf', pdfFile, ...
    'values', values, 'orders', orders);

end

function fig = plot_combined_local(KList, eigErr, dgErr)
% Plot combined convergence curves.
cfg = style_local();
eigLabel = '$|\lambda_1-\lambda_{1}^{\mathrm{DG}}|$';
dgLabel = '$\|u_1-u_{1}^{\mathrm{DG}}\|_{\mathrm{DG}}$';

fig = make_figure_local(cfg);
ax = axes(fig);
lines = semilogy(ax, KList(:), [eigErr(:), dgErr(:)], '-');
box(ax, 'on');
style_axes_local(ax, cfg);
style_line_local(lines(1), cfg, cfg.lineColors(1, :), cfg.markers{1}, eigLabel);
style_line_local(lines(2), cfg, cfg.lineColors(2, :), cfg.markers{2}, dgLabel);
finish_axes_local(ax, KList, [eigErr(:); dgErr(:)], cfg);
xlabel(ax, '$K$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(ax, 'Error', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
legend(ax, 'show', 'Location', 'northeast', 'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize, 'Box', cfg.legend.box);
end

function fig = make_figure_local(cfg)
% Create a figure with the saved size.
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
set(fig, 'ToolBar', 'none', 'MenuBar', 'none');
end

function style_line_local(lineObj, cfg, color, marker, labelText)
% Apply line styling.
set(lineObj, ...
    'LineWidth', cfg.lw, ...
    'Marker', marker, ...
    'MarkerSize', cfg.ms, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color, ...
    'Color', color, ...
    'DisplayName', labelText);
end

function style_axes_local(ax, cfg)
% Apply axis styling.
ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
set(ax, ...
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

function finish_axes_local(ax, KList, err, cfg)
% Finalize axis labels and limits.
yMin = min(err(err > 0));
yMax = max(err);
set(ax, ...
    'XTick', KList, ...
    'XLim', [min(KList) - cfg.padX, max(KList) + cfg.padX], ...
    'YLim', [yMin / (1 + cfg.padYLow), yMax * (1 + cfg.padYHigh)]);
end

function cfg = style_local()
% Return style values for this plot.
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
cfg.markers = {'o','s','^','d','x','+'};
cfg.lw = 1.8;
cfg.ms = 8;
cfg.padX = 0.25;
cfg.padYLow = 2;
cfg.padYHigh = 15;
end

function [eigDecay, dgDecay, eigReduction, dgReduction] = adjacent_rates_local(KList, eigErr, dgErr)
% Compute rates between adjacent data points.
dK = diff(KList(:));
eigReduction = eigErr(1:end-1) ./ eigErr(2:end);
dgReduction = dgErr(1:end-1) ./ dgErr(2:end);
eigDecay = log(eigReduction) ./ dK;
dgDecay = log(dgReduction) ./ dK;
end

function lambda = read_lambda_local(runFile)
% Read eigenvalue data from a MAT file.
S = load(runFile, 'run');
lambda = real(S.run.lambda(1));
end
