function plot_h_convergence()
% Plot h-convergence data.

clc; close all;

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ---------------- user params ----------------
Example    = 'Example_3';
Nc_fixed   = 20;
p_list     = [1 2];
nElem_list = [4 8 16 32 64];
lambda_ref = 5.408249091018272;

%% ---------------- unified style params ----------------
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
cfg.legend.location = 'southeast';
cfg.legend.box      = 'off';

%% ---------------- plot params ----------------
lw = 1.8;
ms = 8;

colors  = [223 122 094;
    060 064 091] / 255;
markers = {'o','s'};

slope_color  = [033 158 188] / 255;
slope_order  = 2;
slope_factor = 0.7 * 10^(-1/4);

xpad = 1.2;
ypad = 1.6;

%% ---------------- paths ----------------
resultRoot = fullfile(pwd, 'result', Example, sprintf('Nc_%02d', Nc_fixed));
plotDir    = fullfile(resultRoot, 'plots_h');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

%% ---------------- read data ----------------
data = struct('p', {}, 'nElem', {}, 'h', {}, 'err', {});

for ip = 1:numel(p_list)
    pdeg = p_list(ip);
    csvFile = fullfile(resultRoot, sprintf('p_%d', pdeg), 'summary.csv');
    assert(exist(csvFile, 'file') == 2, 'Missing h-convergence summary: %s', csvFile);

    T = readtable(csvFile);
    vars = T.Properties.VariableNames;
    assert(all(ismember({'nElem', 'h', 'lambda1'}, vars)), ...
        'Summary file must contain nElem, h, and lambda1: %s', csvFile);

    idx = ismember(T.nElem, nElem_list);
    assert(any(idx), 'No requested mesh sizes found in %s.', csvFile);

    nElem = T.nElem(idx);
    h = T.h(idx);
    err = abs(T.lambda1(idx) - lambda_ref);

    assert(all(isfinite(nElem) & isfinite(h) & isfinite(err) & h > 0 & err > 0), ...
        'Invalid h-convergence data in %s.', csvFile);

    [h, id] = sort(h, 'descend');
    data(end+1).p     = pdeg; %#ok<AGROW>
    data(end).nElem   = nElem(id);
    data(end).h       = h;
    data(end).err     = err(id);

    fprintf('\n========================================\n');
    fprintf('Example 3, Nc = %d, p = %d\n', Nc_fixed, pdeg);
    fprintf('----------------------------------------\n');
    for k = 1:numel(h)-1
        ord = log(err(id(k+1))/err(id(k))) / log(h(k+1)/h(k));
        fprintf('nElem pair (%d -> %d): h = %.4e -> %.4e, order = %.6f\n', ...
            nElem(id(k)), nElem(id(k+1)), h(k), h(k+1), ord);
    end
end

assert(~isempty(data), 'No h-convergence data were loaded.');

%% ---------------- reference slope line ----------------
href = [data(end).h(1); data(end).h(end)];
eref = slope_factor * data(end).err(end) * (href / data(end).h(end)).^slope_order;

%% ---------------- plot ----------------
fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);

ax = axes(fig);
ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];

hold(ax, 'on');
hh = gobjects(1, numel(data) + 1);

for k = 1:numel(data)
    hh(k) = loglog(ax, data(k).h, data(k).err, 'LineWidth', lw);
    hh(k).LineStyle       = '-';
    hh(k).Color           = colors(k,:);
    hh(k).Marker          = markers{k};
    hh(k).MarkerSize      = ms;
    hh(k).MarkerFaceColor = 'w';
    hh(k).MarkerEdgeColor = colors(k,:);
    hh(k).LineWidth       = lw;
    hh(k).DisplayName     = sprintf('$p=%d$', data(k).p);
end

hh(end) = loglog(ax, href, eref, 'LineWidth', lw);
hh(end).LineStyle   = '--';
hh(end).Color       = slope_color;
hh(end).LineWidth   = lw;
hh(end).Marker      = 'none';
hh(end).DisplayName = sprintf('$\\mathrm{Slope}=%d$', slope_order);

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
xlabel(ax, '$h$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
ylabel(ax, '$|\lambda_1-\lambda_{1}^{\mathrm{DG}}|$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);

xAll = [href(:); vertcat(data.h)];
yAll = [eref(:); vertcat(data.err)];
ax.XLim = [min(xAll)/xpad, max(xAll)*xpad];
ax.YLim = [min(yAll)/ypad, max(max(yAll)*ypad, 1e-3)];

ax.XTick = [1e-2 1e-1];
ax.XTickLabel = {'$10^{-2}$','$10^{-1}$'};

% -------- y tick labels in exponential form --------
emin = floor(log10(ax.YLim(1)));
emax = ceil(log10(ax.YLim(2)));
emax = max(emax, -3);
tickExp = emin:emax;
ax.YTick = 10.^tickExp;
ax.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), tickExp, 'UniformOutput', false);

lgd = legend(ax, 'Location', cfg.legend.location, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize);
lgd.Box = cfg.legend.box;

%% ---------------- save ----------------
baseName = 'h';

exportgraphics(fig, fullfile(plotDir, [baseName, '.pdf']), 'ContentType', 'vector');
fprintf('[SAVED] %s\n', fullfile(plotDir, [baseName, '.pdf']));

end
