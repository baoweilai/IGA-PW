function plot_scaled_errors()
%Plot scaled error curves.

clc; close all; format short e;

%% ---------------- global style ----------------
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ---------------- user params ----------------
p_list     = [2];
alpha_list = [1.0, 2.0];

savePNG = true;
savePDF = true;
pngDPI  = 600;

%% ====================== style params ======================
cfg = struct();

% 1) figure
cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

% 2) layout
cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

% 3) curves
cfg.curve.lw       = 1.8;
cfg.curve.ms       = 8;
cfg.curve.colL2    = [223 122 094] / 255;
cfg.curve.colDG    = [060 064 091] / 255;
cfg.curve.colSlope = [033 158 188] / 255;
cfg.curve.mkL2     = 'o';
cfg.curve.mkDG     = 's';
cfg.curve.lsData   = '-';
cfg.curve.lsSlope  = '--';

% 4) axes
cfg.axes.fontSize   = 11;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xScale     = 'log';
cfg.axes.yScale     = 'log';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.box        = 'on';
cfg.axes.labelSize  = 13;

% 5) reference slope line
cfg.slope.lw        = 1.8;
cfg.slope.factor_p1 = [0.18, 0.07, 0.08, 0.10];
cfg.slope.factor_p2 = [0.10, 0.05, 0.07, 0.07];

% 6) padding
cfg.pad.x = 1.2;
cfg.pad.y = 1.8;

% 7) manual y ticks
cfg.ytick      = 10.^[-3 -2 -1];
cfg.yticklabel = {'$10^{-3}$','$10^{-2}$','$10^{-1}$'};

% 8) legend
cfg.legend.loc      = 'southeast';
cfg.legend.box      = 'off';
cfg.legend.fontSize = 11;

%% ---------------- paths ----------------
resultRoot = fullfile(pwd, 'result', 'scaled_errors');
outDir     = fullfile(resultRoot, 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- main loop ----------------
for ip = 1:numel(p_list)
    pdeg = p_list(ip);
    csvFile = fullfile(resultRoot, sprintf('p_%d', pdeg), 'summary.csv');

    assert(exist(csvFile, 'file') == 2, 'Missing scaled-error data file: %s', csvFile);

    T = readtable(csvFile);

    for ia = 1:numel(alpha_list)
        alpha = alpha_list(ia);
        gamma = (3 - alpha) / (2 * alpha);
        slope_ref = 1 + gamma;

        Ti = T(abs(T.alpha - alpha) < 1e-12, :);
        assert(~isempty(Ti), 'Missing scaled-error data for p=%d, alpha=%.1f.', pdeg, alpha);

        Ti = sortrows(Ti, 'h');
        x  = Ti.h(:);
        y1 = Ti.u1_L2(:);
        y2 = Ti.h_gamma_DG(:);

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

        set(ax, ...
            'FontSize',   cfg.axes.fontSize, ...
            'LineWidth',  cfg.axes.lineWidth, ...
            'TickDir',    cfg.axes.tickDir, ...
            'XScale',     cfg.axes.xScale, ...
            'YScale',     cfg.axes.yScale, ...
            'XMinorTick', cfg.axes.xMinorTick, ...
            'YMinorTick', cfg.axes.yMinorTick, ...
            'Box',        cfg.axes.box);

        grid(ax, 'off');
        ax.XMinorGrid = 'off';
        ax.YMinorGrid = 'off';
        ax.XRuler.MinorTick = 'off';
        ax.YRuler.MinorTick = 'off';

        h1 = plot(ax, x, y1, cfg.curve.lsData, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.colL2, ...
            'Marker', cfg.curve.mkL2, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.colL2);

        h2 = plot(ax, x, y2, cfg.curve.lsData, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.colDG, ...
            'Marker', cfg.curve.mkDG, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.colDG);

        % reference slope line
        x0 = x(round(numel(x)/2));
        if pdeg == 1
            fac = cfg.slope.factor_p1(ia);
        else
            fac = cfg.slope.factor_p2(ia);
        end
        yRef = fac * max([y1; y2]) * (x / x0).^slope_ref;

        h3 = plot(ax, x, yRef, cfg.curve.lsSlope, ...
            'Color', cfg.curve.colSlope, ...
            'LineWidth', cfg.slope.lw);

        xlabel(ax, '$h$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);
        ylabel(ax, '$\mathrm{Error}$', 'Interpreter', 'latex', 'FontSize', cfg.axes.labelSize);

        % x-axis expand
        xt = unique(x);
        xt(abs(xt - 0.03333333) < 1e-6 | abs(xt - 0.033) < 1e-6) = [];

        axisSpec.XTick = xt;
        axisSpec.XLim  = [min(x)/cfg.pad.x, max(x)*cfg.pad.x];

        % y-axis expand
        yAll = [y1; y2; yRef];
        yAll = yAll(isfinite(yAll) & yAll > 0);
        axisSpec.YLim = [min(yAll)/cfg.pad.y, max(yAll)*cfg.pad.y];

        % manual y-ticks within current y-limits
        yt = cfg.ytick(cfg.ytick >= axisSpec.YLim(1) & cfg.ytick <= axisSpec.YLim(2));
        if ~isempty(yt)
            axisSpec.YTick = yt;
            axisSpec.YTickLabel = cfg.yticklabel(ismember(cfg.ytick, yt));
        else
            axisSpec.YTick = [];
            axisSpec.YTickLabel = {};
        end

        ax.XTick = axisSpec.XTick;
        ax.XLim  = axisSpec.XLim;
        ax.YLim  = axisSpec.YLim;

        if ~isempty(axisSpec.YTick)
            ax.YTick = axisSpec.YTick;
            ax.YTickLabel = axisSpec.YTickLabel;
        end

        lgd = legend(ax, [h1, h2, h3], ...
            {'$\|u_1-u_1^{\mathrm{DG}}\|_{L^2}$', ...
            '$h^\gamma \|u_1-u_1^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
            sprintf('$\\mathrm{Slope} = %.2f$', slope_ref)}, ...
            'Interpreter', 'latex', ...
            'Location',    cfg.legend.loc, ...
            'FontSize',    cfg.legend.fontSize);
        lgd.Box = cfg.legend.box;

        baseFile = sprintf('scaled_%d', ia);

        pngFile = fullfile(outDir, [baseFile, '.png']);
        pdfFile = fullfile(outDir, [baseFile, '.pdf']);

        if savePNG
            exportgraphics(fig, pngFile, 'Resolution', pngDPI);
            fprintf('[SAVE] %s\n', pngFile);
        end
        if savePDF
            exportgraphics(fig, pdfFile, 'ContentType', 'vector');
            fprintf('[SAVE] %s\n', pdfFile);
        end
    end
end

end
