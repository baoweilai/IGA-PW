function plot_scaled_errors()
% Plot scaled error curves.

clc; close all; format short e;

%% ---------------- global style ----------------
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ---------------- user params ----------------
p_list     = 2;
alpha_list = [1.0, 2.0];

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
cfg.curve.lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;
cfg.curve.colSlope = [033 158 188] / 255;
cfg.curve.mkL2     = 'o';
cfg.curve.mkDG     = 's';
cfg.curve.mkL2u3   = '^';
cfg.curve.mkDGu3   = 'd';
cfg.curve.lsU1     = '-';
cfg.curve.lsU3     = '-';
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
cfg.slope.minGap    = 5;

% 6) padding
cfg.pad.x = 1.2;
cfg.ypad  = [1e-8, 1e-6];
cfg.ytop  = [1, 2];

% 7) manual y ticks for scale = 1 and 2
cfg.ytickExponent = {[-8, -6, -4, -2, 0], [-6, -4, -2, 0]};

% 8) two-column in-axes legend (Slope on the left)
cfg.legend.xL1      = 0.30;
cfg.legend.xL2      = 0.38;
cfg.legend.xLT      = 0.40;
cfg.legend.xR1      = 0.64;
cfg.legend.xR2      = 0.72;
cfg.legend.xRT      = 0.74;
cfg.legend.rowY     = [0.30, 0.18, 0.06];
cfg.legend.fontSize = 10;  % fixed at 10 pt

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
        y3 = Ti.u3_L2(:);
        y4 = Ti.u3_h_gamma_DG(:);

        fig = figure( ...
            'Color',    cfg.fig.bgColor, ...
            'Units',    'inches', ...
            'Position', [1 1 cfg.fig.width cfg.fig.height], ...
            'Renderer', cfg.fig.renderer);

        ax = axes('Parent', fig);
        ax.Toolbar.Visible = 'off';
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

        h1 = plot(ax, x, y1, cfg.curve.lsU1, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.lineColors(1, :), ...
            'Marker', cfg.curve.mkL2, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.lineColors(1, :));

        h2 = plot(ax, x, y2, cfg.curve.lsU1, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.lineColors(2, :), ...
            'Marker', cfg.curve.mkDG, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.lineColors(2, :));

        h3 = plot(ax, x, y3, cfg.curve.lsU3, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.lineColors(3, :), ...
            'Marker', cfg.curve.mkL2u3, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.lineColors(3, :));

        h4 = plot(ax, x, y4, cfg.curve.lsU3, ...
            'LineWidth', cfg.curve.lw, ...
            'Color', cfg.curve.lineColors(4, :), ...
            'Marker', cfg.curve.mkDGu3, ...
            'MarkerSize', cfg.curve.ms, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cfg.curve.lineColors(4, :));

        % reference slope line
        x0 = x(round(numel(x)/2));
        yEnvelope = max([y1, y2, y3, y4], [], 2);
        yScale = max(yEnvelope);
        slopeShape = (x / x0).^slope_ref;
        fac = cfg.slope.minGap * max(yEnvelope ./ (yScale * slopeShape));
        yRef = fac * yScale * slopeShape;

        h5 = plot(ax, x, yRef, cfg.curve.lsSlope, ...
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
        yAll = [y1; y2; y3; y4; yRef];
        yAll = yAll(isfinite(yAll) & yAll > 0);
        axisSpec.YLim = [cfg.ypad(ia), cfg.ytop(ia)];
        assert(min(yAll) >= axisSpec.YLim(1) && max(yAll) <= axisSpec.YLim(2), ...
            'Scaled-error curves exceed the requested y-limits for scale index %d.', ia);

        % manual y-ticks
        tickExponent = cfg.ytickExponent{ia};
        axisSpec.YTick = 10.^tickExponent;
        axisSpec.YTickLabel = arrayfun(@(e) sprintf('$10^{%d}$', e), ...
            tickExponent, 'UniformOutput', false);

        ax.XTick = axisSpec.XTick;
        ax.XLim  = axisSpec.XLim;
        ax.YLim  = axisSpec.YLim;
        ax.YTick = axisSpec.YTick;
        ax.YTickLabel = axisSpec.YTickLabel;

        draw_scaled_legend_in_axes(ax, [h1, h2, h3, h4], h5, slope_ref, cfg);

        baseFile = sprintf('scaled_%d', ia);

        pdfFile = fullfile(outDir, [baseFile, '.pdf']);

        exportgraphics(fig, pdfFile, 'ContentType', 'vector');
        fprintf('[SAVE] %s\n', pdfFile);
    end
end

function draw_scaled_legend_in_axes(ax, hError, hSlope, slopeRef, cfg)
% Draw a two-column legend inside the axes.

labels = { ...
    '$\|u_1-u_1^{\mathrm{DG}}\|_{L^2}$', ...
    '$h^\beta \|u_1-u_1^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
    '$\|u_3-u_3^{\mathrm{DG}}\|_{L^2}$', ...
    '$h^\beta \|u_3-u_3^{\mathrm{DG}}\|_{\mathrm{DG}}$'};

leftEntries = [1, 3];
rightEntries = [2, 4];

for j = 1:2
    idx = leftEntries(j);
    draw_one_fake_entry_in_axes(ax, hError(idx), '-', ...
        cfg.legend.xL1, cfg.legend.xL2, cfg.legend.xLT, ...
        cfg.legend.rowY(j + 1), labels{idx}, cfg.legend.fontSize);
end

for j = 1:2
    idx = rightEntries(j);
    draw_one_fake_entry_in_axes(ax, hError(idx), '-', ...
        cfg.legend.xR1, cfg.legend.xR2, cfg.legend.xRT, ...
        cfg.legend.rowY(j + 1), labels{idx}, cfg.legend.fontSize);
end

slopeLabel = sprintf('$\\mathrm{Slope} = %.2f$', slopeRef);
draw_one_fake_entry_in_axes(ax, hSlope, '--', ...
    cfg.legend.xR1, cfg.legend.xR2, cfg.legend.xRT, ...
    cfg.legend.rowY(1), slopeLabel, cfg.legend.fontSize);

end

function draw_one_fake_entry_in_axes(ax, hLine, legendLineStyle, x1n, x2n, xtn, yn, labelStr, fontSize)
% Draw one fake-legend entry using normalized axes coordinates.

c  = get(hLine, 'Color');
lw = get(hLine, 'LineWidth');

mk  = get(hLine, 'Marker');
ms  = get(hLine, 'MarkerSize');
mfc = get(hLine, 'MarkerFaceColor');
mec = get(hLine, 'MarkerEdgeColor');

[x1, y1] = axes_norm_to_data(ax, x1n, yn);
[x2, y2] = axes_norm_to_data(ax, x2n, yn);
[xt, yt] = axes_norm_to_data(ax, xtn, yn);

xm = sqrt(x1 * x2);
ym = y1;

line(ax, [x1 x2], [y1 y2], ...
    'LineStyle', legendLineStyle, ...
    'Color', c, ...
    'LineWidth', lw, ...
    'Marker', 'none', ...
    'Clipping', 'off');

if ~strcmp(mk, 'none')
    line(ax, xm, ym, ...
        'LineStyle', 'none', ...
        'Color', c, ...
        'Marker', mk, ...
        'MarkerSize', ms, ...
        'MarkerFaceColor', mfc, ...
        'MarkerEdgeColor', mec, ...
        'LineWidth', lw, ...
        'Clipping', 'off');
end

text(ax, xt, yt, labelStr, ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', ...
    'Margin', 0.5, ...
    'Clipping', 'off');

end

function [x, y] = axes_norm_to_data(ax, xn, yn)
% Convert normalized axes coordinates to data coordinates.

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

end
