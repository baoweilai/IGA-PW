function plot_method_fields(resultDir)
% Plot method fields.

assert(exist('resultDir', 'var') == 1, 'plot_method_fields requires resultDir.');

dataFile = fullfile(resultDir, 'fields.mat');
if ~exist(dataFile, 'file')
    error('Cannot find plotting data: %s', dataFile);
end

S = load(dataFile);
Xg = S.Xg;
Yg = S.Yg;
errorFields = S.errorFields;
assert(iscell(errorFields), 'errorFields must be stored as a cell array.');

methodFileLabels = {'PW', 'IGA', 'IGA-PW'};
methodFileNames = {'pw_error', 'iga_error', 'iga_pw_error'};
maxErrors = S.maxErrors;

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 11.0 3.7], ...
    'Visible', 'off');
tl = tiledlayout(fig, 1, numel(errorFields), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(errorFields)
    ax = nexttile(tl, i);
    render_error_surface_local(ax, Xg, Yg, errorFields{i}, maxErrors(i), methodFileLabels{i});
end
export_pdf_local(fig, fullfile(resultDir, 'error_fields'), 600);
close_figure_if_valid_local(fig);

for i = 1:numel(errorFields)
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.15 3.75], ...
        'Visible', 'off');
    ax = axes(fig, 'Position', [0.03 0.02 0.76 0.96]);
    render_error_surface_local(ax, Xg, Yg, errorFields{i}, maxErrors(i), methodFileLabels{i});
    tighten_single_error_layout_local(ax);
    export_pdf_local(fig, fullfile(resultDir, methodFileNames{i}), 600);
    close_figure_if_valid_local(fig);
end
end

function tighten_single_error_layout_local(ax)
% Set the single-error figure layout.
ax.Position = [0.165 0.070 0.765 0.885];
end

function render_error_surface_local(ax, Xg, Yg, E, cmax, methodLabel)
% Render the error surface.
stride = max(1, ceil(max(size(E)) / 180));
rows = 1:stride:size(E, 1);
cols = 1:stride:size(E, 2);
surf(ax, Xg(rows, cols), Yg(rows, cols), E(rows, cols), E(rows, cols), ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
view(ax, [-45 30]);
axis(ax, 'tight');
zlim(ax, sanitize_color_limits_local([0 cmax]));
colormap(ax, example1_reference_colormap_local(256));
clim(ax, sanitize_color_limits_local([0 cmax]));
set(ax, 'XTick', linspace(-2, 2, 5), 'YTick', linspace(-2, 2, 5));
set_axis_error_ticks_local(ax, cmax, methodLabel);
xlabel(ax, '');
ylabel(ax, '');
zlabel(ax, '');
set_axes3d_style_local(ax);
disable_axes_toolbar_local(ax);
end

function set_axis_error_ticks_local(ax, cmax, methodLabel)
% Set ticks for the error axis.
[ticks, tickLabels] = nice_error_ticks_local(cmax, methodLabel);
ax.ZTick = ticks;
ax.ZTickLabel = tickLabels;
ax.TickLabelInterpreter = 'latex';
end

function [ticks, tickLabels] = nice_error_ticks_local(cmax, methodLabel)
% Choose ticks for the error axis.
assert(isfinite(cmax) && cmax > 0, 'Error color limit must be positive.');
ticks = linspace(0, cmax, 4);
tickLabels = format_error_tick_labels_local(ticks, methodLabel);
end

function tickLabels = format_error_tick_labels_local(ticks, methodLabel)
% Format error tick labels.
tickLabels = cell(size(ticks));
if strcmpi(methodLabel, 'IGA-PW')
    fmt = '%.3f';
else
    fmt = '%.2f';
end
for k = 1:numel(ticks)
    if abs(ticks(k)) < 1e-12
        tickLabels{k} = '0';
    else
        tickLabels{k} = sprintf(fmt, ticks(k));
    end
end
end

function cmap = example1_reference_colormap_local(n)
% Build the shared reference colormap.
anchors255 = [
    55, 105, 105;
    140, 190, 170;
    248, 242, 232;
    236, 170, 145;
    196,  85,  60
    ];
cmap = colormap_from_anchors_local(anchors255, n);
end

function cmap = colormap_from_anchors_local(anchors255, n)
% Build a colormap from anchor colors.
A = anchors255 / 255;
x = linspace(0, 1, size(A, 1));
xi = linspace(0, 1, n);
cmap = zeros(n, 3);
for j = 1:3
    cmap(:, j) = interp1(x, A(:, j), xi, 'pchip');
end
cmap = min(max(cmap, 0), 1);
end

function lims = sanitize_color_limits_local(lims)
% Validate color limits.
lims = lims(:).';
assert(numel(lims) == 2 && all(isfinite(lims)) && lims(2) > lims(1), ...
    'Color limits must be finite and increasing.');
end

function set_axes3d_style_local(ax)
% Apply 3-D axes style settings.
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 0.5, ...
    'GridAlpha', 0.15, 'GridColor', [0.15 0.15 0.15], ...
    'XMinorTick', 'off', 'YMinorTick', 'off');
end

function disable_axes_toolbar_local(ax)
% Hide the axes toolbar.
disableDefaultInteractivity(ax);
if isprop(ax, 'Toolbar')
    ax.Toolbar.Visible = 'off';
end
end

function export_pdf_local(fig, baseName, dpi)
% Export a rasterized PDF.
drawnow;
exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'image', ...
    'Resolution', dpi, 'BackgroundColor', 'white');
end

function close_figure_if_valid_local(fig)
% Close the figure when it exists.
close(fig);
end
