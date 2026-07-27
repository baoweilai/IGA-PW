function figs = plot_helium_fields()
% Plot density and potential fields.

% Resolve workflow paths and paper styling.
rootDir = fileparts(fileparts(mfilename('fullpath')));
add_workflow_paths(fullfile(rootDir, 'model', 'helium_fields'), ...
    {'config', 'core', 'operators', 'solver'});
resultRoot = fullfile(rootDir, 'data', 'result');
plotRoot = fullfile(resultRoot, 'plots');
cfg = paper_style();
apply_paper_plot_style();
cfg.fig.width = 4.8; cfg.fig.height = 3.0;
cfg.fig.renderer = 'painters'; cfg.fig.bgColor = 'w';
cfg.layout.left = 0.14; cfg.layout.right = 0.04;
cfg.layout.bottom = 0.16; cfg.layout.top = 0.08;
cfg.axes.fontSize = 10; cfg.axes.lineWidth = 1.0;
cfg.axes.tickDir = 'out'; cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off'; cfg.axes.labelSize = 12;
cfg.legend.fontSize = 11;

% Load the reference and IGA-PW field data.
refRunFile = fullfile(resultRoot, 'REFERENCE', 'K_30', 'p_2', ...
    'nelem_32', 'run.mat');
cmpRunFile = fullfile(resultRoot, 'IGA', 'K_20', 'p_2', ...
    'nelem_08', 'run.mat');
if ~exist(refRunFile, 'file')
    error('Reference run file not found: %s. Run run_reference first.', refRunFile);
end
if ~exist(cmpRunFile, 'file')
    error('Comparison run file not found: %s. Run run_iga(2,8,20) first.', cmpRunFile);
end
zoomRadius = 0.05;
zoomInsetYPadRatio = 0.05;  % Pad the joint Reference/IGA-PW range in each inset.
zoomInsetYWindowRatio = {[], [], []};  % Relative inset y-window widths.
zoomInsetYLim = {[], [], []};  % Explicit inset y-limits.
zoomInsetYTicks = {[], [], []};
zoomInsetZeroGapFraction = {[], [], 1/3};
zoomInsetPosition = [ ...
    0.20, 0.58, 0.36, 0.31; ...
    0.20, 0.58, 0.36, 0.31; ...
    0.13, 0.13, 0.36, 0.31];  % [left bottom width height].
zoomInsetTickDecimals = [3, 2, 4];
zoomInsetTickExponent = {[], [], []};
mainProfileTickExponent = {[], [], []};
mainProfileYTickDecimals = [2, 2, 2];

ensure_dir_local(plotRoot);
refRun = load_run_local(refRunFile);
cmpRun = load_run_local(cmpRunFile);

[xRef, Lref, mRef] = build_grid_local(refRun.run);
[xCmp, Lcmp, mCmp] = build_grid_local(cmpRun.run);
validate_shared_grid_local(xRef, xCmp, Lref, Lcmp, mRef, mCmp);
xmid = xCmp(:);
domainRadius = Lcmp / 2;
orbitalDenseInnerN = 401;

% Define field styles and output containers.
fieldNames = {'orbital', 'density', 'hartree'};
fieldTags = {'orbital', 'density', 'hartree'};
refLabel = 'Reference';
cmpLabel = 'IGA-PW';
fieldStyles = build_field_styles_local();
profileStyle = example1_profile_style_local(cfg);
zoomOpt = profile_zoom_options_local(zoomRadius, zoomInsetYPadRatio, ...
    zoomInsetPosition(1, :), zoomInsetTickDecimals);
sliceTicks = {[0.10, 0.18, 0.26, 0.34], [0, 0.05, 0.10, 0.15, 0.20], []};
sliceTickLabels = {{'0.10', '0.18', '0.26', '0.34'}, ...
    {'0', '0.05', '0.10', '0.15', '0.20'}, {}};

figs = struct();
figs.profile = gobjects(numel(fieldNames), 1);
figs.slice = gobjects(numel(fieldNames), 1);
figs.profilePdf = cell(numel(fieldNames), 1);
figs.slicePdf = cell(numel(fieldNames), 1);

% Draw and export each radial profile and central slice.
for iField = 1:numel(fieldNames)
    refGrid = get_field_grid_local(refRun.run, fieldNames{iField});
    cmpGrid = get_field_grid_local(cmpRun.run, fieldNames{iField});

    refLine = extract_axis_line_local(refGrid, xmid, 'x');
    cmpLine = extract_axis_line_local(cmpGrid, xmid, 'x');
    [r, refRadial, cmpRadial] = build_radial_profile_local(xmid, refLine, cmpLine);
    [r, refRadial, cmpRadial] = replace_inner_profile_local( ...
        refRun.run, cmpRun.run, fieldNames{iField}, ...
        r, refRadial, cmpRadial, orbitalDenseInnerN);

    cmpSlice = extract_z0_slice_local(cmpGrid, xmid);
    [sliceLimits, cbTicks, cbTickLabels] = slice_colorbar_spec_local( ...
        cmpSlice, fieldStyles(iField).slice_prct, sliceTicks{iField}, ...
        sliceTickLabels{iField});

    figs.profile(iField) = figure('Color', cfg.fig.bgColor, 'Units', 'inches', ...
        'Position', profileStyle.figPos, 'Renderer', cfg.fig.renderer);
    set(figs.profile(iField), 'DefaultAxesFontName', profileStyle.fontName);
    set(figs.profile(iField), 'DefaultTextFontName', profileStyle.fontName);
    ax1 = axes('Parent', figs.profile(iField));
    set_tight_profile_layout_local(ax1, profileStyle);
    hProfile = plot_profile_pair_local(ax1, r, real(refRadial), real(cmpRadial), profileStyle, ...
        refLabel, cmpLabel);
    grid(ax1, 'off');
    xlim(ax1, [0, domainRadius]);
    apply_profile_y_padding_local(ax1, [real(refRadial(:)); real(cmpRadial(:))], ...
        profileStyle.padY);
    apply_main_profile_tick_labels_local(ax1, mainProfileTickExponent{iField});
    set_main_profile_five_ticks_local(ax1, mainProfileYTickDecimals(iField));
    xlabel(ax1, '');
    ylabel(ax1, '');
    title(ax1, '');
    hide_toolbar_local(ax1);
    zoomOptField = zoomOpt;
    zoomOptField.position = zoomInsetPosition(iField, :);
    zoomOptField.tickExponent = zoomInsetTickExponent{iField};
    zoomOptField.tickDecimals = zoomInsetTickDecimals(iField);
    zoomOptField.yWindowRatio = zoomInsetYWindowRatio{iField};
    zoomOptField.yLim = zoomInsetYLim{iField};
    zoomOptField.yTicks = zoomInsetYTicks{iField};
    zoomOptField.zeroGapFraction = zoomInsetZeroGapFraction{iField};
    draw_profile_zoom_inset_local(ax1, r, real(refRadial), real(cmpRadial), ...
        zoomOptField, profileStyle);
    lgd = legend(ax1, hProfile([2 1]), {cmpLabel, refLabel}, ...
        'Location', 'northeast');
    lgd.Box = 'off';
    lgd.FontName = profileStyle.fontName;
    lgd.FontSize = cfg.legend.fontSize;
    lgd.Interpreter = 'latex';
    figs.profilePdf{iField} = export_named_figure_local( ...
        figs.profile(iField), plotRoot, ['helium_profile_' fieldTags{iField}], ...
        'vector');

    figs.slice(iField) = figure('Color', cfg.fig.bgColor, 'Units', 'inches', ...
        'Position', profileStyle.sliceFigPos, 'Renderer', cfg.fig.renderer);
    set(figs.slice(iField), 'DefaultAxesFontName', cfg.fontName);
    set(figs.slice(iField), 'DefaultTextFontName', cfg.fontName);
    ax3 = axes('Parent', figs.slice(iField));
    imagesc(ax3, xmid, xmid, real(cmpSlice).');
    set(ax3, 'YDir', 'normal');
    axis(ax3, 'equal');
    axis(ax3, 'tight');
    xlim(ax3, [-domainRadius, domainRadius]);
    ylim(ax3, [-domainRadius, domainRadius]);
    xlabel(ax3, '');
    ylabel(ax3, '');
    title(ax3, '');
    colormap(ax3, cfg.errorMap);
    clim(ax3, sliceLimits);
    cb3 = colorbar(ax3);
    cb3.Label.String = '';
    cb3.FontSize = cfg.legend.fontSize;
    set_colorbar_ticks_local(cb3, cbTicks, cbTickLabels);
    style_line_axes_local(ax3, cfg);
    set_tight_slice_layout_local(ax3, cb3, profileStyle);
    set(ax3, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    axis(ax3, 'off');
    grid(ax3, 'off');
    hide_toolbar_local(ax3);
    figs.slicePdf{iField} = export_named_figure_local( ...
        figs.slice(iField), plotRoot, ['helium_slice_' fieldTags{iField}], ...
        'image');
end

end

function out = load_run_local(runFile)
% Load one saved run file.
if ~exist(runFile, 'file')
    error('Run file not found: %s', runFile);
end
out = load(runFile, 'run');
if ~isfield(out, 'run')
    error('run is missing in %s', runFile);
end
end

function [xmid, L, mFFT] = build_grid_local(run)
% Build the plotting grid.
if ~isfield(run, 'meta') || ~isfield(run.meta, 'L') || ~isfield(run.meta, 'grid_mFFT')
    error('run.meta.L or run.meta.grid_mFFT is missing.');
end

L = run.meta.L;
mFFT = run.meta.grid_mFFT;
edges = linspace(-L / 2, L / 2, mFFT + 1);
xmid = 0.5 * (edges(1:end-1) + edges(2:end));
end

function validate_shared_grid_local(xRef, xCmp, Lref, Lcmp, mRef, mCmp)
% Validate that field data share one grid.
if mRef ~= mCmp || abs(Lref - Lcmp) > 1e-12 || numel(xRef) ~= numel(xCmp) || ...
        any(abs(xRef(:) - xCmp(:)) > 1e-12)
    error('Reference and comparison runs do not share the same Cartesian grid.');
end
end

function fieldGrid = get_field_grid_local(run, fieldName)
% Return one field grid from saved data.
switch lower(strtrim(fieldName))
    case 'orbital'
        if ~isfield(run, 'meta') || ~isfield(run.meta, 'uGrid') || isempty(run.meta.uGrid)
            error('run.meta.uGrid is missing.');
        end
        fieldGrid = abs(run.meta.uGrid);

    case 'density'
        if ~isfield(run, 'rhoGrid') || isempty(run.rhoGrid)
            error('run.rhoGrid is missing.');
        end
        fieldGrid = real(run.rhoGrid);

    case 'hartree'
        if ~isfield(run, 'rhoGrid') || isempty(run.rhoGrid)
            error('run.rhoGrid is missing; cannot recompute corrected Hartree potential.');
        end
        if ~isfield(run, 'meta') || ~isfield(run.meta, 'L')
            error('run.meta.L is missing; cannot recompute corrected Hartree potential.');
        end
        fieldGrid = solve_poisson_fft_zero_mode_local(real(run.rhoGrid), run.meta.L);

    otherwise
        error('Unsupported field name: %s', fieldName);
end
end

function VH = solve_poisson_fft_zero_mode_local(rho, L)
% Solve the zero-mode Poisson equation by FFT.
m = size(rho, 1);
q = ifftshift(-floor(m / 2):ceil(m / 2) - 1);
[Qx, Qy, Qz] = ndgrid(q, q, q);
G2 = ((2 * pi / L) ^ 2) * (Qx .^ 2 + Qy .^ 2 + Qz .^ 2);

rhoHat = fftn(rho);
VHat = zeros(size(rhoHat));
mask = (G2 > 0);
VHat(mask) = 4 * pi * rhoHat(mask) ./ G2(mask);
VHat(~mask) = 0;

VH = ifftn(VHat, 'symmetric');
end

function lineData = extract_axis_line_local(fieldGrid, xmid, axisName)
% Extract a line through the selected axis.
[i0, i1, w0, w1] = interpolation_weights_for_point_local(xmid, 0);

switch lower(axisName)
    case 'x'
        lineData = w0 * w0 * fieldGrid(:, i0, i0) + ...
            w0 * w1 * fieldGrid(:, i0, i1) + ...
            w1 * w0 * fieldGrid(:, i1, i0) + ...
            w1 * w1 * fieldGrid(:, i1, i1);

    case 'y'
        lineData = w0 * w0 * squeeze(fieldGrid(i0, :, i0)) + ...
            w0 * w1 * squeeze(fieldGrid(i0, :, i1)) + ...
            w1 * w0 * squeeze(fieldGrid(i1, :, i0)) + ...
            w1 * w1 * squeeze(fieldGrid(i1, :, i1));

    case 'z'
        lineData = w0 * w0 * squeeze(fieldGrid(i0, i0, :)) + ...
            w0 * w1 * squeeze(fieldGrid(i0, i1, :)) + ...
            w1 * w0 * squeeze(fieldGrid(i1, i0, :)) + ...
            w1 * w1 * squeeze(fieldGrid(i1, i1, :));

    otherwise
        error('Unsupported axis name: %s', axisName);
end

lineData = lineData(:);
end

function sliceData = extract_z0_slice_local(fieldGrid, xmid)
% Extract the z = 0 slice.
[i0, i1, w0, w1] = interpolation_weights_for_point_local(xmid, 0);
sliceData = w0 * fieldGrid(:, :, i0) + w1 * fieldGrid(:, :, i1);
end

function [r, refRadial, cmpRadial] = build_radial_profile_local(xmid, refLine, cmpLine)
% Build a radial profile from grid data.
posMask = xmid > 0;
negMask = xmid < 0;
[i0, i1, w0, w1] = interpolation_weights_for_point_local(xmid, 0);

r = xmid(posMask);
refPos = refLine(posMask);
cmpPos = cmpLine(posMask);
refNeg = flipud(refLine(negMask));
cmpNeg = flipud(cmpLine(negMask));

n = min([numel(r), numel(refNeg), numel(cmpNeg)]);
refZero = w0 * refLine(i0) + w1 * refLine(i1);
cmpZero = w0 * cmpLine(i0) + w1 * cmpLine(i1);
r = [0; r(1:n)];
refRadial = [refZero; 0.5 * (refPos(1:n) + refNeg(1:n))];
cmpRadial = [cmpZero; 0.5 * (cmpPos(1:n) + cmpNeg(1:n))];
end

function [r, refRadial, cmpRadial] = replace_inner_profile_local( ...
    refRun, cmpRun, fieldName, rGrid, refGridRadial, cmpGridRadial, nDense)
% Dense IGA/FFT evaluation produces smooth profiles inside the IGA radius.
innerRadius = min(get_inner_radius_local(refRun), get_inner_radius_local(cmpRun));
if ~isfinite(innerRadius) || innerRadius <= 0 || nDense < 2
    r = rGrid;
    refRadial = refGridRadial;
    cmpRadial = cmpGridRadial;
    return;
end

rDense = linspace(0, innerRadius, nDense).';
refDense = eval_dense_profile_radial_local(refRun, fieldName, rDense);
cmpDense = eval_dense_profile_radial_local(cmpRun, fieldName, rDense);

outerMask = rGrid > innerRadius + 10 * eps(max(1, innerRadius));
r = [rDense; rGrid(outerMask)];
refRadial = [refDense; refGridRadial(outerMask)];
cmpRadial = [cmpDense; cmpGridRadial(outerMask)];
end

function values = eval_dense_profile_radial_local(run, fieldName, r)
% Evaluate one radial profile from saved coefficients rather than coarse grid samples.
switch lower(strtrim(fieldName))
    case 'orbital'
        values = eval_inner_orbital_radial_local(run, r);
    case 'density'
        values = eval_inner_density_radial_local(run, r);
    case 'hartree'
        values = eval_hartree_radial_fft_local(run, r);
    otherwise
        error('Unsupported dense profile field: %s', fieldName);
end
end

function a = get_inner_radius_local(run)
% Return the IGA subdomain radius.
if ~isfield(run, 'meta') || ~isfield(run.meta, 'a')
    error('run.meta.a is missing; cannot evaluate the inner orbital profile.');
end
a = double(run.meta.a);
end

function orbital = eval_inner_orbital_radial_local(run, r)
% Evaluate the saved IGA orbital on the x-axis.
if ~isfield(run, 'uh') || isempty(run.uh)
    error('run.uh is missing; cannot evaluate the inner orbital profile.');
end
if ~isfield(run, 'n_dofs_nurbs') || isempty(run.n_dofs_nurbs)
    error('run.n_dofs_nurbs is missing; cannot split the IGA coefficients.');
end
if ~isfield(run, 'nurbs_refine') || isempty(run.nurbs_refine)
    error('run.nurbs_refine is missing; cannot evaluate the IGA orbital.');
end

r = r(:);
z = zeros(size(r));
nI = double(run.n_dofs_nurbs);
cI = run.uh(1:nI, 1);
a = get_inner_radius_local(run);
uPos = eval_iga_orbital_on_points_local(run.nurbs_refine, cI, r, z, z, a);
uNeg = eval_iga_orbital_on_points_local(run.nurbs_refine, cI, -r, z, z, a);
orbital = 0.5 * (abs(uPos(:)) + abs(uNeg(:)));
end

function density = eval_inner_density_radial_local(run, r)
% Evaluate the saved IGA density on the x-axis.
if ~isfield(run, 'uh') || isempty(run.uh)
    error('run.uh is missing; cannot evaluate the inner density profile.');
end
if ~isfield(run, 'n_dofs_nurbs') || isempty(run.n_dofs_nurbs)
    error('run.n_dofs_nurbs is missing; cannot split the IGA coefficients.');
end
if ~isfield(run, 'nurbs_refine') || isempty(run.nurbs_refine)
    error('run.nurbs_refine is missing; cannot evaluate the IGA density.');
end

r = r(:);
z = zeros(size(r));
nI = double(run.n_dofs_nurbs);
cI = run.uh(1:nI, 1);
a = get_inner_radius_local(run);
uPos = eval_iga_orbital_on_points_local(run.nurbs_refine, cI, r, z, z, a);
uNeg = eval_iga_orbital_on_points_local(run.nurbs_refine, cI, -r, z, z, a);
density = 0.5 * (2 * abs(uPos(:)) .^ 2 + 2 * abs(uNeg(:)) .^ 2);
end

function hartree = eval_hartree_radial_fft_local(run, r)
% Evaluate Hartree potential on the x-axis by Fourier interpolation.
if ~isfield(run, 'rhoGrid') || isempty(run.rhoGrid)
    error('run.rhoGrid is missing; cannot evaluate the Hartree profile.');
end
if ~isfield(run, 'meta') || ~isfield(run.meta, 'L')
    error('run.meta.L is missing; cannot evaluate the Hartree profile.');
end

r = r(:);
hartreePos = eval_hartree_fft_axis_points_local(real(run.rhoGrid), run.meta.L, r);
hartreeNeg = eval_hartree_fft_axis_points_local(real(run.rhoGrid), run.meta.L, -r);
hartree = 0.5 * (hartreePos(:) + hartreeNeg(:));
end

function values = eval_hartree_fft_axis_points_local(rho, L, xq)
% Evaluate the zero-mode-corrected Hartree potential along the x-axis.
m = size(rho, 1);
q = ifftshift(-floor(m / 2):ceil(m / 2) - 1);
[Qx, Qy, Qz] = ndgrid(q, q, q);
G2 = ((2 * pi / L) ^ 2) * (Qx .^ 2 + Qy .^ 2 + Qz .^ 2);

rhoHat = fftn(rho);
VHat = zeros(size(rhoHat));
mask = (G2 > 0);
VHat(mask) = 4 * pi * rhoHat(mask) ./ G2(mask);
VHat(~mask) = 0;

dx = L / m;
x0 = -L / 2 + dx / 2;
gridIndexAtZero = (0 - x0) / dx;
phaseZero = exp(2i * pi * gridIndexAtZero * q / m);
planeCoeff = squeeze(sum(sum(VHat .* reshape(phaseZero, 1, [], 1) .* ...
    reshape(phaseZero, 1, 1, []), 3), 2));

xIndex = (xq(:) - x0) / dx;
phaseX = exp(2i * pi * xIndex * q / m);
values = real(phaseX * planeCoeff(:) / m ^ 3);
end

function val = eval_iga_orbital_on_points_local(nurbs_refine, coeff, X, Y, Z, a)
% Evaluate the tensor-product B-spline orbital coefficients.
U = nurbs_refine.Ubar;
V = nurbs_refine.Vbar;
W = nurbs_refine.Wbar;
pu = double(nurbs_refine.pu);
pv = double(nurbs_refine.pv);
pw = double(nurbs_refine.pw);
m = double(nurbs_refine.m);
n = double(nurbs_refine.n);

X = X(:);
Y = Y(:);
Z = Z(:);
val = zeros(numel(X), 1);

for q = 1:numel(X)
    u = min(max((X(q) + a) / (2 * a), 0), 1);
    v = min(max((Y(q) + a) / (2 * a), 0), 1);
    w = min(max((Z(q) + a) / (2 * a), 0), 1);

    ispan = findspan(U, pu, u);
    jspan = findspan(V, pv, v);
    kspan = findspan(W, pw, w);

    Nu = bspbasisDers(U, pu, u, 1); Nu = Nu(1, :)';
    Nv = bspbasisDers(V, pv, v, 1); Nv = Nv(1, :)';
    Nw = bspbasisDers(W, pw, w, 1); Nw = Nw(1, :)';

    ii = ispan-pu:ispan;
    jj = jspan-pv:jspan;
    kk = kspan-pw:kspan;

    s = 0;
    for kz = 1:(pw + 1)
        for jy = 1:(pv + 1)
            for ix = 1:(pu + 1)
                gid = ii(ix) + (jj(jy) - 1) * m + (kk(kz) - 1) * m * n;
                s = s + coeff(gid) * Nu(ix) * Nv(jy) * Nw(kz);
            end
        end
    end
    val(q) = s;
end
end

function [idx0, idx1, w0, w1] = interpolation_weights_for_point_local(x, xq)
% Compute interpolation weights at one point.
tol = 10 * eps(max(1, max(abs(x))));
idxExact = find(abs(x - xq) <= tol, 1, 'first');

if ~isempty(idxExact)
    idx0 = idxExact;
    idx1 = idxExact;
    w0 = 1;
    w1 = 0;
    return;
end

idx1 = find(x > xq, 1, 'first');
idx0 = find(x < xq, 1, 'last');

if isempty(idx0)
    idx0 = idx1;
    w0 = 1;
    w1 = 0;
    return;
end

if isempty(idx1)
    idx1 = idx0;
    w0 = 1;
    w1 = 0;
    return;
end

w1 = (xq - x(idx0)) / (x(idx1) - x(idx0));
w0 = 1 - w1;
end

function limits = robust_colormap_limits_local(data, prct)
% Choose robust color limits.
data = real(data(:));
data = data(isfinite(data));
if isempty(data)
    limits = [0, 1];
    return;
end
if all(data >= 0)
    upper = prctile(data, prct(2));
    lower = prctile(data, prct(1));
    if lower <= 0.15 * upper
        lower = 0;
    end
    limits = [lower, upper];
else
    q = prctile(data, prct);
    m = max(abs(q));
    limits = [-m, m];
end
limits = sanitize_limits(limits);
end

function [limits, ticks, tickLabels] = slice_colorbar_spec_local(data, prct, fixedTicks, fixedLabels)
% Return colorbar settings for a slice plot.
if ~isempty(fixedTicks)
    ticks = fixedTicks(:).';
    limits = [min(ticks), max(ticks)];
    tickLabels = fixedLabels;
    return;
end

limits = robust_colormap_limits_local(data, prct);
ticks = linspace(limits(1), limits(2), 5);
tickLabels = arrayfun(@(x) sprintf('%.3g', x), ticks, 'UniformOutput', false);
end

function set_colorbar_ticks_local(cb, ticks, ~)
% Set colorbar tick labels.
cb.Ticks = ticks;
cb.TickLabels = arrayfun(@format_colorbar_tick_local, ticks, ...
    'UniformOutput', false);
cb.TickLabelInterpreter = 'latex';
end

function label = format_colorbar_tick_local(x)
% Format one colorbar tick label.
if abs(x) < 10 * eps(max(1, abs(x)))
    label = '0';
else
    label = sprintf('%.2f', x);
end
end

function style = example1_profile_style_local(cfg)
% Return profile-line styles.
style = struct();
style.figPos = [1, 1, cfg.fig.width, cfg.fig.height];
style.fontName = cfg.fontName;
style.fontSize = cfg.axes.fontSize;
style.labelSize = cfg.axes.labelSize;
style.legendSize = cfg.legend.fontSize;
style.lineColors = [ ...
    223, 122, 094;
    060, 064, 091;
    130, 178, 154;
    242, 204, 142] / 255;
style.markers = {'o', 's', '^', 'd', 'x', '+'};
style.lineWidth = 1.8;
style.markerSize = 8;
style.markerCount = 9;
style.axesLineWidth = cfg.axes.lineWidth;
style.tickDir = cfg.axes.tickDir;
style.padY = 0.05;
style.zoomColor = [0.66, 0.66, 0.66];
style.zoomLineWidth = 0.30;
style.insetLineWidth = 1.2;
style.insetMarkerSize = 0;
style.insetFontSize = 7.0;
style.insetAxesLineWidth = 0.9;
style.profileAxesPosition = [cfg.layout.left, cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];
style.sliceFigPos = [1, 1, 3.05, 2.8];
style.sliceAxesPosition = [0.015, 0.045, 0.835, 0.910];
style.sliceColorbarPosition = [0.890, 0.120, 0.045, 0.780];
if isfield(cfg, 'fontName') && ~isempty(cfg.fontName)
    style.fontName = cfg.fontName;
end
end

function zoomOpt = profile_zoom_options_local(zoomRadius, yPadRatio, insetPosition, tickDecimals)
% Return zoom-inset options.
zoomOpt = struct();
zoomOpt.radius = zoomRadius;
zoomOpt.yPadRatio = yPadRatio;
zoomOpt.yWindowRatio = [];
zoomOpt.yLim = [];
zoomOpt.yTickCount = 3;
zoomOpt.tickDecimals = tickDecimals;
zoomOpt.tickExponent = [];
zoomOpt.position = insetPosition;
zoomOpt.yTicks = [];
zoomOpt.zeroGapFraction = [];
end

function apply_profile_y_padding_local(ax, yData, padRatio)
% Pad profile-axis limits.
yData = yData(isfinite(yData));
if isempty(yData)
    return;
end
yMin = min(yData);
yMax = max(yData);
ySpan = yMax - yMin;
if ySpan <= eps(max(1, max(abs(yData))))
    ySpan = max(1, max(abs(yData))) * 1e-4;
end
ylim(ax, [yMin - padRatio * ySpan, yMax + padRatio * ySpan]);
end

function set_tight_profile_layout_local(ax, style)
% Apply the profile figure layout.
ax.Units = 'normalized';
ax.Position = style.profileAxesPosition;
ax.LooseInset = [0, 0, 0, 0];
end

function set_tight_slice_layout_local(ax, cb, style)
% Apply the slice figure layout.
ax.Units = 'normalized';
cb.Units = 'normalized';
ax.Position = style.sliceAxesPosition;
cb.Position = style.sliceColorbarPosition;
ax.LooseInset = [0, 0, 0, 0];
end

function apply_main_profile_tick_labels_local(ax, tickExponent)
% Set main profile tick labels.
if isempty(tickExponent)
    return;
end
ticks = (1:5) * 10^tickExponent;
ax.YTick = ticks;
ax.YTickLabel = arrayfun(@(x) sprintf('%d', x), 1:5, ...
    'UniformOutput', false);
yl = ylim(ax);
ylim(ax, [min(yl(1), ticks(1)), max(yl(2), ticks(end))]);
text(ax, 0.02, 0.98, sprintf('$\\times 10^{%d}$', tickExponent), ...
    'Units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', ax.FontName, 'FontSize', ax.FontSize, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'Color', 'k', 'Clipping', 'off');
end

function set_main_profile_five_ticks_local(ax, yTickDecimals)
% Set five ticks on the main profile axis.
ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 5);
ax.YTick = linspace(ax.YLim(1), ax.YLim(2), 5);
ax.XTickLabel = arrayfun(@(x) sprintf('%.2f', x), ax.XTick, ...
    'UniformOutput', false);
fmt = sprintf('%%.%df', yTickDecimals);
ax.YTickLabel = arrayfun(@(x) sprintf(fmt, x), ax.YTick, ...
    'UniformOutput', false);
end

function pdfFile = export_named_figure_local(fig, plotRoot, baseName, contentType)
% Export a figure with a fixed filename.
pdfFile = fullfile(plotRoot, [baseName '.pdf']);
drawnow;
set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
exportgraphics(fig, pdfFile, 'ContentType', contentType, 'Resolution', 600);
end

function ensure_dir_local(pathstr)
% Create a directory when it is missing.
if ~exist(pathstr, 'dir')
    mkdir(pathstr);
end
end

function hide_toolbar_local(ax)
% Hide interactive figure toolbars.
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end

function hLines = plot_profile_pair_local(ax, r, refRadial, cmpRadial, style, refLabel, cmpLabel)
% Plot two radial profiles together.
markerIdx = marker_indices_by_x_local(r, style.markerCount);
hLines = plot(ax, r, [refRadial(:), cmpRadial(:)], '-');
set(hLines(1), 'LineWidth', style.lineWidth, ...
    'Color', style.lineColors(2, :), 'Marker', 'o', ...
    'MarkerSize', style.markerSize, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.lineColors(2, :), ...
    'MarkerIndices', markerIdx, 'DisplayName', refLabel);
set(hLines(2), 'LineWidth', style.lineWidth, ...
    'Color', style.lineColors(1, :), 'Marker', 's', ...
    'MarkerSize', style.markerSize, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.lineColors(1, :), ...
    'MarkerIndices', markerIdx, 'DisplayName', cmpLabel);
style_profile_axes_local(ax, style);
end

function markerIdx = marker_indices_by_x_local(x, markerCount)
% Place profile markers at approximately uniform x-locations.
x = x(:);
valid = find(isfinite(x));
if isempty(valid)
    markerIdx = [];
    return;
end

xValid = x(valid);
xMin = min(xValid);
xMax = max(xValid);
if xMax <= xMin || markerCount <= 1
    markerIdx = valid(1);
    return;
end

targets = linspace(xMin, xMax, min(markerCount, numel(valid)));
markerIdx = zeros(size(targets));
for i = 1:numel(targets)
    [~, j] = min(abs(xValid - targets(i)));
    markerIdx(i) = valid(j);
end
markerIdx = unique(markerIdx, 'stable');
end

function draw_profile_zoom_inset_local(axMain, r, refRadial, cmpRadial, zoomOpt, style)
% Draw the profile zoom inset.
% Determine the inset data range and vertical limits.
mask = r >= 0 & r <= zoomOpt.radius;
if nnz(mask) < 2
    return;
end

yZoom = [refRadial(mask); cmpRadial(mask)];
yZoom = yZoom(isfinite(yZoom));
if isempty(yZoom)
    return;
end
ySpan = max(yZoom) - min(yZoom);
if ySpan <= eps(max(abs(yZoom)))
    ySpan = max(1, max(abs(yZoom))) * 1e-4;
end
yAtZero = [refRadial(1), cmpRadial(1)];
yAtZero = yAtZero(isfinite(yAtZero));
if isempty(yAtZero)
    yCenter = median(yZoom);
else
    yCenter = mean(yAtZero);
end
halfSpan = max(abs(yZoom - yCenter));
if halfSpan <= eps(max(abs(yZoom)))
    halfSpan = 0.5 * ySpan;
end
halfSpan = halfSpan * (1 + zoomOpt.yPadRatio);
yLim = [yCenter - halfSpan, yCenter + halfSpan];
if ~isempty(zoomOpt.yLim)
    yLim = zoomOpt.yLim;
elseif ~isempty(zoomOpt.zeroGapFraction) && zoomOpt.zeroGapFraction > 0 && ...
        numel(yAtZero) == 2
    zeroGap = abs(yAtZero(1) - yAtZero(2));
    if zeroGap > 0
        halfSpan = zeroGap / (2 * zoomOpt.zeroGapFraction);
        yLim = [yCenter - halfSpan, yCenter + halfSpan];
    end
elseif ~isempty(zoomOpt.yWindowRatio) && zoomOpt.yWindowRatio > 0 && ...
        zoomOpt.yWindowRatio < 1
    halfSpan = 0.5 * ySpan * zoomOpt.yWindowRatio;
    yLim = [yCenter - halfSpan, yCenter + halfSpan];
end
xInsetLim = inset_x_limits_local(r, refRadial, cmpRadial, yLim, zoomOpt.radius);

% Mark the enlarged region on the main axes.
mainYLim = ylim(axMain);
boxYLim = [max(yLim(1), mainYLim(1)), min(yLim(2), mainYLim(2))];
if boxYLim(1) >= boxYLim(2)
    boxYLim = yLim;
end
line(axMain, [0 xInsetLim(2) xInsetLim(2) 0 0], ...
    [boxYLim(1) boxYLim(1) boxYLim(2) boxYLim(2) boxYLim(1)], ...
    'Color', style.zoomColor, 'LineWidth', style.zoomLineWidth, ...
    'LineStyle', '-', 'HandleVisibility', 'off', 'Clipping', 'on');

% Place the inset relative to the main plotting area.
fig = ancestor(axMain, 'figure');
drawnow;
oldUnits = axMain.Units;
axMain.Units = 'normalized';
pos = axMain.Position;
axMain.Units = oldUnits;
insetRel = zoomOpt.position;
axInset = axes('Parent', fig, 'Position', ...
    [pos(1) + insetRel(1) * pos(3), pos(2) + insetRel(2) * pos(4), ...
    insetRel(3) * pos(3), insetRel(4) * pos(4)], 'Color', 'w');

% Draw both profiles and their sample markers.
markerIdx = find(r >= xInsetLim(1) & r <= xInsetLim(2));
markerIdx = markerIdx(unique(round(linspace(1, numel(markerIdx), min(6, numel(markerIdx))))));
hInset = plot(axInset, r, [refRadial(:), cmpRadial(:)], '-');
set(hInset(1), 'LineWidth', style.insetLineWidth, ...
    'Color', style.lineColors(2, :));
set(hInset(2), 'LineWidth', style.insetLineWidth, ...
    'Color', style.lineColors(1, :));
if style.insetMarkerSize > 0
    set(hInset(1), 'Marker', 'o', ...
        'MarkerSize', style.insetMarkerSize, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.lineColors(2, :), ...
        'MarkerIndices', markerIdx);
    set(hInset(2), 'Marker', 's', ...
        'MarkerSize', style.insetMarkerSize, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.lineColors(1, :), ...
        'MarkerIndices', markerIdx);
end
% Format the inset and connect it to the marked region.
if isempty(zoomOpt.yTicks)
    ticks = linspace(yLim(1), yLim(2), zoomOpt.yTickCount);
else
    ticks = zoomOpt.yTicks;
end
set(axInset, 'XLim', xInsetLim, 'YLim', yLim, ...
    'YTick', ticks, 'YTickLabel', format_zoom_ticks_local(ticks, zoomOpt), ...
    'FontName', style.fontName, 'FontSize', style.insetFontSize, ...
    'LineWidth', style.insetAxesLineWidth, 'TickDir', style.tickDir, ...
    'TickLabelInterpreter', 'latex', 'Box', 'on', ...
    'XMinorTick', 'off', 'YMinorTick', 'off');
grid(axInset, 'off');
axInset.XMinorGrid = 'off';
axInset.YMinorGrid = 'off';
title(axInset, '');
xlabel(axInset, '');
ylabel(axInset, '');
hide_toolbar_local(axInset);
draw_zoom_connectors_local(axMain, axInset, xInsetLim, boxYLim, style);
uistack(axInset, 'top');
end

function xLim = inset_x_limits_local(r, refRadial, cmpRadial, yLim, radius)
% Return x-limits for the inset.
idx = find(r >= 0 & r <= radius & isfinite(r) & ...
    isfinite(refRadial(:)) & isfinite(cmpRadial(:)));
if numel(idx) < 2
    xLim = [0, radius];
    return;
end
x = r(idx);
xMax = radius;
xMax = min(xMax, first_y_boundary_crossing_local(x, refRadial(idx), yLim));
xMax = min(xMax, first_y_boundary_crossing_local(x, cmpRadial(idx), yLim));
if ~isfinite(xMax) || xMax <= 0
    xMax = radius;
end
xLim = [0, xMax];
end

function xCross = first_y_boundary_crossing_local(x, y, yLim)
% Find the first crossing of a y-boundary.
xCross = x(end);
if all(y >= yLim(1) & y <= yLim(2))
    return;
end
for j = 2:numel(x)
    if y(j) < yLim(1) || y(j) > yLim(2)
        if y(j) < yLim(1)
            yBound = yLim(1);
        else
            yBound = yLim(2);
        end
        dy = y(j) - y(j - 1);
        if abs(dy) <= eps(max(1, max(abs(y))))
            xCross = x(j - 1);
        else
            alpha = (yBound - y(j - 1)) / dy;
            alpha = min(max(alpha, 0), 1);
            xCross = x(j - 1) + alpha * (x(j) - x(j - 1));
        end
        return;
    end
end
end

function labels = format_zoom_ticks_local(ticks, zoomOpt)
% Format zoom-inset ticks.
if ~isempty(zoomOpt.tickExponent)
    ticks = ticks ./ 10^zoomOpt.tickExponent;
end
fmt = sprintf('%%.%df', zoomOpt.tickDecimals);
labels = arrayfun(@(x) sprintf(fmt, x), ticks, 'UniformOutput', false);
end

function draw_zoom_connectors_local(axMain, axInset, xLim, yLim, style)
% Draw connectors for the zoom inset.
fig = ancestor(axMain, 'figure');
[xLeft, yTop] = data_to_fig_norm_local(axMain, xLim(1), yLim(2));
[xRight, ~] = data_to_fig_norm_local(axMain, xLim(2), yLim(2));

oldUnits = axInset.Units;
axInset.Units = 'normalized';
insetPos = axInset.Position;
axInset.Units = oldUnits;

annotation(fig, 'line', [insetPos(1), xLeft], [insetPos(2), yTop], ...
    'Color', style.zoomColor, 'LineStyle', '-', 'LineWidth', style.zoomLineWidth);
annotation(fig, 'line', [insetPos(1) + insetPos(3), xRight], ...
    [insetPos(2), yTop], 'Color', style.zoomColor, ...
    'LineStyle', '-', 'LineWidth', style.zoomLineWidth);
end

function [xf, yf] = data_to_fig_norm_local(ax, x, y)
% Map data coordinates to figure coordinates.
[xn, yn] = data_to_axes_norm_local(ax, x, y);
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;
ax.Units = oldUnits;
xf = pos(1) + xn * pos(3);
yf = pos(2) + yn * pos(4);
end

function [xn, yn] = data_to_axes_norm_local(ax, x, y)
% Map data coordinates to axes coordinates.
xLim = ax.XLim;
yLim = ax.YLim;
xn = (x - xLim(1)) / (xLim(2) - xLim(1));
yn = (y - yLim(1)) / (yLim(2) - yLim(1));
xn = min(max(xn, 0), 1);
yn = min(max(yn, 0), 1);
end

function style_profile_axes_local(ax, style)
% Apply profile-axis styling.
set(ax, 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', style.axesLineWidth, 'TickDir', style.tickDir, ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'TickLabelInterpreter', 'latex', 'Box', 'on', 'Layer', 'top');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
end

function style_line_axes_local(ax, cfg)
% Apply line-axis styling.
set(ax, 'FontName', cfg.fontName, 'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, 'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', cfg.axes.xMinorTick, 'YMinorTick', cfg.axes.yMinorTick, ...
    'TickLabelInterpreter', 'latex', 'Box', 'on', 'Layer', 'top');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
end

function styles = build_field_styles_local()
% Return line and colormap settings for field plots.
styles(1) = struct( ...
    'slice_prct', [1, 99]);

styles(2) = struct( ...
    'slice_prct', [1, 99.5]);

styles(3) = struct( ...
    'slice_prct', [0.5, 99.8]);
end
