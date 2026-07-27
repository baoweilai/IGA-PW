function plot_eigen_error()
% Plot eigen error.

clc; close all; format short e;

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');

prepare_example_2_plot_context(mfilename('fullpath'));

% Set current, reference, quadrature, and plot-range parameters.
Example     = 'Example_2';

Nc_fixed    = 30;
pdeg        = 2;
refine_list = 2:6;

ref_Nc      = 45;
ref_pdeg    = 3;
ref_refine  = 8;

u_list      = [1 2 3];

L           = 4;
a           = 0.2;

patchCenters = [-1, 0;
    1, 0];

Csigma      = 10;
dx_in       = 5e-3;
dx_out      = 5e-3;
chunkSize   = 20000;

slope1_factor = 0.05;
slope2_factor = 0.0003;

padX = 0.4;
padY = 0.4;

% Define figure, axes, and curve styles.
cfg = struct();

cfg.fig.width    = 5.039;
cfg.fig.height   = 3.0045;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

cfg.layout.left   = 0.18;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

cfg.axes.fontSize   = 11;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';

cfg.curve.lineWidth       = 1.8;
cfg.curve.markerSize      = 8;
cfg.curve.markerFaceColor = 'w';

cfg.label.fontSize = 12;

% Set fixed axis limits and ticks.
axisSpec.XTick      = [];
axisSpec.XTickLabel = {};
axisSpec.XLim       = [];

axisSpec.YTick      = 10.^[-6 -4 -2 0];
axisSpec.YTickLabel = {'$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$'};
axisSpec.YLim       = [1e-7, 1e0];

% Set manual in-axes legend positions.
legendSpec.xL1 = 0.40;
legendSpec.xL2 = 0.50;
legendSpec.xLT = 0.52;

legendSpec.xR1 = 0.70;
legendSpec.xR2 = 0.88;
legendSpec.xRT = 0.90;

legendSpec.rowY      = [0.4, 0.3, 0.2, 0.1];
legendSpec.slopeRows = [3, 4];

legendSpec.slope1Label = '$\mathrm{Slope} = 1$';
legendSpec.slope2Label = '$\mathrm{Slope} = 2$';
legendSpec.fontSize    = 10;

% Define curve and reference-slope colors.
lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154] / 255;

slope1Color = [239 065 067] / 255;
slope2Color = [033 158 188] / 255;

% Resolve the h-convergence workflow and output paths.
activate_example_workflow('h_convergence', ...
    {'nurbs', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});

ExampleRoot = fullfile(pwd, 'result', Example);
assert(isfolder(ExampleRoot), 'Missing Example 2 result directory: %s', ExampleRoot);

dxTag  = make_dx_tag(dx_in, dx_out);

outDir = fullfile(ExampleRoot, 'eigen_error', ...
    sprintf('Nc_%d', Nc_fixed), dxTag);

if ~exist(outDir,'dir')
    mkdir(outDir);
end

errCsvFile  = fullfile(outDir, 'errors.csv');
rateCsvFile = fullfile(outDir, 'rates.csv');

pdfFile = fullfile(outDir, 'eigen_error.pdf');

% Align current plane-wave indices with the reference basis.
[p_vec_curr, n_pw_curr] = generate_p_vec(Nc_fixed);
[p_vec_ref,  n_pw_ref ] = generate_p_vec(ref_Nc);

[tf_k, loc_k] = ismember(p_vec_curr, p_vec_ref, 'rows');

% Load or compute eigenfunction errors for the selected degree.
fprintf('\n==================== p = %d ====================\n', pdeg);

requiredVars = {'refine','dx', ...
    'u1_L2','u1_DG', ...
    'u2_L2','u2_DG', ...
    'u3_L2','u3_DG', ...
    'ref_Nc','ref_pdeg','ref_refine'};

needRecompute = true;
T = [];

if exist(errCsvFile, 'file')

    T_all = readtable(errCsvFile);

    if all(ismember(requiredVars, T_all.Properties.VariableNames))
        refMask = T_all.ref_Nc == ref_Nc ...
            & T_all.ref_pdeg == ref_pdeg ...
            & T_all.ref_refine == ref_refine;
        T = T_all(ismember(T_all.refine, refine_list) & refMask, :);

        refine_in_cache = unique(T.refine(:)).';
        refine_need     = unique(refine_list(:)).';

        if ~isempty(T) && all(ismember(refine_need, refine_in_cache))
            needRecompute = false;
        end
    end
end

if needRecompute

    ref_runMat = fullfile(ExampleRoot, sprintf('Nc_%d', ref_Nc), ...
        sprintf('p_%d', ref_pdeg), sprintf('refine_%02d', ref_refine), 'run.mat');

    if ~exist(ref_runMat, 'file')
        error('Missing reference run file: %s', ref_runMat);
    end

    S = load(ref_runMat, 'run');
    assert(isfield(S, 'run'), 'Reference MAT file does not contain run.');
    Rref = S.run;

    ref = extract_run_data_example2(Rref, n_pw_ref);

    fprintf('[REF ] Nc=%d p=%d refine=%d nNURBS=%d (n1=%d, n2=%d)\n', ...
        ref_Nc, ref_pdeg, ref_refine, ref.nNURBS, ref.n_dofs_1, ref.n_dofs_2);

    runs = load_runs_over_refine_example2(ExampleRoot, Nc_fixed, pdeg, refine_list, n_pw_curr);

    nR = numel(runs);

    dx_axis       = zeros(nR,1);
    eL2           = zeros(nR, numel(u_list));
    eDG           = zeros(nR, numel(u_list));
    refine_loaded = zeros(nR,1);

    for i = 1:nR
        rr = runs(i);
        refine_loaded(i) = rr.refine;

        h_param = estimate_h_parametric_example2(rr);
        h_phys  = 2 * a * h_param;
        dx_axis(i) = h_phys;

        sigma = Csigma * (Nc_fixed + 1/max(h_phys,1e-14));

        for ju = 1:numel(u_list)
            uk = u_list(ju);

            u_curr = rr.uh(:,uk);
            u_ref  = ref.uh(:,uk);

            curr_NURBS = u_curr(1:rr.nNURBS);
            curr_PW    = u_curr(rr.nNURBS+1 : rr.nNURBS+n_pw_curr);

            ref_NURBS  = u_ref(1:ref.nNURBS);
            ref_PW     = u_ref(ref.nNURBS+1 : ref.nNURBS+n_pw_ref);

            curr_PW_padded = zeros(size(ref_PW));
            curr_PW_padded(loc_k(tf_k)) = curr_PW(tf_k);

            alpha = ref_PW' * curr_PW_padded;
            if abs(alpha) > 1e-14
                phase_factor = exp(-1i * angle(alpha));
            else
                phase_factor = 1;
            end

            curr_NURBS = curr_NURBS * phase_factor;
            curr_PW    = curr_PW    * phase_factor;

            err_IGA_L2 = iga_L2_error_example2( ...
                rr, curr_NURBS, ...
                ref, ref_NURBS, ...
                patchCenters, a, dx_in);

            err_PW_L2 = pw_L2_error_example2( ...
                curr_PW, p_vec_curr, ...
                ref_PW, p_vec_ref, ...
                L, patchCenters, a, dx_out, chunkSize);

            eL2(i,ju) = sqrt(err_IGA_L2^2 + err_PW_L2^2);

            eDG(i,ju) = compute_DG_error_example2( ...
                rr, curr_NURBS, curr_PW, ...
                ref, ref_NURBS, ref_PW, ...
                p_vec_curr, p_vec_ref, L, patchCenters, a, ...
                dx_in, dx_out, chunkSize, sigma);

            fprintf('  refine=%02d  u_%d:  L2=%.3e  DG=%.3e\n', ...
                rr.refine, uk, eL2(i,ju), eDG(i,ju));
        end
    end

    [dx_axis, ord] = sort(dx_axis, 'ascend');
    refine_loaded  = refine_loaded(ord);
    eL2 = eL2(ord,:);
    eDG = eDG(ord,:);

    eL2(eL2 <= 0) = eps;
    eDG(eDG <= 0) = eps;

    T = table( ...
        refine_loaded, dx_axis, ...
        eL2(:,1), eDG(:,1), ...
        eL2(:,2), eDG(:,2), ...
        eL2(:,3), eDG(:,3), ...
        repmat(ref_Nc, nR, 1), ...
        repmat(ref_pdeg, nR, 1), ...
        repmat(ref_refine, nR, 1), ...
        'VariableNames', requiredVars);

    writetable(T, errCsvFile);
end

% Filter the requested refinements and assemble plotted errors.
T = T(ismember(T.refine, refine_list), :);
T = sortrows(T, 'refine');

if isempty(T)
    error('p=%d has missing or invalid run data.', pdeg);
end

Y = [T.u1_DG, T.u2_DG, T.u3_DG];

x_axis = T.dx;
Y(Y <= 0) = eps;

T_rate = build_local_rate_table(T, pdeg);
writetable(T_rate, rateCsvFile);

curveLabels = { ...
    '$i=1$', ...
    '$i=2$', ...
    '$i=3$'};

fig = figure('Color', cfg.fig.bgColor, ...
    'Units','inches', ...
    'Position',[1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);

ax = axes(fig);
hold(ax,'on');
box(ax,'on');

ax.Units = 'normalized';
ax.Position = [ ...
    cfg.layout.left, ...
    cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];

set(ax, ...
    'FontSize', cfg.axes.fontSize, ...
    'TickDir', cfg.axes.tickDir, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'XScale','log', ...
    'YScale','log', ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick);

grid(ax,'off');

mk = {'s','d','>'};
hCurve = gobjects(numel(curveLabels),1);

for k = 1:numel(curveLabels)
    col = lineColors(k,:);
    hCurve(k) = plot(ax, x_axis, Y(:,k), '-', ...
        'LineWidth', cfg.curve.lineWidth, ...
        'Color', col, ...
        'Marker', mk{k}, ...
        'MarkerSize', cfg.curve.markerSize, ...
        'MarkerFaceColor', cfg.curve.markerFaceColor, ...
        'MarkerEdgeColor', col);
end

xGrid = unique(x_axis(:));
anchor = choose_anchor_index(xGrid);
h0 = xGrid(anchor);
yMax = max(Y(:));

yref1 = (slope1_factor * yMax) * (xGrid ./ h0).^1;
hSlope1 = plot(ax, xGrid, yref1, '--', ...
    'Color', slope1Color, ...
    'LineWidth', 1.8);

yref2 = (slope2_factor * yMax) * (xGrid ./ h0).^2;
hSlope2 = plot(ax, xGrid, yref2, '--', ...
    'Color', slope2Color, ...
    'LineWidth', 1.8);

xlabel(ax, '$h$', ...
    'FontSize', cfg.label.fontSize, ...
    'Interpreter', 'latex');

ylabel(ax, '$\|u_i-u_{i}^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
    'FontSize', cfg.label.fontSize, ...
    'Interpreter', 'latex');

apply_manual_log_axes(ax, x_axis, Y, yref1, yref2, axisSpec, padX, padY);
draw_fake_legend_in_axes(ax, hCurve, hSlope1, hSlope2, curveLabels, legendSpec);

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('[SAVE] %s\n', pdfFile);

end

function ref = extract_run_data_example2(R, n_pw_basis)
% Extract the Example 2 reference fields from one saved run.

ref = struct();

ref.uh = R.uh;

ref.nurbs_original_1 = R.nurbs_original_1;
ref.nurbs_refine_1   = R.nurbs_refine_1;
ref.nurbs_original_2 = R.nurbs_original_2;
ref.nurbs_refine_2   = R.nurbs_refine_2;
requiredFields = {'n_dofs_1', 'n_dofs_2', 'n_dofs_nurbs', 'n_pw_basis'};
assert(all(isfield(R, requiredFields)), ...
    'Reference run is missing dimension metadata.');
ref.n_dofs_1 = double(R.n_dofs_1);
ref.n_dofs_2 = double(R.n_dofs_2);
ref.nNURBS = double(R.n_dofs_nurbs);
ref.n_pw_basis = double(R.n_pw_basis);
assert(ref.nNURBS == ref.n_dofs_1 + ref.n_dofs_2, ...
    'Reference NURBS dimensions are inconsistent.');
assert(ref.n_pw_basis == n_pw_basis, ...
    'Reference plane-wave dimension is inconsistent.');

end

function T_rate = build_local_rate_table(T, pdeg)
% Build local rate table.

metricNames = {'u1_L2','u1_DG','u2_L2','u2_DG','u3_L2','u3_DG'};
n = height(T);

if n < 2
    T_rate = cell2table(cell(0,8), ...
        'VariableNames', {'p','metric','refine_from','refine_to','dx_from','dx_to','err_to','rate'});
    return;
end

Rows = cell((n-1) * numel(metricNames), 8);
rc = 0;

for j = 1:numel(metricNames)
    metric = metricNames{j};
    e = T.(metric);

    for i = 1:n-1
        rc = rc + 1;
        rate = log(e(i+1)/e(i)) / log(T.dx(i+1)/T.dx(i));

        Rows(rc,:) = { ...
            pdeg, ...
            metric, ...
            T.refine(i), ...
            T.refine(i+1), ...
            T.dx(i), ...
            T.dx(i+1), ...
            round_sig(e(i+1), 10), ...
            round_sig(rate, 10)};
    end
end

T_rate = cell2table(Rows, ...
    'VariableNames', {'p','metric','refine_from','refine_to','dx_from','dx_to','err_to','rate'});

end

function apply_manual_log_axes(ax, x_axis, Y, yref1, yref2, spec, padX, padY)
% Apply manual log axes.

xAll = x_axis(:);
yAll = [Y(:); yref1(:); yref2(:)];
yAll = yAll(isfinite(yAll) & yAll > 0);

if isempty(spec.XLim)
    xmin = min(xAll);
    xmax = max(xAll);

    if xmin == xmax
        ax.XLim = [xmin/1.5, xmax*1.5];
    else
        ax.XLim = [xmin/(1+padX), xmax*(1+padX)];
    end
else
    ax.XLim = spec.XLim;
end

if ~isempty(spec.XTick)
    ax.XTick = spec.XTick;
end

if ~isempty(spec.XTickLabel)
    ax.XTickLabel = spec.XTickLabel;
end

if isempty(spec.YLim)
    ymin = min(yAll);
    ymax = max(yAll);

    if ymin == ymax
        ax.YLim = [ymin/10, ymax*10];
    else
        ax.YLim = [ymin/(1+padY), ymax*(1+padY)];
    end
else
    ax.YLim = spec.YLim;
end

if ~isempty(spec.YTick)
    ax.YTick = spec.YTick;
end

if ~isempty(spec.YTickLabel)
    ax.YTickLabel = spec.YTickLabel;
end

end

function idx = choose_anchor_index(xGrid)
% Choose the point used for slope annotation.

n = numel(xGrid);

if n <= 2
    idx = 1;
else
    idx = round(n/2);
end

end

function draw_fake_legend_in_axes(ax, hCurve, hSlope1, hSlope2, curveLabels, spec)
% Draw an in-axes legend using proxy curves.

for j = 1:numel(curveLabels)
    draw_one_fake_entry_in_axes(ax, hCurve(j), ...
        spec.xR1, spec.xR2, spec.xRT, spec.rowY(j+1), curveLabels{j}, spec.fontSize);
end

draw_one_fake_entry_in_axes(ax, hSlope1, ...
    spec.xL1, spec.xL2, spec.xLT, spec.rowY(spec.slopeRows(1)), spec.slope1Label, spec.fontSize);

draw_one_fake_entry_in_axes(ax, hSlope2, ...
    spec.xL1, spec.xL2, spec.xLT, spec.rowY(spec.slopeRows(2)), spec.slope2Label, spec.fontSize);

end

function draw_one_fake_entry_in_axes(ax, hLine, x1n, x2n, xtn, yn, labelStr, fontSize)
% Draw one manual legend entry inside the axes.

% Read the source line and marker styling.
c  = get(hLine, 'Color');
ls = get(hLine, 'LineStyle');
lw = get(hLine, 'LineWidth');

mk  = 'none';
ms  = 8;
mfc = 'none';
mec = c;

if isprop(hLine, 'Marker')
    mk = get(hLine, 'Marker');
end

if isprop(hLine, 'MarkerSize')
    ms = get(hLine, 'MarkerSize');
end

if isprop(hLine, 'MarkerFaceColor')
    mfc = get(hLine, 'MarkerFaceColor');
end

if isprop(hLine, 'MarkerEdgeColor')
    mec = get(hLine, 'MarkerEdgeColor');
end

% Convert normalized legend positions to data coordinates.
[x1, y1] = axes_norm_to_data(ax, x1n, yn);
[x2, ~ ] = axes_norm_to_data(ax, x2n, yn);
[xt, yt] = axes_norm_to_data(ax, xtn, yn);

xm = 0.49 * (x1 + x2);

% Draw the line sample, marker, and label.
line(ax, [x1 x2], [y1 y1], ...
    'LineStyle', ls, ...
    'Color', c, ...
    'LineWidth', lw, ...
    'Marker', 'none', ...
    'Clipping', 'off');

line(ax, xm, y1, ...
    'LineStyle', 'none', ...
    'Color', c, ...
    'Marker', mk, ...
    'MarkerSize', ms, ...
    'MarkerFaceColor', mfc, ...
    'MarkerEdgeColor', mec, ...
    'LineWidth', lw, ...
    'Clipping', 'off');

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

function err = iga_L2_error_example2(currRun, cNURBS, refRun, rNURBS, patchCenters, a, dx)
% Compute the IGA-region L2 error against the reference grid.

[X, Y, patchId, wA] = grid_points_inner_example2(patchCenters, a, dx);

vC = iga_eval_val_example2(currRun, cNURBS, X, Y, patchId, patchCenters, a);
vR = iga_eval_val_example2(refRun,  rNURBS, X, Y, patchId, patchCenters, a);

dv = vC - vR;
err = sqrt(sum(abs(dv).^2) * wA);

end

function err = pw_L2_error_example2(cPW, p_curr, rPW, p_ref, L, patchCenters, a, dx, chunkSize)
% Compute the plane-wave-region L2 error against the reference grid.

dy = dx;
x = -L/2 + dx/2 : dx : L/2 - dx/2;
y = x;
[Xg,Yg] = meshgrid(x,y);

mask = true(size(Xg));

for k = 1:size(patchCenters,1)
    xc = patchCenters(k,1);
    yc = patchCenters(k,2);

    mask = mask & ~((Xg >= xc-a) & (Xg <= xc+a) & ...
        (Yg >= yc-a) & (Yg <= yc+a));
end

X = Xg(mask);
Y = Yg(mask);
wA = dx*dy;

acc = 0;
nPts = numel(X);
k = 1;

while k <= nPts
    k2 = min(nPts, k+chunkSize-1);
    Xc = X(k:k2);
    Yc = Y(k:k2);
    F = [Xc(:)'; Yc(:)'];

    expoC = exp((1i*2*pi/L) * (p_curr * F));
    expoR = exp((1i*2*pi/L) * (p_ref  * F));

    vC = (cPW.' * expoC) / L;
    vR = (rPW.' * expoR) / L;

    dv = vC - vR;
    acc = acc + sum(abs(dv).^2);

    k = k2 + 1;
end

err = sqrt(acc * wA);

end

function errDG = compute_DG_error_example2( ...
currRun, cNURBS, cPW, ...
    refRun,  rNURBS, rPW, ...
    p_vec_curr, p_vec_ref, L, patchCenters, a, ...
    dx_in, dx_out, chunkSize, sigma)

% Compute the Example 2 DG eigenfunction error.
[xi, yi, patchId, wA_in] = grid_points_inner_example2(patchCenters, a, dx_in);

[vC, vxC, vyC] = iga_eval_val_grad_example2(currRun, cNURBS, xi, yi, patchId, patchCenters, a);
[vR, vxR, vyR] = iga_eval_val_grad_example2(refRun,  rNURBS, xi, yi, patchId, patchCenters, a);

de  = vC - vR;
dex = vxC - vxR;
dey = vyC - vyR;

H1_in = sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2) * wA_in;

[xo, yo, wA_out] = grid_points_outer_example2(L, patchCenters, a, dx_out);

H1_out = 0;
nPts = numel(xo);
k = 1;

while k <= nPts
    k2 = min(nPts, k+chunkSize-1);
    X = xo(k:k2);
    Y = yo(k:k2);

    [vC, vxC, vyC] = pw_eval_val_grad(cPW, p_vec_curr, X, Y, L);
    [vR, vxR, vyR] = pw_eval_val_grad(rPW, p_vec_ref,  X, Y, L);

    de  = vC - vR;
    dex = vxC - vxR;
    dey = vyC - vyR;

    H1_out = H1_out + sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2);

    k = k2 + 1;
end

H1_out = H1_out * wA_out;

[xg, yg, patchIdG, wL] = boundary_points_example2(patchCenters, a, dx_in);

[vCout, ~, ~] = pw_eval_val_grad(cPW, p_vec_curr, xg, yg, L);
[vRout, ~, ~] = pw_eval_val_grad(rPW, p_vec_ref,  xg, yg, L);

[vCin,  ~, ~] = iga_eval_val_grad_example2(currRun, cNURBS, xg, yg, patchIdG, patchCenters, a);
[vRin,  ~, ~] = iga_eval_val_grad_example2(refRun,  rNURBS, xg, yg, patchIdG, patchCenters, a);

jump = (vCout - vRout) - (vCin - vRin);
J2   = sum(abs(jump).^2) * wL;

errDG = sqrt(H1_in + H1_out + sigma * J2);

end

function [val, gx, gy] = pw_eval_val_grad(coeff, p_vec, X, Y, L)
% Evaluate the field value and gradient.

F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));

val = (coeff.' * expo) / L;

fac = (1i * 2*pi / L) / L;
gx  = (((coeff .* p_vec(:,1)).' * expo)) * fac;
gy  = (((coeff .* p_vec(:,2)).' * expo)) * fac;

val = val(:);
gx  = gx(:);
gy  = gy(:);

end

function val = iga_eval_val_example2(runData, coeff, X, Y, patchId, patchCenters, a)
% Evaluate the Example 2 field value.

[val, ~, ~] = iga_eval_val_grad_example2(runData, coeff, X, Y, patchId, patchCenters, a);

end

function [val, gx, gy] = iga_eval_val_grad_example2(runData, coeff, X, Y, patchId, patchCenters, a)
% Evaluate the Example 2 field value and gradient.

n1 = runData.n_dofs_1;
n2 = runData.n_dofs_2;

if numel(coeff) ~= n1 + n2
    error('Inconsistent numerical data.');
end

coeff1 = coeff(1:n1);
coeff2 = coeff(n1+1:n1+n2);

val = zeros(numel(X),1);
gx  = zeros(numel(X),1);
gy  = zeros(numel(X),1);

for k = 1:numel(X)
    pid = patchId(k);

    if pid == 1
        nurbs = runData.nurbs_refine_1;
        cLoc  = coeff1;
    elseif pid == 2
        nurbs = runData.nurbs_refine_2;
        cLoc  = coeff2;
    else
        error('Invalid patchId: %d', pid);
    end

    xc = patchCenters(pid,1);
    yc = patchCenters(pid,2);

    [val(k), gx(k), gy(k)] = iga_eval_on_one_patch( ...
        nurbs, cLoc, X(k), Y(k), xc, yc, a);
end

end

function [v, gx, gy] = iga_eval_on_one_patch(nurbs, coeff, X, Y, xc, yc, a)
% Evaluate one IGA patch at the requested Cartesian points.

pu = nurbs.pu;
pv = nurbs.pv;
U  = nurbs.Ubar(:).';
V  = nurbs.Vbar(:).';

mU = length(U) - pu - 1;
nV = length(V) - pv - 1;

u = (X - (xc-a)) / (2*a);
vpar = (Y - (yc-a)) / (2*a);

u = max(0, min(1, u));
vpar = max(0, min(1, vpar));

spanU = findspan_local(mU-1, pu, u, U);
spanV = findspan_local(nV-1, pv, vpar, V);

[Nu, dNu] = bspline_basis_and_der1(U, pu, u, spanU);
[Nv, dNv] = bspline_basis_and_der1(V, pv, vpar, spanV);

s  = 0 + 1i*0;
su = 0 + 1i*0;
sv = 0 + 1i*0;

for j1 = (spanV - pv) : spanV
    lv = j1 - (spanV - pv) + 1;

    for i1 = (spanU - pu) : spanU
        lu = i1 - (spanU - pu) + 1;
        row = i1 + (j1-1)*mU;
        c = coeff(row);

        s  = s  + c * Nu(lu)  * Nv(lv);
        su = su + c * dNu(lu) * Nv(lv);
        sv = sv + c * Nu(lu)  * dNv(lv);
    end
end

v  = s;
gx = su / (2*a);
gy = sv / (2*a);

end

function [N, dN] = bspline_basis_and_der1(U, p, u, span)
% Evaluate basis values and first derivatives.

ndu = zeros(p+1, p+1);
left = zeros(1, p+1);
right = zeros(1, p+1);

ndu(1,1) = 1.0;

for j = 1:p
    left(j+1)  = u - U(span+1-j);
    right(j+1) = U(span+j) - u;
    saved = 0.0;

    for r = 0:(j-1)
        ndu(j+1, r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1, j) / ndu(j+1, r+1);
        ndu(r+1, j+1) = saved + right(r+2) * temp;
        saved = left(j-r+1) * temp;
    end

    ndu(j+1, j+1) = saved;
end

N = ndu(1:p+1, p+1).';

ders1 = zeros(1, p+1);

for r = 0:p
    d = 0.0;
    pk = p-1;

    if r >= 1
        d = d + (1.0/ndu(pk+2, r)) * ndu(r, pk+1);
    end

    if r <= p-1
        d = d - (1.0/ndu(pk+2, r+1)) * ndu(r+1, pk+1);
    end

    ders1(r+1) = d;
end

dN = ders1 * p;

end

function span = findspan_local(n, p, u, U)
% Find the active knot span for a parameter value.

if u >= U(n+2)
    span = n+1;
    return;
end

if u <= U(p+1)
    span = p+1;
    return;
end

low = p+1;
high = n+2;
mid = floor((low+high)/2);

while (u < U(mid) || u >= U(mid+1))
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end

    mid = floor((low+high)/2);
end

span = mid;

end

function [X, Y, patchId, wA] = grid_points_inner_example2(patchCenters, a, dx)
% Build inner-region quadrature points for Example 2.

X = [];
Y = [];
patchId = [];

for p = 1:size(patchCenters,1)
    xc = patchCenters(p,1);
    yc = patchCenters(p,2);

    x = (xc-a) + dx/2 : dx : (xc+a) - dx/2;
    y = (yc-a) + dx/2 : dx : (yc+a) - dx/2;
    [Xg, Yg] = meshgrid(x,y);

    X = [X; Xg(:)]; %#ok<AGROW>
    Y = [Y; Yg(:)]; %#ok<AGROW>
    patchId = [patchId; p*ones(numel(Xg),1)]; %#ok<AGROW>
end

wA = dx*dx;

end

function [X, Y, wA] = grid_points_outer_example2(L, patchCenters, a, dx)
% Build outer-region quadrature points for Example 2.

x = -L/2 + dx/2 : dx : L/2 - dx/2;
y = x;
[Xg, Yg] = meshgrid(x,y);

mask = true(size(Xg));

for p = 1:size(patchCenters,1)
    xc = patchCenters(p,1);
    yc = patchCenters(p,2);

    mask = mask & ~((Xg >= xc-a) & (Xg <= xc+a) & ...
        (Yg >= yc-a) & (Yg <= yc+a));
end

X = Xg(mask);
Y = Yg(mask);
wA = dx*dx;

end

function [X, Y, patchId, wL] = boundary_points_example2(patchCenters, a, ds)
% Build interface quadrature points for Example 2.

X = [];
Y = [];
patchId = [];

for p = 1:size(patchCenters,1)
    xc = patchCenters(p,1);
    yc = patchCenters(p,2);

    t = -a + ds/2 : ds : a - ds/2;

    xb = xc + t;  yb = (yc-a)*ones(size(t));
    xt = xc + t;  yt = (yc+a)*ones(size(t));
    yl = yc + t;  xl = (xc-a)*ones(size(t));
    yr = yc + t;  xr = (xc+a)*ones(size(t));

    Xp = [xb, xt, xl, xr].';
    Yp = [yb, yt, yl, yr].';

    X = [X; Xp]; %#ok<AGROW>
    Y = [Y; Yp]; %#ok<AGROW>
    patchId = [patchId; p*ones(numel(Xp),1)]; %#ok<AGROW>
end

wL = ds;

end

function h = estimate_h_parametric_example2(runData)
% Estimate the Example 2 parametric mesh size.

h1 = estimate_h_parametric_one(runData.nurbs_refine_1);
h2 = estimate_h_parametric_one(runData.nurbs_refine_2);
h = max(h1, h2);

end

function h = estimate_h_parametric_one(nurbs)
% Estimate one parametric mesh size.

Uu = unique(nurbs.Ubar(:).');
Vv = unique(nurbs.Vbar(:).');
h  = max(max(diff(Uu)), max(diff(Vv)));

end

function [p_vec, n_pw_basis] = generate_p_vec(Nc)
% Build the two-dimensional plane-wave index set for one cutoff.

N = floor(Nc);
p_vec = zeros((2*N+1)^2, 2);
n_pw_basis = 0;

for ii = -N:N
    m = floor(sqrt(N^2 - ii^2));

    for jj = -m:m
        n_pw_basis = n_pw_basis + 1;
        p_vec(n_pw_basis, :) = [ii, jj];
    end
end

p_vec = p_vec(1:n_pw_basis, :);

end

function runs = load_runs_over_refine_example2(ExampleRoot, Nc_fixed, pdeg, refine_list, n_pw_curr)
% Load the retained Example 2 runs over all refinement levels.

template = struct('refine', 0, ...
    'nurbs_original_1', [], 'nurbs_refine_1', [], ...
    'nurbs_original_2', [], 'nurbs_refine_2', [], ...
    'n_dofs_1', 0, 'n_dofs_2', 0, 'nNURBS', 0, ...
    'uh', []);
runs = repmat(template, numel(refine_list), 1);

for i = 1:numel(refine_list)
    rf = refine_list(i);

    runMat = fullfile(ExampleRoot, sprintf('Nc_%d',Nc_fixed), sprintf('p_%d',pdeg), ...
        sprintf('refine_%02d',rf), 'run.mat');

    assert(isfile(runMat), 'Missing run file: %s', runMat);
    S = load(runMat, 'run');
    assert(isfield(S, 'run'), 'MAT file does not contain run: %s', runMat);
    R = S.run;

    rr.refine = rf;

    rr.nurbs_original_1 = R.nurbs_original_1;
    rr.nurbs_refine_1   = R.nurbs_refine_1;
    rr.nurbs_original_2 = R.nurbs_original_2;
    rr.nurbs_refine_2   = R.nurbs_refine_2;
    rr.uh               = R.uh;

    requiredFields = {'n_dofs_1', 'n_dofs_2', ...
        'n_dofs_nurbs', 'n_pw_basis'};
    assert(all(isfield(R, requiredFields)), ...
        'Run is missing dimension metadata: %s', runMat);
    rr.n_dofs_1 = double(R.n_dofs_1);
    rr.n_dofs_2 = double(R.n_dofs_2);
    rr.nNURBS = double(R.n_dofs_nurbs);
    assert(rr.nNURBS == rr.n_dofs_1 + rr.n_dofs_2, ...
        'refine=%d has inconsistent NURBS degrees of freedom.', rf);
    assert(double(R.n_pw_basis) == n_pw_curr, ...
        'refine=%d has an inconsistent plane-wave dimension.', rf);
    runs(i) = rr;
end

end

function prepare_example_2_plot_context(scriptFile)
% Prepare paths relative to this example script.

figuresDir = fileparts(scriptFile);
exampleDir = fileparts(figuresDir);
projectDir = fileparts(fileparts(exampleDir));
utilsDir   = fullfile(projectDir, 'src', 'utils');
dataDir    = fullfile(exampleDir, 'data');

assert(isfolder(utilsDir), 'Missing utility directory: %s', utilsDir);
assert(isfolder(dataDir), 'Missing data directory: %s', dataDir);

addpath(utilsDir, '-begin');
add_example_paths(exampleDir);
cd(dataDir);

end

function dxTag = make_dx_tag(dx_in, dx_out)
% Encode the inner and outer grid spacings in a filename tag.

sx = dx_to_str(dx_in);
sy = dx_to_str(dx_out);

if abs(dx_in - dx_out) < eps(max(dx_in, dx_out))
    dxTag = sprintf('dx_%s', sx);
else
    dxTag = sprintf('dxIn_%s_dxOut_%s', sx, sy);
end

end

function s = dx_to_str(dx)
% Convert a numeric value to compact text.

s = sprintf('%.1e', dx);
s = strrep(s, '.', 'p');
s = strrep(s, '+', '');

end

function y = round_sig(x, nSig)
% Round a value to significant digits.

y = x;
mask = isfinite(x) & (x ~= 0);

ax = abs(x(mask));
p  = floor(log10(ax));
scale = 10.^(nSig-1-p);

y(mask) = round(x(mask).*scale)./scale;

end
