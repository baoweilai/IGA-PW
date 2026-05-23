function plot_eigen_error()
%Plot eigen error.

clc; close all; format short e;

%% --------------------- Global interpreter -------------------------------
set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');

%% Section
Example     = 'Example_1';
Nc_fixed    = 30;
p_list      = [1 2];
refine_list = 2:6;

ref_Nc      = 45;
ref_pdeg    = 3;
ref_refine  = 8;

u_list      = [1 2];

L           = 4;
a           = 0.2;

Csigma      = 10;
dx_in       = 5e-3;
dx_out      = 5e-3;
chunkSize   = 20000;

slope1_factor = [0.07,   0.05];
slope2_factor = [0.0008, 0.001];

padX = 0.4;
padY = 10;

savePNG = true;
savePDF = true;
pngDPI  = 600;

keepFigureOpen = true;
forceRecompute = true;

%% Section
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
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';

cfg.curve.lineWidth       = 2.0;
cfg.curve.markerSize      = 7;
cfg.curve.markerFaceColor = 'w';

cfg.label.fontSize = 12;

%% Section
% ---------- p = 1 ----------
axisSpec(1).XTick      = [];
axisSpec(1).XTickLabel = {};
axisSpec(1).XLim       = [];

axisSpec(1).YTick      = 10.^[-6 -4 -2 0];
axisSpec(1).YTickLabel = {'$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$'};
axisSpec(1).YLim       = [1e-6, 1e-0];

% ---------- p = 2 ----------
axisSpec(2).XTick      = [];
axisSpec(2).XTickLabel = {};
axisSpec(2).XLim       = [];

axisSpec(2).YTick      = 10.^[-6 -4 -2 0];
axisSpec(2).YTickLabel = {'$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$'};
axisSpec(2).YLim       = [1e-6, 1e-0];

%% Section
legendSpec(1).xL1         = 0.37;
legendSpec(1).xL2         = 0.47;
legendSpec(1).xLT         = 0.48;
legendSpec(1).xR1         = 0.67;
legendSpec(1).xR2         = 0.83;
legendSpec(1).xRT         = 0.85;
legendSpec(1).rowY        = [0.4, 0.3, 0.2, 0.1];
legendSpec(1).slopeRows   = [3, 4];
legendSpec(1).slope1Label = '$\mathrm{Slope} = 1$';
legendSpec(1).slope2Label = '$\mathrm{Slope} = 2$';

legendSpec(2).xL1         = 0.37;
legendSpec(2).xL2         = 0.47;
legendSpec(2).xLT         = 0.48;
legendSpec(2).xR1         = 0.67;
legendSpec(2).xR2         = 0.83;
legendSpec(2).xRT         = 0.84;
legendSpec(2).rowY        = [0.4, 0.3, 0.2, 0.1];
legendSpec(2).slopeRows   = [3, 4];
legendSpec(2).slope1Label = '$\mathrm{Slope} = 1$';
legendSpec(2).slope2Label = '$\mathrm{Slope} = 2$';

%% Section
lineColors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;

slope1Color = [033 158 188] / 255;
slope2Color = [239 065 067] / 255;

%% Section
activate_example_workflow('h_convergence', ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core', 'solver'});

ExampleRoot = find_example_root(pwd, Example);

dxTag  = make_dx_tag(dx_in, dx_out);
refTag = make_refine_tag(refine_list);

outDir = fullfile(ExampleRoot, 'eigen_error', ...
    sprintf('Nc_%d', Nc_fixed), dxTag);

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% Section
[p_vec_curr, n_pw_curr] = generate_p_vec(Nc_fixed);
[p_vec_ref,  n_pw_ref ] = generate_p_vec(ref_Nc);

[tf_k, loc_k] = ismember(p_vec_curr, p_vec_ref, 'rows');

%% Section
refLoaded = false;
ref = struct();

%% Section
for ip = 1:numel(p_list)
    pdeg = p_list(ip);
    fprintf('\n==================== p = %d ====================\n', pdeg);

    pOutDir = fullfile(outDir, sprintf('p_%d', pdeg));
    if ~exist(pOutDir, 'dir'), mkdir(pOutDir); end
    errCsvFile   = fullfile(pOutDir, 'errors.csv');
    rateCsvFile  = fullfile(pOutDir, 'rates.csv');
    pngFile      = fullfile(pOutDir, 'eigen_error.png');
    pdfFile      = fullfile(pOutDir, 'eigen_error.pdf');

    needRecompute = true;
    T = [];

    if ~forceRecompute && exist(errCsvFile, 'file')

        T_all = readtable(errCsvFile);
        T = T_all(ismember(T_all.refine, refine_list), :);

        if isempty(T)
            needRecompute = true;
        else
            refine_in_cache = unique(T.refine(:)).';
            refine_need     = unique(refine_list(:)).';
            if all(ismember(refine_need, refine_in_cache))
                needRecompute = false;
            else
                needRecompute = true;
            end
        end
    end

    if needRecompute

        if ~refLoaded
            ref_runMat = fullfile(ExampleRoot, sprintf('Nc_%d', ref_Nc), ...
                sprintf('p_%d', ref_pdeg), sprintf('refine_%02d', ref_refine), 'run.mat');

            if ~exist(ref_runMat, 'file')
                error('Missing reference run file: %s', ref_runMat);
            end

            S = load(ref_runMat);
            if isfield(S,'run')
                Rref = S.run;
            else
                Rref = S;
            end

            ref.nurbs_refine = Rref.nurbs_refine;
            ref.uh           = Rref.uh;

            if isfield(Rref,'n_dofs_nurbs') && ~isempty(Rref.n_dofs_nurbs)
                ref.nNURBS = double(Rref.n_dofs_nurbs);
            else
                ref.nNURBS = size(Rref.uh,1) - n_pw_ref;
            end

            fprintf('[REF ] Nc=%d p=%d refine=%d nNURBS=%d\n', ...
                ref_Nc, ref_pdeg, ref_refine, ref.nNURBS);

            refLoaded = true;
        end

        runs = load_runs_over_refine(ExampleRoot, Nc_fixed, pdeg, refine_list, n_pw_curr);
        if isempty(runs)
            error('p=%d has missing or invalid run data.', pdeg);
            continue;
        end

        nR = numel(runs);

        dx_axis = zeros(nR,1);
        eL2     = zeros(nR, numel(u_list));
        eDG     = zeros(nR, numel(u_list));
        refine_loaded = zeros(nR,1);

        for i = 1:nR
            rr = runs(i);

            refine_loaded(i) = rr.refine;

            h_param = estimate_h_parametric(rr.nurbs_refine);
            h_phys = 2 * a * h_param;
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

                err_IGA_L2 = Compute_L2_Error( ...
                    rr.nurbs_original, rr.nurbs_refine, curr_NURBS, ...
                    ref.nurbs_refine, ref_NURBS);

                err_PW_L2 = pw_L2_error(curr_PW, p_vec_curr, ref_PW, p_vec_ref, ...
                    L, a, dx_out, chunkSize);

                eL2(i,ju) = sqrt(err_IGA_L2^2 + err_PW_L2^2);

                eDG(i,ju) = compute_DG_error( ...
                    rr.nurbs_refine, curr_NURBS, curr_PW, ...
                    ref.nurbs_refine, ref_NURBS, ref_PW, ...
                    p_vec_curr, p_vec_ref, L, a, ...
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
            eL2(:,1), eDG(:,1), eL2(:,2), eDG(:,2), ...
            'VariableNames', {'refine','dx','u1_L2','u1_DG','u2_L2','u2_DG'});

        writetable(T, errCsvFile);
    end

    T = T(ismember(T.refine, refine_list), :);
    T = sortrows(T, 'refine');

    if isempty(T)
        error('p=%d has missing or invalid run data.', pdeg);
        continue;
    end

    Y = [T.u1_L2, T.u1_DG, T.u2_L2, T.u2_DG];
    x_axis = T.dx;
    Y(Y <= 0) = eps;

    T_rate = build_local_rate_table(T, pdeg);
    writetable(T_rate, rateCsvFile);

    curveLabels = { ...
        '$u_1,\ L^2$', ...
        '$u_1,\ \mathrm{DG}$', ...
        '$u_2,\ L^2$', ...
        '$u_2,\ \mathrm{DG}$'};

    %% --------------------- Plot ---------------------------------------
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
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

    mk = {'o','s','^','d'};
    hCurve = gobjects(4,1);

    for k = 1:4
        col = lineColors(k,:);
        hCurve(k) = plot(ax, x_axis, Y(:,k), '-', ...
            'LineWidth', cfg.curve.lineWidth, ...
            'Color', col, ...
            'Marker', mk{k}, ...
            'MarkerSize', cfg.curve.markerSize, ...
            'MarkerFaceColor', cfg.curve.markerFaceColor, ...
            'MarkerEdgeColor', col);
    end

    %% Section
    xGrid = unique(x_axis(:));
    anchor = choose_anchor_index(xGrid);
    h0 = xGrid(anchor);
    yMax = max(Y(:));

    yref1 = (slope1_factor(ip) * yMax) * (xGrid ./ h0).^1;
    hSlope1 = plot(ax, xGrid, yref1, '--', ...
        'Color', slope1Color, ...
        'LineWidth', 1.8);

    yref2 = (slope2_factor(ip) * yMax) * (xGrid ./ h0).^2;
    hSlope2 = plot(ax, xGrid, yref2, '--', ...
        'Color', slope2Color, ...
        'LineWidth', 1.8);

    xlabel(ax, '$h$', ...
        'FontSize', cfg.label.fontSize, ...
        'Interpreter', 'latex');

    ylabel(ax, '$\|u_i-u_{i}^{\mathrm{DG}}\|$', ...
        'FontSize', cfg.label.fontSize, ...
        'Interpreter', 'latex');

    apply_manual_log_axes(ax, x_axis, Y, yref1, yref2, axisSpec(pdeg), padX, padY);
    draw_fake_legend_in_axes(ax, hCurve, hSlope1, hSlope2, curveLabels, legendSpec(pdeg));

    if savePNG
        exportgraphics(fig, pngFile, 'Resolution', pngDPI);
        fprintf('[SAVE] %s\n', pngFile);
    end

    if savePDF
        exportgraphics(fig, pdfFile, 'ContentType', 'vector');
        fprintf('[SAVE] %s\n', pdfFile);
    end

    if ~keepFigureOpen
        close(fig);
    end
end

end

function T_rate = build_local_rate_table(T, pdeg)
%Build local rate table.

metricNames = {'u1_L2','u1_DG','u2_L2','u2_DG'};
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
%Apply manual log axes.

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
%Choose the point used for slope annotation.
n = numel(xGrid);
if n <= 2
    idx = 1;
else
    idx = round(n/2);
end
end

function draw_fake_legend_in_axes(ax, hCurve, hSlope1, hSlope2, curveLabels, spec)
%Draw fake legend in axes.

xL1 = spec.xL1;  xL2 = spec.xL2;  xLT = spec.xLT;
xR1 = spec.xR1;  xR2 = spec.xR2;  xRT = spec.xRT;
rowY = spec.rowY;

for j = 1:4
    draw_one_fake_entry_in_axes(ax, hCurve(j), xR1, xR2, xRT, rowY(j), curveLabels{j});
end

draw_one_fake_entry_in_axes(ax, hSlope1, xL1, xL2, xLT, rowY(spec.slopeRows(1)), spec.slope1Label);
draw_one_fake_entry_in_axes(ax, hSlope2, xL1, xL2, xLT, rowY(spec.slopeRows(2)), spec.slope2Label);

end

function draw_one_fake_entry_in_axes(ax, hLine, x1n, x2n, xtn, yn, labelStr)
%Draw one fake entry in axes.

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

[x1, y1] = axes_norm_to_data(ax, x1n, yn);
[x2, y2] = axes_norm_to_data(ax, x2n, yn);
[xt, yt] = axes_norm_to_data(ax, xtn, yn);

xm = 0.49 * (x1 + x2);
ym = y1;

line(ax, [x1 x2], [y1 y2], ...
    'LineStyle', ls, ...
    'Color', c, ...
    'LineWidth', lw, ...
    'Marker', 'none', ...
    'Clipping', 'off');

line(ax, xm, ym, ...
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
    'FontSize', 11, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', ...
    'Margin', 0.5, ...
    'Clipping', 'off');

end

function [x, y] = axes_norm_to_data(ax, xn, yn)
%Compute norm to data.

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

function err = pw_L2_error(cPW, p_curr, rPW, p_ref, L, a, dx, chunkSize)
%Compute L2 error.
dy = dx;
x = -L/2 + dx/2 : dx : L/2 - dx/2;
y = x;
[Xg,Yg] = meshgrid(x,y);
mask = ~(Xg >= -a & Xg <= a & Yg >= -a & Yg <= a);
X = Xg(mask); Y = Yg(mask);
wA = dx*dy;

acc = 0;
nPts = numel(X);
k = 1;
while k <= nPts
    k2 = min(nPts, k+chunkSize-1);
    Xc = X(k:k2); Yc = Y(k:k2);
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

function errDG = compute_DG_error( ...
nurbs_curr, cNURBS, cPW, ...
    nurbs_ref,  rNURBS, rPW, ...
    p_vec_curr, p_vec_ref, L, a, ...
    dx_in, dx_out, chunkSize, sigma)

[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
[vC, vxC, vyC]  = iga_eval_val_grad(nurbs_curr, cNURBS, xi, yi, a);
[vR, vxR, vyR]  = iga_eval_val_grad(nurbs_ref,  rNURBS, xi, yi, a);
de  = vC - vR;  dex = vxC - vxR;  dey = vyC - vyR;
H1_in = sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2) * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
H1_out = 0;
nPts = numel(xo);
k = 1;
while k <= nPts
    k2 = min(nPts, k+chunkSize-1);
    X = xo(k:k2); Y = yo(k:k2);

    [vC, vxC, vyC] = pw_eval_val_grad(cPW, p_vec_curr, X, Y, L);
    [vR, vxR, vyR] = pw_eval_val_grad(rPW, p_vec_ref,  X, Y, L);

    de  = vC - vR;  dex = vxC - vxR;  dey = vyC - vyR;
    H1_out = H1_out + sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2);
    k = k2 + 1;
end
H1_out = H1_out * wA_out;

[xg, yg, wL] = boundary_points_square(a, dx_in);
[vCout, ~, ~] = pw_eval_val_grad(cPW, p_vec_curr, xg, yg, L);
[vRout, ~, ~] = pw_eval_val_grad(rPW, p_vec_ref,  xg, yg, L);
[vCin,  ~, ~] = iga_eval_val_grad(nurbs_curr, cNURBS, xg, yg, a);
[vRin,  ~, ~] = iga_eval_val_grad(nurbs_ref,  rNURBS, xg, yg, a);

jump = (vCout - vRout) - (vCin - vRin);
J2   = sum(abs(jump).^2) * wL;

errDG = sqrt(H1_in + H1_out + sigma * J2);
end

function [val, gx, gy] = pw_eval_val_grad(coeff, p_vec, X, Y, L)
%Evaluate the field value and gradient.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));
val = (coeff.' * expo) / L;

fac = (1i * 2*pi / L) / L;
gx  = (((coeff .* p_vec(:,1)).' * expo)) * fac;
gy  = (((coeff .* p_vec(:,2)).' * expo)) * fac;

val = val(:); gx = gx(:); gy = gy(:);
end

function [val, gx, gy] = iga_eval_val_grad(nurbs, coeff, X, Y, a)
%Evaluate the field value and gradient.
pu = nurbs.pu;  pv = nurbs.pv;
U  = nurbs.Ubar(:).';
V  = nurbs.Vbar(:).';
mU = length(U) - pu - 1;
nV = length(V) - pv - 1;

val = zeros(numel(X),1);
gx  = zeros(numel(X),1);
gy  = zeros(numel(X),1);

for k = 1:numel(X)
    u = (X(k) + a) / (2*a);
    v = (Y(k) + a) / (2*a);
    u = max(0, min(1, u));
    v = max(0, min(1, v));

    spanU = findspan_local(mU-1, pu, u, U);
    spanV = findspan_local(nV-1, pv, v, V);

    [Nu, dNu] = bspline_basis_and_der1(U, pu, u, spanU);
    [Nv, dNv] = bspline_basis_and_der1(V, pv, v, spanV);

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

    val(k) = s;
    gx(k)  = su / (2*a);
    gy(k)  = sv / (2*a);
end
end

function [N, dN] = bspline_basis_and_der1(U, p, u, span)
%Evaluate basis values and first derivatives.
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
%Locate an index or object used by the computation.
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

%% ============================ sampling grids ============================
function [X, Y, wA] = grid_points_square(xmin, xmax, dx)
%Build midpoint quadrature points in the square.
x = xmin + dx/2 : dx : xmax - dx/2;
y = x;
[Xg, Yg] = meshgrid(x,y);
X = Xg(:);
Y = Yg(:);
wA = dx*dx;
end

function [X, Y, wA] = grid_points_outer(L, a, dx)
%Build midpoint quadrature points outside the square.
x = -L/2 + dx/2 : dx : L/2 - dx/2;
y = x;
[Xg, Yg] = meshgrid(x,y);
mask = ~(Xg >= -a & Xg <= a & Yg >= -a & Yg <= a);
X = Xg(mask);
Y = Yg(mask);
wA = dx*dx;
end

function [X, Y, wL] = boundary_points_square(a, ds)
%Build midpoint quadrature points on the square boundary.
t  = -a + ds/2 : ds : a - ds/2;
xb = t;  yb = -a*ones(size(t));
xt = t;  yt =  a*ones(size(t));
yl = t;  xl = -a*ones(size(t));
yr = t;  xr =  a*ones(size(t));

X = [xb, xt, xl, xr].';
Y = [yb, yt, yl, yr].';
wL = ds;
end

%% ============================ misc helpers =============================
function h = estimate_h_parametric(nurbs_refine)
%Estimate the parametric mesh size.
Uu = unique(nurbs_refine.Ubar(:).');
Vv = unique(nurbs_refine.Vbar(:).');
h  = max(max(diff(Uu)), max(diff(Vv)));
end

function [p_vec, n_pw_basis] = generate_p_vec(Nc)
%Generate p vec.
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

function runs = load_runs_over_refine(ExampleRoot, Nc_fixed, pdeg, refine_list, n_pw_curr)
%Load runs over refine.
runs = struct('refine',{},'nurbs_original',{},'nurbs_refine',{},'uh',{},'nNURBS',{});

for i = 1:numel(refine_list)
    rf = refine_list(i);
    runMat = fullfile(ExampleRoot, sprintf('Nc_%d',Nc_fixed), sprintf('p_%d',pdeg), ...
        sprintf('refine_%02d',rf), 'run.mat');

    if ~exist(runMat,'file')
        continue;
    end

    S = load(runMat);
    if isfield(S,'run')
        R = S.run;
    else
        R = S;
    end

    rr.refine = rf;
    rr.nurbs_original = R.nurbs_original;
    rr.nurbs_refine   = R.nurbs_refine;
    rr.uh             = R.uh;

    if isfield(R,'n_dofs_nurbs') && ~isempty(R.n_dofs_nurbs)
        rr.nNURBS = double(R.n_dofs_nurbs);
    else
        rr.nNURBS = size(R.uh,1) - n_pw_curr;
    end

    runs(end+1) = rr; %#ok<AGROW>
end
end

function root = find_example_root(startDir, Example)
%Locate an index or object used by the computation.
root = startDir;
for k = 1:10
    if ~isempty(dir(fullfile(root,'Nc_*')))
        return;
    end

    cand = fullfile(root,'result',Example);
    if exist(cand,'dir') && ~isempty(dir(fullfile(cand,'Nc_*')))
        root = cand;
        return;
    end

    parent = fileparts(root);
    if strcmp(parent,root)
        break;
    end
    root = parent;
end

error('Cannot locate data directory for %s.', Example);
end

function dxTag = make_dx_tag(dx_in, dx_out)
%Build dx tag.
sx = dx_to_str(dx_in);
sy = dx_to_str(dx_out);

if abs(dx_in - dx_out) < eps(max(dx_in, dx_out))
    dxTag = sprintf('dx_%s', sx);
else
    dxTag = sprintf('dxIn_%s_dxOut_%s', sx, sy);
end
end

function refTag = make_refine_tag(refine_list)
%Build refine tag.
refine_list = unique(refine_list(:).');
refine_str = sprintf('%d_', refine_list);
refine_str(end) = [];
refTag = ['r_' refine_str];
end

function s = dx_to_str(dx)
%Compute to str.
s = sprintf('%.1e', dx);
s = strrep(s, '.', 'p');
s = strrep(s, '+', '');
end

function y = round_sig(x, nSig)
%Round a value to significant digits.
y = x;
mask = isfinite(x) & (x ~= 0);
ax = abs(x(mask));
p  = floor(log10(ax));
scale = 10.^(nSig-1-p);
y(mask) = round(x(mask).*scale)./scale;
end
