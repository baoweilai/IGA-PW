function plot_dg_error()
% Plot DG error.

clc; close all;

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');

%% ---------------- user params ----------------
Example    = 'Example_3';
Nc_fixed   = 20;
p_list     = [1 2];
nElem_list = [4 8 16 32 64];

% -------- reference --------
ref_Nc     = 40;
ref_p      = 2;
ref_refine = 7;

L         = 4;
a         = 0.2;
beta      = 100;

dx_in     = 5e-3;
dx_out    = 5e-3;
chunkSize = 20000;

%% ---------------- plot params ----------------
cfg = struct();

cfg.fig.width    = 4.8;
cfg.fig.height   = 3.0;
cfg.fig.renderer = 'painters';
cfg.fig.bgColor  = 'w';

cfg.layout.left   = 0.14;
cfg.layout.right  = 0.04;
cfg.layout.bottom = 0.16;
cfg.layout.top    = 0.08;

cfg.curve.colors = [ ...
    223 122 094;
    060 064 091] / 255;

cfg.curve.markers = {'s','o'};
cfg.curve.lineWidth = 2.0;
cfg.curve.markerSize = 7;
cfg.curve.markerFaceColor = 'w';

cfg.axes.fontSize   = 11;
cfg.axes.lineWidth  = 1.0;
cfg.axes.tickDir    = 'out';
cfg.axes.xScale     = 'log';
cfg.axes.yScale     = 'log';
cfg.axes.xMinorTick = 'off';
cfg.axes.yMinorTick = 'off';
cfg.axes.xlabel     = '$h$';
cfg.axes.ylabel     = '$\|u_1-u_{1}^{\mathrm{DG}}\|_{\mathrm{DG}}$';
cfg.axes.labelSize  = 12;

cfg.range.padX = 0.2;
cfg.range.padY = 0.6;

cfg.axis.XTick      = [1e-2 1e-1];
cfg.axis.XTickLabel = {'$10^{-2}$','$10^{-1}$'};
cfg.axis.XLim       = [];
cfg.axis.YTick      = [1e-2 1e-1];
cfg.axis.YTickLabel = {'$10^{-2}$','$10^{-1}$'};
cfg.axis.YLim       = [];

cfg.slope.color1    = [033 158 188] / 255;
cfg.slope.lineStyle = '--';
cfg.slope.lineWidth = 1.8;
cfg.slope.order1    = 1;
cfg.slope.factor1   = 0.5;

cfg.legend.location = 'southeast';
cfg.legend.fontSize = 11;
cfg.legend.box      = 'off';

%% ---------------- paths ----------------
resultRoot = fullfile(pwd, 'result', Example, sprintf('Nc_%02d', Nc_fixed));
if ~exist(resultRoot, 'dir')
    error('Directory not found: %s', resultRoot);
end

outDir = fullfile(resultRoot, 'eigen_error');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ---------------- reference ----------------
refRunMat = fullfile(pwd, 'result', Example, ...
    sprintf('refine_%02d', ref_refine), ...
    sprintf('p_%d', ref_p), ...
    sprintf('Nc_%02d', ref_Nc), ...
    'run.mat');

if ~exist(refRunMat, 'file')
    error('Reference run.mat not found: %s', refRunMat);
end

ref = load_run_data(refRunMat);
ref.u  = ref.uh(:,1);
ref.uI = ref.u(1:ref.nNURBS);
ref.uA = ref.u(ref.nNURBS+1 : ref.nNURBS + size(ref.k_pw,1));

fprintf('[REF] %s\n', refRunMat);

%% ---------------- collect plot data for all p ----------------
plotData = cell(numel(p_list),1);

%% ---------------- main loop: cache per p ----------------
for ip = 1:numel(p_list)
    pdeg = p_list(ip);

    pOutDir = fullfile(outDir, sprintf('p_%d', pdeg));
    if ~exist(pOutDir, 'dir'), mkdir(pOutDir); end
    csvFile = fullfile(pOutDir, 'errors.csv');
    matFile = fullfile(pOutDir, 'errors.mat');

    T_old = table();
    if exist(csvFile, 'file')
        T_old = readtable(csvFile);
        fprintf('[CACHE] %s\n', csvFile);
    end

    rows = table();
    for nElem = nElem_list(:).'
        idx = [];
        if ~isempty(T_old) && ismember('nElem', T_old.Properties.VariableNames)
            idx = find(T_old.nElem == nElem, 1, 'last');
        end

        if ~isempty(idx) && all(ismember({'nElem','h','u1_DG'}, T_old.Properties.VariableNames)) ...
                && isfinite(T_old.h(idx)) && isfinite(T_old.u1_DG(idx))
            one = T_old(idx, {'nElem','h','u1_DG'});
            rows = [rows; one]; %#ok<AGROW>
            fprintf('[CACHE] p=%d, nElem=%02d: h=%.4e, DG=%.3e\n', ...
                pdeg, one.nElem, one.h, one.u1_DG);
            continue;
        end

        runMat = fullfile(resultRoot, ...
            sprintf('p_%d', pdeg), ...
            sprintf('nElem_%02d', nElem), ...
            'run.mat');

        assert(exist(runMat, 'file') == 2, 'Missing run data: %s', runMat);

        rr = load_run_data(runMat);
        u  = rr.uh(:,1);

        uI = u(1:rr.nNURBS);
        uA = u(rr.nNURBS+1 : rr.nNURBS + size(rr.k_pw,1));

        [tf, loc] = ismember(rr.k_pw, ref.k_pw, 'rows');
        if ~all(tf)
            error('Current PW index set is not contained in reference PW index set.');
        end

        uA_pad = zeros(size(ref.uA));
        uA_pad(loc(tf)) = uA(tf);

        alpha = ref.uA' * uA_pad;
        phase = 1;
        if abs(alpha) > 1e-14
            phase = exp(-1i * angle(alpha));
        end

        uI = uI * phase;
        uA = uA * phase;

        h = rr.h;
        assert(isfinite(h) && h > 0, 'Mesh size h must be positive.');
        sigma = beta * (Nc_fixed + 1 / h);

        errDG = compute_DG_error( ...
            rr.nurbs_refine, uI, uA, rr.k_pw, ...
            ref.nurbs_refine, ref.uI, ref.uA, ref.k_pw, ...
            L, a, dx_in, dx_out, chunkSize, sigma);

        one = table(nElem, h, errDG, ...
            'VariableNames', {'nElem','h','u1_DG'});
        rows = [rows; one]; %#ok<AGROW>

        fprintf('[CALC ] p=%d, nElem=%02d: h=%.4e, DG=%.3e\n', ...
            pdeg, nElem, h, errDG);
    end

    assert(~isempty(rows), 'No valid runs found for p = %d.', pdeg);

    [~, ia] = unique(rows.nElem, 'last');
    rows = rows(sort(ia), :);
    rows = sortrows(rows, 'h', 'descend');

    writetable(rows, csvFile);
    T = rows;
    save(matFile, 'T');

    h_desc   = rows.h;
    eDG_desc = rows.u1_DG;
    n_list   = rows.nElem;

    assert(all(eDG_desc > 0), 'DG errors must be positive.');

    fprintf('\n========================================\n');
    fprintf('Example 3, Nc = %d, p = %d\n', Nc_fixed, pdeg);
    fprintf('----------------------------------------\n');
    for k = 1:numel(h_desc)-1
        ordDG = log(eDG_desc(k+1)/eDG_desc(k)) / log(h_desc(k+1)/h_desc(k));
        fprintf('nElem pair (%d -> %d): h = %.4e -> %.4e, DG order = %.6f\n', ...
            n_list(k), n_list(k+1), h_desc(k), h_desc(k+1), ordDG);
    end

    [h_plot, idx] = sort(rows.h, 'ascend');
    plotData{ip}.p  = pdeg;
    plotData{ip}.h  = h_plot;
    assert(all(rows.u1_DG(idx) > 0), 'DG errors must be positive.');
    plotData{ip}.DG = rows.u1_DG(idx);
end

%% ---------------- combined plot for p=1 and p=2 ----------------
validMask = ~cellfun(@isempty, plotData);
plotData = plotData(validMask);

assert(~isempty(plotData), 'No valid plot data found.');

pdfFile = fullfile(outDir, 'dg_error.pdf');

fig = figure('Color', cfg.fig.bgColor, ...
    'Units', 'inches', ...
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
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'Box', 'on', ...
    'XScale', cfg.axes.xScale, ...
    'YScale', cfg.axes.yScale, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick);

hList = gobjects(0);
leg = {};
xall = [];
yall = [];

for ip = 1:numel(plotData)
    pdeg = plotData{ip}.p;
    h    = plotData{ip}.h;
    eDG  = plotData{ip}.DG;

    h2 = plot(ax, h, eDG, '-', ...
        'LineWidth', cfg.curve.lineWidth, ...
        'Color', cfg.curve.colors(ip,:), ...
        'Marker', cfg.curve.markers{ip}, ...
        'MarkerSize', cfg.curve.markerSize, ...
        'MarkerFaceColor', cfg.curve.markerFaceColor, ...
        'MarkerEdgeColor', cfg.curve.colors(ip,:));
    hList(end+1) = h2; %#ok<AGROW>
    leg{end+1} = sprintf('$p=%d$', pdeg); %#ok<AGROW>

    xall = [xall; h(:)]; %#ok<AGROW>
    yall = [yall; eDG(:)]; %#ok<AGROW>
end

% 1st-order reference line
xGrid = unique(sort(xall(:)));
anchor = choose_anchor_index(xGrid);
h0 = xGrid(anchor);
yMax = max(yall(:));
yref1 = (cfg.slope.factor1 * yMax) * (xGrid ./ h0).^cfg.slope.order1;

hSlope = plot(ax, xGrid, yref1, cfg.slope.lineStyle, ...
    'Color', cfg.slope.color1, ...
    'LineWidth', cfg.slope.lineWidth);
hList(end+1) = hSlope;
leg{end+1} = '$\mathrm{Slope} = 1$';

xlabel(ax, cfg.axes.xlabel, 'FontSize', cfg.axes.labelSize, 'Interpreter', 'latex');
ylabel(ax, cfg.axes.ylabel, 'FontSize', cfg.axes.labelSize, 'Interpreter', 'latex');

xall = xall(isfinite(xall) & xall > 0);
yall = [yall; yref1(:)];
yall = yall(isfinite(yall) & yall > 0);

if isempty(cfg.axis.XLim)
    xmin = min(xall); xmax = max(xall);
    if xmin == xmax
        ax.XLim = [xmin/1.5, xmax*1.5];
    else
        ax.XLim = [xmin/(1+cfg.range.padX), xmax*(1+cfg.range.padX)];
    end
else
    ax.XLim = cfg.axis.XLim;
end

if isempty(cfg.axis.YLim)
    ymin = min(yall); ymax = max(yall);
    if ymin == ymax
        ax.YLim = [ymin/10, ymax*10];
    else
        ax.YLim = [ymin/(1+cfg.range.padY), ymax*(1+cfg.range.padY)];
    end
else
    ax.YLim = cfg.axis.YLim;
end

if ~isempty(cfg.axis.XTick), ax.XTick = cfg.axis.XTick; end
if ~isempty(cfg.axis.XTickLabel), ax.XTickLabel = cfg.axis.XTickLabel; end
if ~isempty(cfg.axis.YTick), ax.YTick = cfg.axis.YTick; end
if ~isempty(cfg.axis.YTickLabel), ax.YTickLabel = cfg.axis.YTickLabel; end

lgd = legend(ax, hList, leg, ...
    'Location', cfg.legend.location, ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.legend.fontSize);
lgd.Box = cfg.legend.box;

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('[SAVE] %s\n', pdfFile);
end

function rr = load_run_data(runMat)
% Load one saved DG-PW run.
S = load(runMat);
assert(isfield(S, 'run'), 'Run file must contain a run struct: %s', runMat);
R = S.run;

rr.uh = R.uh;
rr.k_pw = R.k_pw;
rr.nNURBS = double(R.n_dofs_nurbs);
rr.nurbs_refine = R.nurbs_refine;
rr.h = R.meta.h;
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

function errDG = compute_DG_error( ...
nurbs_curr, uI_curr, uA_curr, k_curr, ...
    nurbs_ref,  uI_ref,  uA_ref,  k_ref, ...
    L, a, dx_in, dx_out, chunkSize, sigma)

% Compute the DG eigenfunction error.
[xi, yi, wA_in] = grid_points_square(-a, a, dx_in);
[vI_curr, gxI_curr, gyI_curr] = iga_eval_val_grad(nurbs_curr, uI_curr, xi, yi, a);
[vI_ref,  gxI_ref,  gyI_ref ] = iga_eval_val_grad(nurbs_ref,  uI_ref,  xi, yi, a);

de  = vI_curr - vI_ref;
dex = gxI_curr - gxI_ref;
dey = gyI_curr - gyI_ref;
H1_in = sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2) * wA_in;

[xo, yo, wA_out] = grid_points_outer(L, a, dx_out);
H1_out = 0;
nPts = numel(xo);

k = 1;
while k <= nPts
    k2 = min(nPts, k + chunkSize - 1);
    X  = xo(k:k2);
    Y  = yo(k:k2);

    [vA_curr, gxA_curr, gyA_curr] = pw_eval_val_grad(uA_curr, k_curr, X, Y, L);
    [vA_ref,  gxA_ref,  gyA_ref ] = pw_eval_val_grad(uA_ref,  k_ref,  X, Y, L);

    de  = vA_curr - vA_ref;
    dex = gxA_curr - gxA_ref;
    dey = gyA_curr - gyA_ref;
    H1_out = H1_out + sum(abs(de).^2 + abs(dex).^2 + abs(dey).^2);
    k = k2 + 1;
end
H1_out = H1_out * wA_out;

[xg, yg, wL] = boundary_points_square(a, dx_in);
vA_curr = pw_eval_val(uA_curr, k_curr, xg, yg, L);
vA_ref  = pw_eval_val(uA_ref,  k_ref,  xg, yg, L);
vI_curr = iga_eval_val(nurbs_curr, uI_curr, xg, yg, a);
vI_ref  = iga_eval_val(nurbs_ref,  uI_ref,  xg, yg, a);

jump = (vA_curr - vA_ref) - (vI_curr - vI_ref);
J2 = sum(abs(jump).^2) * wL;

errDG = sqrt(H1_in + H1_out + sigma * J2);
end

function val = pw_eval_val(coeff, p_vec, X, Y, L)
% Evaluate the field value.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));
val = (coeff.' * expo) / L;
val = val(:);
end

function [val, gx, gy] = pw_eval_val_grad(coeff, p_vec, X, Y, L)
% Evaluate the field value and gradient.
F = [X(:)'; Y(:)'];
expo = exp((1i * 2*pi / L) * (p_vec * F));

val = (coeff.' * expo) / L;
fac = (1i * 2*pi / L) / L;
gx  = ((coeff .* p_vec(:,1)).' * expo) * fac;
gy  = ((coeff .* p_vec(:,2)).' * expo) * fac;

val = val(:);
gx  = gx(:);
gy  = gy(:);
end

function val = iga_eval_val(nurbs, coeff, X, Y, a)
% Evaluate the field value.
[val, ~, ~] = iga_eval_val_grad(nurbs, coeff, X, Y, a);
end

function [val, gx, gy] = iga_eval_val_grad(nurbs, coeff, X, Y, a)
% Evaluate the field value and gradient.
pu = nurbs.pu;
pv = nurbs.pv;
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

    s  = 0; su = 0; sv = 0;
    for j = (spanV-pv):spanV
        lv = j - (spanV-pv) + 1;
        for i = (spanU-pu):spanU
            lu = i - (spanU-pu) + 1;
            row = i + (j-1) * mU;
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
% Evaluate basis values and first derivatives.
ndu = zeros(p+1, p+1);
left = zeros(1, p+1);
right = zeros(1, p+1);

ndu(1,1) = 1.0;
for j = 1:p
    left(j+1)  = u - U(span+1-j);
    right(j+1) = U(span+j) - u;
    saved = 0.0;
    for r = 0:j-1
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
    pk = p - 1;
    if r >= 1
        d = d + ndu(r, pk+1) / ndu(pk+2, r);
    end
    if r <= p-1
        d = d - ndu(r+1, pk+1) / ndu(pk+2, r+1);
    end
    ders1(r+1) = d;
end
dN = ders1 * p;
end

function span = findspan_local(n, p, u, U)
% Find the active knot span for a parameter value.
if u >= U(n+2)
    span = n + 1;
    return;
end
if u <= U(p+1)
    span = p + 1;
    return;
end

low = p + 1;
high = n + 2;
mid = floor((low + high) / 2);

while (u < U(mid) || u >= U(mid+1))
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end

span = mid;
end

function [X, Y, wA] = grid_points_square(xmin, xmax, dx)
% Build midpoint quadrature points in the square.
x = xmin + dx/2 : dx : xmax - dx/2;
[Xg, Yg] = meshgrid(x, x);
X = Xg(:);
Y = Yg(:);
wA = dx * dx;
end

function [X, Y, wA] = grid_points_outer(L, a, dx)
% Build midpoint quadrature points outside the square.
x = -L/2 + dx/2 : dx : L/2 - dx/2;
[Xg, Yg] = meshgrid(x, x);
mask = ~(Xg >= -a & Xg <= a & Yg >= -a & Yg <= a);
X = Xg(mask);
Y = Yg(mask);
wA = dx * dx;
end

function [X, Y, wL] = boundary_points_square(a, ds)
% Build midpoint quadrature points on the square boundary.
t = -a + ds/2 : ds : a - ds/2;

xb = t;  yb = -a*ones(size(t));
xt = t;  yt =  a*ones(size(t));
xl = -a*ones(size(t)); yl = t;
xr =  a*ones(size(t)); yr = t;

X = [xb, xt, xl, xr].';
Y = [yb, yt, yl, yr].';
wL = ds;
end
