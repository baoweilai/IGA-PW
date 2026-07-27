function out = plot_cutoff_dg()
% Plot cutoff DG.

clc; close all; format short e;

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

%% ---------------- user params ----------------
Example    = 'Example_3';
refine     = 6;
pdeg       = 1;

Nc_select  = 2:8;   % [] means use all Nc

ref_refine = 7;
ref_p      = 2;
ref_Nc     = 40;

L          = 4;
a          = 0.2;
beta       = 20;

dx_in      = 5e-3;
dx_out     = 5e-3;
chunkSize  = 20000;

%% ---------------- plotting params ----------------
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
cfg.legend.location = 'northeast';

cfg.lineColors = [ ...
    223 122 094;
    060 064 091] / 255;

cfg.lineWidth  = 1.8;
cfg.markerSize = 8;

xPad = 0.10;
yPad = 10;

%% ---------------- paths ----------------
resultRoot = fullfile(pwd, 'result', Example);
if ~exist(resultRoot, 'dir')
    error('Directory not found: %s', resultRoot);
end

plotDir = fullfile(resultRoot, 'plots_lambda_DG_Nc', sprintf('refine_%02d', refine));
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

caseDir = fullfile(resultRoot, sprintf('refine_%02d', refine), sprintf('p_%d', pdeg));
if ~exist(caseDir, 'dir')
    error('Case directory not found: %s', caseDir);
end

%% ---------------- reference for eigenvalue and DG error ----------------
refRunMat = fullfile(resultRoot, ...
    sprintf('refine_%02d', ref_refine), ...
    sprintf('p_%d', ref_p), ...
    sprintf('Nc_%02d', ref_Nc), ...
    'run.mat');

if ~exist(refRunMat, 'file')
    error('Reference run.mat not found: %s', refRunMat);
end

ref = load_run_data(refRunMat);
lambda_ref = ref.lambda(1);
ref.u  = ref.uh(:,1);
ref.uI = ref.u(1:ref.nNURBS);
ref.uA = ref.u(ref.nNURBS+1 : ref.nNURBS + size(ref.k_pw,1));

fprintf('[REF] %s\n', refRunMat);
fprintf('[REF] lambda_ref = %.12f\n', lambda_ref);

%% ---------------- collect eigenvalue error ----------------
csvFile = fullfile(caseDir, 'summary.csv');
if ~exist(csvFile, 'file')
    error('Missing file: %s', csvFile);
end

Tlam = readtable(csvFile);

Nc_lam  = Tlam.Nc(:);
lam     = Tlam.lambda1(:);

valid = isfinite(Nc_lam) & isfinite(lam);
Nc_lam = Nc_lam(valid);
lam    = lam(valid);

[Nc_lam, idx] = sort(Nc_lam, 'ascend');
lam = lam(idx);

if ~isempty(Nc_select)
    keep = ismember(Nc_lam, Nc_select);
    Nc_lam = Nc_lam(keep);
    lam    = lam(keep);
end

errLambda = abs(lam - lambda_ref);
keep = isfinite(errLambda) & errLambda > 0;
Nc_lam = Nc_lam(keep);
errLambda = errLambda(keep);

if isempty(Nc_lam)
    error('No valid eigenvalue error data found.');
end

%% ---------------- collect DG error with cache ----------------
Nc_list = scan_existing_Nc(caseDir);
if ~isempty(Nc_select)
    Nc_list = Nc_list(ismember(Nc_list, Nc_select));
end
if isempty(Nc_list)
    error('No valid Nc results found in %s', caseDir);
end

baseName = 'cutoff_dg';
errCsv   = fullfile(plotDir, [baseName, '.csv']);
errMat   = fullfile(plotDir, [baseName, '.mat']);
pdfFile  = fullfile(plotDir, [baseName, '.pdf']);

if exist(errCsv, 'file')
    T_old = readtable(errCsv);
    fprintf('[CACHE] Found %s\n', errCsv);
    cacheVars = T_old.Properties.VariableNames;
    hasRefColumns = all(ismember({'ref_refine','ref_p','ref_Nc','lambda_ref'}, cacheVars));
    if ~hasRefColumns ...
            || any(T_old.ref_refine ~= ref_refine) ...
            || any(T_old.ref_p ~= ref_p) ...
            || any(T_old.ref_Nc ~= ref_Nc) ...
            || any(abs(T_old.lambda_ref - lambda_ref) > 1e-10)
        fprintf('[CACHE] Ignoring stale DG cache for a different reference.\n');
        T_old = table();
    end
else
    T_old = table();
end

Nc_dg = [];
errDG = [];

for i = 1:numel(Nc_list)
    Nc = Nc_list(i);

    idx = [];
    if ~isempty(T_old) && ismember('Nc', T_old.Properties.VariableNames)
        idx = find(T_old.Nc == Nc, 1, 'last');
    end

    if ~isempty(idx) && ...
            ismember('u1_DG', T_old.Properties.VariableNames) && ...
            isfinite(T_old.u1_DG(idx))

        Nc_dg(end+1,1)  = Nc; %#ok<AGROW>
        errDG(end+1,1)  = T_old.u1_DG(idx); %#ok<AGROW>

        fprintf('[CACHE] p=%d, refine=%d, Nc=%d: DG=%.3e\n', ...
            pdeg, refine, Nc, errDG(end));
        continue;
    end

    runMat = fullfile(caseDir, sprintf('Nc_%02d', Nc), 'run.mat');
    assert(exist(runMat, 'file') == 2, 'Missing run data: %s', runMat);

    rr = load_run_data(runMat);
    u  = rr.uh(:,1);

    uI = u(1:rr.nNURBS);
    uA = u(rr.nNURBS+1 : rr.nNURBS + size(rr.k_pw,1));

    [tf, loc] = ismember(rr.k_pw, ref.k_pw, 'rows');
    if ~all(tf)
        error('Current PW modes are not contained in reference PW modes.');
    end

    uA_pad = zeros(size(ref.uA));
    uA_pad(loc(tf)) = uA(tf);

    alpha = ref.uA' * uA_pad;
    if abs(alpha) > 1e-14
        phase = exp(-1i * angle(alpha));
    else
        phase = 1;
    end

    uI = uI * phase;
    uA = uA * phase;

    assert(isfinite(rr.h) && rr.h > 0, 'Mesh size h must be positive.');
    sigma = beta * (Nc + 1 / rr.h);

    err = compute_DG_error( ...
        rr.nurbs_refine, uI, uA, rr.k_pw, ...
        ref.nurbs_refine, ref.uI, ref.uA, ref.k_pw, ...
        L, a, dx_in, dx_out, chunkSize, sigma);

    Nc_dg(end+1,1) = Nc; %#ok<AGROW>
    errDG(end+1,1) = err; %#ok<AGROW>

    fprintf('[CALC ] p=%d, refine=%d, Nc=%d: DG=%.3e\n', ...
        pdeg, refine, Nc, err);
end

if isempty(Nc_dg)
    error('No valid DG data found.');
end

[Nc_dg, ord] = sort(Nc_dg, 'ascend');
errDG = errDG(ord);

assert(all(errDG > 0) && all(errLambda > 0), 'Plot errors must be positive.');

%% ---------------- align common Nc ----------------
Nc_common = intersect(Nc_lam, Nc_dg, 'stable');
if isempty(Nc_common)
    error('No common Nc values between eigenvalue error and DG error.');
end

[~, ia] = ismember(Nc_common, Nc_lam);
[~, ib] = ismember(Nc_common, Nc_dg);

errLambda_plot = errLambda(ia);
errDG_plot     = errDG(ib);

T = table(Nc_common, errLambda_plot, errDG_plot, ...
    repmat(ref_refine, numel(Nc_common), 1), ...
    repmat(ref_p, numel(Nc_common), 1), ...
    repmat(ref_Nc, numel(Nc_common), 1), ...
    repmat(lambda_ref, numel(Nc_common), 1), ...
    'VariableNames', {'Nc','lambda1_err','u1_DG', ...
    'ref_refine','ref_p','ref_Nc','lambda_ref'});
writetable(T, errCsv);
save(errMat, 'T');

fprintf('[SAVE] %s\n', errCsv);
fprintf('[SAVE] %s\n', errMat);

%% ---------------- plot ----------------
fig = figure('Color', cfg.fig.bgColor, 'Units', 'inches', ...
    'Position', [1 1 cfg.fig.width cfg.fig.height], ...
    'Renderer', cfg.fig.renderer);
ax = axes(fig);

hLines = semilogy(ax, Nc_common, errLambda_plot, '-o', Nc_common, errDG_plot, '-s');
style_two_lines(hLines, cfg);

ax.Units = 'normalized';
ax.Position = [cfg.layout.left, cfg.layout.bottom, ...
    1 - cfg.layout.left - cfg.layout.right, ...
    1 - cfg.layout.bottom - cfg.layout.top];

set(ax, ...
    'FontSize', cfg.axes.fontSize, ...
    'LineWidth', cfg.axes.lineWidth, ...
    'TickDir', cfg.axes.tickDir, ...
    'XMinorTick', cfg.axes.xMinorTick, ...
    'YMinorTick', cfg.axes.yMinorTick, ...
    'Box', 'on');

grid(ax, 'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.XRuler.MinorTick = 'off';
ax.YRuler.MinorTick = 'off';

ax.XTick = Nc_common(:).';
ax.XLim = [Nc_common(1), Nc_common(end)];

xlabel(ax, '$K$', ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);
ylabel(ax, 'Error', ...
    'Interpreter', 'latex', ...
    'FontSize', cfg.axes.labelSize);

lgd = legend(ax, hLines, ...
    {'$|\lambda_1-\lambda_{1}^{\mathrm{DG}}|$', ...
    '$\|u_1-u_{1}^{\mathrm{DG}}\|_{\mathrm{DG}}$'}, ...
    'Location', cfg.legend.location, ...
    'Interpreter', 'latex');
lgd.Box = 'off';
lgd.FontSize = cfg.legend.fontSize;

apply_semilogy_padding(ax, xPad, yPad);

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('[SAVE] %s\n', pdfFile);

out = struct('Example', Example, 'refine', refine, 'pdeg', pdeg, ...
    'Nc', Nc_common, 'lambda_ref', lambda_ref, ...
    'errLambda', errLambda_plot, 'errDG', errDG_plot, ...
    'csv', errCsv, 'mat', errMat, 'pdf', pdfFile);

end

function style_two_lines(hLines, cfg)
% Apply the shared style to both plotted curves.
marks = {'o', 's'};
cols = cfg.lineColors(1:2, :);

for i = 1:numel(hLines)
    set(hLines(i), ...
        'LineWidth', cfg.lineWidth, ...
        'Color', cols(i, :), ...
        'Marker', marks{i}, ...
        'MarkerSize', cfg.markerSize, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cols(i, :));
end
end

function apply_semilogy_padding(ax, xPad, yPad)
% Apply semilogy padding.
x = ax.XLim;
dx = x(2) - x(1);
if dx <= 0
    dx = 1;
end
ax.XLim = [x(1) - xPad * dx, x(2) + xPad * dx];

h = findobj(ax, 'Type', 'line');
yAll = [];
for i = 1:numel(h)
    y = h(i).YData;
    y = y(isfinite(y) & y > 0);
    yAll = [yAll, y]; %#ok<AGROW>
end

assert(~isempty(yAll), 'No positive Y data available for axis limits.');

ax.YScale = 'log';
ax.YLim = [min(yAll) / yPad, max(yAll) * yPad];
end

function rr = load_run_data(runMat)
% Load one saved DG-PW run.
S = load(runMat);
assert(isfield(S, 'run'), 'Run file must contain a run struct: %s', runMat);
R = S.run;

rr.uh = R.uh;
rr.lambda = R.lambda(:).';
rr.k_pw = R.k_pw;
rr.nNURBS = double(R.n_dofs_nurbs);
rr.nurbs_refine = R.nurbs_refine;
rr.h = R.meta.h;
end

function Nc_list = scan_existing_Nc(caseDir)
% Read the available cutoff values from completed case folders.
dd = dir(fullfile(caseDir, 'Nc_*'));
Nc_list = [];

for i = 1:numel(dd)
    if ~dd(i).isdir
        continue;
    end
    tok = regexp(dd(i).name, '^Nc_(\d+)$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    runMat = fullfile(caseDir, dd(i).name, 'run.mat');
    if exist(runMat, 'file')
        Nc_list(end+1,1) = str2double(tok{1}); %#ok<AGROW>
    end
end

Nc_list = unique(sort(Nc_list));
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

    s  = 0;
    su = 0;
    sv = 0;

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
