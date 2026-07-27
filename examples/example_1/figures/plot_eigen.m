function plot_eigen()
% Plot the first four eigenfunction DG errors for p=1 and p=2.

% Initialize paths and fixed comparison parameters.
clc; close all; format short e;

figDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(figDir);
pList = [1 2];
rList = 2:6;
Nc = 30;
dx = 5e-3;

set(groot, 'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');
addpath(fullfile(exampleDir, 'model', 'h_convergence', 'nurbs'));

root = fullfile(exampleDir, 'data', 'result', 'Example_1');
out = fullfile(root, 'eigen', sprintf('Nc_%d', Nc));
if ~exist(out, 'dir'), mkdir(out); end

refFile = fullfile(root, 'Nc_45', 'p_3', 'refine_08', 'run.mat');
assert(exist(refFile, 'file') == 2, 'Missing reference run: %s', refFile);
% Sample the reference solution and its outer-domain Gram matrix.
ref = load_run(refFile);
q = quadrature(ref.meta.L, ref.meta.a, dx);
refSample = reference_samples(ref, q);
refOuter = pw_l2_gram(ref.pw, ref.k, ref.pw, ref.k, ...
    ref.meta.L, ref.meta.a, dx);

% Evaluate aligned DG errors over both polynomial degrees.
D = struct([]);
for ip = 1:numel(pList)
    p = pList(ip);
    D(ip).p = p;
    D(ip).r = rList(:);
    D(ip).h = zeros(numel(rList), 1);
    D(ip).e = zeros(numel(rList), 4);

    for ir = 1:numel(rList)
        r = rList(ir);
        file = fullfile(root, sprintf('Nc_%d', Nc), sprintf('p_%d', p), ...
            sprintf('refine_%02d', r), 'run.mat');
        assert(exist(file, 'file') == 2, 'Missing run: %s', file);
        cur = load_run(file);

        sigma34 = cur.meta.beta * (cur.meta.Nc + 1 / cur.meta.h);
        assert(abs(sigma34 - cur.meta.sigma) <= 1e-12 * cur.meta.sigma, ...
            'Inconsistent penalty in %s.', file);
        sigma12 = 10 * (cur.meta.Nc + 1 / cur.meta.h);

        [G, curSample] = build_grams(cur, ref, refOuter, q, dx);
        A = align_modes(cur, ref, G);
        D(ip).h(ir) = cur.meta.h;
        D(ip).e(ir, :) = dg_errors(cur, ref, refSample, curSample, ...
            q, dx, [sigma12 sigma12 sigma34 sigma34], A);
        fprintf('[p=%d r=%d] h=%.4e DG=%s\n', ...
            p, r, cur.meta.h, mat2str(D(ip).e(ir, :), 6));
    end

    [D(ip).h, order] = sort(D(ip).h);
    D(ip).r = D(ip).r(order);
    D(ip).e = D(ip).e(order, :);
    T = table(D(ip).r, D(ip).h, D(ip).e(:, 1), D(ip).e(:, 2), ...
        D(ip).e(:, 3), D(ip).e(:, 4), ...
        'VariableNames', {'refine', 'h', 'u1_DG', 'u2_DG', 'u3_DG', 'u4_DG'});
    writetable(T, fullfile(out, sprintf('p%d.csv', p)));
end

% Export one convergence figure for each polynomial degree.
draw_plots(D, out, plot_style());
fprintf('[DONE] %s\n', out);
end

function rr = load_run(file)
% Load the four saved eigenvectors.
S = load(file, 'run');
run = S.run;
nNURBS = double(run.n_dofs_nurbs);
rr.iga = run.uh(1:nNURBS, 1:4);
rr.pw = run.uh(nNURBS + 1:end, 1:4);
rr.k = run.k_pw;
rr.nurbs_original = run.nurbs_original;
rr.nurbs = run.nurbs_refine;
rr.meta = run.meta;
end

function [grams, sample] = build_grams(cur, ref, refOuterL2, quad, dx)
% Build the L2 Gram matrices and aligned-error samples.
[MrrIn, MhhIn, MrhIn] = iga_l2_grams(cur, ref);
MhhOut = pw_l2_gram( ...
    cur.pw, cur.k, cur.pw, cur.k, cur.meta.L, cur.meta.a, dx);
MrhOut = pw_l2_gram( ...
    ref.pw, ref.k, cur.pw, cur.k, cur.meta.L, cur.meta.a, dx);

[vi, gxi, gyi] = iga_eval_physical(cur.nurbs, cur.iga, quad.x_in, quad.y_in, cur.meta.a);
sample.vi = vi;
sample.gxi = gxi;
sample.gyi = gyi;

[vib, ~, ~] = iga_eval_physical(cur.nurbs, cur.iga, quad.x_bnd, quad.y_bnd, cur.meta.a);
[vob, ~, ~] = pw_eval(cur.pw, cur.k, quad.x_bnd, quad.y_bnd, cur.meta.L);
sample.jump = vob - vib;

grams.Mrr = hermitian(MrrIn + refOuterL2);
grams.Mhh = hermitian(MhhIn + MhhOut);
grams.Mrh = MrhIn + MrhOut;
end

function align = align_modes(cur, ref, grams)
% Simple modes use phase alignment, while modes 3 and 4 use L2-Procrustes alignment.
nEig = size(cur.iga, 2);
assert(nEig == 4, 'The integrated plot requires four eigenfunctions.');

align.ref = eye(nEig);
align.cur = eye(nEig);
for j = 1:2
    align.cur(j, j) = scalar_pw_phase(cur, ref, j);
end

block = 3:4;
blockGrams.Mrr = grams.Mrr(block, block);
blockGrams.Mhh = grams.Mhh(block, block);
blockGrams.Mrh = grams.Mrh(block, block);
blockAlign = l2_alignment(blockGrams);
align.ref(block, block) = eye(2) / blockAlign.Rr;
align.cur(block, block) = (eye(2) / blockAlign.Rh) * blockAlign.Q;
end

function phase = scalar_pw_phase(cur, ref, eigIndex)
% Match the saved plane-wave phases for simple modes.
[isPresent, refIndex] = ismember(cur.k, ref.k, 'rows');
assert(all(isPresent), 'The reference plane-wave basis does not contain the current basis.');
alpha = ref.pw(refIndex, eigIndex)' * cur.pw(:, eigIndex);
if abs(alpha) > 1e-14
    phase = exp(-1i * angle(alpha));
else
    phase = 1;
end
end

function align = l2_alignment(grams)
% Compute the L2-Procrustes rotation from normalized Gram matrices.
Rr = chol(hermitian(grams.Mrr));
Rh = chol(hermitian(grams.Mhh));
C = Rr' \ (grams.Mrh / Rh);
[L, ~, R] = svd(C, 'econ');
align.Rr = Rr;
align.Rh = Rh;
align.Q = R * L';
end

function err = dg_errors(cur, ref, refSample, curSample, quad, dx, sigma, align)
% Integrate aligned volume and jump residuals.
Hout = pw_outer_h1(cur, ref, dx, align);
vi = refSample.vi * align.ref - curSample.vi * align.cur;
gxi = refSample.gxi * align.ref - curSample.gxi * align.cur;
gyi = refSample.gyi * align.ref - curSample.gyi * align.cur;
jump = refSample.jump * align.ref - curSample.jump * align.cur;
EDG = hermitian(h1_gram(vi, gxi, gyi, quad.w_in) + Hout);
jumpEnergy = real(diag(jump' * jump * quad.w_bnd))';
err = sqrt(max(real(diag(EDG))' + sigma .* jumpEnergy, 0));
end

function H = pw_outer_h1(cur, ref, dx, align)
% Integrate the aligned outer plane-wave H1 residual.
coeff = [ref.pw * align.ref; -cur.pw * align.cur];
k = [ref.k; cur.k];
x = -ref.meta.L / 2 + dx / 2:dx:ref.meta.L / 2 - dx / 2;
[val, gx, gy] = pw_eval_tensor(coeff, k, x, x, ref.meta.L);
[X, Y] = ndgrid(x, x);
outer = ~(X >= -ref.meta.a & X <= ref.meta.a & ...
    Y >= -ref.meta.a & Y <= ref.meta.a);
nEig = size(coeff, 2);
val = reshape(val, [], nEig);
gx = reshape(gx, [], nEig);
gy = reshape(gy, [], nEig);
outer = outer(:);
val = val(outer, :);
gx = gx(outer, :);
gy = gy(outer, :);
H = (val' * val + gx' * gx + gy' * gy) * dx^2;
end

function [val, gx, gy] = pw_eval_tensor(coeff, k, x, y, L)
% Evaluate a plane-wave expansion on a tensor grid without coefficient padding.
alpha = 1i * 2 * pi / L;
kx = unique(k(:, 1));
Ex = exp(alpha * (x(:) * kx'));
nEig = size(coeff, 2);
F = zeros(numel(kx), numel(y), nEig);
Fy = zeros(numel(kx), numel(y), nEig);

for i = 1:numel(kx)
    idx = k(:, 1) == kx(i);
    phaseY = exp(alpha * (k(idx, 2) * y(:)'));
    for j = 1:nEig
        F(i, :, j) = sum(coeff(idx, j) .* phaseY, 1) / L;
        Fy(i, :, j) = sum((alpha * k(idx, 2) .* coeff(idx, j)) .* phaseY, 1) / L;
    end
end

val = zeros(numel(x), numel(y), nEig);
gx = zeros(numel(x), numel(y), nEig);
gy = zeros(numel(x), numel(y), nEig);
for j = 1:nEig
    val(:, :, j) = Ex * F(:, :, j);
    gx(:, :, j) = Ex * ((alpha * kx) .* F(:, :, j));
    gy(:, :, j) = Ex * Fy(:, :, j);
end
end

function [Mrr, Mhh, Mrh] = iga_l2_grams(cur, ref)
% Integrate inner-domain L2 Gram matrices on the current IGA element grid.
nurbs = cur.nurbs;
orig = cur.nurbs_original;
[gp_u, gw_u] = grule(nurbs.pu + 5);
[gp_v, gw_v] = grule(nurbs.pv + 5);
nEig = size(cur.iga, 2);
Mrr = zeros(nEig);
Mhh = zeros(nEig);
Mrh = zeros(nEig);

for e = 1:nurbs.NoEs
    ue = nurbs.Coordinate(e, 1:2);
    ve = nurbs.Coordinate(e, 3:4);
    uJ = (ue(2) - ue(1)) / 2;
    vJ = (ve(2) - ve(1)) / 2;
    for j = 1:numel(gp_v)
        v = vJ * gp_v(j) + (ve(1) + ve(2)) / 2;
        for i = 1:numel(gp_u)
            u = uJ * gp_u(i) + (ue(1) + ue(2)) / 2;
            [~, ~, ~, ~, DF] = NurbsSurfaceDers( ...
                orig.ConPts, orig.knotU, orig.knotV, orig.weights, ...
                orig.pu, u, orig.pv, v);
            w = abs(det(DF)) * gw_u(i) * uJ * gw_v(j) * vJ;
            vh = iga_eval_parametric(cur.nurbs, cur.iga, u, v);
            vr = iga_eval_parametric(ref.nurbs, ref.iga, u, v);
            Mrr = Mrr + vr' * vr * w;
            Mhh = Mhh + vh' * vh * w;
            Mrh = Mrh + vr' * vh * w;
        end
    end
end
end

function sample = reference_samples(ref, quad)
% Build reference volume fields and actual interface jumps.
[sample.vi, sample.gxi, sample.gyi] = iga_eval_physical( ...
    ref.nurbs, ref.iga, quad.x_in, quad.y_in, ref.meta.a);
[vib, ~, ~] = iga_eval_physical( ...
    ref.nurbs, ref.iga, quad.x_bnd, quad.y_bnd, ref.meta.a);
[vob, ~, ~] = pw_eval(ref.pw, ref.k, ...
    quad.x_bnd, quad.y_bnd, ref.meta.L);
sample.jump = vob - vib;
end

function M = pw_l2_gram(Ca, ka, Cb, kb, L, a, dx)
% Integrate the outer-domain plane-wave L2 Gram matrix.
qN = max(abs(ka(:))) + max(abs(kb(:)));
Mker = outer_midpoint_kernel(L, a, dx, qN);
M = zeros(size(Ca, 2), size(Cb, 2));

for i = 1:size(ka, 1)
    qx = kb(:, 1) - ka(i, 1);
    qy = kb(:, 2) - ka(i, 2);
    idx = sub2ind(size(Mker), qx + qN + 1, qy + qN + 1);
    rowM = reshape(Mker(idx), 1, []);
    M = M + Ca(i, :)' * (rowM * Cb);
end
end

function Mker = outer_midpoint_kernel(L, a, dx, qN)
% Build the separable Fourier kernel of the existing outer midpoint quadrature.
q = (-qN:qN)';
x = -a + dx / 2:dx:a - dx / 2;
F = dx * sum(exp((1i * 2 * pi / L) * (q * x)), 2);
Fy = conj(F);
Mker = -(F * Fy') / L^2;
Mker(qN + 1, qN + 1) = Mker(qN + 1, qN + 1) + 1;
end

function G = h1_gram(v, gx, gy, weight)
% Form the H1 Gram matrix from the sampled fields.
G = (v' * v + gx' * gx + gy' * gy) * weight;
end

function [val, gx, gy] = iga_eval_physical(nurbs, coeff, X, Y, a)
% Evaluate complex IGA values and physical gradients at common points.
nPts = numel(X);
nEig = size(coeff, 2);
val = zeros(nPts, nEig);
gx = zeros(nPts, nEig);
gy = zeros(nPts, nEig);
for k = 1:nPts
    u = max(0, min(1, (X(k) + a) / (2 * a)));
    v = max(0, min(1, (Y(k) + a) / (2 * a)));
    [val(k, :), du, dv] = iga_eval_parametric(nurbs, coeff, u, v);
    gx(k, :) = du / (2 * a);
    gy(k, :) = dv / (2 * a);
end
end

function [val, du, dv] = iga_eval_parametric(nurbs, coeff, u, v)
% Evaluate complex IGA values and parametric gradients at one point.
i = findspan(nurbs.Ubar, nurbs.pu, u);
j = findspan(nurbs.Vbar, nurbs.pv, v);
iu = i - nurbs.pu:i;
jv = j - nurbs.pv:j;
row = zeros((nurbs.pu + 1) * (nurbs.pv + 1), 1);
q = 1;
for jj = jv
    for ii = iu
        row(q) = ii + (jj - 1) * nurbs.m;
        q = q + 1;
    end
end

Nu = bspbasisDers(nurbs.Ubar, nurbs.pu, u, 1);
Nv = bspbasisDers(nurbs.Vbar, nurbs.pv, v, 1);
B = Nu(1, :)' * Nv(1, :);
Bu = Nu(2, :)' * Nv(1, :);
Bv = Nu(1, :)' * Nv(2, :);
val = sum(coeff(row, :) .* B(:), 1);
du = sum(coeff(row, :) .* Bu(:), 1);
dv = sum(coeff(row, :) .* Bv(:), 1);
end

function [val, gx, gy] = pw_eval(coeff, k, X, Y, L)
% Evaluate complex PW values and gradients without cross-space coefficient padding.
F = [X(:)'; Y(:)'];
E = exp((1i * 2 * pi / L) * (k * F));
alpha = 1i * 2 * pi / L;
nPts = numel(X);
nEig = size(coeff, 2);
val = zeros(nPts, nEig);
gx = zeros(nPts, nEig);
gy = zeros(nPts, nEig);
for j = 1:nEig
    s0 = sum(coeff(:, j) .* E, 1) / L;
    sx = sum((coeff(:, j) .* k(:, 1)) .* E, 1) * alpha / L;
    sy = sum((coeff(:, j) .* k(:, 2)) .* E, 1) * alpha / L;
    val(:, j) = reshape(s0, [], 1);
    gx(:, j) = reshape(sx, [], 1);
    gy(:, j) = reshape(sy, [], 1);
end
end

function quad = quadrature(L, a, dx)
% Build the existing midpoint volume and interface rules.
x = -a + dx / 2:dx:a - dx / 2;
[X, Y] = meshgrid(x, x);
quad.x_in = X(:);
quad.y_in = Y(:);
quad.w_in = dx^2;

t = x;
quad.x_bnd = [t, t, -a * ones(size(t)), a * ones(size(t))]';
quad.y_bnd = [-a * ones(size(t)), a * ones(size(t)), t, t]';
quad.w_bnd = dx;
quad.L = L;
end

function A = hermitian(A)
% Apply Hermitian symmetrization before factorizations and eigensolves.
A = (A + A') / 2;
end

function style = plot_style()
% Return the separate-p style used by plot_h_convergence.
style.fig.width = 4.8;
style.fig.height = 3.0;
style.fig.renderer = 'painters';
style.fig.bgColor = 'w';
style.layout.left = 0.14;
style.layout.right = 0.04;
style.layout.bottom = 0.16;
style.layout.top = 0.08;
style.axes.fontSize = 11;
style.axes.labelSize = 13;
style.axes.lineWidth = 1.0;
style.axes.tickDir = 'out';
style.axes.padX = 0.4;
style.axes.padYLow = 2;
style.axes.padYHigh = 1;
style.curve.colors = [ ...
    223 122 094;
    060 064 091;
    130 178 154;
    242 204 142] / 255;
style.curve.slopeColor = [033 158 188] / 255;
style.curve.extraSlopeColor = [239 065 067] / 255;
style.curve.markers = {'o', 's', '^', 'd'};
style.curve.lineWidth = 1.8;
style.curve.markerSize = 8;
style.curve.markerFaceColor = 'w';
style.slope.lineWidth = 1.8;
style.legend.fontSize = 10;
style.legend.xL1 = 0.40;
style.legend.xL2 = 0.50;
style.legend.xLT = 0.52;
style.legend.xR1 = 0.70;
style.legend.xR2 = 0.88;
style.legend.xRT = 0.90;
style.legend.rowY = [0.4 0.3 0.2 0.1];
end

function draw_plots(D, outDir, style)
% Plot four eigenfunction DG errors in one figure for each polynomial degree.
assert(numel(D) == 2, 'The integrated plot requires p=1 and p=2.');

% Create and validate one panel for each polynomial degree.
for ip = 1:numel(D)
    h = D(ip).h;
    Y = D(ip).e;
    p = D(ip).p;
    assert(size(Y, 2) == 4 && all(Y(:) > 0), ...
        'Each p figure requires four positive eigenfunction errors.');

    % Create the axes and draw the four computed errors.
    fig = figure('Color', style.fig.bgColor, ...
        'Units', 'inches', 'Position', [1 1 style.fig.width style.fig.height], ...
        'Renderer', style.fig.renderer);
    ax = axes(fig, 'Position', [style.layout.left, style.layout.bottom, ...
        1 - style.layout.left - style.layout.right, ...
        1 - style.layout.bottom - style.layout.top]);
    hold(ax, 'on');

    hCurve = gobjects(4, 1);
    for j = 1:4
        hCurve(j) = plot(ax, h, Y(:, j), '-', ...
            'LineWidth', style.curve.lineWidth, ...
            'Color', style.curve.colors(j, :), ...
            'Marker', style.curve.markers{j}, ...
            'MarkerSize', style.curve.markerSize, ...
            'MarkerFaceColor', style.curve.markerFaceColor, ...
        'MarkerEdgeColor', style.curve.colors(j, :));
    end

    % Construct the reference slopes for the selected polynomial degree.
    anchor = max(2, min(numel(h) - 1, round(numel(h) / 2)));
    if p == 1
        lowerCoefficient = max(Y(:, 4) ./ h);
        upperCoefficient = min(Y(:, 2) ./ h);
        assert(lowerCoefficient < upperCoefficient, ...
            'The slope-1 guide cannot be placed between i=2 and i=4.');
        slopeCoefficient = sqrt(lowerCoefficient * upperCoefficient);
        slopeLine1 = slopeCoefficient * h;
        hSlope1 = plot(ax, h, slopeLine1, '--', ...
            'LineWidth', style.slope.lineWidth, 'Color', style.curve.slopeColor);
        hSlope2 = gobjects(0);
        slopeLines = slopeLine1;
    else
        upperAtAnchor = Y(anchor, 2);
        lowerAtAnchor = Y(anchor, 3);
        slope1AtAnchor = upperAtAnchor^(11 / 15) * lowerAtAnchor^(4 / 15);
        slope2AtAnchor = upperAtAnchor^(4 / 15) * lowerAtAnchor^(11 / 15);
        slopeLine1 = slope1AtAnchor * (h / h(anchor));
        slopeLine2 = slope2AtAnchor * (h / h(anchor)).^2;
        assert(all(slopeLine1 < Y(:, 2) & slopeLine1 > Y(:, 3)) && ...
            all(slopeLine2 < Y(:, 2) & slopeLine2 > Y(:, 3)), ...
            'Both reference guides must remain between i=2 and i=3.');
        hSlope1 = plot(ax, h, slopeLine1, '--', ...
            'LineWidth', style.slope.lineWidth, ...
            'Color', style.curve.slopeColor);
        hSlope2 = plot(ax, h, slopeLine2, '--', ...
            'LineWidth', style.slope.lineWidth, ...
            'Color', style.curve.extraSlopeColor);
        slopeLines = [slopeLine1(:); slopeLine2(:)];
    end

    % Format, export, and close the completed figure.
    set_axes(ax, style);
    xlabel(ax, '$h$', 'Interpreter', 'latex', 'FontSize', style.axes.labelSize);
    ylabel(ax, '$\|u_i-u_{i}^{\mathrm{DG}}\|_{\mathrm{DG}}$', ...
        'Interpreter', 'latex', 'FontSize', style.axes.labelSize);
    set_limits(ax, h, [Y(:); slopeLines(:)], style);
    draw_legend(ax, hSlope1, hSlope2, hCurve, p, style);

    base = fullfile(outDir, sprintf('p%d', p));
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    close(fig);
    fprintf('[SAVE] %s\n', base);
end
end

function set_axes(ax, style)
% Apply the Example 1 log-log axes style.
set(ax, 'XScale', 'log', 'YScale', 'log', ...
    'FontSize', style.axes.fontSize, ...
    'LineWidth', style.axes.lineWidth, 'TickDir', style.axes.tickDir, ...
    'Box', 'on', 'XMinorTick', 'off', 'YMinorTick', 'off');
grid(ax, 'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
end

function set_limits(ax, x, y, style)
% Apply the logarithmic limits and ticks.
y = y(y > 0);
yLim = [min(y) / (1 + style.axes.padYLow), ...
    max(y) * (1 + style.axes.padYHigh)];
tickPower = floor(log10(yLim(1))):ceil(log10(yLim(2)));
set(ax, 'XLim', [min(x) / (1 + style.axes.padX), max(x) * (1 + style.axes.padX)], ...
    'YLim', yLim, 'YTick', 10.^tickPower);
end

function draw_legend(ax, hSlope1, hSlope2, hCurve, p, style)
% The legend geometry and size match plot_h_convergence.
for j = 1:4
    curveLabel = sprintf('$i=%d$', j);
    legend_entry(ax, hCurve(j), style.legend.xR1, style.legend.xR2, ...
        style.legend.xRT, style.legend.rowY(j), ...
        curveLabel, style.legend.fontSize);
end
if p == 1
    legend_entry(ax, hSlope1, style.legend.xL1, style.legend.xL2, ...
        style.legend.xLT, style.legend.rowY(4), ...
        '$\mathrm{Slope} = 1$', style.legend.fontSize);
else
    legend_entry(ax, hSlope1, style.legend.xL1, style.legend.xL2, ...
        style.legend.xLT, style.legend.rowY(3), ...
        '$\mathrm{Slope} = 1$', style.legend.fontSize);
    legend_entry(ax, hSlope2, style.legend.xL1, style.legend.xL2, ...
        style.legend.xLT, style.legend.rowY(4), ...
        '$\mathrm{Slope} = 2$', style.legend.fontSize);
end
end

function legend_entry(ax, hLine, x1n, x2n, xtn, yn, label, fontSize)
% Draw one legend entry in normalized axes coordinates.
[x1, y] = axes_point(ax, x1n, yn);
[x2, ~] = axes_point(ax, x2n, yn);
[xt, ~] = axes_point(ax, xtn, yn);
c = get(hLine, 'Color');
ls = get(hLine, 'LineStyle');
lw = get(hLine, 'LineWidth');
mk = get(hLine, 'Marker');
ms = get(hLine, 'MarkerSize');
mfc = get(hLine, 'MarkerFaceColor');
mec = get(hLine, 'MarkerEdgeColor');

line(ax, [x1 x2], [y y], 'Color', c, 'LineStyle', ls, ...
    'LineWidth', lw, 'Clipping', 'off');
line(ax, 0.49 * (x1 + x2), y, 'Color', c, 'LineStyle', 'none', ...
    'Marker', mk, 'MarkerSize', ms, 'MarkerFaceColor', mfc, ...
    'MarkerEdgeColor', mec, 'LineWidth', lw, 'Clipping', 'off');
text(ax, xt, y, label, 'Interpreter', 'latex', 'FontSize', fontSize, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', 'Margin', 0.5, 'Clipping', 'off');
end

function [x, y] = axes_point(ax, xn, yn)
% Convert normalized axes coordinates to logarithmic data coordinates.
lx = log10(ax.XLim);
ly = log10(ax.YLim);
x = 10^(lx(1) + xn * (lx(2) - lx(1)));
y = 10^(ly(1) + yn * (ly(2) - ly(1)));
end
