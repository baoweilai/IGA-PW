function [errDG, sigma] = cutoff_eigenfunction_errors(runFiles, refFile, eigBlocks, dxIn, dxOut)
% Compute aligned DG errors for cutoff-convergence eigenfunctions.

eigCols = unique([eigBlocks{:}], 'stable');
ref = load_eigen_run(refFile, eigCols);
quad = midpoint_quadrature(ref.meta.L, ref.meta.a, dxIn);
refSample = reference_samples(ref, quad, dxOut);

nP = size(runFiles, 1);
nNc = size(runFiles, 2);
errDG = zeros(nP, nNc, numel(eigCols));
sigma = zeros(nP, nNc);

for ip = 1:nP
    for iNc = 1:nNc
        runFile = runFiles{ip, iNc};
        cur = load_eigen_run(runFile, eigCols);
        assert(cur.meta.L == ref.meta.L && cur.meta.a == ref.meta.a, ...
            'Current and reference domains do not match: %s', runFile);

        sigma(ip, iNc) = cur.meta.beta * (1 / cur.meta.h + ...
            cur.meta.Nc * 2 * pi / cur.meta.L);
        assert(abs(sigma(ip, iNc) - cur.meta.sigma) <= 1e-12 * cur.meta.sigma, ...
            'Stored penalty coefficient is inconsistent in %s.', runFile);

        [grams, curSample] = build_samples(cur, ref, refSample, quad, dxOut);
        for ib = 1:numel(eigBlocks)
            [found, block] = ismember(eigBlocks{ib}, eigCols);
            assert(all(found), 'An eigenspace block is absent from eigCols.');

            blockGrams.Mrr = grams.Mrr(block, block);
            blockGrams.Mhh = grams.Mhh(block, block);
            blockGrams.Mrh = grams.Mrh(block, block);
            align = l2_alignment(blockGrams);
            blockError = aligned_dg_error( ...
                refSample, curSample, quad, dxOut, sigma(ip, iNc), align, block);
            errDG(ip, iNc, block) = reshape(blockError, 1, 1, []);
        end

        fprintf('[EIGFUN] p=%d Nc=%d DG=%s\n', cur.meta.pu, cur.meta.Nc, ...
            mat2str(reshape(errDG(ip, iNc, :), 1, []), 8));
    end
end
end

function rr = load_eigen_run(runFile, eigCols)
% Load only the fields needed for eigenfunction postprocessing.
S = load(runFile, 'run');
run = S.run;
required = {'uh','k_pw','n_dofs_nurbs','nurbs_original','nurbs_refine','meta'};
assert(all(isfield(run, required)), 'Incomplete run structure in %s.', runFile);
assert(size(run.uh, 2) >= max(eigCols), 'Missing eigenvectors in %s.', runFile);

nNURBS = double(run.n_dofs_nurbs);
rr.iga = run.uh(1:nNURBS, eigCols);
rr.pw = run.uh(nNURBS + 1:end, eigCols);
rr.k_pw = run.k_pw;
rr.nurbs_original = run.nurbs_original;
rr.nurbs = run.nurbs_refine;
rr.meta = run.meta;
end

function sample = reference_samples(ref, quad, dxOut)
% Build reference volume fields, jumps, and outer PW fields once.
[sample.vi, sample.gxi, sample.gyi] = iga_eval_physical( ...
    ref.nurbs, ref.iga, quad.x_in, quad.y_in, ref.meta.a);
[vib, ~, ~] = iga_eval_physical( ...
    ref.nurbs, ref.iga, quad.x_bnd, quad.y_bnd, ref.meta.a);
[vob, ~, ~] = pw_eval(ref.pw, ref.k_pw, ...
    quad.x_bnd, quad.y_bnd, ref.meta.L);
sample.jump = vob - vib;
[sample.vo, sample.gxo, sample.gyo] = pw_outer_samples( ...
    ref.pw, ref.k_pw, ref.meta.L, ref.meta.a, dxOut);
end

function [grams, sample] = build_samples(cur, ref, refSample, quad, dxOut)
% Build L2 Gram matrices and common DG samples for one cutoff case.
[MrrIn, MhhIn, MrhIn] = iga_l2_grams(cur, ref);

[sample.vi, sample.gxi, sample.gyi] = iga_eval_physical( ...
    cur.nurbs, cur.iga, quad.x_in, quad.y_in, cur.meta.a);
[vib, ~, ~] = iga_eval_physical( ...
    cur.nurbs, cur.iga, quad.x_bnd, quad.y_bnd, cur.meta.a);
[vob, ~, ~] = pw_eval(cur.pw, cur.k_pw, ...
    quad.x_bnd, quad.y_bnd, cur.meta.L);
sample.jump = vob - vib;
[sample.vo, sample.gxo, sample.gyo] = pw_outer_samples( ...
    cur.pw, cur.k_pw, cur.meta.L, cur.meta.a, dxOut);

grams.Mrr = hermitian(MrrIn + refSample.vo' * refSample.vo * dxOut^2);
grams.Mhh = hermitian(MhhIn + sample.vo' * sample.vo * dxOut^2);
grams.Mrh = MrhIn + refSample.vo' * sample.vo * dxOut^2;
end

function align = l2_alignment(grams)
% Compute the L2-Procrustes rotation on each simple or repeated eigenspace.
Rr = chol(hermitian(grams.Mrr));
Rh = chol(hermitian(grams.Mhh));
C = Rr' \ (grams.Mrh / Rh);
[L, ~, R] = svd(C, 'econ');
align.Rr = Rr;
align.Rh = Rh;
align.Q = R * L';
end

function errDG = aligned_dg_error(ref, cur, quad, dxOut, sigma, align, block)
% Integrate the aligned inner, outer, and interface-jump residuals.
vi = ref.vi(:, block) / align.Rr - (cur.vi(:, block) / align.Rh) * align.Q;
gxi = ref.gxi(:, block) / align.Rr - (cur.gxi(:, block) / align.Rh) * align.Q;
gyi = ref.gyi(:, block) / align.Rr - (cur.gyi(:, block) / align.Rh) * align.Q;
vo = ref.vo(:, block) / align.Rr - (cur.vo(:, block) / align.Rh) * align.Q;
gxo = ref.gxo(:, block) / align.Rr - (cur.gxo(:, block) / align.Rh) * align.Q;
gyo = ref.gyo(:, block) / align.Rr - (cur.gyo(:, block) / align.Rh) * align.Q;
jump = ref.jump(:, block) / align.Rr - (cur.jump(:, block) / align.Rh) * align.Q;

EDG = hermitian( ...
    gram_h1(vi, gxi, gyi, quad.w_in) + ...
    gram_h1(vo, gxo, gyo, dxOut^2) + ...
    sigma * (jump' * jump * quad.w_bnd));
errDG = sqrt(max(real(diag(EDG)), 0))';
end

function [Mrr, Mhh, Mrh] = iga_l2_grams(cur, ref)
% Integrate inner-domain L2 Gram matrices on the current IGA element grid.
nurbs = cur.nurbs;
orig = cur.nurbs_original;
[gp_u, gw_u] = grule(nurbs.pu + 5);
[gp_v, gw_v] = grule(nurbs.pv + 5);
nFun = size(cur.iga, 2);
Mrr = zeros(nFun);
Mhh = zeros(nFun);
Mrh = zeros(nFun);
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

function [val, gx, gy] = pw_outer_samples(coeff, k_pw, L, a, dx)
% Evaluate a PW expansion on the common outer-domain midpoint grid.
x = -L / 2 + dx / 2:dx:L / 2 - dx / 2;
[val, gx, gy] = pw_eval_tensor(coeff, k_pw, x, x, L);
[X, Y] = ndgrid(x, x);
outer = ~(X >= -a & X <= a & Y >= -a & Y <= a);
nFun = size(coeff, 2);
val = reshape(val, [], nFun);
gx = reshape(gx, [], nFun);
gy = reshape(gy, [], nFun);
outer = outer(:);
val = val(outer, :);
gx = gx(outer, :);
gy = gy(outer, :);
end

function [val, gx, gy] = pw_eval_tensor(coeff, k_pw, x, y, L)
% Evaluate a plane-wave expansion efficiently on a tensor grid.
alpha = 1i * 2 * pi / L;
kx = unique(k_pw(:, 1));
Ex = exp(alpha * (x(:) * kx'));
nFun = size(coeff, 2);
F = zeros(numel(kx), numel(y), nFun);
Fy = zeros(numel(kx), numel(y), nFun);

for i = 1:numel(kx)
    idx = k_pw(:, 1) == kx(i);
    phaseY = exp(alpha * (k_pw(idx, 2) * y(:)'));
    for j = 1:nFun
        F(i, :, j) = sum(coeff(idx, j) .* phaseY, 1) / L;
        Fy(i, :, j) = sum((alpha * k_pw(idx, 2) .* coeff(idx, j)) .* phaseY, 1) / L;
    end
end

val = zeros(numel(x), numel(y), nFun);
gx = zeros(numel(x), numel(y), nFun);
gy = zeros(numel(x), numel(y), nFun);
for j = 1:nFun
    val(:, :, j) = Ex * F(:, :, j);
    gx(:, :, j) = Ex * ((alpha * kx) .* F(:, :, j));
    gy(:, :, j) = Ex * Fy(:, :, j);
end
end

function [val, gx, gy] = iga_eval_physical(nurbs, coeff, X, Y, a)
% Evaluate complex IGA values and physical gradients at common points.
nPts = numel(X);
nFun = size(coeff, 2);
val = zeros(nPts, nFun);
gx = zeros(nPts, nFun);
gy = zeros(nPts, nFun);
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

function [val, gx, gy] = pw_eval(coeff, k_pw, X, Y, L)
% Evaluate complex PW values and gradients at arbitrary common points.
F = [X(:)'; Y(:)'];
E = exp((1i * 2 * pi / L) * (k_pw * F));
alpha = 1i * 2 * pi / L;
nPts = numel(X);
nFun = size(coeff, 2);
val = zeros(nPts, nFun);
gx = zeros(nPts, nFun);
gy = zeros(nPts, nFun);
for j = 1:nFun
    val(:, j) = reshape(sum(coeff(:, j) .* E, 1) / L, [], 1);
    gx(:, j) = reshape(sum((coeff(:, j) .* k_pw(:, 1)) .* E, 1) * alpha / L, [], 1);
    gy(:, j) = reshape(sum((coeff(:, j) .* k_pw(:, 2)) .* E, 1) * alpha / L, [], 1);
end
end

function G = gram_h1(v, gx, gy, weight)
% Form the complex H1 Gram matrix of the residual samples.
G = (v' * v + gx' * gx + gy' * gy) * weight;
end

function quad = midpoint_quadrature(L, a, dx)
% Build the midpoint rules used for volume and interface DG terms.
assert(abs(round(L / dx) - L / dx) < 1e-12, 'dx must divide L.');
assert(abs(round(2 * a / dx) - 2 * a / dx) < 1e-12, 'dx must divide 2a.');
x = -a + dx / 2:dx:a - dx / 2;
[X, Y] = meshgrid(x, x);
quad.x_in = X(:);
quad.y_in = Y(:);
quad.w_in = dx^2;
quad.x_bnd = [x, x, -a * ones(size(x)), a * ones(size(x))]';
quad.y_bnd = [-a * ones(size(x)), a * ones(size(x)), x, x]';
quad.w_bnd = dx;
end

function A = hermitian(A)
% Suppress roundoff-level loss of Hermitian symmetry.
A = (A + A') / 2;
end
