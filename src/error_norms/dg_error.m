function errDG = dg_error(referenceRunFile, caseRunFile, opt)
% Compute the DG energy error.
arguments
    referenceRunFile
    caseRunFile
    opt struct
end
if strcmp(referenceRunFile, caseRunFile), errDG = 0; return; end

Sr = load(referenceRunFile, 'run');
Sc = load(caseRunFile, 'run');
R = unpack_run_local(Sr.run, referenceRunFile);
C = unpack_run_local(Sc.run, caseRunFile);
C = align_phase_local(R, C);

[innerErr, hCase] = inner_iga_error_local(R, C, opt.innerGridN);
outerErr = outer_pw_error_local(R, C, opt.outerGridN);
jumpErr = face_jump_error_local(R, C, opt.faceGridN, opt.chunkSize);
sigma = opt.Csigma * (max(1, max(abs(C.k(:)))) + 1 / max(hCase, eps));
errDG = sqrt(max(real(innerErr + outerErr + sigma * jumpErr), 0));
end

function R = unpack_run_local(run, runFile)
% Extract fields from a saved run structure.
if ~isfield(run, 'uh') || isempty(run.uh)
    error('The run file must contain saved eigenvectors: run.uh. Missing in: %s', runFile);
end
R.uh = run.uh(:, 1);
R.k = double(run.k_pw);
R.nI = run.n_dofs_nurbs;
R.cI = R.uh(1:R.nI);
R.cP = R.uh(R.nI + 1:end);
R.nurbs = run.nurbs_refine;
R.L = run.meta.L;
R.a = run.meta.a;
end

function C = align_phase_local(R, C)
% Align coefficient columns with their reference modes.
[~, ia, ib] = intersect(R.k, C.k, 'rows');
alpha = sum(conj(R.cP(ia)) .* C.cP(ib));
if abs(alpha) > 0
    q = exp(-1i * angle(alpha));
    C.cI = C.cI * q;
    C.cP = C.cP * q;
end
end

function [val, hPhys] = inner_iga_error_local(R, C, n)
% Integrate the IGA-region DG error.
[X, Y, Z, w] = cube_points_local(-R.a, R.a, n);
[vR, gxR, gyR, gzR] = iga_val_grad_local(R.nurbs, R.cI, X, Y, Z, R.a);
[vC, gxC, gyC, gzC] = iga_val_grad_local(C.nurbs, C.cI, X, Y, Z, C.a);
val = sum(abs(vC - vR).^2 + abs(gxC - gxR).^2 + ...
    abs(gyC - gyR).^2 + abs(gzC - gzR).^2) * w;
hPhys = estimate_h_local(C.nurbs, C.a);
end

function val = outer_pw_error_local(R, C, nBase)
% Integrate the plane-wave-region DG error.
m = max([nBase, 2 * max(abs(R.k(:))) + 1, 2 * max(abs(C.k(:))) + 1, 8]);
[vR, gxR, gyR, gzR, xmid] = pw_grid_grad_local(R, m);
[vC, gxC, gyC, gzC] = pw_grid_grad_local(C, m);
[X, Y, Z] = ndgrid(xmid, xmid, xmid);
mask = ~(abs(X) <= R.a & abs(Y) <= R.a & abs(Z) <= R.a);
w = (R.L / m) ^ 3;
val = sum(abs(vC(mask) - vR(mask)).^2 + abs(gxC(mask) - gxR(mask)).^2 + ...
    abs(gyC(mask) - gyR(mask)).^2 + abs(gzC(mask) - gzR(mask)).^2) * w;
end

function val = face_jump_error_local(R, C, n, chunkSize)
% Integrate the interface jump error.
[X, Y, Z, w] = face_points_local(R.a, n);
[pR, ~, ~, ~] = pw_points_grad_local(R.cP, R.k, X, Y, Z, R.L, chunkSize);
[pC, ~, ~, ~] = pw_points_grad_local(C.cP, C.k, X, Y, Z, C.L, chunkSize);
iR = iga_val_grad_local(R.nurbs, R.cI, X, Y, Z, R.a);
iC = iga_val_grad_local(C.nurbs, C.cI, X, Y, Z, C.a);
val = sum(abs((pC - pR) - (iC - iR)).^2) * w;
end

function [X, Y, Z, w] = cube_points_local(xmin, xmax, n)
% Create midpoint quadrature points in a cube.
dx = (xmax - xmin) / n;
x = (xmin + dx / 2 : dx : xmax - dx / 2).';
[X, Y, Z] = ndgrid(x, x, x);
X = X(:); Y = Y(:); Z = Z(:);
w = dx ^ 3;
end

function [X, Y, Z, w] = face_points_local(a, n)
% Create midpoint quadrature points on cube faces.
ds = 2 * a / n;
t = (-a + ds / 2 : ds : a - ds / 2).';
[A, B] = ndgrid(t, t);
o = ones(numel(A), 1);
X = [a*o; -a*o; A(:); A(:); A(:); A(:)];
Y = [A(:); A(:); a*o; -a*o; B(:); B(:)];
Z = [B(:); B(:); B(:); B(:); a*o; -a*o];
w = ds ^ 2;
end

function [u, gx, gy, gz, xmid] = pw_grid_grad_local(R, m)
% Evaluate a plane-wave field and gradient on a grid.
raw = zeros(m, m, m);
rawX = raw; rawY = raw; rawZ = raw;
idx = sub2ind([m, m, m], mod(R.k(:, 1), m) + 1, mod(R.k(:, 2), m) + 1, mod(R.k(:, 3), m) + 1);
phase = exp(-1i * pi * sum(R.k, 2) * (1 - 1 / m)) / sqrt(R.L ^ 3);
base = R.cP(:) .* phase;
fac = 1i * 2 * pi / R.L;
raw(idx) = base;
rawX(idx) = base .* (fac * R.k(:, 1));
rawY(idx) = base .* (fac * R.k(:, 2));
rawZ(idx) = base .* (fac * R.k(:, 3));
u = ifftn(raw) * m ^ 3;
gx = ifftn(rawX) * m ^ 3;
gy = ifftn(rawY) * m ^ 3;
gz = ifftn(rawZ) * m ^ 3;
dx = R.L / m;
xmid = (-R.L / 2 + dx / 2 : dx : R.L / 2 - dx / 2).';
end

function [v, gx, gy, gz] = pw_points_grad_local(c, k, X, Y, Z, L, chunkSize)
% Evaluate a plane-wave field and gradient at points.
v = zeros(numel(X), 1); gx = v; gy = v; gz = v;
fac = 1i * 2 * pi / L;
scale = 1 / sqrt(L ^ 3);
for s = 1:chunkSize:numel(X)
    e = min(numel(X), s + chunkSize - 1);
    F = [X(s:e).'; Y(s:e).'; Z(s:e).'];
    E = exp(fac * (k * F));
    v(s:e) = (c.' * E).' * scale;
    gx(s:e) = ((c .* k(:, 1)).' * E).' * fac * scale;
    gy(s:e) = ((c .* k(:, 2)).' * E).' * fac * scale;
    gz(s:e) = ((c .* k(:, 3)).' * E).' * fac * scale;
end
end

function [v, gx, gy, gz] = iga_val_grad_local(nurbs, c, X, Y, Z, a)
% Evaluate the IGA field and gradient at points.
U = nurbs.Ubar; V = nurbs.Vbar; W = nurbs.Wbar;
pu = nurbs.pu; pv = nurbs.pv; pw = nurbs.pw;
m = nurbs.m; nn = nurbs.n;
v = zeros(numel(X), 1); gx = v; gy = v; gz = v;
for q = 1:numel(X)
    u = min(max((X(q) + a) / (2 * a), 0), 1);
    vv = min(max((Y(q) + a) / (2 * a), 0), 1);
    w = min(max((Z(q) + a) / (2 * a), 0), 1);
    is = findspan(U, pu, u); js = findspan(V, pv, vv); ks = findspan(W, pw, w);
    Bu = bspbasisDers(U, pu, u, 1); Bv = bspbasisDers(V, pv, vv, 1); Bw = bspbasisDers(W, pw, w, 1);
    iu = is-pu:is; jv = js-pv:js; kw = ks-pw:ks;
    sq = 0; sx = 0; sy = 0; sz = 0;
    for kz = 1:pw + 1
        for jy = 1:pv + 1
            for ix = 1:pu + 1
                id = iu(ix) + (jv(jy) - 1) * m + (kw(kz) - 1) * m * nn;
                cc = c(id);
                sq = sq + cc * Bu(1, ix) * Bv(1, jy) * Bw(1, kz);
                sx = sx + cc * Bu(2, ix) * Bv(1, jy) * Bw(1, kz);
                sy = sy + cc * Bu(1, ix) * Bv(2, jy) * Bw(1, kz);
                sz = sz + cc * Bu(1, ix) * Bv(1, jy) * Bw(2, kz);
            end
        end
    end
    v(q) = sq; gx(q) = sx / (2 * a); gy(q) = sy / (2 * a); gz(q) = sz / (2 * a);
end
end

function h = estimate_h_local(nurbs, a)
% Estimate the physical mesh size.
h = 2 * a * max([max(diff(unique(nurbs.Ubar))), ...
    max(diff(unique(nurbs.Vbar))), max(diff(unique(nurbs.Wbar)))]);
end
