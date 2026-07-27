function [E, sigma] = cutoff_eigerr(files, refFile, ids, dx)
% Aligned DG eigenfunction errors for the Example 2 cutoff sweep.

ref = load_ref(refFile, ids);
Q = quad(ref.meta.L, ref.meta.a, dx);
R = sample(ref, Q, dx);
E = zeros(numel(files), numel(ids));
sigma = zeros(numel(files), 1);

for k = 1:numel(files)
    cur = load_cur(files{k}, ids);
    sigma(k) = cur.meta.sigma;
    expected = cur.meta.beta * (2 * cur.meta.a / cur.meta.h + cur.meta.Nc);
    assert(abs(sigma(k) - expected) <= 1e-12 * expected, ...
        'Invalid penalty coefficient in %s.', files{k});

    [G, C] = data(cur, ref, R, Q, dx);
    for j = 1:numel(ids)
        nr = sqrt(real(G.rr(j, j)));
        nh = sqrt(real(G.hh(j, j)));
        z = G.rh(j, j) / (nr * nh);
        assert(nr > 0 && nh > 0 && abs(z) > 0, ...
            'Invalid L2 alignment in %s.', files{k});
        q = conj(z) / abs(z);

        vi = R.vi(:, j) / nr - C.vi(:, j) * q / nh;
        gxi = R.gxi(:, j) / nr - C.gxi(:, j) * q / nh;
        gyi = R.gyi(:, j) / nr - C.gyi(:, j) * q / nh;
        vo = R.vo(:, j) / nr - C.vo(:, j) * q / nh;
        gxo = R.gxo(:, j) / nr - C.gxo(:, j) * q / nh;
        gyo = R.gyo(:, j) / nr - C.gyo(:, j) * q / nh;
        jump = R.jump(:, j) / nr - C.jump(:, j) * q / nh;

        d = sum(abs(vi).^2 + abs(gxi).^2 + abs(gyi).^2) * Q.wi + ...
            sum(abs(vo).^2 + abs(gxo).^2 + abs(gyo).^2) * dx^2 + ...
            sigma(k) * sum(abs(jump).^2) * Q.wb;
        assert(isfinite(d) && d > 0, 'Invalid DG error in %s.', files{k});
        E(k, j) = sqrt(d);
    end

    fprintf('[EIGFUN] K=%d DG=%s\n', cur.meta.Nc, mat2str(E(k, :), 8));
end
end

function r = load_cur(file, ids)
% Load a current eigenfunction case.
S = load(file, 'run');
x = S.run;
n1 = double(x.n_dofs_1);
n2 = double(x.n_dofs_2);
r.c = {x.uh(1:n1, ids), x.uh(n1 + (1:n2), ids)};
r.pw = x.uh(n1 + n2 + 1:end, ids);
r.k = x.k_pw;
r.o = {x.nurbs_original(1), x.nurbs_original(2)};
r.n = {x.nurbs_refine(1), x.nurbs_refine(2)};
r.meta = x.meta;
end

function r = load_ref(file, ids)
% Load the reference eigenfunction case.
S = load(file, 'run');
x = S.run;
n1 = double(x.n_dofs_1);
n2 = double(x.n_dofs_2);
r.c = {x.uh(1:n1, ids), x.uh(n1 + (1:n2), ids)};
r.pw = x.uh(n1 + n2 + 1:end, ids);
r.k = x.k_pw;
r.o = {x.nurbs_original_1, x.nurbs_original_2};
r.n = {x.nurbs_refine_1, x.nurbs_refine_2};
r.meta = x.meta;
end

function S = sample(r, Q, dx)
% Sample a field on the volume and interface quadrature points.
[S.vi, S.gxi, S.gyi] = iga_eval(r, Q.xi, Q.yi, Q.pi, Q.c, Q.a);
vib = iga_eval(r, Q.xb, Q.yb, Q.pb, Q.c, Q.a);
vob = pw_eval(r.pw, r.k, Q.xb, Q.yb, Q.L);
S.jump = vob - vib;
[S.vo, S.gxo, S.gyo] = pw_outer(r.pw, r.k, Q.L, Q.c, Q.a, dx);
end

function [G, C] = data(cur, ref, R, Q, dx)
% Build the inner and outer Gram data for one comparison.
[rr, hh, rh] = iga_gram(cur, ref);
C = sample(cur, Q, dx);
G.rr = rr + R.vo' * R.vo * dx^2;
G.hh = hh + C.vo' * C.vo * dx^2;
G.rh = rh + R.vo' * C.vo * dx^2;
end

function [rr, hh, rh] = iga_gram(cur, ref)
% Assemble the IGA Gram matrices for current and reference fields.
m = size(cur.c{1}, 2);
rr = zeros(m);
hh = zeros(m);
rh = zeros(m);

for p = 1:2
    n = cur.n{p};
    o = cur.o{p};
    [gu, wu] = grule(n.pu + 5);
    [gv, wv] = grule(n.pv + 5);
    for e = 1:n.NoEs
        ue = n.Coordinate(e, 1:2);
        ve = n.Coordinate(e, 3:4);
        ju = diff(ue) / 2;
        jv = diff(ve) / 2;
        for b = 1:numel(gv)
            v = jv * gv(b) + sum(ve) / 2;
            for a = 1:numel(gu)
                u = ju * gu(a) + sum(ue) / 2;
                [~, ~, ~, ~, DF] = NurbsSurfaceDers( ...
                    o.ConPts, o.knotU, o.knotV, o.weights, o.pu, u, o.pv, v);
                w = abs(det(DF)) * wu(a) * ju * wv(b) * jv;
                vh = iga_point(n, cur.c{p}, u, v);
                vr = iga_point(ref.n{p}, ref.c{p}, u, v);
                rr = rr + vr' * vr * w;
                hh = hh + vh' * vh * w;
                rh = rh + vr' * vh * w;
            end
        end
    end
end
end

function [v, gx, gy] = iga_eval(r, X, Y, pid, centers, a)
% Evaluate the IGA field and gradient at physical points.
m = size(r.c{1}, 2);
v = zeros(numel(X), m);
gx = zeros(numel(X), m);
gy = zeros(numel(X), m);
for k = 1:numel(X)
    p = pid(k);
    u = (X(k) - centers(p, 1) + a) / (2 * a);
    t = (Y(k) - centers(p, 2) + a) / (2 * a);
    [v(k, :), du, dv] = iga_point(r.n{p}, r.c{p}, u, t);
    gx(k, :) = du / (2 * a);
    gy(k, :) = dv / (2 * a);
end
end

function [v, du, dv] = iga_point(n, c, u, t)
% Evaluate the IGA expansion at one parametric point.
i = findspan(n.Ubar, n.pu, u);
j = findspan(n.Vbar, n.pv, t);
rows = zeros((n.pu + 1) * (n.pv + 1), 1);
k = 1;
for jj = j - n.pv:j
    for ii = i - n.pu:i
        rows(k) = ii + (jj - 1) * n.m;
        k = k + 1;
    end
end
Nu = bspbasisDers(n.Ubar, n.pu, u, 1);
Nv = bspbasisDers(n.Vbar, n.pv, t, 1);
B = Nu(1, :)' * Nv(1, :);
Bu = Nu(2, :)' * Nv(1, :);
Bv = Nu(1, :)' * Nv(2, :);
v = sum(c(rows, :) .* B(:), 1);
du = sum(c(rows, :) .* Bu(:), 1);
dv = sum(c(rows, :) .* Bv(:), 1);
end

function [v, gx, gy] = pw_outer(c, k, L, centers, a, dx)
% Evaluate the plane-wave field and gradient on the outer grid.
x = -L / 2 + dx / 2:dx:L / 2 - dx / 2;
[v, gx, gy] = pw_tensor(c, k, x, L);
[X, Y] = ndgrid(x, x);
keep = true(size(X));
for p = 1:2
    keep = keep & ~(abs(X - centers(p, 1)) <= a & abs(Y - centers(p, 2)) <= a);
end
m = size(c, 2);
keep = keep(:);
v = reshape(v, [], m); v = v(keep, :);
gx = reshape(gx, [], m); gx = gx(keep, :);
gy = reshape(gy, [], m); gy = gy(keep, :);
end

function [v, gx, gy] = pw_tensor(c, k, x, L)
% Evaluate a plane-wave expansion on a tensor grid.
alpha = 1i * 2 * pi / L;
kx = unique(k(:, 1));
Ex = exp(alpha * (x(:) * kx'));
m = size(c, 2);
F = zeros(numel(kx), numel(x), m);
Fy = zeros(numel(kx), numel(x), m);
for a = 1:numel(kx)
    id = k(:, 1) == kx(a);
    Ey = exp(alpha * (k(id, 2) * x(:)'));
    for j = 1:m
        F(a, :, j) = sum(c(id, j) .* Ey, 1) / L;
        Fy(a, :, j) = sum((alpha * k(id, 2) .* c(id, j)) .* Ey, 1) / L;
    end
end
v = zeros(numel(x), numel(x), m);
gx = zeros(numel(x), numel(x), m);
gy = zeros(numel(x), numel(x), m);
for j = 1:m
    v(:, :, j) = Ex * F(:, :, j);
    gx(:, :, j) = Ex * ((alpha * kx) .* F(:, :, j));
    gy(:, :, j) = Ex * Fy(:, :, j);
end
end

function [v, gx, gy] = pw_eval(c, k, X, Y, L)
% Evaluate a plane-wave field and gradient at physical points.
alpha = 1i * 2 * pi / L;
P = exp(alpha * (k * [X(:)'; Y(:)']));
v = (c.' * P).' / L;
gx = ((c .* (alpha * k(:, 1))).' * P).' / L;
gy = ((c .* (alpha * k(:, 2))).' * P).' / L;
end

function Q = quad(L, a, dx)
% Build the volume and interface quadrature points.
Q.L = L;
Q.a = a;
Q.c = [-1 0; 1 0];
t = -a + dx / 2:dx:a - dx / 2;
n = numel(t);
Q.xi = zeros(2 * n^2, 1);
Q.yi = zeros(2 * n^2, 1);
Q.pi = zeros(2 * n^2, 1);
Q.xb = zeros(8 * n, 1);
Q.yb = zeros(8 * n, 1);
Q.pb = zeros(8 * n, 1);
for p = 1:2
    [X, Y] = meshgrid(Q.c(p, 1) + t, Q.c(p, 2) + t);
    ii = (p - 1) * n^2 + (1:n^2);
    Q.xi(ii) = X(:);
    Q.yi(ii) = Y(:);
    Q.pi(ii) = p;

    x = [Q.c(p, 1) + t, Q.c(p, 1) + t, ...
        (Q.c(p, 1) - a) * ones(1, n), (Q.c(p, 1) + a) * ones(1, n)]';
    y = [(Q.c(p, 2) - a) * ones(1, n), (Q.c(p, 2) + a) * ones(1, n), ...
        Q.c(p, 2) + t, Q.c(p, 2) + t]';
    ib = (p - 1) * 4 * n + (1:4 * n);
    Q.xb(ib) = x;
    Q.yb(ib) = y;
    Q.pb(ib) = p;
end
Q.wi = dx^2;
Q.wb = dx;
end
