function [A, M] = generate_A_M_NURBS_2D(g, k, L, nq, a)
% Assemble the Example 1 NURBS volume matrices on the affine square.

[gp, gw] = grule(nq);
u = basis_data(g.UBreaks, g.Ubar, g.pu, gp, gw, a);
v = basis_data(g.VBreaks, g.Vbar, g.pv, gp, gw, a);

nd = (g.pu + 1) * (g.pv + 1);
ne = g.NoEs;
nz = ne * nd^2;
ii = zeros(nz, 1);
jj = zeros(nz, 1);
av = zeros(nz, 1);
mv = zeros(nz, 1);
pos = 1;

for jv = 1:g.vNoEs
    for iu = 1:g.uNoEs
        e = iu + (jv - 1) * g.uNoEs;
        [B, Bx, By, w, x, y] = element_data(u{iu}, v{jv}, a);
        V = potential(k, L, x, y);

        Bw = B .* w.';
        Mloc = Bw * B';
        Kloc = 0.5 * ((Bx .* w.') * Bx' + (By .* w.') * By');
        Aloc = Kloc + (B .* (w .* V).') * B';

        id = pos:(pos + nd^2 - 1);
        dofs = g.Element(e, :).';
        ii(id) = repmat(dofs, nd, 1);
        jj(id) = repelem(dofs, nd, 1);
        av(id) = Aloc(:);
        mv(id) = Mloc(:);
        pos = pos + nd^2;
    end
end

A = sparse(ii, jj, av, g.n_dofs, g.n_dofs);
M = sparse(ii, jj, mv, g.n_dofs, g.n_dofs);
end

function data = basis_data(breaks, knot, p, gp, gw, a)
% Precompute one-dimensional basis values and physical derivatives.

data = cell(numel(breaks) - 1, 1);
for e = 1:numel(data)
    ab = breaks(e:e+1);
    q = ((ab(2) - ab(1)) * gp + ab(1) + ab(2)) / 2;
    N = zeros(p + 1, numel(gp));
    D = N;
    for j = 1:numel(gp)
        ders = bspbasisDers(knot, p, q(j), 1);
        N(:, j) = ders(1, :).';
        D(:, j) = ders(2, :).' / (2 * a);
    end
    data{e} = struct('N', N, 'D', D, ...
        'w', (ab(2) - ab(1)) * gw(:) / 2, ...
        'x', -a + 2 * a * q(:));
end
end

function [B, Bx, By, w, x, y] = element_data(u, v, a)
% Build tensor-product basis operators for one element.

nu = numel(u.x);
nv = numel(v.x);
nd = size(u.N, 1) * size(v.N, 1);
B = zeros(nd, nu * nv);
Bx = B;
By = B;

for j = 1:nv
    id = (j - 1) * nu + (1:nu);
    B(:, id) = kron(v.N(:, j), u.N);
    Bx(:, id) = kron(v.N(:, j), u.D);
    By(:, id) = kron(v.D(:, j), u.N);
end

w = 4 * a^2 * kron(v.w, u.w);
x = repmat(u.x, nv, 1);
y = kron(v.x, ones(nu, 1));
end

function V = potential(k, L, x, y)
% Evaluate the Example 1 radial potential in one batch.

mu = 5;
G = (2 * pi / L) * k;
gn = hypot(G(:, 1), G(:, 2));
id = gn > 0;
c = erfc(gn(id) / (2 * mu)) ./ gn(id);
r = hypot(x, y);
V = -erfc(mu * r) ./ r ...
    - (2 * pi / L^2) * (exp(1i * ([x, y] * G(id, :).')) * c) ...
    + 2 * mu / sqrt(pi);
end
