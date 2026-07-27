function [P, S] = assemble_DG_square_interface_fast(g, k, pw, L, a, n)
% Assemble the four Example 1 interface edges in one fast pass.

alpha = 2 * pi / L;
kx = k(:, 1);
ky = k(:, 2);
Kx = interval_ft(kx - kx.', a, alpha);
Ky = interval_ft(ky - ky.', a, alpha);

P = sparse(n, n);
S = sparse(n, n);
Pww = complex(zeros(size(k, 1)));
Sww = complex(zeros(size(k, 1)));

dv = bspbasisDers(g.Vbar, g.pv, 0, 1);
[P, S, Pww, Sww] = add_edge(P, S, Pww, Sww, ...
    g.UBreaks, g.Ubar, g.pu, dv(1, 1:2), dv(2, 1:2), ...
    g.bottom_dofs_2_layers, -a, [0; -1], k, pw, L, a, alpha, Kx, 1);

dv = bspbasisDers(g.Vbar, g.pv, 1, 1);
[P, S, Pww, Sww] = add_edge(P, S, Pww, Sww, ...
    g.UBreaks, g.Ubar, g.pu, dv(1, end-1:end), dv(2, end-1:end), ...
    g.top_dofs_2_layers, a, [0; 1], k, pw, L, a, alpha, Kx, 1);

du = bspbasisDers(g.Ubar, g.pu, 0, 1);
[P, S, Pww, Sww] = add_edge(P, S, Pww, Sww, ...
    g.VBreaks, g.Vbar, g.pv, du(1, 1:2), du(2, 1:2), ...
    g.left_dofs_2_layers, -a, [-1; 0], k, pw, L, a, alpha, Ky, 2);

du = bspbasisDers(g.Ubar, g.pu, 1, 1);
[P, S, Pww, Sww] = add_edge(P, S, Pww, Sww, ...
    g.VBreaks, g.Vbar, g.pv, du(1, end-1:end), du(2, end-1:end), ...
    g.right_dofs_2_layers, a, [1; 0], k, pw, L, a, alpha, Ky, 2);

P(pw, pw) = P(pw, pw) + sparse(Pww);
S(pw, pw) = S(pw, pw) + sparse(Sww);
end

function [P, S, Pww, Sww] = add_edge(P, S, Pww, Sww, ...
breaks, knot, p, f, df, dofs, c, normal, k, pw, L, a, alpha, K, dir)
% Accumulate one straight interface edge.

[E, F, Pm, Sm] = edge_data( ...
    breaks, knot, p, f, df, dofs, c, normal, k, L, a, alpha, K, dir);

P(dofs, dofs) = P(dofs, dofs) + E.pp;
P(dofs, pw) = P(dofs, pw) + E.pm;
P(pw, dofs) = P(pw, dofs) + E.mp;
S(dofs, dofs) = S(dofs, dofs) + F.pp;
S(dofs, pw) = S(dofs, pw) + F.pm;
S(pw, dofs) = S(pw, dofs) + F.mp;
Pww = Pww + Pm;
Sww = Sww + Sm;
end

function [E, F, Pww, Sww] = edge_data( ...
breaks, knot, p, f, df, dofs, c, normal, k, L, a, alpha, K, dir)
% Build the IGA and plane-wave blocks on one edge.

% Prepare quadrature and sparse trace triplets.
[gp, gw] = grule(10 * p + 5);
ne = numel(breaks) - 1;
nq = ne * numel(gp);
nl = 2 * (p + 1);
row = zeros(nq * nl, 1);
col = row;
bv = row;
av = row;
w = zeros(nq, 1);
t = w;
m = numel(dofs) / 2;

% Select the tangential basis ordering and normal direction.
if dir == 1
    map = [(1:m).', (m+1:2*m).'];
    pick = @(id) reshape(map(id, :), 1, []);
    shape = @(z, q) kron(q(:), z);
    normal_plus = normal(2);
else
    map = [(1:2:2*m).', (2:2:2*m).'];
    pick = @(id) reshape(map(id, :).', 1, []);
    shape = @(z, q) kron(z, q(:));
    normal_plus = normal(1);
end

pos = 1;
iq = 1;

% Sample the IGA traces and normal derivatives.
for e = 1:ne
    ab = breaks(e:e+1);
    q = ((ab(2) - ab(1)) * gp + ab(1) + ab(2)) / 2;
    wt = (ab(2) - ab(1)) * a * gw(:);
    span = findspan(knot, p, ab(1));
    active = pick(span-p:span);

    for j = 1:numel(gp)
        z = bspbasisDers(knot, p, q(j), 1);
        z = z(1, :).';
        id = pos:(pos + nl - 1);
        row(id) = iq;
        col(id) = active;
        bv(id) = shape(z, f);
        av(id) = normal_plus * shape(z, df) / (4 * a);
        w(iq) = wt(j);
        t(iq) = -a + 2 * a * q(j);
        pos = pos + nl;
        iq = iq + 1;
    end
end

% Assemble the weighted IGA trace operators.
Bp = sparse(row, col, bv, nq, numel(dofs));
Ap = sparse(row, col, av, nq, numel(dofs));
WBp = w .* Bp;
WAp = w .* Ap;

if dir == 1
    kt = k(:, 1);
    kn = k(:, 2);
else
    kt = k(:, 2);
    kn = k(:, 1);
end

% Evaluate the plane-wave traces and normal derivatives.
Bm = exp(-1i * alpha * (t * kt.' + c * kn.')) / L;
Wc = w .* conj(Bm);
mult = (1i * alpha / 2) * (k * normal).';

% Contract the quadrature data into the four DG blocks.
E.pp = Bp.' * WBp;
E.pm = -Bp.' * Wc;
E.mp = -Bm.' * WBp;
F.pp = Bp.' * WAp;
F.pm = Bp.' * (Wc .* mult);
F.mp = -Bm.' * WAp;

phase = exp(-1i * alpha * c * kn.');
Pww = ((phase.' * conj(phase)) / L^2) .* K;
Sww = -Pww .* mult;
end

function F = interval_ft(q, a, alpha)
% Evaluate the Fourier integral on [-a,a].

F = zeros(size(q));
z = q == 0;
F(z) = 2 * a;
q = q(~z);
F(~z) = (exp(1i * alpha * q * a) - exp(-1i * alpha * q * a)) ...
    ./ (1i * alpha * q);
end
