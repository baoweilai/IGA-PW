function [rhoGrid, VHGrid, aux] = build_scf_potentials_3D( ...
uh, nI, nurbs_refine, pwData, gridCache, opts)
% Build SCF densities and potentials on the Cartesian midpoint grid.
arguments
uh
nI
nurbs_refine
pwData
gridCache
opts struct
end

rho_override = opts.rhoGrid;
skip_poisson = opts.skip_poisson;
iga_grid_eval_matrix = opts.iga_grid_eval_matrix;

if ~isempty(rho_override)
    rhoGrid = rho_override;
    uGrid = [];
else
    m = gridCache.mFFT;
    uPW = eval_pw_on_grid(uh(nI+1:end), pwData, m);

    uGrid = uPW;
    if gridCache.n_inner > 0
        if ~isempty(iga_grid_eval_matrix)
            uIGA = iga_grid_eval_matrix * uh(1:nI);
        else
            uIGA = eval_iga_on_points( ...
                nurbs_refine, uh(1:nI), ...
                gridCache.x_inner, gridCache.y_inner, gridCache.z_inner, gridCache.a);
        end
        uGrid(gridCache.inner_idx, gridCache.inner_idx, gridCache.inner_idx) = ...
            reshape(uIGA, gridCache.n_inner, gridCache.n_inner, gridCache.n_inner);
    end

    rhoGrid = 2 * abs(uGrid) .^ 2;
end

if skip_poisson
    VHGrid = zeros(size(rhoGrid));
else
    VHGrid = solve_poisson_fft_zero_mode(rhoGrid, gridCache);
end

vol = gridCache.dx ^ 3;
aux.int_rho_vh = real(sum(rhoGrid(:) .* VHGrid(:))) * vol;
aux.total_charge = real(sum(rhoGrid(:))) * vol;
aux.uGrid = uGrid;
end

function u = eval_pw_on_grid(cP, pwData, m)
% Evaluate the plane-wave expansion on the Cartesian grid.
L = pwData.L;
Omega = L ^ 3;
k_pw = pwData.k_pw;

ii = mod(k_pw(:, 1), m) + 1;
jj = mod(k_pw(:, 2), m) + 1;
kk = mod(k_pw(:, 3), m) + 1;

phase = exp(-1i * pi * (sum(k_pw, 2)) * (1 - 1 / m)) / sqrt(Omega);

raw = zeros(m, m, m);
raw(sub2ind([m, m, m], ii, jj, kk)) = cP(:) .* phase;

u = ifftn(raw) * (m ^ 3);
end

function val = eval_iga_on_points(nurbs_refine, coeff, X, Y, Z, a)
% Evaluate the IGA expansion at Cartesian points.
U = nurbs_refine.Ubar; V = nurbs_refine.Vbar; W = nurbs_refine.Wbar;
pu = nurbs_refine.pu; pv = nurbs_refine.pv; pw = nurbs_refine.pw;
m = nurbs_refine.m; n = nurbs_refine.n;

nq = numel(X);
val = zeros(nq, 1);

for q = 1:nq
    u = min(max((X(q) + a) / (2 * a), 0), 1);
    v = min(max((Y(q) + a) / (2 * a), 0), 1);
    w = min(max((Z(q) + a) / (2 * a), 0), 1);

    ispan = findspan(U, pu, u);
    jspan = findspan(V, pv, v);
    kspan = findspan(W, pw, w);

    Nu = bspbasisDers(U, pu, u, 1); Nu = Nu(1, :)';
    Nv = bspbasisDers(V, pv, v, 1); Nv = Nv(1, :)';
    Nw = bspbasisDers(W, pw, w, 1); Nw = Nw(1, :)';

    ii = ispan-pu : ispan;
    jj = jspan-pv : jspan;
    kk = kspan-pw : kspan;

    s = 0;
    for kz = 1:(pw+1)
        for jy = 1:(pv+1)
            for ix = 1:(pu+1)
                gid = ii(ix) + (jj(jy)-1) * m + (kk(kz)-1) * m * n;
                s = s + coeff(gid) * Nu(ix) * Nv(jy) * Nw(kz);
            end
        end
    end
    val(q) = s;
end
end

function VH = solve_poisson_fft_zero_mode(rho, gridCache)
% Solve the periodic Poisson equation with the zero Fourier mode removed.

m = size(rho, 1);
L = gridCache.L;
q = ifftshift(-floor(m / 2):ceil(m / 2) - 1);
[Qx, Qy, Qz] = ndgrid(q, q, q);

Gx = (2 * pi / L) * Qx;
Gy = (2 * pi / L) * Qy;
Gz = (2 * pi / L) * Qz;
G2 = Gx .^ 2 + Gy .^ 2 + Gz .^ 2;

rhoHat = fftn(rho);
VHat = zeros(size(rhoHat));
mask = (G2 > 0);
VHat(mask) = 4 * pi * rhoHat(mask) ./ G2(mask);
VHat(~mask) = 0;

VH = ifftn(VHat, 'symmetric');
end
