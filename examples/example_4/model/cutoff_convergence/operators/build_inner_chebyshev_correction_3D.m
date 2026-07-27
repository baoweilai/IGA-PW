function [VinCheb, info] = build_inner_chebyshev_correction_3D( ...
L, inner_domains, qvals, inner_cheb_n, sample_fun)
% Build inner-cube Fourier coefficients by tensor-product Chebyshev expansion.
% VinCheb(qx,qy,qz) = int_{[-a,a]^3} V(x,y,z) exp(-i alpha q.r) dxdydz / Omega
% The input sample_fun accepts array inputs X,Y,Z of identical size and
arguments
L
inner_domains
qvals
inner_cheb_n
sample_fun
end

if isvector(inner_domains) && numel(inner_domains) == 6
    inner_domains = reshape(inner_domains, 1, 6);
end

xL = inner_domains(1, 1); xR = inner_domains(1, 2);
yL = inner_domains(1, 3); yR = inner_domains(1, 4);
zL = inner_domains(1, 5); zR = inner_domains(1, 6);

ax = 0.5 * (xR - xL);
ay = 0.5 * (yR - yL);
az = 0.5 * (zR - zL);
cx = 0.5 * (xL + xR);
cy = 0.5 * (yL + yR);
cz = 0.5 * (zL + zR);

if abs(cx) > 1e-12 || abs(cy) > 1e-12 || abs(cz) > 1e-12 || ...
        abs(ax - ay) > 1e-12 || abs(ax - az) > 1e-12
    error('Chebyshev inner correction currently expects a centered cube [-a,a]^3.');
end

a = ax;
n = inner_cheb_n;
alpha = 2 * pi / L;
Omega = L ^ 3;

theta = pi * ((0:n-1) + 0.5) / n;
x1d = a * cos(theta);
[X, Y, Z] = ndgrid(x1d, x1d, x1d);
Vg = sample_fun(X, Y, Z);

ell = 0:n-1;
Phi = cos(theta(:) * ell);
scale = (2 / n) * ones(n, 1);
scale(1) = 1 / n;
Tproj = bsxfun(@times, Phi.', scale);

C = mode_product_local(Vg, Tproj, 1);
C = mode_product_local(C, Tproj, 2);
C = mode_product_local(C, Tproj, 3);

nq = max(240, 4 * n);
[gp, gw] = grule(nq);
xq = a * gp(:);
wq = a * gw(:);

theta_q = acos(max(min(gp(:), 1), -1));
Tquad = cos((0:n-1).' * theta_q.');
phase = exp(-1i * alpha * (xq * qvals(:).'));
Q = Tquad * bsxfun(@times, wq, phase);

VinCheb = mode_product_local(C, Q.', 1);
VinCheb = mode_product_local(VinCheb, Q.', 2);
VinCheb = mode_product_local(VinCheb, Q.', 3) / Omega;

info = struct();
info.inner_cheb_n = n;
info.inner_quad_n = nq;
info.x_nodes = x1d(:);
info.coeff_norm = norm(C(:));
end

function Y = mode_product_local(X, A, dim)
% Evaluate a tensor-product mode.
perm = [dim, 1:dim-1, dim+1:ndims(X)];
Xp = permute(X, perm);
sz = size(Xp);
Xmat = reshape(Xp, sz(1), []);
Ymat = A * Xmat;
Yp = reshape(Ymat, [size(A, 1), sz(2:end)]);
Y = ipermute(Yp, perm);
end
