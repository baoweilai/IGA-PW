function [pwData, timing] = generate_A_M_PW_3D( ...
L, Nc, inner_domains, p_Vr, n_pw_Vr, fft_grid_n, opts)
% Assemble 3-D plane-wave operators with FFT-Chebyshev correction.

arguments
L
Nc
inner_domains
p_Vr
n_pw_Vr
fft_grid_n
opts struct
end

inner_cheb_n = opts.inner_cheb_n;
nuclear_charge = opts.nuclear_charge;

t_total = tic;

if isvector(inner_domains) && numel(inner_domains)==6
    inner_domains = reshape(inner_domains, 1, 6);
end

xL = inner_domains(1,1); xR = inner_domains(1,2);
yL = inner_domains(1,3); yR = inner_domains(1,4);
zL = inner_domains(1,5); zR = inner_domains(1,6);

% ---------------- PW basis ----------------
t_build_p = tic;
[k_pw, n_pw] = build_pw_ball(Nc);
timing.t_build_p = toc(t_build_p);

N  = floor(Nc);
K  = 2*N + 1;
qN = 2*N;
Kq = 2*qN + 1;

pwData = struct();
pwData.L             = L;
pwData.Nc            = Nc;
pwData.N             = N;
pwData.K             = K;
pwData.Kq            = Kq;
pwData.k_pw          = k_pw;
pwData.n_pw          = n_pw;
pwData.inner_cheb_n = inner_cheb_n;
pwData.pw_fft_grid_n = fft_grid_n;
pwData.nuclear_charge = nuclear_charge;

maskCube = false(K,K,K);
linBall  = zeros(n_pw,1);

for r = 1:n_pw
    ii = k_pw(r,1) + N + 1;
    jj = k_pw(r,2) + N + 1;
    kk = k_pw(r,3) + N + 1;
    maskCube(ii,jj,kk) = true;
    linBall(r) = sub2ind([K,K,K], ii, jj, kk);
end

pwData.maskCube = maskCube;
pwData.linBall  = linBall;

kvals = -N:N;
[kxCube, kyCube, kzCube] = ndgrid(kvals, kvals, kvals);
pwData.kxCube = kxCube;
pwData.kyCube = kyCube;
pwData.kzCube = kzCube;

alpha = 2*pi/L;
Omega = L^3;
pwData.alpha = alpha;
pwData.Omega = Omega;

% ---------------- Uker ----------------
t_build_U = tic;
qvals = -qN:qN;
[Qx, Qy, Qz] = ndgrid(qvals, qvals, qvals);

Fx = interval_ft_general(Qx, xL, xR, alpha);
Fy = interval_ft_general(Qy, yL, yR, alpha);
Fz = interval_ft_general(Qz, zL, zR, alpha);

Uker = (Fx .* Fy .* Fz) / Omega;
Uker = hermitize_kernel_3d(Uker);

timing.t_build_U = toc(t_build_U);

% ---------------- Vker = full FFTN - inner correction ----------------
t_build_V = tic;
[Vker, VinInner, info] = build_Vker_fullFFT_innerCorrection( ...
    L, N, inner_domains, p_Vr, n_pw_Vr, fft_grid_n, ...
    inner_cheb_n, nuclear_charge);
Vker = hermitize_kernel_3d(Vker);
timing.t_build_V = toc(t_build_V);

pwData.Uker     = Uker;
pwData.Vker     = Vker;
pwData.VinInner = VinInner;
pwData.VinCheb  = VinInner;

convN = K + Kq - 1;
pwData.convN = convN;

pwData.Ufft = fftn(ifftshift(embed_center(Uker, convN)));
pwData.Vfft = fftn(ifftshift(embed_center(Vker, convN)));

vol_in = (xR-xL) * (yR-yL) * (zR-zL);
massDiag0 = 1 - vol_in / Omega;
V0 = real(Vker(qN+1, qN+1, qN+1));

k2 = k_pw(:,1).^2 + k_pw(:,2).^2 + k_pw(:,3).^2;
pwData.massDiag  = massDiag0 * ones(n_pw,1);
pwData.stiffDiag = 0.5 * alpha^2 * k2 * massDiag0 + V0;

if abs(xL + xR) < 1e-14 && abs(yL + yR) < 1e-14 && abs(zL + zR) < 1e-14 ...
        && abs((xR-xL)-(yR-yL)) < 1e-14 && abs((xR-xL)-(zR-zL)) < 1e-14
    a = xR;
    pwData.facePenaltyDiagApprox = 24 * a^2 / Omega;
else
    area_box = 2 * ((yR-yL)*(zR-zL) + (xR-xL)*(zR-zL) + (xR-xL)*(yR-yL));
    pwData.facePenaltyDiagApprox = area_box / Omega;
end

timing.mFFT          = info.mFFT;
timing.dx            = info.dx;
timing.inner_cheb_n  = info.inner_cheb_n;
timing.inner_quad_n  = info.inner_quad_n;
timing.t_full_fft    = info.t_full_fft;
timing.t_inner_cheb  = info.t_inner_cheb;
timing.t_total       = toc(t_total);

end

function F = interval_ft_general(q, aL, aR, alpha)
% Integrate one Fourier mode over an arbitrary interval.
F = zeros(size(q));
idx0 = (q == 0);
F(idx0) = aR - aL;

qq = q(~idx0);
F(~idx0) = (exp(1i * alpha * qq * aR) - exp(1i * alpha * qq * aL)) ./ (1i * alpha * qq);
end

function [Vker, VinInner, info] = build_Vker_fullFFT_innerCorrection( ...
L, N, inner_domains, p_Vr, n_pw_Vr, fft_grid_n, ...
    inner_cheb_n, nuclear_charge)
% Compute outer potential coefficients from FFT data and Chebyshev correction.
% Prepare the FFT grid and smoothed-potential parameters.
mFFT = fft_grid_n;
assert(mFFT >= 4 * N + 1, 'pw_fft_grid_n is too small for this cutoff.');

dx = L / mFFT;
params = make_hatV_params_3D(L, inner_domains, p_Vr, n_pw_Vr, nuclear_charge);

% Transform the full-cell potential and extract requested frequencies.
x1d = -L/2 + dx/2 + (0:mFFT-1)*dx;
t_full_fft = tic;
vext = build_hatV_grid_vectorized_local(x1d, params);

Vfull_raw = fftn(vext) / (mFFT^3);

qN = 2*N;
Kq = 2*qN + 1;
qvals = -qN:qN;

Vfull = zeros(Kq, Kq, Kq);
for i = 1:Kq
    qx = qvals(i);
    ii = mod(qx, mFFT) + 1;
    phx = exp(1i*pi*qx*(1 - 1/mFFT));

    for j = 1:Kq
        qy = qvals(j);
        jj = mod(qy, mFFT) + 1;
        phy = exp(1i*pi*qy*(1 - 1/mFFT));

        for k = 1:Kq
            qz = qvals(k);
            kk = mod(qz, mFFT) + 1;
            phz = exp(1i*pi*qz*(1 - 1/mFFT));

            Vfull(i,j,k) = Vfull_raw(ii,jj,kk) * phx * phy * phz;
        end
    end
end
t_full_fft = toc(t_full_fft);

% Compute and subtract the inner-domain Chebyshev correction.
t_inner_cheb = tic;
[VinInner, chebInfo] = build_inner_chebyshev_correction_3D( ...
    L, inner_domains, qvals, inner_cheb_n, ...
    @(X, Y, Z) hatV_Example1_3D_batch(X, Y, Z, params));
t_inner_cheb = toc(t_inner_cheb);

Vker = Vfull - VinInner;

% Package the grid and timing data.
info = struct();
info.mFFT          = mFFT;
info.dx            = dx;
info.inner_cheb_n  = inner_cheb_n;
info.inner_quad_n  = chebInfo.inner_quad_n;
info.t_full_fft    = t_full_fft;
info.t_inner_cheb  = t_inner_cheb;
end

function params = make_hatV_params_3D(L, inner_domains, p_Vr, n_pw_Vr, nuclear_charge)
% Build parameters for the 3-D smooth potential.
xL = inner_domains(1,1); xR = inner_domains(1,2);
a = 0.5*(xR - xL);

a_c = 0.95 * a;
b   = 0.50 * a_c;

g0 = Vr_3D_Example_1(p_Vr, L, n_pw_Vr, b, 0, 0, nuclear_charge);

params = struct();
params.L       = L;
params.a       = a;
params.a_c     = a_c;
params.b       = b;
params.p_Vr    = p_Vr;
params.n_pw_Vr = n_pw_Vr;
params.nuclear_charge = nuclear_charge;
params.g0      = g0;
end

function vext = build_hatV_grid_vectorized_local(x1d, params)
% Evaluate the smooth potential on a grid.
mFFT = numel(x1d);
vext = zeros(mFFT, mFFT, mFFT);
[Y, Z] = ndgrid(x1d, x1d);
YZ2 = Y .^ 2 + Z .^ 2;

for ix = 1:mFFT
    x = x1d(ix);
    vext(ix, :, :) = reshape(hatV_slice_local(x, Y, Z, YZ2, params), 1, mFFT, mFFT);
end
end

function Vh = hatV_slice_local(x, Y, Z, YZ2, params)
% Evaluate one slice of the smooth potential.
r = sqrt(x ^ 2 + YZ2);
Vorig = Vr_3D_Example_1_batch_local(x, Y, Z, r, params);

if all(r(:) <= params.b)
    Vh = params.g0 * ones(size(r));
    return;
end
if all(r(:) >= params.a_c)
    Vh = Vorig;
    return;
end

eta = zeros(size(r));
eta(r <= params.b) = 1;
mask_mid = (r > params.b) & (r < params.a_c);
if any(mask_mid(:))
    t = (r(mask_mid) - params.b) / (params.a_c - params.b);
    eta(mask_mid) = 1 - smooth_step_theta(t);
end

Vh = Vorig;
mask_eta = eta ~= 0;
Vh(mask_eta) = (1 - eta(mask_eta)) .* Vorig(mask_eta) + eta(mask_eta) * params.g0;
end

function Vh = hatV_Example1_3D_batch(X, Y, Z, params)
% Evaluate the reciprocal-space potential term on an array.
r = sqrt(X .^ 2 + Y .^ 2 + Z .^ 2);
Vorig = Vr_3D_Example_1_batch_local(X, Y, Z, r, params);

Vh = Vorig;
mask_inner = (r <= params.b);
if any(mask_inner(:))
    Vh(mask_inner) = params.g0;
end

mask_mid = (r > params.b) & (r < params.a_c);
if any(mask_mid(:))
    t = (r(mask_mid) - params.b) / (params.a_c - params.b);
    eta = 1 - smooth_step_theta(t);
    Vh(mask_mid) = (1 - eta) .* Vorig(mask_mid) + eta * params.g0;
end
end

function Vr = Vr_3D_Example_1_batch_local(X, Y, Z, r, params)
% Evaluate the 3-D Example 1 potential on an array.
alpha = 5;
Omega = params.L ^ 3;
r_safe = max(r, 1e-14);
Vr = -params.nuclear_charge * erfc(alpha * r_safe) ./ r_safe;

G = (2 * pi / params.L) * params.p_Vr(1:params.n_pw_Vr, :);
if isempty(G)
    Vr = real(Vr + params.nuclear_charge * 2 * alpha / sqrt(pi));
    return;
end

Gnorm2 = G(:, 1) .^ 2 + G(:, 2) .^ 2 + G(:, 3) .^ 2;
mask = Gnorm2 > 0;
if any(mask)
    Gm = G(mask, :);
    coeff = -params.nuclear_charge * (4 * pi / Omega) * ...
        exp(-Gnorm2(mask) / (4 * alpha ^ 2)) ./ Gnorm2(mask);
    phases = bsxfun(@times, X(:), Gm(:, 1).') + ...
        bsxfun(@times, Y(:), Gm(:, 2).') + ...
        bsxfun(@times, Z(:), Gm(:, 3).');
    recip = exp(1i * phases) * coeff;
    Vr = Vr + reshape(recip, size(r));
end

Vr = real(Vr + params.nuclear_charge * 2 * alpha / sqrt(pi));
end

function th = smooth_step_theta(t)
% Evaluate the smooth step function.
th = sfun(t) ./ (sfun(t) + sfun(1 - t));
end

function y = sfun(t)
% Evaluate the transition polynomial.
y = zeros(size(t));
idx = (t > 0);
y(idx) = exp(-1 ./ t(idx));
end

function A = embed_center(X, convN)
% Embed a smaller array at the grid center.
A = zeros(convN, convN, convN);
K = size(X,1);
i0 = floor((convN - K)/2) + 1;
A(i0:i0+K-1, i0:i0+K-1, i0:i0+K-1) = X;
end

function Ksym = hermitize_kernel_3d(K)
% Enforce Hermitian symmetry on a centered 3-D kernel.
Ksym = 0.5 * (K + conj(K(end:-1:1, end:-1:1, end:-1:1)));
end
