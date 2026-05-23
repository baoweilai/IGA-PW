function [H, M, A, timing] = generate_A_M_PW_2D(L, Nc, inner_domains, p_Vr, n_pw_Vr, opts)
%Assemble 2-D plane-wave matrices with FFT-Chebyshev inner correction.

assert(isfield(opts, 'inner_cheb_n'), 'Missing opts.inner_cheb_n.');
assert(isfield(opts, 'pw_fft_grid_n'), 'Missing opts.pw_fft_grid_n.');

t_total = tic;

inner_domains = reshape(inner_domains, [], 4);
inner_cheb_n = opts.inner_cheb_n;
fft_grid_n = opts.pw_fft_grid_n;

t_build_p = tic;
[k_pw, n_pw] = build_pw_disk_local(Nc);
timing.t_build_p = toc(t_build_p);

N = floor(Nc);
qN = 2 * N;
qvals = -qN:qN;
alpha = 2 * pi / L;
Omega = L ^ 2;

t_build_U = tic;
Mker = build_outer_mass_kernel_2D(L, inner_domains, qvals);
timing.t_build_U = toc(t_build_U);

t_build_V = tic;
[Vker, VinInner, info] = build_Vker_fft_chebyshev_2D( ...
    L, N, inner_domains, p_Vr, n_pw_Vr, fft_grid_n, inner_cheb_n);
timing.t_build_V = toc(t_build_V);

t_fill = tic;
[H, M, A] = fill_pw_matrices_from_kernels(k_pw, alpha, Mker, Vker, qN);
timing.t_fill_matrices = toc(t_fill);

timing.inner_cheb_n = inner_cheb_n;
timing.inner_quad_n = info.inner_quad_n;
timing.pw_fft_grid_n = fft_grid_n;
timing.mFFT = info.mFFT;
timing.dx = info.dx;
timing.t_full_fft = info.t_full_fft;
timing.t_inner_cheb = info.t_inner_cheb;
timing.t_total = toc(t_total);
timing.n_pw = n_pw;
timing.mass_center = real(Mker(qN + 1, qN + 1));
timing.vker_center = real(Vker(qN + 1, qN + 1));
timing.vinner_norm = norm(VinInner(:));
end

function [H, M, A] = fill_pw_matrices_from_kernels(k_pw, alpha, Mker, Vker, qN)
%Fill dense PW matrices from Fourier kernels.
n_pw = size(k_pw, 1);
px = k_pw(:, 1).';
py = k_pw(:, 2).';

H = complex(zeros(n_pw, n_pw));
M = complex(zeros(n_pw, n_pw));
A = complex(zeros(n_pw, n_pw));

for i = 1:n_pw
    dkx = px - k_pw(i, 1);
    dky = py - k_pw(i, 2);
    idx = sub2ind(size(Mker), dkx + qN + 1, dky + qN + 1);

    Mi = reshape(Mker(idx), 1, []);
    Ai = 0.5 * alpha ^ 2 * (k_pw(i, 1) * px + k_pw(i, 2) * py) .* Mi;
    Vi = reshape(Vker(idx), 1, []);

    M(i, :) = Mi;
    A(i, :) = Ai;
    H(i, :) = Ai + Vi;
end

M = 0.5 * (M + M');
A = 0.5 * (A + A');
H = 0.5 * (H + H');
end

function Mker = build_outer_mass_kernel_2D(L, inner_domains, qvals)
%Compute outer-domain Fourier coefficients for the mass matrix.
alpha = 2 * pi / L;
Omega = L ^ 2;
qN = (numel(qvals) - 1) / 2;

[Qx, Qy] = ndgrid(qvals, qvals);
Uker = complex(zeros(size(Qx)));

for rr = 1:size(inner_domains, 1)
    rect = inner_domains(rr, :);
    Fx = interval_ft_general(Qx, rect(1), rect(2), alpha);
    Fy = interval_ft_general(Qy, rect(3), rect(4), alpha);
    Uker = Uker + Fx .* Fy;
end

Mker = -Uker / Omega;
Mker(qN + 1, qN + 1) = Mker(qN + 1, qN + 1) + 1;
Mker = hermitize_kernel_2d(Mker);
end

function [Vker, VinInner, info] = build_Vker_fft_chebyshev_2D( ...
    L, N, inner_domains, p_Vr, n_pw_Vr, fft_grid_n, inner_cheb_n)
%Compute outer potential coefficients from FFT data and Chebyshev correction.
mFFT = fft_grid_n;
assert(mFFT >= 4 * N + 1, 'pw_fft_grid_n is too small for this cutoff.');

dx = L / mFFT;
params = make_hatV_params(L, inner_domains, p_Vr, n_pw_Vr);

x1d = -L / 2 + dx / 2 + (0:mFFT-1) * dx;
[X, Y] = ndgrid(x1d, x1d);

t_full_fft = tic;
vext = arrayfun(@(x, y) hatV_Example1(x, y, params), X, Y);
Vfull_raw = ifftn(vext);

qN = 2 * N;
Kq = 2 * qN + 1;
qvals = -qN:qN;
Vfull = complex(zeros(Kq, Kq));

for i = 1:Kq
    qx = qvals(i);
    ii = mod(qx, mFFT) + 1;
    phx = exp(1i * pi * qx * (1 / mFFT + 1));

    for j = 1:Kq
        qy = qvals(j);
        jj = mod(qy, mFFT) + 1;
        phy = exp(1i * pi * qy * (1 / mFFT + 1));
        Vfull(i, j) = Vfull_raw(ii, jj) * phx * phy;
    end
end
info.t_full_fft = toc(t_full_fft);

t_inner_cheb = tic;
sample_fun = @(Xin, Yin) arrayfun(@(x, y) hatV_Example1(x, y, params), Xin, Yin);
[VinInner, chebInfo] = build_inner_chebyshev_correction_2D( ...
    L, inner_domains, qvals, inner_cheb_n, sample_fun);
info.t_inner_cheb = toc(t_inner_cheb);

Vker = hermitize_kernel_2d(Vfull - VinInner);

info.mFFT = mFFT;
info.dx = dx;
info.inner_quad_n = chebInfo.inner_quad_n;
end

function [VinCheb, info] = build_inner_chebyshev_correction_2D( ...
    L, inner_domains, qvals, inner_cheb_n, sample_fun)
%Integrate inner-patch potential coefficients with tensor Chebyshev data.
n = inner_cheb_n;
alpha = 2 * pi / L;
Omega = L ^ 2;
Kq = numel(qvals);

theta = pi * ((0:n-1) + 0.5) / n;
ell = 0:n-1;
Phi = cos(theta(:) * ell);
scale = (2 / n) * ones(n, 1);
scale(1) = 1 / n;
Tproj = bsxfun(@times, Phi.', scale);

nq = max(160, 4 * n);
[gp, gw] = grule(nq);
theta_q = acos(max(min(gp(:), 1), -1));
Tquad = cos((0:n-1).' * theta_q.');

VinCheb = complex(zeros(Kq, Kq));

for rr = 1:size(inner_domains, 1)
    rect = inner_domains(rr, :);
    ax = 0.5 * (rect(2) - rect(1));
    ay = 0.5 * (rect(4) - rect(3));
    cx = 0.5 * (rect(1) + rect(2));
    cy = 0.5 * (rect(3) + rect(4));

    x_nodes = cx + ax * cos(theta);
    y_nodes = cy + ay * cos(theta);
    [X, Y] = ndgrid(x_nodes, y_nodes);
    Vg = sample_fun(X, Y);
    C = Tproj * Vg * Tproj.';

    xq = cx + ax * gp(:);
    yq = cy + ay * gp(:);
    wx = ax * gw(:);
    wy = ay * gw(:);

    phaseX = exp(1i * alpha * (xq * qvals(:).'));
    phaseY = exp(1i * alpha * (yq * qvals(:).'));
    Qx = Tquad * bsxfun(@times, wx, phaseX);
    Qy = Tquad * bsxfun(@times, wy, phaseY);

    VinCheb = VinCheb + Qx.' * C * Qy / Omega;
end

info = struct();
info.inner_cheb_n = n;
info.inner_quad_n = nq;
end

function params = make_hatV_params(L, inner_domains, p_Vr, n_pw_Vr)
%Build smooth-potential parameters for the one-center example.
a = min(0.5 * (inner_domains(1, 2) - inner_domains(1, 1)), ...
    0.5 * (inner_domains(1, 4) - inner_domains(1, 3)));
a_c = 0.95 * a;
b = 0.50 * a_c;
g0 = Vr_2D_Example_1(p_Vr, L, n_pw_Vr, b, 0);

params = struct();
params.L = L;
params.a = a;
params.a_c = a_c;
params.b = b;
params.p_Vr = p_Vr;
params.n_pw_Vr = n_pw_Vr;
params.g0 = g0;
end

function Vh = hatV_Example1(x, y, params)
%Evaluate the smoothed one-center potential.
r = hypot(x, y);

if r <= params.b
    Vh = params.g0;
    return;
elseif r >= params.a_c
    Vh = Vr_2D_Example_1(params.p_Vr, params.L, params.n_pw_Vr, x, y);
    return;
end

eta = cutoff_eta(r, params.b, params.a_c);
Vorig = Vr_2D_Example_1(params.p_Vr, params.L, params.n_pw_Vr, x, y);
Vh = (1 - eta) * Vorig + eta * params.g0;
end

function eta = cutoff_eta(r, b, a_c)
%Evaluate the smooth cutoff function.
if r <= b
    eta = 1;
elseif r >= a_c
    eta = 0;
else
    t = (r - b) / (a_c - b);
    eta = 1 - smooth_step_theta(t);
end
end

function th = smooth_step_theta(t)
%Evaluate the smooth step function.
th = sfun(t) ./ (sfun(t) + sfun(1 - t));
end

function y = sfun(t)
%Evaluate the transition function.
y = zeros(size(t));
idx = (t > 0);
y(idx) = exp(-1 ./ t(idx));
end

function [k_list, n_basis] = build_pw_disk_local(Nc)
%Build the plane-wave disk basis.
N = floor(Nc);
k_list = zeros((2 * N + 1) ^ 2, 2);
n_basis = 0;
for k1 = -N:N
    m = floor(sqrt(N ^ 2 - k1 ^ 2));
    for k2 = -m:m
        n_basis = n_basis + 1;
        k_list(n_basis, :) = [k1, k2];
    end
end
k_list = k_list(1:n_basis, :);
end

function F = interval_ft_general(q, aL, aR, alpha)
%Integrate one-dimensional Fourier modes on an interval.
F = zeros(size(q));
idx0 = (q == 0);
F(idx0) = aR - aL;

qq = q(~idx0);
F(~idx0) = (exp(1i * alpha * qq * aR) - exp(1i * alpha * qq * aL)) ./ (1i * alpha * qq);
end

function Ksym = hermitize_kernel_2d(K)
%Enforce conjugate symmetry on a centered Fourier kernel.
Ksym = 0.5 * (K + conj(K(end:-1:1, end:-1:1)));
end
