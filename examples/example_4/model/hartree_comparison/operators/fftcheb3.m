function [pwData, timing] = fftcheb3(cfg, K, innerDomains, pVr, nPwVr)
% Build the three-dimensional plane-wave data with FFT-Chebyshev contraction.
arguments
    cfg struct
    K (1,1) double
    innerDomains double
    pVr double
    nPwVr (1,1) double
end

% Read the cell, grid, and smoothing parameters.
tTotal = tic;

L = cfg.Lcell;
mFFT = cfg.fftGridN;
nCheb = cfg.chebDegree;
mu = cfg.mu;
charge = cfg.charge;

assert(K >= 0 && K == floor(K), 'K must be a nonnegative integer.');
assert(L > 0, 'cfg.Lcell must be positive.');
assert(mFFT == 300, 'cfg.fftGridN must be 300.');
assert(nCheb == 80, 'cfg.chebDegree must be 80.');
assert(mu == 5, 'cfg.mu must be 5.');
assert(charge == 2, 'cfg.charge must be 2.');
assert(nPwVr >= 1 && nPwVr == floor(nPwVr), 'nPwVr must be a positive integer.');
assert(size(pVr, 2) == 3 && size(pVr, 1) >= nPwVr, ...
    'pVr must contain at least nPwVr three-dimensional modes.');

innerDomains = reshape(innerDomains, 1, []);
assert(numel(innerDomains) == 6, 'innerDomains must contain six bounds.');
xL = innerDomains(1); xR = innerDomains(2);
yL = innerDomains(3); yR = innerDomains(4);
zL = innerDomains(5); zR = innerDomains(6);
ax = 0.5 * (xR - xL);
ay = 0.5 * (yR - yL);
az = 0.5 * (zR - zL);
cx = 0.5 * (xL + xR);
cy = 0.5 * (yL + yR);
cz = 0.5 * (zL + zR);
assert(ax > 0 && abs(ax - ay) < 1e-14 && abs(ax - az) < 1e-14, ...
    'innerDomains must define a cube.');
assert(max(abs([cx, cy, cz])) < 1e-14, ...
    'innerDomains must define a centered cube.');

% Build the plane-wave basis and the outer-domain mass kernel.
N = K;
cubeN = 2 * N + 1;
qMax = 2 * N;
qvals = -qMax:qMax;
qCubeN = numel(qvals);
assert(mFFT >= 4 * N + 1, 'cfg.fftGridN is too small for K.');

tPw = tic;
[kPw, nPw] = build_pw_ball_local(N);
timing.t_build_p = toc(tPw);

subs = kPw + N + 1;
linBall = sub2ind([cubeN, cubeN, cubeN], subs(:, 1), subs(:, 2), subs(:, 3));
maskCube = false(cubeN, cubeN, cubeN);
maskCube(linBall) = true;
[kxCube, kyCube, kzCube] = ndgrid(-N:N, -N:N, -N:N);

waveScale = 2 * pi / L;
Omega = L ^ 3;

tU = tic;
[Qx, Qy, Qz] = ndgrid(qvals, qvals, qvals);
Fx = interval_ft_local(Qx, xL, xR, waveScale);
Fy = interval_ft_local(Qy, yL, yR, waveScale);
Fz = interval_ft_local(Qz, zL, zR, waveScale);
Uker = hermitize_kernel_local((Fx .* Fy .* Fz) / Omega);
timing.t_build_U = toc(tU);

params = smooth_params_local(L, ax, pVr(1:nPwVr, :), mu, charge);

% Sample the smoothed potential and transform it by FFT.
tV = tic;
dx = L / mFFT;
xmid = -L / 2 + dx / 2 + (0:mFFT-1) * dx;

tFull = tic;
vext = smooth_grid_local(xmid, params);
VfullRaw = fftn(vext) / (mFFT ^ 3);
freqIdx = mod(qvals, mFFT) + 1;
midPhase = exp(1i * pi * qvals * (1 - 1 / mFFT));
Vfull = VfullRaw(freqIdx, freqIdx, freqIdx) ...
    .* reshape(midPhase, qCubeN, 1, 1) ...
    .* reshape(midPhase, 1, qCubeN, 1) ...
    .* reshape(midPhase, 1, 1, qCubeN);
timing.t_full_fft = toc(tFull);
clear vext VfullRaw

% Subtract the inner-domain Chebyshev contribution.
tInner = tic;
[VinInner, innerInfo] = inner_cheb_local( ...
    L, ax, qvals, nCheb, @(X, Y, Z) smooth_batch_local(X, Y, Z, params));
timing.t_inner_cheb = toc(tInner);

Vker = hermitize_kernel_local(Vfull - VinInner);
timing.t_build_V = toc(tV);

% Build the convolution kernels and diagonal approximations.
convN = cubeN + qCubeN - 1;
Ufft = fftn(ifftshift(embed_center_local(Uker, convN)));
Vfft = fftn(ifftshift(embed_center_local(Vker, convN)));

volumeInner = (xR - xL) * (yR - yL) * (zR - zL);
mass0 = 1 - volumeInner / Omega;
V0 = real(Vker(qMax + 1, qMax + 1, qMax + 1));
k2 = sum(kPw .^ 2, 2);

% Package the plane-wave operators and timing data.
pwData = struct();
pwData.L = L;
pwData.Nc = K;
pwData.N = N;
pwData.K = cubeN;
pwData.Kq = qCubeN;
pwData.k_pw = kPw;
pwData.n_pw = nPw;
pwData.inner_cheb_n = nCheb;
pwData.pw_fft_grid_n = mFFT;
pwData.nuclear_charge = charge;
pwData.mu = mu;
pwData.maskCube = maskCube;
pwData.linBall = linBall;
pwData.kxCube = kxCube;
pwData.kyCube = kyCube;
pwData.kzCube = kzCube;
pwData.alpha = waveScale;
pwData.Omega = Omega;
pwData.Uker = Uker;
pwData.Vker = Vker;
pwData.VinInner = VinInner;
pwData.VinCheb = VinInner;
pwData.convN = convN;
pwData.Ufft = Ufft;
pwData.Vfft = Vfft;
pwData.massDiag = mass0 * ones(nPw, 1);
pwData.stiffDiag = 0.5 * waveScale ^ 2 * k2 * mass0 + V0;
pwData.facePenaltyDiagApprox = 24 * ax ^ 2 / Omega;

timing.mFFT = mFFT;
timing.dx = dx;
timing.inner_cheb_n = nCheb;
timing.inner_quad_n = innerInfo.inner_quad_n;
timing.t_total = toc(tTotal);
end

function [kPw, nPw] = build_pw_ball_local(N)
% Enumerate the plane-wave modes inside a three-dimensional cutoff ball.
kPw = zeros((2 * N + 1) ^ 3, 3);
nPw = 0;
for k1 = -N:N
    for k2 = -N:N
        rem2 = N ^ 2 - k1 ^ 2 - k2 ^ 2;
        if rem2 >= 0
            m = floor(sqrt(rem2));
            kz = (-m:m).';
            count = numel(kz);
            ids = nPw + (1:count);
            kPw(ids, 1) = k1;
            kPw(ids, 2) = k2;
            kPw(ids, 3) = kz;
            nPw = nPw + count;
        end
    end
end
kPw = kPw(1:nPw, :);
end

function F = interval_ft_local(q, left, right, waveScale)
% Evaluate the Fourier transform of an interval indicator.
F = zeros(size(q));
zeroMask = (q == 0);
F(zeroMask) = right - left;
qq = q(~zeroMask);
F(~zeroMask) = (exp(1i * waveScale * qq * right) ...
    - exp(1i * waveScale * qq * left)) ./ (1i * waveScale * qq);
end

function params = smooth_params_local(L, a, pVr, mu, charge)
% Build the potential-smoothing parameters.
params = struct();
params.L = L;
params.Omega = L ^ 3;
params.a = a;
params.ac = 0.95 * a;
params.b = 0.50 * params.ac;
params.mu = mu;
params.charge = charge;
params.G = (2 * pi / L) * pVr;
g2 = sum(params.G .^ 2, 2);
params.recipCoeff = zeros(size(g2));
nonzeroMask = (g2 > 0);
params.recipCoeff(nonzeroMask) = -charge * (4 * pi / params.Omega) ...
    * exp(-g2(nonzeroMask) / (4 * mu ^ 2)) ./ g2(nonzeroMask);
params.g0 = raw_batch_local(params.b, 0, 0, params);
end

function Vhat = smooth_grid_local(xmid, params)
% Evaluate the smoothed potential on the full grid.
mFFT = numel(xmid);
Vhat = zeros(mFFT, mFFT, mFFT);
[Y, Z] = ndgrid(xmid, xmid);
YZ2 = Y .^ 2 + Z .^ 2;
nModes = size(params.G, 1);
phaseYZ = cell(nModes, 1);
for ig = 1:nModes
    if params.recipCoeff(ig) ~= 0
        phaseYZ{ig} = exp(1i * (params.G(ig, 2) * Y + params.G(ig, 3) * Z));
    end
end

for ix = 1:mFFT
    x = xmid(ix);
    r = sqrt(x ^ 2 + YZ2);
    safeR = max(r, 1e-14);
    Vorig = -params.charge * erfc(params.mu * safeR) ./ safeR;
    for ig = 1:nModes
        if params.recipCoeff(ig) ~= 0
            Vorig = Vorig + params.recipCoeff(ig) ...
                * exp(1i * params.G(ig, 1) * x) * phaseYZ{ig};
        end
    end
    Vorig = real(Vorig + params.charge * 2 * params.mu / sqrt(pi));
    Vhat(ix, :, :) = reshape(blend_slice_local(Vorig, r, params), 1, mFFT, mFFT);
end
end

function Vhat = blend_slice_local(Vorig, r, params)
% Blend the smooth and raw potential values on one grid slice.
Vhat = Vorig;
innerMask = (r <= params.b);
Vhat(innerMask) = params.g0;
middleMask = (r > params.b) & (r < params.ac);
if any(middleMask(:))
    t = (r(middleMask) - params.b) / (params.ac - params.b);
    eta = 1 - smooth_theta_local(t);
    Vhat(middleMask) = (1 - eta) .* Vorig(middleMask) + eta * params.g0;
end
end

function Vhat = smooth_batch_local(X, Y, Z, params)
% Evaluate the smoothed periodic potential in batches.
r = sqrt(X .^ 2 + Y .^ 2 + Z .^ 2);
Vorig = raw_batch_local(X, Y, Z, params);
Vhat = Vorig;
innerMask = (r <= params.b);
Vhat(innerMask) = params.g0;
middleMask = (r > params.b) & (r < params.ac);
if any(middleMask(:))
    t = (r(middleMask) - params.b) / (params.ac - params.b);
    eta = 1 - smooth_theta_local(t);
    Vhat(middleMask) = (1 - eta) .* Vorig(middleMask) + eta * params.g0;
end
end

function V = raw_batch_local(X, Y, Z, params)
% Evaluate the unmodified periodic potential in batches.
r = sqrt(X .^ 2 + Y .^ 2 + Z .^ 2);
safeR = max(r, 1e-14);
V = -params.charge * erfc(params.mu * safeR) ./ safeR;
nonzeroMask = (params.recipCoeff ~= 0);
G = params.G(nonzeroMask, :);
coeff = params.recipCoeff(nonzeroMask);
if ~isempty(coeff)
    phases = X(:) * G(:, 1).' + Y(:) * G(:, 2).' + Z(:) * G(:, 3).';
    reciprocalPart = exp(1i * phases) * coeff;
    V = V + reshape(reciprocalPart, size(r));
end
V = real(V + params.charge * 2 * params.mu / sqrt(pi));
end

function theta = smooth_theta_local(t)
% Evaluate the normalized smooth transition.
left = smooth_seed_local(t);
right = smooth_seed_local(1 - t);
theta = left ./ (left + right);
end

function y = smooth_seed_local(t)
% Evaluate the potential-smoothing seed function.
y = zeros(size(t));
positiveMask = (t > 0);
y(positiveMask) = exp(-1 ./ t(positiveMask));
end

function [VinInner, info] = inner_cheb_local(L, a, qvals, nCheb, sampleFun)
% Compute the inner-domain Chebyshev correction.
waveScale = 2 * pi / L;
Omega = L ^ 3;
angles = pi * ((0:nCheb-1) + 0.5) / nCheb;
nodes = a * cos(angles);
[X, Y, Z] = ndgrid(nodes, nodes, nodes);
Vg = sampleFun(X, Y, Z);

degrees = 0:nCheb-1;
Phi = cos(angles(:) * degrees);
scale = (2 / nCheb) * ones(nCheb, 1);
scale(1) = 1 / nCheb;
projection = bsxfun(@times, Phi.', scale);
C = mode_product_local(Vg, projection, 1);
C = mode_product_local(C, projection, 2);
C = mode_product_local(C, projection, 3);

nQuad = max(240, 4 * nCheb);
[gp, gw] = grule(nQuad);
xq = a * gp(:);
wq = a * gw(:);
quadAngles = acos(max(min(gp(:), 1), -1));
Tquad = cos((0:nCheb-1).' * quadAngles.');
phase = exp(-1i * waveScale * (xq * qvals(:).'));
Q = Tquad * bsxfun(@times, wq, phase);

VinInner = mode_product_local(C, Q.', 1);
VinInner = mode_product_local(VinInner, Q.', 2);
VinInner = mode_product_local(VinInner, Q.', 3) / Omega;

info = struct();
info.inner_quad_n = nQuad;
end

function Y = mode_product_local(X, A, dim)
% Multiply a tensor along one mode.
order = [dim, 1:dim-1, dim+1:ndims(X)];
Xp = permute(X, order);
shape = size(Xp);
Xmat = reshape(Xp, shape(1), []);
Ymat = A * Xmat;
Yp = reshape(Ymat, [size(A, 1), shape(2:end)]);
Y = ipermute(Yp, order);
end

function A = embed_center_local(X, targetN)
% Embed an array in the center of a larger grid.
A = zeros(targetN, targetN, targetN);
sourceN = size(X, 1);
first = floor((targetN - sourceN) / 2) + 1;
A(first:first+sourceN-1, first:first+sourceN-1, first:first+sourceN-1) = X;
end

function Kherm = hermitize_kernel_local(Kraw)
% Enforce Hermitian symmetry in a Fourier kernel.
Kherm = 0.5 * (Kraw + conj(Kraw(end:-1:1, end:-1:1, end:-1:1)));
end
