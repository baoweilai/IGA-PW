function T = run_vout_table_3d()
% Compare the selected 3-D Vout cases.

% Set the retained cases and build the quadrature reference.
testDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(testDir));
addpath(fullfile(rootDir, 'src', 'nurbs'), '-begin');

L = 4;
a = 0.2;
nuclearCharge = 2;
Nc = 20;
maskedGridN = 1000;
chebGridN = 300;
chebN = 100;
nRef = 192;
nRuns = 3;

[pVr, nPwVr] = build_pw_ball_local(1);
qvals = -2 * Nc:2 * Nc;
innerDomains = [-a, a, -a, a, -a, a];
paramsRef = make_params_local(L, a, pVr, nPwVr, nuclearCharge);
Vout_ref = reference_outer_gauss_local(L, a, qvals, nRef, paramsRef);

t_total_masked_runs = zeros(nRuns, 1);
t_total_cheb_runs = zeros(nRuns, 1);
t_cheb_runs = zeros(nRuns, 1);
L2_error_masked_runs = zeros(nRuns, 1);
L2_error_cheb_runs = zeros(nRuns, 1);

% Run each method three times.
for runIdx = 1:nRuns
    fprintf('Run %d/%d\n', runIdx, nRuns);
    [Vout_masked, t_total_masked_runs(runIdx)] = masked_fft_outer_local( ...
        L, a, qvals, maskedGridN, pVr, nPwVr, nuclearCharge);
    L2_error_masked_runs(runIdx) = norm(Vout_masked(:) - Vout_ref(:), 2);
    clear Vout_masked;

    [Vout_cheb, t_total_cheb_runs(runIdx), t_cheb_runs(runIdx)] = ...
        fft_cheb_outer_local(L, innerDomains, qvals, chebGridN, chebN, ...
        pVr, nPwVr, nuclearCharge);
    L2_error_cheb_runs(runIdx) = norm(Vout_cheb(:) - Vout_ref(:), 2);
    clear Vout_cheb;
end

% Average the timings and errors.
t_total_masked = mean(t_total_masked_runs);
t_total_cheb = mean(t_total_cheb_runs);
t_cheb = mean(t_cheb_runs);
L2_error_masked = mean(L2_error_masked_runs);
L2_error_cheb = mean(L2_error_cheb_runs);

% Build, display, and save the comparison table.
Dimension = ["3D"; "3D"];
Method = ["Masked FFT"; "FFT-Chebyshev"];
FFT_grid = ["1000^3"; "300^3"];
n = ["--"; string(chebN)];
t_total_s = [t_total_masked; t_total_cheb];
t_cheb_s = ["--"; compose("%.15g", t_cheb)];
t_cheb_fraction = ["--"; compose("%.2f%%", 100 * t_cheb / t_total_cheb)];
L2_error = [L2_error_masked; L2_error_cheb];
T = table(Dimension, Method, FFT_grid, n, t_total_s, t_cheb_s, ...
    t_cheb_fraction, L2_error, ...
    'VariableNames', {'Dimension', 'Method', 'FFT grid', 'n', ...
    'T_total (s)', 'T_Cheb (s)', 'T_Cheb/T_total', 'L2 error'});

fprintf('%-10s %-20s %-14s %-7s %-16s %-14s %-17s %s\n', ...
    'Dimension', 'Method', 'FFT grid', 'n', 'T_total (s)', ...
    'T_Cheb (s)', 'T_Cheb/T_total', 'L2 error');
fprintf('%-10s %-20s %-14s %-7s %-16.6f %-14s %-17s %.10e\n', ...
    '3D', 'Masked FFT', '1000^3', '--', ...
    t_total_masked, '--', '--', L2_error_masked);
fprintf('%-10s %-20s %-14s %-7d %-16.6f %-14.6f %-16.2f%% %.10e\n', ...
    '3D', 'FFT-Chebyshev', '300^3', chebN, t_total_cheb, t_cheb, ...
    100 * t_cheb / t_total_cheb, L2_error_cheb);

writetable(T, fullfile(testDir, 'results_3d.csv'));
end

function Vout = reference_outer_gauss_local(L, a, qvals, nRef, params)
% Compute the outer-domain Gauss reference coefficients.
halfL = L / 2;
boxes = [
    -halfL, -a,    -halfL, halfL, -halfL, halfL
     a,      halfL, -halfL, halfL, -halfL, halfL
    -a,      a,     -halfL, -a,    -halfL, halfL
    -a,      a,      a,      halfL, -halfL, halfL
    -a,      a,     -a,      a,     -halfL, -a
    -a,      a,     -a,      a,      a,      halfL
    ];

Kq = numel(qvals);
Vout = zeros(Kq, Kq, Kq);
for ib = 1:size(boxes, 1)
    Vout = Vout + integrate_box_gauss_local( ...
        L, qvals, nRef, boxes(ib, :), params);
end
Vout = Vout / L ^ 3;
end

function Vbox = integrate_box_gauss_local(L, qvals, n, box, params)
% Integrate the potential over one box with tensor Gauss quadrature.
alpha = 2 * pi / L;
[gp, gw] = grule(n);
x = map_gauss_local(gp, box(1), box(2));
y = map_gauss_local(gp, box(3), box(4));
z = map_gauss_local(gp, box(5), box(6));
wx = gw(:) * (box(2) - box(1)) / 2;
wy = gw(:) * (box(4) - box(3)) / 2;
wz = gw(:) * (box(6) - box(5)) / 2;

Kq = numel(qvals);
Vg = zeros(n, n, n);
[Y, Z] = ndgrid(y, z);
YZ2 = Y .^ 2 + Z .^ 2;
for ix = 1:n
    Vg(ix, :, :) = reshape(periodic_coulomb_slice_local( ...
        x(ix), Y, Z, YZ2, params), 1, n, n);
end

W3 = reshape(wx, n, 1, 1) .* reshape(wy, 1, n, 1) .* ...
    reshape(wz, 1, 1, n);
T = Vg .* W3;
Ex = exp(-1i * alpha * (x(:) * qvals(:).'));
Ey = exp(-1i * alpha * (y(:) * qvals(:).'));
Ez = exp(-1i * alpha * (z(:) * qvals(:).'));
A = Ex.' * reshape(T, n, n * n);
A = reshape(A, Kq, n, n);
B = zeros(Kq, Kq, n);
for iz = 1:n
    B(:, :, iz) = A(:, :, iz) * Ey;
end
Vbox = reshape(reshape(B, Kq * Kq, n) * Ez, Kq, Kq, Kq);
end

function [Vout, tTotal] = masked_fft_outer_local( ...
L, a, qvals, fftGridN, pVr, nPwVr, nuclearCharge)
% Compute the outer potential coefficients with masked FFT.
tAll = tic;
params = make_params_local(L, a, pVr, nPwVr, nuclearCharge);
m = fftGridN;
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vgrid = zeros(m, m, m);
[Y, Z] = ndgrid(x1d, x1d);
innerYZ = abs(Y) <= a & abs(Z) <= a;
Ey = exp(1i * (x1d(:) * params.G(:, 2).'));
Ez = exp(1i * (x1d(:) * params.G(:, 3).'));
for ix = 1:m
    [slice, ~] = periodic_potential_separable_slice_local( ...
        x1d(ix), x1d, x1d, Ey, Ez, params);
    if abs(x1d(ix)) <= a
        slice(innerYZ) = 0;
    end
    vgrid(ix, :, :) = reshape(slice, 1, m, m);
end
Vraw = fftn(vgrid) / m ^ 3;
clear vgrid;
Vout = extract_fft_coeffs_local(Vraw, qvals);
tTotal = toc(tAll);
end

function [Vout, tTotal, tCheb] = fft_cheb_outer_local( ...
L, innerDomains, qvals, fftGridN, chebN, pVr, nPwVr, nuclearCharge)
% Compute the outer potential coefficients with FFT-Chebyshev.
tAll = tic;
params = make_params_local(L, innerDomains(2), ...
    pVr, nPwVr, nuclearCharge);
m = fftGridN;
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vgrid = zeros(m, m, m);
Ey = exp(1i * (x1d(:) * params.G(:, 2).'));
Ez = exp(1i * (x1d(:) * params.G(:, 3).'));
for ix = 1:m
    [Vbase, r] = periodic_potential_separable_slice_local( ...
        x1d(ix), x1d, x1d, Ey, Ez, params);
    slice = patched_potential_local(r, Vbase, params);
    vgrid(ix, :, :) = reshape(slice, 1, m, m);
end
Vraw = fftn(vgrid) / m ^ 3;
clear vgrid;
Vfull = extract_fft_coeffs_local(Vraw, qvals);
clear Vraw;

[Vin, tCheb] = inner_chebyshev_correction_local( ...
    L, innerDomains, qvals, chebN, params);
Vout = Vfull - Vin;
tTotal = toc(tAll);
end

function [VinCheb, tCheb] = inner_chebyshev_correction_local( ...
L, innerDomains, qvals, n, params)
% Compute the inner-domain Chebyshev correction and timing.
tChebStart = tic;
a = 0.5 * (innerDomains(2) - innerDomains(1));
alpha = 2 * pi / L;
Omega = L ^ 3;

theta = pi * ((0:n - 1) + 0.5) / n;
x1d = a * cos(theta);
Vg = zeros(n, n, n);
Ey = exp(1i * (x1d(:) * params.G(:, 2).'));
Ez = exp(1i * (x1d(:) * params.G(:, 3).'));
for ix = 1:n
    [Vbase, r] = periodic_potential_separable_slice_local( ...
        x1d(ix), x1d, x1d, Ey, Ez, params);
    slice = patched_potential_local(r, Vbase, params);
    Vg(ix, :, :) = reshape(slice, 1, n, n);
end
ell = 0:n - 1;
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
Tquad = cos((0:n - 1).' * theta_q.');
phase = exp(-1i * alpha * (xq * qvals(:).'));
Q = Tquad * bsxfun(@times, wq, phase);

VinCheb = mode_product_local(C, Q.', 1);
VinCheb = mode_product_local(VinCheb, Q.', 2);
VinCheb = mode_product_local(VinCheb, Q.', 3) / Omega;
tCheb = toc(tChebStart);
end

function V = extract_fft_coeffs_local(Vraw, qvals)
% Extract the requested Fourier coefficients from a three-dimensional FFT grid.
m = size(Vraw, 1);
idx = mod(qvals, m) + 1;
phase = exp(1i * pi * qvals * (1 - 1 / m));
P = reshape(phase, [], 1, 1) .* reshape(phase, 1, [], 1) .* ...
    reshape(phase, 1, 1, []);
V = Vraw(idx, idx, idx) .* P;
end

function params = make_params_local(L, a, pVr, nPwVr, nuclearCharge)
% Build the potential parameters.
params = struct();
params.L = L;
params.a = a;
params.a_c = 0.95 * a;
params.b = 0.5 * params.a_c;
params.p_Vr = pVr;
params.n_pw_Vr = nPwVr;
params.nuclear_charge = nuclearCharge;
params.ewald_alpha = 5;
params.G = (2 * pi / L) * pVr(1:nPwVr, :);
Gnorm2 = sum(params.G .^ 2, 2);
nonzero = Gnorm2 > 0;
params.G = params.G(nonzero, :);
Gnorm2 = Gnorm2(nonzero);
params.recip_coeff = -nuclearCharge * (4 * pi / L ^ 3) * ...
    exp(-Gnorm2 / (4 * params.ewald_alpha ^ 2)) ./ Gnorm2;
params.constant = nuclearCharge * 2 * params.ewald_alpha / sqrt(pi);
params.g0 = periodic_coulomb_batch_local( ...
    params.b, 0, 0, params.b, params);
end

function [Vbase, r] = periodic_potential_separable_slice_local( ...
x, y, z, Ey, Ez, params)
% Evaluate the separable periodic potential on one three-dimensional grid slice.
r = sqrt(x ^ 2 + y(:) .^ 2 + (z(:).') .^ 2);
weightedCoeff = params.recip_coeff .* exp(1i * x * params.G(:, 1));
recip = real(bsxfun(@times, Ey, weightedCoeff.') * Ez.');
rSafe = max(r, 1e-14);
Vbase = -params.nuclear_charge * ...
    erfc(params.ewald_alpha * rSafe) ./ rSafe + recip + params.constant;
end

function Vh = patched_potential_local(r, Vbase, params)
% Evaluate the smoothly patched potential.
Vh = Vbase;
maskInner = r <= params.b;
Vh(maskInner) = params.g0;
maskMid = (r > params.b) & (r < params.a_c);
t = (r(maskMid) - params.b) / (params.a_c - params.b);
eta = 1 - smooth_step_theta_local(t);
Vh(maskMid) = (1 - eta) .* Vbase(maskMid) + eta * params.g0;
end

function Vr = periodic_coulomb_slice_local(x, Y, Z, YZ2, params)
% Evaluate the periodic Coulomb potential on one grid slice.
r = sqrt(x ^ 2 + YZ2);
Vr = periodic_coulomb_batch_local(x + zeros(size(Y)), Y, Z, r, params);
end

function Vr = periodic_coulomb_batch_local(X, Y, Z, r, params)
% Evaluate the periodic Coulomb potential at batches of points.
rSafe = max(r, 1e-14);
Vr = -params.nuclear_charge * ...
    erfc(params.ewald_alpha * rSafe) ./ rSafe;
phases = bsxfun(@times, X(:), params.G(:, 1).') + ...
    bsxfun(@times, Y(:), params.G(:, 2).') + ...
    bsxfun(@times, Z(:), params.G(:, 3).');
recip = exp(1i * phases) * params.recip_coeff;
Vr = real(Vr + reshape(recip, size(r)) + params.constant);
end

function Y = mode_product_local(X, A, dim)
% Multiply a tensor along one mode.
perm = [dim, 1:dim - 1, dim + 1:ndims(X)];
Xp = permute(X, perm);
sz = size(Xp);
Xmat = reshape(Xp, sz(1), []);
Ymat = A * Xmat;
Yp = reshape(Ymat, [size(A, 1), sz(2:end)]);
Y = ipermute(Yp, perm);
end

function th = smooth_step_theta_local(t)
% Evaluate the normalized smooth transition.
th = sfun_local(t) ./ (sfun_local(t) + sfun_local(1 - t));
end

function y = sfun_local(t)
% Evaluate the smooth cutoff seed function.
y = zeros(size(t));
idx = t > 0;
y(idx) = exp(-1 ./ t(idx));
end

function x = map_gauss_local(gp, xL, xR)
% Map Gauss nodes to a physical interval.
x = ((xR - xL) * gp(:) + xL + xR) / 2;
end

function [kList, nBasis] = build_pw_ball_local(Nc)
% Enumerate the plane-wave modes inside a three-dimensional cutoff ball.
N = floor(Nc);
kList = zeros((2 * N + 1) ^ 3, 3);
nBasis = 0;
for k1 = -N:N
    for k2 = -N:N
        radius2 = N ^ 2 - k1 ^ 2 - k2 ^ 2;
        if radius2 < 0
            continue;
        end
        m = floor(sqrt(radius2));
        for k3 = -m:m
            nBasis = nBasis + 1;
            kList(nBasis, :) = [k1, k2, k3];
        end
    end
end
kList = kList(1:nBasis, :);
end
