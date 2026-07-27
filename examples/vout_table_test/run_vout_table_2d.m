function T = run_vout_table_2d()
% Compare the selected 2-D Vout cases.

% Set the retained cases and build the quadrature reference.
testDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(testDir));
addpath(fullfile(rootDir, 'src', 'nurbs'), '-begin');

L = 4;
a = 0.2;
Nc = 25;
maskedGridN = 20000;
chebGridNs = [300, 600, 300];
chebNs = [50, 50, 100];
nRef = 192;
nRuns = 3;
nCases = numel(chebNs);

[pVr, nPwVr] = build_pw_disk_local(2);
qvals = -2 * Nc:2 * Nc;
paramsRef = make_params_local(L, a, pVr, nPwVr);
Vout_ref = reference_outer_gauss_local(L, a, qvals, nRef, paramsRef);

t_total_masked_runs = zeros(nRuns, 1);
t_total_cheb_runs = zeros(nRuns, nCases);
t_cheb_runs = zeros(nRuns, nCases);
L2_error_masked_runs = zeros(nRuns, 1);
L2_error_cheb_runs = zeros(nRuns, nCases);

% Run each method three times.
for runIdx = 1:nRuns
    fprintf('Run %d/%d\n', runIdx, nRuns);
    [Vout_masked, t_total_masked_runs(runIdx)] = masked_fft_outer_local( ...
        L, a, qvals, maskedGridN, pVr, nPwVr);
    L2_error_masked_runs(runIdx) = norm(Vout_masked(:) - Vout_ref(:), 2);
    clear Vout_masked;

    for caseIdx = 1:nCases
        [Vout_cheb, t_total_cheb_runs(runIdx, caseIdx), ...
            t_cheb_runs(runIdx, caseIdx)] = fft_cheb_outer_local( ...
            L, a, qvals, chebGridNs(caseIdx), chebNs(caseIdx), pVr, nPwVr);
        L2_error_cheb_runs(runIdx, caseIdx) = ...
            norm(Vout_cheb(:) - Vout_ref(:), 2);
        clear Vout_cheb;
    end
end

% Average the timings and errors.
t_total_masked = mean(t_total_masked_runs);
t_total_cheb = mean(t_total_cheb_runs, 1).';
t_cheb = mean(t_cheb_runs, 1).';
L2_error_masked = mean(L2_error_masked_runs);
L2_error_cheb = mean(L2_error_cheb_runs, 1).';

% Build, display, and save the comparison table.
Dimension = repmat("2D", nCases + 1, 1);
Method = ["Masked FFT"; repmat("FFT-Chebyshev", nCases, 1)];
FFT_grid = [string(maskedGridN) + "^2"; string(chebGridNs(:)) + "^2"];
n = ["--"; string(chebNs(:))];
t_total_s = [t_total_masked; t_total_cheb(:)];
t_cheb_s = ["--"; compose("%.15g", t_cheb)];
t_cheb_fraction = ["--"; compose("%.2f%%", 100 * t_cheb ./ t_total_cheb)];
L2_error = [L2_error_masked; L2_error_cheb(:)];
T = table(Dimension, Method, FFT_grid, n, t_total_s, t_cheb_s, ...
    t_cheb_fraction, L2_error, ...
    'VariableNames', {'Dimension', 'Method', 'FFT grid', 'n', ...
    'T_total (s)', 'T_Cheb (s)', 'T_Cheb/T_total', 'L2 error'});

fprintf('%-10s %-20s %-14s %-7s %-16s %-14s %-17s %s\n', ...
    'Dimension', 'Method', 'FFT grid', 'n', 'T_total (s)', ...
    'T_Cheb (s)', 'T_Cheb/T_total', 'L2 error');
fprintf('%-10s %-20s %-14s %-7s %-16.6f %-14s %-17s %.10e\n', ...
    '2D', 'Masked FFT', char(FFT_grid(1)), '--', ...
    t_total_masked, '--', '--', L2_error_masked);
for caseIdx = 1:nCases
    fprintf('%-10s %-20s %-14s %-7d %-16.6f %-14.6f %-16.2f%% %.10e\n', ...
        '2D', 'FFT-Chebyshev', char(FFT_grid(caseIdx + 1)), ...
        chebNs(caseIdx), t_total_cheb(caseIdx), t_cheb(caseIdx), ...
        100 * t_cheb(caseIdx) / t_total_cheb(caseIdx), ...
        L2_error_cheb(caseIdx));
end

writetable(T, fullfile(testDir, 'results_2d.csv'));
end

function Vout = reference_outer_gauss_local(L, a, qvals, nRef, params)
% Compute the outer-domain Gauss reference coefficients.
halfL = L / 2;
rectangles = [
    -halfL, -a,     -halfL, halfL
     a,      halfL, -halfL, halfL
    -a,      a,     -halfL, -a
    -a,      a,      a,      halfL
    ];

Kq = numel(qvals);
Vout = complex(zeros(Kq, Kq));
for ir = 1:size(rectangles, 1)
    Vout = Vout + integrate_rectangle_gauss_local( ...
        L, qvals, nRef, rectangles(ir, :), params);
end
Vout = Vout / L ^ 2;
end

function Vrect = integrate_rectangle_gauss_local(L, qvals, n, rect, params)
% Integrate the potential over one rectangle with tensor Gauss quadrature.
alpha = 2 * pi / L;
[gp, gw] = grule(n);
x = map_gauss_local(gp, rect(1), rect(2));
y = map_gauss_local(gp, rect(3), rect(4));
wx = gw(:) * (rect(2) - rect(1)) / 2;
wy = gw(:) * (rect(4) - rect(3)) / 2;
[X, Y] = ndgrid(x, y);
Vg = periodic_coulomb_batch_local(X, Y, params);
W2 = wx * wy.';
Ex = exp(1i * alpha * (x(:) * qvals(:).'));
Ey = exp(1i * alpha * (y(:) * qvals(:).'));
Vrect = Ex.' * (Vg .* W2) * Ey;
end

function [Vout, tTotal] = masked_fft_outer_local( ...
L, a, qvals, fftGridN, pVr, nPwVr)
% Compute the outer potential coefficients with masked FFT.
tAll = tic;
params = make_params_local(L, a, pVr, nPwVr);
m = fftGridN;
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vgrid = zeros(m, m);

G = params.G;
Ex = exp(1i * (x1d(:) * G(:, 1).'));
Ey = exp(1i * (x1d(:) * G(:, 2).'));
innerY = abs(x1d) <= a;

rowBlockN = 128;
for firstRow = 1:rowBlockN:m
    rows = firstRow:min(firstRow + rowBlockN - 1, m);
    [vgrid(rows, :), ~] = periodic_potential_separable_local( ...
        x1d(rows), x1d, Ex(rows, :), Ey, params);
    innerRows = rows(abs(x1d(rows)) <= a);
    vgrid(innerRows, innerY) = 0;
end
clear Ex Ey;

Vraw = ifftn(vgrid);
clear vgrid;
Vout = extract_ifft_coeffs_local(Vraw, qvals);
tTotal = toc(tAll);
end

function [Vout, tTotal, tCheb] = fft_cheb_outer_local( ...
L, a, qvals, fftGridN, chebN, pVr, nPwVr)
% Compute the outer potential coefficients with FFT-Chebyshev.
tAll = tic;
params = make_params_local(L, a, pVr, nPwVr);
m = fftGridN;
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vext = zeros(m, m);

G = params.G;
Ex = exp(1i * (x1d(:) * G(:, 1).'));
Ey = exp(1i * (x1d(:) * G(:, 2).'));
rowBlockN = 128;
for firstRow = 1:rowBlockN:m
    rows = firstRow:min(firstRow + rowBlockN - 1, m);
    [Vbase, r] = periodic_potential_separable_local( ...
        x1d(rows), x1d, Ex(rows, :), Ey, params);
    vext(rows, :) = patched_potential_local(r, Vbase, params);
end
clear Ex Ey;

Vraw = ifftn(vext);
clear vext;
Vfull = extract_ifft_coeffs_local(Vraw, qvals);
clear Vraw;

[Vin, tCheb] = inner_chebyshev_correction_local( ...
    L, [-a, a, -a, a], qvals, chebN, params);
Vout = hermitize_kernel_2d_local(Vfull - Vin);
tTotal = toc(tAll);
end

function [VinCheb, tCheb] = inner_chebyshev_correction_local( ...
L, innerDomains, qvals, n, params)
% Compute the inner-domain Chebyshev correction and timing.
tChebStart = tic;
alpha = 2 * pi / L;
Omega = L ^ 2;
Kq = numel(qvals);

theta = pi * ((0:n - 1) + 0.5) / n;
ell = 0:n - 1;
Phi = cos(theta(:) * ell);
scale = (2 / n) * ones(n, 1);
scale(1) = 1 / n;
Tproj = bsxfun(@times, Phi.', scale);

nq = max(160, 4 * n);
[gp, gw] = grule(nq);
theta_q = acos(max(min(gp(:), 1), -1));
Tquad = cos((0:n - 1).' * theta_q.');

VinCheb = complex(zeros(Kq, Kq));
for rr = 1:size(innerDomains, 1)
    rect = innerDomains(rr, :);
    ax = 0.5 * (rect(2) - rect(1));
    ay = 0.5 * (rect(4) - rect(3));
    cx = 0.5 * (rect(1) + rect(2));
    cy = 0.5 * (rect(3) + rect(4));

    xNodes = cx + ax * cos(theta);
    yNodes = cy + ay * cos(theta);
    ExNodes = exp(1i * (xNodes(:) * params.G(:, 1).'));
    EyNodes = exp(1i * (yNodes(:) * params.G(:, 2).'));
    [Vbase, r] = periodic_potential_separable_local( ...
        xNodes, yNodes, ExNodes, EyNodes, params);
    Vg = patched_potential_local(r, Vbase, params);
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
tCheb = toc(tChebStart);
end

function V = extract_ifft_coeffs_local(Vraw, qvals)
% Extract the requested Fourier coefficients from a two-dimensional inverse FFT grid.
m = size(Vraw, 1);
idx = mod(qvals, m) + 1;
phase = exp(1i * pi * qvals * (1 / m + 1));
V = Vraw(idx, idx) .* (phase(:) * phase(:).');
end

function params = make_params_local(L, a, pVr, nPwVr)
% Build the potential parameters.
params = struct();
params.L = L;
params.a = a;
params.a_c = 0.95 * a;
params.b = 0.5 * params.a_c;
params.p_Vr = pVr;
params.n_pw_Vr = nPwVr;
params.ewald_alpha = 5;
params.G = (2 * pi / L) * pVr(1:nPwVr, :);
Gnorm = sqrt(sum(params.G .^ 2, 2));
nonzero = Gnorm > 0;
params.G = params.G(nonzero, :);
Gnorm = Gnorm(nonzero);
params.recip_coeff = -(2 * pi / L ^ 2) * ...
    erfc(Gnorm / (2 * params.ewald_alpha)) ./ Gnorm;
params.constant = 2 * params.ewald_alpha / sqrt(pi);
params.g0 = real(periodic_coulomb_batch_local(params.b, 0, params));
end

function [Vbase, r] = periodic_potential_separable_local( ...
x, y, Ex, Ey, params)
% Evaluate the separable periodic potential on a two-dimensional tensor grid.
r = hypot(x(:), y(:).');
recip = real(bsxfun(@times, Ex, params.recip_coeff.') * Ey.');
Vbase = -erfc(params.ewald_alpha * r) ./ r + ...
    recip + params.constant;
end

function Vh = patched_potential_local(r, Vbase, params)
% Evaluate the smoothly patched potential.
Vh = Vbase;
innerMask = r <= params.b;
transitionMask = r > params.b & r < params.a_c;
Vh(innerMask) = params.g0;

t = (r(transitionMask) - params.b) / (params.a_c - params.b);
eta = 1 - smooth_step_theta_local(t);
Vh(transitionMask) = (1 - eta) .* Vbase(transitionMask) + ...
    eta .* params.g0;
end

function Vr = periodic_coulomb_batch_local(X, Y, params)
% Evaluate the periodic Coulomb potential at batches of points.
r = hypot(X, Y);
phases = bsxfun(@times, X(:), params.G(:, 1).') + ...
    bsxfun(@times, Y(:), params.G(:, 2).');
recip = exp(1i * phases) * params.recip_coeff;
Vr = -erfc(params.ewald_alpha * r) ./ r + ...
    reshape(recip, size(r)) + params.constant;
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

function [kList, nBasis] = build_pw_disk_local(Nc)
% Enumerate the plane-wave modes inside a two-dimensional cutoff disk.
N = floor(Nc);
kList = zeros((2 * N + 1) ^ 2, 2);
nBasis = 0;
for k1 = -N:N
    m = floor(sqrt(N ^ 2 - k1 ^ 2));
    for k2 = -m:m
        nBasis = nBasis + 1;
        kList(nBasis, :) = [k1, k2];
    end
end
kList = kList(1:nBasis, :);
end

function Ksym = hermitize_kernel_2d_local(K)
% Enforce Hermitian symmetry in a two-dimensional Fourier kernel.
Ksym = 0.5 * (K + conj(K(end:-1:1, end:-1:1)));
end
