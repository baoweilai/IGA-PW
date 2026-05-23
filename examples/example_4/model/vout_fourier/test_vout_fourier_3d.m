function out = test_vout_fourier_3d(userCfg)
%Compare 3-D outer-potential quadrature methods.
arguments
    userCfg struct
end

activate_example_workflow('vout_fourier', {'config', 'core', 'operators', 'solver'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
rootDir = fullfile(exampleDir, 'data');
cfg = default_config(struct());

Nc = userCfg.Nc;
fftGridN = userCfg.fft_grid_n;
chebN = userCfg.cheb_n;
nRef = userCfg.n_ref;
tag = userCfg.tag;
methodsToRun = string(userCfg.methods);
fftPrecision = char(userCfg.fft_precision);
memoryPeakFactor = userCfg.memory_peak_factor;
memoryFreeFraction = userCfg.memory_free_fraction;

L = cfg.L;
a = cfg.a;
nuclearCharge = cfg.nuclear_charge;

[pVr, nPwVr] = build_pw_ball(1);
params = make_vout_params_local(L, a, pVr, nPwVr, nuclearCharge);
qvals = -2 * Nc:2 * Nc;
innerDomains = [-a, a, -a, a, -a, a];

outDir = fullfile(rootDir, 'result', 'VoutFourier');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

preflight_fft_memory_local(fftGridN, qvals, fftPrecision, ...
    memoryPeakFactor, memoryFreeFraction);

tRef = tic;
Vref = reference_outer_gauss_local(L, a, qvals, nRef, params);
referenceTime = toc(tRef);
refNorm = norm(Vref(:));
refSym = symmetry_residual_local(Vref);

methodNames = cell(numel(methodsToRun), 1);
Vmethods = cell(numel(methodsToRun), 1);
timings = cell(numel(methodsToRun), 1);
for im = 1:numel(methodsToRun)
    methodKey = lower(strtrim(methodsToRun(im)));
    switch methodKey
        case lower("Masked FFT")
            [Vmethods{im}, timings{im}] = masked_fft_outer_local( ...
                L, a, qvals, fftGridN, params, fftPrecision);
            methodNames{im} = 'Masked FFT';
        case lower("FFT-Chebyshev")
            [Vmethods{im}, timings{im}] = fft_cheb_outer_local( ...
                L, innerDomains, qvals, fftGridN, chebN, params, fftPrecision);
            methodNames{im} = 'FFT-Chebyshev';
        otherwise
            error('Unknown method: %s', methodsToRun(im));
    end
end

tableOut = build_result_table_local( ...
    methodNames, Nc, fftGridN, chebN, nRef, ...
    Vref, Vmethods, timings, referenceTime);

csvFile = fullfile(outDir, sprintf('vout_fourier_error_%s_Nc%d_m%d_cheb%d_nref%d.csv', ...
    tag, Nc, fftGridN, chebN, nRef));
mdFile = fullfile(outDir, sprintf('vout_fourier_error_%s_Nc%d_m%d_cheb%d_nref%d.md', ...
    tag, Nc, fftGridN, chebN, nRef));
writetable(tableOut, csvFile);
write_markdown_local(mdFile, tableOut, L, a, qvals, referenceTime, refSym);

out = struct('table', tableOut, 'csv', csvFile, 'markdown', mdFile);
end

function T = build_result_table_local(methods, Nc, fftGridN, chebN, nRef, ...
Vref, Vmethods, timings, referenceTime)
nMethod = numel(methods);
method = strings(nMethod, 1);
NcCol = zeros(nMethod, 1);
fftGridCol = zeros(nMethod, 1);
chebCol = zeros(nMethod, 1);
innerQuadCol = zeros(nMethod, 1);
errAbs2 = zeros(nMethod, 1);
errRel2 = zeros(nMethod, 1);
maxAbsErr = zeros(nMethod, 1);
cpuTotal = zeros(nMethod, 1);
timeFftSampling = zeros(nMethod, 1);
timeInnerCorrection = zeros(nMethod, 1);
symmetryResidual = zeros(nMethod, 1);
referenceN = nRef * ones(nMethod, 1);
referenceCpuTime = referenceTime * ones(nMethod, 1);
refNorm = norm(Vref(:));

for i = 1:nMethod
    V = Vmethods{i};
    method(i) = methods{i};
    NcCol(i) = Nc;
    fftGridCol(i) = fftGridN;
    if isfield(timings{i}, 'cheb_n')
        chebCol(i) = timings{i}.cheb_n;
    end
    if isfield(timings{i}, 'inner_quad_n')
        innerQuadCol(i) = timings{i}.inner_quad_n;
    end
    E = V - Vref;
    errAbs2(i) = norm(E(:));
    errRel2(i) = errAbs2(i) / refNorm;
    maxAbsErr(i) = max(abs(E(:)));
    cpuTotal(i) = timings{i}.total;
    timeFftSampling(i) = timings{i}.fft_sampling;
    timeInnerCorrection(i) = timings{i}.inner_correction;
    symmetryResidual(i) = symmetry_residual_local(V);
end

T = table(method, NcCol, fftGridCol, chebCol, innerQuadCol, ...
    errAbs2, errRel2, maxAbsErr, cpuTotal, timeFftSampling, ...
    timeInnerCorrection, symmetryResidual, referenceN, referenceCpuTime, ...
    'VariableNames', {'method', 'Nc', 'fft_grid_n', 'cheb_n', ...
    'inner_quad_n', 'err_abs_2', 'err_rel_2', 'max_abs_error', ...
    'CPU_time_total', 'time_fft_sampling', 'time_inner_correction', ...
    'symmetry_residual', 'reference_n', 'reference_CPU_time'});
end

function Vout = reference_outer_gauss_local(L, a, qvals, nRef, params)
%Integrate the outer potential by Gaussian quadrature.
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
    Vout = Vout + integrate_box_gauss_local(L, qvals, nRef, boxes(ib, :), params);
end
Vout = Vout / L ^ 3;
end

function Vbox = integrate_box_gauss_local(L, qvals, n, box, params)
%Integrate over a box by tensor-product Gauss quadrature.
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
    Vg(ix, :, :) = reshape(periodic_coulomb_slice_local(x(ix), Y, Z, YZ2, params), 1, n, n);
end

W3 = reshape(wx, n, 1, 1) .* reshape(wy, 1, n, 1) .* reshape(wz, 1, 1, n);
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

function [Vout, timing] = masked_fft_outer_local(L, a, qvals, fftGridN, params, fftPrecision)
%Compute the masked outer potential with FFT data.
tTotal = tic;
tFft = tic;
m = adjusted_fft_grid_local(fftGridN, qvals, fftPrecision);
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vgrid = zeros(m, m, m, fftPrecision);
[Y, Z] = ndgrid(x1d, x1d);
YZ2 = Y .^ 2 + Z .^ 2;
innerYZ = abs(Y) <= a & abs(Z) <= a;
for ix = 1:m
    slice = patched_potential_slice_local(x1d(ix), Y, Z, YZ2, params);
    if abs(x1d(ix)) <= a
        slice(innerYZ) = 0;
    end
    vgrid(ix, :, :) = cast(reshape(slice, 1, m, m), fftPrecision);
end
Vraw = fftn(vgrid) / m ^ 3;
clear vgrid;
Vout = double(extract_fft_coeffs_local(Vraw, qvals));
timing.fft_sampling = toc(tFft);
timing.inner_correction = 0;
timing.total = toc(tTotal);
timing.m = m;
timing.fft_precision = fftPrecision;
end

function [Vout, timing] = fft_cheb_outer_local(L, innerDomains, qvals, fftGridN, chebN, params, fftPrecision)
%Compute the outer potential with FFT-Chebyshev data.
tTotal = tic;
tFft = tic;
m = adjusted_fft_grid_local(fftGridN, qvals, fftPrecision);
dx = L / m;
x1d = -L / 2 + dx / 2 + (0:m - 1) * dx;
vgrid = zeros(m, m, m, fftPrecision);
[Y, Z] = ndgrid(x1d, x1d);
YZ2 = Y .^ 2 + Z .^ 2;
for ix = 1:m
    vgrid(ix, :, :) = cast(reshape(patched_potential_slice_local(x1d(ix), Y, Z, YZ2, params), 1, m, m), fftPrecision);
end
Vraw = fftn(vgrid) / m ^ 3;
clear vgrid;
Vfull = double(extract_fft_coeffs_local(Vraw, qvals));
timing.fft_sampling = toc(tFft);

tInner = tic;
Vin = build_inner_chebyshev_correction_3D( ...
    L, innerDomains, qvals, chebN, ...
    @(X, Y, Z) patched_potential_batch_local(X, Y, Z, params));
timing.inner_correction = toc(tInner);
timing.total = toc(tTotal);
timing.m = m;
timing.cheb_n = chebN;
timing.inner_quad_n = max(240, 4 * chebN);
timing.fft_precision = fftPrecision;
Vout = Vfull - Vin;
end

function V = extract_fft_coeffs_local(Vraw, qvals)
%Extract Fourier coefficients from grid values.
m = size(Vraw, 1);
Kq = numel(qvals);
V = zeros(Kq, Kq, Kq);
phase = exp(1i * pi * qvals * (1 - 1 / m));
for i = 1:Kq
    ii = mod(qvals(i), m) + 1;
    for j = 1:Kq
        jj = mod(qvals(j), m) + 1;
        for k = 1:Kq
            kk = mod(qvals(k), m) + 1;
            V(i, j, k) = Vraw(ii, jj, kk) * phase(i) * phase(j) * phase(k);
        end
    end
end
end

function m = adjusted_fft_grid_local(fftGridN, qvals, fftPrecision)
%Choose an FFT grid size for the requested cutoff.
m = max(round(fftGridN), 2 * max(abs(qvals)) + 1);
if mod(m, 2) ~= 0
    m = m + 1;
end
m = max(m, 8);
if strcmpi(fftPrecision, 'single')
    bytesReal = 4 * m ^ 3;
    bytesComplex = 8 * m ^ 3;
else
    bytesReal = 8 * m ^ 3;
    bytesComplex = 16 * m ^ 3;
end
end

function preflight_fft_memory_local(fftGridN, qvals, fftPrecision, ...
memoryPeakFactor, memoryFreeFraction)
m = max(round(fftGridN), 2 * max(abs(qvals)) + 1);
if mod(m, 2) ~= 0
    m = m + 1;
end
points = double(m) ^ 3;
if strcmpi(fftPrecision, 'single')
    bytesReal = 4 * points;
    bytesComplex = 8 * points;
else
    bytesReal = 8 * points;
    bytesComplex = 16 * points;
end
estimatedPeak = memoryPeakFactor * (bytesReal + bytesComplex);
freeBytes = get_free_memory_bytes_local();
if estimatedPeak > memoryFreeFraction * freeBytes
    error(['Requested %d^3 %s FFT is not safe: estimated peak %.3f GiB ', ...
        'exceeds %.1f%% of free physical memory %.3f GiB. Stop before allocation.'], ...
        m, fftPrecision, estimatedPeak / 1024 ^ 3, ...
        100 * memoryFreeFraction, freeBytes / 1024 ^ 3);
end
end

function freeBytes = get_free_memory_bytes_local()
%Return the available system memory.
assert(ispc, 'The memory preflight currently expects Windows.');
[status, txt] = system('powershell -NoProfile -Command "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"');
assert(status == 0, 'Failed to query free physical memory.');
freeKb = str2double(strtrim(txt));
assert(isfinite(freeKb), 'Free physical memory query returned a nonnumeric value.');
freeBytes = freeKb * 1024;
end

function params = make_vout_params_local(L, a, pVr, nPwVr, nuclearCharge)
%Build parameters for the Vout test.
params = struct();
params.L = L;
params.a = a;
params.a_c = 0.95 * a;
params.b = 0.5 * params.a_c;
params.p_Vr = pVr;
params.n_pw_Vr = nPwVr;
params.nuclear_charge = nuclearCharge;
params.g0 = Vr_3D_Example_1(pVr, L, nPwVr, params.b, 0, 0, nuclearCharge);
end

function Vh = patched_potential_batch_local(X, Y, Z, params)
%Evaluate the patched potential in batches.
r = sqrt(X .^ 2 + Y .^ 2 + Z .^ 2);
Vorig = periodic_coulomb_batch_local(X, Y, Z, r, params);
Vh = Vorig;

maskInner = r <= params.b;
if any(maskInner(:))
    Vh(maskInner) = params.g0;
end

maskMid = (r > params.b) & (r < params.a_c);
if any(maskMid(:))
    t = (r(maskMid) - params.b) / (params.a_c - params.b);
    eta = 1 - smooth_step_theta_local(t);
    Vh(maskMid) = (1 - eta) .* Vorig(maskMid) + eta * params.g0;
end
end

function Vh = patched_potential_slice_local(x, Y, Z, YZ2, params)
%Evaluate one slice of the patched potential.
r = sqrt(x ^ 2 + YZ2);
Vorig = periodic_coulomb_slice_local(x, Y, Z, YZ2, params);
Vh = Vorig;

maskInner = r <= params.b;
if any(maskInner(:))
    Vh(maskInner) = params.g0;
end

maskMid = (r > params.b) & (r < params.a_c);
if any(maskMid(:))
    t = (r(maskMid) - params.b) / (params.a_c - params.b);
    eta = 1 - smooth_step_theta_local(t);
    Vh(maskMid) = (1 - eta) .* Vorig(maskMid) + eta * params.g0;
end
end

function Vr = periodic_coulomb_slice_local(x, Y, Z, YZ2, params)
%Evaluate one periodic Coulomb slice.
r = sqrt(x ^ 2 + YZ2);
Vr = periodic_coulomb_batch_local(x + zeros(size(Y)), Y, Z, r, params);
end

function Vr = periodic_coulomb_batch_local(X, Y, Z, r, params)
%Evaluate the periodic Coulomb potential in batches.
ewaldAlpha = 5;
Omega = params.L ^ 3;
rSafe = max(r, 1e-14);
Vr = -params.nuclear_charge * erfc(ewaldAlpha * rSafe) ./ rSafe;

G = (2 * pi / params.L) * params.p_Vr(1:params.n_pw_Vr, :);
Gnorm2 = G(:, 1) .^ 2 + G(:, 2) .^ 2 + G(:, 3) .^ 2;
mask = Gnorm2 > 0;
if any(mask)
    Gm = G(mask, :);
    coeff = -params.nuclear_charge * (4 * pi / Omega) * ...
        exp(-Gnorm2(mask) / (4 * ewaldAlpha ^ 2)) ./ Gnorm2(mask);
    phases = bsxfun(@times, X(:), Gm(:, 1).') + ...
        bsxfun(@times, Y(:), Gm(:, 2).') + ...
        bsxfun(@times, Z(:), Gm(:, 3).');
    recip = exp(1i * phases) * coeff;
    Vr = Vr + reshape(recip, size(r));
end

Vr = real(Vr + params.nuclear_charge * 2 * ewaldAlpha / sqrt(pi));
end

function x = map_gauss_local(gp, xL, xR)
%Map Gauss points to a physical interval.
x = ((xR - xL) * gp(:) + xL + xR) / 2;
end

function th = smooth_step_theta_local(t)
%Evaluate the smooth step function.
th = sfun_local(t) ./ (sfun_local(t) + sfun_local(1 - t));
end

function y = sfun_local(t)
%Evaluate the transition polynomial.
y = zeros(size(t));
idx = t > 0;
y(idx) = exp(-1 ./ t(idx));
end

function r = symmetry_residual_local(V)
%Measure symmetry residuals in the computed potential.
den = norm(V(:));
if den == 0
    r = 0;
else
    r = norm(V(:) - reshape(conj(V(end:-1:1, end:-1:1, end:-1:1)), [], 1)) / den;
end
end

function write_markdown_local(mdFile, T, L, a, qvals, referenceTime, refSym)
%Write the Vout comparison report.
fid = fopen(mdFile, 'w');
fprintf(fid, '# Vout Fourier Coefficient Test\n\n');
fprintf(fid, '- Box: `[-L/2,L/2]^3`, `L = %.16g`\n', L);
fprintf(fid, '- Inner cube: `[-a,a]^3`, `a = %.16g`\n', a);
fprintf(fid, '- q range: `%d:%d`, tensor size `%d^3`\n', min(qvals), max(qvals), numel(qvals));
fprintf(fid, '- Sign convention: `exp(-i*kappa_q dot r)` for reference and methods.\n');
fprintf(fid, '- Reference: direct Gauss-Legendre on six outer boxes.\n');
fprintf(fid, '- Reference CPU time: `%.6g s`, symmetry residual: `%.6e`.\n\n', referenceTime, refSym);
fprintf(fid, '| method | err_abs_2 | err_rel_2 | max_abs_error | CPU_time_total | time_fft_sampling | time_inner_correction | symmetry_residual |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(T)
    fprintf(fid, '| %s | %.8e | %.8e | %.8e | %.6g | %.6g | %.6g | %.8e |\n', ...
        T.method(i), T.err_abs_2(i), T.err_rel_2(i), T.max_abs_error(i), ...
        T.CPU_time_total(i), T.time_fft_sampling(i), ...
        T.time_inner_correction(i), T.symmetry_residual(i));
end
fclose(fid);
end
