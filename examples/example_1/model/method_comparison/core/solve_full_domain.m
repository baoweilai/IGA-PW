function [result, meta] = solve_full_domain(mode, cfg, opts)
%Compute fulldomain singlemethod.

format long;
tWall = tic;

mode = lower(string(mode));
L = get_opt_local(opts, 'L', 4);
Example = 'Example_1';
nEigenvalues = get_opt_local(opts, 'n_eigenvalues', 1);
solveMode = "none";
lambdaRef = get_opt_local(opts, 'lambda_ref', 4.969971740613);

ops = struct();
ops.tol = get_opt_local(opts, 'primme_tol', 1e-8);
ops.maxit = get_opt_local(opts, 'primme_maxit', 1e7);
ops.reportLevel = get_opt_local(opts, 'primme_reportLevel', 0);

seed = get_opt_local(opts, 'seed', 20260404);
rng(seed, 'twister');

tBuild = tic;
switch mode
    case "pw"
        [Mat, M, staticMeta] = build_pw_problem_local(cfg, opts, L);
    case "iga"
        [Mat, M, staticMeta] = build_iga_problem_local(cfg, opts, L, Example);
    otherwise
        error('Unsupported mode: %s', mode);
end
timeBuild = toc(tBuild);

n = size(Mat, 1);
if isreal(Mat) && isreal(M)
    v0 = randn(n, 1);
else
    v0 = randn(n, 1) + 1i * randn(n, 1);
end
v0 = v0 / max(norm(v0), eps);
ops.v0 = v0;

target = get_opt_local(opts, 'primme_target', lambdaRef);
method = get_opt_local(opts, 'primme_method', 'DEFAULT_MIN_MATVECS');

tSolve = tic;
[uhSolve, D, rnorms, ~] = solve_generalized_eigs_local( ...
    Mat, M, nEigenvalues, target, ops, method);
timeSolve = toc(tSolve);

lambda = sort(real(diag(D)), 'ascend');
[~, perm] = sort(real(diag(D)), 'ascend');
uhSolve = uhSolve(:, perm);

for i = 1:nEigenvalues
    ni = real(uhSolve(:, i)' * M * uhSolve(:, i));
    uhSolve(:, i) = uhSolve(:, i) ./ sqrt(max(abs(ni), eps));
end

uh = uhSolve;
if mode == "iga" && isfield(staticMeta, 'full_to_reduced') && ~isempty(staticMeta.full_to_reduced)
    uh = expand_periodic_coeff_local(uhSolve, staticMeta.full_to_reduced);
end

meta = staticMeta;
meta.mode = char(mode);
meta.solve_mode = char(solveMode);
meta.prec_type = 'None';
meta.lambda_ref = lambdaRef;
meta.time_build_s = timeBuild;
meta.time_solve_s = timeSolve;
meta.time_wall_s = toc(tWall);
meta.primme_rnorms = rnorms;
meta.rng_seed = seed;

result = struct();
result.lambda = lambda(:).';
result.err1 = abs(result.lambda(1) - lambdaRef);
result.uh = uh;
result.n_dofs_total = size(Mat, 1);
result.n_dofs_eval = size(uh, 1);
result.time_build_s = timeBuild;
result.time_solve_s = timeSolve;
result.time_wall_s = meta.time_wall_s;

if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end
    run = struct();
    run.mode = char(mode);
    run.lambda = result.lambda;
    run.err1 = result.err1;
    run.uh = uh;
    run.n_dofs_total = result.n_dofs_total;
    run.n_dofs_eval = result.n_dofs_eval;
    run.time_build_s = result.time_build_s;
    run.time_solve_s = result.time_solve_s;
    run.time_wall_s = result.time_wall_s;
    run.meta = meta;
    if mode == "pw"
        run.k_pw = staticMeta.k_pw;
    else
        run.nurbs_refine = staticMeta.nurbs_refine;
        if isfield(staticMeta, 'full_to_reduced')
            run.uh_reduced = uhSolve;
        end
    end
    save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
end
end

function [Mat, M, meta] = build_pw_problem_local(cfg, opts, L)
%Build PW problem.
Nc = cfg.Nc;
Nvr = get_opt_local(opts, 'potential_Nc_pw', 2);
sampleM = get_opt_local(opts, 'pw_potential_grid_m', 1024);

[kVr, nVr] = build_pw_disk_local(Nvr);
[Mat, M, kPw, cacheHit] = get_full_pw_matrices_cached_local(L, Nc, kVr, nVr, sampleM, opts);

Mat = 0.5 * (Mat + Mat');
M = 0.5 * (M + M');

meta = struct();
meta.Nc = Nc;
meta.potential_Nc = Nvr;
meta.L = L;
meta.k_pw = kPw;
meta.periodic_bc = true;
meta.n_dofs_total = size(Mat, 1);
meta.cache_hit = cacheHit;
end

function [Mat, M, meta] = build_iga_problem_local(cfg, opts, L, Example)
%Build IGA problem.
pdeg = cfg.pdeg;
refine = cfg.refine;
t = pdeg - 1;
pu0 = 1;
pv0 = 1;

[knotUe, knotVe] = IGADegreeElevSurface([0 0 1 1], [0 0 1 1], t);
pu = pu0 + t;
pv = pv0 + t;

nurbsOriginal = make_rect_patch_local([-L/2, L/2, -L/2, L/2]);
nurbsOriginal.pu = pu0;
nurbsOriginal.pv = pv0;
nurbsOriginal.knotU = [0 0 1 1];
nurbsOriginal.knotV = [0 0 1 1];

nurbsRefine = IGA_2D_Grid(knotUe, knotVe, pu, pv, refine);

Nvr = get_opt_local(opts, 'potential_Nc_iga', 2);
[kVr, nVr] = build_pw_disk_local(Nvr);
nGp = get_opt_local(opts, 'n_gp', 10);

cacheOK = get_opt_local(opts, 'use_nurbs_cache', true) ...
    && isfield(opts, 'cacheNurbsRoot') && ~isempty(opts.cacheNurbsRoot);
if cacheOK
    cacheFile = fullfile(opts.cacheNurbsRoot, sprintf( ...
        'NURBS_EX1_FULL_L_%g_p_%d_ref_%d_NVr_%d_ngp_%d.mat', ...
        L, pdeg, refine, nVr, nGp));
else
    cacheFile = '';
end

if cacheOK && exist(cacheFile, 'file')
    S = load(cacheFile, 'Mat', 'M');
    Mat = S.Mat;
    M = S.M;
    cacheHit = true;
else
    [Mat, M] = generate_A_M_NURBS_2D( ...
        nurbsOriginal, nurbsRefine, kVr, nVr, L, nGp, Example, opts);
    cacheHit = false;
    if cacheOK
        if ~exist(opts.cacheNurbsRoot, 'dir'), mkdir(opts.cacheNurbsRoot); end
        save(cacheFile, 'Mat', 'M', '-v7.3');
    end
end

Mat = 0.5 * (Mat + Mat');
M = 0.5 * (M + M');

usePeriodicBC = get_opt_local(opts, 'periodic_bc_iga', true);
fullToReduced = [];
if usePeriodicBC
    [Mat, M, fullToReduced] = apply_periodic_reduction_local(Mat, M, nurbsRefine);
end

meta = struct();
meta.pdeg = pdeg;
meta.refine = refine;
meta.h = L / (2 ^ refine);
meta.L = L;
meta.nurbs_refine = nurbsRefine;
meta.periodic_bc = usePeriodicBC;
meta.n_dofs_full = nurbsRefine.n_dofs_domains;
meta.n_dofs_total = size(Mat, 1);
meta.cache_hit = cacheHit;
if usePeriodicBC
    meta.full_to_reduced = fullToReduced;
end
end

function [Mat, M, kPw, cacheHit] = get_full_pw_matrices_cached_local(L, Nc, pVr, nVr, sampleM, opts)
%Return full PW matrices cached.
cacheOK = get_opt_local(opts, 'use_pw_cache', true) ...
    && isfield(opts, 'cacheRoot') && ~isempty(opts.cacheRoot);

if cacheOK
    cacheFile = fullfile(opts.cacheRoot, sprintf( ...
        'PW_FULL_EX1_L_%g_Nc_%d_NVr_%d_m_%d.mat', L, Nc, nVr, sampleM));
    if exist(cacheFile, 'file')
        S = load(cacheFile, 'Mat', 'M', 'kPw');
        Mat = S.Mat;
        M = S.M;
        kPw = S.kPw;
        cacheHit = true;
        return;
    end
else
    cacheFile = '';
end

[kPw, nBasis] = build_pw_disk_local(Nc);
Vfft = sample_full_potential_fft_local(L, pVr, sampleM);

M = speye(nBasis);
Mat = zeros(nBasis, nBasis);
alpha = 2 * pi / L;

for ii = 1:nBasis
    kii = kPw(ii, :);
    kinetic = 0.5 * (kii * kii') * alpha ^ 2;
    for jj = 1:nBasis
        dk = kPw(jj, :) - kii;
        if all(dk == 0)
            Mat(ii, jj) = kinetic;
        end
        dkIdx = mod(dk, sampleM) + 1;
        Vcoef = Vfft(dkIdx(2), dkIdx(1)) ...
            * exp(1i * pi * dk(1) * (1 / sampleM + 1)) ...
            * exp(1i * pi * dk(2) * (1 / sampleM + 1));
        Mat(ii, jj) = Mat(ii, jj) + Vcoef;
    end
end

Mat = sparse(real(Mat));
cacheHit = false;
if cacheOK
    if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
    save(cacheFile, 'Mat', 'M', 'kPw', '-v7.3');
end
end

function Vfft = sample_full_potential_fft_local(L, pVr, sampleM)
%Sample full potential FFT.
dx = L / sampleM;
x = -L / 2 + dx / 2 : dx : L / 2 - dx / 2;
[X, Y] = meshgrid(x, x);
V = example1_potential_grid_local(pVr, L, X, Y);
Vfft = ifftn(V);
end

function Vr = example1_potential_grid_local(p, L, X, Y)
%Compute potential grid.
alpha = 5;
p = p * 2 * pi / L;

r = hypot(X, Y);
rSafe = max(r, eps);

s = zeros(size(X));
for i = 1:size(p, 1)
    pNorm = norm(p(i, :));
    if pNorm > 0
        fac = erfc(pNorm / (2 * alpha)) / pNorm;
        s = s + fac * exp(1i * (p(i, 1) * X + p(i, 2) * Y));
    end
end

Vr = -erfc(alpha * rSafe) ./ rSafe - 2 * pi * s / (L * L) ...
    + 2 * alpha / sqrt(pi);
Vr = real(Vr);
end

function [uh, D, rnorms, stats] = solve_generalized_eigs_local( ...
Mat, M, nEigs, target, ops, method)
%Solve the full-domain eigenproblem without preconditioning.
[uh, D, rnorms, stats] = call_primme_noprec_local(Mat, M, nEigs, target, ops, method);
end

function [uh, D, rnorms, stats] = call_primme_noprec_local(Mat, M, nEigs, target, ops, method)
%Call PRIMME noprec.
if is_identity_matrix_local(M)
    [uh, D, rnorms, stats] = call_primme_standard_noprec_local(Mat, nEigs, target, ops, method);
    return;
end

if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, nEigs, target, ops, method);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, nEigs, target, ops);
end
end

function [uh, D, rnorms, stats] = call_primme_standard_noprec_local(Mat, nEigs, target, ops, method)
%Call PRIMME standard noprec.
if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, nEigs, target, ops, method);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, nEigs, target, ops);
end
end

function tf = is_identity_matrix_local(M)
%Compute identity matrix.
n = size(M, 1);
tf = size(M, 2) == n && nnz(M - speye(n)) == 0;
end

function nurbsOriginal = make_rect_patch_local(rect)
%Build rect patch.
x1 = rect(1); x2 = rect(2);
y1 = rect(3); y2 = rect(4);

ConPts = zeros(2, 2, 2);
ConPts(:, :, 1) = [x1 x1; x2 x2];
ConPts(:, :, 2) = [y1 y2; y1 y2];

nurbsOriginal = struct();
nurbsOriginal.ConPts = ConPts;
nurbsOriginal.weights = [1 1; 1 1];
end

function [kList, nBasis] = build_pw_disk_local(Nc)
%Build PW disk.
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

function [MatRed, MRed, fullToReduced] = apply_periodic_reduction_local(MatFull, MFull, nurbsRefine)
%Apply periodic reduction.
fullToReduced = build_periodic_trace_map_local(nurbsRefine.m, nurbsRefine.n);
nFull = numel(fullToReduced);
nRed = max(fullToReduced);
R = sparse((1:nFull).', fullToReduced(:), 1, nFull, nRed);
MatRed = R' * MatFull * R;
MRed = R' * MFull * R;
MatRed = 0.5 * (MatRed + MatRed');
MRed = 0.5 * (MRed + MRed');
end

function fullToReduced = build_periodic_trace_map_local(m, n)
%Build periodic trace map.
rep = reshape(1:(m * n), m, n);
rep(m, :) = rep(1, :);
rep(:, n) = rep(:, 1);
[~, ~, fullToReduced] = unique(rep(:), 'stable');
fullToReduced = double(fullToReduced(:));
end

function uhFull = expand_periodic_coeff_local(uhRed, fullToReduced)
%Compute periodic coeff.
uhFull = uhRed(fullToReduced, :);
end

function value = get_opt_local(opts, name, defaultValue)
%Return one option value.
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
