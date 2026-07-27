function run_method_data(resultDir)
% Compute Example 3 method-comparison error fields.

assert(exist('resultDir', 'var') == 1, 'run_method_data requires resultDir.');
ensure_dir_local(resultDir);
exampleDir = fileparts(fileparts(resultDir));

% Set the three method cases and the reference case.
cfgPW = struct('Nc', 30);
cfgIGA = struct('pdeg', 1, 't', 0, 'nElem', 64);
cfgHybrid = struct('Nc', 10, 'pdeg', 1, 't', 0, 'nElem', 40);
cfgRef = struct('Nc', 40, 'pdeg', 2, 'refine', 7);
parameters = struct('pw', cfgPW, 'iga', cfgIGA, ...
    'iga_pw', cfgHybrid, 'reference', cfgRef);

opts = make_common_opts_local(resultDir);

fprintf('\n============================================================\n');
fprintf('[RUN ] Example 3 PW / IGA / IGA-PW absolute-error fields\n');
fprintf('[CASE] PW: Nc=%d\n', cfgPW.Nc);
fprintf('[CASE] IGA: p=%d, nElem=%d\n', cfgIGA.pdeg, cfgIGA.nElem);
fprintf('[CASE] IGA-PW: Nc=%d, p=%d, nElem=%d\n', ...
    cfgHybrid.Nc, cfgHybrid.pdeg, cfgHybrid.nElem);

% Load or compute the reference and three method solutions.
refRun = load_reference_run_local(exampleDir, cfgRef);
pwRun = load_or_run_pw_case_local(fullfile(resultDir, 'pw'), cfgPW, opts);
igaRun = load_or_run_iga_case_local(fullfile(resultDir, 'iga'), cfgIGA, opts);
igapwRun = run_hybrid_case_local(fullfile(resultDir, 'iga_pw'), cfgHybrid, opts);

% Evaluate and align all solutions on the common grid.
gridN = 401;
[Xg, Yg] = meshgrid(linspace(-2, 2, gridN), linspace(-2, 2, gridN));
dxg = Xg(1, 2) - Xg(1, 1);

UrefRaw = evaluate_run_on_grid_local(refRun, 'igapw', Xg, Yg);
[Uref, ~] = normalize_and_align_local(UrefRaw, [], dxg);

UpwRaw = evaluate_run_on_grid_local(pwRun, 'pw', Xg, Yg);
[Upw, ~] = normalize_and_align_local(UpwRaw, Uref, dxg);

UigaRaw = evaluate_run_on_grid_local(igaRun, 'iga', Xg, Yg);
[Uiga, ~] = normalize_and_align_local(UigaRaw, Uref, dxg);

UigapwRaw = evaluate_run_on_grid_local(igapwRun, 'igapw', Xg, Yg);
[Uigapw, ~] = normalize_and_align_local(UigapwRaw, Uref, dxg);

methodLabels = {'PW', 'IGA', 'IGA-PW'};
sourceIters = [get_scf_iterations_local(pwRun), get_scf_iterations_local(igaRun), ...
    get_scf_iterations_local(igapwRun)];

errorFields = {abs(Upw - Uref), abs(Uiga - Uref), abs(Uigapw - Uref)};
maxErrors = [max(errorFields{1}(:)), max(errorFields{2}(:)), max(errorFields{3}(:))];
l2Errors = [l2_error_from_grid_local(errorFields{1}, dxg), ...
    l2_error_from_grid_local(errorFields{2}, dxg), ...
    l2_error_from_grid_local(errorFields{3}, dxg)];
refNorm = sqrt(sum(abs(Uref(:)).^2) * dxg * dxg);
relL2Errors = l2Errors ./ max(refNorm, eps);
globalMax = max(maxErrors);

% Save the aligned error fields and error measures.
save(fullfile(resultDir, 'fields.mat'), 'Xg', 'Yg', 'Uref', 'methodLabels', ...
    'errorFields', 'sourceIters', 'maxErrors', 'l2Errors', ...
    'relL2Errors', 'globalMax', 'parameters', '-v7.3');

fprintf('[DONE] fields.mat saved: %s\n', fullfile(resultDir, 'fields.mat'));
fprintf('       L2 errors: PW=%.6e, IGA=%.6e, IGA-PW=%.6e\n', ...
    l2Errors(1), l2Errors(2), l2Errors(3));
fprintf('============================================================\n\n');
end

function opts = make_common_opts_local(resultDir)
% Build the common solver options.
opts = struct();
opts.Example = 'Example_3';
opts.beta = 100;
opts.n_gp = 30;
opts.inner_cheb_n = 150;
opts.pw_fft_grid_n = 550;
opts.primme_tol = 1e-9;
opts.primme_maxit = 1e8;
opts.primme_method = 'DEFAULT_MIN_TIME';
opts.primme_reportLevel = 0;
opts.block_targetShift = 0.0;
opts.eps_diag = 1e-12;
opts.iface_reg = 1e-12;
opts.scf_maxit = 70;
opts.scf_pw_grid_m = 600;
opts.scf_tol_lambda = 1e-7;
opts.scf_mixing = 0.9;
opts.scf_track_n_eigs = 1;
opts.use_pw_cache = true;
opts.use_nurbs_cache = true;
opts.cacheRoot = fullfile(resultDir, 'cache_pw');
opts.cacheNurbsRoot = fullfile(resultDir, 'cache_nurbs');
ensure_dir_local(opts.cacheRoot);
ensure_dir_local(opts.cacheNurbsRoot);
end

function run = load_reference_run_local(exampleDir, cfg)
% Load the saved reference run.
runFile = fullfile(exampleDir, 'data', 'result', 'Example_3', ...
    sprintf('refine_%02d', cfg.refine), sprintf('p_%d', cfg.pdeg), ...
    sprintf('Nc_%02d', cfg.Nc), 'run.mat');
assert(exist(runFile, 'file') == 2, 'Missing reference run: %s', runFile);
S = load(runFile, 'run');
run = S.run;
assert(isfield(run, 'uh') && isfield(run, 'k_pw') && isfield(run, 'nurbs_refine'), ...
    'Reference run lacks field data: %s', runFile);
end

function run = run_hybrid_case_local(outDir, cfg, opts)
% Compute and load the IGA-PW case.
ensure_dir_local(outDir);
fprintf('[CASE] IGA-PW solve started: Nc=%d, p=%d, nElem=%d\n', ...
    cfg.Nc, cfg.pdeg, cfg.nElem);
[~, ~, ~, run] = solve_iga_pw_dg(cfg.nElem, cfg.t, cfg.Nc, 1, opts);
run.mode = 'iga_pw';
save(fullfile(outDir, 'run.mat'), 'run', '-v7.3');
fprintf('[CASE] IGA-PW solve finished.\n');
end

function run = load_or_run_pw_case_local(outDir, cfg, opts)
% Load or compute the full-domain plane-wave case.
runFile = fullfile(outDir, 'run.mat');
if exist(runFile, 'file')
    S = load(runFile, 'run');
    run = S.run;
    if isfield(run, 'mode') && strcmpi(run.mode, 'pw') && ...
            isfield(run, 'meta') && isfield(run.meta, 'Nc') && ...
            double(run.meta.Nc) == double(cfg.Nc)
        fprintf('[CASE] PW loaded: Nc=%d\n', cfg.Nc);
        return;
    end
end
run = run_pw_case_local(outDir, cfg, opts);
end

function run = load_or_run_iga_case_local(outDir, cfg, opts)
% Load or compute the full-domain IGA case.
runFile = fullfile(outDir, 'run.mat');
if exist(runFile, 'file')
    S = load(runFile, 'run');
    run = S.run;
    if isfield(run, 'mode') && strcmpi(run.mode, 'iga') && ...
            isfield(run, 'meta') && isfield(run.meta, 'pdeg') && ...
            isfield(run.meta, 'nElem') && ...
            double(run.meta.pdeg) == double(cfg.pdeg) && ...
            double(run.meta.nElem) == double(cfg.nElem)
        fprintf('[CASE] IGA loaded: p=%d, nElem=%d\n', cfg.pdeg, cfg.nElem);
        return;
    end
end
run = run_iga_case_local(outDir, cfg, opts);
end

function run = run_pw_case_local(outDir, cfg, opts)
% Compute the full-domain plane-wave case.
ensure_dir_local(outDir);
fprintf('[CASE] PW solve started: Nc=%d\n', cfg.Nc);
L = 4;
N_Vr = 1;
[k_Vr, n_pw_Vr] = build_pw_disk_local(N_Vr);
[Mat0, M, k_pw] = build_full_pw_problem_cached_local(L, cfg.Nc, k_Vr, n_pw_Vr, opts);
buildN = @(u) assemble_nonlinear_pw_from_grid( ...
    k_pw, sample_rho_pw_full_grid_local(u, k_pw, L, opts.scf_pw_grid_m), L);
fixPhase = @(u) fix_global_phase_single_pw_local(u, k_pw);
[lambda1, u1, info] = solve_scf_full_local(Mat0, M, buildN, fixPhase, opts);

run = struct();
run.mode = 'pw';
run.lambda = lambda1;
run.uh = u1;
run.k_pw = k_pw;
run.n_dofs_total = numel(u1);
run.n_dofs_eval = numel(u1);
run.meta = struct('mode', 'pw', 'Nc', cfg.Nc, 'L', L, ...
    'n_pw_basis', size(k_pw, 1), 'scf_iterations', info.n_iters, ...
    'scf_abslambda', info.abslambda, 'scf_converged', info.converged);
save(fullfile(outDir, 'run.mat'), 'run', '-v7.3');
fprintf('[CASE] PW solve finished.\n');
end

function run = run_iga_case_local(outDir, cfg, opts)
% Compute the full-domain IGA case.
ensure_dir_local(outDir);
fprintf('[CASE] IGA solve started: p=%d, nElem=%d\n', cfg.pdeg, cfg.nElem);
L = 4;
pu0 = 1;
pv0 = 1;
pu = pu0 + cfg.t;
pv = pv0 + cfg.t;

nurbsBase = make_rect_patch_local([-L/2, L/2, -L/2, L/2]);
nurbsBase.pu = pu0;
nurbsBase.pv = pv0;
nurbsBase.knotU = [0 0 1 1];
nurbsBase.knotV = [0 0 1 1];
nurbs_refine = IGA_2D_Grid_nElem([], [], pu, pv, cfg.nElem);

N_Vr = 1;
[k_Vr, n_pw_Vr] = build_pw_disk_local(N_Vr);
[MatFull, MFull] = build_full_iga_problem_cached_local( ...
    nurbsBase, nurbs_refine, cfg.nElem, cfg.t, k_Vr, n_pw_Vr, L, opts);
[Mat0, M, fullToReduced, R] = apply_periodic_reduction_local(MatFull, MFull, nurbs_refine);

buildN = @(u) R' * assemble_nonlinear_nurbs( ...
    nurbsBase, nurbs_refine, ...
    sample_rho_nurbs(expand_periodic_coeff_local(u, fullToReduced), ...
    nurbsBase, nurbs_refine, opts.n_gp), opts.n_gp) * R;
fixPhase = @(u) u;
[lambda1, uRed, info] = solve_scf_full_local(Mat0, M, buildN, fixPhase, opts);
uFull = expand_periodic_coeff_local(uRed, fullToReduced);

run = struct();
run.mode = 'iga';
run.lambda = lambda1;
run.uh = uFull;
run.uh_reduced = uRed;
run.nurbs_base = nurbsBase;
run.nurbs_refine = nurbs_refine;
run.n_dofs_total = numel(uRed);
run.n_dofs_eval = numel(uFull);
run.meta = struct('mode', 'iga', 'pdeg', cfg.pdeg, 'nElem', cfg.nElem, ...
    'h', L / cfg.nElem, 'L', L, 'periodic_bc', true, ...
    'n_dofs_full', numel(uFull), 'n_dofs_reduced', numel(uRed), ...
    'scf_iterations', info.n_iters, 'scf_abslambda', info.abslambda, ...
    'scf_converged', info.converged);
save(fullfile(outDir, 'run.mat'), 'run', '-v7.3');
fprintf('[CASE] IGA solve finished.\n');
end

function [lambda1, u1, info] = solve_scf_full_local(Mat0, M, buildN, fixPhase, opts)
% Solve the self-consistent full-domain eigenproblem.
ops = struct();
ops.tol = opts.primme_tol;
ops.maxit = opts.primme_maxit;
ops.reportLevel = opts.primme_reportLevel;

nTrack = min(opts.scf_track_n_eigs, size(Mat0, 1) - 1);
nTrack = max(nTrack, 1);
[U0, D0] = solve_lowest_local(Mat0, M, nTrack, ops, opts);
uPrev = normalize_in_M_local(fixPhase(U0(:, 1)), M);
lambdaPrev = real(D0(1, 1));

lambdaHist = zeros(opts.scf_maxit, 1);
info = struct('n_iters', 0, 'abslambda', inf, 'converged', false);
scfMixing = opts.scf_mixing;

for it = 1:opts.scf_maxit
    MatCur = Mat0 + buildN(uPrev);
    MatCur = 0.5 * (MatCur + MatCur');
    [Ucand, Dcand] = solve_lowest_local(MatCur, M, nTrack, ops, opts);
    [uRaw, lambdaCand] = select_branch_by_overlap_local(Ucand, Dcand, uPrev, M);
    uRaw = align_phase_local(fixPhase(uRaw), uPrev, M);
    uRaw = normalize_in_M_local(uRaw, M);

    uMix = scfMixing * uRaw + (1 - scfMixing) * uPrev;
    uMix = align_phase_local(uMix, uPrev, M);
    uMix = normalize_in_M_local(uMix, M);

    abslambda = abs(lambdaCand - lambdaPrev);
    lambdaHist(it) = lambdaCand;
    info.n_iters = it;
    info.abslambda = abslambda;

    if it >= 3
        if abs(lambdaHist(it) - lambdaHist(it - 2)) < 5e-5 && ...
                abs(lambdaHist(it) - lambdaHist(it - 1)) > 2e-4
            scfMixing = max(0.02, 0.5 * scfMixing);
        end
    end

    uPrev = uMix;
    lambdaPrev = lambdaCand;

    if abslambda < opts.scf_tol_lambda
        info.converged = true;
        break;
    end
end

MatFinal = Mat0 + buildN(uPrev);
MatFinal = 0.5 * (MatFinal + MatFinal');
[Uf, Df] = solve_lowest_local(MatFinal, M, nTrack, ops, opts);
[u1, lambda1] = select_branch_by_overlap_local(Uf, Df, uPrev, M);
u1 = normalize_in_M_local(align_phase_local(fixPhase(u1), uPrev, M), M);
lambda1 = real(lambda1);
end

function [uh, D] = solve_lowest_local(Mat, M, nEigs, ops, opts)
% Compute the lowest generalized eigenpairs.
d = abs(diag(Mat));
d(d < opts.eps_diag) = 1;
Pfun = @(X) bsxfun(@rdivide, X, d);
[uh, D] = primme_eigs(Mat, M, nEigs, 'SA', ops, opts.primme_method, Pfun);
[lam, perm] = sort(real(diag(D)), 'ascend');
uh = uh(:, perm);
D = diag(lam);
end

function [uSel, lambdaSel] = select_branch_by_overlap_local(Ucand, Dcand, uPrev, M)
% Select the eigenvector branch with the largest mass overlap.
lams = real(diag(Dcand));
overlaps = zeros(size(Ucand, 2), 1);
for j = 1:size(Ucand, 2)
    overlaps(j) = abs(uPrev' * M * Ucand(:, j));
end
[~, idx] = max(overlaps);
uSel = Ucand(:, idx);
lambdaSel = lams(idx);
end

function [Mat, M, k_pw] = build_full_pw_problem_cached_local(L, Nc, pVr, nPwVr, opts)
% Load or assemble the cached full-domain plane-wave matrices.
cacheFile = fullfile(opts.cacheRoot, sprintf( ...
    'PW_FULL_EX3_L_%g_Nc_%d_NVr_%d_m_%d.mat', L, Nc, nPwVr, opts.scf_pw_grid_m));
if opts.use_pw_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'Mat', 'M', 'k_pw');
    Mat = S.Mat;
    M = S.M;
    k_pw = S.k_pw;
    return;
end

[k_pw, nBasis] = build_pw_disk_local(Nc);
Vfft = sample_full_potential_fft_local(L, pVr, nPwVr, opts.scf_pw_grid_m);
alpha = 2 * pi / L;
M = speye(nBasis);
Mat = complex(zeros(nBasis, nBasis));
for ii = 1:nBasis
    kii = k_pw(ii, :);
    kin = 0.5 * alpha ^ 2 * (kii * kii');
    for jj = 1:nBasis
        dk = k_pw(jj, :) - kii;
        if all(dk == 0)
            Mat(ii, jj) = kin;
        end
        idx = mod(dk, opts.scf_pw_grid_m) + 1;
        phase = exp(1i * pi * dk(1) * (1 / opts.scf_pw_grid_m + 1)) * ...
            exp(1i * pi * dk(2) * (1 / opts.scf_pw_grid_m + 1));
        Mat(ii, jj) = Mat(ii, jj) + Vfft(idx(1), idx(2)) * phase;
    end
end
Mat = real(0.5 * (Mat + Mat'));
save(cacheFile, 'Mat', 'M', 'k_pw', '-v7.3');
end

function [Mat, M] = build_full_iga_problem_cached_local( ...
nurbsBase, nurbs_refine, nElem, t, pVr, nPwVr, L, opts)
% Load or assemble the cached full-domain IGA matrices.
cacheFile = fullfile(opts.cacheNurbsRoot, sprintf( ...
    'NURBS_FULL_EX3_nElem_%d_t_%d_nd_%d_NVr_%d_ngp_%d.mat', ...
    nElem, t, nurbs_refine.n_dofs_domains, nPwVr, opts.n_gp));
if opts.use_nurbs_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'Mat', 'M');
    Mat = S.Mat;
    M = S.M;
    return;
end
[Mat, M] = generate_A_M_NURBS_2D( ...
    nurbsBase, nurbs_refine, pVr, nPwVr, L, opts.n_gp, 'Example_3');
save(cacheFile, 'Mat', 'M', '-v7.3');
end

function Vfft = sample_full_potential_fft_local(L, pVr, nPwVr, sampleM)
% Sample the external potential on the FFT grid.
dx = L / sampleM;
x = -L / 2 + dx / 2 + (0:sampleM - 1) * dx;
[X, Y] = ndgrid(x, x);
V = arrayfun(@(xv, yv) Vr_2D_Example_1(pVr, L, nPwVr, xv, yv), X, Y);
Vfft = ifftn(V);
end

function rho = sample_rho_pw_full_grid_local(cPw, kPw, L, m)
% Evaluate the plane-wave density on the full periodic grid.
omegaArea = L * L;
C = complex(zeros(m, m));
for s = 1:size(kPw, 1)
    k1 = kPw(s, 1);
    k2 = kPw(s, 2);
    i1 = mod(k1, m) + 1;
    i2 = mod(k2, m) + 1;
    phaseMid = (-1) ^ (k1 + k2) * exp(1i * pi * (k1 + k2) / m);
    C(i1, i2) = C(i1, i2) + cPw(s) * phaseMid;
end
uGrid = (m * m / sqrt(omegaArea)) * ifftn(C);
rho = real(abs(uGrid) .^ 2);
end

function Ugrid = evaluate_run_on_grid_local(run, kind, Xg, Yg)
% Evaluate a saved method solution on the common grid.
L = 4;
coeff = run.uh(:, 1);
switch lower(kind)
    case 'pw'
        Ugrid = reshape(pw_eval_val_chunked_local(coeff, run.k_pw, Xg(:), Yg(:), L, 20000), size(Xg));
    case 'iga'
        Ugrid = reshape(iga_eval_on_one_patch_local( ...
            run.nurbs_refine, coeff, Xg(:), Yg(:), 0.0, 0.0, 2.0), size(Xg));
    case 'igapw'
        a = 0.2;
        nI = run.n_dofs_nurbs;
        coeffI = coeff(1:nI);
        coeffPw = coeff(nI + (1:run.n_pw_basis));
        x = Xg(:);
        y = Yg(:);
        maskInner = abs(x) <= a & abs(y) <= a;
        val = zeros(size(x));
        if any(~maskInner)
            val(~maskInner) = pw_eval_val_chunked_local( ...
                coeffPw, run.k_pw, x(~maskInner), y(~maskInner), L, 20000);
        end
        if any(maskInner)
            val(maskInner) = iga_eval_on_one_patch_local( ...
                run.nurbs_refine, coeffI, x(maskInner), y(maskInner), 0.0, 0.0, a);
        end
        Ugrid = reshape(val, size(Xg));
    otherwise
        error('Unknown field kind: %s', kind);
end
end

function val = pw_eval_val_chunked_local(coeff, pVec, X, Y, L, chunkSize)
% Evaluate plane-wave values in point chunks.
nPts = numel(X);
val = zeros(nPts, 1);
for i1 = 1:chunkSize:nPts
    i2 = min(i1 + chunkSize - 1, nPts);
    F = [X(i1:i2).'; Y(i1:i2).'];
    expo = exp((1i * 2 * pi / L) * (pVec * F));
    val(i1:i2) = ((coeff.' * expo) / L).';
end
end

function val = iga_eval_on_one_patch_local(nurbs, coeff, X, Y, xc, yc, a)
% Evaluate an IGA field on one rectangular patch.
pu = nurbs.pu;
pv = nurbs.pv;
U = nurbs.Ubar(:).';
V = nurbs.Vbar(:).';
mU = length(U) - pu - 1;
nV = length(V) - pv - 1;
val = zeros(numel(X), 1);
for k = 1:numel(X)
    u = max(0, min(1, (X(k) - (xc - a)) / (2 * a)));
    v = max(0, min(1, (Y(k) - (yc - a)) / (2 * a)));
    spanU = findspan_local(mU - 1, pu, u, U);
    spanV = findspan_local(nV - 1, pv, v, V);
    Nu = bspline_basis_local(U, pu, u, spanU);
    Nv = bspline_basis_local(V, pv, v, spanV);
    s = 0.0 + 0.0i;
    for j1 = (spanV - pv):spanV
        lv = j1 - (spanV - pv) + 1;
        for i1 = (spanU - pu):spanU
            lu = i1 - (spanU - pu) + 1;
            row = i1 + (j1 - 1) * mU;
            s = s + coeff(row) * Nu(lu) * Nv(lv);
        end
    end
    val(k) = s;
end
end

function [MatRed, MRed, fullToReduced, R] = apply_periodic_reduction_local(MatFull, MFull, nurbs)
% Reduce the full IGA matrices by periodic trace identification.
fullToReduced = build_periodic_trace_map_local(nurbs.m, nurbs.n);
nFull = numel(fullToReduced);
nRed = max(fullToReduced);
R = sparse((1:nFull).', fullToReduced(:), 1, nFull, nRed);
MatRed = R' * MatFull * R;
MRed = R' * MFull * R;
MatRed = 0.5 * (MatRed + MatRed');
MRed = 0.5 * (MRed + MRed');
end

function fullToReduced = build_periodic_trace_map_local(m, n)
% Map full IGA boundary degrees of freedom to periodic reduced degrees of freedom.
rep = reshape(1:(m * n), m, n);
rep(m, :) = rep(1, :);
rep(:, n) = rep(:, 1);
[~, ~, fullToReduced] = unique(rep(:), 'stable');
fullToReduced = double(fullToReduced(:));
end

function uFull = expand_periodic_coeff_local(uRed, fullToReduced)
% Expand reduced periodic coefficients to the full IGA basis.
uFull = uRed(fullToReduced, :);
end

function [Uout, meta] = normalize_and_align_local(Uin, Uref, dx)
% Normalize a grid field and align its phase.
u = Uin;
normBefore = sqrt(max(sum(abs(u(:)).^2) * dx * dx, eps));
u = u ./ normBefore;
phaseApplied = 0;
if isempty(Uref)
    [~, idx] = max(abs(u(:)));
    phaseApplied = angle(u(idx));
    u = u * exp(-1i * phaseApplied);
    if real(u(idx)) < 0
        u = -u;
        phaseApplied = phaseApplied + pi;
    end
else
    innerVal = sum(u(:) .* conj(Uref(:))) * dx * dx;
    if abs(innerVal) > 0
        phaseApplied = angle(innerVal);
        u = u * exp(-1i * phaseApplied);
    end
    if real(sum(u(:) .* conj(Uref(:))) * dx * dx) < 0
        u = -u;
        phaseApplied = phaseApplied + pi;
    end
end
Uout = u;
meta = struct('norm_before', normBefore, 'phase_applied', phaseApplied);
end

function val = l2_error_from_grid_local(fieldVals, dx)
% Compute the discrete L2 error from grid values.
val = sqrt(sum(abs(fieldVals(:)) .^ 2) * dx * dx);
end

function n = get_scf_iterations_local(run)
% Read the SCF iteration count from a saved run.
n = double(run.meta.scf_iterations);
end

function u = normalize_in_M_local(u, M)
% Normalize a vector in the mass inner product.
nu = real(u' * M * u);
u = u / sqrt(abs(nu));
end

function uAligned = align_phase_local(uAligned, uTarget, M)
% Align a vector phase with a target in the mass inner product.
alpha = uTarget' * M * uAligned;
if abs(alpha) > 0
    uAligned = exp(-1i * angle(alpha)) * uAligned;
end
end

function u = fix_global_phase_single_pw_local(u, kPw)
% Fix the global phase from the dominant plane-wave coefficient.
[~, idx] = max(abs(u));
k0 = kPw(idx, :);
idxNeg = find(kPw(:, 1) == -k0(1) & kPw(:, 2) == -k0(2), 1);
cK = u(idx);
if ~isempty(idxNeg) && idxNeg ~= idx
    cMk = u(idxNeg);
    theta = 0.5 * (angle(cK) + angle(cMk));
else
    theta = angle(cK);
end
u = exp(-1i * theta) * u;
if real(u(idx)) < 0
    u = -u;
end
end

function N = bspline_basis_local(U, p, u, span)
% Evaluate the nonzero B-spline basis functions.
ndu = zeros(p + 1, p + 1);
left = zeros(1, p + 1);
right = zeros(1, p + 1);
ndu(1, 1) = 1.0;
for j = 1:p
    left(j + 1) = u - U(span + 1 - j);
    right(j + 1) = U(span + j) - u;
    saved = 0.0;
    for r = 0:(j - 1)
        ndu(j + 1, r + 1) = right(r + 2) + left(j - r + 1);
        temp = ndu(r + 1, j) / ndu(j + 1, r + 1);
        ndu(r + 1, j + 1) = saved + right(r + 2) * temp;
        saved = left(j - r + 1) * temp;
    end
    ndu(j + 1, j + 1) = saved;
end
N = ndu(1:p + 1, p + 1).';
end

function span = findspan_local(n, p, u, U)
% Find the active B-spline knot span.
if u >= U(n + 2)
    span = n + 1;
    return;
end
if u <= U(p + 1)
    span = p + 1;
    return;
end
low = p + 1;
high = n + 2;
mid = floor((low + high) / 2);
while (u < U(mid) || u >= U(mid + 1))
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end
span = mid;
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

function nurbsBase = make_rect_patch_local(rect)
% Build a rectangular NURBS patch.
x1 = rect(1);
x2 = rect(2);
y1 = rect(3);
y2 = rect(4);
ConPts = zeros(2, 2, 2);
ConPts(:, :, 1) = [x1 x1; x2 x2];
ConPts(:, :, 2) = [y1 y2; y1 y2];
nurbsBase = struct();
nurbsBase.ConPts = ConPts;
nurbsBase.weights = [1 1; 1 1];
end

function ensure_dir_local(pathName)
% Create a directory when it does not exist.
if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end
