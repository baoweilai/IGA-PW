function [result, meta] = solve_full_domain(mode, cfg, opts)
% Solve one full-domain IGA or plane-wave problem.

% Read the method and eigensolver settings.
format long;
t_total = tic;

mode = lower(string(mode));
L = 4;
Example = 'Example_2';
n_eigenvalues = get_opt_local(opts, 'n_eigenvalues', 2);
solve_mode = lower(string(get_opt_local(opts, 'solve_mode', 'none')));

ops = struct();
ops.tol = get_opt_local(opts, 'primme_tol', 1e-12);
ops.maxit = get_opt_local(opts, 'primme_maxit', 2e8);
ops.reportLevel = get_opt_local(opts, 'primme_reportLevel', 0);
primme_method = get_opt_local(opts, 'primme_method', 'DEFAULT_MIN_TIME');

% Assemble the selected full-domain matrix pair.
switch mode
    case "pw"
        [Mat, M, static_meta] = build_pw_problem_local(cfg, opts, L);
    case "iga"
        [Mat, M, static_meta] = build_iga_problem_local(cfg, opts, L, Example);
    otherwise
        error('Unsupported mode: %s', mode);
end

[uh_solve, D, rnorms, ~] = solve_generalized_eigs_local( ...
    Mat, M, solve_mode, n_eigenvalues, ops, primme_method, opts);
lambda = sort(real(diag(D)), 'ascend');
[~, perm] = sort(real(diag(D)), 'ascend');
uh_solve = uh_solve(:, perm);

for i = 1:n_eigenvalues
    ni = real(uh_solve(:, i)' * M * uh_solve(:, i));
    uh_solve(:, i) = uh_solve(:, i) ./ sqrt(max(abs(ni), eps));
end

uh = uh_solve;
if mode == "iga" && isfield(static_meta, 'full_to_reduced') && ~isempty(static_meta.full_to_reduced)
    uh = expand_periodic_coeff_local(uh_solve, static_meta.full_to_reduced);
end

if mode == "pw"
    uh = fix_global_phase_pw_local(uh, static_meta.k_pw);
end

% Package the normalized eigenpairs and solver data.
meta = static_meta;
meta.mode = char(mode);
meta.solve_mode = char(solve_mode);
meta.prec_type = solve_mode_to_label_local(solve_mode);
meta.time_total = toc(t_total);
meta.primme_rnorms = rnorms;

result = struct();
result.lambda = lambda(:).';
result.uh = uh;
result.n_dofs_total = size(Mat, 1);
result.n_dofs_eval = size(uh, 1);

% Save the retained full-domain result.
if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end
    run = struct();
    run.mode = char(mode);
    run.lambda = result.lambda;
    run.uh = uh;
    run.n_dofs_total = result.n_dofs_total;
    run.n_dofs_eval = result.n_dofs_eval;
    run.meta = meta;
    if mode == "pw"
        run.k_pw = static_meta.k_pw;
    else
        run.nurbs_refine = static_meta.nurbs_refine;
        if isfield(static_meta, 'full_to_reduced')
            run.uh_reduced = uh_solve;
        end
    end
    save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
end

end

function [Mat, M, meta] = build_pw_problem_local(cfg, opts, L)
% Assemble the full-domain plane-wave eigenproblem.
Nc = cfg.Nc;
N_Vr = get_opt_local(opts, 'potential_Nc_pw', 2);
sample_m = get_opt_local(opts, 'pw_potential_grid_m', 2048);

[k_Vr, n_pw_Vr] = build_pw_disk_local(N_Vr);
[Mat, M, k_pw] = get_full_pw_matrices_cached_local(L, Nc, k_Vr, n_pw_Vr, sample_m, opts);

Mat = 0.5 * (Mat + Mat');
M = 0.5 * (M + M');

meta = struct();
meta.Nc = Nc;
meta.potential_Nc = N_Vr;
meta.L = L;
meta.k_pw = k_pw;
meta.periodic_bc = true;
end

function [Mat, M, meta] = build_iga_problem_local(cfg, opts, L, Example)
% Build the full-domain IGA matrix pair.
% Construct the elevated and refined rectangular patch.
pdeg = cfg.pdeg;
refine = cfg.refine;
t = pdeg - 1;
pu0 = 1;
pv0 = 1;

[knotUe, knotVe] = IGADegreeElevSurface([0 0 1 1], [0 0 1 1], t);
pu = pu0 + t;
pv = pv0 + t;

nurbs_original = make_rect_patch_local([-L/2, L/2, -L/2, L/2]);
nurbs_original.pu = pu0;
nurbs_original.pv = pv0;
nurbs_original.knotU = [0 0 1 1];
nurbs_original.knotV = [0 0 1 1];

nurbs_refine = IGA_2D_Grid(knotUe, knotVe, pu, pv, refine);

N_Vr = get_opt_local(opts, 'potential_Nc_iga', 2);
[k_Vr, n_pw_Vr] = build_pw_disk_local(N_Vr);
n_gp = get_opt_local(opts, 'n_gp', 10);
cache_ok = isfield(opts, 'cacheNurbsRoot') && ~isempty(opts.cacheNurbsRoot);
if cache_ok
    cacheFile = fullfile(opts.cacheNurbsRoot, sprintf( ...
        'NURBS_EX2_FULL_nu_%d_nv_%d_pu_%d_pv_%d_nd_%d_NVr_%d_ngp_%d.mat', ...
        numel(nurbs_refine.Ubar), numel(nurbs_refine.Vbar), ...
        nurbs_refine.pu, nurbs_refine.pv, nurbs_refine.n_dofs_domains, n_pw_Vr, n_gp));
else
    cacheFile = '';
end

% Load or assemble the NURBS matrices.
if cache_ok && exist(cacheFile, 'file')
    S = load(cacheFile, 'Mat', 'M');
    Mat = S.Mat;
    M = S.M;
else
    [Mat, M] = generate_A_M_NURBS_2D(nurbs_original, nurbs_refine, k_Vr, n_pw_Vr, L, n_gp, Example);
    if cache_ok
        if ~exist(opts.cacheNurbsRoot, 'dir'), mkdir(opts.cacheNurbsRoot); end
        save(cacheFile, 'Mat', 'M', '-v7.3');
    end
end

Mat = 0.5 * (Mat + Mat');
M = 0.5 * (M + M');

% Apply periodic identification and package the mesh data.
[Mat, M, full_to_reduced] = ...
    apply_periodic_reduction_local(Mat, M, nurbs_refine);

meta = struct();
meta.pdeg = pdeg;
meta.refine = refine;
meta.h = L / (2^refine);
meta.L = L;
meta.nurbs_refine = nurbs_refine;
meta.periodic_bc = true;
meta.n_dofs_full = nurbs_refine.n_dofs_domains;
meta.n_dofs_reduced = size(Mat, 1);
meta.full_to_reduced = full_to_reduced;
end

function [Mat, M, k_pw] = get_full_pw_matrices_cached_local(L, Nc, p_Vr, n_pw_Vr, sample_m, opts)
% Load or assemble cached full-domain plane-wave matrices.
cache_ok = get_opt_local(opts, 'use_pw_cache', true) ...
    && isfield(opts, 'cacheRoot') && ~isempty(opts.cacheRoot);

if cache_ok
    cacheFile = fullfile(opts.cacheRoot, sprintf('PW_FULL_EX2_L_%g_Nc_%d_NVr_%d_m_%d.mat', L, Nc, n_pw_Vr, sample_m));
    if exist(cacheFile, 'file')
        S = load(cacheFile, 'Mat', 'M', 'k_pw');
        Mat = S.Mat;
        M = S.M;
        k_pw = S.k_pw;
        return;
    end
end

[k_pw, n_basis] = build_pw_disk_local(Nc);
Vfft = sample_full_potential_fft_local(L, p_Vr, sample_m);

M = speye(n_basis);
Mat = zeros(n_basis, n_basis);

for ii = 1:n_basis
    kii = k_pw(ii, :);
    kin = 0.5 * (kii * kii') * (2 * pi / L)^2;
    for jj = 1:n_basis
        dk = k_pw(jj, :) - kii;
        if all(dk == 0)
            Mat(ii, jj) = kin;
        end
        dk_idx = mod(dk, sample_m) + 1;
        % meshgrid stores x along columns and y along rows, so the FFT
        % The coefficient lookup uses (y-index, x-index).
        Vcoef = Vfft(dk_idx(2), dk_idx(1)) ...
            * exp(1i * pi * dk(1) * (1 / sample_m + 1)) ...
            * exp(1i * pi * dk(2) * (1 / sample_m + 1));
        Mat(ii, jj) = Mat(ii, jj) + Vcoef;
    end
end

Mat = sparse(real(Mat));
if cache_ok
    if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
    save(cacheFile, 'Mat', 'M', 'k_pw', '-v7.3');
end
end

function Vfft = sample_full_potential_fft_local(L, p_Vr, sample_m)
% Sample full potential FFT.
dx = L / sample_m;
x = -L / 2 + dx / 2 : dx : L / 2 - dx / 2;
[X, Y] = meshgrid(x, x);

V = example2_potential_grid_local(p_Vr, L, X, Y);
Vfft = ifftn(V);
end

function Vr = example2_potential_grid_local(p, L, X, Y)
% Compute potential grid.
alpha = 5;
p = p * 2 * pi / L;

r1x = X + 1;
r1y = Y;
r2x = X - 1;
r2y = Y;
r1 = hypot(r1x, r1y);
r2 = hypot(r2x, r2y);

s1 = zeros(size(X));
s2 = zeros(size(X));
for i = 1:size(p, 1)
    p_norm = norm(p(i, :));
    if p_norm > 0
        fac = erfc(p_norm / (2 * alpha)) / p_norm;
        s1 = s1 + fac * exp(1i * (p(i, 1) * r1x + p(i, 2) * r1y));
        s2 = s2 + fac * exp(1i * (p(i, 1) * r2x + p(i, 2) * r2y));
    end
end

Vr1 = -2 * pi * s1 / (L * L) - erfc(alpha * r1) ./ r1 + 2 * alpha / sqrt(pi);
Vr2 = -2 * pi * s2 / (L * L) - erfc(alpha * r2) ./ r2 + 2 * alpha / sqrt(pi);
Vr = real(Vr1 + Vr2);
end

function [uh, D, rnorms, stats] = solve_generalized_eigs_local(Mat, M, solve_mode, n_eigs, ops, method, opts)
% Solve the selected generalized eigenproblem.
Pfun = build_selected_preconditioner_local(Mat, solve_mode, opts);
if isempty(Pfun)
    [uh, D, rnorms, stats] = call_primme_noprec_local(Mat, M, n_eigs, ops, method);
else
    [uh, D, rnorms, stats] = call_primme_with_prec_local(Mat, M, n_eigs, ops, method, Pfun);
end
end

function uh = fix_global_phase_pw_local(uh, k_pw)
% Fix the global phase using the zero plane-wave coefficient.
for i = 1:size(uh, 2)
    [~, idx] = max(abs(uh(:, i)));
    k0 = k_pw(idx, :);
    idx_neg = find(k_pw(:, 1) == -k0(1) & k_pw(:, 2) == -k0(2), 1);
    ck = uh(idx, i);

    if ~isempty(idx_neg) && idx_neg ~= idx
        cmk = uh(idx_neg, i);
        theta = 0.5 * (angle(ck) + angle(cmk));
        uh(:, i) = exp(-1i * theta) * uh(:, i);
    else
        theta = angle(ck);
        uh(:, i) = exp(-1i * theta) * uh(:, i);
        if real(uh(idx, i)) < 0
            uh(:, i) = -uh(:, i);
        end
    end
end
end

function nurbs_original = make_rect_patch_local(rect)
% Build an affine rectangular NURBS patch.
x1 = rect(1); x2 = rect(2);
y1 = rect(3); y2 = rect(4);

ConPts = zeros(2, 2, 2);
ConPts(:, :, 1) = [x1 x1; x2 x2];
ConPts(:, :, 2) = [y1 y2; y1 y2];

nurbs_original = struct();
nurbs_original.ConPts = ConPts;
nurbs_original.weights = [1 1; 1 1];
end

function [k_list, n_basis] = build_pw_disk_local(Nc)
% Build the circular plane-wave index set.
N = floor(Nc);
k_list = zeros((2 * N + 1)^2, 2);
n_basis = 0;
for k1 = -N:N
    m = floor(sqrt(N^2 - k1^2));
    for k2 = -m:m
        n_basis = n_basis + 1;
        k_list(n_basis, :) = [k1, k2];
    end
end
k_list = k_list(1:n_basis, :);
end

function [MatRed, MRed, full_to_reduced] = apply_periodic_reduction_local(MatFull, MFull, nurbs_refine)
% Apply periodic reduction.
full_to_reduced = build_periodic_trace_map_local(nurbs_refine.m, nurbs_refine.n);
nFull = numel(full_to_reduced);
nRed = max(full_to_reduced);
R = sparse((1:nFull).', full_to_reduced(:), 1, nFull, nRed);
MatRed = R' * MatFull * R;
MRed = R' * MFull * R;
MatRed = 0.5 * (MatRed + MatRed');
MRed = 0.5 * (MRed + MRed');
end

function full_to_reduced = build_periodic_trace_map_local(m, n)
% Build periodic trace map.
rep = reshape(1:(m * n), m, n);
rep(m, :) = rep(1, :);
rep(:, n) = rep(:, 1);
[~, ~, full_to_reduced] = unique(rep(:), 'stable');
full_to_reduced = double(full_to_reduced(:));
end

function uh_full = expand_periodic_coeff_local(uh_red, full_to_reduced)
% Expand reduced periodic coefficients to the full basis.
uh_full = uh_red(full_to_reduced, :);
end

function [uh, D, rnorms, stats] = call_primme_noprec_local(Mat, M, n_eigs, ops, method)
% Solve the generalized eigenproblem with PRIMME without preconditioning.
if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, method);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops);
end
end

function [uh, D, rnorms, stats] = call_primme_with_prec_local(Mat, M, n_eigs, ops, method, Pfun)
% Solve with PRIMME using the supplied preconditioner.
if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, method, Pfun);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, [], Pfun);
end
end

function Pfun = build_selected_preconditioner_local(Ktau, solve_mode, opts)
% Build selected preconditioner.
solve_mode = lower(string(solve_mode));
switch solve_mode
    case "none"
        Pfun = [];
    case "purediag"
        eps_diag = get_opt_local(opts, 'eps_diag', 1e-12);
        d = abs(diag(Ktau));
        d(d < eps_diag) = 1;
        dinv = 1 ./ d;
        Pfun = @(X) bsxfun(@times, X, dinv);
    otherwise
        error('Unsupported solve_mode for full-domain solver: %s', solve_mode);
end
end

function label = solve_mode_to_label_local(solve_mode)
% Convert the solver-mode identifier to a display label.
switch lower(string(solve_mode))
    case "none"
        label = 'None';
    case "purediag"
        label = 'PureDiag-Jacobi';
    otherwise
        label = char(solve_mode);
end
end

function value = get_opt_local(opts, name, default)
% Read one option with a default value.
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = default;
end
end
