function [result, n_dofs_total, meta] = solve_iga_pw_dg(Refinement, t, Nc, n_eigenvalues, opts)
% Solve the Example 1 IGA-PW-DG problem.

format long;

% Read the discretization and eigensolver settings.
beta = opts.beta;
lambda_ref = opts.lambda_ref;
tau_shift = opts.tau_shift;
eps_diag = opts.eps_diag;
iface_reg = opts.iface_reg;

L = opts.L;
a = opts.a;
inner_domains_coordinates = [-a, a, -a, a];

% Build the inner NURBS patch and plane-wave bases.
nu = 2;
nv = 2;
ConPts = zeros(nu, nv, 2);
x = [-a, a];
y = [-a, a];
ConPts(:, :, 1) = [x(1) x(1); x(2) x(2)];
ConPts(:, :, 2) = [y(1) y(2); y(1) y(2)];

weights = [1 1; 1 1];
knotU = [0 0 1 1];
knotV = [0 0 1 1];

pu0 = 1;
pv0 = 1;
nurbs_original = struct();
nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.pu = pu0;
nurbs_original.pv = pv0;
nurbs_original.knotU = knotU;
nurbs_original.knotV = knotV;

[knotUe, knotVe] = IGADegreeElevSurface(knotU, knotV, t);
pu = pu0 + t;
pv = pv0 + t;

nurbs_refine = IGA_2D_Grid(knotUe, knotVe, pu, pv, Refinement);
h = 2 * a / (2 ^ Refinement);
n_dofs_nurbs = nurbs_refine.n_dofs_domains;

[k_pw, n_pw_basis] = build_pw_disk_local(Nc);
pw_dofs_indices = n_dofs_nurbs + (1:n_pw_basis);
n_dofs_total = n_dofs_nurbs + n_pw_basis;

N_Vr = 2;
[k_Vr, n_pw_Vr] = build_pw_disk_local(N_Vr);

idxI = 1:n_dofs_nurbs;
idxP = pw_dofs_indices;
use_square_dg_fast = logical(opts.use_square_dg_fast);

H = sparse(n_dofs_total, n_dofs_total);
M = sparse(n_dofs_total, n_dofs_total);

% Assemble the uncoupled NURBS and plane-wave volume blocks.
t_nurbs = tic;
[H_nurbs, M_nurbs, nurbs_build_meta] = generate_A_M_NURBS_2D( ...
    nurbs_original, nurbs_refine, k_Vr, n_pw_Vr, L, opts.n_gp, opts.Example, opts);
time_nurbs_build = toc(t_nurbs);

H(idxI, idxI) = H_nurbs;
M(idxI, idxI) = M_nurbs;

t_pw = tic;
[H_pw, M_pw, pw_timing, pw_cache] = get_pw_matrices_cached( ...
    L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);
time_pw_build = toc(t_pw);

H(idxP, idxP) = H_pw;
M(idxP, idxP) = M_pw;

% Assemble the DG interface terms and shifted eigenproblem.
t_dg = tic;
p = k_pw;
if use_square_dg_fast
    [P, S, dg_fast_meta] = assemble_DG_square_interface_fast( ...
        nurbs_refine, p, pw_dofs_indices, L, a, n_dofs_total);
    dg_assembly_method = 'square_fast';
else
    [P_Bottom, S_Bottom] = IGA_DG_Bottom_Edge_Assemble( ...
        nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
    [P_Top, S_Top] = IGA_DG_Top_Edge_Assemble( ...
        nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
    [P_Left, S_Left] = IGA_DG_Left_Edge_Assemble( ...
        nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
    [P_Right, S_Right] = IGA_DG_Right_Edge_Assemble( ...
        nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
    P = P_Bottom + P_Top + P_Left + P_Right;
    S = S_Bottom + S_Top + S_Left + S_Right;
    dg_fast_meta = struct();
    dg_assembly_method = 'legacy_edge_loops';
end

sigma = beta * (1 / h + Nc);
diagP = full(diag(P));

Mat = H - 0.5 * S - 0.5 * S' + sigma * P;
Mat = 0.5 * (Mat + Mat');
M = 0.5 * (M + M');
Atau = Mat - tau_shift * M;
time_dg_assembly = toc(t_dg);

time_shared_core = time_nurbs_build + time_pw_build + time_dg_assembly;
time_shared_build = time_shared_core;

% Package assembly data and configure the eigensolver.
meta = struct();
meta.Refinement = Refinement;
meta.Nc = Nc;
meta.L = L;
meta.a = a;
meta.h = h;
meta.beta = beta;
meta.sigma = sigma;
meta.lambda_ref = lambda_ref;
meta.n_dofs_nurbs = n_dofs_nurbs;
meta.n_pw_basis = n_pw_basis;
meta.time_nurbs_build = time_nurbs_build;
meta.time_pw_build = time_pw_build;
meta.time_dg_assembly = time_dg_assembly;
meta.time_shared_core = time_shared_core;
meta.time_shared_build = time_shared_build;
meta.enabled_variants = {'InterfaceBlock'};
meta.use_square_dg_fast = use_square_dg_fast;
meta.dg_assembly_method = dg_assembly_method;
meta.dg_fast_meta = dg_fast_meta;
meta.nurbs_assembly_method = nurbs_build_meta.method;
meta.nurbs_build_meta = nurbs_build_meta;
meta.pw_timing = pw_timing;
meta.pw_cache = pw_cache;

ops = struct();
ops.tol = opts.primme_tol;
ops.maxit = opts.primme_maxit;
ops.reportLevel = opts.primme_reportLevel;

primme_method = opts.primme_method;
fprintf('[PRIMME] method = %s\n', string_or_empty_local(primme_method));

primme_target = opts.primme_target;
fprintf('[PRIMME] target = %s\n', target_to_string_local(primme_target));

seed = 20260404;
rng(seed, 'twister');
v0 = randn(n_dofs_total, 1) + 1i * randn(n_dofs_total, 1);
v0 = v0 / max(norm(v0), eps);
ops.v0 = v0;
meta.rng_seed = seed;
fprintf('[PRIMME] shared v0 seed             = %d\n', seed);

result = struct();

% Build the interface preconditioner and solve the retained eigenpairs.
fprintf('\n---------------- InterfaceBlock (PRIMME) -------------------\n');

[bd_ifb, ifb, time_build_ifb] = ...
    build_interface_block_prec_data_local( ...
    Atau, P, H_pw, M_pw, sigma, tau_shift, diagP, idxI, idxP, ...
    eps_diag, iface_reg);

res = struct();
res.prec_type = 'InterfaceBlock';
res.time_build_prec = time_build_ifb;

precfun = @(x) apply_interface_block_prec_local(x, bd_ifb, ifb);

t_step = tic;
[uh, D, rnorms, ~] = call_primme_with_prec_local( ...
    Mat, M, n_eigenvalues, primme_target, ops, primme_method, precfun);
res.time_primme = toc(t_step);
res.time_total = res.time_build_prec + res.time_primme;

[res.lambda, res.uh] = postprocess_eigs_local( ...
    uh, D, M, n_eigenvalues, n_dofs_nurbs, k_pw);
res.err1 = abs(res.lambda(1) - lambda_ref);
res.rnorms = rnorms;

fprintf('[RESULT] InterfaceBlock: lambda1 = %.12f, err1 = %.3e, total = %.6f\n', ...
    res.lambda(1), res.err1, res.time_total);

result.InterfaceBlock = res;

meta.time_total = time_shared_build + sum_variant_times_local(result);

% Save the retained eigenpairs and optional operators.
if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir')
        mkdir(opts.outDir);
    end

    run = struct();
    run.result = result;
    run.meta = meta;
    run.n_dofs_total = n_dofs_total;
    run.n_dofs_nurbs = n_dofs_nurbs;
    run.n_pw_basis = n_pw_basis;
    run.pw_timing = pw_timing;
    run.pw_cache = pw_cache;

    if isfield(opts, 'save_matrices') && opts.save_matrices
        run.M = M;
        if isfield(opts, 'save_mat') && opts.save_mat
            run.Mat = Mat;
        end
        run.P = P;
        run.S = S;
    end

    if isfield(opts, 'save_pw_index') && opts.save_pw_index
        run.k_pw = k_pw;
    end

    if isfield(opts, 'save_nurbs') && opts.save_nurbs
        run.nurbs_original = nurbs_original;
        run.nurbs_refine = nurbs_refine;
    end

    if isfield(opts, 'save_eigenvectors') && opts.save_eigenvectors
        save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
    else
        save(fullfile(opts.outDir, 'run.mat'), 'run');
    end
end
end

function [bd, time_build_prec] = build_blockdiag_data_local( ...
Atau, H_pw, M_pw, sigma, tau_shift, diagP, idxI, idxP, eps_diag)
% Build the diagonal data used by the interface block.
t_build = tic;

dII = abs(diag(Atau(idxI, idxI)));
dII(abs(dII) < eps_diag) = 1;

dPP = abs(diag(H_pw) + sigma * diagP(idxP) - tau_shift * diag(M_pw));
dPP(abs(dPP) < eps_diag) = 1;

bd = struct();
bd.idxI = idxI;
bd.idxP = idxP;
bd.dII = dII;
bd.dPP = dPP;
bd.dall = [dII; dPP];
bd.dinv = 1 ./ bd.dall;

time_build_prec = toc(t_build);
end

function [bd, ifb, time_build_prec] = build_interface_block_prec_data_local( ...
    Atau, P, H_pw, M_pw, sigma, tau_shift, diagP, idxI, idxP, ...
    eps_diag, iface_reg)
% Build the TB-DG interface-block preconditioner.
[bd, time_build_bd] = build_blockdiag_data_local( ...
    Atau, H_pw, M_pw, sigma, tau_shift, diagP, idxI, idxP, eps_diag);

t_build = tic;
gamma = find(sum(abs(P), 2) ~= 0);
eta = setdiff((1:size(Atau, 1))', gamma);
Ag = 0.5 * (Atau(gamma, gamma) + Atau(gamma, gamma)');
delta = iface_reg * max(1, norm(Ag, 1));
Dg = decomposition(Ag + delta * speye(size(Ag)), 'chol');

ifb = struct();
ifb.gamma = gamma;
ifb.eta = eta;
ifb.solve = Dg;
ifb.Ag = Ag;
ifb.delta = delta;

time_build_prec = time_build_bd + toc(t_build);
end

function total = sum_variant_times_local(result)
% Sum the recorded timing components for one solver variant.
total = 0;
fn = fieldnames(result);
for i = 1:numel(fn)
    entry = result.(fn{i});
    if isstruct(entry) && isfield(entry, 'time_total')
        total = total + entry.time_total;
    end
end
end

function [k_list, n_basis] = build_pw_disk_local(Nc)
% Build the circular plane-wave index set.
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

function [H_pw, M_pw, pw_timing, cache_meta] = get_pw_matrices_cached( ...
L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts)
% Load or assemble cached plane-wave matrices.
cache_ok = isfield(opts, 'use_pw_cache') && opts.use_pw_cache ...
    && isfield(opts, 'cacheRoot') && ~isempty(opts.cacheRoot);

cache_meta = struct();
cache_meta.cache_ok = cache_ok;
cache_meta.cache_hit = false;
cache_meta.cache_file = '';
cache_meta.cache_key = '';
cache_meta.cache_version = 'PW2D_KERNEL_CHEB_V1';
cache_meta.t_cache_load = 0;

if cache_ok
    cache_key = build_pw_cache_key_local(L, Nc, n_pw_Vr, opts, cache_meta.cache_version);
    cache_meta.cache_key = cache_key;
    cache_meta.cache_file = fullfile(opts.cacheRoot, [cache_key '.mat']);

    if exist(cache_meta.cache_file, 'file')
        t_load = tic;
        S = load(cache_meta.cache_file, 'H_pw', 'M_pw', 'pw_timing', 'pw_meta');
        cache_meta.t_cache_load = toc(t_load);
        cache_meta.cache_hit = true;
        H_pw = S.H_pw;
        M_pw = S.M_pw;
        pw_timing = S.pw_timing;
        pw_timing.cache_hit = true;
        pw_timing.t_cache_load = cache_meta.t_cache_load;
        if isfield(S, 'pw_meta')
            cache_meta = merge_struct_local(cache_meta, S.pw_meta);
        end
        return;
    end
end

[H_pw, M_pw, ~, pw_timing] = generate_A_M_PW_2D( ...
    L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);
pw_timing.cache_hit = false;
pw_timing.t_cache_load = 0;

if cache_ok
    if ~exist(opts.cacheRoot, 'dir')
        mkdir(opts.cacheRoot);
    end
    pw_meta = cache_meta;
    save(cache_meta.cache_file, 'H_pw', 'M_pw', 'pw_timing', 'pw_meta', '-v7.3');
end
end

function key = build_pw_cache_key_local(L, Nc, n_pw_Vr, opts, version)
% Build a cache key from the plane-wave discretization settings.
inner_cheb_n = opts.inner_cheb_n;
fft_grid_n = opts.pw_fft_grid_n;

key = sprintf('%s_L_%s_Nc_%d_NVr_%d_cheb_%d_fft_%d', ...
    version, sanitize_number_local(L), Nc, n_pw_Vr, inner_cheb_n, fft_grid_n);
end

function txt = sanitize_number_local(x)
% Convert a numeric value to a path tag.
txt = strrep(num2str(x, '%.6g'), '.', 'p');
txt = strrep(txt, '-', 'm');
end

function out = merge_struct_local(base, override)
% Merge override fields into the base structure.
out = base;
if ~isstruct(override)
    return;
end

fn = fieldnames(override);
for i = 1:numel(fn)
    out.(fn{i}) = override.(fn{i});
end
end

function y = apply_blockdiag_prec_local(x, bd)
% Apply the stored block-diagonal preconditioner.
y = bsxfun(@times, x, bd.dinv);
end

function y = apply_interface_block_prec_local(x, bd, ifb)
% Apply the interface-block preconditioner.
y = apply_blockdiag_prec_local(x, bd);
y(ifb.gamma, :) = ifb.solve \ x(ifb.gamma, :);
end

function [uh, D, rnorms, stats] = call_primme_with_prec_local( ...
Mat, M, n_eigs, target, ops, method, precfun)
% Call PRIMME with a preconditioner.
assert(~isempty(precfun), 'call_primme_with_prec_local requires precfun.');

if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, target, ops, method, precfun);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, target, ops, [], precfun);
end
end

function [lambda, uh] = postprocess_eigs_local( ...
uh, D, M, n_eigenvalues, n_dofs_nurbs, k_pw)
% Sort and normalize eigenpairs.
lambda = diag(D);
[lambda, perm] = sort(real(lambda), 'ascend');
uh = uh(:, perm);

for i = 1:n_eigenvalues
    ni = real(uh(:, i)' * M * uh(:, i));
    uh(:, i) = uh(:, i) / sqrt(abs(ni));
end

uh_pw = uh(n_dofs_nurbs+1:end, :);

for i = 1:n_eigenvalues
    [~, idx] = max(abs(uh_pw(:, i)));
    k0 = k_pw(idx, :);
    idx_neg = find(k_pw(:, 1) == -k0(1) & k_pw(:, 2) == -k0(2), 1);

    c_k = uh_pw(idx, i);
    if ~isempty(idx_neg) && idx_neg ~= idx
        c_mk = uh_pw(idx_neg, i);
        theta = 0.5 * (angle(c_k) + angle(c_mk));
        uh(:, i) = exp(-1i * theta) * uh(:, i);
    else
        theta = angle(c_k);
        uh(:, i) = exp(-1i * theta) * uh(:, i);
        uh_pw2 = uh(n_dofs_nurbs+1:end, i);
        [~, idx2] = max(abs(uh_pw2));
        if real(uh_pw2(idx2)) < 0
            uh(:, i) = -uh(:, i);
        end
    end
end

lambda = lambda(:).';
end

function s = string_or_empty_local(x)
% Convert a value to text or return an empty string.
if isempty(x)
    s = '';
else
    s = char(string(x));
end
end

function s = target_to_string_local(x)
% Convert an eigenvalue target to text.
if isnumeric(x)
    s = char(string(x));
elseif isstring(x) || ischar(x)
    s = char(string(x));
else
    s = '<custom>';
end
end
