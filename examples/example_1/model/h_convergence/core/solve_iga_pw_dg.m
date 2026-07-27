function [lambda, n_dofs_total, meta] = solve_iga_pw_dg(Refinement, t, Nc, n_eigenvalues, opts)
% Solve the Example 1 IGA-PW-DG problem.

format long;
t_total = tic;

% ---------------- Parameters ----------------
beta = opts.beta;

% Example 1 geometry
L  = 4;       % Omega = [-L/2, L/2]^2
a  = 0.2;     % inner domain = [-a,a]^2
inner_domains_coordinates = [-a, a, -a, a];  % [x1 x2 y1 y2]

% ---------------- Build NURBS patch (single patch) ----------------
nu = 2; nv = 2;

ConPts = zeros(nu, nv, 2);
x = [-a, a];
y = [-a, a];

ConPts(:,:,1) = [x(1) x(1); x(2) x(2)];
ConPts(:,:,2) = [y(1) y(2); y(1) y(2)];

weights = [1 1; 1 1];
knotU   = [0 0 1 1];
knotV   = [0 0 1 1];

pu0 = 1; pv0 = 1;
nurbs_original = struct();
nurbs_original.ConPts  = ConPts;
nurbs_original.weights = weights;
nurbs_original.pu      = pu0;
nurbs_original.pv      = pv0;
nurbs_original.knotU   = knotU;
nurbs_original.knotV   = knotV;

% degree elevation: final degrees
[knotUe, knotVe] = IGADegreeElevSurface(knotU, knotV, t);
pu = pu0 + t;
pv = pv0 + t;

% refined grid
nurbs_refine = IGA_2D_Grid(knotUe, knotVe, pu, pv, Refinement);

% physical h
h = 2*a/(2^Refinement);

% total dofs in inner domain
n_dofs_nurbs = nurbs_refine.n_dofs_domains;

% ---------------- Plane-wave basis (|k| <= Nc) ----------------
[k_pw, n_pw_basis] = build_pw_disk(Nc);
pw_dofs_indices = n_dofs_nurbs + (1:n_pw_basis);
n_dofs_total = n_dofs_nurbs + n_pw_basis;

% For periodized Vr truncation
N_Vr = 2;
[k_Vr, n_pw_Vr] = build_pw_disk(N_Vr);

% ---------------- Assemble matrices ----------------
H = sparse(n_dofs_total, n_dofs_total);
M = sparse(n_dofs_total, n_dofs_total);

n_gp = opts.n_gp;

% Assemble NURBS block
t_nurbs = tic;
[H_nurbs, M_nurbs] = generate_A_M_NURBS_2D(nurbs_refine, k_Vr, L, n_gp, a);
time_nurbs = toc(t_nurbs);
fprintf('[IGA] volume assembly time = %.4f s\n', time_nurbs);
H(1:n_dofs_nurbs, 1:n_dofs_nurbs) = H_nurbs;
M(1:n_dofs_nurbs, 1:n_dofs_nurbs) = M_nurbs;

% Assemble PW block (cache by Nc)
t_pw = tic;
[H_pw, M_pw] = get_pw_matrices_cached(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);
time_pw = toc(t_pw);
fprintf('[PW ] get_pw_matrices_cached time = %.4f s\n', time_pw);

H(pw_dofs_indices, pw_dofs_indices) = H_pw;
M(pw_dofs_indices, pw_dofs_indices) = M_pw;

% Interface assembly (single patch)
t_dg = tic;
[P, S] = assemble_DG_square_interface_fast( ...
    nurbs_refine, k_pw, pw_dofs_indices, L, a, n_dofs_total);
time_dg = toc(t_dg);
fprintf('[DG ] interface assembly time = %.4f s\n', time_dg);

% DG penalty coefficient
sigma = beta * (1/h + Nc);

Mat = H - 0.5*S - 0.5*S' + sigma * P;

% Symmetrize
Mat = 0.5 * (Mat + Mat');
M   = 0.5 * (M   + M');

% ---------------- meta info ----------------
meta = struct();
meta.prec_type      = 'InterfaceBlock';
meta.pattern_tol    = opts.pattern_tol;
meta.pu             = pu;
meta.pv             = pv;
meta.Refinement     = Refinement;
meta.Nc             = Nc;
meta.L              = L;
meta.a              = a;
meta.h              = h;
meta.beta           = beta;
meta.sigma          = sigma;
meta.n_dofs_nurbs   = n_dofs_nurbs;
meta.n_pw_basis     = n_pw_basis;
meta.time_nurbs     = time_nurbs;
meta.time_pw        = time_pw;
meta.time_dg        = time_dg;

% ---------------- PRIMME options ----------------
ops = struct();
ops.tol         = opts.primme_tol;
ops.maxit       = opts.primme_maxit;
ops.reportLevel = opts.primme_reportLevel;

primme_method = [];
if isfield(opts,'primme_method')
    primme_method = opts.primme_method;
end

% ---------------- Solve generalized eigenproblem by PRIMME ----------------
t_build_prec = tic;
[Pfun, precInfo] = build_interface_block_preconditioner( ...
    Mat, P, pw_dofs_indices, opts.pattern_tol, opts.eps_diag, opts.iface_reg);
meta.time_build_prec = toc(t_build_prec);
meta.n_gamma = precInfo.n_gamma;
meta.n_gamma_iga = precInfo.n_gamma_iga;
meta.n_eta = precInfo.n_eta;
meta.penalty_pattern_dropped = precInfo.penalty_pattern_dropped;
meta.trace_pattern_dropped = precInfo.trace_pattern_dropped;
fprintf(['[PREC] patternTol=%.1e nGamma=%d nGammaIGA=%d ' ...
    'dropped(P)=%d dropped(trace)=%d\n'], ...
    opts.pattern_tol, precInfo.n_gamma, precInfo.n_gamma_iga, ...
    precInfo.penalty_pattern_dropped, precInfo.trace_pattern_dropped);

t_solve = tic;
[uh, D, rnorms, stats] = call_primme_with_prec(Mat, M, n_eigenvalues, ops, primme_method, Pfun);
meta.time_eigs = toc(t_solve);

if isstruct(stats)
    if isfield(stats,'numOuterIterations'), meta.primme_numOuterIterations = stats.numOuterIterations; end
    if isfield(stats,'estimateLargestSVal'),meta.primme_estLargestSVal = stats.estimateLargestSVal; end
    if isfield(stats,'elapsedTime'),        meta.primme_elapsed = stats.elapsedTime; end
end
meta.primme_rnorms = rnorms;

fprintf('[PRIMME] Nc=%d refine=%d target=SA method=%s\n', ...
    Nc, Refinement, string_or_empty(primme_method));
fprintf('        interfaceblock : time = %.6f s\n', meta.time_eigs);

lambda = diag(D);
[lambda, perm] = sort(real(lambda), 'ascend');
uh = uh(:, perm);

% ---------------- Normalize eigenvectors in M-inner product ----------------
for i = 1:n_eigenvalues
    ni = real(uh(:,i)' * M * uh(:,i));
    uh(:,i) = uh(:,i) / sqrt(abs(ni));
end

% ---------------- Fix global phase (robust) ----------------
uh_pw = uh(n_dofs_nurbs+1:end, :);

for i = 1:n_eigenvalues
    [~, idx] = max(abs(uh_pw(:,i)));
    k0 = k_pw(idx,:);
    idx_neg = find(k_pw(:,1)==-k0(1) & k_pw(:,2)==-k0(2), 1);

    c_k  = uh_pw(idx, i);

    if ~isempty(idx_neg) && idx_neg ~= idx
        c_mk = uh_pw(idx_neg, i);
        theta = 0.5 * (angle(c_k) + angle(c_mk));
        uh(:,i) = exp(-1i*theta) * uh(:,i);
    else
        theta = angle(c_k);
        uh(:,i) = exp(-1i*theta) * uh(:,i);

        uh_pw2 = uh(n_dofs_nurbs+1:end, i);
        [~, idx2] = max(abs(uh_pw2));
        if real(uh_pw2(idx2)) < 0
            uh(:,i) = -uh(:,i);
        end
    end
end

meta.time_total = toc(t_total);

% ---------------- Save run data ----------------
if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end

    run = struct();
    run.lambda        = lambda(:).';
    run.n_dofs_total  = n_dofs_total;
    run.n_dofs_nurbs  = n_dofs_nurbs;
    run.n_pw_basis    = n_pw_basis;
    run.meta          = meta;

    if isfield(opts,'save_matrices') && opts.save_matrices
        run.M = M;
        if isfield(opts,'save_mat') && opts.save_mat
            run.Mat = Mat;
        end
    end

    if isfield(opts,'save_pw_index') && opts.save_pw_index
        run.k_pw = k_pw;
    end

    if isfield(opts,'save_nurbs') && opts.save_nurbs
        run.nurbs_original = nurbs_original;
        run.nurbs_refine   = nurbs_refine;
    end

    if isfield(opts,'save_eigenvectors') && opts.save_eigenvectors
        run.uh = uh;
        save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
    else
        save(fullfile(opts.outDir, 'run.mat'), 'run');
    end
end

end


function [k_list, n_basis] = build_pw_disk(Nc)
% Build the plane-wave disk basis.
N = floor(Nc);
k_list = zeros((2*N+1)^2, 2);
n_basis = 0;
for k1 = -N:N
    m = floor(sqrt(N^2 - k1^2));
    for k2 = -m:m
        n_basis = n_basis + 1;
        k_list(n_basis,:) = [k1, k2];
    end
end
k_list = k_list(1:n_basis,:);
end

function [H_pw, M_pw] = get_pw_matrices_cached(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts)
% Load or assemble cached plane-wave matrices.
cache_ok = isfield(opts,'use_pw_cache') && opts.use_pw_cache ...
    && isfield(opts,'cacheRoot') && ~isempty(opts.cacheRoot);

if cache_ok
    cacheFile = fullfile(opts.cacheRoot, sprintf( ...
        'PW_L_%g_Nc_%d_NVr_%d_fft_%d_cheb_%d.mat', ...
        L, Nc, n_pw_Vr, opts.pw_fft_grid_n, opts.inner_cheb_n));
    if exist(cacheFile, 'file')
        S = load(cacheFile, 'H_pw', 'M_pw');
        H_pw = S.H_pw; M_pw = S.M_pw;
        return;
    end
end

[H_pw, M_pw, ~] = generate_A_M_PW_2D(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);

if cache_ok
    if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
    save(cacheFile, 'H_pw', 'M_pw', '-v7.3');
end
end

function [uh, D, rnorms, stats] = call_primme_with_prec(Mat, M, n_eigs, ops, method, Pfun)
% Call PRIMME with the TB-DG preconditioner.
assert(~isempty(Pfun), 'call_primme_with_prec requires Pfun.');
if ~isempty(method)
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, method, Pfun);
else
    [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, [], Pfun);
end
end

function [Pfun, info] = build_interface_block_preconditioner( ...
    Ktau, P, pw_dofs_indices, pattern_tol, eps_diag, iface_reg)
% Build the TB-DG interface-block preconditioner.
d = abs(diag(Ktau));
d(d < eps_diag) = 1;

[ii, jj, vv] = find(P);
keep = abs(vv) >= pattern_tol;
Ppat = sparse(ii(keep), jj(keep), vv(keep), size(P, 1), size(P, 2));
[ii, jj] = find(Ppat);
gamma = unique([pw_dofs_indices(:); ii; jj]);

Ag = 0.5 * (Ktau(gamma, gamma) + Ktau(gamma, gamma)');
[ia, ja, va] = find(Ag);
keepTrace = abs(va) >= pattern_tol;
Ag = sparse(ia(keepTrace), ja(keepTrace), va(keepTrace), ...
    numel(gamma), numel(gamma));
delta = iface_reg * max(1, norm(Ag, 1));
Dg = decomposition(Ag + delta * speye(size(Ag)), 'chol');
dinv = 1 ./ d;
Pfun = @(X) apply_interface_block_preconditioner(X, dinv, gamma, Dg);

info.n_gamma = numel(gamma);
info.n_gamma_iga = numel(setdiff(gamma, pw_dofs_indices(:)));
info.n_eta = size(P, 1) - info.n_gamma;
info.penalty_pattern_dropped = nnz(~keep);
info.trace_pattern_dropped = nnz(~keepTrace);
end

function Z = apply_interface_block_preconditioner(X, dinv, gamma, Dg)
% Apply the TB-DG interface-block preconditioner.
Z = bsxfun(@times, X, dinv);
Z(gamma, :) = Dg \ X(gamma, :);
end

function s = string_or_empty(x)
% Convert a text value to a string.
if isempty(x)
    s = '';
else
    s = char(string(x));
end
end
