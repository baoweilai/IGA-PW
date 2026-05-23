function [lambda, n_dofs_total, meta] = solve_iga_pw_dg(nElem, t, Nc, ~, opts)
%Solve the Example 3 IGA-PW-DG problem.

format long;
t_total = tic;

beta = opts.beta;

L = 4;
a = 0.2;
inner_domains_coordinates = [-a, a, -a, a];

nu = 2; nv = 2;
ConPts = zeros(nu, nv, 2);
x = [-a, a];
y = [-a, a];
ConPts(:,:,1) = [x(1) x(1); x(2) x(2)];
ConPts(:,:,2) = [y(1) y(2); y(1) y(2)];

weights = [1 1; 1 1];
pu0 = 1; pv0 = 1;

nurbs_original = struct();
nurbs_original.ConPts  = ConPts;
nurbs_original.weights = weights;
nurbs_original.pu      = pu0;
nurbs_original.pv      = pv0;
nurbs_original.knotU   = [0 0 1 1];
nurbs_original.knotV   = [0 0 1 1];

pu = pu0 + t;
pv = pv0 + t;

nurbs_refine = IGA_2D_Grid_nElem([], [], pu, pv, nElem);
h = 2 * a / nElem;
n_dofs_nurbs = nurbs_refine.n_dofs_domains;

[k_pw, n_pw_basis] = build_pw_disk(Nc);
pw_dofs_indices = n_dofs_nurbs + (1:n_pw_basis);
n_dofs_total = n_dofs_nurbs + n_pw_basis;

N_Vr = 1;
[k_Vr, n_pw_Vr] = build_pw_disk(N_Vr);

H = sparse(n_dofs_total, n_dofs_total);
M = sparse(n_dofs_total, n_dofs_total);

n_gp = opts.n_gp;

t_nurbs = tic;
[H_nurbs, M_nurbs] = get_nurbs_matrices_cached( ...
    nurbs_original, nurbs_refine, nElem, t, k_Vr, n_pw_Vr, L, n_gp, opts);
time_nurbs = toc(t_nurbs);
fprintf('[IGA] get_nurbs_matrices_cached time = %.4f s\n', time_nurbs);

H(1:n_dofs_nurbs, 1:n_dofs_nurbs) = H_nurbs;
M(1:n_dofs_nurbs, 1:n_dofs_nurbs) = M_nurbs;

t_pw = tic;
[H_pw, M_pw] = get_pw_matrices_cached(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);
time_pw = toc(t_pw);
fprintf('[PW ] get_pw_matrices_cached time = %.4f s\n', time_pw);

H(pw_dofs_indices, pw_dofs_indices) = H_pw;
M(pw_dofs_indices, pw_dofs_indices) = M_pw;

p = k_pw;

[P_Bottom, S_Bottom] = IGA_DG_Bottom_Edge_Assemble(nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
[P_Top,    S_Top   ] = IGA_DG_Top_Edge_Assemble   (nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
[P_Left,   S_Left  ] = IGA_DG_Left_Edge_Assemble  (nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);
[P_Right,  S_Right ] = IGA_DG_Right_Edge_Assemble (nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs_total);

P = P_Bottom + P_Top + P_Left + P_Right;
S = S_Bottom + S_Top + S_Left + S_Right;

sigma = beta * (1 / h + Nc);

Mat0 = H - 0.5 * S - 0.5 * S' + sigma * P;
Mat0 = 0.5 * (Mat0 + Mat0');
M    = 0.5 * (M + M');

ops = struct();
ops.tol         = opts.primme_tol;
ops.maxit       = opts.primme_maxit;
ops.reportLevel = opts.primme_reportLevel;
primme_method   = opts.primme_method;
targetShift     = opts.block_targetShift;

[lambda1_scf, u1_scf, scfInfo] = solve_scf_groundstate( ...
    Mat0, M, P, n_dofs_nurbs, pw_dofs_indices, ...
    nurbs_original, nurbs_refine, k_pw, L, a, Nc, n_gp, ...
    targetShift, ops, primme_method, opts);

lambda = lambda1_scf;

meta = struct();
meta.prec_type      = 'InterfaceBlock';
meta.nElem          = nElem;
meta.pu             = pu;
meta.pv             = pv;
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

meta.lambda1_scf    = lambda1_scf;
meta.scf_iterations = scfInfo.n_iters;
meta.scf_abslambda  = scfInfo.abslambda;
meta.scf_converged  = scfInfo.converged;
meta.time_build_prec = scfInfo.time_build_prec_total;
meta.time_eigs      = scfInfo.time_eigs_total;
meta.time_total     = toc(t_total);

fprintf('[SCF-GS] lambda1 = %.12f, iters = %d, abslambda = %.3e, conv = %d\n', ...
    lambda1_scf, meta.scf_iterations, meta.scf_abslambda, meta.scf_converged);

if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end

    run = struct();
    run.lambda        = lambda(:).';
    run.n_dofs_total  = n_dofs_total;
    run.n_dofs_nurbs  = n_dofs_nurbs;
    run.n_pw_basis    = n_pw_basis;
    run.meta          = meta;

    if opts.save_matrices
        run.M = M;
        if opts.save_mat
            run.Mat0    = Mat0;
            run.Mat     = scfInfo.Mat_final;
            run.H_nurbs = H_nurbs;
            run.M_nurbs = M_nurbs;
        end
    end

    if opts.save_pw_index
        run.k_pw = k_pw;
    end

    if opts.save_nurbs
        run.nurbs_original = nurbs_original;
        run.nurbs_refine   = nurbs_refine;
    end

    if opts.save_eigenvectors
        run.uh     = u1_scf;
        run.u1_scf = u1_scf;
        save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
    else
        save(fullfile(opts.outDir, 'run.mat'), 'run');
    end
end
end

function [k_list, n_basis] = build_pw_disk(Nc)
%Build the plane-wave disk basis.
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

function [H_nurbs, M_nurbs] = get_nurbs_matrices_cached( ...
    nurbs_original, nurbs_refine, nElem, t, k_Vr, n_pw_Vr, L, n_gp, opts)

% Load or assemble cached NURBS matrices.

cacheFile = fullfile(opts.cacheNurbsRoot, ...
    sprintf('NURBS_%s_nElem_%02d_t_%d_ngp_%d_L_%g_NVr_%d.mat', ...
    opts.Example, nElem, t, n_gp, L, n_pw_Vr));

if opts.use_nurbs_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'H_nurbs', 'M_nurbs');
    H_nurbs = S.H_nurbs;
    M_nurbs = S.M_nurbs;
    return;
end

[H_nurbs, M_nurbs] = generate_A_M_NURBS_2D( ...
    nurbs_original, nurbs_refine, k_Vr, n_pw_Vr, L, n_gp, opts.Example);

if opts.use_nurbs_cache
    if ~exist(opts.cacheNurbsRoot, 'dir'), mkdir(opts.cacheNurbsRoot); end
    save(cacheFile, 'H_nurbs', 'M_nurbs', '-v7.3');
end
end

function [H_pw, M_pw] = get_pw_matrices_cached(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts)
%Load or assemble cached plane-wave matrices.
cacheFile = fullfile(opts.cacheRoot, ...
    sprintf('PW_%s_L_%g_Nc_%d_NVr_%d_fft_%d_cheb_%d.mat', ...
    opts.Example, L, Nc, n_pw_Vr, opts.pw_fft_grid_n, opts.inner_cheb_n));

if opts.use_pw_cache && exist(cacheFile, 'file')
    S = load(cacheFile, 'H_pw', 'M_pw');
    H_pw = S.H_pw;
    M_pw = S.M_pw;
    return;
end

[H_pw, M_pw] = generate_A_M_PW_2D(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);

if opts.use_pw_cache
    if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
    save(cacheFile, 'H_pw', 'M_pw', '-v7.3');
end
end
