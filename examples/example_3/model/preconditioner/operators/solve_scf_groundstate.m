function [lambda1_final, u1_final, info] = solve_scf_groundstate( ...
Mat0, M, n_dofs_nurbs, pw_dofs_indices, ...
    nurbs_original, nurbs_refine, k_pw, L, a, Nc, n_gp, ...
    targetShift, ops, primme_method, idxI, idxA, opts)

scf_maxit      = opts.scf_maxit;
scf_tol_lambda = opts.scf_tol_lambda;
scf_mixing     = opts.scf_mixing;
track_n_eigs   = opts.scf_track_n_eigs;

assert(scf_maxit >= 1, 'SCF maximum iteration count must be positive.');
assert(track_n_eigs >= 1, 'SCF tracked eigenvalue count must be positive.');
assert(isfield(opts, 'scf_pw_grid_m') && isscalar(opts.scf_pw_grid_m) && opts.scf_pw_grid_m > 0, ...
    'opts.scf_pw_grid_m must be provided.');
assert(isfinite(scf_mixing) && scf_mixing >= 0 && scf_mixing <= 1, ...
    'SCF mixing parameter must be in [0, 1].');
m_pw = opts.scf_pw_grid_m;

info = struct();
info.n_iters                  = 0;
info.converged                = false;
info.track_n_eigs             = track_n_eigs;
info.last_branch_id           = 1;
info.time_eigs_total          = 0;
info.time_build_prec_total    = 0;

lambda_hist = zeros(scf_maxit, 1);

n_track_init = min(track_n_eigs, size(Mat0,1) - 1);
if n_track_init < 1
    n_track_init = 1;
end

[U0, D0, ~, ~, solveInfo0, ~] = solve_generalized_eigs_with_blockdiag_jacobi( ...
    Mat0, M, n_track_init, targetShift, ops, primme_method, idxI, idxA, n_dofs_nurbs, k_pw);

u_old = U0(:,1);
u_old = normalize_in_M(u_old, M);
lambda_old = real(D0(1,1));

info.time_eigs_total       = info.time_eigs_total + solveInfo0.time_eigs;
info.time_build_prec_total = info.time_build_prec_total + solveInfo0.time_build_prec;

for it = 1:scf_maxit
    [Nmat, ~] = build_nonlinear_operator( ...
        u_old, n_dofs_nurbs, pw_dofs_indices, ...
        nurbs_original, nurbs_refine, n_gp, ...
        k_pw, L, a, m_pw);

    Mat_cur = Mat0 + Nmat;
    Mat_cur = 0.5 * (Mat_cur + Mat_cur');

    n_track = min(track_n_eigs, size(Mat_cur,1) - 1);
    if n_track < 1
        n_track = 1;
    end

    [Ucand, Dcand, ~, ~, solveInfo1, ~] = solve_generalized_eigs_with_blockdiag_jacobi( ...
        Mat_cur, M, n_track, targetShift, ops, primme_method, idxI, idxA, n_dofs_nurbs, k_pw);

    [u_raw, lambda_new, branch_id] = select_branch_by_overlap(Ucand, Dcand, u_old, M);

    u_raw = align_phase_to_previous(u_raw, u_old, M);
    u_raw = normalize_in_M(u_raw, M);

    if scf_mixing >= 1
        u_mix = u_raw;
    elseif scf_mixing <= 0
        u_mix = u_old;
    else
        u_mix = scf_mixing * u_raw + (1 - scf_mixing) * u_old;
    end

    u_mix = align_phase_to_previous(u_mix, u_old, M);
    u_mix = normalize_in_M(u_mix, M);

    abslambda = abs(lambda_new - lambda_old);

    info.time_eigs_total       = info.time_eigs_total + solveInfo1.time_eigs;
    info.time_build_prec_total = info.time_build_prec_total + solveInfo1.time_build_prec;

    lambda_hist(it) = lambda_new;

    info.n_iters        = it;
    info.abslambda      = abslambda;
    info.last_branch_id = branch_id;

    if it >= 3
        if abs(lambda_hist(it) - lambda_hist(it-2)) < 5e-5 && ...
                abs(lambda_hist(it) - lambda_hist(it-1)) > 2e-4
            scf_mixing = max(0.02, 0.5 * scf_mixing);
        end
    end

    u_old = u_mix;
    lambda_old = lambda_new;

    if abslambda < scf_tol_lambda
        info.converged = true;
        break;
    end
end

[Nmat_final, ~] = build_nonlinear_operator( ...
    u_old, n_dofs_nurbs, pw_dofs_indices, ...
    nurbs_original, nurbs_refine, n_gp, ...
    k_pw, L, a, m_pw);

Mat_final = Mat0 + Nmat_final;
Mat_final = 0.5 * (Mat_final + Mat_final');

n_track_final = min(track_n_eigs, size(Mat_final,1) - 1);
if n_track_final < 1
    n_track_final = 1;
end

[Uf, Df, ~, ~, solveInfof, ~] = solve_generalized_eigs_with_blockdiag_jacobi( ...
    Mat_final, M, n_track_final, targetShift, ops, primme_method, idxI, idxA, n_dofs_nurbs, k_pw);

[u1_final, lambda1_final, branch_id_final] = select_branch_by_overlap(Uf, Df, u_old, M);
u1_final = align_phase_to_previous(u1_final, u_old, M);
u1_final = normalize_in_M(u1_final, M);

info.time_eigs_total       = info.time_eigs_total + solveInfof.time_eigs;
info.time_build_prec_total = info.time_build_prec_total + solveInfof.time_build_prec;

info.Mat_final = Mat_final;

info.last_branch_id = branch_id_final;

end

function [u_sel, lambda_sel, idx_sel] = select_branch_by_overlap(Ucand, Dcand, u_prev, M)
%Select the eigenvector branch by overlap.
lams = real(diag(Dcand));
nCand = size(Ucand, 2);
overlaps = zeros(nCand, 1);
for j = 1:nCand
    overlaps(j) = abs(u_prev' * M * Ucand(:,j));
end
[~, idx_sel] = max(overlaps);
u_sel = Ucand(:, idx_sel);
lambda_sel = lams(idx_sel);
end

function [uh, D, rnorms, stats, solveInfo, Pfun] = solve_generalized_eigs_with_blockdiag_jacobi( ...
Mat, M, n_eigs, targetShift, ops, primme_method, idxI, idxA, n_dofs_nurbs, k_pw)

solveInfo = struct();

Ktau = Mat - targetShift * M;

t_build = tic;
[Pfun, dI, dA] = build_block_diag_jacobi_prec(Ktau, idxI, idxA);
solveInfo.time_build_prec = toc(t_build);

t_solve = tic;
[uh, D, rnorms, stats] = call_primme_local(Mat, M, n_eigs, targetShift, ops, primme_method, Pfun);
solveInfo.time_eigs = toc(t_solve);

lam = diag(D);
[lam, perm] = sort(real(lam), 'ascend');
uh = uh(:, perm);
D  = diag(lam);

for j = 1:size(uh,2)
    uh(:,j) = normalize_in_M(uh(:,j), M);
    uh(:,j) = fix_global_phase_single_pw(uh(:,j), n_dofs_nurbs, k_pw);
end

Pfun = @(X) apply_block_diag_jacobi_prec(X, dI, dA);
end

function [Pfun, dI, dA] = build_block_diag_jacobi_prec(K, idxI, idxA)
%Build the block-diagonal Jacobi preconditioner.
dI = full(diag(K(idxI, idxI)));
dA = full(diag(K(idxA, idxA)));

regI = 1e-12 * max(1, norm(dI, inf));
regA = 1e-12 * max(1, norm(dA, inf));

dI(abs(dI) < regI) = regI;
dA(abs(dA) < regA) = regA;

Pfun = @(X) apply_block_diag_jacobi_prec(X, dI, dA);
end

function Z = apply_block_diag_jacobi_prec(X, dI, dA)
%Apply the block-diagonal Jacobi preconditioner.
nI = numel(dI);
X1 = X(1:nI, :);
X2 = X(nI+1:end, :);

Z1 = bsxfun(@rdivide, X1, dI);
Z2 = bsxfun(@rdivide, X2, dA);

Z = [Z1; Z2];
end

function [uh, D, rnorms, stats] = call_primme_local(Mat, M, n_eigs, targetShift, ops, primme_method, Pfun)
%Call PRIMME for the local eigenproblem.
[uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, targetShift, ops, primme_method, Pfun);
end

function u = normalize_in_M(u, M)
%Normalize an eigenvector in the mass inner product.
nu = real(u' * M * u);
u = u / sqrt(abs(nu));
end

function u = fix_global_phase_single_pw(u, n_dofs_nurbs, k_pw)
%Fix the global phase using one plane-wave coefficient.
u_pw = u(n_dofs_nurbs+1:end);
[~, idx] = max(abs(u_pw));
k0 = k_pw(idx, :);
idx_neg = find(k_pw(:,1) == -k0(1) & k_pw(:,2) == -k0(2), 1);

c_k = u_pw(idx);

if ~isempty(idx_neg) && idx_neg ~= idx
    c_mk = u_pw(idx_neg);
    theta = 0.5 * (angle(c_k) + angle(c_mk));
    u = exp(-1i * theta) * u;
else
    theta = angle(c_k);
    u = exp(-1i * theta) * u;
    u_pw2 = u(n_dofs_nurbs+1:end);
    [~, idx2] = max(abs(u_pw2));
    if real(u_pw2(idx2)) < 0
        u = -u;
    end
end
end

function u_new = align_phase_to_previous(u_new, u_old, M)
%Align the eigenvector phase to the previous iterate.
alpha = u_old' * M * u_new;
if abs(alpha) > 0
    u_new = exp(-1i * angle(alpha)) * u_new;
end
end
