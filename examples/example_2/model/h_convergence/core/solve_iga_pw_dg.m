function [lambda, n_dofs_total, meta] = solve_iga_pw_dg(Refinement, t, Nc, n_eigenvalues, opts)
% Solve the Example 2 IGA-PW-DG problem.

format long;
t_total = tic;

% ---------------- Parameters ----------------
beta = opts.beta;

% Example 2 geometry
L = 4;
a = 0.2;
% Two inner squares:
% [-1.2,-0.8]x[-0.2,0.2] and [0.8,1.2]x[-0.2,0.2]
inner_domains_coordinates = [ -1.2, -0.8, -0.2, 0.2;
    0.8,  1.2, -0.2, 0.2 ];

% ---------------- Build NURBS patches (two patches) ----------------
pu0 = 1; pv0 = 1;
[knotUe, knotVe] = IGADegreeElevSurface([0 0 1 1], [0 0 1 1], t);
pu = pu0 + t;
pv = pv0 + t;

% Patch 1
nurbs_original_1 = make_rect_patch(inner_domains_coordinates(1,:));
nurbs_original_1.pu    = pu0;
nurbs_original_1.pv    = pv0;
nurbs_original_1.knotU = [0 0 1 1];
nurbs_original_1.knotV = [0 0 1 1];
nurbs_refine_1 = IGA_2D_Grid(knotUe, knotVe, pu, pv, Refinement);
n_dofs_1 = nurbs_refine_1.n_dofs_domains;

% Patch 2
nurbs_original_2 = make_rect_patch(inner_domains_coordinates(2,:));
nurbs_original_2.pu    = pu0;
nurbs_original_2.pv    = pv0;
nurbs_original_2.knotU = [0 0 1 1];
nurbs_original_2.knotV = [0 0 1 1];
nurbs_refine_2 = IGA_2D_Grid(knotUe, knotVe, pu, pv, Refinement);
n_dofs_2 = nurbs_refine_2.n_dofs_domains;

n_dofs_nurbs = n_dofs_1 + n_dofs_2;

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

% ---- NURBS blocks (two patches, block diagonal) ----
[H1, M1] = generate_A_M_NURBS_2D(nurbs_original_1, nurbs_refine_1, k_Vr, n_pw_Vr, L, n_gp, opts.Example);
[H2, M2] = generate_A_M_NURBS_2D(nurbs_original_2, nurbs_refine_2, k_Vr, n_pw_Vr, L, n_gp, opts.Example);

idx1 = 1:n_dofs_1;
idx2 = n_dofs_1 + (1:n_dofs_2);

H(idx1, idx1) = H(idx1, idx1) + H1;
M(idx1, idx1) = M(idx1, idx1) + M1;

H(idx2, idx2) = H(idx2, idx2) + H2;
M(idx2, idx2) = M(idx2, idx2) + M2;

% ---- PW block ----
t_pw = tic;
[H_pw, M_pw] = get_pw_matrices_cached_example2(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts);
time_pw = toc(t_pw);
fprintf('[PW ] get_pw_matrices_cached time = %.6f s\n', time_pw);

H(pw_dofs_indices, pw_dofs_indices) = H_pw;
M(pw_dofs_indices, pw_dofs_indices) = M_pw;

% ---------------- Interface assembly (two patches) ----------------
S = sparse(n_dofs_total, n_dofs_total);
P = sparse(n_dofs_total, n_dofs_total);

p = k_pw;

% Patch 1: offset = 0
off1 = 0;
[P_B1, S_B1] = IGA_DG_Bottom_Edge_Assemble_offset(nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1);
[P_T1, S_T1] = IGA_DG_Top_Edge_Assemble_offset   (nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1);
[P_L1, S_L1] = IGA_DG_Left_Edge_Assemble_offset  (nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1);
[P_R1, S_R1] = IGA_DG_Right_Edge_Assemble_offset (nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1);

% Patch 2: offset = n_dofs_1
off2 = n_dofs_1;
[P_B2, S_B2] = IGA_DG_Bottom_Edge_Assemble_offset(nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2);
[P_T2, S_T2] = IGA_DG_Top_Edge_Assemble_offset   (nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2);
[P_L2, S_L2] = IGA_DG_Left_Edge_Assemble_offset  (nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2);
[P_R2, S_R2] = IGA_DG_Right_Edge_Assemble_offset (nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2);

P = P + P_B1 + P_T1 + P_L1 + P_R1 + P_B2 + P_T2 + P_L2 + P_R2;
S = S + S_B1 + S_T1 + S_L1 + S_R1 + S_B2 + S_T2 + S_L2 + S_R2;

% ---------------- DG penalty coefficient ----------------
h = 0.4 / (2^Refinement);
sigma = beta * (2*a/h + Nc);

Mat = H - 0.5*S - 0.5*S' + sigma * P;

% Symmetrize
Mat = 0.5 * (Mat + Mat');
M   = 0.5 * (M   + M');

% ---------------- meta info ----------------
meta = struct();
meta.prec_type      = 'InterfaceBlock';
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
meta.n_dofs_1       = n_dofs_1;
meta.n_dofs_2       = n_dofs_2;
meta.n_pw_basis     = n_pw_basis;
meta.time_pw        = time_pw;

% ---------------- Solve generalized eigenproblem by PRIMME ----------------
ops = struct();
ops.tol         = opts.primme_tol;
ops.maxit       = opts.primme_maxit;
ops.reportLevel = opts.primme_reportLevel;

t_build_prec = tic;
Pfun = build_interface_block_preconditioner(Mat, P, opts.eps_diag, opts.iface_reg);
meta.time_build_prec = toc(t_build_prec);

t_solve = tic;
[uh, D, rnorms, stats] = call_primme_with_prec(Mat, M, n_eigenvalues, 'SA', ops, opts.primme_method, Pfun);
meta.time_eigs = toc(t_solve);

if isstruct(stats)
    if isfield(stats,'numOuterIterations'), meta.primme_numOuterIterations = stats.numOuterIterations; end
    if isfield(stats,'estimateLargestSVal'),meta.primme_estLargestSVal = stats.estimateLargestSVal; end
    if isfield(stats,'elapsedTime'),        meta.primme_elapsed = stats.elapsedTime; end
end
meta.primme_rnorms = rnorms;

fprintf('[PRIMME] Nc=%d refine=%d target=SA method=%s\n', ...
    Nc, Refinement, string_or_empty(opts.primme_method));
fprintf('        interfaceblock : time = %.6f s\n', meta.time_eigs);

lambda = diag(D);
[lambda, perm] = sort(real(lambda), 'ascend');
uh = uh(:, perm);

fprintf('lambda(interfaceblock) = [');
for k = 1:n_eigenvalues
    if k < n_eigenvalues
        fprintf('%.12f ', lambda(k));
    else
        fprintf('%.12f', lambda(k));
    end
end
fprintf(']\n');

% ---------------- Normalize eigenvectors ----------------
for i = 1:n_eigenvalues
    ni = real(uh(:,i)'*M*uh(:,i));
    uh(:,i) = uh(:,i)./sqrt(abs(ni));
end

% ---------------- Phase alignment (robust) ----------------
uh_pw = uh(n_dofs_nurbs+1:end, :);

for i = 1:n_eigenvalues
    [~, idx] = max(abs(uh_pw(:,i)));
    k0 = k_pw(idx,:);
    idx_neg = find(k_pw(:,1)==-k0(1) & k_pw(:,2)==-k0(2), 1);

    c_k = uh_pw(idx, i);

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
    run.n_dofs_1      = n_dofs_1;
    run.n_dofs_2      = n_dofs_2;
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
        run.nurbs_original_1 = nurbs_original_1;
        run.nurbs_refine_1   = nurbs_refine_1;
        run.nurbs_original_2 = nurbs_original_2;
        run.nurbs_refine_2   = nurbs_refine_2;
    end

    if isfield(opts,'save_eigenvectors') && opts.save_eigenvectors
        run.uh = uh;
        save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
    else
        save(fullfile(opts.outDir, 'run.mat'), 'run');
    end
end

end


function nurbs_original = make_rect_patch(rect)
% Build an affine rectangular NURBS patch.
x1 = rect(1); x2 = rect(2);
y1 = rect(3); y2 = rect(4);

nu = 2; nv = 2;
ConPts = zeros(nu, nv, 2);

ConPts(:,:,1) = [x1 x1; x2 x2];
ConPts(:,:,2) = [y1 y2; y1 y2];

nurbs_original = struct();
nurbs_original.ConPts  = ConPts;
nurbs_original.weights = [1 1; 1 1];
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

function [H_pw, M_pw] = get_pw_matrices_cached_example2(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts)
% Load or assemble the Example 2 plane-wave matrices.
cache_ok = isfield(opts,'use_pw_cache') && opts.use_pw_cache ...
    && isfield(opts,'cacheRoot') && ~isempty(opts.cacheRoot);

if cache_ok
    cacheFile = fullfile(opts.cacheRoot, sprintf( ...
        'PW_EX2_L_%g_Nc_%d_NVr_%d_fft_%d_cheb_%d.mat', ...
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

function [uh, D, rnorms, stats] = call_primme_with_prec(Mat, M, n_eigs, targetMode, ops, method, Pfun)
% Call PRIMME with the TB-DG preconditioner.
assert(~isempty(Pfun), 'call_primme_with_prec requires Pfun.');

if strcmpi(targetMode, 'SA')
    if ~isempty(method)
        [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, method, Pfun);
    else
        [uh, D, rnorms, stats] = primme_eigs(Mat, M, n_eigs, 'SA', ops, [], Pfun);
    end
else
    error('Only targetMode = SA is used in this version.');
end
end

function Pfun = build_interface_block_preconditioner(Ktau, P, eps_diag, iface_reg)
% Build the TB-DG interface-block preconditioner.
d = abs(diag(Ktau));
d(d < eps_diag) = 1;
gamma = find(sum(abs(P), 2) ~= 0);
Ag = 0.5 * (Ktau(gamma, gamma) + Ktau(gamma, gamma)');
delta = iface_reg * max(1, norm(Ag, 1));
Dg = decomposition(Ag + delta * speye(size(Ag)), 'chol');
dinv = 1 ./ d;
Pfun = @(X) apply_interface_block_preconditioner(X, dinv, gamma, Dg);
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

% Assemble the four translated IGA-PW boundary interfaces.

function [P,S] = IGA_DG_Bottom_Edge_Assemble_offset(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, offset)
% Assemble the bottom-edge jump and average operators.

% Read the geometry and refined-basis data.

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;

knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

S = sparse(n_dofs,n_dofs);
P = S;

UBreaks          = nurbs_refine.UBreaks;
uNoEs            = length(UBreaks) - 1;
bottom_edge_dofs = nurbs_refine.bottom_edge_dofs;

% Prepare bottom-edge segments, quadrature, and plane-wave traces.
bottom_edge_node = zeros(uNoEs,2);
for i=1:uNoEs
    bottom_edge_node(i,:) = [UBreaks(i),UBreaks(i+1)];
end

[gp,gw] = grule(10*pu+5);
n_gp = length(gp);

v_bottom_inner = 0;

edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index,1);
basis_grad_minus = zeros(n_pw_basis,2);
basis_minus      = zeros(n_pw_basis,1);
Omega_area       = L*L;

for e=1:uNoEs
    edge_dofs_plus  = bottom_edge_dofs(e,:) + offset;
    edge_dofs       = [edge_dofs_plus,edge_dofs_minus];
    n_edge_dofs     = length(edge_dofs);

    ue = bottom_edge_node(e,:);
    a  = ue(1); b = ue(2);
    J1 = (b-a)/2;

    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        u  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v_bottom_inner);

        tau    = DF_plus(:,1);
        ds     = norm(tau);
        normal = [tau(2);-tau(1)]/ds;

        Jacobi = J1*gw(i)*ds;

        Uders_plus  = bspbasisDers(knotU,pu,u,1);
        Nu_plus     = Uders_plus(1,:)';          DNu_plus = Uders_plus(2,:)';
        Vders_plus  = bspbasisDers(knotV,pv,v_bottom_inner,1);
        Nv_plus     = Vders_plus(1,1:2);         DNv_plus = Vders_plus(2,1:2);

        basis_plus      = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus      = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus      = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;

        for k = 1:n_pw_basis
            basis_minus(k)        = exp(-1i * 2*pi/L *pw_index(k,:)*F)/sqrt(Omega_area);
            basis_grad_minus(k,:) = (1i*2*pi/L)*pw_index(k,:)*basis_minus(k)';
        end

        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*(edge_jump')*Jacobi;

        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2;
        edge_average    = edge_average.';
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;
    end

    P(edge_dofs,edge_dofs) = P(edge_dofs,edge_dofs) + edge_jump_Ae;
    S(edge_dofs,edge_dofs) = S(edge_dofs,edge_dofs) + edge_average_Ae;
end
end

function [P,S] = IGA_DG_Top_Edge_Assemble_offset(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, offset)
% Assemble the top-edge jump and average operators.

% Read the geometry and refined-basis data.

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;

knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

S = sparse(n_dofs,n_dofs);
P = S;

UBreaks        = nurbs_refine.UBreaks;
uNoEs          = length(UBreaks) - 1;
top_edge_dofs  = nurbs_refine.top_edge_dofs;

% Prepare top-edge segments, quadrature, and plane-wave traces.
top_edge_node = zeros(uNoEs,2);
for i=1:uNoEs
    top_edge_node(i,:) = [UBreaks(i),UBreaks(i+1)];
end

[gp,gw] = grule(10*pu+5);
n_gp = length(gp);

v_top_inner = 1;

edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index,1);
basis_grad_minus = zeros(n_pw_basis,2);
basis_minus      = zeros(n_pw_basis,1);
Omega_area       = L*L;

for e=1:uNoEs
    edge_dofs_plus  = top_edge_dofs(e,:) + offset;
    edge_dofs       = [edge_dofs_plus,edge_dofs_minus];
    n_edge_dofs     = length(edge_dofs);

    ue = top_edge_node(e,:);
    a  = ue(1); b = ue(2);
    J1 = (b-a)/2;

    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        u  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u,knotV_o,pv_o,v_top_inner);

        tau    = DF_plus(:,1);
        ds     = norm(tau);
        normal = [-tau(2);tau(1)]/ds;

        Jacobi = J1*gw(i)*ds;

        Uders_plus  = bspbasisDers(knotU,pu,u,1);
        Nu_plus     = Uders_plus(1,:)';              DNu_plus = Uders_plus(2,:)';
        Vders_plus  = bspbasisDers(knotV,pv,v_top_inner,1);
        Nv_plus     = Vders_plus(1,end-1:end);       DNv_plus = Vders_plus(2,end-1:end);

        basis_plus      = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus      = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus      = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;

        for k = 1:n_pw_basis
            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_area);
            basis_grad_minus(k,:) = (1i*2*pi/L)*pw_index(k,:)*basis_minus(k)';
        end

        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;

        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2;
        edge_average    = edge_average.';
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;
    end

    P(edge_dofs,edge_dofs) = P(edge_dofs,edge_dofs) + edge_jump_Ae;
    S(edge_dofs,edge_dofs) = S(edge_dofs,edge_dofs) + edge_average_Ae;
end
end

function [P,S] = IGA_DG_Left_Edge_Assemble_offset(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, offset)
% Assemble the left-edge jump and average operators.

% Read the geometry and refined-basis data.

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;

knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

S = sparse(n_dofs,n_dofs);
P = S;

VBreaks         = nurbs_refine.VBreaks;
vNoEs           = length(VBreaks) - 1;
left_edge_dofs  = nurbs_refine.left_edge_dofs;

% Prepare left-edge segments, quadrature, and plane-wave traces.
left_edge_node = zeros(vNoEs,2);
for i=1:vNoEs
    left_edge_node(i,:) = [VBreaks(i),VBreaks(i+1)];
end

[gp,gw] = grule(10*pv+5);
n_gp = length(gp);

u_left_inner = 0;

edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index,1);
basis_grad_minus = zeros(n_pw_basis,2);
basis_minus      = zeros(n_pw_basis,1);
Omega_area       = L*L;

for e = 1:vNoEs
    edge_dofs_plus  = left_edge_dofs(e,:) + offset;
    edge_dofs       = [edge_dofs_plus,edge_dofs_minus];
    n_edge_dofs     = length(edge_dofs);

    ve = left_edge_node(e,:);
    a  = ve(1); b = ve(2);
    J1 = (b-a)/2;

    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        v  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u_left_inner,knotV_o,pv_o,v);

        tau    = DF_plus(:,2);
        ds     = norm(tau);
        normal = [-tau(2);tau(1)]/ds;

        Jacobi = J1*gw(i)*ds;

        Uders_plus  = bspbasisDers(knotU,pu,u_left_inner,1);
        Nu_plus     = Uders_plus(1,1:2)';  DNu_plus = Uders_plus(2,1:2)';
        Vders_plus  = bspbasisDers(knotV,pv,v,1);
        Nv_plus     = Vders_plus(1,:);     DNv_plus = Vders_plus(2,:);

        basis_plus      = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus      = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus      = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;

        for k = 1:n_pw_basis
            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_area);
            basis_grad_minus(k,:) = (1i*2*pi/L)*pw_index(k,:)*basis_minus(k)';
        end

        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;

        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2;
        edge_average    = edge_average.';
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;
    end

    P(edge_dofs,edge_dofs) = P(edge_dofs,edge_dofs) + edge_jump_Ae;
    S(edge_dofs,edge_dofs) = S(edge_dofs,edge_dofs) + edge_average_Ae;
end
end

function [P,S] = IGA_DG_Right_Edge_Assemble_offset(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, offset)
% Assemble the right-edge jump and average operators.

% Read the geometry and refined-basis data.

ConPts_o   = nurbs_original.ConPts;
weights_o  = nurbs_original.weights;
knotU_o    = nurbs_original.knotU;
knotV_o    = nurbs_original.knotV;
pu_o       = nurbs_original.pu;
pv_o       = nurbs_original.pv;

knotU      = nurbs_refine.Ubar;
knotV      = nurbs_refine.Vbar;
pu         = nurbs_refine.pu;
pv         = nurbs_refine.pv;

S = sparse(n_dofs,n_dofs);
P = S;

VBreaks          = nurbs_refine.VBreaks;
vNoEs            = length(VBreaks) - 1;
right_edge_dofs  = nurbs_refine.right_edge_dofs;

% Prepare right-edge segments, quadrature, and plane-wave traces.
right_edge_node = zeros(vNoEs,2);
for i=1:vNoEs
    right_edge_node(i,:) = [VBreaks(i),VBreaks(i+1)];
end

[gp,gw] = grule(10*pv+5);
n_gp = length(gp);

u_right_inner = 1;

edge_dofs_minus  = plane_wave_dofs_index;
n_pw_basis       = size(pw_index,1);
basis_grad_minus = zeros(n_pw_basis,2);
basis_minus      = zeros(n_pw_basis,1);
Omega_area       = L*L;

for e = 1:vNoEs
    edge_dofs_plus  = right_edge_dofs(e,:) + offset;
    edge_dofs       = [edge_dofs_plus,edge_dofs_minus];
    n_edge_dofs     = length(edge_dofs);

    ve = right_edge_node(e,:);
    a  = ve(1); b = ve(2);
    J1 = (b-a)/2;

    edge_jump_Ae    = zeros(n_edge_dofs,n_edge_dofs);
    edge_average_Ae = zeros(n_edge_dofs,n_edge_dofs);

    for i=1:n_gp
        v  = ((b-a)*gp(i) +a+b)/2;
        [F,DF_plus] = NurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u_right_inner,knotV_o,pv_o,v);

        tau    = DF_plus(:,2);
        ds     = norm(tau);
        normal = [tau(2); -tau(1)]/ds;

        Jacobi = J1*gw(i)*ds;

        Uders_plus  = bspbasisDers(knotU,pu,u_right_inner,1);
        Nu_plus     = Uders_plus(1,end-1:end)';   DNu_plus = Uders_plus(2,end-1:end)';
        Vders_plus  = bspbasisDers(knotV,pv,v,1);
        Nv_plus     = Vders_plus(1,:);            DNv_plus = Vders_plus(2,:);

        basis_plus      = Nu_plus*Nv_plus;    basis_plus = basis_plus(:);
        DNu_v_plus      = DNu_plus*Nv_plus;   DNu_v_plus = DNu_v_plus(:);
        DNv_u_plus      = Nu_plus*DNv_plus;   DNv_u_plus = DNv_u_plus(:);
        basis_grad_plus = [DNu_v_plus,DNv_u_plus]/DF_plus;

        for k = 1:n_pw_basis
            basis_minus(k)        = exp(-1i*2*pi/L*pw_index(k,:)*F)/sqrt(Omega_area);
            basis_grad_minus(k,:) = (1i*2*pi/L)*pw_index(k,:)*basis_minus(k)';
        end

        edge_jump    = [basis_plus; - basis_minus];
        edge_jump_Ae = edge_jump_Ae + edge_jump*edge_jump'*Jacobi;

        edge_average    = [basis_grad_plus; basis_grad_minus]*normal/2;
        edge_average    = edge_average.';
        edge_average_Ae = edge_average_Ae + edge_jump*edge_average*Jacobi;
    end

    P(edge_dofs,edge_dofs) = P(edge_dofs,edge_dofs) + edge_jump_Ae;
    S(edge_dofs,edge_dofs) = S(edge_dofs,edge_dofs) + edge_average_Ae;
end
end
