function [lambda, n_dofs_total, meta] = solve_iga_pw_dg(Refinement, t, Nc, n_eigenvalues, opts)
%Solve the Example 2 IGA-PW-DG problem.

format long;
t_total = tic;
solve_mode = "interfaceblock";

% ---------------- Parameters ----------------
beta = opts.beta;
DIM  = 2; %#ok<NASGU>

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
fprintf('[IGA] assembling / loading patch 1 NURBS block ...\n');
t_nurbs1 = tic;
[H1, M1] = get_nurbs_block_cached_example2( ...
    'patch1', nurbs_original_1, nurbs_refine_1, k_Vr, n_pw_Vr, L, n_gp, opts);
fprintf('[IGA] patch 1 done in %.6f s\n', toc(t_nurbs1));

fprintf('[IGA] assembling / loading patch 2 NURBS block ...\n');
t_nurbs2 = tic;
[H2, M2] = get_nurbs_block_cached_example2( ...
    'patch2', nurbs_original_2, nurbs_refine_2, k_Vr, n_pw_Vr, L, n_gp, opts);
fprintf('[IGA] patch 2 done in %.6f s\n', toc(t_nurbs2));

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
fprintf('[DG ] assembling interface terms ...\n');
t_iface = tic;
S = sparse(n_dofs_total, n_dofs_total);
P = sparse(n_dofs_total, n_dofs_total);

p = k_pw;

% Patch 1: offset = 0
off1 = 0;
[P_B1, S_B1] = get_interface_edge_cached_example2('patch1_bottom', @IGA_DG_Bottom_Edge_Assemble_offset, ...
    nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1, Nc, opts);
[P_T1, S_T1] = get_interface_edge_cached_example2('patch1_top', @IGA_DG_Top_Edge_Assemble_offset, ...
    nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1, Nc, opts);
[P_L1, S_L1] = get_interface_edge_cached_example2('patch1_left', @IGA_DG_Left_Edge_Assemble_offset, ...
    nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1, Nc, opts);
[P_R1, S_R1] = get_interface_edge_cached_example2('patch1_right', @IGA_DG_Right_Edge_Assemble_offset, ...
    nurbs_original_1, nurbs_refine_1, p, pw_dofs_indices, L, n_dofs_total, off1, Nc, opts);

% Patch 2: offset = n_dofs_1
off2 = n_dofs_1;
[P_B2, S_B2] = get_interface_edge_cached_example2('patch2_bottom', @IGA_DG_Bottom_Edge_Assemble_offset, ...
    nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2, Nc, opts);
[P_T2, S_T2] = get_interface_edge_cached_example2('patch2_top', @IGA_DG_Top_Edge_Assemble_offset, ...
    nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2, Nc, opts);
[P_L2, S_L2] = get_interface_edge_cached_example2('patch2_left', @IGA_DG_Left_Edge_Assemble_offset, ...
    nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2, Nc, opts);
[P_R2, S_R2] = get_interface_edge_cached_example2('patch2_right', @IGA_DG_Right_Edge_Assemble_offset, ...
    nurbs_original_2, nurbs_refine_2, p, pw_dofs_indices, L, n_dofs_total, off2, Nc, opts);

P = P + P_B1 + P_T1 + P_L1 + P_R1 + P_B2 + P_T2 + P_L2 + P_R2;
S = S + S_B1 + S_T1 + S_L1 + S_R1 + S_B2 + S_T2 + S_L2 + S_R2;
fprintf('[DG ] interface assembly done in %.6f s\n', toc(t_iface));

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
meta.solve_mode     = char(solve_mode);
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
Pfun = build_interface_block_preconditioner_local(Mat, P, opts);
meta.time_build_prec = toc(t_build_prec);
fprintf('[PREC] %s build done in %.6f s\n', meta.prec_type, meta.time_build_prec);

t_solve = tic;
[uh, D, rnorms, ~] = call_primme_with_prec(Mat, M, n_eigenvalues, 'SA', ops, opts.primme_method, Pfun);
meta.time_eigs = toc(t_solve);

meta.primme_rnorms = rnorms;

fprintf('[PRIMME] Nc=%d refine=%d target=SA method=%s\n', ...
    Nc, Refinement, string_or_empty(opts.primme_method));
fprintf('        %s : build = %.6f s, solve = %.6f s\n', ...
    meta.prec_type, meta.time_build_prec, meta.time_eigs);

lambda = diag(D);
[lambda, perm] = sort(real(lambda), 'ascend');
uh = uh(:, perm);

fprintf('lambda(%s) = [', meta.prec_type);
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
%Build rect patch.
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

function [H_pw, M_pw] = get_pw_matrices_cached_example2(L, Nc, inner_domains_coordinates, k_Vr, n_pw_Vr, opts)
%Return PW matrices cached example2.
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

function [H_blk, M_blk] = get_nurbs_block_cached_example2(tag, nurbs_original, nurbs_refine, k_Vr, n_pw_Vr, L, n_gp, opts)
%Return NURBS block cached example2.
cache_ok = isfield(opts, 'cacheNurbsRoot') && ~isempty(opts.cacheNurbsRoot);

if cache_ok
    cacheFile = fullfile(opts.cacheNurbsRoot, sprintf( ...
        'NURBS_EX2_%s_nu_%d_nv_%d_pu_%d_pv_%d_nd_%d_NVr_%d_ngp_%d.mat', ...
        tag, numel(nurbs_refine.Ubar), numel(nurbs_refine.Vbar), ...
        nurbs_refine.pu, nurbs_refine.pv, nurbs_refine.n_dofs_domains, n_pw_Vr, n_gp));
    if exist(cacheFile, 'file')
        S = load(cacheFile, 'H_blk', 'M_blk');
        H_blk = S.H_blk;
        M_blk = S.M_blk;
        return;
    end
end

[H_blk, M_blk] = generate_A_M_NURBS_2D(nurbs_original, nurbs_refine, k_Vr, n_pw_Vr, L, n_gp, opts.Example);

if cache_ok
    if ~exist(opts.cacheNurbsRoot, 'dir'), mkdir(opts.cacheNurbsRoot); end
    save(cacheFile, 'H_blk', 'M_blk', '-v7.3');
end
end

function [P_edge, S_edge] = get_interface_edge_cached_example2(tag, edgefun, ...
nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs, offset, Nc, opts)
%Load or assemble one cached interface edge.

cache_ok = isfield(opts, 'cacheInterfaceRoot') && ~isempty(opts.cacheInterfaceRoot);

if cache_ok
    cacheFile = fullfile(opts.cacheInterfaceRoot, sprintf( ...
        'IFACE_EX2_%s_nu_%d_nv_%d_pu_%d_pv_%d_nd_%d_Nc_%d.mat', ...
        tag, numel(nurbs_refine.Ubar), numel(nurbs_refine.Vbar), ...
        nurbs_refine.pu, nurbs_refine.pv, nurbs_refine.n_dofs_domains, Nc));
    if exist(cacheFile, 'file')
        fprintf('[DG ] %s loaded from cache.\n', tag);
        S = load(cacheFile, 'P_edge', 'S_edge');
        P_edge = S.P_edge;
        S_edge = S.S_edge;
        return;
    end
else
    cacheFile = '';
end

fprintf('[DG ] %s computing ...\n', tag);
t_edge = tic;
[P_edge, S_edge] = edgefun(nurbs_original, nurbs_refine, p, pw_dofs_indices, L, n_dofs, offset);
fprintf('[DG ] %s done in %.6f s\n', tag, toc(t_edge));

if cache_ok
    if ~exist(opts.cacheInterfaceRoot, 'dir'), mkdir(opts.cacheInterfaceRoot); end
    save(cacheFile, 'P_edge', 'S_edge', '-v7.3');
end
end

function [uh, D, rnorms, stats] = call_primme_with_prec(Mat, M, n_eigs, targetMode, ops, method, Pfun)
%Call PRIMME with prec.
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

function Pfun = build_interface_block_preconditioner_local(Ktau, P, opts)
%Build the TB-DG interface-block preconditioner.
eps_diag = get_opt_local(opts, 'eps_diag', 1e-12);
iface_reg = get_opt_local(opts, 'iface_reg', 1e-12);
d = abs(diag(Ktau));
d(d < eps_diag) = 1;
gamma = find(sum(abs(P), 2) ~= 0);
Ag = 0.5 * (Ktau(gamma, gamma) + Ktau(gamma, gamma)');
delta = iface_reg * max(1, norm(Ag, 1));
Dg = decomposition(Ag + delta * speye(size(Ag)), 'chol');
dinv = 1 ./ d;
Pfun = @(X) apply_interface_block_prec_local(X, dinv, gamma, Dg);
end

function Z = apply_interface_block_prec_local(X, dinv, gamma, Dg)
%Apply the interface-block preconditioner.
Z = bsxfun(@times, X, dinv);
Z(gamma, :) = Dg \ X(gamma, :);
end

function s = string_or_empty(x)
%Convert a text value to a string.
if isempty(x)
    s = '';
else
    s = char(string(x));
end
end

function value = get_opt_local(opts, name, default)
%Return one option value.
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = default;
end
end

% Offset versions of your 4 edge assembly functions
% Only difference from original: edge_dofs_plus = local_edge_dofs + offset

function [P,S] = IGA_DG_Bottom_Edge_Assemble_offset(nurbs_original, nurbs_refine, pw_index, plane_wave_dofs_index, L, n_dofs, offset)
%Assemble matrices or interface terms for the method.

DIM = 2; %#ok<NASGU>

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
pw_grad_scale    = 1i * 2*pi / L;
pw_val_scale     = 1 / sqrt(Omega_area);

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

        phase_arg         = pw_index * F;
        basis_minus       = exp(-1i * 2*pi/L * phase_arg) * pw_val_scale;
        basis_grad_minus  = pw_grad_scale * bsxfun(@times, pw_index, basis_minus);

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
%Assemble matrices or interface terms for the method.

DIM = 2; %#ok<NASGU>

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
pw_grad_scale    = 1i * 2*pi / L;
pw_val_scale     = 1 / sqrt(Omega_area);

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

        phase_arg         = pw_index * F;
        basis_minus       = exp(-1i * 2*pi/L * phase_arg) * pw_val_scale;
        basis_grad_minus  = pw_grad_scale * bsxfun(@times, pw_index, basis_minus);

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
%Assemble matrices or interface terms for the method.

DIM = 2; %#ok<NASGU>

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
pw_grad_scale    = 1i * 2*pi / L;
pw_val_scale     = 1 / sqrt(Omega_area);

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

        phase_arg         = pw_index * F;
        basis_minus       = exp(-1i * 2*pi/L * phase_arg) * pw_val_scale;
        basis_grad_minus  = pw_grad_scale * bsxfun(@times, pw_index, basis_minus);

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
%Assemble matrices or interface terms for the method.

DIM = 2; %#ok<NASGU>

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
pw_grad_scale    = 1i * 2*pi / L;
pw_val_scale     = 1 / sqrt(Omega_area);

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

        phase_arg         = pw_index * F;
        basis_minus       = exp(-1i * 2*pi/L * phase_arg) * pw_val_scale;
        basis_grad_minus  = pw_grad_scale * bsxfun(@times, pw_index, basis_minus);

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
