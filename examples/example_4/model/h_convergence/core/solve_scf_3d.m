function [lambda, n_dofs_total, meta] = solve_scf_3d(Refinement, t, Nc, opts)
%Solve the 3-D SCF IGA-PW-DG problem.
arguments
    Refinement
    t
    Nc
    opts struct
end

format long;
t_total = tic;

assert(exist('primme_eigs', 'file') == 2, ...
    'primme_eigs not found. Please compile primme_mex and add PRIMME/Matlab.');

pdeg = 1 + t;
[refine_mode, refine_value, refine_tag] = parse_refinement_input(Refinement);

beta = opts.beta;
L = opts.L;
a = opts.a;
nuclearCharge = opts.nuclear_charge;
n_gp = opts.n_gp;
inner_cheb_n = opts.inner_cheb_n;
targetShift = opts.block_targetShift;
scf_maxit = opts.scf_maxit;
scf_tol_eig = opts.scf_tol_eig;
scf_tol_rho = opts.scf_tol_rho;
scf_beta = opts.scf_beta;
scf_stopping_rule = opts.scf_stopping_rule;
if isstring(scf_stopping_rule)
    scf_stopping_rule = char(scf_stopping_rule);
end
scf_stopping_rule = lower(scf_stopping_rule);
referenceLambda = opts.reference_lambda;
referenceEnergy = opts.reference_energy;
useExchangeCorrelation = logical(opts.use_exchange_correlation);
eps_iface = opts.eps_iface;
iface_explicit_gamma_max = opts.iface_explicit_gamma_max;
use_direct_interface_gamma = logical(opts.use_direct_interface_gamma);
use_tensor_api = logical(opts.use_tensor_api);
use_tensor_face_data = logical(opts.use_tensor_face_data);
use_iga_grid_eval_cache = logical(opts.use_iga_grid_eval_cache);
iface_direct_trace_entry_max = opts.iface_direct_trace_entry_max;
preconditionerType = opts.preconditioner_type;
if isstring(preconditionerType)
    preconditionerType = char(preconditionerType);
end

primme_tol = opts.primme_tol;
primme_maxit = opts.primme_maxit;
primme_method = opts.primme_method;
primme_reportLevel = opts.primme_reportLevel;

inner_box = [-a, a, -a, a, -a, a];

% 1) Build inner IGA patch
nu = 2; nv = 2; nw = 2;
ConPts = zeros(nu, nv, nw, 3);

xs = [-a, a];
ys = [-a, a];
zs = [-a, a];

for i = 1:nu
    for j = 1:nv
        for k = 1:nw
            ConPts(i, j, k, 1) = xs(i);
            ConPts(i, j, k, 2) = ys(j);
            ConPts(i, j, k, 3) = zs(k);
        end
    end
end

weights = ones(nu, nv, nw);
knotU = [0 0 1 1];
knotV = [0 0 1 1];
knotW = [0 0 1 1];

nurbs_original = struct();
nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.pu = 1;
nurbs_original.pv = 1;
nurbs_original.pw = 1;
nurbs_original.knotU = knotU;
nurbs_original.knotV = knotV;
nurbs_original.knotW = knotW;

[knotUe, knotVe, knotWe] = IGADegreeElevVolume(knotU, knotV, knotW, t);
pu = 1 + t; pv = 1 + t; pw = 1 + t;

nurbs_refine = IGA_3D_Grid(knotUe, knotVe, knotWe, pu, pv, pw, Refinement);

hx = 2 * a / nurbs_refine.uNoEs;
hy = 2 * a / nurbs_refine.vNoEs;
hz = 2 * a / nurbs_refine.wNoEs;
hmin = min([hx, hy, hz]);
hmax = max([hx, hy, hz]);

n_dofs_nurbs = nurbs_refine.n_dofs_domains;

% 2) PW basis and static operators
[k_pw, n_pw_basis] = build_pw_ball(Nc);
n_dofs_total = n_dofs_nurbs + n_pw_basis;

N_Vr = 1;
[k_Vr, n_pw_Vr] = build_pw_ball(N_Vr);

t_static_build = tic;
t_kin_mass = tic;
[K_nurbs_kin, M_nurbs, nurbs_kin_mass_meta] = generate_K_M_NURBS_3D( ...
    nurbs_original, nurbs_refine, n_gp, opts);
time_build_nurbs_kin_mass = toc(t_kin_mass);
K_nurbs_kin = 0.5 * (K_nurbs_kin + K_nurbs_kin');
M_nurbs = 0.5 * (M_nurbs + M_nurbs');

t_pw_static = tic;
[pwDataStatic, ~, pwCacheFile, pwCacheHit] = get_pw_operator_cached_3D_local( ...
    L, Nc, inner_box, k_Vr, n_pw_Vr, opts);
time_build_pw_static = toc(t_pw_static);

t_face = tic;
if use_tensor_face_data
    faceData = build_face_data_3D(nurbs_refine, a, L, floor(Nc));
else
    faceData = build_face_data_3D(nurbs_refine, a, L, floor(Nc), 'legacy');
end
if ~use_tensor_api
    for iface = 1:numel(faceData)
        faceData{iface}.use_tensor_api = false;
    end
end
time_build_face_data = toc(t_face);

P_ii = sparse(n_dofs_nurbs, n_dofs_nurbs);
S_ii = sparse(n_dofs_nurbs, n_dofs_nurbs);

for iface = 1:numel(faceData)
    F = faceData{iface};
    W = spdiags(F.w, 0, numel(F.w), numel(F.w));
    P_ii = P_ii + F.TI' * W * F.TI;
    S_ii = S_ii + F.TI' * W * F.GI;
end

sigma = beta * (1 / hmin + Nc);

% 3) External potential grid and Hartree-grid cache
hartree_grid_n = get_hartree_grid_n(opts);
gridCache = build_midpoint_grid_cache_3D(L, hartree_grid_n, a);
t_iga_grid_eval = tic;
if use_iga_grid_eval_cache
    igaGridEval = build_iga_grid_eval_matrix_3D( ...
        nurbs_refine, gridCache.x_inner, gridCache.y_inner, gridCache.z_inner, a);
else
    igaGridEval = [];
end
time_build_iga_grid_eval = toc(t_iga_grid_eval);
t_vext = tic;
[VextGrid, vextCacheFile, vextCacheHit] = get_vext_grid_cached_3D( ...
    gridCache, k_Vr, n_pw_Vr, opts);
time_build_vext_grid = toc(t_vext);
t_vext_iga = tic;
[Vext_ii, nurbs_vext_meta] = assemble_NURBS_potential_interp_3D( ...
    nurbs_original, nurbs_refine, VextGrid, L, n_gp, opts);
time_build_vext_iga = toc(t_vext_iga);
time_build_static_operator = toc(t_static_build);

K_ii_static = K_nurbs_kin + Vext_ii - 0.25 * S_ii - 0.25 * S_ii' + sigma * P_ii;
K_ii_static = 0.5 * (K_ii_static + K_ii_static');
Bfun_static = @(x) apply_global_B(x, M_nurbs, n_dofs_nurbs, pwDataStatic);
AstaticFun = @(x) apply_global_A(x, K_ii_static, n_dofs_nurbs, sigma, pwDataStatic, faceData);

% 4) Initial guess
uh_prev = build_initial_guess_local(pwDataStatic, n_dofs_nurbs, n_dofs_total);
uh_prev = normalize_in_B(uh_prev, Bfun_static);
lambda_prev = real(uh_prev' * AstaticFun(uh_prev));
rhoOnlyOpts = struct( ...
    'rhoGrid', [], ...
    'skip_poisson', true, ...
    'use_exchange_correlation', useExchangeCorrelation, ...
    'iga_grid_eval_matrix', igaGridEval);
[rho_old, ~, ~, ~, ~] = build_scf_potentials_3D( ...
    uh_prev, n_dofs_nurbs, nurbs_refine, pwDataStatic, gridCache, rhoOnlyOpts);

lambda_hist = zeros(scf_maxit, 1);
lambda_rel_hist = zeros(scf_maxit, 1);
energy_hist = zeros(scf_maxit, 1);
hartree_energy_hist = zeros(scf_maxit, 1);
rho_res_hist = zeros(scf_maxit, 1);
charge_hist = zeros(scf_maxit, 1);
time_scf_potential_hist = zeros(scf_maxit, 1);
time_hartree_iga_hist = zeros(scf_maxit, 1);
time_pw_update_hist = zeros(scf_maxit, 1);
time_prec_build_hist = zeros(scf_maxit, 1);
time_primme_solve_hist = zeros(scf_maxit, 1);
scf_converged = false;
preconditionerName = 'interface_block';
preconditionerMeta = struct();
rhoGrid = [];
VHGrid = [];
VxcGrid = [];
epsxcGrid = [];
aux = struct();
uGrid = [];
lambda_rel_change = 0;
rho_res = 0;

% 5) SCF loop
for it = 1:scf_maxit
    rho_input = rho_old;
    t_scf_pot = tic;
    [~, VHGridIn, VxcGridIn, ~, ~] = build_scf_potentials_3D( ...
        [], n_dofs_nurbs, nurbs_refine, pwDataStatic, gridCache, ...
        struct('rhoGrid', rho_input, ...
        'skip_poisson', false, ...
        'use_exchange_correlation', useExchangeCorrelation, ...
        'iga_grid_eval_matrix', []));
    time_scf_potential_hist(it) = toc(t_scf_pot);

    VscfGrid = VHGridIn + VxcGridIn;
    t_hartree_iga = tic;
    Vhartree_ii = assemble_NURBS_potential_interp_3D( ...
        nurbs_original, nurbs_refine, VscfGrid, L, n_gp, opts);
    time_hartree_iga_hist(it) = toc(t_hartree_iga);

    K_ii = K_ii_static + Vhartree_ii;
    K_ii = 0.5 * (K_ii + K_ii');

    pwUpdateOpts = struct();
    pwUpdateOpts.inner_cheb_n = inner_cheb_n;
    pwUpdateOpts.combine_with_static = true;
    t_pw_update = tic;
    pwDataIter = update_pw_potential_from_grid_3D( ...
        VscfGrid, L, inner_box, pwDataStatic, pwUpdateOpts);
    time_pw_update_hist(it) = toc(t_pw_update);

    Afun = @(x) apply_global_A(x, K_ii, n_dofs_nurbs, sigma, pwDataIter, faceData);
    Bfun = @(x) apply_global_B(x, M_nurbs, n_dofs_nurbs, pwDataStatic);

    t_prec_build = tic;
    [Pfun, preconditionerName, preconditionerMeta] = build_preconditioner_local( ...
        preconditionerType, K_ii, M_nurbs, targetShift, pwDataIter, sigma, ...
        n_dofs_nurbs, n_dofs_total, faceData, Afun, Bfun, ...
        eps_iface, iface_explicit_gamma_max, use_direct_interface_gamma, ...
        iface_direct_trace_entry_max);
    time_prec_build_hist(it) = toc(t_prec_build);

    ops = struct();
    ops.tol = primme_tol;
    ops.maxit = primme_maxit;
    ops.reportLevel = primme_reportLevel;
    ops.display = 0;
    ops.isreal = false;
    ops.isdouble = true;
    ops.ishermitian = true;
    ops.v0 = uh_prev;

    t_primme = tic;
    [uh_raw, D] = call_primme_operator( ...
        Afun, Bfun, n_dofs_total, 1, targetShift, ops, primme_method, Pfun);
    time_primme_solve_hist(it) = toc(t_primme);

    lambda_raw = real(D(1, 1));
    uh_raw = uh_raw(:, 1);

    uh_raw = normalize_in_B(uh_raw, Bfun_static);
    uh_raw = align_phase(uh_raw, uh_prev, Bfun_static);

    [rho_state, VHGridState, VxcGridState, epsxcGridState, aux_state] = build_scf_potentials_3D( ...
        uh_raw, n_dofs_nurbs, nurbs_refine, pwDataStatic, gridCache, ...
        struct('rhoGrid', [], ...
        'skip_poisson', false, ...
        'use_exchange_correlation', useExchangeCorrelation, ...
        'iga_grid_eval_matrix', igaGridEval));
    rho_res = norm(rho_state(:) - rho_input(:)) / max(norm(rho_state(:)), 1);
    lambda_rel_change = abs(lambda_raw - lambda_prev) / max(abs(lambda_raw), 1);

    if useExchangeCorrelation
        Etot = 2 * lambda_raw - 0.5 * aux_state.int_rho_vh - aux_state.int_rho_vxc + aux_state.Exc;
    else
        Etot = compute_total_hartree_energy(uh_raw, AstaticFun, aux_state.int_rho_vh);
    end

    lambda_hist(it) = lambda_raw;
    lambda_rel_hist(it) = lambda_rel_change;
    energy_hist(it) = Etot;
    hartree_energy_hist(it) = 0.5 * aux_state.int_rho_vh;
    rho_res_hist(it) = rho_res;
    charge_hist(it) = aux_state.total_charge;

    uh_prev = uh_raw;
    lambda_prev = lambda_raw;
    rhoGrid = rho_state;
    VHGrid = VHGridState;
    VxcGrid = VxcGridState;
    epsxcGrid = epsxcGridState;
    aux = aux_state;
    uGrid = aux_state.uGrid;

    switch scf_stopping_rule
        case {'lambda_only', 'lambda'}
            scf_should_stop = it >= 2 && lambda_rel_change < scf_tol_eig;
        case {'lambda_and_rho', 'both'}
            scf_should_stop = it >= 2 && ...
                lambda_rel_change < scf_tol_eig && rho_res < scf_tol_rho;
        case {'rho_only', 'rho'}
            scf_should_stop = it >= 2 && rho_res < scf_tol_rho;
        otherwise
            error('Unsupported scf_stopping_rule: %s', scf_stopping_rule);
    end

    if scf_should_stop
        scf_converged = true;
        break;
    end

    rho_old = (1 - scf_beta) * rho_input + scf_beta * rho_state;
end

% 6) Final postprocess
if isempty(rhoGrid)
    [rhoGrid, VHGrid, VxcGrid, epsxcGrid, aux] = build_scf_potentials_3D( ...
        uh_prev, n_dofs_nurbs, nurbs_refine, pwDataStatic, gridCache, ...
        struct('rhoGrid', [], ...
        'skip_poisson', false, ...
        'use_exchange_correlation', useExchangeCorrelation, ...
        'iga_grid_eval_matrix', igaGridEval));
    uGrid = aux.uGrid;
end

if useExchangeCorrelation
    Etot = 2 * lambda_prev - 0.5 * aux.int_rho_vh - aux.int_rho_vxc + aux.Exc;
else
    Etot = compute_total_hartree_energy(uh_prev, AstaticFun, aux.int_rho_vh);
end
uCenter = evaluate_hybrid_state_3D( ...
    uh_prev, n_dofs_nurbs, nurbs_refine, pwDataStatic, a, 0, 0, 0);
lambda = lambda_prev;
lambdaError = compute_reference_error(lambda, referenceLambda);
energyError = compute_reference_error(Etot, referenceEnergy);
scfIters = it;

meta = struct();
meta.Example = opts.Example;
meta.example_title = 'Example 5 (3D Hartree model of a helium atom)';
meta.run_mode = opts.runMode;
meta.Nc = Nc;
meta.t = t;
meta.pdeg = pdeg;
meta.refine_mode = refine_mode;
meta.refine_value = refine_value;
meta.refine_tag = refine_tag;
meta.L = L;
meta.a = a;
meta.hmin = hmin;
meta.hmax = hmax;
meta.sigma = sigma;
meta.n_dofs_nurbs = n_dofs_nurbs;
meta.n_pw_basis = n_pw_basis;
meta.n_dofs_total = n_dofs_total;
meta.scf_iters = scfIters;
meta.scf_converged = scf_converged;
meta.lambda_hist = lambda_hist(1:scfIters);
meta.lambda_rel_hist = lambda_rel_hist(1:scfIters);
meta.energy_hist = energy_hist(1:scfIters);
meta.hartree_energy_hist = hartree_energy_hist(1:scfIters);
meta.rho_res_hist = rho_res_hist(1:scfIters);
meta.charge_hist = charge_hist(1:scfIters);
meta.scf_beta = scf_beta;
meta.scf_tol_eig = scf_tol_eig;
meta.scf_tol_rho = scf_tol_rho;
meta.scf_stopping_rule = scf_stopping_rule;
meta.primme_tol = primme_tol;
meta.primme_maxit = primme_maxit;
meta.primme_method = primme_method;
meta.primme_reportLevel = primme_reportLevel;
meta.final_lambda_rel_change = lambda_rel_change;
meta.final_density_residual = rho_res;
meta.reference_lambda = referenceLambda;
meta.lambda_error = lambdaError;
meta.energy_total = Etot;
meta.hartree_energy = 0.5 * aux.int_rho_vh;
meta.u_center = uCenter;
meta.reference_energy = referenceEnergy;
meta.energy_error = energyError;
meta.init_guess_mode = opts.init_guess_mode;
meta.rhoGrid = rhoGrid;
meta.VHGrid = VHGrid;
meta.VxcGrid = VxcGrid;
meta.epsxcGrid = epsxcGrid;
meta.uGrid = uGrid;
meta.grid_mFFT = gridCache.mFFT;
meta.hartree_grid_n = gridCache.mFFT;
meta.global_fft_grid_n = gridCache.mFFT;
meta.inner_cheb_n = pwDataStatic.inner_cheb_n;
meta.pw_fft_grid_n = pwDataStatic.pw_fft_grid_n;
meta.use_direct_interface_gamma = use_direct_interface_gamma;
meta.use_tensor_api = use_tensor_api;
meta.use_tensor_face_data = use_tensor_face_data;
meta.use_iga_grid_eval_cache = use_iga_grid_eval_cache;
meta.iface_direct_trace_entry_max = iface_direct_trace_entry_max;
meta.use_exchange_correlation = useExchangeCorrelation;
meta.nuclear_charge = nuclearCharge;
meta.alpha = 5;
meta.hartree_zero_mode = 0;
meta.rho_definition = 'rho = 2|u|^2';
meta.pw_cache_file = pwCacheFile;
meta.pw_cache_hit = pwCacheHit;
meta.vext_cache_file = vextCacheFile;
meta.vext_cache_hit = vextCacheHit;
meta.preconditioner_name = preconditionerName;
meta.preconditioner_type = preconditionerType;
meta.preconditioner_info = preconditionerMeta;
meta.time_build_operator = time_build_static_operator;
meta.time_build_static_operator = time_build_static_operator;
meta.time_build_nurbs_kin_mass = time_build_nurbs_kin_mass;
meta.time_build_pw_static = time_build_pw_static;
meta.time_build_face_data = time_build_face_data;
meta.time_build_iga_grid_eval = time_build_iga_grid_eval;
meta.time_build_vext_grid = time_build_vext_grid;
meta.time_build_vext_iga = time_build_vext_iga;
meta.nurbs_kin_mass_info = nurbs_kin_mass_meta;
meta.nurbs_vext_info = nurbs_vext_meta;
meta.time_scf_potential_hist = time_scf_potential_hist(1:scfIters);
meta.time_hartree_iga_hist = time_hartree_iga_hist(1:scfIters);
meta.time_pw_update_hist = time_pw_update_hist(1:scfIters);
meta.time_prec_build_hist = time_prec_build_hist(1:scfIters);
meta.time_primme_solve_hist = time_primme_solve_hist(1:scfIters);
meta.time_scf_potential_total = sum(meta.time_scf_potential_hist);
meta.time_hartree_iga_total = sum(meta.time_hartree_iga_hist);
meta.time_pw_update_total = sum(meta.time_pw_update_hist);
meta.time_build_prec_total = sum(meta.time_prec_build_hist);
meta.time_eigs = sum(meta.time_primme_solve_hist);
meta.time_total = toc(t_total);
meta.density_diagnostics = struct();
meta.static_blocks = struct( ...
    'mass_cached', true, ...
    'kinetic_cached', true, ...
    'dg_interface_cached', true, ...
    'external_potential_cached', true, ...
    'dynamic_updates', {{'IGA/IGA volume Hartree block', 'PW/PW potential block'}}, ...
    'mixed_block_fixed', true);

% 7) Save
if isfield(opts, 'outDir') && ~isempty(opts.outDir)
    if ~exist(opts.outDir, 'dir'), mkdir(opts.outDir); end

    if opts.save_density_diagnostics
        meta.density_diagnostics = save_density_diagnostics( ...
            opts.outDir, uh_prev, n_dofs_nurbs, nurbs_refine, pwDataStatic, ...
            rhoGrid, VHGrid, gridCache, a, meta.lambda_hist, meta.rho_res_hist);
    end

    run = struct();
    run.lambda = lambda;
    run.n_dofs_total = n_dofs_total;
    run.n_dofs_nurbs = n_dofs_nurbs;
    run.n_pw_basis = n_pw_basis;
    run.meta = meta;
    run.rhoGrid = rhoGrid;
    run.VHGrid = VHGrid;
    run.scf_history = struct( ...
        'lambda', meta.lambda_hist, ...
        'lambda_rel', meta.lambda_rel_hist, ...
        'energy', meta.energy_hist, ...
        'hartree_energy', meta.hartree_energy_hist, ...
        'rho_res', meta.rho_res_hist, ...
        'charge', meta.charge_hist);

    if opts.save_eigenvectors
        run.uh = uh_prev;
    end
    if opts.save_pw_index
        run.k_pw = k_pw;
    end
    if opts.save_nurbs
        run.nurbs_original = nurbs_original;
        run.nurbs_refine = nurbs_refine;
    end

    save(fullfile(opts.outDir, 'run.mat'), 'run', '-v7.3');
end
end

function y = apply_global_B(x, M_nurbs, nI, pwData)
%Apply global b.
xI = x(1:nI, :); xP = x(nI+1:end, :);
y = [M_nurbs * xI; PW3D_apply_mass(xP, pwData)];
end

function y = apply_global_A(x, K_ii, nI, sigma, pwData, faceData)
%Apply global a.
xI = x(1:nI, :); xP = x(nI+1:end, :);
yI = K_ii * xI;
[yP, yI_face] = PW3D_apply_stiff(xP, xI, pwData, faceData, sigma);
y = [yI + yI_face; yP];
end

function [Pfun, preconditionerName, info] = build_preconditioner_local( ...
preconditionerType, K_ii, M_nurbs, targetShift, pwDataIter, sigma, nI, nTotal, faceData, Afun, Bfun, eps_iface, gamma_max, ...
    use_direct_gamma, trace_entry_max)
%Build the requested preconditioner.
switch lower(preconditionerType)
    case 'interface_block'
        [Pfun, preconditionerName, info] = build_interface_block_prec_local( ...
            K_ii, M_nurbs, targetShift, nI, nTotal, faceData, Afun, Bfun, ...
            pwDataIter, sigma, eps_iface, gamma_max, use_direct_gamma, trace_entry_max);
    case 'blockdiag_jacobi'
        [Pfun, preconditionerName, info] = build_blockdiag_prec_local( ...
            K_ii, M_nurbs, targetShift, pwDataIter, sigma, nI);
    otherwise
        error('Unsupported preconditioner_type: %s', preconditionerType);
end
end

function [Pfun, preconditionerName, info] = build_blockdiag_prec_local( ...
    K_ii, M_nurbs, targetShift, pwDataIter, sigma, nI)
%Build the explicit block-diagonal Jacobi preconditioner.

D_ii = abs(diag(K_ii - targetShift * M_nurbs));
D_ii(D_ii < 1e-12) = 1;

D_pw = abs(pwDataIter.stiffDiag + sigma * pwDataIter.facePenaltyDiagApprox ...
    - targetShift * pwDataIter.massDiag);
D_pw(D_pw < 1e-12) = 1;

Pfun = @(x) apply_blockdiag_jacobi(x, D_ii, D_pw, nI);
preconditionerName = 'blockdiag_jacobi';
info = struct();
info.status = 'connected';
info.type = preconditionerName;
info.n_gamma = 0;
info.n_eta = nI;
end

function y = apply_blockdiag_jacobi(x, D_ii, D_pw, nI)
%Apply the block-diagonal Jacobi preconditioner.
xI = x(1:nI, :) ./ D_ii;
xP = x(nI+1:end, :) ./ D_pw;
y = [xI; xP];
end

function [Pfun, preconditionerName, info] = build_interface_block_prec_local( ...
K_ii, M_nurbs, targetShift, nI, nTotal, faceData, Afun, Bfun, ...
    pwData, sigma, eps_iface, gamma_max, use_direct_gamma, trace_entry_max)
%Build the interface block.

[gamma, eta, interfaceIga] = identify_interface_index_sets_local(nI, nTotal, faceData);

info = struct();
info.status = 'prepared';
info.type = 'interface_block';
info.n_gamma = numel(gamma);
info.n_eta = numel(eta);
info.interface_iga_dofs = interfaceIga;
info.gamma_build_method = 'probe';
info.time_build_gamma = 0;
info.trace_entry_max = trace_entry_max;

assert(~isempty(gamma), 'The interface-block preconditioner has no gamma index set.');
assert(numel(gamma) <= gamma_max, ...
    'The interface-block gamma set exceeds iface_explicit_gamma_max.');

diagAtauI = diag(K_ii - targetShift * M_nurbs);
diagEta = diagAtauI(eta);

t_gamma = tic;
if use_direct_gamma
    [Ag, directInfo] = build_interface_gamma_block_direct_local( ...
        K_ii, M_nurbs, targetShift, pwData, sigma, faceData, interfaceIga, trace_entry_max);
    info.gamma_build_method = 'direct_explicit';
    info.direct_gamma_info = directInfo;
else
    Ag = extract_gamma_block_local(Afun, Bfun, gamma, targetShift, nTotal);
end
info.time_build_gamma = toc(t_gamma);

prec = build_interface_block_prec_3D(Ag, gamma, eta, diagEta, nTotal, ...
    struct('eps_iface', eps_iface));

Pfun = prec.apply;
preconditionerName = 'interface_block';
info.status = 'connected';
info.factor_type = prec.factor_type;
info.delta = prec.delta;
info.gamma = gamma;
info.eta = eta;
end

function [gamma, eta, interfaceIga] = identify_interface_index_sets_local(nI, nTotal, faceData)
%Compute interface index sets.
interfaceMask = false(nI, 1);
for iface = 1:numel(faceData)
    cols = any(abs(faceData{iface}.TI) > 0, 1);
    interfaceMask = interfaceMask | cols(:);
end

interfaceIga = find(interfaceMask);
pwIdx = (nI+1:nTotal).';
gamma = unique([interfaceIga; pwIdx]);
eta = setdiff((1:nI).', interfaceIga);
end

function Ag = extract_gamma_block_local(Afun, Bfun, gamma, targetShift, nTotal)
%Extract gamma block.
ng = numel(gamma);
Ag = zeros(ng, ng);
blockSize = min(32, ng);

for j0 = 1:blockSize:ng
    j1 = min(j0 + blockSize - 1, ng);
    cols = j0:j1;
    E = sparse(gamma(cols), 1:numel(cols), 1, nTotal, numel(cols));
    Y = Afun(E) - targetShift * Bfun(E);
    Ag(:, cols) = full(Y(gamma, :));
end

Ag = sparse(0.5 * (Ag + Ag'));
end

function [Ag, info] = build_interface_gamma_block_direct_local( ...
K_ii, M_nurbs, targetShift, pwData, sigma, faceData, interfaceIga, trace_entry_max)
%Build the interface gamma block directly.

info = struct();
info.status = 'prepared';
info.method = 'direct_explicit';
info.n_interface_iga = numel(interfaceIga);
info.n_pw = pwData.n_pw;
info.max_face_trace_entries = 0;
info.total_face_quadrature = 0;
info.n_trace_blocks = 0;

nPw = pwData.n_pw;
nGammaIga = numel(interfaceIga);

for iface = 1:numel(faceData)
    nFaceEntries = numel(faceData{iface}.w) * nPw;
    info.max_face_trace_entries = max(info.max_face_trace_entries, nFaceEntries);
    info.total_face_quadrature = info.total_face_quadrature + numel(faceData{iface}.w);
end

Aii = K_ii(interfaceIga, interfaceIga) - targetShift * M_nurbs(interfaceIga, interfaceIga);
[App, Mpp] = assemble_pw_volume_blocks_direct_local(pwData);
App = App - targetShift * Mpp;
Api = complex(zeros(nPw, nGammaIga));

for iface = 1:numel(faceData)
    F = faceData{iface};

    App = App + assemble_pw_face_block_direct_local(F, pwData, sigma);

    if isfield(F, 'tensor') && logical(F.use_tensor_api)
        [ApiFace, nGroups] = assemble_iga_pw_face_coupling_tensor_local( ...
            F, pwData, sigma, interfaceIga);
        Api = Api + ApiFace;
        info.n_trace_blocks = info.n_trace_blocks + nGroups;
    else
        TIg = F.TI(:, interfaceIga);
        GIg = F.GI(:, interfaceIga);
        W = F.w(:);

        WTIg = spdiags(W, 0, numel(W), numel(W)) * TIg;
        WGIg = spdiags(W, 0, numel(W), numel(W)) * GIg;

        blockSize = max(1, min(nPw, floor(trace_entry_max / max(1, numel(W)))));
        for j0 = 1:blockSize:nPw
            j1 = min(j0 + blockSize - 1, nPw);
            cols = j0:j1;
            [Tm, normalFactor] = build_face_trace_mat_block_local(F, pwData, cols);

            TmWGIg = Tm' * WGIg;
            TmWTIg = Tm' * WTIg;
            Api(cols, :) = Api(cols, :) ...
                + 0.25 * TmWGIg ...
                - 0.25 * bsxfun(@times, conj(normalFactor(:)), TmWTIg) ...
                - sigma * TmWTIg;
            info.n_trace_blocks = info.n_trace_blocks + 1;
        end
    end
end

Aip = Api';
Ag = [Aii, Aip; Api, App];
Ag = sparse(0.5 * (Ag + Ag'));
info.status = 'direct_ok';
end

function [App, Mpp] = assemble_pw_volume_blocks_direct_local(pwData)
%Assemble matrices or interface terms for the method.
k = pwData.k_pw;
nPw = size(k, 1);
qN = floor((size(pwData.Uker, 1) - 1) / 2);

diffx = k(:, 1) - k(:, 1).';
diffy = k(:, 2) - k(:, 2).';
diffz = k(:, 3) - k(:, 3).';

lin = sub2ind(size(pwData.Uker), diffx + qN + 1, diffy + qN + 1, diffz + qN + 1);
Umat = pwData.Uker(lin);
Vmat = pwData.Vker(lin);

k2 = k(:, 1) .^ 2 + k(:, 2) .^ 2 + k(:, 3) .^ 2;
dotK = k(:, 1) * k(:, 1).' + k(:, 2) * k(:, 2).' + k(:, 3) * k(:, 3).';

App = 0.5 * pwData.alpha ^ 2 * (diag(k2) - dotK .* Umat) + Vmat;
Mpp = eye(nPw) - Umat;

App = 0.5 * (App + App');
Mpp = 0.5 * (Mpp + Mpp');
end

function AppFace = assemble_pw_face_block_direct_local(F, pwData, sigma)
%Assemble matrices or interface terms for the method.
alpha = pwData.alpha;
N = pwData.N;
Omega = pwData.Omega;
k = pwData.k_pw;

switch F.type
    case 'x'
        tang1 = k(:, 2) + N + 1;
        tang2 = k(:, 3) + N + 1;
        kn = k(:, 1);
    case 'y'
        tang1 = k(:, 1) + N + 1;
        tang2 = k(:, 3) + N + 1;
        kn = k(:, 2);
    case 'z'
        tang1 = k(:, 1) + N + 1;
        tang2 = k(:, 2) + N + 1;
        kn = k(:, 3);
    otherwise
        error('Unsupported face type: %s', F.type);
end

Wmat = reshape(F.w(:), F.nq1, F.nq2);
w1 = Wmat(:, 1);
w2 = Wmat(1, :).' / Wmat(1, 1);
G1 = F.E1' * bsxfun(@times, F.E1, w1);
G2 = F.E2' * bsxfun(@times, F.E2, w2);

phase = exp(1i * alpha * F.fixedCoord * (kn.' - kn));
Pface = phase .* G1(tang1, tang1) .* G2(tang2, tang2) / Omega;

normalFactor = 1i * alpha * F.normalSign * kn;
AppFace = Pface .* (sigma + 0.25 * (normalFactor.' + conj(normalFactor)));
end

function [ApiFace, nGroups] = assemble_iga_pw_face_coupling_tensor_local( ...
F, pwData, sigma, interfaceIga)
%Assemble the IGA-PW face coupling.
alpha = pwData.alpha;
N = pwData.N;
OmegaSqrt = sqrt(pwData.Omega);
k = pwData.k_pw;
nPw = pwData.n_pw;
nGammaIga = numel(interfaceIga);
T = F.tensor;

W1B1 = spdiags(T.w1, 0, numel(T.w1), numel(T.w1)) * T.B1;
W2B2 = spdiags(T.w2, 0, numel(T.w2), numel(T.w2)) * T.B2;
C1 = F.E1' * W1B1;
C2 = F.E2' * W2B2;

switch F.type
    case 'x'
        tang1 = k(:, 2) + N + 1;
        tang2 = k(:, 3) + N + 1;
        kn = k(:, 1);
    case 'y'
        tang1 = k(:, 1) + N + 1;
        tang2 = k(:, 3) + N + 1;
        kn = k(:, 2);
    case 'z'
        tang1 = k(:, 1) + N + 1;
        tang2 = k(:, 2) + N + 1;
        kn = k(:, 3);
    otherwise
        error('Unsupported face type: %s', F.type);
end

ApiFace = complex(zeros(nPw, nGammaIga));
pairs = unique([tang1, tang2], 'rows');
nGroups = size(pairs, 1);

for ig = 1:nGroups
    i1 = pairs(ig, 1);
    i2 = pairs(ig, 2);
    rows = find(tang1 == i1 & tang2 == i2);

    switch F.type
        case 'x'
            baseTI = kron(C2(i2, :), kron(C1(i1, :), T.Bfix));
            baseGI = kron(C2(i2, :), kron(C1(i1, :), T.Dfix));
        case 'y'
            baseTI = kron(C2(i2, :), kron(T.Bfix, C1(i1, :)));
            baseGI = kron(C2(i2, :), kron(T.Dfix, C1(i1, :)));
        case 'z'
            baseTI = kron(T.Bfix, kron(C2(i2, :), C1(i1, :)));
            baseGI = kron(T.Dfix, kron(C2(i2, :), C1(i1, :)));
    end

    baseTI = full(baseTI(interfaceIga));
    baseGI = full(baseGI(interfaceIga));

    phase = exp(1i * alpha * F.fixedCoord * kn(rows)) / OmegaSqrt;
    coeff = conj(phase(:));
    normalFactor = 1i * alpha * F.normalSign * kn(rows);
    coeffTI = -(0.25 * conj(normalFactor(:)) + sigma) .* coeff;
    coeffGI = 0.25 * coeff;

    ApiFace(rows, :) = coeffGI * baseGI + coeffTI * baseTI;
end
end

function [Tm, normalFactor] = build_face_trace_mat_block_local(F, pwData, cols)
%Build face trace MAT block.
nq = numel(F.w);
alpha = pwData.alpha;
N = pwData.N;
OmegaSqrt = sqrt(pwData.Omega);
k = pwData.k_pw;

cols = cols(:).';
Tm = complex(zeros(nq, numel(cols)));

switch F.type
    case 'x'
        tang1 = k(cols, 2) + N + 1;
        tang2 = k(cols, 3) + N + 1;
        kn = k(cols, 1);
    case 'y'
        tang1 = k(cols, 1) + N + 1;
        tang2 = k(cols, 3) + N + 1;
        kn = k(cols, 2);
    case 'z'
        tang1 = k(cols, 1) + N + 1;
        tang2 = k(cols, 2) + N + 1;
        kn = k(cols, 3);
    otherwise
        error('Unsupported face type: %s', F.type);
end

phaseFix = exp(1i * alpha * F.fixedCoord * kn) / OmegaSqrt;
normalFactor = 1i * alpha * F.normalSign * kn;

for j = 1:numel(cols)
    traceCol = phaseFix(j) * kron(F.E2(:, tang2(j)), F.E1(:, tang1(j)));
    Tm(:, j) = traceCol;
end
end

function u = normalize_in_B(u, Bfun)
%Normalize in b.
nu = real(u' * Bfun(u));
u = u / sqrt(nu);
end

function u = align_phase(u, uref, Bfun)
%Align phase.
ov = uref' * Bfun(u);
if abs(ov) > 1e-14
    u = exp(-1i * angle(ov)) * u;
end
end

function [mode_name, mode_value, mode_tag] = parse_refinement_input(Refinement)
%Compute refinement input.
if isnumeric(Refinement) && isscalar(Refinement)
    mode_name = 'dyadic';
    mode_value = Refinement;
    mode_tag = sprintf('refine_%02d', Refinement);
elseif isstruct(Refinement) && isfield(Refinement, 'mode') && strcmpi(Refinement.mode, 'nelem')
    mode_name = 'nelem';
    mode_value = Refinement.value;
    if isscalar(mode_value)
        mode_tag = sprintf('nelem_%02d', mode_value);
    else
        mode_tag = sprintf('nelem_%02d_%02d_%02d', mode_value(1), mode_value(2), mode_value(3));
    end
else
    error('Unsupported Refinement input.');
end
end

function [uh, D, rnorms, stats, hist] = call_primme_operator( ...
Afun, Bfun, n, n_eigs, targetShift, ops, method, Pfun, saveHistory)
%Call PRIMME with operator handles.
arguments
Afun
Bfun
n
n_eigs
targetShift
ops
method
Pfun = []
saveHistory = false
end
use_prec = ~isempty(Pfun);
hist = [];

if saveHistory
    ops.reportLevel = max(ops.reportLevel, 2);
    ops.display = 0;
end

if use_prec
    if saveHistory
        [uh, D, rnorms, stats, hist] = primme_eigs( ...
            Afun, Bfun, n, n_eigs, targetShift, ops, method, Pfun);
    else
        [uh, D, rnorms, stats] = primme_eigs( ...
            Afun, Bfun, n, n_eigs, targetShift, ops, method, Pfun);
    end
else
    if saveHistory
        [uh, D, rnorms, stats, hist] = primme_eigs( ...
            Afun, Bfun, n, n_eigs, targetShift, ops, method);
    else
        [uh, D, rnorms, stats] = primme_eigs( ...
            Afun, Bfun, n, n_eigs, targetShift, ops, method);
    end
end
end

function uh0 = build_initial_guess_local(pwDataStatic, nI, nTotal)
%Build the constant hybrid initial state.
uh0 = build_constant_hybrid_initial_guess(pwDataStatic, nI, nTotal);
end

function uh0 = build_constant_hybrid_initial_guess(pwDataStatic, nI, nTotal)
%Build constant hybrid initial guess.
uh0 = zeros(nTotal, 1);
uh0(1:nI) = 1;

k0 = find(all(pwDataStatic.k_pw == 0, 2), 1, 'first');
if isempty(k0)
    error('The PW basis does not contain the zero mode needed for the constant initial guess.');
end

uh0(nI + k0) = sqrt(pwDataStatic.Omega);
end

function [pwData, timing, cacheFile, cacheHit] = get_pw_operator_cached_3D_local( ...
L, Nc, inner_domains, k_Vr, n_pw_Vr, opts)
%Return cached PW operators.

cacheHit = false;
cacheFile = '';
fft_grid_n = get_pw_fft_grid_n(opts, Nc);
inner_cheb_n = opts.inner_cheb_n;
nuclear_charge = opts.nuclear_charge;

if ~(isfield(opts, 'use_pw_cache') && opts.use_pw_cache && ...
        isfield(opts, 'cacheRoot') && ~isempty(opts.cacheRoot))
    [pwData, timing] = generate_A_M_PW_3D( ...
        L, Nc, inner_domains, k_Vr, n_pw_Vr, fft_grid_n, opts);
    return;
end

xL = inner_domains(1); xR = inner_domains(2);
yL = inner_domains(3); yR = inner_domains(4);
zL = inner_domains(5); zR = inner_domains(6);

cacheFile = fullfile(opts.cacheRoot, ...
    sprintf(['PW3D_OPERATOR_%s_L_%g_Nc_%d_NVr_%d_Z_%g_ic_%d_', ...
    'mFFT_%d_x_%g_%g_y_%g_%g_z_%g_%g.mat'], ...
    opts.Example, L, Nc, n_pw_Vr, ...
    nuclear_charge, inner_cheb_n, fft_grid_n, ...
    xL, xR, yL, yR, zL, zR));

if exist(cacheFile, 'file')
    S = load(cacheFile, 'pwData', 'timing');
    pwData = S.pwData;
    timing = S.timing;
    cacheHit = true;
    return;
end

[pwData, timing] = generate_A_M_PW_3D( ...
    L, Nc, inner_domains, k_Vr, n_pw_Vr, fft_grid_n, opts);
if ~exist(opts.cacheRoot, 'dir'), mkdir(opts.cacheRoot); end
save(cacheFile, 'pwData', 'timing', '-v7.3');
end

function [VextGrid, cacheFile, cacheHit] = get_vext_grid_cached_3D(gridCache, k_Vr, n_pw_Vr, opts)
%Return vext grid cached 3D.
cacheHit = false;
cacheFile = '';
L = gridCache.L;
mFFT = gridCache.mFFT;
nuclear_charge = opts.nuclear_charge;

if ~(isfield(opts, 'use_vext_cache') && opts.use_vext_cache && ...
        isfield(opts, 'cacheRoot') && ~isempty(opts.cacheRoot))
    VextGrid = build_vext_grid_dg_pw_iga_0414_3D(gridCache, k_Vr, n_pw_Vr, nuclear_charge);
    return;
end

cacheFile = fullfile(opts.cacheRoot, ...
    sprintf('VextGrid_DG0414_L_%g_mFFT_%d_NVr_%d_Z_%g_alpha_5.mat', L, mFFT, n_pw_Vr, nuclear_charge));

if exist(cacheFile, 'file')
    S = load(cacheFile, 'VextGrid');
    VextGrid = S.VextGrid;
    cacheHit = true;
    return;
end

VextGrid = build_vext_grid_dg_pw_iga_0414_3D(gridCache, k_Vr, n_pw_Vr, nuclear_charge);
save(cacheFile, 'VextGrid', '-v7.3');
end

function VextGrid = build_vext_grid_dg_pw_iga_0414_3D(gridCache, k_Vr, n_pw_Vr, nuclear_charge)
%Build vext grid DG PW IGA 0414 3D.
alpha = 5;
L = gridCache.L;
Omega = L ^ 3;
xmid = gridCache.xmid;
mFFT = gridCache.mFFT;

VextGrid = zeros(mFFT, mFFT, mFFT);
[Y, Z] = ndgrid(xmid, xmid);
YZ2 = Y .^ 2 + Z .^ 2;

G = (2 * pi / L) * k_Vr;
recip_coeff = zeros(n_pw_Vr, 1);
phase_yz = cell(n_pw_Vr, 1);
for ig = 1:n_pw_Vr
    Gnorm2 = G(ig, 1) ^ 2 + G(ig, 2) ^ 2 + G(ig, 3) ^ 2;
    if Gnorm2 > 0
        recip_coeff(ig) = -nuclear_charge * (4 * pi / Omega) * exp(-Gnorm2 / (4 * alpha ^ 2)) / Gnorm2;
        phase_yz{ig} = exp(1i * (G(ig, 2) * Y + G(ig, 3) * Z));
    end
end

for ix = 1:mFFT
    r = sqrt(xmid(ix) ^ 2 + YZ2);
    r = max(r, 1e-14);
    Vslice = -nuclear_charge * erfc(alpha * r) ./ r;

    for ig = 1:n_pw_Vr
        if recip_coeff(ig) ~= 0
            phase_x = exp(1i * G(ig, 1) * xmid(ix));
            Vslice = Vslice + recip_coeff(ig) * phase_x * phase_yz{ig};
        end
    end

    VextGrid(ix, :, :) = real(Vslice + nuclear_charge * 2 * alpha / sqrt(pi));
end
end

function gridCache = build_midpoint_grid_cache_3D(L, mFFT, a)
%Build midpoint grid cache 3D.
dx = L / mFFT;
xmid = -L / 2 + dx / 2 + (0:mFFT-1) * dx;
inner_idx = find(abs(xmid) <= a + 1e-12);
x_inner = xmid(inner_idx);
[Xin, Yin, Zin] = ndgrid(x_inner, x_inner, x_inner);

gridCache = struct();
gridCache.L = L;
gridCache.a = a;
gridCache.dx = dx;
gridCache.mFFT = mFFT;
gridCache.xmid = xmid;
gridCache.inner_idx = inner_idx;
gridCache.n_inner = numel(inner_idx);
gridCache.x_inner = Xin(:);
gridCache.y_inner = Yin(:);
gridCache.z_inner = Zin(:);
end

function fft_grid_n = get_pw_fft_grid_n(opts, Nc)
%Return PW FFT grid n.
fft_grid_n = opts.pw_fft_grid_n;
fft_grid_n = max(fft_grid_n, 4 * floor(Nc) + 1);
if mod(fft_grid_n, 2) ~= 0
    fft_grid_n = fft_grid_n + 1;
end
fft_grid_n = max(fft_grid_n, 8);
end

function hartree_grid_n = get_hartree_grid_n(opts)
%Return hartree grid n.
hartree_grid_n = opts.hartree_grid_n;
hartree_grid_n = max(round(hartree_grid_n), 8);
if mod(hartree_grid_n, 2) ~= 0
    hartree_grid_n = hartree_grid_n + 1;
end
end

function err = compute_reference_error(value, referenceValue)
%Compute reference error.
if ~isempty(referenceValue) && isfinite(referenceValue)
    err = abs(value - referenceValue);
else
    err = [];
end
end

function Etot = compute_total_hartree_energy(uh, AstaticFun, int_rho_vh)
%Compute total hartree energy.
singleParticle = real(uh' * AstaticFun(uh));
Etot = 2 * singleParticle + 0.5 * int_rho_vh;
end

function info = save_density_diagnostics( ...
outDir, uh, nI, nurbs_refine, pwData, rhoGrid, VHGrid, gridCache, a, lambdaHist, rhoResHist)
%Save density diagnostics.
xmid = gridCache.xmid;
uLine = evaluate_hybrid_state_3D( ...
    uh, nI, nurbs_refine, pwData, a, xmid(:), zeros(numel(xmid), 1), zeros(numel(xmid), 1));
rhoX = 2 * abs(uLine(:)) .^ 2;
vhX = extract_x_axis_field(VHGrid, xmid);
rhoZ0 = extract_z0_slice_field(rhoGrid, xmid);
iter = (1:numel(lambdaHist)).';

diagData = struct();
diagData.iter = iter;
diagData.lambda = lambdaHist(:);
diagData.rho_res = rhoResHist(:);
diagData.x = xmid(:);
diagData.u_x_axis = uLine(:);
diagData.rho_x_axis = real(rhoX(:));
diagData.VH_x_axis = real(vhX(:));
diagData.slice_x = xmid(:);
diagData.slice_y = xmid(:);
diagData.rho_z0 = real(rhoZ0);

dataFile = fullfile(outDir, 'density_diagnostics.mat');
scfFile = fullfile(outDir, 'scf_convergence.png');
lineFile = fullfile(outDir, 'line_profile_x.png');
sliceFile = fullfile(outDir, 'rho_z0_slice.png');

save(dataFile, 'diagData', '-v7.3');

fig1 = figure('Visible', 'off', 'Color', 'w');
tl = tiledlayout(fig1, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax1 = nexttile(tl);
semilogy(ax1, iter, max(real(rhoResHist(:)), eps), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel(ax1, 'SCF iteration', 'Interpreter', 'none');
ylabel(ax1, 'relative density residual', 'Interpreter', 'none');
title(ax1, 'SCF residual', 'Interpreter', 'none');
grid(ax1, 'on');

ax2 = nexttile(tl);
plot(ax2, iter, real(lambdaHist(:)), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel(ax2, 'SCF iteration', 'Interpreter', 'none');
ylabel(ax2, '\lambda_1', 'Interpreter', 'tex');
title(ax2, 'Lowest eigenvalue', 'Interpreter', 'none');
grid(ax2, 'on');
exportgraphics(fig1, scfFile, 'Resolution', 200);
close(fig1);

fig2 = figure('Visible', 'off', 'Color', 'w');
yyaxis left;
plot(xmid, real(uLine), 'LineWidth', 1.5);
hold on;
plot(xmid, real(rhoX), 'LineWidth', 1.5);
xlabel('x', 'Interpreter', 'none');
ylabel('u(x,0,0), \rho(x,0,0)', 'Interpreter', 'tex');
yyaxis right;
plot(xmid, real(vhX), '--', 'LineWidth', 1.2);
ylabel('V_H(x,0,0)', 'Interpreter', 'tex');
title('Helium line profile on the x-axis', 'Interpreter', 'none');
legend({'u(x,0,0)', '\rho(x,0,0)', 'V_H(x,0,0)'}, ...
    'Location', 'best', 'Interpreter', 'tex');
grid on;
exportgraphics(fig2, lineFile, 'Resolution', 200);
close(fig2);

fig3 = figure('Visible', 'off', 'Color', 'w');
imagesc(xmid, xmid, real(rhoZ0).');
set(gca, 'YDir', 'normal');
axis equal tight;
xlabel('x', 'Interpreter', 'none');
ylabel('y', 'Interpreter', 'none');
title('Helium density slice on z = 0', 'Interpreter', 'none');
colorbar;
exportgraphics(fig3, sliceFile, 'Resolution', 200);
close(fig3);

info = struct();
info.data_file = dataFile;
info.scf_convergence_plot = scfFile;
info.x_axis_plot = lineFile;
info.z0_slice_plot = sliceFile;
info.figure_files = {scfFile, lineFile, sliceFile};
end

function val = evaluate_hybrid_state_3D(uh, nI, nurbs_refine, pwData, a, X, Y, Z)
%Evaluate the 3-D hybrid solution state.
Xv = X(:);
Yv = Y(:);
Zv = Z(:);

val = eval_pw_on_points_local(uh(nI+1:end), pwData, Xv, Yv, Zv);
innerMask = abs(Xv) <= a + 1e-12 & abs(Yv) <= a + 1e-12 & abs(Zv) <= a + 1e-12;
if any(innerMask)
    val(innerMask) = eval_iga_on_points_local( ...
        nurbs_refine, uh(1:nI), Xv(innerMask), Yv(innerMask), Zv(innerMask), a);
end

val = reshape(val, size(X));
end

function val = eval_pw_on_points_local(cP, pwData, X, Y, Z)
%Evaluate PW on points.
k_pw = pwData.k_pw;
alpha = pwData.alpha;
Omega = pwData.Omega;
nq = numel(X);
val = zeros(nq, 1);
batchSize = 256;

for q0 = 1:batchSize:nq
    q1 = min(q0 + batchSize - 1, nq);
    idx = q0:q1;
    phase = exp(1i * alpha * ( ...
        X(idx) * k_pw(:, 1).' + ...
        Y(idx) * k_pw(:, 2).' + ...
        Z(idx) * k_pw(:, 3).')) / sqrt(Omega);
    val(idx) = phase * cP(:);
end
end

function val = eval_iga_on_points_local(nurbs_refine, coeff, X, Y, Z, a)
%Evaluate IGA on points.
U = nurbs_refine.Ubar;
V = nurbs_refine.Vbar;
W = nurbs_refine.Wbar;
pu = nurbs_refine.pu;
pv = nurbs_refine.pv;
pw = nurbs_refine.pw;
m = nurbs_refine.m;
n = nurbs_refine.n;

nq = numel(X);
val = zeros(nq, 1);

for q = 1:nq
    u = min(max((X(q) + a) / (2 * a), 0), 1);
    v = min(max((Y(q) + a) / (2 * a), 0), 1);
    w = min(max((Z(q) + a) / (2 * a), 0), 1);

    ispan = findspan(U, pu, u);
    jspan = findspan(V, pv, v);
    kspan = findspan(W, pw, w);

    Nu = bspbasisDers(U, pu, u, 1); Nu = Nu(1, :)';
    Nv = bspbasisDers(V, pv, v, 1); Nv = Nv(1, :)';
    Nw = bspbasisDers(W, pw, w, 1); Nw = Nw(1, :)';

    ii = ispan-pu:ispan;
    jj = jspan-pv:jspan;
    kk = kspan-pw:kspan;

    s = 0;
    for kz = 1:(pw + 1)
        for jy = 1:(pv + 1)
            for ix = 1:(pu + 1)
                gid = ii(ix) + (jj(jy) - 1) * m + (kk(kz) - 1) * m * n;
                s = s + coeff(gid) * Nu(ix) * Nv(jy) * Nw(kz);
            end
        end
    end
    val(q) = s;
end
end

function rhoX = extract_x_axis_field(rhoGrid, xmid)
%Extract x axis field.
[iy0, iy1, wy0, wy1] = interpolation_weights_for_point(xmid, 0);
[iz0, iz1, wz0, wz1] = interpolation_weights_for_point(xmid, 0);

rhoX = wy0 * wz0 * rhoGrid(:, iy0, iz0);
rhoX = rhoX + wy0 * wz1 * rhoGrid(:, iy0, iz1);
rhoX = rhoX + wy1 * wz0 * rhoGrid(:, iy1, iz0);
rhoX = rhoX + wy1 * wz1 * rhoGrid(:, iy1, iz1);
rhoX = rhoX(:);
end

function rhoZ0 = extract_z0_slice_field(rhoGrid, xmid)
%Extract z0 slice field.
[iz0, iz1, wz0, wz1] = interpolation_weights_for_point(xmid, 0);
rhoZ0 = wz0 * rhoGrid(:, :, iz0) + wz1 * rhoGrid(:, :, iz1);
end

function [idx0, idx1, w0, w1] = interpolation_weights_for_point(x, xq)
%Compute weights for point.
tol = 10 * eps(max(1, max(abs(x))));
idxExact = find(abs(x - xq) <= tol, 1, 'first');

if ~isempty(idxExact)
    idx0 = idxExact;
    idx1 = idxExact;
    w0 = 1;
    w1 = 0;
    return;
end

idx1 = find(x > xq, 1, 'first');
idx0 = find(x < xq, 1, 'last');

if isempty(idx0)
    idx0 = idx1;
    w0 = 1;
    w1 = 0;
    return;
end

if isempty(idx1)
    idx1 = idx0;
    w0 = 1;
    w1 = 0;
    return;
end

w1 = (xq - x(idx0)) / (x(idx1) - x(idx0));
w0 = 1 - w1;
end
