function cfg = default_config(userCfg)
%Return fixed model parameters.
arguments
    userCfg struct = struct()
end

cfg = struct();
cfg.smoke = false;
cfg.force = false;

cfg.reference = struct('p', 2, 'Nc', 45, 'Nelement', 32);
cfg.pw = struct('fixed_p', 1, 'fixed_Nelement', 12, 'Nc_list', 4:8);
cfg.iga = struct('fixed_Nc', 20, 'p_list', [1 2], ...
    'Nelement_list', [2 4 8 12]);
cfg.energy = struct('fixed_p', [2], ...
    'K_list', [2, 4, 8, 12, 16], ...
    'refine_list', [1, 2, 3, 4, 5]);

cfg.reference_lambda = [];
cfg.reference_energy = [];

cfg.inner_cheb_n = 80;
cfg.pw_fft_grid_n = 300;
cfg.hartree_grid_n = 300;
cfg.scf_maxit = 40;
cfg.scf_tol_eig = 1e-7;
cfg.scf_tol_rho = 1e-7;
cfg.scf_stopping_rule = 'lambda_and_rho';
cfg.scf_beta = 0.8;
cfg.primme_tol = 1e-9;
cfg.preconditioner_type = 'blockdiag_jacobi';
cfg.iface_explicit_gamma_max = 12000;
cfg.state_error_grid_n = 100;

cfg = merge_struct_local(cfg, userCfg);
end

function A = merge_struct_local(A, B)
%Merge struct.
if ~isstruct(B)
    return;
end
fn = fieldnames(B);
for i = 1:numel(fn)
    key = fn{i};
    if isfield(A, key) && isstruct(A.(key)) && isstruct(B.(key))
        A.(key) = merge_struct_local(A.(key), B.(key));
    else
        A.(key) = B.(key);
    end
end
end
