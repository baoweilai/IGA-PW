function out = run_energy_data()
% Generate energy-convergence data.

% Load the energy configuration and reference values.
activate_example_workflow('energy', {'config', 'core', 'operators', 'solver'});
cfg = default_config(struct());
[referenceLambda, referenceEnergy] = load_reference_values_local();
if isempty(cfg.reference_lambda)
    cfg.reference_lambda = referenceLambda;
end
if isempty(cfg.reference_energy)
    cfg.reference_energy = referenceEnergy;
end
H = example_helpers(cfg);

energyCfg = H.cfg.energy;
K_list = reshape(energyCfg.K_list, 1, []);
refine_list = reshape(energyCfg.refine_list, 1, []);
p_list = reshape(energyCfg.p_list, 1, []);
if numel(K_list) ~= numel(refine_list)
    error('energy.K_list and energy.refine_list must have the same length.');
end

nCase = numel(K_list);
nP = numel(p_list);
lambda_list = zeros(nP, nCase);
energy_total = zeros(nP, nCase);
time_total = zeros(nP, nCase);
partial_file = fullfile(H.cfg.resultRoot, 'energy_runs_partial.mat');

% Run every polynomial degree and refinement case.
for ip = 1:nP
    pdeg = p_list(ip);
    for i = 1:nCase
        caseOut = H.run_refinement_case( ...
            'energy', pdeg, K_list(i), refine_list(i));
        lambda_list(ip, i) = caseOut.lambda;

        assert(isfield(caseOut.run, 'meta') ...
            && isfield(caseOut.run.meta, 'energy_total'), ...
            'Missing energy_total in %s.', caseOut.runFile);
        meta = caseOut.run.meta;
        energy_total(ip, i) = meta.energy_total;
        time_total(ip, i) = meta.time_total;
        validate_case_local(meta, H.cfg, pdeg, K_list(i), refine_list(i));
        save(partial_file, 'p_list', 'K_list', 'refine_list', ...
            'lambda_list', 'energy_total', 'time_total', ...
            'referenceLambda', 'referenceEnergy');
        clear caseOut meta
    end
end

% Compute the reference errors and save the completed dataset.
reference_lambda = H.cfg.reference_lambda;
reference_energy = H.cfg.reference_energy;
assert(isfinite(reference_lambda), 'Reference eigenvalue is not finite.');
assert(isfinite(reference_energy), 'Reference energy is not finite.');
err_lambda = abs(lambda_list - reference_lambda);
err_energy = abs(energy_total - reference_energy);

result_file = fullfile(H.cfg.resultRoot, 'energy_runs.mat');
save(result_file, 'p_list', 'K_list', 'refine_list', ...
    'lambda_list', 'energy_total', 'time_total', ...
    'reference_lambda', 'reference_energy', ...
    'err_lambda', 'err_energy');
if isfile(partial_file)
    delete(partial_file);
end

out = load(result_file);
out.result_file = result_file;
end

function [referenceLambda, referenceEnergy] = load_reference_values_local()
% Read the strict production reference scalars without loading its field grids.

workflowDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(fileparts(workflowDir));
refFile = fullfile(exampleDir, 'data', 'result', ...
    'REFERENCE', 'K_30', 'p_2', ...
    'nelem_32', 'run.mat');
assert(isfile(refFile), 'Reference run file does not exist: %s', refFile);
referenceLambda = double(h5read(refFile, '/run/lambda'));
referenceEnergy = double(h5read(refFile, '/run/meta/energy_total'));
referenceLambda = real(referenceLambda(1));
referenceEnergy = real(referenceEnergy(1));

referenceDof = double(h5read(refFile, '/run/n_dofs_total'));
referenceK = double(h5read(refFile, '/run/meta/Nc'));
referenceP = double(h5read(refFile, '/run/meta/pdeg'));
referenceNelem = double(h5read(refFile, '/run/meta/refine_value'));
referenceConverged = double(h5read(refFile, '/run/meta/scf_converged'));
referenceIters = double(h5read(refFile, '/run/meta/scf_iters'));
referenceAbs = double(h5read(refFile, ...
    '/run/meta/final_lambda_abs_change'));
referenceRho = double(h5read(refFile, ...
    '/run/meta/final_density_residual'));
referenceInnerCheb = double(h5read(refFile, '/run/meta/inner_cheb_n'));
referencePwGrid = double(h5read(refFile, '/run/meta/pw_fft_grid_n'));
referenceHartreeGrid = double(h5read(refFile, ...
    '/run/meta/hartree_grid_n'));
referenceTolEig = double(h5read(refFile, '/run/meta/scf_tol_eig'));
referenceTolRho = double(h5read(refFile, '/run/meta/scf_tol_rho'));
referencePrimme = double(h5read(refFile, '/run/meta/primme_tol'));
referenceNGamma = double(h5read(refFile, ...
    '/run/meta/preconditioner_info/n_gamma'));

assert(referenceDof == 152385 && referenceK == 30 ...
    && referenceP == 2 && referenceNelem == 32, ...
    'The production reference identity or DOF is inconsistent.');
assert(referenceConverged == 1 && referenceIters == 14 ...
    && referenceAbs < 1e-10 && referenceRho < 1e-10, ...
    'The production reference did not meet strict SCF convergence.');
assert(referenceInnerCheb == 80 && referencePwGrid == 300 ...
    && referenceHartreeGrid == 300, ...
    'The production reference integration grids are inconsistent.');
assert(referenceTolEig == 1e-10 && referenceTolRho == 1e-10 ...
    && referencePrimme == 1e-12, ...
    'The production reference tolerances are inconsistent.');
assert(referenceNGamma == 0, ...
    'The production reference must use n_gamma=0.');
end

function validate_case_local(meta, cfg, pdeg, K, refine)
% Validate one strict energy point before proceeding to the next one.

assert(meta.scf_converged, ...
    'Energy case K=%d p=%d refine=%d did not converge.', K, pdeg, refine);
assert(meta.Nc == K && meta.pdeg == pdeg ...
    && strcmp(meta.refine_mode, 'dyadic') && meta.refine_value == refine, ...
    'Saved energy-case identity is inconsistent.');
assert(meta.scf_iters <= cfg.scf_maxit, 'SCF iteration count is invalid.');
assert(meta.final_lambda_abs_change < cfg.scf_tol_eig, ...
    'Final absolute eigenvalue change does not meet the strict tolerance.');
assert(meta.final_density_residual < cfg.scf_tol_rho, ...
    'Final density residual does not meet the strict tolerance.');
assert(meta.scf_tol_eig == cfg.scf_tol_eig ...
    && meta.scf_tol_rho == cfg.scf_tol_rho ...
    && meta.primme_tol == cfg.primme_tol, ...
    'Saved solver tolerances do not match the strict configuration.');
assert(strcmp(meta.scf_stopping_rule, 'lambda_and_rho'), ...
    'Energy cases must use the lambda-and-density stopping rule.');
assert(strcmp(meta.preconditioner_type, 'blockdiag_jacobi'), ...
    'Energy cases must use blockdiag_jacobi.');
assert(meta.preconditioner_info.n_gamma == 0, ...
    'Energy cases must use n_gamma=0.');
assert(strcmp(meta.vext_assembly, 'iga_gauss_direct_ewald'), ...
    'Energy cases must use direct-Gauss external potential assembly.');
assert(meta.inner_cheb_n == cfg.inner_cheb_n ...
    && meta.pw_fft_grid_n == cfg.pw_fft_grid_n ...
    && meta.hartree_grid_n == cfg.hartree_grid_n ...
    && meta.grid_mFFT == cfg.hartree_grid_n ...
    && meta.global_fft_grid_n == cfg.hartree_grid_n, ...
    'Saved energy-case integration grids do not match the strict configuration.');
assert(abs(meta.reference_lambda - cfg.reference_lambda) ...
    <= 10 * eps(max(1, abs(cfg.reference_lambda))) ...
    && abs(meta.reference_energy - cfg.reference_energy) ...
    <= 10 * eps(max(1, abs(cfg.reference_energy))), ...
    'Saved energy-case reference values are inconsistent.');
assert(isfinite(meta.energy_total) && isfinite(meta.time_total) ...
    && meta.time_total > 0, ...
    'Saved energy or runtime is invalid.');
expectedH = 0.4 / (2 ^ refine);
assert(abs(meta.hmin - expectedH) <= 1e-12 * max(1, expectedH), ...
    'Saved hmin is inconsistent with h=0.4/2^r.');
end
