function out = run_energy_data(forceFlag, userCfg)
%Generate energy-convergence data.
assert(exist('forceFlag', 'var') == 1, 'run_energy_data requires forceFlag.');
assert(exist('userCfg', 'var') == 1, 'run_energy_data requires userCfg.');

activate_example_workflow('energy', {'config', 'core', 'operators', 'solver'});
cfg = default_config(userCfg);
cfg.force = logical(forceFlag);
H = example_helpers(cfg);

energyCfg = H.cfg.energy;
K_list = reshape(energyCfg.K_list, 1, []);
refine_list = reshape(energyCfg.refine_list, 1, []);
if numel(K_list) ~= numel(refine_list)
    error('energy.K_list and energy.refine_list must have the same length.');
end

pdeg = energyCfg.fixed_p;
nCase = numel(K_list);
lambda_list = zeros(1, nCase);
energy_total = zeros(1, nCase);

for i = 1:nCase
    caseOut = H.run_refinement_case('energy', pdeg, K_list(i), refine_list(i));
    lambda_list(i) = caseOut.lambda;

    assert(isfield(caseOut.run, 'meta') && isfield(caseOut.run.meta, 'energy_total'), ...
        'Missing energy_total in %s.', caseOut.runFile);
    energy_total(i) = caseOut.run.meta.energy_total;
end

reference_energy = H.cfg.reference_energy;
if isempty(reference_energy)
    reference_energy = load_reference_energy_local(H.cfg.resultRoot);
end
assert(isfinite(reference_energy), 'Reference energy is not finite.');
err_energy = abs(energy_total - reference_energy);

result_file = fullfile(H.cfg.resultRoot, 'energy_runs.mat');
save(result_file, 'K_list', 'refine_list', 'pdeg', ...
    'lambda_list', 'energy_total', ...
    'reference_energy', 'err_energy', '-v7.3');

out = load(result_file);
out.result_file = result_file;
end

function referenceEnergy = load_reference_energy_local(resultRoot)
%Read the reference total energy.

refFile = fullfile(resultRoot, 'REFERENCE', 'K_45', 'p_2', ...
    'nelem_32', 'run.mat');
referenceEnergy = double(h5read(refFile, '/run/meta/energy_total'));
end
