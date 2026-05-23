function out = run_reference_solution(forceFlag, userCfg)
%Generate the reference solution data.
assert(exist('forceFlag', 'var') == 1, 'run_reference_solution requires forceFlag.');
assert(exist('userCfg', 'var') == 1, 'run_reference_solution requires userCfg.');

activate_example_workflow('reference', {'config', 'core', 'operators', 'solver'});
cfg = default_config(userCfg);
cfg.force = logical(forceFlag);
H = example_helpers(cfg);
ref = H.cfg.reference;
caseOut = H.run_case('reference', ref.p, ref.Nc, ref.Nelement);

lambda_ref = caseOut.lambda;
uh_ref = caseOut.run.uh;
rho_ref = caseOut.run.rhoGrid;
scf_history = caseOut.run.scf_history;

p_ref = ref.p;
Nc_ref = ref.Nc;
Nelement_ref = ref.Nelement;
refRunFile = caseOut.runFile;
config = H.cfg;
save(H.cfg.reference_result_file, 'lambda_ref', 'uh_ref', 'rho_ref', ...
    'scf_history', 'config', 'p_ref', 'Nc_ref', 'Nelement_ref', ...
    'refRunFile', '-v7.3');

out = load(H.cfg.reference_result_file);
out.run_file = refRunFile;
end
