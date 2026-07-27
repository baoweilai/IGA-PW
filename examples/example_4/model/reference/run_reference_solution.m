function out = run_reference_solution()
% Generate the reference solution data.

% Load the reference configuration and solve its case.
activate_example_workflow('reference', {'config', 'core', 'operators', 'solver'});
cfg = default_config(struct());
H = example_helpers(cfg);
ref = H.cfg.reference;
caseOut = H.run_case('reference', ref.p, ref.Nc, ref.Nelement);

lambda_ref = caseOut.lambda;

% Save the reference values, metadata, and run location.
p_ref = ref.p;
Nc_ref = ref.Nc;
Nelement_ref = ref.Nelement;
refRunFile = caseOut.runFile;
config = H.cfg;
save(H.cfg.reference_result_file, 'lambda_ref', 'config', ...
    'p_ref', 'Nc_ref', 'Nelement_ref', ...
    'refRunFile', '-v7.3');

out = load(H.cfg.reference_result_file);
out.run_file = refRunFile;
end
