function out = run_cutoff_data()
% Generate cutoff-convergence data.

% Load the cutoff configuration and helper functions.
activate_example_workflow('cutoff_convergence', {'config', 'core', 'operators', 'solver'});
cfg = default_config(struct());
H = example_helpers(cfg);
pw = H.cfg.pw;

Nc_list = reshape(pw.Nc_list, 1, []);
nCase = numel(Nc_list);
lambda_list = zeros(1, nCase);

% Run every retained plane-wave cutoff case.
for i = 1:nCase
    caseOut = H.run_case('pw', pw.fixed_p, Nc_list(i), pw.fixed_Nelement);
    lambda_list(i) = caseOut.lambda;
end

% Compare with the reference and save the result.
lambda_ref = H.cfg.reference_lambda;
reference_run_file = H.case_run_file('reference', ...
    H.cfg.reference.p, H.cfg.reference.Nc, H.cfg.reference.Nelement);
if isempty(lambda_ref)
    S = load(H.cfg.reference_result_file, 'lambda_ref');
    assert(isfield(S, 'lambda_ref'), 'Missing lambda_ref in %s.', H.cfg.reference_result_file);
    lambda_ref = S.lambda_ref;
end
assert(isfinite(lambda_ref), 'Reference eigenvalue is not finite.');
err_lambda = abs(lambda_list - lambda_ref);

fixed_p = pw.fixed_p;
fixed_Nelement = pw.fixed_Nelement;
result_file = H.cfg.pw_result_file;
save(result_file, 'Nc_list', 'fixed_p', 'fixed_Nelement', ...
    'lambda_list', 'lambda_ref', 'err_lambda', ...
    'reference_run_file', '-v7.3');

out = load(result_file);
out.result_file = result_file;
end
