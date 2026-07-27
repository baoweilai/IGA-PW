function out = run_pw(Nc, pdeg, nElem)
% Run one full-domain plane-wave reference case.

assert(exist('Nc', 'var') == 1, 'run_pw requires Nc.');
assert(exist('pdeg', 'var') == 1, 'run_pw requires pdeg.');
assert(exist('nElem', 'var') == 1, 'run_pw requires nElem.');

% Load the workflow and run the requested plane-wave case.
activate_example_workflow('cutoff_convergence', {'config', 'core', 'operators', 'solver'});
cfg = default_config(struct());
H = example_helpers(cfg);
out = H.run_case('pw', pdeg, Nc, nElem);
end
