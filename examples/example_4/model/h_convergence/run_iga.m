function out = run_iga(pdeg, nElem, Nc)
% Run the IGA-only convergence calculation.

assert(exist('pdeg', 'var') == 1, 'run_iga requires pdeg.');
assert(exist('nElem', 'var') == 1, 'run_iga requires nElem.');
assert(exist('Nc', 'var') == 1, 'run_iga requires Nc.');

% Load the workflow and run the requested IGA case.
activate_example_workflow('h_convergence', {'config', 'core', 'operators', 'solver'});
cfg = default_config(struct());
H = example_helpers(cfg);
out = H.run_case('iga', pdeg, Nc, nElem);
end
