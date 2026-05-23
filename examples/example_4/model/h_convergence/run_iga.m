function out = run_iga(pdeg, nElem, Nc, forceFlag, userCfg)
%Run the IGA-only convergence calculation.

assert(exist('pdeg', 'var') == 1, 'run_iga requires pdeg.');
assert(exist('nElem', 'var') == 1, 'run_iga requires nElem.');
assert(exist('Nc', 'var') == 1, 'run_iga requires Nc.');
assert(exist('forceFlag', 'var') == 1, 'run_iga requires forceFlag.');
assert(exist('userCfg', 'var') == 1, 'run_iga requires userCfg.');

activate_example_workflow('h_convergence', {'config', 'core', 'operators', 'solver'});
cfg = default_config(userCfg);

cfg.force = logical(forceFlag);
H = example_helpers(cfg);
out = H.run_case('iga', pdeg, Nc, nElem);
end
