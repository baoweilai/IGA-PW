function out = run_pw(Nc, pdeg, nElem, forceFlag, userCfg)
%Run PW.

assert(exist('Nc', 'var') == 1, 'run_pw requires Nc.');
assert(exist('pdeg', 'var') == 1, 'run_pw requires pdeg.');
assert(exist('nElem', 'var') == 1, 'run_pw requires nElem.');
assert(exist('forceFlag', 'var') == 1, 'run_pw requires forceFlag.');
assert(exist('userCfg', 'var') == 1, 'run_pw requires userCfg.');

activate_example_workflow('cutoff_convergence', {'config', 'core', 'operators', 'solver'});
cfg = default_config(userCfg);

cfg.force = logical(forceFlag);
H = example_helpers(cfg);
out = H.run_case('pw', pdeg, Nc, nElem);
end
