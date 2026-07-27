function out = run_h_data(pList)
% Generate h-convergence data.
assert(exist('pList', 'var') == 1, 'run_h_data requires pList.');

% Load the workflow and set the retained mesh cases.
clc;
activate_example_workflow('h_convergence', {'config', 'core', 'operators', 'solver'});
nElemList = [2 4 8 16];
Nc = 20;

% Run every requested degree and mesh size.
out = cell(numel(pList), numel(nElemList));
for ip = 1:numel(pList)
    for in = 1:numel(nElemList)
        out{ip, in} = run_iga(pList(ip), nElemList(in), Nc);
    end
end
end
