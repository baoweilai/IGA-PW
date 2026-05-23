function out = run_h_data(pList)
%Generate h-convergence data.
assert(exist('pList', 'var') == 1, 'run_h_data requires pList.');

clc;
activate_example_workflow('h_convergence', {'config', 'core', 'operators', 'solver'});
cfg = struct('force', true, 'primme_tol', 1e-9, ...
    'pw_fft_grid_n', 300, 'hartree_grid_n', 300, ...
    'save_eigenvectors', true, 'save_nurbs', true, 'save_pw_index', true);
nElemList = [2 4 8 12];
Nc = 20;

out = cell(numel(pList), numel(nElemList));
for ip = 1:numel(pList)
    for in = 1:numel(nElemList)
        out{ip, in} = run_iga(pList(ip), nElemList(in), Nc, true, cfg);
    end
end
end
