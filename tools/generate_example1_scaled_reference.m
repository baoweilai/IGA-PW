function generate_example1_scaled_reference()
% Generate the reference solution used by Example 1 scaled-error plots.

projectDir = fileparts(fileparts(mfilename('fullpath')));
add_required_paths(projectDir);

exampleDir = fullfile(projectDir, 'examples', 'example_1');
add_example_paths(exampleDir);
add_workflow_paths(fullfile(exampleDir, 'model', 'scaled_errors'), ...
    {'nurbs', 'dg', 'iga', 'assembly', 'operators', 'error_norms', 'core'});

dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

outDir = fullfile(dataDir, 'result', 'Example_1_verify49', ...
    'reference', 'p_3', 'refine_06_Nc_48');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cacheRoot = fullfile(dataDir, 'result', 'scaled_errors', 'cache_pw');
if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

opts = struct();
opts.Example = 'Example_1';
opts.beta = 20;
opts.n_gp = 10;
opts.inner_cheb_n = 200;
opts.inner_quad_n = 1000;
opts.pw_fft_grid_n = 256;
opts.primme_tol = 1e-12;
opts.primme_maxit = 5e7;
opts.primme_method = 'DEFAULT_MIN_TIME';
opts.primme_reportLevel = 0;
opts.eps_diag = 1e-12;
opts.iface_reg = 1e-12;
opts.save_eigenvectors = true;
opts.save_nurbs = true;
opts.save_pw_index = true;
opts.use_pw_cache = true;
opts.cacheRoot = cacheRoot;
opts.outDir = outDir;

solve_iga_pw_dg(6, 2, 48, 1, opts);
end

function add_required_paths(projectDir)
% Add project folders needed to build the reference solution.

requiredDirs = {
    fullfile(projectDir, 'src', 'assembly')
    fullfile(projectDir, 'src', 'dg')
    fullfile(projectDir, 'src', 'iga')
    fullfile(projectDir, 'src', 'pw')
    fullfile(projectDir, 'src', 'nurbs')
    fullfile(projectDir, 'src', 'error_norms')
    fullfile(projectDir, 'src', 'postprocess')
    fullfile(projectDir, 'src', 'plotting')
    fullfile(projectDir, 'src', 'utils')
    fullfile(projectDir, 'external', 'primme', 'Matlab')
    };

for k = 1:numel(requiredDirs)
    assert(isfolder(requiredDirs{k}), 'Missing project directory: %s', requiredDirs{k});
    addpath(requiredDirs{k});
end
end
