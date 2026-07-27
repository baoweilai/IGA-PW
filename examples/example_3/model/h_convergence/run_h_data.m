function run_h_data()
% Generate h-convergence data.

clc; close all;

% Set project paths and the workflow data directory.
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
projectDir = fileparts(fileparts(exampleDir));
utilsDir = fullfile(projectDir, 'src', 'utils');
assert(isfolder(utilsDir), 'Missing project utils directory: %s', utilsDir);
addpath(utilsDir, '-begin');
primmeDir = fullfile(projectDir, 'external', 'primme', 'Matlab');
assert(isfolder(primmeDir), 'Missing PRIMME directory: %s', primmeDir);
addpath(primmeDir, '-begin');

activate_example_workflow('h_convergence', ...
    {'nurbs', 'iga', 'assembly', 'operators', 'core'});
dataDir = fullfile(exampleDir, 'data');
startDir = pwd;
cleanupObj = onCleanup(@() cd(startDir));
cd(dataDir);

% Set the h-convergence and SCF parameters.
Example    = 'Example_3';
Nc         = 20;
nElem_list = [4 8 16 32 64];
t_list     = 0:1;   % p = 1, 2

beta = 100;
n_gp = 30;
inner_cheb_n = 150;
pw_fft_grid_n = 550;

primme_tol         = 1e-9;
primme_maxit       = 1e8;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
block_targetShift  = 0.0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

scf_maxit        = 60;
scf_pw_grid_m    = 500;
scf_tol_lambda   = 1e-7;
scf_mixing       = 0.9;
scf_track_n_eigs = 1;

save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;

use_pw_cache    = true;
use_nurbs_cache = true;

resultRoot = fullfile(pwd, 'result', Example, sprintf('Nc_%02d', Nc));
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cachePwRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cachePwRoot, 'dir'), mkdir(cachePwRoot); end

cacheNurbsRoot = fullfile(resultRoot, 'cache_nurbs');
if ~exist(cacheNurbsRoot, 'dir'), mkdir(cacheNurbsRoot); end

% Solve or load each degree and mesh case.
for t = t_list
    pdeg = 1 + t;
    caseDir = fullfile(resultRoot, sprintf('p_%d', pdeg));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    ncase = numel(nElem_list);
    h_list       = zeros(ncase, 1);
    lambda1_list = zeros(ncase, 1);
    dof_list     = zeros(ncase, 1);

    fprintf('\n============================================================\n');
    fprintf('[INFO] Example = %s, Nc = %d, p = %d\n', Example, Nc, pdeg);
    fprintf('[INFO] nElem_list = %s\n', mat2str(nElem_list));

    for ii = 1:ncase
        nElem = nElem_list(ii);
        runDir = fullfile(caseDir, sprintf('nElem_%02d', nElem));
        if ~exist(runDir, 'dir'), mkdir(runDir); end
        runMat = fullfile(runDir, 'run.mat');

        if exist(runMat, 'file')
            S = load(runMat, 'run');
            run = S.run;
            meta = run.meta;

            fprintf('[LOAD] p=%d, Nc=%d, nElem=%d\n', pdeg, Nc, nElem);

            lambda1_list(ii) = run.lambda(1);
            dof_list(ii)     = run.n_dofs_total;
            h_list(ii)        = meta.h;
            continue;
        end

        opts = struct();
        opts.Example            = Example;
        opts.beta               = beta;
        opts.n_gp               = n_gp;
        opts.inner_cheb_n       = inner_cheb_n;
        opts.pw_fft_grid_n      = pw_fft_grid_n;

        opts.primme_tol         = primme_tol;
        opts.primme_maxit       = primme_maxit;
        opts.primme_method      = primme_method;
        opts.primme_reportLevel = primme_reportLevel;
        opts.block_targetShift  = block_targetShift;
        opts.eps_diag           = eps_diag;
        opts.iface_reg          = iface_reg;

        opts.scf_maxit          = scf_maxit;
        opts.scf_pw_grid_m      = scf_pw_grid_m;
        opts.scf_tol_lambda     = scf_tol_lambda;
        opts.scf_mixing         = scf_mixing;
        opts.scf_track_n_eigs   = scf_track_n_eigs;

        opts.save_eigenvectors  = save_eigenvectors;
        opts.save_nurbs         = save_nurbs;
        opts.save_pw_index      = save_pw_index;

        opts.use_pw_cache       = use_pw_cache;
        opts.use_nurbs_cache    = use_nurbs_cache;
        opts.cacheRoot          = cachePwRoot;
        opts.cacheNurbsRoot     = cacheNurbsRoot;

        opts.outDir             = runDir;
        opts.save_matrices      = true;
        opts.save_mat           = true;

        fprintf('\n------------------------------------------------------------\n');
        fprintf('[RUN ] p=%d, Nc=%d, nElem=%d\n', pdeg, Nc, nElem);

        [lambda, ndof, meta] = solve_iga_pw_dg(nElem, t, Nc, 1, opts);

        lambda1_list(ii)  = lambda(1);
        dof_list(ii)      = ndof;
        h_list(ii)        = meta.h;

        fprintf('lambda1 = %.12f\n', lambda1_list(ii));
        fprintf('DOF     = %d\n', ndof);
        fprintf('h       = %.8f\n', h_list(ii));
        fprintf('time    = %.6fs, eigs = %.6fs\n', meta.time_total, meta.time_eigs);
    end

    % Save the h-convergence summary.
    T = table(nElem_list(:), h_list, dof_list, lambda1_list, ...
        'VariableNames', {'nElem','h','dof','lambda1'});

    writetable(T, fullfile(caseDir, 'summary.csv'));

    fprintf('\n[SAVED] %s\n', fullfile(caseDir, 'summary.csv'));
end
end
