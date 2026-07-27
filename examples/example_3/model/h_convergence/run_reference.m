function run_reference()
% Build the fine-grid solution used by the DG error curve.

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

dataDir = fullfile(exampleDir, 'data');
startDir = pwd;
cleanupObj = onCleanup(@() cd(startDir));
cd(dataDir);

add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
    {'nurbs', 'iga', 'assembly', 'operators', 'core'});

% Set the fine-grid reference parameters.
Example = 'Example_3';
Refinement = 7;
t = 1;
Nc = 40;
n_eigenvalues = 1;

beta = 100;
n_gp = 30;
inner_cheb_n = 150;
pw_fft_grid_n = 550;

primme_tol = 1e-9;
primme_maxit = 1e8;
primme_method = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
block_targetShift = 0.0;
eps_diag = 1e-12;
iface_reg = 1e-12;

scf_maxit = 60;
scf_pw_grid_m = 600;
scf_tol_lambda = 1e-7;
scf_mixing = 0.9;
scf_track_n_eigs = 1;

save_eigenvectors = true;
save_nurbs = true;
save_pw_index = true;

use_pw_cache = true;
use_nurbs_cache = true;

resultRoot = fullfile(pwd, 'result', Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cachePwRoot = fullfile(resultRoot, 'cache_pw_reference');
if ~exist(cachePwRoot, 'dir'), mkdir(cachePwRoot); end

cacheNurbsRoot = fullfile(resultRoot, 'cache_nurbs_reference');
if ~exist(cacheNurbsRoot, 'dir'), mkdir(cacheNurbsRoot); end

pdeg = 1 + t;
runDir = fullfile(resultRoot, sprintf('refine_%02d', Refinement), ...
    sprintf('p_%d', pdeg), sprintf('Nc_%02d', Nc));
if ~exist(runDir, 'dir'), mkdir(runDir); end

% Load a complete reference run when it is available.
runMat = fullfile(runDir, 'run.mat');
if exist(runMat, 'file')
    S = load(runMat, 'run');
    assert(isfield(S, 'run'), 'Reference file does not contain run: %s', runMat);
    assert(isfield(S.run, 'uh'), 'Reference file does not contain eigenvectors: %s', runMat);
    assert(isfield(S.run, 'nurbs_refine'), 'Reference file does not contain NURBS data: %s', runMat);
    assert(isfield(S.run, 'k_pw'), 'Reference file does not contain plane-wave indices: %s', runMat);
    return;
end

% Solve and save the reference case.
opts = struct();
opts.Example = Example;
opts.beta = beta;
opts.n_gp = n_gp;
opts.inner_cheb_n = inner_cheb_n;
opts.pw_fft_grid_n = pw_fft_grid_n;

opts.primme_tol = primme_tol;
opts.primme_maxit = primme_maxit;
opts.primme_method = primme_method;
opts.primme_reportLevel = primme_reportLevel;
opts.block_targetShift = block_targetShift;
opts.eps_diag = eps_diag;
opts.iface_reg = iface_reg;

opts.scf_maxit = scf_maxit;
opts.scf_pw_grid_m = scf_pw_grid_m;
opts.scf_tol_lambda = scf_tol_lambda;
opts.scf_mixing = scf_mixing;
opts.scf_track_n_eigs = scf_track_n_eigs;

opts.save_eigenvectors = save_eigenvectors;
opts.save_nurbs = save_nurbs;
opts.save_pw_index = save_pw_index;

opts.use_pw_cache = use_pw_cache;
opts.use_nurbs_cache = use_nurbs_cache;
opts.cacheRoot = cachePwRoot;
opts.cacheNurbsRoot = cacheNurbsRoot;

opts.outDir = runDir;
opts.save_matrices = true;
opts.save_mat = true;

solve_iga_pw_dg(Refinement, t, Nc, n_eigenvalues, opts);
end
