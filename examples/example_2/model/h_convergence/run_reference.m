function run_reference(Nc, pdeg, refines)
%Generate the fine-grid reference data.

clc; close all;

activate_example_workflow('h_convergence', ...
    {'nurbs', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core', 'solver'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

Example = 'Example_2';
if nargin < 1 || isempty(Nc)
    Nc = 45;
end
if nargin < 2 || isempty(pdeg)
    pdeg = 3;
end
if nargin < 3 || isempty(refines)
    refines = 8;
end
t = pdeg - 1;
n_eigenvalues = 4;

beta = 20;
n_gp = 10;
inner_cheb_n = 100;
pw_fft_grid_n = 500;

primme_tol         = 1e-12;
primme_maxit       = 2e8;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;
use_pw_cache      = true;

resultRoot = fullfile(pwd, 'result', Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

caseDir = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg));
if ~exist(caseDir, 'dir'), mkdir(caseDir); end

summaryCsv = fullfile(caseDir, 'summary.csv');
n_ref = numel(refines);
Lambda_h = zeros(n_ref, n_eigenvalues);
n_dofs = zeros(n_ref, 1);

fprintf('[INFO] caseDir = %s\n', caseDir);
fprintf('[INFO] Refines to record = %s\n', mat2str(refines.'));

for ii = 1:n_ref
    refine = refines(ii);
    runDir = fullfile(caseDir, sprintf('refine_%02d', refine));
    if ~exist(runDir, 'dir'), mkdir(runDir); end
    runMat = fullfile(runDir, 'run.mat');

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
    opts.eps_diag           = eps_diag;
    opts.iface_reg          = iface_reg;
    opts.save_eigenvectors  = save_eigenvectors;
    opts.save_nurbs         = save_nurbs;
    opts.save_pw_index      = save_pw_index;
    opts.use_pw_cache       = use_pw_cache;
    opts.cacheRoot          = cacheRoot;
    opts.outDir             = runDir;
    opts.save_matrices      = true;
    opts.save_mat           = true;

    fprintf('\n============================================================\n');
    fprintf('[RUN ] Nc=%d, p=%d, refine=%d\n', Nc, pdeg, refine);

    [lambda, ndof] = solve_iga_pw_dg(refine, t, Nc, n_eigenvalues, opts);

    assert(exist(runMat, 'file') == 2, 'Missing reference run file: %s.', runMat);
    assert(numel(lambda) >= n_eigenvalues, 'The solver returned too few eigenvalues.');
    Lambda_h(ii,:) = lambda(1:n_eigenvalues).';
    n_dofs(ii) = ndof;
end

Tnew = table(refines(:), n_dofs, 'VariableNames', {'refine','dof'});
for k = 1:n_eigenvalues
    Tnew.(sprintf('lambda%d', k)) = Lambda_h(:, k);
end

if exist(summaryCsv, 'file') == 2
    Told = readtable(summaryCsv);
    Told(ismember(Told.refine, Tnew.refine), :) = [];
    T = sortrows([Told; Tnew], 'refine');
else
    T = Tnew;
end

writetable(T, summaryCsv);
fprintf('\n[SAVED] %s\n', summaryCsv);
end
