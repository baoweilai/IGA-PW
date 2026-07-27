function run_cutoff_data()
% Generate Example 1 cutoff-convergence data.

clc; close all;

% Set workflow paths, parameters, and cache directories.
activate_example_workflow('cutoff_convergence', ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
rootDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(rootDir);

Example = 'Example_1';

Nc_list         = [2 4 6 8 10 12 14];
t_list          = 0:1;
Refinement_list = [1:4 7];
n_eigenvalues   = 6;

beta = 20;
n_gp = 20;
inner_cheb_n = 200;
pw_fft_grid_n = 256;

primme_tol         = 1e-13;
primme_maxit       = 2e8;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

block_targetShift  = 0.0;

save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;
use_pw_cache      = true;

resultRoot = fullfile(rootDir, 'result', Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cacheRootPW = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRootPW, 'dir'), mkdir(cacheRootPW); end

cacheRootRef = fullfile(resultRoot, 'cache_refine');
if ~exist(cacheRootRef, 'dir'), mkdir(cacheRootRef); end

% Solve each polynomial degree, refinement, and cutoff case.
for t = t_list
    pdeg = 1 + t;

    cachePdegDir = fullfile(cacheRootRef, sprintf('p_%d', pdeg));
    if ~exist(cachePdegDir, 'dir'), mkdir(cachePdegDir); end

    refList = sort(unique(Refinement_list(:)));

    for rr = 1:numel(refList)
        refine = refList(rr);

        refineCacheDir  = fullfile(cachePdegDir, sprintf('refine_%02d', refine));
        if ~exist(refineCacheDir, 'dir'), mkdir(refineCacheDir); end
        refineCacheFile = fullfile(refineCacheDir, 'refine_cache.mat');

        fprintf('\n============================================================\n');
        fprintf('[INFO] Example=%s, p=%d, refine=%d\n', Example, pdeg, refine);

        for Nc = Nc_list
            caseDir = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg));
            if ~exist(caseDir, 'dir'), mkdir(caseDir); end

            runDir = fullfile(caseDir, sprintf('refine_%02d', refine));
            if ~exist(runDir, 'dir'), mkdir(runDir); end
            runMat = fullfile(runDir, 'run.mat');

            if exist(runMat, 'file')
                assert_valid_runmat(runMat, n_eigenvalues);
                fprintf('[LOAD] Nc=%d p=%d refine=%d\n', Nc, pdeg, refine);
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
            opts.eps_diag           = eps_diag;
            opts.iface_reg          = iface_reg;
            opts.block_targetShift  = block_targetShift;
            opts.save_eigenvectors  = save_eigenvectors;
            opts.save_nurbs         = save_nurbs;
            opts.save_pw_index      = save_pw_index;
            opts.use_pw_cache       = use_pw_cache;
            opts.cacheRoot          = cacheRootPW;
            opts.outDir             = runDir;
            opts.refineCacheFile    = refineCacheFile;
            opts.save_matrices      = true;
            opts.save_mat           = true;

            fprintf('\n------------------------------------------------------------\n');
            fprintf('[RUN ] Nc=%d, p=%d, refine=%d\n', Nc, pdeg, refine);

            solve_iga_pw_dg(refine, t, Nc, n_eigenvalues, opts);
        end

        % Collect and save the cutoff results for this refinement.
        summaryCsv = fullfile(refineCacheDir, 'summary.csv');
        LamNc = zeros(numel(Nc_list), n_eigenvalues);
        DofNc = zeros(numel(Nc_list), 1);

        for a = 1:numel(Nc_list)
            Nc = Nc_list(a);
            runMat  = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg), ...
                sprintf('refine_%02d', refine), 'run.mat');
            run = load_run(runMat, n_eigenvalues);
            LamNc(a,:) = run.lambda(1:n_eigenvalues);
            DofNc(a)   = run.n_dofs_total;
        end

        T = table(Nc_list(:), DofNc, 'VariableNames', {'Nc','dof'});
        for k = 1:n_eigenvalues
            T.(sprintf('lambda%d', k)) = LamNc(:, k);
        end
        writetable(T, summaryCsv);

        fprintf('\n[SAVED] %s\n', summaryCsv);
    end
end

end

function run = load_run(runMat, n_eigs)
% Read one completed cutoff case.

S = load(runMat, 'run');
assert(isfield(S, 'run'), 'Missing run structure in %s.', runMat);
run = S.run;
assert(isfield(run, 'lambda') && numel(run.lambda) >= n_eigs, ...
    'Missing eigenvalues in %s.', runMat);
assert(isfield(run, 'n_dofs_total'), 'Missing total DOFs in %s.', runMat);
end

function assert_valid_runmat(runMat, n_eigs)
% Check the cached run file used by cutoff data.

run = load_run(runMat, n_eigs);
assert(isfield(run, 'meta') && isfield(run.meta, 'prec_type'), ...
    'Missing preconditioner metadata in %s.', runMat);
assert(strcmp(run.meta.prec_type, 'InterfaceBlock'), ...
    'Unexpected preconditioner type in %s.', runMat);
end
