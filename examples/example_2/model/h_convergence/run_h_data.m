function run_h_data(refines)
% Generate Example 2 h-convergence data.

clc; close all;

% Set workflow paths, solver parameters, and output directories.
activate_example_workflow('h_convergence', ...
    {'nurbs', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

Example = 'Example_2';

caseSpecs = [ ...
    struct('Nc', 30, 't', 1, 'refines', refines)];
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

% Solve or load each configured refinement case.
for icase = 1:numel(caseSpecs)
    Nc = caseSpecs(icase).Nc;
    t = caseSpecs(icase).t;
    pdeg = 1 + t;

    caseDir = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    summaryCsv = fullfile(caseDir, 'summary.csv');

    refList = caseSpecs(icase).refines(:);
    n_ref = numel(refList);

    Lambda_h = zeros(n_ref, n_eigenvalues);
    n_dofs   = zeros(n_ref, 1);

    fprintf('[INFO] caseDir = %s\n', caseDir);
    fprintf('[INFO] Refines to record = %s\n', mat2str(refList.'));

    for ii = 1:n_ref
        refine = refList(ii);

        runDir = fullfile(caseDir, sprintf('refine_%02d', refine));
        if ~exist(runDir, 'dir'), mkdir(runDir); end
        runMat = fullfile(runDir, 'run.mat');

        if exist(runMat, 'file')
            S = load(runMat, 'run');
            assert(isfield(S, 'run'), 'Missing run structure in %s.', runMat);
            assert(isfield(S.run, 'lambda') && numel(S.run.lambda) >= n_eigenvalues, ...
                'Missing eigenvalues in %s.', runMat);
            assert(isfield(S.run, 'n_dofs_total'), 'Missing total DOFs in %s.', runMat);

            fprintf('[LOAD] Nc=%d p=%d refine=%d\n', Nc, pdeg, refine);
            Lambda_h(ii,:) = S.run.lambda(1:n_eigenvalues);
            n_dofs(ii)     = S.run.n_dofs_total;
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

        assert(numel(lambda) >= n_eigenvalues, 'The solver returned too few eigenvalues.');
        Lambda_h(ii,:) = lambda(1:n_eigenvalues).';
        n_dofs(ii)     = ndof;
    end

    % Merge and save the refinement summary.
    Tnew = table(refList, n_dofs, 'VariableNames', {'refine','dof'});
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

end
