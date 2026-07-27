function run_h_data(caseSpecs)
% Generate Example 1 h-convergence data.

clc; close all;

% Set workflow paths, solver parameters, and output directories.
activate_example_workflow('h_convergence', ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

Example = 'Example_1';
n_eigenvalues = 4;

beta = 20;
n_gp = 10;
inner_cheb_n = 200;
pw_fft_grid_n = 500;

primme_tol         = 1e-12;
primme_maxit       = 5e7;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;
pattern_tol        = 1e-12;

save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;
use_pw_cache      = true;

resultRoot = fullfile(pwd, 'result', Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

defaultCacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(defaultCacheRoot, 'dir'), mkdir(defaultCacheRoot); end

% Solve or load each configured refinement case.
for icase = 1:numel(caseSpecs)
    Nc = caseSpecs(icase).Nc;
    t = caseSpecs(icase).t;
    pdeg = 1 + t;
    cacheRoot = defaultCacheRoot;
    if isfield(caseSpecs(icase), 'cacheRoot') && ~isempty(caseSpecs(icase).cacheRoot)
        cacheRoot = caseSpecs(icase).cacheRoot;
    end
    if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

    caseDir = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    summaryCsv = fullfile(caseDir, 'summary.csv');

    refList = unique([caseSpecs(icase).refines(:); scan_existing_refines(caseDir)]);
    refList = sort(refList(:));
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
        opts.pattern_tol        = pattern_tol;
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

    % Save the refinement summary for this case.
    T = table(refList, n_dofs, 'VariableNames', {'refine','dof'});
    for k = 1:n_eigenvalues
        T.(sprintf('lambda%d', k)) = Lambda_h(:, k);
    end

    writetable(T, summaryCsv);
    fprintf('\n[SAVED] %s\n', summaryCsv);
end

end

function refList = scan_existing_refines(caseDir)
% Find completed refinement folders.

refList = [];
dd = dir(fullfile(caseDir, 'refine_*'));

for i = 1:numel(dd)
    if ~dd(i).isdir, continue; end
    tok = regexp(dd(i).name, '^refine_(\d+)$', 'tokens', 'once');
    if isempty(tok), continue; end

    runMat = fullfile(caseDir, dd(i).name, 'run.mat');
    if exist(runMat, 'file')
        refList(end+1,1) = str2double(tok{1}); %#ok<AGROW>
    end
end

refList = unique(sort(refList));
end
