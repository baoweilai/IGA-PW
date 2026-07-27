function run_cutoff_data()
% Generate Example 2 cutoff-convergence data.

clc; close all;

% Set workflow paths, parameters, and cache directories.
activate_example_workflow('cutoff_convergence', ...
    {'nurbs', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
rootDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(rootDir);

Example = 'Example_2';

Nc_list         = [2 4 6 8 10 12 14 16];
t_list          = 0:1;
Refinement_list = 7;
n_eigenvalues   = 4;

beta = 20;
n_gp = 10;
inner_cheb_n = 48;
pw_fft_grid_n = 256;

primme_tol         = 1e-12;
primme_maxit       = 2e8;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

resultRoot = fullfile(rootDir, 'result', Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cacheRootPW  = fullfile(resultRoot, 'cache_pw');
if ~exist(cacheRootPW, 'dir'), mkdir(cacheRootPW); end

cacheRootRef = fullfile(resultRoot, 'cache_refine');
if ~exist(cacheRootRef, 'dir'), mkdir(cacheRootRef); end

% Solve each polynomial degree, refinement, and cutoff case.
for t = t_list
    pdeg = 1 + t;
    refList = sort(unique(Refinement_list(:)));

    for rr = 1:numel(refList)
        refine = refList(rr);

        nNc = numel(Nc_list);
        Lambda_Nc = zeros(nNc, n_eigenvalues);
        n_dofs    = zeros(nNc, 1);

        summaryDir = fullfile(cacheRootRef, sprintf('p_%d', pdeg), sprintf('refine_%02d', refine));
        if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
        summaryCsv = fullfile(summaryDir, 'summary.csv');

        fprintf('\n============================================================\n');
        fprintf('[INFO] Example=%s, p=%d, refine=%d\n', Example, pdeg, refine);

        for ii = 1:nNc
            Nc = Nc_list(ii);
            caseDir = fullfile(resultRoot, sprintf('Nc_%d', Nc), sprintf('p_%d', pdeg));
            if ~exist(caseDir, 'dir'), mkdir(caseDir); end

            runDir = fullfile(caseDir, sprintf('refine_%02d', refine));
            if ~exist(runDir, 'dir'), mkdir(runDir); end
            runMat = fullfile(runDir, 'run.mat');

            if exist(runMat, 'file')
                run = load_run(runMat, n_eigenvalues);
                fprintf('[LOAD] Nc=%d p=%d refine=%d\n', Nc, pdeg, refine);
            else
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
                opts.use_pw_cache       = true;
                opts.cacheRoot          = cacheRootPW;
                opts.outDir             = runDir;
                opts.save_eigenvectors  = true;
                opts.save_nurbs         = true;
                opts.save_pw_index      = true;

                fprintf('\n------------------------------------------------------------\n');
                fprintf('[RUN ] Nc=%d, p=%d, refine=%d\n', Nc, pdeg, refine);

                solve_iga_pw_dg(refine, t, Nc, n_eigenvalues, opts);
                run = load_run(runMat, n_eigenvalues);
            end

            Lambda_Nc(ii,:) = run.lambda(1:n_eigenvalues);
            n_dofs(ii)      = run.n_dofs_total;
        end

        % Save the cutoff summary for this refinement.
        T = table(Nc_list(:), n_dofs, 'VariableNames', {'Nc','dof'});
        for k = 1:n_eigenvalues
            T.(sprintf('lambda%d', k)) = Lambda_Nc(:, k);
        end

        writetable(T, summaryCsv);
        fprintf('[SAVED] %s\n', summaryCsv);
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
