function run_cutoff_data()
%Generate cutoff-convergence data.

clc; close all;

activate_example_workflow('cutoff_convergence', ...
    {'nurbs', 'iga', 'assembly', 'operators', 'core', 'solver'});
exampleDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(exampleDir, 'data');
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

Example = 'Example_3';

fixed_Refinement_list = 6;
fixed_t_list          = 0;
Nc_list               = 2:8;

n_eigenvalues = 1;
beta = 20;
n_gp = 10;
inner_cheb_n = 80;
pw_fft_grid_n = 128;

primme_tol         = 1e-8;
primme_maxit       = 1e8;
primme_method      = 'DEFAULT_MIN_TIME';
primme_reportLevel = 0;

block_targetShift  = 0.0;
eps_diag           = 1e-12;
iface_reg          = 1e-12;

scf_maxit        = 60;
scf_pw_grid_m    = 500;
scf_tol_lambda   = 1e-6;
scf_mixing       = 0.9;
scf_track_n_eigs = 2;

save_eigenvectors = true;
save_nurbs        = true;
save_pw_index     = true;

use_pw_cache    = true;
use_nurbs_cache = true;

resultBase = fullfile(pwd, 'result');
if ~exist(resultBase, 'dir'), mkdir(resultBase); end

resultRoot = fullfile(resultBase, Example);
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

cachePwRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(cachePwRoot, 'dir'), mkdir(cachePwRoot); end

cacheNurbsRoot = fullfile(resultRoot, 'cache_nurbs');
if ~exist(cacheNurbsRoot, 'dir'), mkdir(cacheNurbsRoot); end

for fixed_Refinement = fixed_Refinement_list
    for fixed_t = fixed_t_list

        pdeg = 1 + fixed_t;

        caseDir = fullfile(resultRoot, sprintf('refine_%02d', fixed_Refinement), sprintf('p_%d', pdeg));
        if ~exist(caseDir, 'dir'), mkdir(caseDir); end

        summaryCsv = fullfile(caseDir, 'summary.csv');

        Nc_from_disk = scan_existing_Nc(caseDir);
        Nc_all = unique([Nc_list(:); Nc_from_disk(:)]);
        Nc_all = sort(Nc_all(:));
        n_Nc = numel(Nc_all);

        fprintf('\n============================================================\n');
        fprintf('[INFO] caseDir = %s\n', caseDir);
        fprintf('[INFO] Fixed IGA: refine=%d, p=%d\n', fixed_Refinement, pdeg);
        fprintf('[INFO] Nc to record (merged) = %s\n', mat2str(Nc_all.'));

        Lambda_h = zeros(n_Nc, n_eigenvalues);
        n_dofs   = zeros(n_Nc, 1);

        for ii = 1:n_Nc
            Nc = Nc_all(ii);

            runDir = fullfile(caseDir, sprintf('Nc_%02d', Nc));
            if ~exist(runDir, 'dir'), mkdir(runDir); end
            runMat = fullfile(runDir, 'run.mat');

            if exist(runMat, 'file')
                S = load(runMat, 'run');
                run = S.run;
                fprintf('[LOAD] refine=%d p=%d Nc=%d (from run.mat)\n', fixed_Refinement, pdeg, Nc);

                Lambda_h(ii,:) = run.lambda(1:n_eigenvalues);
                n_dofs(ii)     = run.n_dofs_total;
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
            opts.scf_tol_lambda     = scf_tol_lambda;
            opts.scf_mixing         = scf_mixing;
            opts.scf_pw_grid_m      = scf_pw_grid_m;
            opts.scf_track_n_eigs   = scf_track_n_eigs;

            opts.save_eigenvectors  = save_eigenvectors;
            opts.save_nurbs         = save_nurbs;
            opts.save_pw_index      = save_pw_index;

            opts.use_pw_cache       = use_pw_cache;
            opts.cacheRoot          = cachePwRoot;
            opts.use_nurbs_cache    = use_nurbs_cache;
            opts.cacheNurbsRoot     = cacheNurbsRoot;

            opts.outDir             = runDir;
            opts.save_matrices      = true;
            opts.save_mat           = true;

            fprintf('\n------------------------------------------------------------\n');
            fprintf('[RUN ] refine=%d, p=%d, Nc=%d\n', fixed_Refinement, pdeg, Nc);

            [lambda, ndof, meta] = solve_iga_pw_dg(fixed_Refinement, fixed_t, Nc, n_eigenvalues, opts);

            Lambda_h(ii,:) = lambda(:).';
            n_dofs(ii)     = ndof;

            fprintf('lambda1 = %.12f\n', lambda(1));
            fprintf('DOF = %d\n', ndof);
            fprintf('time_total = %.6fs, time_eigs = %.6fs\n', meta.time_total, meta.time_eigs);
        end

        T = table(Nc_all, n_dofs, 'VariableNames', {'Nc','dof'});

        for k = 1:n_eigenvalues
            T.(sprintf('lambda%d', k)) = Lambda_h(:, k);
        end

        writetable(T, summaryCsv);

        fprintf('\n[SAVED] %s\n', summaryCsv);
    end
end

end

function NcList = scan_existing_Nc(caseDir)
%Compute existing nc.
NcList = [];
dd = dir(fullfile(caseDir, 'Nc_*'));
for i = 1:numel(dd)
    if ~dd(i).isdir, continue; end
    tok = regexp(dd(i).name, '^Nc_(\d+)$', 'tokens', 'once');
    if isempty(tok), continue; end
    Nc = str2double(tok{1});
    runMat = fullfile(caseDir, dd(i).name, 'run.mat');
    if exist(runMat, 'file')
        NcList(end+1,1) = Nc; %#ok<AGROW>
    end
end
NcList = unique(sort(NcList));
end
