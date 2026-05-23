function H = example_helpers(userCfg)
%Return shared configuration values.
assert(exist('userCfg', 'var') == 1, 'example_helpers requires userCfg.');

workflowDir = fileparts(fileparts(mfilename('fullpath')));
exampleDir = fileparts(fileparts(workflowDir));
setup_paths_local(workflowDir);
cfg = build_cfg_local(fullfile(exampleDir, 'data'), userCfg);

H = struct();
H.cfg = cfg;
H.ensure_reference = @ensure_reference_local;
H.run_case = @run_case_local;
H.case_run_file = @case_run_file_local;
H.run_refinement_case = @run_refinement_case_local;
H.refinement_case_run_file = @refinement_case_run_file_local;

    function ref = ensure_reference_local()
        %Ensure reference.
        ref = run_case_local('reference', cfg.reference.p, cfg.reference.Nc, cfg.reference.Nelement);

        lambda_ref = ref.lambda;
        uh_ref = ref.run.uh;
        rho_ref = ref.run.rhoGrid;
        scf_history = ref.run.scf_history;
        config = cfg;
        p_ref = cfg.reference.p;
        Nc_ref = cfg.reference.Nc;
        Nelement_ref = cfg.reference.Nelement;
        refRunFile = ref.runFile;

        save(cfg.reference_result_file, ...
            'lambda_ref', 'uh_ref', 'rho_ref', 'scf_history', 'config', ...
            'p_ref', 'Nc_ref', 'Nelement_ref', 'refRunFile', '-v7.3');
    end

    function caseOut = run_case_local(studyName, pdeg, Nc, nElem)
        %Run one parameter case.
        runDir = get_case_dir_local(studyName, pdeg, Nc, nElem);
        runFile = fullfile(runDir, 'run.mat');

        if ~exist(runDir, 'dir')
            mkdir(runDir);
        end

        if ~exist(runFile, 'file') || cfg.force
            opts = build_opts_local(runDir, studyName);
            Refinement = struct('mode', 'nelem', 'value', nElem);
            t = max(round(pdeg) - 1, 0);

            solve_scf_3d(Refinement, t, Nc, opts);
        end

        caseOut = load_run_local(runFile);
        caseOut.study = studyName;
        caseOut.p = pdeg;
        caseOut.Nc = Nc;
        caseOut.Nelement = nElem;
    end

    function caseOut = run_refinement_case_local(studyName, pdeg, Nc, refineValue)
        %Run refinement case.
        runDir = get_refinement_case_dir_local(studyName, pdeg, Nc, refineValue);
        runFile = fullfile(runDir, 'run.mat');

        if ~exist(runDir, 'dir')
            mkdir(runDir);
        end

        if ~exist(runFile, 'file') || cfg.force
            opts = build_opts_local(runDir, studyName);
            Refinement = refineValue;
            t = max(round(pdeg) - 1, 0);

            solve_scf_3d(Refinement, t, Nc, opts);
        end

        caseOut = load_run_local(runFile);
        caseOut.study = studyName;
        caseOut.p = pdeg;
        caseOut.K = Nc;
        caseOut.Nc = Nc;
        caseOut.refine = refineValue;
    end

    function runFile = case_run_file_local(studyName, pdeg, Nc, nElem)
        %Build the MAT-file path for one run.
        runFile = fullfile(get_case_dir_local(studyName, pdeg, Nc, nElem), 'run.mat');
    end

    function runFile = refinement_case_run_file_local(studyName, pdeg, Nc, refineValue)
        %Build the MAT-file path for one case.
        runFile = fullfile(get_refinement_case_dir_local(studyName, pdeg, Nc, refineValue), ...
            'run.mat');
    end

    function runDir = get_case_dir_local(studyName, pdeg, Nc, nElem)
        %Return case dir.
        runDir = fullfile(cfg.resultRoot, upper(studyName), ...
            sprintf('K_%d', Nc), sprintf('p_%d', pdeg), ...
            sprintf('nelem_%02d', nElem));
    end

    function runDir = get_refinement_case_dir_local(studyName, pdeg, Nc, refineValue)
        %Return refinement case dir.
        runDir = fullfile(cfg.resultRoot, lower(studyName), sprintf('K_%d', Nc), ...
            sprintf('p_%d', pdeg), sprintf('refine_%02d', refineValue));
    end

    function caseOut = load_run_local(runFile)
        %Load one saved run file.
        S = load(runFile, 'run');
        caseOut = struct();
        caseOut.runFile = runFile;
        caseOut.run = S.run;
        caseOut.lambda = real(S.run.lambda(1));
    end

    function opts = build_opts_local(runDir, studyName)
        %Build solver option values.
        opts = struct();
        opts.Example = cfg.Example;
        opts.runMode = cfg.mode_name;
        opts.reference_lambda = get_reference_lambda_local(studyName);
        opts.reference_energy = cfg.reference_energy;
        opts.beta = cfg.beta;
        opts.n_gp = cfg.n_gp;
        opts.L = cfg.L;
        opts.a = cfg.a;
        opts.nuclear_charge = cfg.nuclear_charge;
        opts.inner_cheb_n = cfg.inner_cheb_n;
        opts.pw_fft_grid_n = cfg.pw_fft_grid_n;
        opts.hartree_grid_n = cfg.hartree_grid_n;
        opts.block_targetShift = 0.0;
        opts.primme_tol = cfg.primme_tol;
        opts.primme_maxit = cfg.primme_maxit;
        opts.primme_method = cfg.primme_method;
        opts.primme_reportLevel = 0;
        opts.scf_maxit = cfg.scf_maxit;
        opts.scf_tol_eig = cfg.scf_tol_eig;
        opts.scf_tol_rho = cfg.scf_tol_rho;
        opts.scf_stopping_rule = cfg.scf_stopping_rule;
        opts.scf_beta = cfg.scf_beta;
        opts.init_guess_mode = cfg.init_guess_mode;
        opts.eps_iface = cfg.eps_iface;
        opts.iface_explicit_gamma_max = cfg.iface_explicit_gamma_max;
        opts.use_direct_interface_gamma = cfg.use_direct_interface_gamma;
        opts.iface_direct_trace_entry_max = cfg.iface_direct_trace_entry_max;
        opts.use_tensor_api = cfg.use_tensor_api;
        opts.use_tensor_face_data = cfg.use_tensor_face_data;
        opts.use_iga_grid_eval_cache = cfg.use_iga_grid_eval_cache;
        opts.use_affine_cube_fast = cfg.use_affine_cube_fast;
        opts.use_exchange_correlation = false;
        opts.use_pw_cache = true;
        opts.use_vext_cache = true;
        opts.save_eigenvectors = cfg.save_eigenvectors;
        opts.save_nurbs = cfg.save_nurbs;
        opts.save_pw_index = cfg.save_pw_index;
        opts.save_density_diagnostics = cfg.save_density_diagnostics;
        opts.cacheRoot = cfg.cacheRoot;
        opts.outDir = runDir;
        opts.smoke = cfg.smoke;
    end

    function lambda_ref = get_reference_lambda_local(studyName)
        %Return reference lambda.
        lambda_ref = [];
        if strcmpi(studyName, 'reference')
            return;
        end
        if isfield(cfg, 'reference_lambda') && ~isempty(cfg.reference_lambda) && isfinite(cfg.reference_lambda)
            lambda_ref = cfg.reference_lambda;
            return;
        end
        if exist(cfg.reference_result_file, 'file')
            S = load(cfg.reference_result_file, 'lambda_ref');
            if isfield(S, 'lambda_ref')
                lambda_ref = S.lambda_ref;
            end
        end
    end
end

function cfg = build_cfg_local(rootDir, userCfg)
%Build the Example 4 configuration struct.
cfg = struct();
cfg.Example = 'Example_5_3D_DG_PW_IGA_0414';
cfg.smoke = userCfg.smoke;
cfg.force = userCfg.force;
cfg.reference_lambda = userCfg.reference_lambda;
cfg.reference_energy = userCfg.reference_energy;
cfg.beta = 20;
cfg.n_gp = 10;
cfg.L = 4;
cfg.a = 0.2;
cfg.nuclear_charge = 2;
cfg.eps_iface = 1e-10;
cfg.use_direct_interface_gamma = true;
cfg.iface_direct_trace_entry_max = 2e7;
cfg.use_tensor_api = true;
cfg.use_tensor_face_data = true;
cfg.use_iga_grid_eval_cache = true;
cfg.use_affine_cube_fast = true;
cfg.primme_method = 'DEFAULT_MIN_MATVECS';
cfg.primme_maxit = 1e8;
cfg.init_guess_mode = 'hybrid_constant';
cfg.scf_stopping_rule = 'lambda_only';
cfg.save_eigenvectors = true;
cfg.save_nurbs = true;
cfg.save_pw_index = true;
cfg.save_density_diagnostics = true;

if cfg.smoke
    cfg.mode_name = 'smoke';
    cfg.reference = struct('p', 2, 'Nc', 4, 'Nelement', 6);
    cfg.pw = struct('fixed_p', 2, 'fixed_Nelement', 6, 'Nc_list', 2:4);
    cfg.iga = struct('fixed_Nc', 4, 'p_list', 1:2, 'Nelement_list', 2:2:4);
    cfg.inner_cheb_n = 24;
    cfg.pw_fft_grid_n = 200;
    cfg.hartree_grid_n = 200;
    cfg.scf_maxit = 20;
    cfg.scf_tol_eig = 1e-6;
    cfg.scf_tol_rho = 1e-4;
    cfg.scf_beta = 0.25;
    cfg.primme_tol = 1e-6;
    cfg.iface_explicit_gamma_max = 12000;
    cfg.state_error_grid_n = 48;
else
    cfg.mode_name = 'formal';
    cfg.reference = struct('p', 2, 'Nc', 45, 'Nelement', 32);
    cfg.pw = struct('fixed_p', 1, 'fixed_Nelement', 12, 'Nc_list', 4:8);
    cfg.iga = struct('fixed_Nc', 20, 'p_list', [1 2], 'Nelement_list', [2, 4, 8, 12]);
    cfg.energy = struct('fixed_p', 2, ...
        'K_list', [2, 4, 8, 12, 16], ...
        'refine_list', [1, 2, 3, 4, 5]);
    cfg.inner_cheb_n = 80;
    cfg.pw_fft_grid_n = 300;
    cfg.hartree_grid_n = 300;
    cfg.scf_maxit = 40;
    cfg.scf_tol_eig = 1e-7;
    cfg.scf_tol_rho = 1e-7;
    cfg.scf_beta = 0.80;
    cfg.scf_stopping_rule = 'lambda_and_rho';
    cfg.primme_tol = 1e-9;
    cfg.iface_explicit_gamma_max = 12000;
    cfg.state_error_grid_n = 100;
end

cfg = merge_struct_local(cfg, userCfg);

cfg.resultRoot = fullfile(rootDir, 'result');
cfg.caseRoot = cfg.resultRoot;
cfg.cacheRoot = fullfile(cfg.resultRoot, 'cache');
cfg.plotRoot = fullfile(cfg.resultRoot, 'plots');
cfg.reference_result_file = fullfile(cfg.resultRoot, 'reference_solution.mat');
cfg.pw_result_file = fullfile(cfg.resultRoot, 'PW', 'pw_convergence.mat');
cfg.iga_result_file = fullfile(cfg.resultRoot, 'IGA', 'iga_convergence.mat');

ensure_dir(cfg.resultRoot);
ensure_dir(cfg.caseRoot);
ensure_dir(cfg.cacheRoot);
ensure_dir(cfg.plotRoot);
end

function setup_paths_local(rootDir)
%Add required model paths.
addpath(fullfile(rootDir, 'config'), '-begin');
addpath(fullfile(rootDir, 'core'), '-begin');
addpath(fullfile(rootDir, 'solver'), '-begin');
addpath(fullfile(rootDir, 'operators'), '-begin');
rehash;
end

function ensure_dir(pathstr)
%Ensure dir.
if ~exist(pathstr, 'dir')
    mkdir(pathstr);
end
end

function out = merge_struct_local(base, override)
%Merge struct.
out = base;
if ~isstruct(override)
    return;
end

fn = fieldnames(override);
for i = 1:numel(fn)
    name = fn{i};
    value = override.(name);
    if isfield(out, name) && isstruct(out.(name)) && isstruct(value)
        out.(name) = merge_struct_local(out.(name), value);
    else
        out.(name) = value;
    end
end
end
