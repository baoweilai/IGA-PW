function results = run_cheb_sensitivity()
% Run the Chebyshev sensitivity study for the first eigenvalue.

clc;
format long g;

% 1. Configure project paths and output directories.
modelDir = fileparts(mfilename('fullpath'));
exampleDir = fileparts(modelDir);
projectDir = fileparts(fileparts(exampleDir));
dataDir = fullfile(exampleDir, 'data');

setup_project_paths_local(projectDir);
add_workflow_paths(fullfile(modelDir, 'h_convergence'), ...
    {'nurbs', 'dg', 'iga', 'assembly', ...
    'operators', 'error_norms', 'core'});

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(dataDir);

resultRoot = fullfile(dataDir, 'result', ...
    'cheb_sensitivity');
cacheRoot = fullfile(resultRoot, 'cache_pw');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end
if ~exist(cacheRoot, 'dir'), mkdir(cacheRoot); end

% 2. Define the reference values and the nine parameter combinations.
% Example 1 reference used by cutoff convergence and full-domain comparison.
lambda1Ref = 4.969971740613;
referenceRefine = 8;
referenceH = 0.4 / (2^referenceRefine);

% The sampled meshes are coarser than the refine-8 reference mesh.
resolutionList = ["coarse_h_K25"; "fine_h_K25"; "finer_h_K35"];
refineList = [4; 6; 7];
KList = [25; 25; 35];
% The n=10 case represents the low-order Chebyshev truncation regime.
innerChebList = [15; 50; 100];

pdeg = 2;
t = pdeg - 1;
nEigenvalues = 1;
pwFftGridN = 500;
innerQuadN = 1000;
ifaceReg = 1e-12;
patternTol = 1e-12;

nCases = numel(resolutionList) * numel(innerChebList);
resolution = strings(nCases, 1);
pColumn = zeros(nCases, 1);
refineColumn = zeros(nCases, 1);
hColumn = zeros(nCases, 1);
referenceHColumn = referenceH * ones(nCases, 1);
KColumn = zeros(nCases, 1);
innerChebColumn = zeros(nCases, 1);
innerQuadColumn = zeros(nCases, 1);
fftColumn = zeros(nCases, 1);
patternTolColumn = patternTol * ones(nCases, 1);
lambda1Column = zeros(nCases, 1);
lambda1RefColumn = lambda1Ref * ones(nCases, 1);
absErrorColumn = zeros(nCases, 1);
dofColumn = zeros(nCases, 1);
timeColumn = zeros(nCases, 1);
sourceColumn = strings(nCases, 1);

% 3. Load each completed case or solve it from the fixed configuration.
row = 0;
for ilevel = 1:numel(resolutionList)
    refine = refineList(ilevel);
    K = KList(ilevel);

    for in = 1:numel(innerChebList)
        innerChebN = innerChebList(in);
        row = row + 1;

        caseTag = sprintf('refine_%02d_K_%03d_cheb_%03d_fft_%03d', ...
            refine, K, innerChebN, pwFftGridN);
        caseFile = fullfile(resultRoot, [caseTag, '.mat']);

        if exist(caseFile, 'file') == 2
            S = load(caseFile, 'caseResult');
            savedCase = S.caseResult;
            assert(savedCase.p == pdeg ...
                && savedCase.refine == refine && savedCase.K == K ...
                && savedCase.inner_cheb_n == innerChebN ...
                && savedCase.inner_quad_n == innerQuadN ...
                && savedCase.pw_fft_grid_n == pwFftGridN ...
                && savedCase.iface_reg == ifaceReg ...
                && savedCase.pattern_tol == patternTol, ...
                'Saved case metadata do not match %s.', caseFile);
            assert(isfield(savedCase, 'lambda') ...
                && numel(savedCase.lambda) >= nEigenvalues, ...
                'Saved case does not contain the first eigenvalue.');
            lambda1Current = savedCase.lambda(1);
            caseResult = build_case_result_local( ...
                resolutionList(ilevel), pdeg, refine, savedCase.h, K, ...
                innerChebN, innerQuadN, pwFftGridN, ifaceReg, patternTol, ...
                lambda1Current, lambda1Ref, savedCase.n_dofs_total, ...
                savedCase.time_total_s);
            save(caseFile, 'caseResult');
            source = "loaded";
        else
            opts = build_solver_options_local( ...
                cacheRoot, innerChebN, pwFftGridN, ifaceReg, patternTol);
            fprintf('\n============================================================\n');
            fprintf('[CASE %d/%d] %s, p=%d, refine=%d, K=%d, n=%d, FFTN=%d\n', ...
                row, nCases, resolutionList(ilevel), pdeg, refine, K, ...
                innerChebN, pwFftGridN);

            [lambda, nDofsTotal, meta] = solve_iga_pw_dg( ...
                refine, t, K, nEigenvalues, opts);
            assert(numel(lambda) >= nEigenvalues, ...
                'The solver did not return the first eigenvalue.');
            lambda = lambda(:);
            lambda1Current = lambda(1);
            assert(isfinite(lambda1Current) ...
                && lambda1Current > 0 ...
                && abs(lambda1Current - lambda1Ref) < 0.1, ...
                'Invalid first eigenvalue for refine=%d, K=%d, n=%d.', ...
                refine, K, innerChebN);

            caseResult = build_case_result_local( ...
                resolutionList(ilevel), pdeg, refine, meta.h, K, ...
                innerChebN, innerQuadN, pwFftGridN, ifaceReg, patternTol, ...
                lambda1Current, lambda1Ref, nDofsTotal, meta.time_total);
            save(caseFile, 'caseResult');
            source = "computed";
        end

        assert(isfinite(caseResult.lambda1) ...
            && caseResult.lambda1 > 0 ...
            && abs(caseResult.lambda1 - lambda1Ref) < 0.1, ...
            'Saved or computed first eigenvalue is invalid.');

        resolution(row) = resolutionList(ilevel);
        pColumn(row) = pdeg;
        refineColumn(row) = refine;
        hColumn(row) = caseResult.h;
        KColumn(row) = K;
        innerChebColumn(row) = innerChebN;
        innerQuadColumn(row) = innerQuadN;
        fftColumn(row) = pwFftGridN;
        lambda1Column(row) = caseResult.lambda1;
        absErrorColumn(row) = caseResult.abs_error;
        dofColumn(row) = caseResult.n_dofs_total;
        timeColumn(row) = caseResult.time_total_s;
        sourceColumn(row) = source;

        results = build_results_table_local(row, resolution, pColumn, ...
            refineColumn, hColumn, referenceHColumn, KColumn, innerChebColumn, ...
            innerQuadColumn, fftColumn, patternTolColumn, lambda1Column, ...
            lambda1RefColumn, absErrorColumn, dofColumn, timeColumn, ...
            sourceColumn);
    end
end

assert(height(results) == 9, 'The final table must contain exactly nine rows.');

% 4. Compute the truncation differences and trend diagnostics.
delta1ToNmax = zeros(height(results), 1);
for ilevel = 1:numel(resolutionList)
    levelRows = results.resolution == resolutionList(ilevel);
    nmaxRow = levelRows & results.inner_cheb_n == innerChebList(end);
    lambda1Nmax = results.lambda1(nmaxRow);
    delta1ToNmax(levelRows) = ...
        abs(results.lambda1(levelRows) - lambda1Nmax);
end
results.delta_lambda1_to_nmax = delta1ToNmax;
results = movevars(results, 'delta_lambda1_to_nmax', 'After', 'abs_error');

fineSmallKSmallN = results(results.resolution == "fine_h_K25" ...
    & results.inner_cheb_n == innerChebList(1), :);
fineLargeKSmallN = results(results.resolution == "finer_h_K35" ...
    & results.inner_cheb_n == innerChebList(1), :);
coarseLargeN = results(results.resolution == "coarse_h_K25" ...
    & results.inner_cheb_n == innerChebList(3), :);
fineSmallKLargeN = results(results.resolution == "fine_h_K25" ...
    & results.inner_cheb_n == innerChebList(3), :);
fineLargeKLargeN = results(results.resolution == "finer_h_K35" ...
    & results.inner_cheb_n == innerChebList(3), :);
fineMedium = results(results.resolution == "finer_h_K35" ...
    & results.inner_cheb_n == innerChebList(2), :);
fineLarge = results(results.resolution == "finer_h_K35" ...
    & results.inner_cheb_n == innerChebList(3), :);

diagnostics = struct();
diagnostics.small_n_fine_resolution_factor = max( ...
    fineSmallKSmallN.abs_error, fineLargeKSmallN.abs_error) / min( ...
    fineSmallKSmallN.abs_error, fineLargeKSmallN.abs_error);
diagnostics.small_n_truncation_delta_factor = max( ...
    fineSmallKSmallN.delta_lambda1_to_nmax, ...
    fineLargeKSmallN.delta_lambda1_to_nmax) / min( ...
    fineSmallKSmallN.delta_lambda1_to_nmax, ...
    fineLargeKSmallN.delta_lambda1_to_nmax);
diagnostics.large_n_h_reduction = ...
    coarseLargeN.abs_error / fineSmallKLargeN.abs_error;
diagnostics.lambda1_joint_refinement_reduction = ...
    fineSmallKLargeN.abs_error / fineLargeKLargeN.abs_error;
diagnostics.fine_large_n_lambda_gap = abs(fineMedium.lambda1 - fineLarge.lambda1);
diagnostics.fine_large_n_gap_over_error = ...
    diagnostics.fine_large_n_lambda_gap / fineLarge.abs_error;
diagnostics.fine_largeK_reference_error_reduction = ...
    fineLargeKSmallN.abs_error / fineMedium.abs_error;
diagnostics.fine_largeK_cheb_delta_reduction = ...
    fineLargeKSmallN.delta_lambda1_to_nmax / ...
    fineMedium.delta_lambda1_to_nmax;

assert(all(results.pw_fft_grid_n == pwFftGridN), ...
    'FFTN must remain fixed in all nine cases.');
assert(all(results.pattern_tol == patternTol), ...
    'The preconditioner trace threshold must remain fixed at 1e-12.');
assert(all(results.h > results.reference_h), ...
    'Every test mesh must be strictly coarser than the refine-8 reference mesh.');
assert(isequal(unique(results.inner_cheb_n), innerChebList), ...
    'The final table must contain exactly the three requested Chebyshev sizes.');
assert(all(isfinite(results.lambda1)) && all(isfinite(results.abs_error)), ...
    'All first-eigenvalue results and absolute errors must be finite.');
for ilevel = 1:numel(resolutionList)
    assert(nnz(results.resolution == resolutionList(ilevel)) == 3, ...
        'Each h/K level must contain exactly three Chebyshev sizes.');
end
% Record the convergence trends after validating the result table.
trendChecks = struct();
trendChecks.lambda1_small_n_plateau = ...
    diagnostics.small_n_truncation_delta_factor <= 1.2;
trendChecks.lambda1_h_reduction = diagnostics.large_n_h_reduction >= 5;
trendChecks.lambda1_reference_error_drop = ...
    diagnostics.fine_largeK_reference_error_reduction >= 1.2;
trendChecks.lambda1_cheb_error_drop = ...
    diagnostics.fine_largeK_cheb_delta_reduction >= 1e3;
trendChecks.lambda1_large_n_insensitivity = ...
    diagnostics.fine_large_n_gap_over_error <= 1e-3;
diagnostics.trend_checks = trendChecks;
diagnostics.all_trend_checks_pass = all(structfun(@(x) x, trendChecks));

% 5. Save the result table, configuration, and diagnostic values.
save(fullfile(resultRoot, 'results.mat'), ...
    'results', 'diagnostics', 'resolutionList', 'refineList', 'KList', ...
    'innerChebList', 'pwFftGridN', 'innerQuadN', 'ifaceReg', ...
    'patternTol', 'nEigenvalues', 'lambda1Ref', 'referenceRefine', ...
    'referenceH');

fprintf('\n[SAVED] %s\n', fullfile(resultRoot, 'results.mat'));
fprintf('[CHECK] small-n fine-level factor     = %.6g\n', ...
    diagnostics.small_n_fine_resolution_factor);
fprintf('[CHECK] small-n truncation factor     = %.6g\n', ...
    diagnostics.small_n_truncation_delta_factor);
fprintf('[CHECK] large-n h reduction           = %.6g\n', ...
    diagnostics.large_n_h_reduction);
fprintf('[CHECK] lambda1 joint-refine reduction = %.6g\n', ...
    diagnostics.lambda1_joint_refinement_reduction);
fprintf('[CHECK] fine large-K ref-error drop  = %.6g\n', ...
    diagnostics.fine_largeK_reference_error_reduction);
fprintf('[CHECK] fine large-K Cheb reduction  = %.6g\n', ...
    diagnostics.fine_largeK_cheb_delta_reduction);
fprintf('[CHECK] fine medium/large lambda gap  = %.6e\n', ...
    diagnostics.fine_large_n_lambda_gap);
fprintf('[CHECK] gap / fine large-n error      = %.6g\n', ...
    diagnostics.fine_large_n_gap_over_error);
trendNames = fieldnames(trendChecks);
for itrend = 1:numel(trendNames)
    trendName = trendNames{itrend};
    if trendChecks.(trendName)
        trendStatus = 'PASS';
    else
        trendStatus = 'FAIL';
    end
    fprintf('[TREND %s] %s\n', trendStatus, trendName);
end

disp(results(:, {'resolution', 'h', 'K', 'inner_cheb_n', ...
    'pw_fft_grid_n', 'lambda1', 'abs_error', ...
    'delta_lambda1_to_nmax'}));
end

function caseResult = build_case_result_local( ...
    resolution, pdeg, refine, h, K, innerChebN, innerQuadN, ...
    pwFftGridN, ifaceReg, patternTol, lambda1, lambda1Ref, ...
    nDofsTotal, timeTotal)
% Build one result containing only the first eigenvalue.
caseResult = struct();
caseResult.resolution = resolution;
caseResult.p = pdeg;
caseResult.refine = refine;
caseResult.h = h;
caseResult.K = K;
caseResult.inner_cheb_n = innerChebN;
caseResult.inner_quad_n = innerQuadN;
caseResult.pw_fft_grid_n = pwFftGridN;
caseResult.iface_reg = ifaceReg;
caseResult.pattern_tol = patternTol;
caseResult.n_eigenvalues = 1;
caseResult.lambda = lambda1;
caseResult.lambda1 = lambda1;
caseResult.lambda1_ref = lambda1Ref;
caseResult.abs_error = abs(lambda1 - lambda1Ref);
caseResult.n_dofs_total = nDofsTotal;
caseResult.time_total_s = timeTotal;
end

function opts = build_solver_options_local( ...
    cacheRoot, innerChebN, pwFftGridN, ifaceReg, patternTol)
% Build the fixed solver options for one sensitivity case.
opts = struct();
opts.Example = 'Example_1';
opts.beta = 20;
opts.n_gp = 10;
opts.inner_cheb_n = innerChebN;
opts.pw_fft_grid_n = pwFftGridN;
opts.primme_tol = 1e-12;
opts.primme_maxit = 5e7;
opts.primme_method = 'DEFAULT_MIN_TIME';
opts.primme_reportLevel = 0;
opts.eps_diag = 1e-12;
opts.iface_reg = ifaceReg;
opts.pattern_tol = patternTol;
opts.use_pw_cache = true;
opts.cacheRoot = cacheRoot;
end

function results = build_results_table_local(nRows, resolution, pColumn, ...
    refineColumn, hColumn, referenceHColumn, KColumn, innerChebColumn, ...
    innerQuadColumn, ...
    fftColumn, patternTolColumn, lambda1Column, lambda1RefColumn, ...
    absErrorColumn, dofColumn, timeColumn, sourceColumn)
% Build the completed part of the nine-row results table.
idx = (1:nRows).';
results = table(resolution(idx), pColumn(idx), refineColumn(idx), ...
    hColumn(idx), referenceHColumn(idx), KColumn(idx), innerChebColumn(idx), ...
    innerQuadColumn(idx), fftColumn(idx), patternTolColumn(idx), ...
    lambda1Column(idx), lambda1RefColumn(idx), absErrorColumn(idx), ...
    dofColumn(idx), timeColumn(idx), sourceColumn(idx), ...
    'VariableNames', {'resolution', 'p', 'refine', 'h', 'reference_h', 'K', ...
    'inner_cheb_n', 'inner_quad_n', 'pw_fft_grid_n', 'pattern_tol', ...
    'lambda1', 'lambda1_ref', 'abs_error', 'n_dofs_total', ...
    'time_total_s', 'source'});
end

function setup_project_paths_local(projectDir)
% Add project dependencies needed by the standalone model driver.
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
    addpath(requiredDirs{k}, '-begin');
end
end
