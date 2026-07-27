function out = run_hartree_comparison(action)
% Run or validate the retained Hartree comparison cases.

arguments
    action (1,1) string = "run"
end

% Load the workflow configuration and result directory.
activate_example_workflow('hartree_comparison', ...
    {'config', 'core', 'operators', 'solver'});
cfg = default_config();
ensure_directory_local(cfg.resultDir);

% Run the cases or validate the saved results.
switch lower(action)
    case "run"
        out = run_cases_local(cfg);
    case "check"
        out = check_results_local(cfg);
    otherwise
        error('Unknown hartree_comparison action: %s', action);
end
end

function out = run_cases_local(cfg)
% Run four IGA-PW cases and four DG-PW cases.

% Reset workflow outputs and load the shared reference values.
reset_outputs_local(cfg);
[reference, cfg] = load_reference_local(cfg);
save_reference_local(cfg, reference);

igaCsv = fullfile(cfg.resultDir, 'iga.csv');
igaMat = fullfile(cfg.resultDir, 'iga.mat');
dgCsv = fullfile(cfg.resultDir, 'dg.csv');
dgMat = fullfile(cfg.resultDir, 'dg.mat');
comparisonCsv = fullfile(cfg.resultDir, 'comparison.csv');
comparisonMat = fullfile(cfg.resultDir, 'comparison.mat');
logFile = fullfile(cfg.resultDir, 'run.log');

% Record the retained case definitions and run start.
diary(logFile);
fprintf('[comparison] start=%s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('[comparison] IGA cases=%s\n', mat2str(cfg.igaCases));
fprintf('[comparison] DG cases=%s\n', mat2str(cfg.dgCases));

% Solve and store the four IGA-PW cases.
igaResults = empty_iga_results_table_local();
for caseIndex = 1:size(cfg.igaCases, 1)
    parameters = cfg.igaCases(caseIndex, :);
    K = parameters(1);
    p = parameters(2);
    refine = parameters(3);
    fprintf('[comparison] IGA START %d (%d,%d,%d)\n', ...
        caseIndex, K, p, refine);
    point = solve_iga(cfg, K, p, refine);
    validate_iga_point_local(point, cfg);
    save(fullfile(cfg.resultDir, sprintf('iga_%d.mat', caseIndex)), ...
        'point', 'reference', '-v7.3');
    igaResults = [igaResults; iga_row_local(point, cfg)]; %#ok<AGROW>
    writetable(igaResults, igaCsv);
    save(igaMat, 'igaResults', 'reference', 'cfg');
    fprintf(['[comparison] IGA DONE %d lambda=%.16e energy=%.16e ' ...
        'SCF=%d totalTime=%.6f s\n'], caseIndex, point.lambda, ...
        point.energy, point.scfIters, point.totalTime);
    clear point
end

% Solve and store the four DG-PW cases.
dgResults = empty_dg_results_table_local();
for caseIndex = 1:size(cfg.dgCases, 1)
    parameters = cfg.dgCases(caseIndex, :);
    K = parameters(1);
    n_r = parameters(2);
    L_m = parameters(3);
    assert(n_r == L_m, 'DG-PW uses equal n_r and L_m.');
    fprintf('[comparison] DG START %d (%d,%d,%d)\n', ...
        caseIndex, K, n_r, L_m);
    point = solve_dg(cfg, K, n_r);
    validate_dg_point_local(point, cfg);
    fileStem = sprintf('dg_%d', caseIndex);
    writetable(point.stageTimes, fullfile(cfg.resultDir, ...
        [fileStem, '_times.csv']));
    save(fullfile(cfg.resultDir, [fileStem, '.mat']), ...
        'point', 'reference', '-v7.3');
    dgResults = [dgResults; dg_row_local(point, cfg)]; %#ok<AGROW>
    writetable(dgResults, dgCsv);
    save(dgMat, 'dgResults', 'reference', 'cfg');
    fprintf(['[comparison] DG DONE %d lambda=%.16e energy=%.16e ' ...
        'SCF=%d totalTime=%.6f s\n'], caseIndex, point.lambda, ...
        point.energy, point.scfIters, point.totalTime);
    clear point
end

% Combine the two methods and package the workflow outputs.
assert(height(igaResults) == 4, 'The IGA-PW table must contain four rows.');
assert(height(dgResults) == 4, 'The DG-PW table must contain four rows.');
comparison = build_comparison_local(igaResults, dgResults);
writetable(comparison, comparisonCsv);
save(comparisonMat, 'comparison', 'igaResults', 'dgResults', ...
    'reference', 'cfg');
fprintf('[comparison] finish=%s rows=%d\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    height(comparison));
diary off;

out = struct('cfg', cfg, 'reference', reference, ...
    'igaResults', igaResults, 'dgResults', dgResults, ...
    'comparison', comparison, 'csvFile', comparisonCsv, ...
    'matFile', comparisonMat, 'logFile', logFile);
end

function validate_iga_point_local(point, cfg)
% Validate one IGA-PW result.

assert(point.scfConverged, 'The IGA-PW case did not converge.');
assert(point.finalLambdaChange < cfg.scfTolEig ...
    && point.finalDensityChange < cfg.scfTolRho, ...
    'The IGA-PW case did not meet both SCF tolerances.');
assert(point.nDof == point.nPw + point.nLocal, ...
    'The IGA-PW case has inconsistent dimensions.');
end

function validate_dg_point_local(point, cfg)
% Validate one DG-PW result.

assert(point.scfConverged, 'The DG-PW case did not converge.');
assert(point.finalLambdaChange <= cfg.dgScfTolEig ...
    && point.finalDensityChange <= cfg.dgScfTolRho, ...
    'The DG-PW case did not meet both SCF tolerances.');
assert(point.nDof == point.nPw + point.nLocal, ...
    'The DG-PW case has inconsistent dimensions.');
end

function comparison = build_comparison_local(igaResults, dgResults)
% Build the combined eight-row result table.

nIga = height(igaResults);
nDg = height(dgResults);
method = [igaResults.method; dgResults.method];
parameters = strings(nIga + nDg, 1);
sourceFile = strings(nIga + nDg, 1);

for caseIndex = 1:nIga
    parameters(caseIndex) = sprintf('(%d,%d,%d)', ...
        igaResults.K(caseIndex), igaResults.p(caseIndex), ...
        igaResults.refine(caseIndex));
    sourceFile(caseIndex) = sprintf('iga_%d.mat', caseIndex);
end
for caseIndex = 1:nDg
    row = nIga + caseIndex;
    parameters(row) = sprintf('(%d,%d,%d)', ...
        dgResults.K(caseIndex), dgResults.n_r(caseIndex), ...
        dgResults.L_m(caseIndex));
    sourceFile(row) = sprintf('dg_%d.mat', caseIndex);
end

nDof = [igaResults.nDof; dgResults.nDof];
energyError = [igaResults.energyError; dgResults.energyError];
totalTime = [igaResults.totalTime; dgResults.totalTime];
comparison = table(method, parameters, nDof, energyError, totalTime, ...
    sourceFile);
end

function out = check_results_local(cfg)
% Validate the retained code paths and results.

% Check the required workflow directories and files.
requiredDirs = {
    cfg.workflowDir
    cfg.configDir
    cfg.coreDir
    cfg.operatorsDir
    cfg.commonOperators
    cfg.solverDir
    cfg.referenceOperators
    };
for index = 1:numel(requiredDirs)
    assert(isfolder(requiredDirs{index}), ...
        'Missing directory: %s', requiredDirs{index});
end

requiredFiles = {
    cfg.referenceFile
    cfg.dgScfFile
    fullfile(cfg.coreDir, 'solve_iga.m')
    fullfile(cfg.coreDir, 'solve_dg.m')
    fullfile(cfg.operatorsDir, 'build_ewald_potential.m')
    fullfile(cfg.solverDir, 'rfp_chol_mex.c')
    fullfile(cfg.resultDir, 'reference.csv')
    fullfile(cfg.resultDir, 'reference.mat')
    fullfile(cfg.resultDir, 'iga.csv')
    fullfile(cfg.resultDir, 'iga.mat')
    fullfile(cfg.resultDir, 'dg.csv')
    fullfile(cfg.resultDir, 'dg.mat')
    fullfile(cfg.resultDir, 'comparison.csv')
    fullfile(cfg.resultDir, 'comparison.mat')
    fullfile(cfg.resultDir, 'run.log')
    };
for caseIndex = 1:4
    requiredFiles{end + 1} = fullfile(cfg.resultDir, ...
        sprintf('iga_%d.mat', caseIndex)); %#ok<AGROW>
    requiredFiles{end + 1} = fullfile(cfg.resultDir, ...
        sprintf('dg_%d.mat', caseIndex)); %#ok<AGROW>
    requiredFiles{end + 1} = fullfile(cfg.resultDir, ...
        sprintf('dg_%d_times.csv', caseIndex)); %#ok<AGROW>
end
for index = 1:numel(requiredFiles)
    assert(isfile(requiredFiles{index}), ...
        'Missing file: %s', requiredFiles{index});
end
stored = dir(cfg.resultDir);
stored = string({stored(~[stored.isdir]).name}).';
allowed = [string(run_output_names_local()); "table.csv"];
assert(isempty(setdiff(stored, allowed)), ...
    'The result directory contains files outside the eight retained cases.');

% Check that MATLAB resolves the retained workflow implementations.
assert(strcmpi(which('solve_iga'), fullfile(cfg.coreDir, 'solve_iga.m')), ...
    'solve_iga does not resolve to this workflow.');
assert(strcmpi(which('solve_dg'), fullfile(cfg.coreDir, 'solve_dg.m')), ...
    'solve_dg does not resolve to this workflow.');
assert(strcmpi(which('dg_scf'), cfg.dgScfFile), ...
    'dg_scf does not resolve to this workflow.');
assert(strcmpi(which('build_ewald_potential'), ...
    fullfile(cfg.operatorsDir, 'build_ewald_potential.m')), ...
    'build_ewald_potential does not resolve to this workflow.');

% Load the tabular and MAT results for cross-checking.
igaResults = readtable(fullfile(cfg.resultDir, 'iga.csv'), ...
    'TextType', 'string');
dgResults = readtable(fullfile(cfg.resultDir, 'dg.csv'), ...
    'TextType', 'string');
comparison = readtable(fullfile(cfg.resultDir, 'comparison.csv'), ...
    'TextType', 'string');
igaStored = load(fullfile(cfg.resultDir, 'iga.mat'), 'igaResults');
dgStored = load(fullfile(cfg.resultDir, 'dg.mat'), 'dgResults');
comparisonStored = load(fullfile(cfg.resultDir, 'comparison.mat'), ...
    'comparison');

igaTemplate = empty_iga_results_table_local();
dgTemplate = empty_dg_results_table_local();

% Validate schemas, case parameters, and convergence data.
assert(isequal(igaTemplate.Properties.VariableNames, ...
    igaResults.Properties.VariableNames), ...
    'The IGA-PW result schema is inconsistent.');
assert(isequal(dgTemplate.Properties.VariableNames, ...
    dgResults.Properties.VariableNames), ...
    'The DG-PW result schema is inconsistent.');
assert(isequal([igaResults.K, igaResults.p, igaResults.refine], ...
    cfg.igaCases), 'The IGA-PW cases are inconsistent.');
assert(isequal([dgResults.K, dgResults.n_r, dgResults.L_m], ...
    cfg.dgCases), 'The DG-PW cases are inconsistent.');
assert(isequal([igaStored.igaResults.K, igaStored.igaResults.p, ...
    igaStored.igaResults.refine], cfg.igaCases), ...
    'The IGA-PW MAT cases are inconsistent.');
assert(isequal([dgStored.dgResults.K, dgStored.dgResults.n_r, ...
    dgStored.dgResults.L_m], cfg.dgCases), ...
    'The DG-PW MAT cases are inconsistent.');
assert(all(igaResults.scfConverged) && all(dgResults.scfConverged), ...
    'Every retained case must be converged.');
assert(all(igaResults.nDof == igaResults.nPw + igaResults.nLocal), ...
    'The IGA-PW DOF values are inconsistent.');
assert(all(dgResults.nDof == dgResults.nPw + dgResults.nLocal), ...
    'The DG-PW DOF values are inconsistent.');
assert(isscalar(unique([igaResults.energyRef; dgResults.energyRef])), ...
    'Both methods must use one reference energy.');

expected = build_comparison_local(igaResults, dgResults);
assert(height(comparison) == 8, ...
    'The comparison table must contain eight rows.');
assert(isequal(comparison.Properties.VariableNames, ...
    expected.Properties.VariableNames), ...
    'The comparison result schema is inconsistent.');
assert(isequal(string(comparison.method), expected.method), ...
    'The comparison methods are inconsistent.');
assert(isequal(string(comparison.parameters), expected.parameters), ...
    'The comparison parameters are inconsistent.');
assert(isequal(comparison.nDof, expected.nDof), ...
    'The comparison DOF values are inconsistent.');
assert(all(abs(comparison.energyError - expected.energyError) ...
    <= 1e-12 * max(1, abs(expected.energyError))), ...
    'The comparison energy errors are inconsistent.');
assert(all(abs(comparison.totalTime - expected.totalTime) ...
    <= 1e-12 * max(1, abs(expected.totalTime))), ...
    'The comparison times are inconsistent.');
assert(isequal(string(comparison.sourceFile), expected.sourceFile), ...
    'The comparison source files are inconsistent.');
assert(isequal(string(comparisonStored.comparison.parameters), ...
    expected.parameters), 'The comparison MAT parameters are inconsistent.');

% Cross-check each per-case MAT file against the summary tables.
for caseIndex = 1:4
    igaPoint = load(fullfile(cfg.resultDir, ...
        sprintf('iga_%d.mat', caseIndex)), 'point');
    dgPoint = load(fullfile(cfg.resultDir, ...
        sprintf('dg_%d.mat', caseIndex)), 'point');
    assert(isequal([igaPoint.point.K, igaPoint.point.p, ...
        igaPoint.point.refine], cfg.igaCases(caseIndex, :)), ...
        'An IGA-PW case file has inconsistent parameters.');
    assert(isequal([dgPoint.point.K, dgPoint.point.q, dgPoint.point.q], ...
        cfg.dgCases(caseIndex, :)), ...
        'A DG-PW case file has inconsistent parameters.');
    assert(igaPoint.point.nDof == igaResults.nDof(caseIndex), ...
        'An IGA-PW case file has inconsistent DOF.');
    assert(dgPoint.point.nDof == dgResults.nDof(caseIndex), ...
        'A DG-PW case file has inconsistent DOF.');
end

% Package the validation result.
out = struct('status', "PASS", 'cfg', cfg, ...
    'igaRows', height(igaResults), 'dgRows', height(dgResults), ...
    'comparisonRows', height(comparison), ...
    'energyRef', igaResults.energyRef(1));
end

function [reference, cfg] = load_reference_local(cfg)
% Load the shared reference eigenvalue and energy.

assert(isfile(cfg.referenceFile), ...
    'Reference run file does not exist: %s', cfg.referenceFile);
lambdaRef = real(double(h5read(cfg.referenceFile, '/run/lambda')));
energyRef = real(double(h5read( ...
    cfg.referenceFile, '/run/meta/energy_total')));
lambdaRef = lambdaRef(1);
energyRef = energyRef(1);
assert(isfinite(lambdaRef) && isfinite(energyRef), ...
    'The saved reference values must be finite.');

referenceInnerChebN = double(h5read( ...
    cfg.referenceFile, '/run/meta/inner_cheb_n'));
referencePwFftGridN = double(h5read( ...
    cfg.referenceFile, '/run/meta/pw_fft_grid_n'));
referenceHartreeGridN = double(h5read( ...
    cfg.referenceFile, '/run/meta/hartree_grid_n'));
referenceScfTolEig = double(h5read( ...
    cfg.referenceFile, '/run/meta/scf_tol_eig'));
referenceScfTolRho = double(h5read( ...
    cfg.referenceFile, '/run/meta/scf_tol_rho'));
referencePrimmeTol = double(h5read( ...
    cfg.referenceFile, '/run/meta/primme_tol'));
assert(cfg.nGauss == 10, ...
    'IGA-PW direct Gauss quadrature must use 10 points.');
assert(cfg.chebDegree == referenceInnerChebN ...
    && cfg.fftGridN == referencePwFftGridN ...
    && cfg.hartreeGridN == referenceHartreeGridN, ...
    'IGA-PW integration grids do not match the reference run.');
assert(all(isfinite([referenceScfTolEig, referenceScfTolRho, ...
    referencePrimmeTol])) && all([referenceScfTolEig, ...
    referenceScfTolRho, referencePrimmeTol] > 0), ...
    'The saved reference tolerance metadata is invalid.');

cfg.lambdaRef = lambdaRef;
cfg.energyRef = energyRef;
reference = table(lambdaRef, energyRef, cfg.KRef, cfg.pRef, cfg.nelemRef, ...
    string(cfg.referenceRelativeFile), 'VariableNames', ...
    {'lambdaRef', 'energyRef', 'KRef', 'pRef', 'nelemRef', 'sourceFile'});
end

function save_reference_local(cfg, reference)
% Save the reference metadata used by the comparison.

writetable(reference, fullfile(cfg.resultDir, 'reference.csv'));
save(fullfile(cfg.resultDir, 'reference.mat'), 'reference');
end

function results = empty_iga_results_table_local()
% Create the IGA-PW result schema.

results = table( ...
    strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    false(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'method', 'K', 'p', 'refine', 'nelem', 'h', ...
    'nPw', 'nLocal', 'nDof', 'lambda', 'energy', 'lambdaRef', ...
    'energyRef', 'lambdaError', 'energyError', 'scfIters', ...
    'scfConverged', 'finalLambdaChange', 'finalDensityChange', 'totalTime'});
end

function results = empty_dg_results_table_local()
% Create the DG-PW result schema.

results = table( ...
    strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), ...
    false(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'method', 'K', 'n_r', 'L_m', 'nPw', 'nLocal', ...
    'nDof', 'lambda', 'energy', 'lambdaRef', 'energyRef', ...
    'lambdaError', 'energyError', 'scfIters', 'scfConverged', ...
    'finalLambdaChange', 'finalDensityChange', 'totalTime'});
end

function row = iga_row_local(point, cfg)
% Convert one IGA-PW result to a table row.

row = table("IGA-PW", point.K, point.p, point.refine, point.nelem, ...
    point.h, point.nPw, point.nLocal, point.nDof, point.lambda, ...
    point.energy, cfg.lambdaRef, cfg.energyRef, point.lambdaError, ...
    point.energyError, point.scfIters, point.scfConverged, ...
    point.finalLambdaChange, point.finalDensityChange, point.totalTime, ...
    'VariableNames', {'method', 'K', 'p', 'refine', 'nelem', 'h', ...
    'nPw', 'nLocal', 'nDof', 'lambda', 'energy', 'lambdaRef', ...
    'energyRef', 'lambdaError', 'energyError', 'scfIters', ...
    'scfConverged', 'finalLambdaChange', 'finalDensityChange', 'totalTime'});
end

function row = dg_row_local(point, cfg)
% Convert one DG-PW result to a table row.

row = table("DG-PW", point.K, point.q, point.q, ...
    point.nPw, point.nLocal, point.nDof, point.lambda, point.energy, ...
    cfg.lambdaRef, cfg.energyRef, point.lambdaError, point.energyError, ...
    point.scfIters, point.scfConverged, point.finalLambdaChange, ...
    point.finalDensityChange, point.totalTime, 'VariableNames', ...
    {'method', 'K', 'n_r', 'L_m', 'nPw', 'nLocal', 'nDof', 'lambda', ...
    'energy', 'lambdaRef', 'energyRef', 'lambdaError', 'energyError', ...
    'scfIters', 'scfConverged', 'finalLambdaChange', ...
    'finalDensityChange', 'totalTime'});
end

function reset_outputs_local(cfg)
% Remove outputs produced by the eight-case workflow.

names = [run_output_names_local(); {'table.csv'}];
for index = 1:numel(names)
    file = fullfile(cfg.resultDir, names{index});
    if isfile(file)
        delete(file);
    end
end
end

function names = run_output_names_local()
% List files produced by the eight-case computation.

names = {
    'reference.csv'
    'reference.mat'
    'iga.csv'
    'iga.mat'
    'dg.csv'
    'dg.mat'
    'comparison.csv'
    'comparison.mat'
    'run.log'
    };
for caseIndex = 1:4
    names{end + 1} = sprintf('iga_%d.mat', caseIndex); %#ok<AGROW>
    names{end + 1} = sprintf('dg_%d.mat', caseIndex); %#ok<AGROW>
    names{end + 1} = sprintf('dg_%d_times.csv', caseIndex); %#ok<AGROW>
end
end

function ensure_directory_local(directory)
% Create an output directory.

if ~isfolder(directory)
    mkdir(directory);
end
end
