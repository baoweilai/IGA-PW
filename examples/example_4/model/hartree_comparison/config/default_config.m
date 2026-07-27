function cfg = default_config()
% Return the fixed periodic-helium comparison configuration.

% Resolve all paths from the configuration file location.
workflowDir = fileparts(fileparts(mfilename('fullpath')));
modelDir = fileparts(workflowDir);
exampleDir = fileparts(modelDir);
examplesDir = fileparts(exampleDir);
projectDir = fileparts(examplesDir);

cfg = struct();
cfg.workflowName = 'hartree_comparison';
cfg.workflowDir = workflowDir;
cfg.exampleDir = exampleDir;
cfg.projectDir = projectDir;
cfg.configDir = fullfile(workflowDir, 'config');
cfg.coreDir = fullfile(workflowDir, 'core');
cfg.operatorsDir = fullfile(workflowDir, 'operators');
cfg.solverDir = fullfile(workflowDir, 'solver');
cfg.commonOperators = fullfile(modelDir, 'common', 'operators');
cfg.resultDir = fullfile(exampleDir, 'data', 'result', cfg.workflowName);
cfg.referenceOperators = fullfile(modelDir, 'reference', 'operators');
cfg.referenceRelativeFile = fullfile('REFERENCE', 'K_30', 'p_2', ...
    'nelem_32', 'run.mat');
cfg.referenceFile = fullfile(exampleDir, 'data', 'result', ...
    cfg.referenceRelativeFile);
cfg.dgScfFile = fullfile(cfg.coreDir, 'dg_scf.m');

cfg.projectPaths = {
    fullfile(projectDir, 'src', 'assembly')
    fullfile(projectDir, 'src', 'dg')
    fullfile(projectDir, 'src', 'iga')
    fullfile(projectDir, 'src', 'pw')
    fullfile(projectDir, 'src', 'nurbs')
    fullfile(projectDir, 'src', 'plotting')
    fullfile(projectDir, 'src', 'utils')
    fullfile(projectDir, 'external', 'primme', 'Matlab')
    };

% Define the retained IGA-PW and DG-PW cases.
cfg.igaCases = [
    10, 2, 2
    10, 2, 4
    15, 2, 2
    15, 2, 4
    ];
cfg.dgCases = [
    10, 2, 2
    10, 3, 3
    15, 2, 2
    15, 3, 3
    ];

% Set the physical and discretization parameters.
cfg.Lcell = 4;
cfg.mu = 5;
cfg.charge = 2;
cfg.innerRadius = 0.2;

cfg.nGauss = 10;
cfg.chebDegree = 80;
cfg.fftGridN = 300;
cfg.hartreeGridN = 300;
cfg.ewaldCutoff = 1;
cfg.penaltyBeta = 20;
cfg.patternTol = 1e-12;
cfg.targetShift = 0;
cfg.traceFactorization = 'rfp';
cfg.traceBlockRows = 64;

% Set the SCF and eigensolver tolerances.
cfg.scfMaxit = 80;
cfg.scfTolEig = 1e-8;
cfg.scfTolRho = 1e-8;
cfg.scfBeta = 0.8;

cfg.numEvals = 1;
cfg.primmeTol = 1e-9;
cfg.primmeMaxit = 1e8;
cfg.primmeMethod = 'DEFAULT_MIN_MATVECS';
cfg.primmeReportLevel = 0;

cfg.dgPenalty = 10000;
cfg.dgScfTolEig = 1e-8;
cfg.dgScfTolRho = 1e-8;
cfg.dgEigsTol = 1e-9;

% Define the reference discretization.
cfg.KRef = 30;
cfg.pRef = 2;
cfg.nelemRef = 32;
end
