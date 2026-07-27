function out = solve_dg(cfg, K, q)
% Run one DG-PW Hartree-comparison case.

% Activate and verify the nonlinear DG solver.
assert(isfile(cfg.dgScfFile), ...
    'DG-PW nonlinear entry is missing: %s', cfg.dgScfFile);

addpath(cfg.operatorsDir, '-begin');
addpath(cfg.coreDir, '-begin');
rehash;
assert(strcmpi(which('dg_scf'), cfg.dgScfFile), ...
    'MATLAB did not resolve the DG-PW nonlinear entry first.');

% Set the basis sizes and run the nonlinear solve.
nPw = count_pw_ball_local(K);
nLocal = q*(q+1)^2;
nDof = nPw+nLocal;

fprintf('[DG-PW] K=%d N=L=%d nPw=%d nLocal=%d nDof=%d\n', ...
    K, q, nPw, nLocal, nDof);
t0 = tic;
dgOptions = struct( ...
    'scfTolEig', cfg.dgScfTolEig, ...
    'scfTolRho', cfg.dgScfTolRho, ...
    'eigsTol', cfg.dgEigsTol, ...
    'scfMaxit', cfg.scfMaxit, ...
    'mixingNew', cfg.scfBeta);
fprintf(['[DG-PW] tolerances tol_eig=%.1e tol_rho=%.1e ' ...
    'eigs_tol=%.1e maxit=%d\n'], dgOptions.scfTolEig, ...
    dgOptions.scfTolRho, dgOptions.eigsTol, ...
    dgOptions.scfMaxit);
[energyRaw, lambdaRaw, ~, scfIters, finalLambdaChange, ...
    scfConverged, solverMeta] = dg_scf( ...
    cfg.Lcell, cfg.innerRadius, K, q, q, cfg.dgPenalty, dgOptions);
totalTime = toc(t0);

lambda = real(lambdaRaw(1,1));
energy = real(energyRaw(1,1));
clear lambdaRaw energyRaw
assert(isfinite(lambda) && isreal(lambda), ...
    'DG-PW returned an invalid eigenvalue.');
assert(isfinite(energy) && isreal(energy), ...
    'DG-PW returned an invalid Hartree total energy.');
assert(isfinite(totalTime) && totalTime > 0, ...
    'DG-PW returned an invalid total time.');

% Package the convergence and timing results.
out = struct();
out.q = q;
out.K = K;
out.nPw = nPw;
out.nLocal = nLocal;
out.nDof = nDof;
out.lambda = lambda;
out.energy = energy;
out.lambdaError = abs(lambda-cfg.lambdaRef);
out.energyError = abs(energy-cfg.energyRef);
out.scfIters = scfIters;
out.scfConverged = scfConverged;
out.finalLambdaChange = finalLambdaChange;
out.finalDensityChange = solverMeta.finalDensityChange;
out.totalTime = totalTime;
out.solverTime = solverMeta.totalTime;
out.finalEigenResidual = solverMeta.finalEigenResidual;
out.massNormalizationResidual = solverMeta.massNormalizationResidual;
out.stageTimes = solverMeta.stageTimes;
fprintf(['[DG-PW] completed K=%d N=L=%d lambda=%.16e energy=%.16e ' ...
    'totalTime=%.6f s\n'], K, q, lambda, energy, totalTime);
end

function nPw = count_pw_ball_local(K)
% Count the plane-wave modes inside a three-dimensional cutoff ball.
nPw = 0;
N = floor(K);
for ii = -N:N
    jmax = floor(sqrt(K^2-ii^2));
    for jj = -jmax:jmax
        kmax = floor(sqrt(K^2-ii^2-jj^2));
        nPw = nPw+2*kmax+1;
    end
end
end
