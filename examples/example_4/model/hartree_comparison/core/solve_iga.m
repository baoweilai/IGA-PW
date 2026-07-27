function out = solve_iga(cfg, K, p, refine)

% Solve one retained IGA Hartree-comparison case.
arguments
    cfg struct
    K (1,1) double
    p (1,1) double
    refine (1,1) double
end

% Activate the project and workflow directories.
for k = 1:numel(cfg.projectPaths)
    assert(isfolder(cfg.projectPaths{k}), 'Missing project directory: %s', cfg.projectPaths{k});
    addpath(cfg.projectPaths{k});
end
assert(isfolder(cfg.referenceOperators), ...
    'Missing Example 4 reference operators: %s', cfg.referenceOperators);
addpath(cfg.referenceOperators, '-begin');
workflowPaths = {cfg.commonOperators, cfg.operatorsDir, cfg.solverDir, cfg.coreDir};
for k = 1:numel(workflowPaths)
    assert(isfolder(workflowPaths{k}), 'Missing workflow directory: %s', workflowPaths{k});
    addpath(workflowPaths{k}, '-begin');
end
rehash;
assert(exist('primme_eigs', 'file') == 2, 'PRIMME MATLAB interface is unavailable.');
assert(strcmp(cfg.traceFactorization, 'rfp'), ...
    'IGA-PW uses packed/RFP TB-DG.');
assert(cfg.patternTol == 1e-12, 'The fixed pattern threshold must be 1e-12.');
assert(cfg.targetShift == 0, 'The fixed TB-DG shift must be zero.');

% Build the inner IGA patch and the outer plane-wave basis.
format long;
t0 = tic;
nelem = 2 ^ refine;
fprintf('[IGA-PW] K=%d p=%d refine=%d nelem=%d: building static operators\n', ...
    K, p, refine, nelem);

L = cfg.Lcell;
a = cfg.innerRadius;
innerDomains = [-a, a, -a, a, -a, a];
[nurbsOriginal, nurbsRefine] = build_patch_local(p, nelem, a);
nLocal = nurbsRefine.n_dofs_domains;
hmin = 2*a/nelem;

operatorOpts = struct('use_affine_cube_fast', true);
[Kkin, Mlocal] = generate_K_M_NURBS_3D( ...
    nurbsOriginal, nurbsRefine, cfg.nGauss, operatorOpts);
Kkin = 0.5*(Kkin+Kkin');
Mlocal = 0.5*(Mlocal+Mlocal');

[pVr, nPwVr] = build_pw_ball(cfg.ewaldCutoff);
[pwData, fftTiming] = fftcheb3(cfg, K, innerDomains, pVr, nPwVr);
nPw = pwData.n_pw;
nDof = nLocal+nPw;

% Assemble the interface and static Hamiltonian operators.
faceData = build_face_data_3D(nurbsRefine, a, L, floor(K));
[Pii, Sii] = assemble_local_face_terms_local(faceData, nLocal);
sigma = cfg.penaltyBeta*(1/hmin+K);

fprintf('[IGA-PW] K=%d p=%d nelem=%d: assembling the Ewald IGA block at %.3f s\n', ...
    K, p, nelem, toc(t0));
potentialOpts = struct('use_affine_cube_fast', true);
ewald = struct('L', cfg.Lcell, 'mu', cfg.mu, 'charge', cfg.charge);
[Vlocal, potentialMeta] = assemble_vext_direct_3D( ...
    nurbsOriginal, nurbsRefine, pVr, nPwVr, ewald, ...
    cfg.nGauss, potentialOpts);

KlocalStatic = Kkin+Vlocal-0.25*Sii-0.25*Sii'+sigma*Pii;
KlocalStatic = 0.5*(KlocalStatic+KlocalStatic');
AstaticFun = @(x) apply_global_A_local( ...
    x, KlocalStatic, nLocal, sigma, pwData, faceData);
Bfun = @(x) apply_global_B_local(x, Mlocal, nLocal, pwData);

% Prepare the Hartree grid and the initial density.
gridCache = build_midpoint_grid_cache_3D(L, cfg.hartreeGridN, a);
igaGridEval = build_iga_grid_eval_matrix_3D( ...
    nurbsRefine, gridCache.x_inner, gridCache.y_inner, gridCache.z_inner, a);

pwIdx = (nLocal+1:nDof).';
PGamma = sparse(nDof, nDof);
PGamma(1:nLocal, 1:nLocal) = Pii;

uPrev = build_initial_guess_local(pwData, nLocal, nDof, Bfun);
lambdaPrev = real(uPrev' * AstaticFun(uPrev));
rhoInput = evaluate_density_local( ...
    uPrev, nLocal, pwData, gridCache, igaGridEval);

scfConverged = false;
lambda = lambdaPrev;
energy = 0;
lambdaChange = 0;
densityChange = 0;
precInfo = struct();

% Iterate the Hartree potential and ground-state eigenpair.
for it = 1:cfg.scfMaxit
    VHInput = solve_hartree_fft_local(rhoInput, gridCache);
    VhartreeLocal = assemble_NURBS_potential_interp_3D( ...
        nurbsOriginal, nurbsRefine, VHInput, L, cfg.nGauss, potentialOpts);
    KlocalIter = KlocalStatic + VhartreeLocal;
    KlocalIter = 0.5 * (KlocalIter + KlocalIter');
    pwDataIter = fftcheb3_hartree(cfg, VHInput, innerDomains, pwData);
    clear VHInput VhartreeLocal

    Afun = @(x) apply_global_A_local( ...
        x, KlocalIter, nLocal, sigma, pwDataIter, faceData);
    localDiagonal = diag(KlocalIter-cfg.targetShift*Mlocal);
    pwDiagonal = pwDataIter.stiffDiag ...
        +sigma*pwDataIter.facePenaltyDiagApprox ...
        -cfg.targetShift*pwDataIter.massDiag;
    diagShifted = [localDiagonal(:); pwDiagonal(:)];
    caseInfo = struct('K', K, 'p', p, 'nelem', nelem);
    traceFactorBuilder = @(gamma) build_trace_rfp( ...
        gamma, KlocalIter, Mlocal, pwDataIter, faceData, sigma, ...
        nLocal, nDof, cfg.targetShift, cfg.patternTol, cfg.traceBlockRows);

    fprintf(['[IGA-PW] K=%d p=%d refine=%d SCF=%d: ' ...
        'building packed/RFP TB-DG at %.3f s\n'], ...
        K, p, refine, it, toc(t0));
    [Pfun, precInfo] = tbprec_rfp( ...
        PGamma, pwIdx, diagShifted, nDof, ...
        traceFactorBuilder, cfg, caseInfo);

    ops = struct();
    ops.tol = cfg.primmeTol;
    ops.maxit = cfg.primmeMaxit;
    ops.reportLevel = cfg.primmeReportLevel;
    ops.display = 0;
    ops.isreal = false;
    ops.isdouble = true;
    ops.ishermitian = true;
    ops.v0 = uPrev;

    [uRaw, D] = primme_eigs( ...
        Afun, Bfun, nDof, cfg.numEvals, cfg.targetShift, ...
        ops, cfg.primmeMethod, Pfun);
    rfp_chol_mex('free');
    uRaw = uRaw(:, 1);
    uRaw = uRaw / sqrt(real(uRaw' * Bfun(uRaw)));
    phaseOverlap = uPrev' * Bfun(uRaw);
    if abs(phaseOverlap) > 0
        uRaw = uRaw * conj(phaseOverlap) / abs(phaseOverlap);
    end

    lambda = real(D(1,1));
    rhoOutput = evaluate_density_local( ...
        uRaw, nLocal, pwData, gridCache, igaGridEval);
    VHOutput = solve_hartree_fft_local(rhoOutput, gridCache);
    intRhoVH = real(sum(rhoOutput(:) .* VHOutput(:))) * gridCache.dx ^ 3;
    singleParticle = real(uRaw' * AstaticFun(uRaw));
    energy = 2 * singleParticle + 0.5 * intRhoVH;

    lambdaChange = abs(lambda-lambdaPrev);
    densityChange = norm(rhoOutput(:)-rhoInput(:)) ...
        / max(norm(rhoOutput(:)), 1);
    fprintf(['[IGA-PW] K=%d p=%d refine=%d SCF=%d ' ...
        'lambda=%.16e energy=%.16e dlambda=%.3e drho=%.3e\n'], ...
        K, p, refine, it, lambda, energy, lambdaChange, densityChange);

    uPrev = uRaw;
    lambdaPrev = lambda;
    if it >= 2 && lambdaChange < cfg.scfTolEig ...
            && densityChange < cfg.scfTolRho
        scfConverged = true;
        break;
    end
    rhoInput = (1-cfg.scfBeta)*rhoInput + cfg.scfBeta*rhoOutput;
    clear rhoOutput VHOutput pwDataIter KlocalIter
end
tTotal = toc(t0);

% Validate and package the retained case results.
assert(isfinite(lambda) && isreal(lambda), 'IGA-PW returned a non-finite or non-real eigenvalue.');
assert(isfinite(energy) && isreal(energy), 'IGA-PW returned a non-finite or non-real energy.');
assert(isfinite(tTotal) && tTotal > 0, 'IGA-PW returned an invalid total time.');
assert(nDof == nPw+nLocal, 'IGA-PW degree-of-freedom count is inconsistent.');

out = struct();
out.K = K;
out.p = p;
out.refine = refine;
out.nelem = nelem;
out.h = hmin;
out.nPw = nPw;
out.nLocal = nLocal;
out.nDof = nDof;
out.lambda = lambda;
out.energy = energy;
out.lambdaError = abs(lambda-cfg.lambdaRef);
out.energyError = abs(energy-cfg.energyRef);
out.scfIters = it;
out.scfConverged = scfConverged;
out.finalLambdaChange = lambdaChange;
out.finalDensityChange = densityChange;
out.totalTime = tTotal;
out.sigma = sigma;
out.fftGridN = cfg.fftGridN;
out.chebDegree = cfg.chebDegree;
out.innerQuadN = fftTiming.inner_quad_n;
out.nGamma = precInfo.nGamma;
out.nEta = precInfo.nEta;
out.nnzTrace = precInfo.nnzTrace;
out.factorType = precInfo.factorType;
out.potentialAssembly = potentialMeta.method;
out.potentialAssemblyTime = potentialMeta.t_total;
out.potentialReciprocalModes = potentialMeta.n_reciprocal_modes;
fprintf(['[IGA-PW] completed K=%d p=%d refine=%d nelem=%d ' ...
    'lambda=%.16e energy=%.16e totalTime=%.6f s converged=%d\n'], ...
    K, p, refine, nelem, lambda, energy, tTotal, scfConverged);
end

function [nurbsOriginal, nurbsRefine] = build_patch_local(p, nelem, a)
% Build and refine the inner NURBS patch.
base = [-a, a];
[X, Y, Z] = ndgrid(base, base, base);
controlPoints = zeros(2, 2, 2, 3);
controlPoints(:, :, :, 1) = X;
controlPoints(:, :, :, 2) = Y;
controlPoints(:, :, :, 3) = Z;

nurbsOriginal = struct();
nurbsOriginal.ConPts = controlPoints;
nurbsOriginal.weights = ones(2, 2, 2);
nurbsOriginal.pu = 1;
nurbsOriginal.pv = 1;
nurbsOriginal.pw = 1;
nurbsOriginal.knotU = [0, 0, 1, 1];
nurbsOriginal.knotV = [0, 0, 1, 1];
nurbsOriginal.knotW = [0, 0, 1, 1];

t = p-1;
[knotU, knotV, knotW] = IGADegreeElevVolume( ...
    nurbsOriginal.knotU, nurbsOriginal.knotV, nurbsOriginal.knotW, t);
refinement = struct('mode', 'nelem', 'value', nelem);
nurbsRefine = IGA_3D_Grid(knotU, knotV, knotW, p, p, p, refinement);
end

function [Pii, Sii] = assemble_local_face_terms_local(faceData, nLocal)
% Assemble the local IGA face penalty and flux terms.
Pii = sparse(nLocal, nLocal);
Sii = sparse(nLocal, nLocal);
for iface = 1:numel(faceData)
    F = faceData{iface};
    W = spdiags(F.w, 0, numel(F.w), numel(F.w));
    Pii = Pii+F.TI'*W*F.TI;
    Sii = Sii+F.TI'*W*F.GI;
end
end

function y = apply_global_B_local(x, Mlocal, nLocal, pwData)
% Apply the global mass operator.
xLocal = x(1:nLocal, :);
xPw = x(nLocal+1:end, :);
y = [Mlocal*xLocal; PW3D_apply_mass(xPw, pwData)];
end

function y = apply_global_A_local(x, Klocal, nLocal, sigma, pwData, faceData)
% Apply the global Hamiltonian operator.
xLocal = x(1:nLocal, :);
xPw = x(nLocal+1:end, :);
yLocal = Klocal*xLocal;
[yPw, yLocalFace] = PW3D_apply_stiff(xPw, xLocal, pwData, faceData, sigma);
y = [yLocal+yLocalFace; yPw];
end

function v0 = build_initial_guess_local(pwData, nLocal, nDof, Bfun)
% Build and mass-normalize the initial eigenvector guess.
v0 = zeros(nDof, 1);
v0(1:nLocal) = 1;
k0 = find(all(pwData.k_pw == 0, 2), 1, 'first');
assert(~isempty(k0), 'The PW basis does not contain its zero mode.');
v0(nLocal+k0) = sqrt(pwData.Omega);
v0 = v0/sqrt(real(v0'*Bfun(v0)));
end

function rhoGrid = evaluate_density_local( ...
    u, nLocal, pwData, gridCache, igaGridEval)

% Evaluate the mixed IGA-PW density on the FFT grid.
m = gridCache.mFFT;
kPw = pwData.k_pw;
ii = mod(kPw(:, 1), m) + 1;
jj = mod(kPw(:, 2), m) + 1;
kk = mod(kPw(:, 3), m) + 1;
phase = exp(-1i*pi*sum(kPw, 2)*(1-1/m)) / sqrt(pwData.Omega);
raw = complex(zeros(m, m, m));
raw(sub2ind([m, m, m], ii, jj, kk)) = u(nLocal+1:end) .* phase;
uGrid = ifftn(raw) * m ^ 3;
clear raw

innerValues = igaGridEval * u(1:nLocal);
idx = gridCache.inner_idx;
uGrid(idx, idx, idx) = reshape( ...
    innerValues, gridCache.n_inner, gridCache.n_inner, gridCache.n_inner);
rhoGrid = 2 * abs(uGrid) .^ 2;
end

function VHGrid = solve_hartree_fft_local(rhoGrid, gridCache)
% Solve the periodic Hartree potential on the FFT grid.
m = gridCache.mFFT;
q = ifftshift(-floor(m/2):ceil(m/2)-1);
[Qx, Qy, Qz] = ndgrid(q, q, q);
scale = 2*pi/gridCache.L;
G2 = scale^2 * (Qx.^2 + Qy.^2 + Qz.^2);
rhoHat = fftn(rhoGrid);
VHat = complex(zeros(size(rhoHat)));
nonzeroMode = (G2 > 0);
VHat(nonzeroMode) = 4*pi*rhoHat(nonzeroMode) ./ G2(nonzeroMode);
VHat(~nonzeroMode) = 0;
VHGrid = ifftn(VHat, 'symmetric');
end
