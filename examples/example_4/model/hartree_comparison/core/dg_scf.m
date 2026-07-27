function [energy,lambda,phi,iter,error,scfConverged,meta] = ...
    dg_scf( ...
    L, R, N_p, n_r, Lm, sigma, userOptions)
% DG_SCF Solve the nonlinear DG-PW problem.
% The implementation uses the spherical Bessel/harmonic discretization,
% scalar harmonic m-sum, density mixing, and Hartree stopping rule described
% by the comparison workflow.

options = default_options_local(userOptions);

% Set the basis sizes and physical constants.
assert(N_p >= 1 && N_p == floor(N_p), ...
    'The DG-PW solver requires a positive integer K.');
assert(n_r >= 1 && n_r == floor(n_r), ...
    'n_r must be a positive integer.');
assert(Lm >= 0 && Lm == floor(Lm), ...
    'L_m must be a nonnegative integer.');

totalTimer = tic;
N = floor(N_p);
waveScale = 2*pi/L;
volumeOmega = L^3;
volumeInner = 4*pi*R^3/3;
volumeOuter = volumeOmega-volumeInner;
delta = 1e-6;
lmax = Lm;
mu = 5;
charge = 2;

kk = enumerate_pw_ball_local(N_p);
nPw = size(kk, 1);
nAngular = (Lm+1)^2;
nLocal = n_r*nAngular;
n = nPw+nLocal;
zeroPw = find(all(kk == 0, 2), 1);
assert(isscalar(zeroPw), 'The periodic zero mode is missing.');

fprintf(['[DG-PW] K=%d n_r=%d L_m=%d lmax=%d nPw=%d ' ...
    'nLocal=%d nDof=%d\n'], N_p, n_r, Lm, lmax, nPw, nLocal, n);

stageNames = strings(0, 1);
stageSeconds = zeros(0, 1);

% Precompute plane-wave indices, local basis data, and radial kernels.
stageTimer = tic;
pwNormSq = sum(kk.^2, 2);
pwNorm = sqrt(pwNormSq);
pwWaveNorm = waveScale*pwNorm;
[localK, localL, localM, angularL, angularM] = ...
    local_basis_metadata_local(n_r, Lm);

smallSize = 2*N+1;
densitySize = 4*N+1;
rhsSize = 10*N+1;
vhConvSize = 8*N+1;

pwSmallLinear = sub2ind([smallSize, smallSize, smallSize], ...
    kk(:,1)+N+1, kk(:,2)+N+1, kk(:,3)+N+1);
pwDensityLinear = sub2ind([densitySize, densitySize, densitySize], ...
    kk(:,1)+2*N+1, kk(:,2)+2*N+1, kk(:,3)+2*N+1);
pwRhsLinear = sub2ind([rhsSize, rhsSize, rhsSize], ...
    kk(:,1)+5*N+1, kk(:,2)+5*N+1, kk(:,3)+5*N+1);

[densityX, densityY, densityZ] = ndgrid(-2*N:2*N);
densityMask = densityX.^2+densityY.^2+densityZ.^2 ...
    <= (2*N_p)^2;
radialKernel = radial_kernel_local(N, waveScale, R);
kernelFftRhs = fftn(radialKernel, [rhsSize, rhsSize, rhsSize]);
kernelFftVh = fftn(radialKernel, ...
    [vhConvSize, vhConvSize, vhConvSize]);

jValue = zeros(nPw, lmax+1);
djValue = zeros(nPw, lmax+1);
for ell = 0:lmax
    jValue(:,ell+1) = spherical_bessel_original_array_local( ...
        ell, pwWaveNorm*R);
    jPlus = spherical_bessel_original_array_local( ...
        ell, pwWaveNorm*(R+delta));
    jMinus = spherical_bessel_original_array_local( ...
        ell, pwWaveNorm*(R-delta));
    djValue(:,ell+1) = (jPlus-jMinus)/(2*delta);
end

overlapLookup = zeros((2*N)^2+1,1);
overlapLookup(1) = volumeOuter/volumeOmega;
for differenceSquared = 1:(2*N)^2
    differenceWaveNorm = waveScale*sqrt(differenceSquared);
    overlapLookup(differenceSquared+1) = -4*pi*R^2 ...
        *spherical_bessel(1,differenceWaveNorm*R) ...
        /(differenceWaveNorm*volumeOmega);
end

traceYPlus = complex(zeros(nPw,(lmax+1)^2));
nonzeroPw = find(pwNorm > 0);
for ell = 0:lmax
    for emm = -ell:ell
        angularIndex = ell^2+emm+ell+1;
        for index = 1:numel(nonzeroPw)
            q = nonzeroPw(index);
            traceYPlus(q,angularIndex) = spherical_harmonic_xyz( ...
                ell,emm,kk(q,1),kk(q,2),kk(q,3));
        end
    end
end

stageNames(end+1,1) = "basis_and_fft_setup";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

% Assemble the external potential in both subdomains.
stageTimer = tic;
    [vext, Vinner] = build_ewald_potential( ...
    L, R, N, n_r, Lm, mu, charge);
stageNames(end+1,1) = "external_potential";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

H = complex(zeros(n, n));
A = complex(zeros(n, n));
M = zeros(n, n);
pwDifferenceLinear = zeros(nPw,nPw,'uint32');

% Assemble the plane-wave volume matrices.
stageTimer = tic;
for rowFirst = 1:options.pwBlockRows:nPw
    rowLast = min(nPw, rowFirst+options.pwBlockRows-1);
    rows = rowFirst:rowLast;

    dotInteger = kk(rows,:)*kk.';
    differenceSq = pwNormSq(rows)+pwNormSq.'-2*dotInteger;
    differenceSq = max(differenceSq, 0);

    Mblock = reshape(overlapLookup(differenceSq(:)+1), ...
        size(differenceSq));

    kineticH = 0.5*waveScale^2*dotInteger.*Mblock;
    kineticA = waveScale^2*dotInteger.*Mblock;

    dx = kk(:,1).'-kk(rows,1);
    dy = kk(:,2).'-kk(rows,2);
    dz = kk(:,3).'-kk(rows,3);
    vextLinear = sub2ind(size(vext), ...
        dx+2*N+1, dy+2*N+1, dz+2*N+1);
    pwDifferenceLinear(rows,:) = uint32(vextLinear);
    externalBlock = reshape(vext(vextLinear), size(dx));

    H(rows,1:nPw) = kineticH+externalBlock;
    A(rows,1:nPw) = kineticA;
    M(rows,1:nPw) = Mblock;
end

% Add the plane-wave trace and penalty terms.
traceFactor = 0.25*(4*pi)^2*R^2/volumeOmega;
traceTimer = tic;
traceProgressStep = max(1,ceil(nPw/60));
for p = 1:nPw
    if pwNorm(p) == 0
        continue;
    end
    for q = p:nPw
        if pwNorm(q) == 0
            continue;
        end
        sumHarmonicH = 0;
        sumHarmonicA = 0;
        for ell = 0:lmax
            harmonicSum = 0;
            for emm = -ell:ell
                angularIndex = ell^2+emm+ell+1;
                ylm1 = spherical_harmonic_xyz(ell,emm, ...
                    -kk(p,1),-kk(p,2),-kk(p,3));
                harmonicSum = harmonicSum ...
                    +ylm1'*traceYPlus(q,angularIndex);
            end
            jp = jValue(p,ell+1);
            jq = jValue(q,ell+1);
            djp = djValue(p,ell+1);
            djq = djValue(q,ell+1);
            sumHarmonicH = sumHarmonicH+(-1)^ell*harmonicSum ...
                *(jp*djq+djp*jq+4*sigma*jp*jq);
            sumHarmonicA = sumHarmonicA+(-1)^ell*harmonicSum ...
                *(2*jp*djq+2*djp*jq+4*sigma*jp*jq);
        end
        traceValueH = traceFactor*sumHarmonicH;
        traceValueA = traceFactor*sumHarmonicA;
        H(p,q) = H(p,q)+traceValueH;
        A(p,q) = A(p,q)+traceValueA;
        if q ~= p
            H(q,p) = H(q,p)+conj(traceValueH);
            A(q,p) = A(q,p)+conj(traceValueA);
        end
    end
    if mod(p,traceProgressStep) == 0 || p == nPw
        pairFraction = p*(2*nPw-p+1)/(nPw*(nPw+1));
        traceElapsed = toc(traceTimer);
        traceEta = traceElapsed*(1-pairFraction)/pairFraction;
        fprintf(['[DG-PW] harmonic m-sum p=%d/%d ' ...
            'pairs=%.1f%% elapsed=%.1f s eta=%.1f s\n'], ...
            p,nPw,100*pairFraction,traceElapsed,traceEta);
    end
end

zeroTraceH = pi*R^2/volumeOmega ...
    *(djValue(nonzeroPw,1)+4*sigma*jValue(nonzeroPw,1));
zeroTraceA = pi*R^2/volumeOmega ...
    *(2*djValue(nonzeroPw,1)+4*sigma*jValue(nonzeroPw,1));
H(nonzeroPw,zeroPw) = H(nonzeroPw,zeroPw)+zeroTraceH;
H(zeroPw,nonzeroPw) = H(zeroPw,nonzeroPw)+zeroTraceH';
A(nonzeroPw,zeroPw) = A(nonzeroPw,zeroPw)+zeroTraceA;
A(zeroPw,nonzeroPw) = A(zeroPw,nonzeroPw)+zeroTraceA';
H(zeroPw,zeroPw) = H(zeroPw,zeroPw) ...
    +4*pi*R^2*sigma/volumeOmega;
A(zeroPw,zeroPw) = A(zeroPw,zeroPw) ...
    +4*pi*R^2*sigma/volumeOmega;

clear vext densityX densityY densityZ overlapLookup traceYPlus
stageNames(end+1,1) = "pw_pw_matrix";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

% Assemble the local radial and angular matrices.
stageTimer = tic;
intR = zeros(n_r, n_r);
intNoR2 = zeros(n_r, n_r);
intDR = zeros(n_r, n_r);
for k1 = 1:n_r
    for k2 = 1:n_r
        intR(k1,k2) = int_overlap_radius(k1,k2,R,options.nIntegral);
        intNoR2(k1,k2) = int_overlap(k1,k2,R,options.nIntegral);
        intDR(k1,k2) = int_overlap_radius_d(k1,k2,R,options.nIntegral);
    end
end

Hlocal = Vinner;
Alocal = zeros(nLocal, nLocal);
Mlocal = zeros(nLocal, nLocal);
for pLocal = 1:nLocal
    k1 = localK(pLocal);
    l1 = localL(pLocal);
    m1 = localM(pLocal);
    for qLocal = 1:nLocal
        k2 = localK(qLocal);
        if l1 == localL(qLocal) && m1 == localM(qLocal)
            Mlocal(pLocal,qLocal) = intR(k1,k2);
            kinetic = l1*(l1+1)*intNoR2(k1,k2)+intDR(k1,k2);
            Hlocal(pLocal,qLocal) = Hlocal(pLocal,qLocal)+0.5*kinetic;
            Alocal(pLocal,qLocal) = kinetic;

            rp = -r_basis(k1,R,R);
            rq = -r_basis(k2,R,R);
            drp = (r_basis(k1,R+delta,R)-r_basis(k1,R-delta,R)) ...
                /(2*delta);
            drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R)) ...
                /(2*delta);
            Hlocal(pLocal,qLocal) = Hlocal(pLocal,qLocal) ...
                +0.25*(rp*drq+drp*rq+4*sigma*rp*rq)*R^2;
            Alocal(pLocal,qLocal) = Alocal(pLocal,qLocal) ...
                +0.25*(2*rp*drq+2*drp*rq+4*sigma*rp*rq)*R^2;
        end
    end
end
localRange = nPw+(1:nLocal);
H(localRange,localRange) = Hlocal;
A(localRange,localRange) = Alocal;
M(localRange,localRange) = Mlocal;
clear Vinner Hlocal Alocal Mlocal
stageNames(end+1,1) = "local_matrix";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

% Assemble the coupling between local and plane-wave bases.
stageTimer = tic;
Yminus = complex(zeros(nPw,nAngular));
nonzeroPw = find(pwNorm > 0);
for angularIndex = 1:nAngular
    ell = angularL(angularIndex);
    emm = angularM(angularIndex);
    for index = 1:numel(nonzeroPw)
        p = nonzeroPw(index);
        Yminus(p,angularIndex) = spherical_harmonic_xyz( ...
            ell,emm,-kk(p,1),-kk(p,2),-kk(p,3));
    end
end

for qLocal = 1:nLocal
    k2 = localK(qLocal);
    ell = localL(qLocal);
    emm = localM(qLocal);
    angularIndex = ell^2+emm+ell+1;
    rq = -r_basis(k2,R,R);
    drq = (r_basis(k2,R+delta,R)-r_basis(k2,R-delta,R))/(2*delta);

    crossH = complex(zeros(nPw,1));
    crossA = complex(zeros(nPw,1));
    nz = pwNorm > 0;
    jp = jValue(nz,ell+1);
    djp = djValue(nz,ell+1);
    angularValue = Yminus(nz,angularIndex);
    prefactor = pi*R^2*(1i)^ell/sqrt(volumeOmega);
    crossH(nz) = prefactor*angularValue ...
        .*(jp*drq+djp*rq+4*sigma*jp*rq);
    crossA(nz) = prefactor*angularValue ...
        .*(2*jp*drq+2*djp*rq+4*sigma*jp*rq);

    if ell == 0
        zeroPrefactor = 0.25*R^2*sqrt(4*pi)/sqrt(volumeOmega);
        crossH(zeroPw) = zeroPrefactor*(drq+4*sigma*rq);
        crossA(zeroPw) = zeroPrefactor*(2*drq+4*sigma*rq);
    end

    qGlobal = nPw+qLocal;
    H(1:nPw,qGlobal) = crossH;
    H(qGlobal,1:nPw) = crossH';
    A(1:nPw,qGlobal) = crossA;
    A(qGlobal,1:nPw) = crossA';
end
clear Yminus jValue djValue
stageNames(end+1,1) = "pw_local_matrix";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

% Solve the initial generalized eigenvalue problem.
stageTimer = tic;
constantMode = zeros(n,1);
constantMode(zeroPw) = sqrt(volumeOmega);
constantMode(nPw+1) = R*sqrt(4*pi);
constantMass = M*constantMode;
constantNorm = real(constantMode'*constantMass);
zeroModeResidual = norm(A*constantMode)/max(norm(constantMass),1);
A = 0.5*(A+A');
A = A+(constantMass*constantMass')/constantNorm;
H = 0.5*(H+H');

eigsOptions = struct('tol',options.eigsTol);
[phi, lambdaMatrix] = eigs(H, M, 1, -100, eigsOptions);
lambdaNew = real(lambdaMatrix(1,1));
c = phi;
nrm = sqrt(real(c'*M*c));
stageNames(end+1,1) = "constant_mode_and_initial_eigs";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));
fprintf('[DG-PW] Hartree constant-mode residual=%.6e\n', ...
    zeroModeResidual);

% Build the local contractions used by the nonlinear terms.
stageTimer = tic;
[rhsLocalMatrix, hartreeLocalMatrix] = local_triple_matrices_local( ...
    n_r, Lm, R, options.nIntegral, angularL, angularM);
stageNames(end+1,1) = "local_triple_tensors";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));

% Factor the Poisson operator used in every SCF iteration.
stageTimer = tic;
[poissonL, cholFlag] = chol(A, 'lower');
if cholFlag == 0
    poissonMethod = "cholesky";
    poissonDecomposition = [];
else
    poissonMethod = "lu";
    poissonL = [];
    poissonDecomposition = decomposition(A, 'lu');
end
clear A
stageNames(end+1,1) = "poisson_factorization";
stageSeconds(end+1,1) = toc(stageTimer);
print_stage_local(stageNames(end), stageSeconds(end));
fprintf('[DG-PW] Poisson factorization=%s\n', ...
    poissonMethod);

lambdaOld = 0;
finalLambdaChange = abs(lambdaNew-lambdaOld);
finalDensityChange = inf;
fIn = complex(zeros(n,1));
iteration = 0;

% Iterate the density, Hartree potential, and ground-state eigenpair.
scfTimer = tic;
Hiter = complex(zeros(n,n));
while (finalLambdaChange > options.scfTolEig ...
        || finalDensityChange > options.scfTolRho) ...
        && iteration < options.scfMaxit
    iterationTimer = tic;

    coefficientCube = complex(zeros(smallSize,smallSize,smallSize));
    coefficientCube(pwSmallLinear) = c(1:nPw)/nrm;
    reversedCoefficient = conj(flip3_local(coefficientCube));
    rhoCube = ifftn( ...
        fftn(coefficientCube,[densitySize,densitySize,densitySize]) ...
         .*fftn(reversedCoefficient,[densitySize,densitySize,densitySize]));
    rhoCube(~densityMask) = 0;

    rhsConvolution = ifftn( ...
        fftn(rhoCube,[rhsSize,rhsSize,rhsSize]).*kernelFftRhs);
    rhoAtPw = rhoCube(pwDensityLinear);
    radialSumAtPw = rhsConvolution(pwRhsLinear);
    fPw = rhoAtPw*4*pi*volumeOuter/(volumeOmega^(3/2)) ...
        -radialSumAtPw*(4*pi)^2/(volumeOmega^(3/2));

    cLocal = c(localRange);
    localDensityMatrix = cLocal*cLocal';
    fLocal = 4*pi/nrm^2 ...
        *(rhsLocalMatrix*localDensityMatrix(:));
    fOut = 2*[fPw; fLocal];
    if iteration == 0
        fIn = fOut;
        finalDensityChange = inf;
    else
        finalDensityChange = norm(fOut-fIn)/max(norm(fOut),1);
    end
    fIn = (1-options.mixingNew)*fIn+options.mixingNew*fOut;

    fProjected = fIn-constantMass*(constantMode'*fIn)/constantNorm;
    if poissonMethod == "cholesky"
        X = poissonL'\(poissonL\fProjected);
    else
        X = poissonDecomposition\fProjected;
    end
    X = X-constantMode*(constantMode'*M*X)/constantNorm;

    Xcube = complex(zeros(smallSize,smallSize,smallSize));
    Xcube(pwSmallLinear) = X(1:nPw);
    Xconvolution = ifftn( ...
        fftn(Xcube,[vhConvSize,vhConvSize,vhConvSize]).*kernelFftVh);
    convolutionSlice = 2*N+1:6*N+1;
    vH = -4*pi/(volumeOmega^(3/2))*flip3_local( ...
        Xconvolution(convolutionSlice,convolutionSlice,convolutionSlice));
    centralSlice = N+1:3*N+1;
    vH(centralSlice,centralSlice,centralSlice) = ...
        vH(centralSlice,centralSlice,centralSlice) ...
        +flip3_local(Xcube)*volumeOuter/(volumeOmega^(3/2));

    Hiter(:) = 0;
    Hiter(1:nPw,1:nPw) = vH(pwDifferenceLinear);
    Hiter(localRange,localRange) = reshape( ...
        hartreeLocalMatrix*X(localRange), nLocal, nLocal);
    Hiter = 0.5*(Hiter+Hiter');
    Hnew = H+Hiter;
    Hnew = 0.5*(Hnew+Hnew');

    [phi, lambdaMatrix] = eigs(Hnew, M, 1, -10, eigsOptions);

    iteration = iteration+1;
    lambdaOld = lambdaNew;
    lambdaNew = real(lambdaMatrix(1,1));
    c = phi;
    nrm = sqrt(real(c'*M*c));
    finalLambdaChange = abs(lambdaNew-lambdaOld);
    iterationTime = toc(iterationTimer);

    fprintf(['[DG-PW] SCF=%d/%d lambda=%.16e abs_dlambda=%.3e ' ...
        'rho=%.3e tol_eig=%.1e tol_rho=%.1e eigs_tol=%.1e ' ...
        'time=%.6f s\n'], iteration, options.scfMaxit, lambdaNew, ...
        finalLambdaChange, finalDensityChange, options.scfTolEig, ...
        options.scfTolRho, options.eigsTol, iterationTime);
end
scfSeconds = toc(scfTimer);
stageNames(end+1,1) = "scf";
stageSeconds(end+1,1) = scfSeconds;
print_stage_local(stageNames(end), stageSeconds(end));

% Compute the final residuals and total energy.
scfConverged = finalLambdaChange <= options.scfTolEig ...
    && finalDensityChange <= options.scfTolRho;
cNormalized = c/nrm;
finalEigenResidual = norm(Hnew*cNormalized-lambdaNew*M*cNormalized) ...
    /max([norm(Hnew*cNormalized), ...
    abs(lambdaNew)*norm(M*cNormalized),1]);
massNormalizationResidual = abs(real(cNormalized'*M*cNormalized)-1);
hartreeOrbital = real(c'*Hiter*c)/real(nrm^2);
energy = 2*lambdaNew-hartreeOrbital;
totalSeconds = toc(totalTimer);

% Package the timing data.
stageTimes = table(stageNames, stageSeconds, ...
    'VariableNames', {'stage','seconds'});

out = struct();
out.K = N_p;
out.n_r = n_r;
out.L_m = Lm;
out.nPw = nPw;
out.nLocal = nLocal;
out.nDof = n;
out.lambda = lambdaNew;
out.energy = energy;
out.scfIters = iteration;
out.scfConverged = scfConverged;
out.finalLambdaChange = finalLambdaChange;
out.finalDensityChange = finalDensityChange;
out.finalEigenResidual = finalEigenResidual;
out.massNormalizationResidual = massNormalizationResidual;
out.totalTime = totalSeconds;
out.zeroModeResidual = zeroModeResidual;
out.poissonMethod = poissonMethod;
out.stageTimes = stageTimes;
out.finalEigenvector = c;
out.options = options;

fprintf(['[DG-PW] completed K=%d n_r=%d L_m=%d ' ...
    'lambda=%.16e energy=%.16e SCF=%d converged=%d totalTime=%.6f s\n'], ...
    N_p, n_r, Lm, out.lambda, out.energy, out.scfIters, ...
    out.scfConverged, out.totalTime);

energy = out.energy;
lambda = out.lambda;
phi = out.finalEigenvector;
iter = out.scfIters;
error = out.finalLambdaChange;
scfConverged = out.scfConverged;
meta = out;
end

function options = default_options_local(userOptions)
% Fill the SCF options with default values.
options = struct();
options.nIntegral = 100;
options.scfTolEig = 1e-10;
options.scfTolRho = 1e-10;
options.eigsTol = 1e-12;
options.scfMaxit = 80;
options.mixingNew = 0.8;
options.pwBlockRows = 128;

names = fieldnames(userOptions);
for index = 1:numel(names)
    options.(names{index}) = userOptions.(names{index});
end
end

function kk = enumerate_pw_ball_local(cutoff)
% Enumerate integer wave vectors inside the cutoff ball.
N = floor(cutoff);
kk = zeros((2*N+1)^3,3);
count = 0;
for ii = -N:N
    jmax = floor(sqrt(cutoff^2-ii^2));
    for jj = -jmax:jmax
        kmax = floor(sqrt(cutoff^2-ii^2-jj^2));
        for kz = -kmax:kmax
            count = count+1;
            kk(count,:) = [ii,jj,kz];
        end
    end
end
kk = kk(1:count,:);
end

function values = spherical_bessel_original_array_local(ell, z)
% Evaluate a spherical Bessel array with the fixed recurrence.
values = zeros(size(z));
for index = 1:numel(z)
    if z(index) == 0
        values(index) = double(ell == 0);
    else
        values(index) = spherical_bessel(ell,z(index));
    end
end
end

function kernel = radial_kernel_local(N, waveScale, R)
% Build the radial Coulomb kernel.
coordinates = -3*N:3*N;
[x,y,z] = ndgrid(coordinates,coordinates,coordinates);
squaredNorm = x.^2+y.^2+z.^2;
uniqueSquaredNorm = unique(squaredNorm(:));
radialLookup = zeros(max(uniqueSquaredNorm)+1,1);
for index = 2:numel(uniqueSquaredNorm)
    squared = uniqueSquaredNorm(index);
    waveNorm = waveScale*sqrt(squared);
    radialLookup(squared+1) = R^2 ...
        *spherical_bessel(1,waveNorm*R)/waveNorm;
end
kernel = radialLookup(squaredNorm+1);
end

function flipped = flip3_local(values)
% Reverse an array along all three dimensions.
flipped = flip(flip(flip(values,1),2),3);
end

function [localK,localL,localM,angularL,angularM] = ...
    local_basis_metadata_local(n_r,Lm)
% Enumerate the local radial and angular basis metadata.
nAngular = (Lm+1)^2;
nLocal = n_r*nAngular;
angularL = zeros(nAngular,1);
angularM = zeros(nAngular,1);
for ell = 0:Lm
    for emm = -ell:ell
        index = ell^2+emm+ell+1;
        angularL(index) = ell;
        angularM(index) = emm;
    end
end

localK = zeros(nLocal,1);
localL = zeros(nLocal,1);
localM = zeros(nLocal,1);
for radial = 1:n_r
    range = (radial-1)*nAngular+(1:nAngular);
    localK(range) = radial;
    localL(range) = angularL;
    localM(range) = angularM;
end
end

function [rhsMatrix,hartreeMatrix] = local_triple_matrices_local( ...
    n_r,Lm,R,nIntegral,angularL,angularM)
% Assemble local nonlinear triple-product matrices.
nAngular = (Lm+1)^2;
nLocal = n_r*nAngular;

% Integrate the radial triple products.
int3r = zeros(n_r,n_r,n_r);
for k1 = 1:n_r
    for k2 = 1:n_r
        for k3 = 1:n_r
            int3r(k1,k2,k3) = int_overlap_radius_3( ...
                k1,k2,k3,R,nIntegral);
        end
    end
end

% Assemble the angular coupling tensors from Wigner symbols.
rhsAngular = zeros(nAngular,nAngular,nAngular);
hartreeAngular = zeros(nAngular,nAngular,nAngular);
for a1 = 1:nAngular
    l1 = angularL(a1);
    m1 = angularM(a1);
    for a2 = 1:nAngular
        l2 = angularL(a2);
        m2 = angularM(a2);
        for a3 = 1:nAngular
            l3 = angularL(a3);
            m3 = angularM(a3);
            if m1-m2-m3 ~= 0
                continue;
            end
            w0 = Wigner3j(l1,l2,l3,0,0,0);
            if w0 == 0
                continue;
            end
            common = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)) ...
                /sqrt(4*pi)*w0;
            rhsAngular(a3,a1,a2) = (-1)^(m2+m3)*common ...
                *Wigner3j(l1,l2,l3,m1,-m2,-m3);
            hartreeAngular(a1,a2,a3) = (-1)^m1*common ...
                *Wigner3j(l1,l2,l3,-m1,m2,m3);
        end
    end
end

% Combine the radial and angular factors in local basis order.
rhsTensor = zeros(nLocal,nLocal,nLocal);
hartreeTensor = zeros(nLocal,nLocal,nLocal);
for k1 = 1:n_r
    pRange = (k1-1)*nAngular+(1:nAngular);
    for k2 = 1:n_r
        qRange = (k2-1)*nAngular+(1:nAngular);
        for k3 = 1:n_r
            iRange = (k3-1)*nAngular+(1:nAngular);
            radial = int3r(k1,k2,k3);
            rhsTensor(iRange,pRange,qRange) = radial*rhsAngular;
            hartreeTensor(pRange,qRange,iRange) = radial*hartreeAngular;
        end
    end
end

rhsMatrix = reshape(rhsTensor,nLocal,nLocal*nLocal);
hartreeMatrix = reshape(hartreeTensor,nLocal*nLocal,nLocal);
end

function print_stage_local(name, seconds)
% Print the elapsed time for one SCF stage.
fprintf('[DG-PW] stage=%s time=%.6f s\n',name,seconds);
end
