function report = build_trace_rfp( ...
    gamma, Klocal, Mlocal, pwData, faceData, sigma, ...
    nLocal, nDof, targetShift, patternTol, blockRows)

% Build the thresholded upper-triangular RFP trace factor.
arguments
    gamma (:,1) double
    Klocal
    Mlocal
    pwData struct
    faceData cell
    sigma (1,1) double
    nLocal (1,1) double
    nDof (1,1) double
    targetShift (1,1) double
    patternTol (1,1) double
    blockRows (1,1) double
end

% Partition the interface indices and prepare the face operators.
assert(targetShift == 0, 'The RFP trace builder requires targetShift=0.');
assert(patternTol == 1e-12, 'The RFP trace builder requires patternTol=1e-12.');
assert(blockRows >= 1 && blockRows == floor(blockRows), ...
    'traceBlockRows must be a positive integer.');

tBuild = tic;
interfaceLocal = gamma(gamma <= nLocal);
pwGlobal = gamma(gamma > nLocal);
expectedPw = (nLocal+1:nDof).';
assert(isequal(pwGlobal, expectedPw), ...
    'TB-DG gamma must contain every PW degree of freedom.');

nIga = numel(interfaceLocal);
nPw = pwData.n_pw;
nGamma = numel(gamma);
assert(nGamma == nIga+nPw, 'The RFP trace dimension is inconsistent.');

rfp_chol_mex('init', nGamma);
faceOps = prepare_face_ops_local(faceData, pwData, sigma);
nnzUpper = 0;
nnzDiagonal = 0;
minAbsKept = inf;

% Insert the local-interface block into the packed upper triangle.
Aii = Klocal(interfaceLocal, interfaceLocal) ...
    -targetShift*Mlocal(interfaceLocal, interfaceLocal);
Aii = full(0.5*(Aii+Aii'));
[Aii, countUpper, countDiagonal, minKept] = ...
    threshold_upper_local(Aii, (1:nIga).', (1:nIga).', patternTol);
rfp_chol_mex('put', (1:nIga).', (1:nIga).', Aii);
nnzUpper = nnzUpper+countUpper;
nnzDiagonal = nnzDiagonal+countDiagonal;
minAbsKept = min(minAbsKept, minKept);
clear Aii

% Insert the mixed local and plane-wave blocks.
for firstRow = 1:blockRows:nPw
    rows = (firstRow:min(firstRow+blockRows-1, nPw)).';
    ApiRows = mixed_rows_local(rows, interfaceLocal, faceOps, pwData, sigma);
    mixedUpper = ApiRows';
    keep = abs(mixedUpper) >= patternTol;
    mixedUpper(~keep) = 0;
    nnzUpper = nnzUpper+nnz(keep);
    if any(keep(:))
        minAbsKept = min(minAbsKept, min(abs(mixedUpper(keep))));
    end
    rfp_chol_mex('put', (1:nIga).', nIga+rows, mixedUpper);
end

% Insert the upper-triangular plane-wave blocks.
for firstRow = 1:blockRows:nPw
    rows = (firstRow:min(firstRow+blockRows-1, nPw)).';
    cols = (firstRow:nPw).';
    AppTile = pw_tile_local(rows, cols, pwData, faceOps);
    [AppTile, countUpper, countDiagonal, minKept] = ...
        threshold_upper_local(AppTile, nIga+rows, nIga+cols, patternTol);
    rfp_chol_mex('put', nIga+rows, nIga+cols, AppTile);
    nnzUpper = nnzUpper+countUpper;
    nnzDiagonal = nnzDiagonal+countDiagonal;
    minAbsKept = min(minAbsKept, minKept);

    if mod(ceil(firstRow/blockRows), 32) == 0 ...
            || rows(end) == nPw
        fprintf('[RFP-TRACE] PW rows=%d/%d elapsed=%.3f s\n', ...
            rows(end), nPw, toc(tBuild));
    end
end

assemblyTime = toc(tBuild);

% Factor the assembled RFP matrix.
tFactor = tic;
pflag = rfp_chol_mex('factor');
factorTime = toc(tFactor);
assert(pflag == 0, 'RFP Cholesky returned pflag=%d.', pflag);
rfpInfo = rfp_chol_mex('info');

% Package the sparsity and factorization statistics.
report = struct();
report.nnzTrace = 2*nnzUpper-nnzDiagonal;
report.nnzUpperTrace = nnzUpper;
report.nnzDiagonalTrace = nnzDiagonal;
report.minAbsKept = minAbsKept;
report.nGammaIga = nIga;
report.nPw = nPw;
report.nGamma = nGamma;
report.assemblyTime = assemblyTime;
report.factorTime = factorTime;
report.rfpElements = rfpInfo.elements;
report.rfpBytes = rfpInfo.bytes;
report.rfpGiB = rfpInfo.bytes/2^30;
report.putCalls = rfpInfo.putCalls;
report.factorized = logical(rfpInfo.factorized);

fprintf(['[RFP-TRACE] completed nGamma=%d nnz=%0.0f ' ...
    'assembly=%.3f s factor=%.3f s storage=%.6f GiB\n'], ...
    nGamma, report.nnzTrace, assemblyTime, factorTime, report.rfpGiB);
end

function faceOps = prepare_face_ops_local(faceData, pwData, sigma)
% Precompute the face operators used by trace-factor blocks.
faceOps = cell(numel(faceData), 1);
k = pwData.k_pw;
N = pwData.N;
alpha = pwData.alpha;
normalModes = (-N:N).';

for iface = 1:numel(faceData)
    F = faceData{iface};
    assert(isfield(F, 'tensor'), 'The fixed tensor face representation is required.');
    Wmat = reshape(F.w(:), F.nq1, F.nq2);
    w1 = Wmat(:, 1);
    w2 = Wmat(1, :).'/Wmat(1, 1);
    T = F.tensor;
    W1B1 = spdiags(T.w1, 0, numel(T.w1), numel(T.w1))*T.B1;
    W2B2 = spdiags(T.w2, 0, numel(T.w2), numel(T.w2))*T.B2;

    op = struct();
    op.G1 = F.E1'*bsxfun(@times, F.E1, w1);
    op.G2 = F.E2'*bsxfun(@times, F.E2, w2);
    op.C1 = F.E1'*W1B1;
    op.C2 = F.E2'*W2B2;
    op.T = T;
    op.type = F.type;
    op.fixedCoord = F.fixedCoord;
    switch F.type
        case 'x'
            op.tang1 = k(:, 2)+N+1;
            op.tang2 = k(:, 3)+N+1;
            op.kn = k(:, 1);
        case 'y'
            op.tang1 = k(:, 1)+N+1;
            op.tang2 = k(:, 3)+N+1;
            op.kn = k(:, 2);
        case 'z'
            op.tang1 = k(:, 1)+N+1;
            op.tang2 = k(:, 2)+N+1;
            op.kn = k(:, 3);
        otherwise
            error('Unsupported face type: %s', F.type);
    end
    op.normalFactor = 1i*alpha*F.normalSign*op.kn;
    normalFactorModes = 1i*alpha*F.normalSign*normalModes;
    op.normalIndex = op.kn+N+1;
    op.phaseForward = exp(1i*alpha*F.fixedCoord ...
        *(normalModes.'-normalModes));
    op.phaseReverse = exp(1i*alpha*F.fixedCoord ...
        *(normalModes-normalModes.'));
    op.coefficientForward = sigma ...
        +0.25*(normalFactorModes.'+conj(normalFactorModes));
    op.coefficientReverse = sigma ...
        +0.25*(normalFactorModes+conj(normalFactorModes).');
    faceOps{iface} = op;
end
end

function ApiRows = mixed_rows_local( ...
    rows, interfaceLocal, faceOps, pwData, sigma)

% Assemble mixed trace rows coupling the local and plane-wave bases.
ApiRows = complex(zeros(numel(rows), numel(interfaceLocal)));
alpha = pwData.alpha;

for iface = 1:numel(faceOps)
    op = faceOps{iface};
    tang1Rows = op.tang1(rows);
    tang2Rows = op.tang2(rows);
    pairs = unique([tang1Rows, tang2Rows], 'rows');
    for ig = 1:size(pairs, 1)
        i1 = pairs(ig, 1);
        i2 = pairs(ig, 2);
        positions = find(tang1Rows == i1 & tang2Rows == i2);
        pwRows = rows(positions);
        switch op.type
            case 'x'
                baseTI = kron(op.C2(i2, :), kron(op.C1(i1, :), op.T.Bfix));
                baseGI = kron(op.C2(i2, :), kron(op.C1(i1, :), op.T.Dfix));
            case 'y'
                baseTI = kron(op.C2(i2, :), kron(op.T.Bfix, op.C1(i1, :)));
                baseGI = kron(op.C2(i2, :), kron(op.T.Dfix, op.C1(i1, :)));
            case 'z'
                baseTI = kron(op.T.Bfix, kron(op.C2(i2, :), op.C1(i1, :)));
                baseGI = kron(op.T.Dfix, kron(op.C2(i2, :), op.C1(i1, :)));
        end
        baseTI = full(baseTI(interfaceLocal));
        baseGI = full(baseGI(interfaceLocal));
        phase = exp(1i*alpha*op.fixedCoord*op.kn(pwRows))/sqrt(pwData.Omega);
        coeff = conj(phase(:));
        normalFactor = op.normalFactor(pwRows);
        coeffTI = -(0.25*conj(normalFactor(:))+sigma).*coeff;
        coeffGI = 0.25*coeff;
        ApiRows(positions, :) = ApiRows(positions, :) ...
            +coeffGI*baseGI+coeffTI*baseTI;
    end
end
end

function blockValues = pw_tile_local(rows, cols, pwData, faceOps)
% Assemble one plane-wave trace block.
k = pwData.k_pw;
alpha = pwData.alpha;
qN = floor((size(pwData.Uker, 1)-1)/2);
kRows = k(rows, :);
kCols = k(cols, :);

diffx = kRows(:, 1)-kCols(:, 1).';
diffy = kRows(:, 2)-kCols(:, 2).';
diffz = kRows(:, 3)-kCols(:, 3).';
forwardIndex = sub2ind(size(pwData.Uker), ...
    diffx+qN+1, diffy+qN+1, diffz+qN+1);
reverseIndex = sub2ind(size(pwData.Uker), ...
    -diffx+qN+1, -diffy+qN+1, -diffz+qN+1);
dotK = kRows(:, 1)*kCols(:, 1).'+kRows(:, 2)*kCols(:, 2).' ...
    +kRows(:, 3)*kCols(:, 3).';

coreForward = -dotK.*pwData.Uker(forwardIndex);
coreReverse = -dotK.*pwData.Uker(reverseIndex);
for row = 1:numel(rows)
    columnPosition = rows(row)-cols(1)+1;
    diagonalTerm = sum(k(rows(row), :).^2);
    coreForward(row, columnPosition) = ...
        coreForward(row, columnPosition)+diagonalTerm;
    coreReverse(row, columnPosition) = ...
        coreReverse(row, columnPosition)+diagonalTerm;
end
rawForward = 0.5*alpha^2*coreForward+pwData.Vker(forwardIndex);
rawReverse = 0.5*alpha^2*coreReverse+pwData.Vker(reverseIndex);

left = 0.5*(rawForward+conj(rawReverse));
rightReverse = 0.5*(rawReverse+conj(rawForward));

for iface = 1:numel(faceOps)
    [faceForward, faceReverse] = face_tiles_local( ...
        rows, cols, faceOps{iface}, pwData);
    left = left+faceForward;
    rightReverse = rightReverse+faceReverse;
end

blockValues = 0.5*(left+conj(rightReverse));
end

function [faceForward, faceReverse] = face_tiles_local( ...
    rows, cols, op, pwData)

% Build the forward and reverse coupling tiles for one face.
Omega = pwData.Omega;
tang1Rows = op.tang1(rows);
tang2Rows = op.tang2(rows);
tang1Cols = op.tang1(cols);
tang2Cols = op.tang2(cols);
normalIndexRows = op.normalIndex(rows);
normalIndexCols = op.normalIndex(cols);

phaseForward = op.phaseForward(normalIndexRows, normalIndexCols);
penaltyForward = phaseForward ...
    .* op.G1(tang1Rows, tang1Cols) ...
    .* op.G2(tang2Rows, tang2Cols)/Omega;
coefficientForward = op.coefficientForward( ...
    normalIndexRows, normalIndexCols);
faceForward = penaltyForward.*coefficientForward;

phaseReverse = op.phaseReverse(normalIndexRows, normalIndexCols);
penaltyReverse = phaseReverse ...
    .* op.G1(tang1Cols, tang1Rows).' ...
    .* op.G2(tang2Cols, tang2Rows).'/Omega;
coefficientReverse = op.coefficientReverse( ...
    normalIndexRows, normalIndexCols);
faceReverse = penaltyReverse.*coefficientReverse;
end

function [block, countUpper, countDiagonal, minKept] = ...
    threshold_upper_local(block, globalRows, globalCols, patternTol)

% Threshold and retain the upper-triangular block entries.
absBlock = abs(block);
keep = absBlock >= patternTol;
block(~keep) = 0;
upperMask = globalRows <= globalCols.';
keptUpper = keep & upperMask;
countUpper = nnz(keptUpper);

diagonalMask = globalRows == globalCols.';
countDiagonal = nnz(keep & diagonalMask);
if any(keptUpper(:))
    minKept = min(absBlock(keptUpper));
else
    minKept = inf;
end
end
