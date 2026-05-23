function [yP, yI_face] = PW3D_apply_stiff(xP, xI, pwData, faceData, sigma)
%Apply the 3-D plane-wave stiffness operator.

[n_pw, nb] = size(xP);
nI = size(xI,1);

yP      = zeros(n_pw, nb);
yI_face = zeros(nI,   nb);

k2cube = pwData.kxCube.^2 + pwData.kyCube.^2 + pwData.kzCube.^2;

for ib = 1:nb
    X  = ball_to_cube_single(xP(:,ib), pwData);
    xIi = xI(:,ib);

    % 1) Volume part on PW block
    Y = 0.5 * pwData.alpha^2 * (k2cube .* X);

    UXx = conv3_same_single(pwData.kxCube .* X, pwData.Ufft, pwData.convN, pwData.K);
    UXy = conv3_same_single(pwData.kyCube .* X, pwData.Ufft, pwData.convN, pwData.K);
    UXz = conv3_same_single(pwData.kzCube .* X, pwData.Ufft, pwData.convN, pwData.K);

    Y = Y - 0.5 * pwData.alpha^2 * ...
        (pwData.kxCube .* UXx + pwData.kyCube .* UXy + pwData.kzCube .* UXz);

    VX = conv3_same_single(X, pwData.Vfft, pwData.convN, pwData.K);
    Y  = Y + VX;

    yP_col      = cube_to_ball_single(Y, pwData);
    yI_face_col = zeros(nI, 1);

    % 2) Face coupling part
    for iface = 1:numel(faceData)
        F = faceData{iface};

        % IGA trial traces
        uT = F.TI * xIi;
        uN = F.GI * xIi;

        % PW trial traces
        vT = face_forward_value_single(X, F, pwData, false);
        vN = face_forward_value_single(X, F, pwData, true);

        % Mixed block B xI
        yP_col = yP_col ...
            + 0.25 * face_adjoint_value_single (F.w .* uN, F, pwData) ...
            - 0.25 * face_adjoint_normal_single(F.w .* uT, F, pwData) ...
            - sigma * face_adjoint_value_single (F.w .* uT, F, pwData);

        % Hermitian counterpart B^* xP
        yI_face_col = yI_face_col ...
            + 0.25 * F.GI' * (F.w .* vT) ...
            - 0.25 * F.TI' * (F.w .* vN) ...
            - sigma * F.TI' * (F.w .* vT);

        % PW-PW face block D xP
        yP_col = yP_col ...
            + 0.25 * face_adjoint_value_single (F.w .* vN, F, pwData) ...
            + 0.25 * face_adjoint_normal_single(F.w .* vT, F, pwData) ...
            + sigma * face_adjoint_value_single (F.w .* vT, F, pwData);
    end

    yP(:,ib)      = yP_col;
    yI_face(:,ib) = yI_face_col;
end

end

function X = ball_to_cube_single(x, pwData)
%Compute to cube single.
X = zeros(pwData.K, pwData.K, pwData.K);
X(pwData.linBall) = x;
end

function x = cube_to_ball_single(X, pwData)
%Compute to ball single.
x = X(pwData.linBall);
end

function Y = conv3_same_single(X, KerFFT, convN, K)
%Compute same single.
Xp = ifftshift(embed_center_single(X, convN));
Yf = ifftn(fftn(Xp) .* KerFFT);
Yf = fftshift(Yf);
Y  = extract_center_single(Yf, K);
end

function A = embed_center_single(X, convN)
%Compute center single.
A = zeros(convN, convN, convN);
K = size(X,1);
i0 = floor((convN - K)/2) + 1;
A(i0:i0+K-1, i0:i0+K-1, i0:i0+K-1) = X;
end

function X = extract_center_single(A, K)
%Extract center single.
i0 = floor((size(A,1) - K)/2) + 1;
X = A(i0:i0+K-1, i0:i0+K-1, i0:i0+K-1);
end

function val = face_forward_value_single(X, F, pwData, withNormal)
%Compute forward value single.

kvals = -pwData.N : pwData.N;
alpha = pwData.alpha;

switch F.type
    case 'x'
        phase_fix = exp(1i * alpha * kvals * F.fixedCoord);
        if withNormal
            phase_fix = (1i * alpha * F.normalSign * kvals) .* phase_fix;
        end
        A = squeeze(sum(reshape(phase_fix, [],1,1) .* X, 1));
        V = F.E1 * A * F.E2.' / sqrt(pwData.Omega);

    case 'y'
        phase_fix = exp(1i * alpha * kvals * F.fixedCoord);
        if withNormal
            phase_fix = (1i * alpha * F.normalSign * kvals) .* phase_fix;
        end
        A = squeeze(sum(reshape(phase_fix, 1,[],1) .* X, 2));
        V = F.E1 * A * F.E2.' / sqrt(pwData.Omega);

    case 'z'
        phase_fix = exp(1i * alpha * kvals * F.fixedCoord);
        if withNormal
            phase_fix = (1i * alpha * F.normalSign * kvals) .* phase_fix;
        end
        A = squeeze(sum(reshape(phase_fix, 1,1,[]) .* X, 3));
        V = F.E1 * A * F.E2.' / sqrt(pwData.Omega);

    otherwise
        error('Unknown face type.');
end

val = V(:);
end

function yP = face_adjoint_value_single(g, F, pwData)
%Accumulate one face value contribution.
R = reshape(g, F.nq1, F.nq2);
Aadj = F.E1' * R * conj(F.E2) / sqrt(pwData.Omega);

kvals = -pwData.N : pwData.N;
alpha = pwData.alpha;

switch F.type
    case 'x'
        phase_fix = exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, [],1,1) .* reshape(Aadj, 1,pwData.K,pwData.K);

    case 'y'
        phase_fix = exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, 1,[],1) .* reshape(Aadj, pwData.K,1,pwData.K);

    case 'z'
        phase_fix = exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, 1,1,[]) .* reshape(Aadj, pwData.K,pwData.K,1);

    otherwise
        error('Unknown face type.');
end

yP = cube_to_ball_single(Xadj, pwData);
end

function yP = face_adjoint_normal_single(g, F, pwData)
%Accumulate one face normal contribution.
R = reshape(g, F.nq1, F.nq2);
Aadj = F.E1' * R * conj(F.E2) / sqrt(pwData.Omega);

kvals = -pwData.N : pwData.N;
alpha = pwData.alpha;

switch F.type
    case 'x'
        phase_fix = (-1i * alpha * F.normalSign * kvals) .* ...
            exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, [],1,1) .* reshape(Aadj, 1,pwData.K,pwData.K);

    case 'y'
        phase_fix = (-1i * alpha * F.normalSign * kvals) .* ...
            exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, 1,[],1) .* reshape(Aadj, pwData.K,1,pwData.K);

    case 'z'
        phase_fix = (-1i * alpha * F.normalSign * kvals) .* ...
            exp(-1i * alpha * kvals * F.fixedCoord);
        Xadj = reshape(phase_fix, 1,1,[]) .* reshape(Aadj, pwData.K,pwData.K,1);

    otherwise
        error('Unknown face type.');
end

yP = cube_to_ball_single(Xadj, pwData);
end
