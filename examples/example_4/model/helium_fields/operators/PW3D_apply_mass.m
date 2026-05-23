function y = PW3D_apply_mass(x, pwData)
%Apply the 3-D plane-wave mass operator.

[n_pw, nb] = size(x);
y = zeros(n_pw, nb);

for ib = 1:nb
    X  = ball_to_cube_single(x(:,ib), pwData);
    UX = conv3_same_single(X, pwData.Ufft, pwData.convN, pwData.K);
    Y  = X - UX;
    y(:,ib) = cube_to_ball_single(Y, pwData);
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
