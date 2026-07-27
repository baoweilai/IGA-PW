function pwDataNew = fftcheb3_hartree(cfg, Vgrid, innerDomains, pwDataStatic)
% Update the PW Hartree block with the current tensorized FFT--Chebyshev path.

assert(cfg.fftGridN == 300, 'The PW FFT grid must remain 300.');
assert(cfg.chebDegree == 80, 'The inner Chebyshev degree must remain 80.');
assert(size(Vgrid, 1) == 300 && size(Vgrid, 2) == 300 ...
    && size(Vgrid, 3) == 300, 'The Hartree grid must be 300-by-300-by-300.');

L = cfg.Lcell;
mFFT = size(Vgrid, 1);
qMax = 2 * pwDataStatic.N;
qvals = -qMax:qMax;
qCount = numel(qvals);

Vraw = fftn(Vgrid) / (mFFT ^ 3);
freqIdx = mod(qvals, mFFT) + 1;
phase = exp(1i * pi * qvals * (1 - 1 / mFFT));
Vfull = Vraw(freqIdx, freqIdx, freqIdx) ...
    .* reshape(phase, qCount, 1, 1) ...
    .* reshape(phase, 1, qCount, 1) ...
    .* reshape(phase, 1, 1, qCount);
clear Vraw

innerDomains = reshape(innerDomains, 1, 6);
a = 0.5 * (innerDomains(2) - innerDomains(1));
VinInner = inner_cheb_local( ...
    L, a, qvals, cfg.chebDegree, ...
    @(X, Y, Z) sample_grid_local(Vgrid, L, X, Y, Z));

Vdelta = hermitize_kernel_local(Vfull - VinInner);
pwDataNew = pwDataStatic;
pwDataNew.Vker = hermitize_kernel_local(pwDataStatic.Vker + Vdelta);
pwDataNew.VinInner = VinInner;
pwDataNew.VinCheb = VinInner;
pwDataNew.Vfft = fftn(ifftshift(embed_center_local( ...
    pwDataNew.Vker, pwDataStatic.convN)));

V0delta = real(Vdelta(qMax + 1, qMax + 1, qMax + 1));
pwDataNew.stiffDiag = pwDataStatic.stiffDiag + V0delta;
end

function Vsample = sample_grid_local(Vgrid, L, X, Y, Z)
% Interpolate a periodic grid at physical points.
mFFT = size(Vgrid, 1);
dx = L / mFFT;
xmid = -L / 2 + dx / 2 + (0:mFFT-1) * dx;
Vsample = interpn(xmid, xmid, xmid, Vgrid, X, Y, Z, 'linear');
end

function VinInner = inner_cheb_local(L, a, qvals, nCheb, sampleFun)
% Compute the inner-domain Chebyshev correction.
waveScale = 2 * pi / L;
Omega = L ^ 3;
angles = pi * ((0:nCheb-1) + 0.5) / nCheb;
nodes = a * cos(angles);
[X, Y, Z] = ndgrid(nodes, nodes, nodes);
Vg = sampleFun(X, Y, Z);

degrees = 0:nCheb-1;
Phi = cos(angles(:) * degrees);
scale = (2 / nCheb) * ones(nCheb, 1);
scale(1) = 1 / nCheb;
projection = bsxfun(@times, Phi.', scale);
C = mode_product_local(Vg, projection, 1);
C = mode_product_local(C, projection, 2);
C = mode_product_local(C, projection, 3);

nQuad = max(240, 4 * nCheb);
[gp, gw] = grule(nQuad);
xq = a * gp(:);
wq = a * gw(:);
quadAngles = acos(max(min(gp(:), 1), -1));
Tquad = cos((0:nCheb-1).' * quadAngles.');
phase = exp(-1i * waveScale * (xq * qvals(:).'));
Q = Tquad * bsxfun(@times, wq, phase);

VinInner = mode_product_local(C, Q.', 1);
VinInner = mode_product_local(VinInner, Q.', 2);
VinInner = mode_product_local(VinInner, Q.', 3) / Omega;
end

function Y = mode_product_local(X, A, dim)
% Multiply a tensor along one mode.
order = [dim, 1:dim-1, dim+1:ndims(X)];
Xp = permute(X, order);
shape = size(Xp);
Xmat = reshape(Xp, shape(1), []);
Ymat = A * Xmat;
Yp = reshape(Ymat, [size(A, 1), shape(2:end)]);
Y = ipermute(Yp, order);
end

function A = embed_center_local(X, targetN)
% Embed an array in the center of a larger grid.
A = zeros(targetN, targetN, targetN);
sourceN = size(X, 1);
first = floor((targetN - sourceN) / 2) + 1;
A(first:first+sourceN-1, first:first+sourceN-1, first:first+sourceN-1) = X;
end

function Kherm = hermitize_kernel_local(Kraw)
% Enforce Hermitian symmetry in a Fourier kernel.
Kherm = 0.5 * (Kraw + conj(Kraw(end:-1:1, end:-1:1, end:-1:1)));
end
