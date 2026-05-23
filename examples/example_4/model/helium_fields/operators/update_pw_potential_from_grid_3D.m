function pwDataNew = update_pw_potential_from_grid_3D( ...
Vgrid, L, inner_domains, pwDataStatic, opts)
%Update the PW potential operator from a grid potential.
arguments
Vgrid
L
inner_domains
pwDataStatic
opts struct
end

inner_cheb_n = opts.inner_cheb_n;
combine_with_static = opts.combine_with_static;

mFFT = size(Vgrid, 1);
qN = 2 * pwDataStatic.N;
Kq = 2 * qN + 1;
qvals = -qN:qN;

% 1) full-domain Fourier coefficients from midpoint FFT grid
Vraw = fftn(Vgrid) / (mFFT ^ 3);

freqIdx = mod(qvals, mFFT) + 1;
phase = exp(1i * pi * qvals * (1 - 1 / mFFT));
Vfull = Vraw(freqIdx, freqIdx, freqIdx) ...
    .* reshape(phase, Kq, 1, 1) ...
    .* reshape(phase, 1, Kq, 1) ...
    .* reshape(phase, 1, 1, Kq);

% 2) inner-domain coefficients by Chebyshev correction
VinInner = build_inner_chebyshev_correction_3D( ...
    L, inner_domains, qvals, inner_cheb_n, ...
    @(X, Y, Z) sample_grid_values(Vgrid, L, X, Y, Z));

Vdelta = Vfull - VinInner;
Vdelta = hermitize_kernel_3d(Vdelta);

% 3) update pwData
pwDataNew = pwDataStatic;
if combine_with_static
    pwDataNew.Vker = hermitize_kernel_3d(pwDataStatic.Vker + Vdelta);
else
    pwDataNew.Vker = Vdelta;
end

pwDataNew.VinInner = VinInner;
pwDataNew.VinCheb = VinInner;

pwDataNew.Vfft = fftn(ifftshift(embed_center(pwDataNew.Vker, pwDataStatic.convN)));
pwDataNew.inner_cheb_n = inner_cheb_n;

V0delta = real(Vdelta(qN + 1, qN + 1, qN + 1));
if combine_with_static
    pwDataNew.stiffDiag = pwDataStatic.stiffDiag + V0delta;
else
    k2 = pwDataStatic.k_pw(:, 1).^2 + pwDataStatic.k_pw(:, 2).^2 + pwDataStatic.k_pw(:, 3).^2;
    mass0 = pwDataStatic.massDiag(1);
    pwDataNew.stiffDiag = 0.5 * pwDataStatic.alpha ^ 2 * k2 * mass0 + real(pwDataNew.Vker(qN + 1, qN + 1, qN + 1));
end
end

function Vsample = sample_grid_values(Vgrid, L, X, Y, Z)
%Sample grid values.
mFFT = size(Vgrid, 1);
dx = L / mFFT;
x1d = -L / 2 + dx / 2 + (0:mFFT-1) * dx;
Vsample = interpn(x1d, x1d, x1d, Vgrid, X, Y, Z, 'linear');
end

function A = embed_center(X, convN)
%Embed a smaller array at the grid center.
A = zeros(convN, convN, convN);
K = size(X, 1);
i0 = floor((convN - K) / 2) + 1;
A(i0:i0+K-1, i0:i0+K-1, i0:i0+K-1) = X;
end

function Ksym = hermitize_kernel_3d(K)
%Compute kernel 3D.
Ksym = 0.5 * (K + conj(K(end:-1:1, end:-1:1, end:-1:1)));
end
