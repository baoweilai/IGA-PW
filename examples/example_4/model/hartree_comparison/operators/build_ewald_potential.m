function [vext,Vinner] = build_ewald_potential( ...
    L,R,N,n_r,Lm,mu,charge)
% BUILD_EWALD_POTENTIAL Assemble the DG-PW Ewald potential.
% The plane-wave and radial-spherical blocks use the same quadrature nodes,
% trapezoidal weights, reciprocal cutoff, and angular expansion.

Omega = L^3;
waveScale = 2*pi/L;
constantShift = 2*charge*mu/sqrt(pi);
unitModes = [ ...
    -1, 0, 0;
     1, 0, 0;
     0,-1, 0;
     0, 1, 0;
     0, 0,-1;
     0, 0, 1];
gNorm = waveScale;
recipCoefficient = -charge*(4*pi/Omega) ...
    *exp(-gNorm^2/(4*mu^2))/gNorm^2;

nExtIntegral = 1000;
radialGrid = linspace(0,R,nExtIntegral+1);
radialWeight = (R/nExtIntegral)*ones(size(radialGrid));
radialWeight([1,end]) = 0.5*R/nExtIntegral;
radialIntegrand = -charge*radialGrid.*erfc(mu*radialGrid) ...
    +constantShift*radialGrid.^2;

coordinates = -2*N:2*N;
[qx,qy,qz] = ndgrid(coordinates,coordinates,coordinates);
qSquared = qx.^2+qy.^2+qz.^2;
qNorm = waveScale*sqrt(qSquared);

fullCoefficient = zeros(size(qNorm));
zero = qSquared == 0;
fullCoefficient(zero) = constantShift-charge*pi/(mu^2*Omega);
nonzero = ~zero;
q = qNorm(nonzero);
fullCoefficient(nonzero) = -charge*(4*pi/Omega) ...
    .*(1-exp(-q.^2/(4*mu^2)))./q.^2;
unitMask = (abs(qx)+abs(qy)+abs(qz) == 1) ...
    & (qSquared == 1);
fullCoefficient(unitMask) = fullCoefficient(unitMask)+recipCoefficient;

uniqueSquared = unique(qSquared(:));
radialLookup = zeros(max(uniqueSquared)+1,1);
for index = 1:numel(uniqueSquared)
    squared = uniqueSquared(index);
    normValue = waveScale*sqrt(squared);
    j0 = spherical_bessel_original_array_local(0,normValue*radialGrid);
    radialLookup(squared+1) = 4*pi ...
        *sum(radialWeight.*radialIntegrand.*j0);
end
innerRadial = radialLookup(qSquared+1);

innerReciprocal = zeros(size(qNorm));
for mode = 1:size(unitModes,1)
    combinedNorm = waveScale*sqrt( ...
        (qx+unitModes(mode,1)).^2 ...
        +(qy+unitModes(mode,2)).^2 ...
        +(qz+unitModes(mode,3)).^2);
    innerReciprocal = innerReciprocal ...
        +recipCoefficient*sphere_fourier_array_local(combinedNorm,R);
end
vext = real(fullCoefficient-(innerRadial+innerReciprocal)/Omega);

nAngular = (Lm+1)^2;
nLocal = n_r*nAngular;
Vinner = complex(zeros(nLocal,nLocal));

radialLocal = zeros(n_r,n_r);
radialBessel = zeros(n_r,n_r,2*Lm+1);
for k1 = 1:n_r
    basis1 = k1*radialGrid.^(k1-1)/R^k1;
    for k2 = 1:n_r
        basis2 = k2*radialGrid.^(k2-1)/R^k2;
        basisProduct = basis1.*basis2;
        radialLocal(k1,k2) = sum( ...
            radialWeight.*basisProduct.*radialIntegrand);
        for ell = 0:2*Lm
            radialBessel(k1,k2,ell+1) = sum( ...
                radialWeight.*radialGrid.^2.*basisProduct ...
                .*spherical_bessel_original_array_local( ...
                ell,gNorm*radialGrid));
        end
    end
end

planeWaveAmplitude = cell(2*Lm+1,1);
for ell = 0:2*Lm
    amplitudes = complex(zeros(2*ell+1,1));
    for emm = -ell:ell
        directionSum = 0;
        for mode = 1:size(unitModes,1)
            directionSum = directionSum+conj(spherical_harmonic_xyz( ...
                ell,emm,unitModes(mode,1),unitModes(mode,2),unitModes(mode,3)));
        end
        amplitudes(emm+ell+1) = recipCoefficient*4*pi ...
            *(1i)^ell*directionSum;
    end
    planeWaveAmplitude{ell+1} = amplitudes;
end

for k1 = 1:n_r
    for l1 = 0:Lm
        for m1 = -l1:l1
            p = (k1-1)*nAngular+l1^2+m1+l1+1;
            for k2 = 1:n_r
                for l2 = 0:Lm
                    for m2 = -l2:l2
                        q = (k2-1)*nAngular+l2^2+m2+l2+1;
                        value = 0;
                        if l1 == l2 && m1 == m2
                            value = value+radialLocal(k1,k2);
                        end
                        for ell = abs(l1-l2):2:(l1+l2)
                            emm = m1-m2;
                            if abs(emm) > ell
                                continue;
                            end
                            gaunt = (-1)^m1 ...
                                *sqrt((2*l1+1)*(2*ell+1)*(2*l2+1)/(4*pi)) ...
                                *Wigner3j(l1,ell,l2,0,0,0) ...
                                *Wigner3j(l1,ell,l2,-m1,emm,m2);
                            amplitude = planeWaveAmplitude{ell+1}(emm+ell+1);
                            value = value+amplitude*gaunt ...
                                *radialBessel(k1,k2,ell+1);
                        end
                        Vinner(p,q) = value;
                    end
                end
            end
        end
    end
end
Vinner = 0.5*(Vinner+Vinner');
end

function values = spherical_bessel_original_array_local(ell,z)
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

function values = sphere_fourier_array_local(waveNorm,R)
% Evaluate the Fourier coefficients of a ball indicator.
values = zeros(size(waveNorm));
zero = waveNorm == 0;
values(zero) = 4*pi*R^3/3;
nonzero = ~zero;
k = waveNorm(nonzero);
values(nonzero) = 4*pi*R^2 ...
    .*spherical_bessel_original_array_local(1,k*R)./k;
end
