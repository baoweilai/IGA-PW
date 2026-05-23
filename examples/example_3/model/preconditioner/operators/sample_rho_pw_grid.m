function [rho_pw_grid, dx] = sample_rho_pw_grid(cA, k_pw, L, a, m)
%Sample the plane-wave density on a grid.

dx = L / m;
Omega_area = L * L;

C = complex(zeros(m, m));
n_pw = size(k_pw, 1);

for s = 1:n_pw
    k1 = k_pw(s,1);
    k2 = k_pw(s,2);

    i1 = mod(k1, m) + 1;
    i2 = mod(k2, m) + 1;

    phase_mid = (-1)^(k1 + k2) * exp(1i * pi * (k1 + k2) / m);
    C(i1, i2) = C(i1, i2) + cA(s) * phase_mid;
end

u_grid = (m * m / sqrt(Omega_area)) * ifftn(C);

x1d = -L/2 + dx/2 + (0:m-1) * dx;
[X, Y] = ndgrid(x1d, x1d);

mask_out = (abs(X) > a) | (abs(Y) > a);

rho_pw_grid = abs(u_grid).^2;
rho_pw_grid(~mask_out) = 0;
rho_pw_grid = real(rho_pw_grid);
end
