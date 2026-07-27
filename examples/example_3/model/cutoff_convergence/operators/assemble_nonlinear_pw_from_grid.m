function N_pw = assemble_nonlinear_pw_from_grid(k_pw, rho_pw_grid, ~)
% Assemble the plane-wave nonlinear operator from sampled density.

m = size(rho_pw_grid, 1);
n_pw = size(k_pw, 1);

rho_fft = ifftn(rho_pw_grid);
N_pw = zeros(n_pw, n_pw);

for ii = 1:n_pw
    ki = k_pw(ii,:);

    for jj = 1:n_pw
        kj = k_pw(jj,:);
        dk = kj - ki;

        i1 = mod(dk(1), m) + 1;
        i2 = mod(dk(2), m) + 1;

        phase = exp(1i*pi*dk(1)*(1 + 1/m)) * exp(1i*pi*dk(2)*(1 + 1/m));
        N_pw(ii,jj) = phase * rho_fft(i1, i2);
    end
end

N_pw = sparse(0.5 * (N_pw + N_pw'));
end
