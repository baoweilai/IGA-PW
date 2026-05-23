function prec = build_interface_block_prec_3D(Ag, gamma, eta, diag_eta, n_total, opts)
%Build the 3-D interface-block preconditioner.
arguments
    Ag
    gamma
    eta
    diag_eta
    n_total
    opts struct
end

eps_iface = opts.eps_iface;

gamma = gamma(:);
eta = eta(:);
diag_eta = diag_eta(:);

diag_eta(abs(diag_eta) < eps_iface) = 1;
Ag = 0.5 * (Ag + Ag');
delta = eps_iface * max(1, norm(Ag, 1));
Areg = Ag + delta * speye(size(Ag, 1));

[solve_gamma, factor_type] = factorize_gamma_block(Areg);

prec = struct();
prec.gamma = gamma;
prec.eta = eta;
prec.diag_eta = diag_eta;
prec.dinv_eta = 1 ./ diag_eta;
prec.delta = delta;
prec.factor_type = factor_type;
prec.solve_gamma = solve_gamma;
prec.n_total = n_total;
prec.apply = @(x) apply_interface_block_prec_local(x, gamma, eta, prec.dinv_eta, solve_gamma, factor_type);
end

function y = apply_interface_block_prec_local(x, gamma, eta, dinv_eta, solve_gamma, factor_type)
%Apply the interface-block preconditioner.
y = zeros(size(x));
y(eta, :) = bsxfun(@times, x(eta, :), dinv_eta);
y(gamma, :) = solve_gamma_block(solve_gamma, factor_type, x(gamma, :));
end

function [solve_gamma, factor_type] = factorize_gamma_block(Areg)
%Factor the interface block.
[R, pflag] = chol(Areg);
if pflag == 0
    solve_gamma = struct('R', R);
    factor_type = 'chol';
else
    [L, D, p] = ldl(Areg, 'vector');
    solve_gamma = struct('L', L, 'D', D, 'p', p);
    factor_type = 'ldl';
end
end

function y = solve_gamma_block(solve_gamma, factor_type, x)
%Apply the interface-block factorization.
switch factor_type
    case 'chol'
        y = solve_gamma.R \ (solve_gamma.R' \ x);

    case 'ldl'
        xp = x(solve_gamma.p, :);
        yp = solve_gamma.L \ xp;
        zp = solve_gamma.D \ yp;
        wp = solve_gamma.L' \ zp;
        y = zeros(size(x));
        y(solve_gamma.p, :) = wp;

    otherwise
        error('Unsupported factor_type: %s', factor_type);
end
end
