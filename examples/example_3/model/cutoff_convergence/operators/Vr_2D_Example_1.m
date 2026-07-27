function Vr = Vr_2D_Example_1(p, L, n_pw, x, y)
% Evaluate the Example 1 radial potential.

% Ewald-type potential in 2D Example 1
% Robust version: regularize the singular term near r = 0

p = p * (2*pi / L);
r = [x; y];
alpha = 5;

r_norm = norm(r);
r_safe = max(r_norm, 1e-12);

s = 0;
for i = 1:n_pw
    p_i = p(i,:);
    p_norm = norm(p_i);

    if p_norm > 0
        s = s + erfc(p_norm / (2*alpha)) * exp(1i * (p_i * r)) / p_norm;
    end
end

Vr = -erfc(alpha * r_safe) / r_safe ...
    - (2*pi / (L*L)) * s ...
    + 2 * alpha / sqrt(pi);

Vr = real(Vr);

end
