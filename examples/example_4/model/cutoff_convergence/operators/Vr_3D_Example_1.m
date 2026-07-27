function Vr = Vr_3D_Example_1(p, L, n_pw, x, y, z, nuclear_charge)
% Evaluate the 3-D Example 1 potential at one point.
arguments
    p
    L
    n_pw
    x
    y
    z
    nuclear_charge = 1
end

% scaled reciprocal vectors G
G = p * (2*pi / L);

r = [x; y; z];
r_norm = norm(r);

alpha = 5;
Omega = L^3;

% Evaluate the real-space Ewald term away from its singular center.
assert(r_norm > 0, 'The Ewald potential is singular at the origin.');

real_part = -erfc(alpha * r_norm) / r_norm;

% ---------------- reciprocal-space term (3D corrected) ----------------
s = 0;
for i = 1:n_pw
    G_norm2 = G(i,1)^2 + G(i,2)^2 + G(i,3)^2;
    if G_norm2 > 0
        phase_i = exp(1i * (G(i,:) * r));
        s = s + exp(-G_norm2 / (4 * alpha^2)) * phase_i / G_norm2;
    end
end

recip_part = -(4*pi / Omega) * s;

% ---------------- additive constant ----------------
% The formula includes the periodic Ewald background constant.
C0 = 2*alpha/sqrt(pi);
% final value
Vr = real(nuclear_charge * (real_part + recip_part + C0));

end
