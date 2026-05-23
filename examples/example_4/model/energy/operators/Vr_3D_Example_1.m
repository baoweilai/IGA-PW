function Vr = Vr_3D_Example_1(p, L, n_pw, x, y, z, nuclear_charge)
%Compute 3D example 1.
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

% ---------------- real-space singular term ----------------
% In practical DG-IGA-PW quadrature, r = 0 should be avoided by the mesh /
% quadrature and/or patched potential usage on the PW side.
% We keep a tiny guard only to avoid literal division-by-zero crash.
if r_norm < 1e-14
    r_norm = 1e-14;
end

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
% Standard 3D periodic Ewald-style background constant.
% If you want a different convention later, this is the place to change.
C0 = 2*alpha/sqrt(pi);
% final value
Vr = real(nuclear_charge * (real_part + recip_part + C0));

end
