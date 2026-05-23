function Vr = Vr_2D_Example_2(p,L,n_pw,x,y)
%Evaluate the Example 2 radial potential.

p      = p*2*pi/L;
pnt_1 = [-1;0];
pnt_2 = [1;0];
r1      = [x;y] - pnt_1;
r1_norm = norm(r1);
r2      = [x;y] - pnt_2;
r2_norm = norm(r2);
alpha  = 5;
s1 = 0;
s2 = 0;
for i = 1:n_pw
    p_norm = norm(p(i,:));
    if p_norm>0
        s1 = s1 + erfc(p_norm/2/alpha)* exp( 1i*p(i,:)*r1 )/p_norm;
        s2 = s2 + erfc(p_norm/2/alpha)* exp( 1i*p(i,:)*r2 )/p_norm;
    end
end

Vr_1 = -2*pi*s1/L/L - erfc(alpha*r1_norm)/r1_norm + 2*alpha/sqrt(pi);
Vr_2 = -2*pi*s2/L/L - erfc(alpha*r2_norm)/r2_norm + 2*alpha/sqrt(pi);
Vr   = Vr_1 + Vr_2;

end
