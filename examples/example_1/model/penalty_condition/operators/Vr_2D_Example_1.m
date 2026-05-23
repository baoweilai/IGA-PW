function Vr = Vr_2D_Example_1(p,L,n_pw,x,y)
%Evaluate the Example 1 radial potential.
p      = p*2*pi/L;
r      = [x;y];
r_norm = norm(r);
alpha  = 5;

s = 0;
for i=1:n_pw
    p_norm = norm(p(i,:));
    if p_norm>0
        s = s + erfc(p_norm/2/alpha)* exp( 1i*p(i,:)*r )/p_norm;
    end
end
Vr = - erfc(alpha*r_norm)/r_norm -2*pi*s/L/L + 2*alpha/sqrt(pi);
end
