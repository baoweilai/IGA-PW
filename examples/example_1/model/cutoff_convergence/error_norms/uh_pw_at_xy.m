function uh_pw_xy = uh_pw_at_xy(uh_pw,p,x,y,L)
%Evaluate the plane-wave solution at 2-D points.

n_pw_basis = size(p,1);
F = [x;y];
uh_pw_xy = 0;

for k = 1:n_pw_basis
    uh_pw_xy = uh_pw_xy + uh_pw(k)*exp(1i*2*pi/L*p(k,:)*F);
end
end
