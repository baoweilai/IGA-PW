% integration y=\int_0^R chi_i(r)*chi_j(r) dr
function y = int_overlap(i,j,R,n_integral)
% Integrate the radial overlap term.
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (r_basis(i,xl,R)*r_basis(j,xl,R) + r_basis(i,xr,R)*r_basis(j,xr,R));
end
y=s;
return;
