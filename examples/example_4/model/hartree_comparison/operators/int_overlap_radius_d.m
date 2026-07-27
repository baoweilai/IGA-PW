% integration y=\int_0^R r^2*chi'_i(r)*chi'_j(r) dr
function y = int_overlap_radius_d(i,j,R,n_integral)
% Integrate the radial derivative overlap.
N = n_integral; % integration discretization
h = R/N;
s = 0;
for k=1:N
    xl = (k-1)*h;
    xr = k*h;
    s = s + 0.5 * h * (xl^2*r_basis_d(i,xl,R)*r_basis_d(j,xl,R) + xr^2*r_basis_d(i,xr,R)*r_basis_d(j,xr,R));
end
y=s;
return;
