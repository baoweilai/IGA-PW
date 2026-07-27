function y = r_basis(i,r,R)
% Evaluate one radial basis function.
% Apply the radius-dependent normalization.
y = r^(i-1) *i/(R^i);
return;
