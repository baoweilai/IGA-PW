function y = r_basis_d(i,r,R)
% Evaluate the derivative of one radial basis function.
% Apply the radius-dependent normalization.
if i == 1
    y = 0;
else
    y = (i-1)*r^(i-2) *i/(R^i);
end
return;
