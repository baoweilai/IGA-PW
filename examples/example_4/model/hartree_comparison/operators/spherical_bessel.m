% spherical Bellsel functions
function y = spherical_bessel(n,z)
% Evaluate a spherical Bessel function.
y = besselj(n+1/2,z)*sqrt(pi/(2*z));
return;
