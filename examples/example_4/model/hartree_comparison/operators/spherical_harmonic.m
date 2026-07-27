function [real, imag] = spherical_harmonic(l,m,theta,phi)
% Evaluate the real and imaginary parts of a spherical harmonic.

pre_part = sqrt(factorial(l-m)/factorial(l+m));
if m >= 0
    legendre_part = legendre(l,cos(theta));
    real = pre_part * legendre_part(m+1) * cos(m*phi);
    imag = pre_part * legendre_part(m+1) * sin(m*phi);
else
    legendre_part = (-1)^m * factorial(l+m)/factorial(l-m) * legendre(l,cos(theta));
    real = pre_part * legendre_part(-m+1) * cos(m*phi);
    imag = pre_part * legendre_part(-m+1) * sin(m*phi);
end

real = sqrt((2*l+1)/(4*pi)) * real;
imag = sqrt((2*l+1)/(4*pi)) * imag;
