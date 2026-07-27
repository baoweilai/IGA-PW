function y = spherical_harmonic_xyz(l,m,x,y,z)
% Evaluate a spherical harmonic from Cartesian coordinates.

% Convert the Cartesian point to spherical angles.
r = (x^2+y^2+z^2)^(1/2);
if r == 0
    error('error r=0');
else
    theta = acos(z/r);
    if y==0
        if x>=0
            phi = pi/2;
        else
            phi = -pi/2;
        end
    elseif y>0
        phi = atan(x/y);
    else
        if x>=0
            phi = atan(x/y)+pi;
        else
            phi = atan(x/y)-pi;
        end
    end
    [real, imag] = spherical_harmonic(l,m,theta,phi);
end
% Combine the real and imaginary components.
y = real + 1i*imag;
