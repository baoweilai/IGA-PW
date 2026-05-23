function [knotUbar,knotVbar,knotWbar]=IGADegreeElevVolume(knotU,knotV,knotW,t)
%Elevate a NURBS volume degree.

% Elevate the degree of B-spline basis functions to (pi +t) for i=1,2 and


if t==0 % t=0 means that we do not elevate the degree of B-splines.
    knotUbar=knotU;knotVbar=knotV; knotWbar=knotW;
end

if t>=1

    knotUbar = [zeros(1,t),knotU,ones(1,t)];
    knotVbar = [zeros(1,t),knotV,ones(1,t)];
    knotWbar = [zeros(1,t),knotW,ones(1,t)];

end


end
