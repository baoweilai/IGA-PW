function [Ubar,Vbar]=IGADegreeElevSurface(U,V,t)
% Elevate a NURBS surface degree.

Ubar = [zeros(1,t),U,ones(1,t)];

Vbar = [zeros(1,t),V,ones(1,t)];
end
