function [knotU,knotV,knotW] =  updateKnotVectors_3D(nurbs_original)
%Update 3-D knot vectors after refinement.

knotU = nurbs_original.knotU;
knotV = nurbs_original.knotV;
knotW = nurbs_original.knotW;


% nurbs_original.Lengths

a = nurbs_original.Lengths;

dx = min(a);

n = a/dx;

n = fix(n);


if n(1)>1
    n(1) = fix( n(1)/2 );
    inner_knots = linspace(0,1,n(1) + 1 );
    inner_knots([1,end]) = [];
    knotU_refine = [knotU,inner_knots];
    knotU        = sort(knotU_refine);
end


if n(2)>1
    n(2) = fix( n(2)/2 );
    inner_knots = linspace(0,1,n(2) + 1 );
    inner_knots([1,end]) = [];
    knotV_refine = [knotV,inner_knots];
    knotV        = sort(knotV_refine);
end


if n(3)>1
    n(3) = fix( n(3)/2 );
    inner_knots = linspace(0,1,n(3) + 1 );
    inner_knots([1,end]) = [];
    knotW_refine = [knotW,inner_knots];
    knotW        = sort(knotW_refine);
end


end
