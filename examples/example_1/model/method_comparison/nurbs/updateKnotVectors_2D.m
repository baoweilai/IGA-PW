function [knotU,knotV] =  updateKnotVectors_2D(knotU,knotV,a,b)
%Update 2-D knot vectors after refinement.
if b> a
n = fix(b/a);
else
    n = fix(a/b);
end

inner_knots = linspace(0,1,n+1);
inner_knots([1,end]) = [];

if a>b
     knotU_refine = [knotU,inner_knots];
     knotU_refine = sort(knotU_refine);
     knotV_refine = knotV;
else
    knotV_refine = [knotV,inner_knots];
    knotV_refine = sort(knotV_refine);
    knotU_refine = knotU;
end

knotU = knotU_refine;
knotV = knotV_refine;

end
