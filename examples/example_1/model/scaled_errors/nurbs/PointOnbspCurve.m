function C=PointOnbspCurve(P,U,p,u)
% Evaluate a point on a B-spline curve.
i=findspan(U,p,u);
N=bsplinebasis(U,p,u);
% Combine the active control points with the nonzero basis values.
C=P(:,i-p:i)*N;
