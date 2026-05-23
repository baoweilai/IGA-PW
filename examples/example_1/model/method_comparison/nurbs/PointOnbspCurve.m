function C=PointOnbspCurve(P,U,p,u)
%Evaluate a point on a B-spline curve.
i=findspan(U,p,u);
N=bsplinebasis(U,p,u);
%----------------- C=zeros(ndim,1);
C=P(:,i-p:i)*N;


%% =======Test================
% P(:,1)=[0;1];P(:,2)=[1;1];P(:,3)=[1;0];
