function [W,DW,C,DC]=NurbsCurve(P,U,w,p,u)
%Evaluate a NURBS curve.
[ndim,~]=size(P);
Pw=WightedConPtsCurve(P,w);
Cw=PointOnbspCurve(Pw,U,p,u);
W=Cw(end);
C=Project(Cw);
Cwder=bspCurveder(Pw,U,p,u);
Ader=Cwder(1:ndim);
DW=Cwder(end);
DC=(Ader-C*DW)/W;
end

%==================Test====================================================
%===== The Nurbs Book P126 Ex4.2;
%    0 1 1];u=0;
