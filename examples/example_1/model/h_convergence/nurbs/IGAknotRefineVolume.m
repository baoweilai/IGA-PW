function [Ubar,Vbar,Wbar,dof]=IGAknotRefineVolume(knotU,pu,knotV,pv,knotW,pw,Refinement)
%Refine NURBS volume knot vectors.


Ubar = knotU;

Vbar = knotV;

Wbar = knotW;


for i=1:Refinement
    UBreks=unique(Ubar);VBreks=unique(Vbar); WBreks=unique(Wbar);
    Xu= (UBreks(1:end-1)+UBreks(2:end))/2;
    Xv= (VBreks(1:end-1)+VBreks(2:end))/2;
    Xw=(WBreks(1:end-1)+WBreks(2:end))/2;
    temp=[Ubar,Xu]; Ubar=temp; Ubar=sort(Ubar);
    temp=[Vbar,Xv];  Vbar=temp; Vbar=sort(Vbar);
    temp=[Wbar,Xw];Wbar=temp;Wbar=sort(Wbar);
end

nu=length(Ubar)-pu-1;nv=length(Vbar)-pv-1; nw=length(Wbar)-pw-1;
dof=nu*nv*nw;
