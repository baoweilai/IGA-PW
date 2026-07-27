function [Ubar,Vbar,dof]=IGAknotRefineSurface(knotU,knotV,pu,pv,Refinement)
% Refine NURBS surface knot vectors.
Ubar = knotU;
Vbar = knotV;
for i=1:Refinement
    UBreks=unique(Ubar);VBreks=unique(Vbar);
    Xu=(UBreks(1:end-1)+UBreks(2:end))/2;
    Xv=(VBreks(1:end-1)+VBreks(2:end))/2;
    temp=[Ubar,Xu];Ubar=temp;Ubar=sort(Ubar); % Ubar=[Ubar,Xu];Ubar=sort(Ubar);
    temp=[Vbar,Xv];Vbar=temp;Vbar=sort(Vbar);% Vbar=[Vbar,Xv]; Vbar=sort(Vbar);
end

mu=length(Ubar)-pu-1;nv=length(Vbar)-pv-1;
dof=mu*nv;
