function Cder=bspCurveder(P,U,p,u)
%Evaluate B-spline curve derivatives.
uspan=findspan(U,p,u);
ndim=size(P,1);
temp=P(:,uspan-p:uspan);
Q=zeros(ndim,p);

%% Section
for  i=1:ndim
    Q(i,:)=p*(temp(i,2:end)-temp(i,1:end-1))./(U(uspan+1:uspan+p)-U(uspan-p+1:uspan));
end
Ubar=U;
Ubar([1,end])=[];
Nu=bsplinebasis(Ubar,p-1,u);
Cder=Q*Nu;
end

%------------------Test----------------
