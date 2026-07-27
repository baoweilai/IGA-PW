function ndu=AllBasisFuns(U,p,u)
% Evaluate all nonzero B-spline basis functions.
i=findspan(U,p,u);

ndu=zeros(p+1,p+1);
ndu(1,1)=1.0;
for j=1:p
    saved=0;
for r=1:j
    ndu(j+1,r)=U(i+r)-U(i-j+r);
     temp=ndu(r,j)/ndu(j+1,r);
ndu(r,j+1)=saved+(U(i+r)-u)*temp;
saved=(u-U(i-j+r))*temp;
end
ndu(j+1,j+1)=saved;
end
