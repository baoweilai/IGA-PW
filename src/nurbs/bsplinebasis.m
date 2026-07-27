function  N=bsplinebasis(U,p,u)
% Evaluate B-spline basis functions.
N=zeros(p+1,1);
N(1)=1.0;
i=findspan(U,p,u);

for j=1:p
    saved=0;
for k=1:j
temp=N(k);
left=U(i+k)-u;  right=u-U(i+k-j);
tmp=temp/(U(i+k)-U(i+k-j));
N(k)=saved+tmp*left;
saved=tmp*right;
end
N(j+1)=saved;
end
end
