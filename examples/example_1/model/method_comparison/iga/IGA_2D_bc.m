function [A,rhs]=IGA_2D_bc(A,rhs,m,n)
%Mark boundary degrees of freedom for the 2-D IGA mesh.
Node_down=1:m;
Node_up=m*(n-1)+(1:m);
Node_left=1+(0:n-1)*m;
Node_right=(1:n)*m;
be=[Node_up,Node_down,Node_left,Node_right];
be=unique(be);
A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
A(be,be)=eye(length(be));
