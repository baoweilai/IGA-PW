function  b=Bin(n,i)
%Compute a binomial coefficient.
numerator=1;
denominator=1;

for j=0:i-1
numerator=numerator*(n-j);
denominator=denominator*(i-j);
end
b=numerator/denominator;
