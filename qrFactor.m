#Michael Henke
#MHC34

function [Q,R] = qrFactor(A)

n=size(A,1); 
Q=eye(n); 
R=A; 
I=eye(n);

for j=1:n-1
  x=R(j:n,j);
  v=-sign(x(1))*norm(x)*eye(n-j+1,1)-x;
  if norm(v)>0,
    v=v/norm(v);
    P=I; 
    P(j:n,j:n)=P(j:n,j:n)-2*v*v';
    R=P*R;
    Q=Q*P;
   end
end

#R is very close to 0

endfunction